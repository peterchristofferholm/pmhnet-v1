from dataclasses import dataclass, field
from typing import List

import numpy as np
import pandas as pd
import psycopg2.extensions
from psycopg2.extras import NamedTupleCursor
from scipy import stats
from torch.utils.data import Dataset

from pmhnet.utils.dataset_creation import make_survarray
from pmhnet.utils.discretization import compute_breaks


class DataAssembler(Dataset):

    def __init__(self, conn, training=True):

        self.conn = conn
        self.conn.cursor_factory = NamedTupleCursor

        cursor = self.conn.cursor()
        cursor.execute(
            """
            SELECT pid FROM surv_master WHERE train = %s ORDER BY pid
            """,
            (training, )
        )
        self.pids = [r.pid for r in cursor]
        self.datasets = []

    def add_labels(self, n_intervals, max_time):

        cursor = self.conn.cursor()
        cursor.execute(
            """
            SELECT time, event FROM surv_master WHERE pid=ANY(%s) ORDER BY pid
            """,
            (self.pids, )
        )
        surv = np.stack(cursor.fetchall())  # shape: (n_pats, 2)

        # break, labels, and identifiers
        self.breaks = compute_breaks(cursor, max_time, n_intervals)
        self.labels = np.apply_along_axis(
            lambda x: make_survarray(x[0], x[1], self.breaks), 1, surv
        )

    def add_features(self, alias, **kwargs):

        try:
            cls = FeatureCollection.get[alias]
        except KeyError as e:
            avail = FeatureCollection.get
            print(f"Data {fts_name} not implemented, try any of {avail}")
            raise e

        dataset = cls(conn=self.conn, pids=self.pids, **kwargs)
        dataset.fetch_data()  # run queries on dbserver
        self.datasets.append(dataset)

    @property
    def colnames(self):
        return [c for d in self.datasets for c in d.colnames]

    @property
    def n_features(self):

        return sum(len(d.colnames) for d in self.datasets)

    def __getitem__(self, idx):

        label = self.labels[idx]
        features = [dataset[idx] for dataset in self.datasets]
        features = np.concatenate(features)

        return label, features

    def __len__(self):

        return len(self.pids)


@dataclass
class FeatureCollection:

    conn : psycopg2.extensions.connection = field(repr=False)
    pids : List[int] = field(repr=False)  # maintain order

    get = {}

    def __init_subclass__(cls):
        cls.get[cls.alias] = cls

    def fetch_data(self, **kwargs):
        raise NotImplementedError

    @property
    def alias(self):
        raise NotImplementedError

    @property
    def colnames(self):
        return self._colnames


@dataclass
class SimpleDiagnoses(FeatureCollection):
    """One-hot encoded icd10 diagnosis codes"""

    horizon: int     # window to look for diagnoses
    code_level: str  # icd10 hierarchy

    alias = "diag-simple"

    def __post_init__(self):

        cursor = self.conn.cursor()
        cursor.execute(
            """
            SELECT DISTINCT code FROM diagnoses WHERE level = %s AND time > %s
            """,
            (self.code_level, -self.horizon)
        )
        codes = np.concatenate(cursor.fetchall())  # tuples to 1d array
        codes.sort()

        self._codes = codes
        self._colnames = [f"{self.alias}_{col}" for col in codes]

    def fetch_data(self):

        cursor = self.conn.cursor()

        # Get idx corresponding to pid and code
        pid2idx = dict(zip(self.pids, range(len(self.pids))))
        code2idx = dict(zip(self._codes, range(self._codes.shape[0])))

        # Initialize empty np.array
        data = np.zeros(
            (len(self.pids), self._codes.shape[0]), dtype=np.float32
        )

        # Record diagnosis "fingerprint"
        cursor.execute(
            """
            SELECT DISTINCT pid, code FROM diagnoses
            WHERE level = %s AND time > %s AND pid = ANY(%s)
            """,
            (self.code_level, -self.horizon, self.pids)
        )
        for pid, code in cursor:
            data[pid2idx[pid], code2idx[code]] = 1

        self._data = data

    def __getitem__(self, idx):

        return self._data[idx, :]


@dataclass
class ClinicalOne(FeatureCollection):
    """Collection of clinical features corresponding to those also found in
    GRACE. Corresponds to the following columns in the database:

    ClinicalOne (
        age      FLOAT8,      -- age at debut
        pulse    INTEGER,     -- heart rate (bpm)
        sys_bp   INTEGER,     -- systolic blood pressure (mmHg)
        crea     NUMERIC(5),  -- serum creatinine (umol/L)
        arrest   BOOLEAN,     -- cardiac arrest at admission (0/1)
        enzymes  BOOLEAN,     -- increased cardiac enzymes (0/1)
        stemi    BOOLEAN,     -- stemi (0/1)
        killip   INTEGER      -- killip class (1/2/3/4)
    )
    """

    remove: tuple = None

    alias = "clnc-1"


    def fetch_data(self):

        # create cursor and fetch data
        cursor = self.conn.cursor()
        cursor.execute(
            """
            SELECT
                pid, age, pulse, sys_bp, crea, arrest, enzymes, stemi, killip
            FROM baseline
            """
        )

        # convert to pandas dataframe for easy handling
        df = pd.DataFrame(cursor.fetchall())

        # min-max normalization
        df.age = (df.age - 15) / (105 - 15)  # min = 15, max = 105

        # z-score normalization
        df.pulse = (df.pulse
            .transform(lambda x: (x - x.mean()) / x.std())
            .fillna(0)
        )

        # z-score normalization
        df.sys_bp = (df.sys_bp
            .transform(lambda x: (x - x.mean()) / x.std())
            .fillna(0)
        )

        # remove underscores to make one-hot variables stand out
        df.rename(columns={"sys_bp" : "sysbp"}, inplace=True)

        # log then z-score normalization
        df.crea = stats.zscore(
            np.log(df.crea.astype("f4") + 1), nan_policy="omit"
        )
        df.crea = np.nan_to_num(df.crea, nan=0.0)

        # convert enzymes and killip into categories
        df.enzymes = df.enzymes.astype(
            pd.CategoricalDtype(categories=[True, False])
        )
        df.killip = df.killip.astype(
            pd.CategoricalDtype(categories=[1, 2, 3, 4])
        )

        # one-hot-encoding of enzymes and killip
        df = pd.get_dummies(df, columns=["enzymes", "killip"], dummy_na=True)

        # limit dataset to only included pids
        df = df.set_index("pid").reindex(self.pids)

        # remove columns if specified
        if self.remove:
            _drop = lambda x: any(x.startswith(y) for y in self.remove)
            drop = [x for x in df.columns if _drop(x)]
            df.drop(labels=drop, axis="columns", inplace=True)

        # set instance variables
        self._data = df.to_numpy(dtype=np.float32)
        self._colnames = [f"{self.alias}_{col}" for col in df.columns]


    def __getitem__(self, idx):

        return self._data[idx, :]


@dataclass
class ClinicalTwo(FeatureCollection):
    """Collection of clinical features available from pats and other sources
    which doesn't fit in any of the other feature collections. Features part of
    this collection is:

    CREATE TABLE pmhnet.clinicaltwo (
        pid integer,
        sex boolean,          -- true = M, false = F
        height numeric(5,0),  -- body height
        weight numeric(5,0),  -- body weight
        dia_bp integer,       -- diastolic blood pressure (mmHg)
        smoking character,    -- yes/no/ex
        vessels character,    -- DIF/1VD/2VD/3VD
        dominance character,  -- R/L/B
        nyha integer,
        ccs integer,
        familiary_ihd text,
        icd_or_pm boolean,
        lvef text,
        ischemia_test text,
        abnormal_ekg text,
        abnormal_qrs_st text
    )
    """

    remove: tuple = None

    alias = "clnc-2"


    def fetch_data(self):

        # create cursor and fetch data
        cursor = self.conn.cursor()
        cursor.execute(
            """
            SELECT
                pid, sex, height::float, weight::float, dia_bp, smoking,
                vessels, dominance, nyha, ccs, familiary_ihd,
                icd_or_pm, lvef, ischemia_test, abnormal_ekg, abnormal_qrs_st
            FROM clinicaltwo
            """
        )

        # convert to pandas dataframe for easy handling
        df = pd.DataFrame(cursor.fetchall())

        # z-score normalization
        df.height = (df.height
            .transform(lambda x: (x - x.mean()) / x.std())
            .fillna(0)
        )

        # z-score normalization
        df.weight = (df.weight
            .transform(lambda x: (x - x.mean()) / x.std())
            .fillna(0)
        )

        # z-score normalization
        df.dia_bp = (df.dia_bp
            .transform(lambda x: (x - x.mean()) / x.std())
            .fillna(0)
        )

        # convert smoking, vessels, dominance, nyha, ccs into categories for
        # proper one-hot encoding in case of missing categories
        df.smoking = df.smoking.astype(
            pd.CategoricalDtype(categories=["no", "ex", "yes"])
        )
        df.vessels = df.vessels.astype(
            pd.CategoricalDtype(categories=["DIF", "1VD", "2VD", "3VD"])
        )
        df.dominance = df.dominance.astype(
            pd.CategoricalDtype(categories=["B", "R", "L"])
        )
        df.nyha = df.nyha.astype(
            pd.CategoricalDtype(categories=[1, 2, 3, 4])
        )
        df.ccs = df.ccs.astype(
            pd.CategoricalDtype(categories=[0, 1, 2, 3, 4])
        )

        # remove underscores to make one-hot encoded variables stand out
        df.rename(columns=lambda x: x.replace("_", "-"), inplace=True)

        # one-hot encoding of categorical columns with NAs
        df = pd.get_dummies(df, dummy_na=True,
            columns=["smoking", "vessels", "dominance", "nyha", "ccs"]
        )

        # one-hot encoding of categorical columns with explicit missing
        df = pd.get_dummies(
            df, dummy_na=False,  # already present
            columns=[
                "lvef", "ischemia-test", "abnormal-ekg", "abnormal-qrs-st",
                "familiary-ihd"
            ]
        )

        # limit dataset to only included pids
        df = df.set_index("pid").reindex(self.pids)

        # order columns alphabetically
        df = df.reindex(sorted(df.columns), axis=1)

        # remove columns if specified
        if self.remove:
            _drop = lambda x: any(x.startswith(y) for y in self.remove)
            drop = [x for x in df.columns if _drop(x)]
            df.drop(labels=drop, axis="columns", inplace=True)

        # set instance variables
        self._data = df.to_numpy(dtype=np.float32)
        self._colnames = [f"{self.alias}_{col}" for col in df.columns]


    def __getitem__(self, idx):

        return self._data[idx, :]


@dataclass
class SimpleBiochemical(FeatureCollection):
    """Collection of biochemical results. Fetches from a table on the dbserver
    with the following structure:

    CREATE TABLE biochem (
        pid        INTEGER REFERENCES surv_master (pid),
        qid        VARCHAR NOT NULL,  -- quantity identifier
        time       FLOAT8 NOT NULL,   -- time of test relative to index
        value      INTEGER NOT NULL   -- flagging system (-1/0/1)
    )
    """

    horizon: int       # Window to look for results
    pattern: str = ""  # Regex pattern to filter qids

    alias = "bioc-1"


    def __post_init__(self):

        cursor = self.conn.cursor()
        cursor.execute(
            """
            SELECT DISTINCT qid FROM biochem WHERE qid ~ %s
            """,
            (self.pattern, )
        )
        qids = np.concatenate(cursor.fetchall())  # tuples to 1d array
        qids.sort()

        self._qids = qids
        self._colnames = [
            f"{self.alias}_{col}_{val}"
                for col in qids for val in ["low", "normal", "high"]
        ]


    def fetch_data(self):

        # Get idx corresponding to pids and codes
        pid2idx = dict(zip(self.pids, range(len(self.pids))))
        qid2idx = dict(zip(self._qids, range(len(self._qids))))

        # Init empty array with shape [n_pid, n_qid, n_cols]
        # Variables are (-1, 0, 1) -> n_cols is 3
        data = np.zeros((len(self.pids), len(self._qids), 3), dtype=np.float32)

        # Create cursor and fetch data
        cursor = self.conn.cursor()
        cursor.execute(
            """
            SELECT pid, qid, value + 1 as value FROM biochem
                WHERE pid = ANY(%s) and time <= %s and qid ~ %s;
            """,
            (self.pids, self.horizon, self.pattern)
        )

        mat = np.eye(3)  # vars: {low, normal, high}
        for row in cursor:
            data[pid2idx[row.pid], qid2idx[row.qid], :] = mat[row.value]

        self._data = data


    def __getitem__(self, idx):

        return self._data[idx].reshape(-1)


@dataclass
class ProcedureCodes(FeatureCollection):
    """Collection of procedure and investigations codes. The codes are binary
    encoded indicating if subjects have received the given code.
    """

    horizon: int  # window to look for codes

    alias = "proc-1"

    def __post_init__(self):

        cursor = self.conn.cursor()

        cursor.execute("SELECT DISTINCT code FROM nomesco")
        codes = np.concatenate(cursor.fetchall())
        codes.sort()  # alphabetically

        self._codes = codes
        self._colnames = [f"{self.alias}_{col}" for col in codes]

    def fetch_data(self):

        # Get idx corresponding to pid and code
        pid2idx = dict(zip(self.pids, range(len(self.pids))))
        code2idx = dict(zip(self._codes, range(self._codes.shape[0])))

        # Initialize empty np.array
        data = np.zeros(
            (len(self.pids), self._codes.shape[0]), dtype=np.float32
        )

        # Record diagnosis "fingerprint"
        cursor = self.conn.cursor()
        cursor.execute(
            """
            SELECT DISTINCT pid, code FROM nomesco
            WHERE time < %s AND pid = ANY(%s)
            """,
            (self.horizon, self.pids)
        )

        for pid, code in cursor:
            data[pid2idx[pid], code2idx[code]] = 1

        self._data = data

    def __getitem__(self, idx):

        return self._data[idx, :]

@dataclass
class PolygenicData(FeatureCollection):
    """Collection of results from different PGSes. The scores are z-score
    normalised and scaled to the interval [-1;1]. Missing values are set to 0.0
    """

    winsorize: float = 0.0

    alias = "pgss-1"

    def fetch_data(self):

        cursor = self.conn.cursor()
        cursor.execute(
            "select pid, name, score from pmhnet.polygenic_scores"
        )
        df = pd.DataFrame(cursor.fetchall())

        # winsorize and z-score normalize
        df["score"] = (df["score"]
            .groupby(df["name"])
            .transform(
                lambda x: stats.mstats.winsorize(x, limits=self.winsorize)
            )
            .groupby(df["name"])
            .transform(lambda x: (x - x.mean()) / x.std())
        )

        # long to wide and add missing observations
        df = (df
            .pivot(index="pid", columns="name", values="score")
            .reindex(self.pids, fill_value=0.0)
            .sort_index()
        )

        # set instance vars
        self._colnames = [f"{self.alias}_{col}" for col in df.columns]
        self._data = df.to_numpy(dtype=np.float32)

    def __getitem__(self, idx):

        return self._data[idx, :]
