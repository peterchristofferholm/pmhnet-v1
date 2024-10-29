import sys
import os
import importlib.util

from scipy.interpolate import interp1d
import numpy as np
import pandas as pd
import optuna
import psycopg2
import torch

###############################################################################
sys.path.append("./")

sid = snakemake.wildcards.sid
schema1 = snakemake.config["dbschema_optuna"]
schema2 = snakemake.config["dbschema_pmhnet"]
params = snakemake.params

spec = importlib.util.spec_from_file_location(
    "hpo", f"scripts/studies/study_{sid}.py"
)
hpo = importlib.util.module_from_spec(spec)
spec.loader.exec_module(hpo)

dsn = "postgresql://{dbuser}:{dbpass}@{dbhost}/{dbname}".format(
    dbuser=snakemake.config["dbuser"],
    dbhost=snakemake.config["dbhost"],
    dbname=snakemake.config["dbname"],
    dbpass=os.environ["DBPASS"]
)
storage = optuna.storages.RDBStorage(
    url=dsn, engine_kwargs={
        "pool_size": 0,  # 0 is no limit
        "connect_args" : {"options" : f"-c search_path={schema1}"}
    }
)
study = optuna.load_study(study_name=sid, storage=storage)
best_trial = study.best_trial

with psycopg2.connect(dsn, options=f"-c search_path={schema2}") as conn:
    dtrain = hpo.define_data(best_trial, conn, params, training=True)
    dtest  = hpo.define_data(best_trial, conn, params, training=False)

X_train = torch.from_numpy(np.vstack([d[1] for d in dtrain]))
X_train = X_train.type(torch.FloatTensor)
X_test  = torch.from_numpy(np.vstack([d[1] for d in dtest]))
X_test  = X_test.type(torch.FloatTensor)

model  = hpo.define_model(best_trial, in_shape=X_train.shape[1], params=params)
model.load_state_dict(torch.load(snakemake.input["wabs"]))

###############################################################################

colnames = pd.Series(dtest.colnames)
fnames   = colnames.str.extract("^([^_]+_[^_]+)")[0].unique()

def get_predictions(fname):
    X_hat = torch.clone(X_test)

    if fname:
        m = colnames.str.match(fname)
        idx = m[m].index.values
        # replace columns with medians
        X_hat[:, idx] = X_train[:, idx].quantile(q=0.5, axis=0)

    with torch.no_grad():
        pred = model(X_hat).detach().numpy()

    time = np.array((7*4*6, 364, 365*2, 365*3, 365*4, 365*5), dtype=np.int)
    func = interp1d(
        dtrain.breaks[1:],
        np.cumprod(pred[:, 0], axis=1),
        kind="linear",
        axis=1
    )
    out = pd.DataFrame(
        func(time),
        columns=["6m", "1y", "2y", "3y", "4y", "5y"],
        index=dtest.pids
    )
    out["fname"] = fname
    return out

get_predictions(None).to_csv(snakemake.output[0], header=True)
for fname in fnames:
    get_predictions(fname).to_csv(snakemake.output[0], mode="a", header=False)
