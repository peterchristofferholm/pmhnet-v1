import numpy as np
import pandas as pd

import os
import shap
import optuna
import torch
import psycopg2
import torch.nn as nn
import importlib.util

torch.device("cpu")

import sys; sys.path.append("./")
from pmhnet.network import HazToSurv
from datetime import datetime
from sqlalchemy import create_engine


sid = snakemake.wildcards.sid
schema1 = snakemake.config["dbschema_optuna"]
schema2 = snakemake.config["dbschema_pmhnet"]
params = snakemake.params

## import hpo script
spec = importlib.util.spec_from_file_location(
    "hpo", f"scripts/studies/study_{sid}.py"
)
hpo = importlib.util.module_from_spec(spec)
spec.loader.exec_module(hpo)

## get hyperparams from optuna
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

## load model and data
with psycopg2.connect(dsn, options=f"-c search_path={schema2}") as conn:
    dtrain = hpo.define_data(best_trial, conn, params, training=True)
    dtest  = hpo.define_data(best_trial, conn, params, training=False)

model  = hpo.define_model(best_trial, in_shape=dtrain.n_features, params=params)
model.load_state_dict(torch.load(snakemake.input[0]))

# The model is patched to produce outputs on the survival scale. This is done
# by adding a cumulative product layer on top of the original output layer.
model.add_module("surv", HazToSurv(dtrain.breaks, params.xout))
model.eval()

## assemble data for shap
X_train = torch.from_numpy(np.vstack([d[1] for d in dtrain]))
X_train = X_train.type(torch.FloatTensor)
X_test = torch.from_numpy(np.vstack([d[1] for d in dtest]))
X_test = X_test.type(torch.FloatTensor)

# GradientSHAP on entire test set
es = shap.GradientExplainer(model, X_train)
attr = es.shap_values(X_test)

df = pd.concat({
    "shap" : pd.DataFrame(attr, index=dtest.pids, columns=dtest.colnames),
    "vals" : pd.DataFrame(
        X_test.detach().numpy(), index=dtest.pids, columns=dtest.colnames
    )
})

# DataFrame from wide to long format
df = (df
    .rename_axis(index=["type", "pid"], columns=["feature"])
    .stack(level="feature")
    .unstack(level="type")
    .reset_index()
)

df["timepoint"] = params.xout

df.to_sql(
    con=create_engine(dsn), name=f"shap_{sid}", schema=schema2,
    method="multi", if_exists="append", index=False
)

with open(snakemake.output[0], "wt") as f:
    f.write(f"timestamp: {datetime.now()}\n")
