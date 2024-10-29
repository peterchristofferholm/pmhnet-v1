import sys; sys.path.append("./")
import numpy as np
import pandas as pd

import os
import optuna
import psycopg2
import pmhnet
import importlib.util

from torch.utils.data import DataLoader
import torch; torch.device("cpu")

# get variables from snakemake
sid = snakemake.wildcards["sid"]
schema1 = snakemake.config["dbschema_optuna"]
schema2 = snakemake.config["dbschema_pmhnet"]
params = snakemake.params

# import functions from HPO study
spec = importlib.util.spec_from_file_location(
    "hpo_study", f"scripts/studies/study_{sid}.py"
)
hpo = importlib.util.module_from_spec(spec)
spec.loader.exec_module(hpo)


# dsn for database connection
dsn = "postgresql://{dbuser}:{dbpass}@{dbhost}/{dbname}".format(
    dbuser=snakemake.config["dbuser"],
    dbhost=snakemake.config["dbhost"],
    dbname=snakemake.config["dbname"],
    dbpass=os.environ["DBPASS"]
)

# setup url for database connection
storage = optuna.storages.RDBStorage(
    url=dsn, engine_kwargs={
        "pool_size": 0,  # 0 is no limit
        "connect_args" : {"options" : f"-c search_path={schema1}"}
    }
)

if snakemake.wildcards["sid"] == "007":
    sid = "006"

# get best trial from HPO study
study = optuna.load_study(study_name=sid, storage=storage)
trial = study.best_trial

if snakemake.wildcards["sid"] == "007":
    trial.params["pgs_limit"] = 0.05

with psycopg2.connect(dsn, options=f"-c search_path={schema2}") as conn:

    dtrain = hpo.define_data(trial, conn, params, training=True)
    dtest  = hpo.define_data(trial, conn, params, training=False)

# save list of features
with open(snakemake.output.features, "wt") as f:
    for feature in dtrain.colnames:
        f.write(feature + "\n")

model  = hpo.define_model(trial, in_shape=dtrain.n_features, params=params)
batchsize = trial.params["batchsize"]

# training of model with full training data
training = pmhnet.training.Trainer(
    model, loss_fn=pmhnet.network.NegLogLik(),
    optimizer_fn=hpo.define_optimizer(trial),
    train_dl=DataLoader(dtrain, batch_size=batchsize),
    test_dl=DataLoader(dtest, batch_size=batchsize),
    max_epochs=trial.last_step - trial.params["patience"]
)

# run training and collect history
hist = [(epoch.epoch, epoch.train_loss, epoch.val_loss) for epoch in training]
np.savetxt(snakemake.output.history, np.stack(hist))

# save model weigths and biases
model.eval()
torch.save(model.state_dict(), snakemake.output["m_dict"])

# compute model predictions
with torch.no_grad():

    breaks = dtrain.breaks[1:]

    # on testing data
    x = torch.from_numpy(np.vstack([d[1] for d in dtest]))
    y_hat = model(x).detach().numpy()
    preds_test = pd.DataFrame(
        np.cumprod(y_hat[:, 0], axis=1), columns=breaks, index=dtest.pids
    )

    # on training data
    x = torch.from_numpy(np.vstack([d[1] for d in dtrain]))
    y_hat = model(x).detach().numpy()
    preds_train = pd.DataFrame(
        np.cumprod(y_hat[:, 0], axis=1), columns=breaks, index=dtrain.pids
    )

# write predictions to disk
preds_test.to_csv(snakemake.output["p_test"], index_label="pid")
preds_train.to_csv(snakemake.output["p_train"], index_label="pid")
