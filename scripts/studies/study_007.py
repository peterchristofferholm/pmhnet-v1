import sys; sys.path.append("./")
import torch
import torch.nn as nn
import torch.nn.functional as F
import torch.optim as optim

torch.device("cpu")

from pmhnet.datasets import DataAssembler

from pmhnet.network import SurvProbability, NegLogLik
from pmhnet.training import CrossValidation
from pmhnet.utils import kfold_loaders, study_exists

from functools import partial
from concurrent.futures import ProcessPoolExecutor
from datetime import datetime

import psycopg2
import optuna
import git
import os
import logging


def define_data(trial, conn, params, training=True):

    data = DataAssembler(conn, training)
    data.add_labels(params.n_intervals, params.max_time)

    data.add_features("clnc-1")
    data.add_features("clnc-2")

    data.add_features("diag-simple",
        horizon=trial.suggest_float("diag_window", 0.5, 20, step=0.5)*365,
        code_level=trial.suggest_categorical(
            "diag_level", ["lvl4", "lvl3", "block", "chap"]
        )
    )

    data.add_features("bioc-1",
        horizon=trial.suggest_int("bioc_window", 0.5, 5, step=0.25)*365,
    )
    data.add_features("proc-1",
        horizon=trial.suggest_float("proc_window", 0.5, 20, step=0.5)*365
    )

    data.add_features("pgss-1",
        winsorize=trial.suggest_float("pgs_limit", 0.0, 0.49, step=0.01)
    )

    return data


def define_model(trial, in_shape, params):

    layers = []
    n_layers = trial.suggest_int("n_layers", 0, 3)

    for i in range(n_layers):

        out_shape = trial.suggest_int(f"n_units_l{i}", 10, 200)
        layers.extend((
            nn.Linear(in_shape, out_shape), nn.ReLU()
        ))
        droprate = trial.suggest_float(f"droprate_l{i}", 0.05, 0.80)
        layers.append(nn.Dropout(droprate))

        in_shape = out_shape

    layers.append(
        SurvProbability(in_shape, params.n_intervals)
    )

    return nn.Sequential(*layers)


def define_optimizer(trial):

    lr = trial.suggest_float("lr", 1e-6, 1, log=True)
    momentum = trial.suggest_float("momentum", 0.0, 1.0)

    return partial(optim.SGD, lr=lr, momentum=momentum)


def objective(trial, conn, params):

    data = define_data(trial, conn, params)
    model = define_model(trial, in_shape=data.n_features, params=params)
    optimizer_fn = define_optimizer(trial)

    # Initialize cv-training class
    kf_splits = kfold_loaders(
        data, k=5, batch_size=params.batchsize, shuffle=True
    )
    patience = trial.suggest_int("patience", 20, 500)

    cv_training = CrossValidation(
        model, kf_splits, NegLogLik(), optimizer_fn,
        patience=patience, max_epochs=params.max_epochs
    )

    for epoch in cv_training:

        # Get loss and epoch number
        loss = epoch.loss
        epoch = epoch.epoch

        # Report intermediate value
        trial.report(loss, epoch)

        # Handle pruning based on intermediate value
        if trial.should_prune():
            raise optuna.TrialPruned()

    return loss


def optimize(args):

    storage, n_trials = args

    study = optuna.load_study(
        storage=storage, study_name=STUDY_NAME,
        sampler=optuna.samplers.TPESampler(
            n_startup_trials=400
        ),
        pruner=optuna.pruners.HyperbandPruner(
            max_resource=params.max_epochs, min_resource=20, reduction_factor=3
        )
    )

    connection = psycopg2.connect(
        DSN, options = f"-c search_path={PMHNET_SCHEMA}"
    )

    with connection as conn:

        _objective = partial(objective, conn=conn, params=params)
        study.optimize(_objective, n_trials=n_trials, n_jobs=1)


if __name__ == "__main__":

    params = snakemake.params

    # get parameters from snakemake config
    OPTUNA_SCHEMA = snakemake.config["dbschema_optuna"]
    PMHNET_SCHEMA = snakemake.config["dbschema_pmhnet"]
    STUDY_NAME = snakemake.wildcards["sid"]

    # setup url for database connection
    DSN = "postgresql://{dbuser}:{dbpass}@{dbhost}/{dbname}".format(
        dbuser=snakemake.config["dbuser"],
        dbhost=snakemake.config["dbhost"],
        dbname=snakemake.config["dbname"],
        dbpass=os.environ["DBPASS"]
    )

    # setup logging
    logging.basicConfig(
        filename=snakemake.output[0], filemode="w", level=logging.INFO,
        format="%(levelname)s | %(message)s"
    )
    logging.info(f"Script started: {datetime.utcnow().isoformat()}")

    # setup optuna storage
    storage = optuna.storages.RDBStorage(
        url=DSN, engine_kwargs={
            "pool_size": 0,  # 0 is no limit
            "connect_args" : {"options" : f"-c search_path={OPTUNA_SCHEMA}"}
        }
    )

    # remove old versions of study if exists and not protected
    if study_exists(STUDY_NAME, storage):
        if STUDY_NAME in snakemake.config["protected"]:
            raise RuntimeError(f"Study {STUDY_NAME} is protected.")
        else:
            optuna.delete_study(STUDY_NAME, storage)

    study = optuna.create_study(
        storage=storage, study_name=STUDY_NAME, direction="minimize"
    )

    study.set_user_attr("contributors", ["pechris"])

    n_jobs = snakemake.threads
    n_trials = snakemake.params.n_trials // n_jobs

    with ProcessPoolExecutor(max_workers=n_jobs) as executor:
        executor.map(optimize, [(storage, n_trials)] * n_jobs)

    is_pruned = lambda t: t.state == optuna.trial.TrialState.PRUNED
    is_complete = lambda t: t.state == optuna.trial.TrialState.COMPLETE
    pruned_trials = [t for t in study.trials if is_pruned(t)]
    complete_trials = [t for t in study.trials if is_complete(t)]

    logging.info("Study statistics: ")
    logging.info(f"  Number of finished trials: {len(study.trials)}")
    logging.info(f"  Number of pruned trials: {len(pruned_trials)}")
    logging.info(f"  Number of complete trials: {len(complete_trials)}")

    logging.info("Best trial:")
    trial = study.best_trial

    logging.info(f"  Value: {trial.value}")
    logging.info("  Params:")
    for key, value in trial.params.items():
        logging.info(f"    {key}: {value}")

    repo = git.Repo(search_parent_directories=True)
    logging.info("Study finished.")
    logging.info(f"last_commit: {repo.head.commit.hexsha}")
    logging.info(f"timestamp: {datetime.now()}")




