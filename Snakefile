stdies = glob_wildcards("scripts/studies/study_{sid}.py").sid

workflowsingularity_args = "-B $HOME/.odbc.ini:/etc/odbc.ini"

## TARGET RULES ###############################################################

rule get_studies:
    input: expand("logs/studies/study_{sid}.log", sid=studies)

rule get_evals:
    input: expand("results/study_{sid}/evaluation", sid=studies)

rule get_shap:
    input: 
        files = [
            f"logs/shap/s{sid}-{tp}.log" 
            for sid in studies for tp in "6m 1y 3y 5y".split()
        ]

## SETUP AND PREPROCESSING ####################################################

container: "/users/singularity/pmhnet.latest.sif"  # singularity image
configfile: "config.yaml"

## FIND BEST HYPERPARAMETERS ##################################################

def get_vmem(wildcards, attempt, threads):
    return int(1024 * threads * (1 + (attempt - 1)/4))  # add 25% 

rule hpo_search:
    output: 
        logfile  = "logs/studies/study_{sid}.log",
        trialobj = "results/study_{sid}/best-trial.pkl"
    resources: mem_mb=get_vmem, tmin=60*24*7, queue="batch"
    benchmark: "logs/benchmarks/study_{sid}"
    threads: 40
    params: 
        n_trials    = 2500,  # number of trials in hpo
        max_epochs  = 1000,  # max number of epochs per trial
        n_intervals = 30,    # number of prediction bins
        max_time    = 365*5, # max prediction time
    script: "scripts/studies/study_{wildcards.sid}.py"

rule hpo_plots:
    input: "logs/studies/study_{sid}.log"
    output: 
        p1 = "plots/hpo/s{sid}-history.pdf",
        p2 = "plots/hpo/s{sid}-sliceplot.pdf"
    resources: mem_mb=1024*5, tmin=30
    script: "scripts/hpo_plots.R"

## FIT MODELS AND EVALUATE ####################################################

rule model_training:
    input: "logs/studies/study_{sid}.log"
    params: **rules.hpo_search.params
    output: 
        history  = "results/study_{sid}/training_history.txt",
        m_dict   = "results/study_{sid}/model_dict.pt",
        p_train  = "results/study_{sid}/train_predictions.csv",
        p_test   = "results/study_{sid}/test_predictions.csv",
        features = "results/study_{sid}/features.txt"
    resources: mem_mb=1024*5, tmin=60*6
    script: "scripts/model_training.py"

def tp2days(wc):
    return {"6m" : 182.5, "1y" : 365, "3y" : 1095, "5y" : 1825}[wc.tp]

rule model_explainability:
    input: "results/study_{sid}/model_dict.pt"
    output: "logs/shap/s{sid}-{tp}.log"
    params: **rules.hpo_search.params, xout = tp2days
    resources: mem_mb=1024*30, tmin=60*2
    script: "scripts/explainability.py"
