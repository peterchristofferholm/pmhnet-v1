rule base_overview_plots:
    input: ".database_timestamp"
    output: 
        p01 = "plots/misc/s1_km-curves.pdf",
        txt = "stats/survival-stats.txt"
    resources: vmem=1024, tmin=10
    script: "../scripts/baseline_plots.R"

rule experiment_removing_data:
    input: wabs  = "results/study_{sid}/model_dict.pt"
    output: "results/study_{sid}/experiments/sparse-predictions.csv"
    params: **rules.hpo_search.params
    resources: mem_mb=1024*5, tmin=60
    script: "../scripts/evaluation/experiment_removing_features.py"

## EVALUATION OF 006 #########################################################

rule m006_performance_plots:
    input: 
        tra = "results/study_006/train_predictions.csv",
        val = "results/study_006/test_predictions.csv"
    output: 
        "results/study_006/210712/01_yhat-histogram.pdf", 
        "results/study_006/210712/02_yhat-vs-km.pdf",
        "results/study_006/210712/03_combined.pdf", 
        "results/study_006/210712/04_calibration.pdf"
    resources: vmem=1024*5, tmin=60
    script: "../scripts/evaluation/210712_performance-measures.R"
