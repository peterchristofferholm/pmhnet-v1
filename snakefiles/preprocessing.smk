# PREPROCESSING RULES ---------------------------------------------------------

configfile: "config.yaml"

rule inclusion_criteria:
    input:
        dia = "data/raw_symlinks/t_diagadms.tsv",
        bds = "data/raw_symlinks/t_person.tsv",
        kag = "data/raw_symlinks/pats-kag.tsv"
    output: 
        pat = "data/interim/patients.tsv",
        dia = "data/interim/diagnoses.tsv", 
        kag = "data/interim/pats-kag.tsv",
        srv = "data/interim/survival.tsv"
    log: "logs/00_define-inclusion.log"  # attrition stats
    params: test_size=5000
    resources: vmem=1024*20, tmin=15
    script: "../scripts/preprocessing/00_inclusion-criteria.R" 

rule prepare_extra_clinical_features:
    input: 
        pids = "data/interim/patients.tsv",
        bioc = "data/raw_symlinks/biochem.tsv"
    output:
        sup1 = "data/interim/crea_and_enz.tsv",
        sup2 = "data/interim/pulse_and_sbp.tsv"
    resources: vmem=1024*150, tmin=60
    script: "../scripts/preprocessing/01_prepare-extra-data.R"

rule populate_database:
    input:
        # schema and definitions
        sql = "scripts/database_schema.sql",
        icd = "data/resources/icd10_definitions.tsv",
        ref = "data/raw_symlinks/refintervals_lookup.tsv",
        # interim data files
        srv = "data/interim/survival.tsv",
        kag = "data/interim/pats-kag.tsv",
        pat = "data/interim/patients.tsv",
        dia = "data/interim/diagnoses.tsv",
        sup1 = "data/interim/pulse_and_sbp.tsv",
        sup2 = "data/interim/crea_and_enz.tsv"
    output: "logs/02_populate-database.log"
    params: 
        uid = os.environ["DBUSER"], 
        pwd = os.environ["DBPASS"]
    resources: vmem=1024*5, tmin=30
    script: "../scripts/preprocessing/02_populate-database.R" 

rule reference_scores:
    input: "logs/02_populate-database.log"
    output: "logs/03_reference-scores.log"
    threads: 12
    resources: vmem=1024*5, tmin=60*10
    script: "../scripts/preprocessing/03_risk-scores.R"

rule prepare_clinical_two:
    input:
        log = "logs/02_populate-database.log",
        kag = "data/interim/pats-kag.tsv",
        pat = "data/interim/patients.tsv"
    output: "logs/04_prepare-clinical-two.log"
    resources: vmem=1024*10, tmin=15
    script: "../scripts/preprocessing/04_prepare-clinical-two.R"

rule prepare_nomesco_codes:
    input: 
        "logs/02_populate-database.log",
        "data/interim/patients.tsv"
    output: "logs/05_nomesco-codes.log"
    resources: vmem=1024*5, tmin=60
    script: "../scripts/preprocessing/05_prepare-nomesco-codes.R"

# currently included PGS
pgs, = glob_wildcards("data/raw_symlinks/pgs/{pgs}.tsv")

rule prepare_pgs:
    input: 
        pgs=expand("data/raw_symlinks/pgs/{pgs}.tsv", pgs=pgs),
        prev="logs/05_nomesco-codes.log"
    output: ".database_timestamp"  # last updated
    resources: vmem=1024*5, tmin=60
    script: "../scripts/preprocessing/06_prepare-pgs-table.R"
