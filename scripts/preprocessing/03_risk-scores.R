library(tidyverse)
library(doParallel)
library(dbplyr)
library(DBI)

# GRACE2.0
source("R/cardrisk.R")

# Connect to postgres database
con <- dbConnect(odbc::odbc(),
    driver = "PostgreSQL", server = "dbserver", database = "database",
    uid = Sys.getenv("uid"), pwd = Sys.getenv("pwd")
)

# Redirect stdout to logging file
log <- file(snakemake@output[[1]], "wt")
sink(log)

## ----------------------------------------------------------------------------

surv <- tbl(con, ident_q("pmhnet.surv_master"))
base <- tbl(con, ident_q("pmhnet.baseline"))

# Move tables to memory
df <- inner_join(surv, base, by = "pid") %>% collect()

# Rename sys_bp, and set killip to 1 if missing
df <- df %>%
    rename(sbp = sys_bp) %>%
    mutate(killip = coalesce(killip, 1))

# Subset of data with no missing grace arguments
df_nm <- df %>%
    filter(across(formalArgs(grace_v2), ~!is.na(.x)))

# For imputation, convert categorical columns to factors. Also, we have to
# convert the tibble to a classical data.frame. The imputation is performed
# with missForest

df_impute <- df %>%
    select(-pid, -train, -time, -event) %>%
    mutate(
        across(
            c(sex, arrest, enzymes, stemi, smoking, vessels, dominance,
              nyha, ccs, killip), as_factor
        )
    ) %>%
    as.data.frame()

registerDoParallel(cores = snakemake@threads)
mf <- missForest::missForest(
    df_impute, maxiter = 15, parallelize = "forests", ntree = 75
)

# Compute grace risk scores for all individuals
grace_risks <- as_tibble(mf$ximp) %>%
    add_column(pid = df$pid, .before = 1) %>%
    mutate(
        crea = crea / 88.42,  # convert to mg/dL
        across(
            c(arrest, enzymes, stemi, killip), ~as.numeric(as.character(.x))
        )
    ) %>%
    rowwise() %>%
    transmute(
        pid, imputed = !pid %in% df_nm$pid,
        grace = list(
            grace_v2(age, pulse, sbp, crea, killip, arrest, enzymes, stemi)
        )
    ) %>%
    unnest(grace) %>%
    transmute(
        pid, imputed,
        time = case_when(
            model == "6m" ~ 183,
            model == "1y" ~ 365,
            model == "3y" ~ 1095
        ),
        type = "GRACE2.0",
        value = prob
    )

grace_risks <- grace_risks %>% select(pid, time, type, imputed, value)

dbWriteTable(con,
    SQL("pmhnet.risk_scores"), grace_risks, append = TRUE, row.names = FALSE
)

print(glue::glue(
    "GRACE2.0 scores computed for {x} of {y} patients",
    x = n_distinct(grace_risks$pid),
    y = n_distinct(df$pid)
))

print(glue::glue(
    "GRACE2.0 variables imputed for {x} of {y} patients",
    x = n_distinct(df$pid) - n_distinct(df_nm$pid),
    y = n_distinct(df$pid)
))

print("Percentage missing in dataset:")
df %>%
    select_at(formalArgs(grace_v2)) %>%
    summarise(across(everything(), ~sum(is.na(.x))/n()))
