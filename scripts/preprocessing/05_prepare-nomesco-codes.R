library(tidyverse)
library(dbplyr)

# Redirect stdout to logging file
log <- file(snakemake@output[[1]], "wt")
sink(log)

# Establish connection to postgres databases
con_bth <- DBI::dbConnect(odbc::odbc(),
    driver = "PostgreSQL", server = "dbserver", database = "bth",
    uid = Sys.getenv("uid"), pwd = Sys.getenv("pwd")
)

con_pmh <- DBI::dbConnect(odbc::odbc(),
    driver = "PostgreSQL", server = "dbserver", database = "database",
    uid = Sys.getenv("uid"), pwd = Sys.getenv("pwd")
)

## Pull nomesco codes from LPR tables in bth database -------------------------

adm <- tbl(con_bth, ident_q("lpr_210312.adm"))

opr <- tbl(con_bth, ident_q("lpr_210312.sksopr")) %>%
    inner_join(adm, by = c("visit_id", "hospital_code")) %>%
    select(
        pid = person_id, code = operation_code, date = operation_date
    ) %>%
    filter(sql("code ~ '^K'"))

ube <- tbl(con_bth, ident_q("lpr_210312.sksube")) %>%
    inner_join(adm, by = c("visit_id", "hospital_code")) %>%
    select(
        pid = person_id, code = intervention_code, date = intervention_date
    ) %>%
    filter(sql("code ~ '^U'"))

## Remove uncommon codes ------------------------------------------------------

patients <- snakemake@input[[2]] %>%
    read_tsv(col_types = cols()) %>%
    select(pid, index_date)

codes <- union_all(opr, ube) %>%
    inner_join(patients, by = "pid", copy = TRUE) %>%
    mutate(
        pid = as.integer(pid),
        time = index_date - date,
        code = sql("substring(code, 1, 6)")
    ) %>%
    filter(time > 0) %>%
    group_by(pid, code) %>%
    filter(row_number(time) == 1L) %>%
    select(pid, code, time) %>%
    ungroup() %>%
    compute()

included_codes <- codes %>%
    count(code, sort = TRUE) %>%
    filter(n > local(nrow(patients) * 0.01))

glue::glue(
    "number of codes included: {n}",
    n = nrow(collect(included_codes))
)

## Move results to pmhnet.nomesco ---------------------------------------------

codes %>%
    semi_join(included_codes, by = "code") %>%
    {DBI::dbWriteTable(
        conn = con_pmh,
        name = DBI::SQL("pmhnet.nomesco"),
        value = as.data.frame(collect(.))
    )}

glue::glue(
    "Table `pmhnet.nomesco` created\nScript finished: {ts}",
    ts = format(Sys.time(), "%a %b %d %X %Y")
)
