library(tidyverse)
library(lubridate)
library(dbplyr)

patients <- read_tsv(snakemake@input$pids, col_types = cols(pid = "i"))
bio <- read_tsv(snakemake@input$bioc, col_types = cols(pid = "i"))

## Get creatinine and enzymes -------------------------------------------------

comps <- c(
    "TROPONIN I, CARDIAC MUSCLE",
    "TROPONIN T, CARDIAC MUSCLE",
    "CREATINE KINASE MB",
    "CREATININIUM"
)

bio <- bio %>%
    semi_join(patients, by = "pid") %>%
    filter(component_db %in% comps)  %>%
    inner_join(patients, by = "pid") %>%
    mutate(
        t_since = time_length(index_date %--% date, "days"),
        across(shown_clean, ~suppressWarnings(as.numeric(.x)))
    ) %>%
    filter(t_since >= -90, t_since <= 21, !is.na(shown_value))

# For creatinine, find a measurement close to the index_date. We first try to
# find a measurement before (or on) the index_date, but if not available, then
# we allow measurements from the future with a 21 day cutoff.

crea <- bio %>%
    filter(component_db == "CREATININIUM", unit_clean == "umol/L") %>%
    arrange((t_since - 1) > 0, abs(t_since)) %>%
    distinct(pid, .keep_all = TRUE) %>%
    transmute(pid, crea = shown_clean)

# For cardiac enzymes, it's more problematic to allow measurements from the
# "future", and therefore have to use a more strict cutoff. To start with, we
# use a threshold of 7 days. In this case, we are not using the most recent
# value, but instead just check if any of the measurements are above the upper
# reference range

enzymes <- bio %>%
    filter(component_db != "CREATININIUM", t_since < 7) %>%
    group_by(pid) %>%
    summarise(enzymes = any(as.logical(FLAG)), .groups = "drop")

# Finally, we join the two tables and write the result to disk

full_join(crea, enzymes, by = "pid") %>%
    write_tsv(snakemake@output$sup1)

## Get systolic_bp and pulse --------------------------------------------------

# Connection to dbserver
con <- DBI::dbConnect(odbc::odbc(),
    driver = "PostgreSQL", server = "server", database = "database",
    UID = Sys.getenv("uid"), PWD = Sys.getenv("pwd")
)

patients_db <- copy_to(con, patients, temporary = TRUE)
txt <- tbl(con, in_schema("epj", "txt"))

# Get all pulse measurements
pulse <- tbl(con, in_schema("rxtractor", "pulse")) %>%
    inner_join(txt, by = c("entryid" = "txt_id")) %>%
    inner_join(patients_db, by = "pid") %>%
    select(pulse, index_date, rekvdt, pid) %>%
    collect()

# For pulse, we find recorded measurements within half year prior to and one
# month after the index_date. As with enzymes and creatinine, we prioritize
# measurements prior to the event, but allow future measurements if otherwise
# missing.

pulse <- pulse %>%
    mutate(
        across(rekvdt, ymd_hm),
        t_since = time_length(index_date %--% rekvdt, "days")
    ) %>%
    filter(t_since %>% between(-182, 31)) %>%
    arrange((t_since - 1) > 0, abs(t_since)) %>%
    distinct(pid, .keep_all = TRUE) %>%
    select(pid, pulse)

# Get all blood pressure measurements
sys_bp <- tbl(con, in_schema("rxtractor", "bloodpressure")) %>%
    inner_join(txt, by = c("entryid" = "txt_id")) %>%
    inner_join(patients_db, by = "pid") %>%
    select(bp_sys=systolic, index_date, rekvdt, pid) %>%
    collect()

# For the systolic blood pressure meaurements, we use the exact same approach
# as with the pulse measurements.

sys_bp <- sys_bp %>%
    mutate(
        across(rekvdt, ymd_hm),
        t_since = time_length(index_date %--% rekvdt, "days")
    ) %>%
    filter(t_since %>% between(-182, 31)) %>%
    arrange((t_since - 1) > 0, abs(t_since)) %>%
    distinct(pid, .keep_all = TRUE) %>%
    select(pid, sys_bp = bp_sys)

# Join the two tables, and write to disk
full_join(pulse, sys_bp, by = "pid") %>%
    write_tsv(snakemake@output$sup2)
