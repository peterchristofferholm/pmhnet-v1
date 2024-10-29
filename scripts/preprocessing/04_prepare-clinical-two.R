library(tidyverse)
library(dbplyr)

pats <- read_tsv(snakemake@input$kag, col_types = cols())
meta <- read_tsv(snakemake@input$pat, col_types = cols())

# Prepare variables not already prepared in 02_populate-database.R
clnc2 <- inner_join(meta, pats, by = c("entryid", "pid")) %>%
  transmute(
    pid,
    familiary_ihd = case_when(
      familiary_ihd == "Ja" ~ "yes",
      familiary_ihd == "Nej" ~ "no",
      TRUE ~ "missing"
    ),
    icd_or_pm = arrhythmia_device %in% c("BIV-ICD", "BIV-PM", "PM", "ICD") %>%
      coalesce(FALSE),
    lvef = if_else(lvef > 60, 60, lvef) %>%
      cut_interval(length = 10) %>%
      fct_explicit_na("missing"),
    ischemia_test = case_when(
      ischemia_test == "Ikke konklusiv" ~ "inconclusive",
      ischemia_test == "Ikke udført" ~ "missing",
      ischemia_test == "Negativ" ~ "negative",
      ischemia_test == "Positiv" ~ "positive",
      is.na(ischemia_test) ~ "missing"
    ),
    abnormal_ekg = case_when(
      abnormal_ekg == "Ja"  ~ "yes",
      abnormal_ekg == "Nej" ~ "no",
      TRUE ~ "missing"
    ),
    abnormal_qrs_st = case_when(
      abnormal_qrs == "Ja"  ~ "yes",
      abnormal_qrs == "Nej" ~ "no",
      TRUE ~ "missing"
    )
  )

# Establish connection to database
con <- DBI::dbConnect(odbc::odbc(),
    driver = "PostgreSQL",
    server = snakemake@config$dbhost, database = snakemake@config$dbname,
    uid = snakemake@params$dbuser, pwd = Sys.getenv("DBPASS")
)

# Combine features and move to new table
clnc2 <- tbl(con, in_schema("pmhnet", "baseline")) %>%
  select(
    pid, sex, height, weight, dia_bp, smoking, vessels, dominance,
    nyha, ccs
  ) %>%
  inner_join(clnc2, by = "pid", copy = TRUE) %>%
  compute(
    name = in_schema("pmhnet", "clinicaltwo"),  # destination
    temporary = FALSE,
    unique_indexes = list("pid")
  )

# Save timestamp to log-file
write_file(as.character(Sys.time()), snakemake@output[[1]])
