library(tidyverse)
library(dbplyr)

# Redirect stdout to logging file
log <- file(snakemake@output[[1]], "wt")
sink(log)

# Establish connection to postgres databases
con <- DBI::dbConnect(odbc::odbc(),
    driver = "PostgreSQL", server = "dbserver", database = "database",
    uid = Sys.getenv("uid"), pwd = Sys.getenv("pwd")
)

# Get all pids in cohort
pids <- con %>%
    tbl(ident_q("pmhnet.baseline")) %>%
    select(pid) %>%  collect()

# Read files and filter to only include pids in the cohort. Due to some
# problems in calculating the PGSes, the direction of all but the "engage"
# PGSes has been flipped and needs to be adjusted accordingly
df <- tibble(path = snakemake@input$pgs) %>%
    mutate(
        name = str_extract(path, "[^/]+(?=.tsv)"),
        data = map(path, read_tsv, col_types = cols(cpr = "i"))
    ) %>%
    unnest(data) %>%
    semi_join(pids, by = c("cpr" = "pid")) %>%
    mutate(
        score = if_else(
            str_ends(name, "ENGAGE"),
            true  = combined_risk_score,
            false = -1 * combined_risk_score
        )
    ) %>%
    select(pid = cpr, name, score)

glue::glue(
    "{n1} of {n2} individuals with available PGSes",
    n1 = n_distinct(df$pid),
    n2 = n_distinct(pids$pid)
)

df %>%
    {DBI::dbWriteTable(
        conn = con,
        name = DBI::SQL("pmhnet.polygenic_scores"),
        value = as.data.frame(collect(.))
    )}

glue::glue(
    "Table `pmhnet.polygenic_scores` created\nScript finished: {ts}",
    ts = format(Sys.time(), "%a %b %d %X %Y")
)
