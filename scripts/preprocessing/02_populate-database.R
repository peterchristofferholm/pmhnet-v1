library(tidyverse)
library(lubridate)
library(DBI)  # communicate with DB
library(rlang)

## Misc. functions used -------------------------------------------------------
swap_if <- function(cond, x, y) {
    names <- c(as_name(ensym(x)), as_name(ensym(y)))
    out_x <- if_else(cond, y, x)
    out_y <- if_else(cond, x, y)
    set_names(tibble(out_x, out_y), names)
}

extract_roman <- function(x) {
    x %>% str_extract("[IV]+") %>% as.roman() %>% as.integer()
}

trim_icd10 <- function(x) {
    if_else(str_length(x) == 4L, str_c(x, "9"), str_sub(x, end = 5L))
}

## Setup and initialization ---------------------------------------------------

# connect to DB server and send creation query from file
con <- dbConnect(odbc::odbc(),
    driver = "PostgreSQL", server = "server", database = "database",
    uid = snakemake@params$uid, pwd = snakemake@params$pwd
)

dbExecute(con, statement = read_file(snakemake@input$sql))

# Interim data files
surv <- read_tsv(snakemake@input$srv, col_types = cols())
pats <- read_tsv(snakemake@input$kag, col_types = cols())
meta <- read_tsv(snakemake@input$pat, col_types = cols())
diag <- read_tsv(snakemake@input$dia, col_types = cols())

## populate pmhnet.surv_master ------------------------------------------------

surv_db <- surv %>%
    select(pid, train, time, event)

dbWriteTable(
    con, SQL("pmhnet.surv_master"), surv_db, append = TRUE, row.names = FALSE
)

## prepare pmhnet.baseline ----------------------------------------------------

pulse_and_sbp <- read_tsv(snakemake@input$sup1, col_types = cols())
crea_and_enz <- read_tsv(snakemake@input$sup2, col_types = cols())

# part1 of baseline
baseline_1 <- meta %>%
    right_join(pats, by = "pid") %>%
    rename(height = height_cm, weight = weight_kg, crea = creatinine) %>%
    mutate(
        sex = (sex == "M"),
        weight = if_else(weight > 500, weight / 10, weight),
        swap_if(height < 150 & weight > 150, weight, height),
        height = if_else(height < 130 | height > 250, NA_real_, height),
        sys_bp = if_else(sys_bp > 300, NA_real_, sys_bp),
        dia_bp = if_else(dia_bp > 300, NA_real_, dia_bp),
        swap_if(dia_bp > sys_bp, sys_bp, dia_bp),
        enzymes = (enzymes == "Positiv"),
        arrest = coalesce(arrest == "Ja", FALSE),
        stemi = coalesce(stemi, FALSE)
    ) %>%
    left_join(crea_and_enz, by = "pid") %>%
    left_join(pulse_and_sbp, by = "pid") %>%
    mutate(
        crea = coalesce(crea.x, crea.y),
        enzymes = coalesce(enzymes.x, enzymes.y),
        pulse = coalesce(pulse.x, pulse.y) %>% as.integer(),
        sys_bp = coalesce(sys_bp.x, sys_bp.y) %>% as.integer(),
        crea = if_else(crea > 400, NA_real_, crea),
        pulse = if_else(pulse > 300, NA_integer_, pulse),
        pulse = if_else(pulse < 30, NA_integer_, pulse)
    ) %>%
    arrange(pid, proc_date) %>%
    select(
        pid, sex, age, starts_with("entryid"), height, weight, pulse, sys_bp,
        dia_bp, crea, arrest, enzymes, stemi
    ) %>%
    group_by(pid) %>%
    fill(height:dia_bp, .direction = "downup") %>%
    ungroup() %>%
    filter(entryid.x == entryid.y) %>%
    distinct() %>%
    select(-starts_with("entryid"))

baseline_2 <- meta %>%
    inner_join(pats, by = c("pid", "entryid")) %>%
    transmute(
        pid,
        smoking = case_when(
            str_starts(smoking_status, "Akt") ~ "yes",
            str_starts(smoking_status, "Ald") ~ "no",
            str_starts(smoking_status, "Eks") ~ "ex",
            TRUE ~ NA_character_
        ),
        vessels = case_when(
            str_starts(native_kar, "1 gebet") ~ "1VD",
            str_starts(native_kar, "2 gebet") ~ "2VD",
            str_starts(native_kar, "3 gebet") ~ "3VD",
            TRUE ~ "DIF"
        ),
        dominance = case_when(
            str_starts(dominance, "Balanceret") ~ "B",
            str_starts(dominance, "Højre")      ~ "R",
            str_starts(dominance, "Venstre")    ~ "L",
            TRUE ~ NA_character_
        ),
        nyha = extract_roman(nyha_class),
        ccs = case_when(
            str_starts(ccs_class, "CCS") ~ extract_roman(ccs_class),
            str_starts(ccs_class, "Ing") ~ 0L,
            TRUE ~ NA_integer_
        ),
        killip = extract_roman(killip_class)
    )

baseline_db <- inner_join(baseline_1, baseline_2, by = "pid")

dbWriteTable(
    con, SQL("pmhnet.baseline"), baseline_db, append = TRUE, row.names = FALSE
)

## populate pmhnet.diagnoses --------------------------------------------------

# lvl3, block & chapter parents for each lvl4 code
icd10_parents <- snakemake@input$icd %>%
    read_tsv(col_types = cols()) %>%
    mutate(icd10_code = str_c("D", icd10)) %>%
    filter(str_length(icd10_code) == 5L) %>%
    transmute(
        code_lvl4 = icd10_code,
        code_lvl3 = str_sub(icd10_code, end = 4L),
        code_block = block,
        code_chap = str_c("Ch", chapter)
    )

# wrangle data into correct shape
diag_db <- diag %>%
    inner_join(select(meta, pid, index_date), by = "pid") %>%
    filter(
        str_starts(icd10_type, "A|B"),  # primary or secondary codes
        disc_dttm < index_date          # only codes before index
    ) %>%
    mutate(
        code_lvl4 = trim_icd10(icd10_code),
        time = time_length(index_date %--% disc_dttm, "days") %>% round(0)
    ) %>%
    select(pid, time, code_lvl4) %>%
    arrange(pid, -time) %>%
    inner_join(icd10_parents, by = "code_lvl4") %>%
    pivot_longer(
        starts_with("code"),  names_to = "level",
        names_prefix = "code_", values_to = "code"
    )

# remove codes not present in more than 0.5% of the population
rare_codes <- diag_db %>%
    distinct(pid, code) %>%
    count(code, sort = TRUE) %>%
    filter(n / nrow(surv) < 0.005)

diag_db <- diag_db %>%
    anti_join(rare_codes, by = "code")

# Populate diagnoses table in database
dbWriteTable(
    con, SQL("pmhnet.diagnoses"), diag_db, append = TRUE, row.names = FALSE
)

## populate pmhnet.biochem ----------------------------------------------------

con_bth <- dbConnect(odbc::odbc(),
    driver = "PostgreSQL", server = "dbserver", database = "database",
    uid = snakemake@params$uid, pwd = snakemake@params$pwd
)

# we need data from bth.biochem.tests
biochem <- tbl(con_bth, dbplyr::ident_q("biochem.tests"))

# and the lookup table for reference intervals
refint <- read_tsv(snakemake@input$ref, col_types = "cddDDcccDD")
refint <- copy_to(con_bth, refint)

# move data to temp. table in bth database for efficient joins
meta_db <- copy_to(con_bth, select(meta, pid, index_date, age, sex))

# get all pre-index tests and add missing reference intervals
bio_step1 <- biochem %>%
    inner_join(meta_db, by = "pid") %>%
    mutate(
        ref_lower = if_else(ref_type == "upper_limit", NaN, ref_lower),
        ref_upper = if_else(ref_type == "lower_limit", NaN, ref_upper)
    )  %>%
    filter(
        index_date >= sql("drawn_date - interval '1 day'"),
        !is.na(result_num)
    ) %>%
    left_join(refint, by = c("quantity_id" = "npu_code", "sex")) %>%
    mutate(
        age = sql(
            "(((age*365.25)::int -(drawn_date - index_date)) || 'd')::interval"
        ),
        age_1 = sql("('P' || age_1)::interval"),  # iso8601 prefix
        age_2 = sql("('P' || age_2)::interval")
    ) %>%
    filter(
        between(drawn_date, gldfrad, gldtild) %>% coalesce(TRUE),
        between(drawn_date, gldfrad_refint, gldtild_refint) %>% coalesce(TRUE),
        between(age, age_1, age_2) %>% coalesce(TRUE)
    )

# get distinct observations + ref_intervals
bio_step2 <- bio_step1 %>%
    transmute(
        pid, db_source,
        qid = quantity_id,
        drawn_date, index_date,
        result = result_num,
        ref_1 = coalesce(ref_lower, refint_1),
        ref_2 = coalesce(ref_upper, refint_2),
        gldfrad_ref = gldfrad_refint,
        gldtild_ref = gldtild_refint
    ) %>%
    distinct() %>%
    compute()

# params for final steps
# TODO: consider adding these as hyperparameters?
cutoff <- nrow(meta) * 0.05

# wrangle columns to match sql-schema specs
bio_step3 <- bio_step2 %>%
    filter(!is.na(ref_1) | !is.na(ref_2)) %>%
    # only keep most recent test
    group_by(pid, qid) %>%
    filter(row_number(desc(drawn_date)) == 1L) %>%
    # remove rare tests
    group_by(qid) %>%
    filter(n() > cutoff) %>%
    # flag result according to refint
    ungroup() %>%
    mutate(
        time = (index_date - drawn_date),
        onesided = sql("ref_1 is null or ref_2 is null"),
        t1 = coalesce(result < ref_1, FALSE),
        t2 = coalesce(result > ref_2, FALSE)
    ) %>%
    # final columns
    transmute(
        pid, qid, time,
        value = if_else(t1 & onesided, -1L, as.integer(t1|t2))
    ) %>%
    collect()

# Populate biochemical table in database
dbWriteTable(
    con, SQL("pmhnet.biochem"), bio_step3, append = TRUE, row.names = FALSE
)

# close connection to bth database
dbDisconnect(con_bth)

## Finalizing -----------------------------------------------------------------
dbDisconnect(con)

# touch database timestamp
write_file(as.character(Sys.time()), snakemake@output[[1]])
