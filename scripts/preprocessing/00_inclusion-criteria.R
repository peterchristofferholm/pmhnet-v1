library(tidyverse)
library(lubridate)
library(magrittr)
library(glue)

# Set seed and redirect stdout to logging file
set.seed(42)
log <- file(snakemake@log[[1]], "wt")
sink(log)

# Get input paths from snakefile
dia_path <- snakemake@input[["dia"]]  # t_diagadms
bds_path <- snakemake@input[["bds"]]  # t_person
kag_path <- snakemake@input[["kag"]]  # pats-kag

# Read in the different files
dia <- read_tsv(dia_path, col_types = "ccDDiiicccc--", na = "NULL")
bds <- read_tsv(bds_path, col_types = "ccDcDcccD", na = "NULL")
kag <- read_tsv(kag_path, col_types = cols(.default = "c"))

# Collapse to one row per entryid
kag_summary <- kag %>%
    mutate_at("stenosegrad_eyeball", parse_integer) %>%
    group_by(entryid) %>%
    summarise(
        amount = str_c(stenosegrad_eyeball, collapse = ","),
        sites = str_c(proceduresite, collapse = ",")
    )

# Dates in the PATS data is in excel format
parse_exceldate <- function(x) {
    as_date(as.integer(x), origin = "1899-12-30")
}

# Select columns of interest and add more informative names
kag <- kag %>%
    transmute(
        pid = bth_pid, entryid, proc_date = parse_exceldate(proceduredato),
        prev_ami = tidligereami_frasetevt_aktuellestemiel_nstemi,
        prev_cabg = tidligerecabg,
        prev_pci = tidligerepci,
        native_kar = koronarpatologi_nativekar_udfyldesforallept_,
        dominance = dominans,
        smoking_status = rygning,
        smoking_amount = antalpakke_r,
        lung_function = lungefunktion,
        weight_kg,
        height_cm,
        nyha_class = nyhaklasse,
        ccs_class = stabilangina_ccsklasse,
        euro_logistic = logisticeuroscore,
        euro_additive = additiveeuroscore,
        pulse = parse_integer(puls),
        sys_bp = parse_integer(sbt),
        dia_bp = parse_integer(dbt),
        creatinine = parse_integer(se_creatinin),  # µmol/L
        killip_class = killipklassevedami %>% str_extract("[IV]+"),
        arrest = genopliverefterinstitio,
        familiary_ihd,
        arrhythmia_device = arytmidevice,
        active_endocarditis = aktivendokardit_antibiot_beh_,
        lvef = ejectionfractionvalueifknown,
        ischemia_test = isk_mitest_3mdr,
        abnormal_ekg = abnormtekg,
        abnormal_qrs = abnormqrskonfig__st_t,
        pulmonal_hypertension = pulmonalhypertensiongl,
        enzymes = troponint_iel_mark_rer %>% na_if("Foreligger ikke"),
        stemi = (samletstelevationimm > 0)
                | str_detect(diagnoseefterkag, regex("I21\\.?(3|[10]B)"))
                | str_detect(aktionsdiagnose_a, regex("I21\\.?(3|[10]B)"))
    ) %>%
    distinct() %>%
    mutate_at(vars(starts_with("prev")), coalesce, "Nej") %>%
    inner_join(kag_summary, by = "entryid") %>%
    arrange(pid, proc_date)

# All diagnoses related to ischemic heart disease
ischemic_dia <- dia %>%
    filter(str_starts(C_DIAG, "DI2[0-5]")) %>%
    filter(str_detect(C_DIAGTYPE, "[ABG]"))

"Distinct patients in NPR with IHD:"
n_distinct(ischemic_dia$PID)

"Distinct patients in Pats"
n_distinct(kag$pid)

###############################################################################
### MAIN INCLUSION STEPS ######################################################

kag_debuts <- kag %>%
    filter(str_starts(native_kar, "[1-3] gebet|Ateromatose")) %T>%
    {print(glue("step1: {n_distinct(.$pid)}, ihd-pats"))} %>%
    filter(pid %in% ischemic_dia$PID) %T>%
    {print(glue("step2: {n_distinct(.$pid)}, with icd10"))} %>%
    filter(proc_date < date("2017-01-01")) %T>%
    {print(glue("step3: {n_distinct(.$pid)}, before 2017"))} %>%
    group_by(pid) %>%
    filter(all(proc_date > date("2006-01-01"))) %>%
    ungroup() %>%
    arrange(pid, proc_date) %>%
    distinct(pid, .keep_all = TRUE) %T>%
    {print(glue("step4: {n_distinct(.$pid)}, after 2006"))} %>%
    filter_at(vars(starts_with("prev")), ~ . == "Nej") %T>%
    {print(glue("step5: {n_distinct(.$pid)}, no prev"))}

patients <- kag_debuts %>%
    left_join(bds, by = c("pid" = "v_pnr_enc")) %>%
    transmute(
        pid, entryid,
        sex = recode(C_KON, "K" = "F"),
        index_type = str_extract(native_kar, "\\d") %>% coalesce("0"),
        index_date = proc_date,
        age = time_length(D_FODDATO %--% index_date, "years") %>% round(2),
        status = recode(C_STATUS, "90" = "dead", "01" = "alive"),
        status_date = coalesce(D_STATUS_HEN_START, END_OF_DATA)
    ) %>%
    filter(status %in% c("dead", "alive")) %T>%
    {print(glue("step6: {nrow(.)}, non-emigrants"))} %>%
    filter(age >= 18, index_date <= status_date) %T>%
    {print(glue("step7: {nrow(.)}, >18 years"))}

"Number of patients excluded in step6 and 7"
print(glue('{n_distinct(kag_debuts$pid) - n_distinct(patients$pid)}'))

"Type of patients included"
patients %>%
    count(index_type)

write_tsv(patients, snakemake@output[["pat"]])

### INITIAL DATA CLEANING #####################################################

# Filter and clean up diagnosis data
dia %>%
    rename_all(str_to_lower) %>%
    semi_join(patients, by = "pid") %>%
    mutate_if(~class(.x) == "integer", ~ coalesce(., 0L)) %>%
    transmute(
        pid, recnum = v_recnum,
        from_dttm = d_inddto + duration(hour = v_indtime, min = v_indminut),
        disc_dttm = d_uddto + duration(hour = v_udtime),
        icd10_code = c_diag, icd10_type = c_diagtype, icd10_extra = c_tildiag
    ) %>%
    filter(str_starts(icd10_code, "D[A-Z][0-9]")) %>%
    write_tsv(snakemake@output[["dia"]])

# Limit PATS data to included patients
kag %>%
    semi_join(patients, by = "pid") %>%
    select(-starts_with("prev")) %>%
    arrange(pid, proc_date) %>%
    write_tsv(snakemake@output[["kag"]])

# Train/test split
test_size = as.integer(snakemake@params[["test_size"]])

# For log:
print(glue("min/max age in data: {range(patients$age)}"))

# Wrangle data for survival analysis
patients %>%
    transmute(
        pid, time = time_length(index_date %--% status_date, "days"),
        event = as.integer(status == "dead"),
        train = as.integer(!row_number() %in% sample(n(), test_size)),
        sex, age
    ) %>%
    write_tsv(snakemake@output[["srv"]])
