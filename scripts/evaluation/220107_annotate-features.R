library(tidyverse)

con <- DBI::dbConnect(
    odbc::odbc(), driver = "PostgreSQL", server = "trans-db-01", 
    database = "pmhnet", uid = Sys.getenv("DBUSER"),
    pwd = Sys.getenv("DBPASS"), bigint = "integer"
)

## prepare classification ######################################################

# mapping between npu code and component names
tbl(con, dbplyr::ident_q("core.biochem")) %>% 
  distinct(code = quantity_id, analyte = component) %>% 
  collect() -> bioc_map

# icd10 code mappings
read_tsv("data/resources/icd10_definitions.tsv") %>% 
  mutate(across(icd10, ~str_c("D", .x))) %>% 
  select(icd10, definition) -> icd10_map

# procedure code mappings
read_tsv("data/resources/proc_definitions.tsv") %>% 
  distinct(code, .keep_all = TRUE) -> proc_map


## categorise features #########################################################

features <- read_table(
    "results/study_006/features.txt", col_names = "feature", col_types = cols()
  ) %>% 
  separate(
    feature, into = c("fgroup", "fname", "fvalue"), sep = "_", fill = "right"
  ) %>% 
  mutate(across(fgroup, recode, "diag-simple" = "diag-1")) %>% 
  group_by(fgroup, fname) %>% 
  summarise(
    values = str_c(fvalue, collapse = "/"), .groups = "drop"
  )

types <- read_tsv("data/scratch/006-input.tsv", col_types = cols()) %>% 
  summarise(across(-pid, ~ n_distinct(.x) > 2)) %>% 
  pivot_longer(everything(), values_to = "continuous", names_to = "fname") %>% 
  mutate(across(fname, str_extract, "(?<=_)[^_]+")) %>% 
  group_by(fname) %>% 
  summarise(
    type = case_when(
        n() > 1 ~ "categorical", continuous ~ "continuous", TRUE ~ "logical"
    ) %>% unique(), 
    .groups = "drop"
  )

features %>% 
  inner_join(types, by = "fname") %>% 
  left_join(bioc_map, by = c("fname" = "code")) %>% 
  left_join(icd10_map, by = c("fname" = "icd10")) %>% 
  left_join(proc_map, by = c("fname" = "code")) %>% 
  mutate(
    definition = coalesce(analyte, definition, text), .keep = "unused"
  ) %>% 
  write_tsv("data/scratch/feature-list.tsv")
