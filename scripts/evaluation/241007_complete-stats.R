library(tidyverse)  
library(dbplyr)

features <- read_lines("data/scratch/006-input.tsv", n_max = 1) |> 
  str_split_1("\t") |>  
  map_chr(\(x) str_split_1(x, "_")[2]) |> 
  unique()

con <- DBI::dbConnect(odbc::odbc(), "transdb", database = "pmhnet")

splits <- tbl(con, in_schema("pmhnet", "surv_master")) |> 
  collect() |> transmute(pid, split = c("test", "train")[as.integer(train) + 1])

###############################################################################
## biochem

bioc <- tbl(con, in_schema("pmhnet", "biochem")) |> 
  filter(time |> between(0, 3650)) |> 
  slice_min(order_by = time,  by = c(pid, qid),  with_ties = FALSE, n = 1) |> 
  select(pid, qid, value) |> 
  pivot_wider(names_from = qid, values_from = value) |> 
  collect() |> 
  mutate(across(-pid, \(x) c("low", "normal", "high")[(x + 2)]))

bioc_stats <- bioc |> 
  right_join(splits)  |> 
  mutate(across(-c(pid, split), \(x) fct_na_value_to_level(x, "missing"))) |> 
  pivot_longer(-c(pid, split)) |> 
  filter(name %in% features) |> 
  count(split, name, value) |> 
  mutate(n = if_else(n < 5, 5, n)) |> 
  mutate(
    p = n * 100 / sum(n), .by = c(split, name)
  ) |> 
  arrange(split, name)
  
###############################################################################
## diagnoses

diag <- tbl(con, in_schema("pmhnet", "diagnoses")) |> 
  filter(time |> between(-3650, 0)) |> 
  filter(level == "lvl3") |> 
  distinct(pid, code) |> 
  collect() 

diag_stats <- diag |> 
  filter(code %in% features) |> 
  mutate(code = str_sub(code, start = 2L)) |> 
  add_column(seen = TRUE) |>
  right_join(splits) |> 
  complete(nesting(pid, split), code, fill = list(seen = FALSE))  |> 
  filter(!is.na(code)) |> 
  count(split, code, seen) |> 
  mutate(n = if_else(n < 5, 5, n)) |> 
  mutate(p = n * 100 / sum(n), .by = c(code, split)) |> 
  arrange(split, code)


###############################################################################
## procedures

proc <- tbl(con, in_schema("pmhnet", "nomesco")) |> 
  filter(time |> between(0, 3650)) |> 
  distinct(pid, code) |> 
  collect()

proc_stats <- proc |> 
  filter(code %in% features) |> 
  add_column(seen = TRUE) |>
  right_join(splits) |> 
  complete(nesting(pid, split), code, fill = list(seen = FALSE))  |> 
  filter(!is.na(code)) |> 
  count(split, code, seen) |> 
  mutate(n = if_else(n < 5, 5, n)) |> 
  mutate(p = n * 100 / sum(n), .by = c(code, split)) |> 
  arrange(split, code)

###############################################################################
## clinical

# note: joining on 'na' values doesn't work in postgres

clin <- inner_join(
  tbl(con, in_schema("pmhnet", "baseline")) |> collect(),
  tbl(con, in_schema("pmhnet", "clinicaltwo")) |> collect()
)

clin_stats_1 <- clin |> 
  inner_join(splits, join_by(pid)) |> 
  summarise(
    across(
      c(age, height, weight, pulse, crea, sys_bp, dia_bp),
      function(x) {
        tibble(
          mean    = mean(x, na.rm = TRUE),
          sd = sd(x, na.rm = TRUE),
          missing = sum(is.na(x)) * 100 / length(x)
        )
      }
    ),
    .by = split
  ) |> 
  pivot_longer(-split) |> 
  unpack(value)
  

clin_stats_2 <- clin |> 
  inner_join(splits, join_by(pid)) |> 
  select(-pid) |> 
  select(-c(age, height, weight, pulse, crea, sys_bp, dia_bp)) |> 
  mutate(
    sex = fct_recode(sex, male = "1", female = "0"),
    across(
      c(arrest, stemi, icd_or_pm, enzymes),
      \(x) x |> 
        fct_recode(yes = "1", no = "0") |> 
        fct_na_value_to_level("missing")
    ),
    across(
      where(is_character) | c(nyha, ccs, killip), 
      \(x) x |> 
        as.character() |> 
        fct_na_value_to_level("missing")
    )
  ) |> 
  pivot_longer(-split) |> 
  count(split, name, value) |> 
  mutate(p = n * 100 / sum(n), .by = c(name, split))


## save

diag_stats   |> write_tsv("results/misc/241007/stats_diagnoses.tsv")
bioc_stats   |> write_tsv("results/misc/241007/stats_biochemical.tsv")
proc_stats   |> write_tsv("results/misc/241007/stats_procedures.tsv")
clin_stats_1 |> write_tsv("results/misc/241007/stats_clinical_num.tsv") 
clin_stats_2 |> write_tsv("results/misc/241007/stats_clinical_fct.tsv")