library(tidyverse)

con <- DBI::dbConnect(odbc::odbc(), "transdb", database = "pmhnet")
DBI::dbExecute(con, "set search_path to pmhnet")

###############################################################################
# compute shap sumstats

shap_sumstats <- tbl(con, "shap_006") |>  
  mutate(
    feature = sql("substring(feature from '^[^_]+_[^_]+')")
  ) |> 
  group_by(timepoint, feature, pid) |> 
  summarise(shap = sum(shap)) |>        # per tau, feature, patient
  summarise(shap = sum(abs(shap))) |>   # per tau, feature
  mutate(shap = shap / sum(shap)) |>    # per tau
  mutate(
    rank = dense_rank(desc(shap)),
    fgroup = sql("substring(feature from '^[^_]+')"),
    fname  = sql("substring(feature from '^[^_]+_([^_]+)')")
  ) |>
  mutate(
    fgroup = if_else(fgroup == "diag-simple", "diag-1", fgroup)
  ) |> 
  select(-feature) |> 
  ungroup() |> collect() |>
  mutate(
    timepoint = timepoint |> as.character() |> 
      recode("183" = "6m", "365" = "1y", "1095" = "3y", "1825" = "5y") |> 
      fct_inorder()
  ) |> 
  select(timepoint, fgroup, fname, shap)

exports_dir <- "/users/projects/pmhnet/exports/221116/"
write_tsv(
  x = shap_sumstats,
  file = fs::path(exports_dir, "shap-sumstats.tsv.xz")
)
  
###############################################################################

colors <- c(
  "bioc-1" = "#ebb391",
  "diag-1" = "#7ba79d",
  "clnc-1" = "#e8d19d",
  "clnc-2" = "#bfd3e4",
  "proc-1" = "#cadcd8"
)

colors |> 
  enframe("fgroup", "color") |> 
  write_tsv(fs::path(exports_dir, "colors.tsv"))