library(tidyverse)
library(survival)
library(patchwork)
library(riskRegression)
library(glue)

###############################################################################
sprintf_transformer <- function(text, envir) {
  m <- regexpr(":.+$", text)
  if (m != -1) {
    format <- substring(regmatches(text, m), 2)
    regmatches(text, m) <- ""
    res <- eval(parse(text = text, keep.source = FALSE), envir)
    do.call(sprintf, list(glue("%{format}"), res))
  } else {
    eval(parse(text = text, keep.source = FALSE), envir)
  }
}

glue_fmt <- function(..., .envir = parent.frame()) {
  glue(..., .transformer = sprintf_transformer, .envir = .envir)
}
###############################################################################

con <- DBI::dbConnect(odbc::odbc(), "transdb", database = "pmhnet")
DBI::dbExecute(con, "set search_path to pmhnet")

# observed survival of the test set
surv <- tbl(con, "surv_master") |> 
  select(pid, train, time, event) |> 
  collect() |> 
  mutate(
    event = event == "1",
    train = train == "1"
  ) |> 
  filter(!train)

# feature attributions (shap)
shap <- tbl(con, "shap_013") |>  
  mutate(
    feature = sql("substring(feature from '^[^_]+_[^_]+')")
  ) |> 
  filter(timepoint == 1825L) |> 
  group_by(feature, pid) |> 
  summarise(shap = sum(shap)) |> 
  summarise(shap = mean(abs(shap)) * 100) |>
  mutate(shap_rank = dense_rank(desc(shap))) |> 
  collect()

# model prediction with different features left out
data <- read_csv(
    "results/study_006/experiments/sparse-predictions.csv",
    show_col_types = FALSE 
  ) |> 
  rename(pid = `...1`) |> 
  transmute(
    pid = as.integer(pid), 
    fname = fname |> coalesce("full"), 
    pred = `5y`
  ) |> 
  nest(model = c(pid, pred)) |> 
  deframe()

###############################################################################
 
predictRisk.data.frame <- function(object, newdata, times, ...) {
  object |> 
    right_join(
      newdata |> select(pid), join_by(pid)
    ) |> 
    arrange(match(pid, newdata$pid)) |> 
    transmute(est = 1 - pred) |> 
    data.matrix()
}

results <- Score(
  object = data,
  summary = c("ipa"),
  formula = Surv(time, event) ~ 1,
  contrasts = list(1:585),
  times = 5 * 365,
  data = surv |> select(pid, time, event)
)

###############################################################################

results[["AUC"]][["contrasts"]]  |> 
  as_tibble() |> 
  mutate(
    feature = model |> as.character(),
    p_adjusted = p.adjust(p, "BH"),
  ) |> 
  left_join(shap, join_by(feature)) |> 
  mutate(across(c(delta.AUC, lower, upper), \(x) x * 100)) |> 
  arrange(delta.AUC) |> 
  transmute(
    feature_missing = feature,
    delta_auc = glue_fmt(
      "{delta.AUC:+.1f} [{lower:+.1f}; {upper:+.1f}]"
    ),
    shap_rank, p_adjusted
  )

results[["Brier"]][["contrasts"]]  |> 
  as_tibble() |> 
  mutate(
    feature = model |> as.character(),
    p_adjusted = p.adjust(p, "BH"),
  ) |> 
  left_join(shap, join_by(feature)) |> 
  mutate(across(c(delta.Brier, lower, upper), \(x) x * 100)) |> 
  arrange(desc(delta.Brier)) |> 
  transmute(
    feature_missing = feature,
    delta_brier = glue_fmt(
      "{delta.Brier:+.2f} [{lower:+.1f}; {upper:+.1f}]"
    ),
    shap_rank, p_adjusted
  )
