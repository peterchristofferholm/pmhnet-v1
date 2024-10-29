library(tidyverse)
library(survival)
library(patchwork)
library(riskRegression)

###############################################################################
theme_set(theme_minimal(base_size = 10, base_family = "Helvetica"))
options(
  ggplot2.discrete.colour = RColorBrewer::brewer.pal(n = 8, name = "Dark2"),
  ggplot2.discrete.fill = RColorBrewer::brewer.pal(n = 8, name = "Dark2")
)

###############################################################################

con <- DBI::dbConnect(odbc::odbc(), "transdb", database = "pmhnet")
DBI::dbExecute(con, "set search_path to pmhnet")

surv <- tbl(con, "surv_master") |> 
  select(pid, train, time, event) |> 
  collect() |> 
  mutate(
    event = event == "1",
    train = train == "1"
  ) |> 
  filter(!train)

###############################################################################

get_predictions <- function(path) {
  path |> 
    read_csv(show_col_types = FALSE) |> 
    rowwise() |> 
    transmute(
      pid = pid,
      model = 
        approxfun(
          x = pick(-pid) |> colnames() |> as.numeric(),
          y = c_across(-pid) 
        ) |> list()
    )
}

predictRisk.data.frame <- function(object, newdata, times, ...) {
  estimates <- object |> 
    semi_join(newdata, by = "pid") |> 
    rowwise() |> 
    reframe(
      pid = pid,
      tau = times,
      est = 1 - model(times) 
    ) |> 
    pivot_wider(names_from = tau, values_from = est)
  
  newdata |> 
    select(pid) |> 
    left_join(estimates, by = "pid") |> 
    select(-pid) |> data.matrix()
}

###############################################################################

fit1 <- get_predictions("results/study_002/test_predictions.csv")
fit2 <- get_predictions("results/study_006/test_predictions.csv")
fit3 <- get_predictions("results/study_013/test_predictions.csv")

res <- Score(
  object = list("002" = fit1, "006" = fit2, "013" = fit3), 
  summary = c("ipa"), plots = c("ROC", "calibration"),
  formula = Surv(time, event) ~ 1,
  times = c(0.5, 1, 3, 5) * 365,
  data = surv |> select(pid, time, event)
)

res[["AUC"]][["score"]] |> 
  as_tibble()

res[["Brier"]][["score"]] |> 
  as_tibble()
