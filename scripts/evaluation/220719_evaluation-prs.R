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
  semi_join(tbl(con, "polygenic_scores"), by = "pid") |> 
  collect() |> 
  mutate(
    event = event == "1",
    train = train == "1"
  ) |> 
  filter(!train)

###############################################################################

get_predictions <- function(path) {
  out <- path |> 
    read_csv(col_types = cols()) |> 
    pivot_longer(
      -pid, names_to = "timepoint", values_to = "estimate", 
      names_transform = list(timepoint = as.numeric)
    ) |> 
    mutate(
      estimate = 1 - estimate,
      pid = as.integer(pid)
    )
  class(out) <- c(class(out), "pmhnet")
  return(out)
}

predictRisk.pmhnet <- function(object, newdata, times, ...) {
  estimates <- newdata |> 
    inner_join(object, by = "pid") |> 
    group_by(pid) |> 
    summarise(
      tau = times,
      val = approx(timepoint, estimate, xout = tau)$y,
      .groups = "drop"
    ) |> 
    pivot_wider(names_from = tau, values_from = val)
  
  newdata |> select(pid) |> 
    left_join(estimates, by = "pid") |> 
    select(-pid) |> data.matrix()
}

###############################################################################

models <- fs::dir_ls(
    glob = "results/study_*/test_predictions.csv", recurse = TRUE
  ) |> 
  enframe("model", "path") |> 
  mutate(across(model, ~ str_match(., "study_(\\d{3})")[,2])) |> 
  filter(as.integer(model) %in% c(1, 8, 9, 10, 11)) |> 
  deframe() |> map(get_predictions)
  
results_extended <- Score(
  object = models, contrasts = list(c(3, 2), c(1, 4)),
  formula = Surv(time, event) ~ 1,
  data = surv, times = seq(0.25, 5, by = 0.125) * 365
)

(p1 <- results_extended |> 
  pluck("AUC", "score") |>
  as_tibble() |> set_names(str_to_lower) |> 
  ggplot(aes(times / 365, auc)) +
  geom_line(aes(color = model)) +
  geom_text(
    data = ~ filter(.x, times == min(times)),
    aes(label = model, color = model),
    nudge_x = - 1/10
  ) +
  scale_y_continuous(
    labels = scales::label_percent(), 
    minor_breaks = scales::breaks_width(0.01)
  ) +
  scale_x_continuous(limits = c(0, 5)) +
  labs(
    x = "Prediction horizon (years)", y = "Time-dependent AUC"
  ) +
  theme(legend.position = "none")
)

ggsave(
  "plots/misc/220719_pgs-discrimination.pdf", plot = p1,
  width = 170, height = 100, scale = .9, unit = "mm"
)

###############################################################################

(p2 <- results_extended$AUC$contrasts |> 
  as_tibble() |> 
  mutate(
    reference = reference |> 
      recode(
        "001" = "Diagnoses", 
        "009" = "ClinicalOne + ClinicalTwo"
      )
  ) |> 
  ggplot(aes(x = times/365, y = delta.AUC, fill = reference)) +
  geom_ribbon(
    aes(ymin = lower, ymax = upper), alpha = .25
  ) +
  geom_line() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  facet_wrap(~reference, ncol = 1) +
  theme(legend.position = "none") +
  scale_x_continuous(limits = c(0, 5)) +
  scale_y_continuous(labels = scales::label_percent()) +
  labs(
    title = "Change in tdAUC after adding PGS data",
    x = "Prediction horizon (years)", 
    y = "delta-tdAUC"
  )
)

ggsave(
  "plots/misc/220719_pgs-delta-auc.pdf", plot = p2,
  width = 170, height = 100, scale = .9, unit = "mm"
)
