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

con <- DBI::dbConnect(odbc::odbc(), "dbserver", database = "pmhnet")
DBI::dbExecute(con, "set search_path to pmhnet")

surv <- tbl(con, "surv_master") %>%
  select(pid, train, time, event) %>%
  collect() %>%
  mutate(
    event = event == "1",
    train = train == "1"
  ) %>%
  filter(!train)

###############################################################################

get_predictions <- function(path) {
  out <- path %>%
    read_csv(col_types = cols()) %>%
    pivot_longer(
      -pid, names_to = "timepoint", values_to = "estimate",
      names_transform = list(timepoint = as.numeric)
    ) %>%
    mutate(
      estimate = 1 - estimate,
      pid = as.integer(pid)
    )
  class(out) <- c(class(out), "pmhnet")
  return(out)
}

predictRisk.pmhnet <- function(object, newdata, times, ...) {
  estimates <- newdata %>%
    inner_join(object, by = "pid") %>%
    group_by(pid) %>%
    summarise(
      tau = times,
      val = approx(timepoint, estimate, xout = tau)$y,
      .groups = "drop"
    ) %>%
    pivot_wider(names_from = tau, values_from = val)

  newdata %>% select(pid) %>%
    left_join(estimates, by = "pid") %>%
    select(-pid) %>% data.matrix()
}

###############################################################################

fit1 <- get_predictions("results/study_002/test_predictions.csv")
fit2 <- get_predictions("results/study_006/test_predictions.csv")

###############################################################################

# GRACE risk score
fit0 <- tbl(con, "risk_scores") %>%
  select(pid, timepoint = time, estimate = value) %>%
  collect() %>%
  mutate(
    timepoint = if_else(timepoint == 183, 182.5, as.numeric(timepoint)),
    estimate = estimate / 100
  ) %>%
  semi_join(fit2, by = "pid")

class(fit0) <- c(class(fit0), "pmhnet")

###############################################################################

results <- Score(
  object = list("GRACE2.0" = fit0, "model_2" = fit1, "model_6" = fit2),
  formula = Surv(time, event) ~ 1,
  data = surv, times = c(365/2, 365*1, 365*3),
  plots = c("ROC", "Calibration")
)

###############################################################################
# 220516_calibration-plot_pmhnet-vs-grace.pdf

get_calibration_data <- function(results, time) {
  plotCalibration(results, cens.method = "local", plot = F, times = time) %>%
    pluck("plotFrames") %>%
    map_dfr(~as_tibble(.x, .name_repair = str_to_lower), .id = "model")
}

p1 <- tibble(tau = c("6m", "1y", "3y"), times = c(365/2, 365*1, 365*3)) %>%
  summarise(
    tau = as_factor(tau),
    val = map(times, get_calibration_data, results=results)
  ) %>%
  unnest(val) %>%
  ggplot(aes(pred, obs, color = model)) +
  geom_abline(slope = 1, linetype = "dashed", alpha = .5) +
  geom_line(size = 0.9) +
  scale_x_continuous(limits = c(0, 1), labels = scales::label_percent()) +
  scale_y_continuous(limits = c(0, 1), labels = scales::label_percent()) +
  facet_wrap(~tau) +
  coord_equal() +
  labs(x = "Predicted risk", y = "Estimated actual risk", color = "") +
  theme(legend.position = "bottom")

ggsave(
  "plots/misc/220516_calibration-plot_pmhnet-vs-grace.pdf", plot = p1,
  width = 190, height = 100, scale = .9, unit = "mm"
)

###############################################################################
# 220516_discrimination-plot_pmhnet-vs-grace.pdf

plot_roc_curves <- function(score_obj) {
  auc <- results %>% pluck("AUC", "score") %>% as_tibble() %>%
    mutate(
      timepoint = as_factor(times) %>%
        lvls_revalue(c("6m", "1y", "3y")),
      across(c(AUC, lower, upper), ~format(round(.x, digits = 2), nsmall = 2)),
      label = glue::glue(
        "{AUC} [{lower};{upper}]"
      ),
      ypos = (seq(0, 0.2, length.out = 3)[as.integer(model)-1])
    )
  results %>% pluck("ROC", "plotframe") %>% as_tibble() %>%
    mutate(
      timepoint = as_factor(times) %>%
        lvls_revalue(c("6m", "1y", "3y"))
    ) %>%
    ggplot(aes(x = FPR, y = TPR, color = model)) +
    geom_line(size = 0.9) +
    geom_label(
      data = auc,
      aes(x = 1, y = ypos, label = label, fill = model),
      color = "white", size = 2.5, nudge_x = -0.2, nudge_y = 0.03
    ) +
    geom_abline(slope = 1, alpha = .5, linetype = "dashed") +
    facet_wrap(~timepoint) +
    scale_x_continuous(limits = c(0, 1), labels = scales::label_percent()) +
    scale_y_continuous(limits = c(0, 1), labels = scales::label_percent()) +
    coord_equal() +
    labs(x = "1 - Specificity", y = "Sensitivity", color = "") +
    theme(legend.position = "none")
}

p2 <- plot_roc_curves(results)

ggsave(
  "plots/misc/220516_discrimination-plot_pmhnet-vs-grace.pdf", plot = p2,
  width = 190, height = 100, scale = .9, unit = "mm"
)

###############################################################################

comparisons <- list("GRACE2.0" = fit0, "model_6" = fit2, "model_2" = fit1) %>%
  map(predictRisk, surv, times = c(365*3)) %>%
  map(as.double) %>%
  map_df(~as_tibble(.x) %>% rownames_to_column("id"), .id = "model") %>%
  pivot_wider(names_from = model, values_from = value)

create_comparison_plot <- function(data, x, y) {
  ggplot(data, aes_(x = substitute(x), y = substitute(y))) +
    stat_bin2d(
      aes(color = after_stat(count)),
      geom = "point", binwidth = 1/50, size = .5, fill = "black"
    ) +
    geom_abline(slope = 1, alpha = .5, linetype = "dashed") +
    geom_density2d(
      contour_var = "ndensity", breaks = c(5/100, 10/100, 20/100), color = "black"
    ) +
    scale_x_continuous(limits = c(0, 1), labels = scales::label_percent()) +
    scale_y_continuous(limits = c(0, 1), labels = scales::label_percent()) +
    scale_color_viridis_c(
      trans = "log10", limits = c(1, 500), oob = scales::squish
    ) +
    coord_equal()
}

p3 <- create_comparison_plot(comparisons, GRACE2.0, model_6) +
  create_comparison_plot(comparisons, model_2,  model_6) +
  plot_layout(guides = "collect")

ggsave(
  "plots/misc/220516_model-comparison-plot.pdf", plot = p3,
  width = 190, height = 100, scale = .9, unit = "mm"
)

###############################################################################
# extended results for all versions of pmhnet

models <- fs::dir_ls(
    glob = "results/study_*/test_predictions.csv", recurse = TRUE
  ) %>%
  enframe("model", "path") %>%
  mutate(across(model, ~ str_match(., "study_(\\d{3})")[,2])) %>%
  filter(as.integer(model) %in% 1:6) %>%
  deframe() %>% map(get_predictions)

results_extended <- Score(
  object = models,
  formula = Surv(time, event) ~ 1,
  data = surv, times = seq(0.5, 5, by = 0.125) * 365
)

(p4 <- results_extended %>%
  pluck("AUC", "score") %>%
  as_tibble() %>% set_names(str_to_lower) %>%
  ggplot(aes(times / 365, auc)) +
  geom_line(aes(color = model)) +
  geom_text(
    data = ~ filter(.x, times == 365/2),
    aes(label = model, color = model),
    nudge_x = - 1/10
  ) +
  scale_y_continuous(
    labels = scales::label_percent(),
    minor_breaks = scales::breaks_width(0.01)
) +
  scale_x_continuous(limits = c(0, 5)) +
  labs(
    title = "Model discrimination over time",
    x = "Prediction horizon (years)", y = "Time-dependent AUC"
  ) +
  theme(legend.position = "none")
)

ggsave(
  "plots/misc/220516_pmhnet-auc-over-time.pdf", plot = p4,
  width = 190, height = 100, scale = .9, unit = "mm"
)

capture.output(print(results_extended, nrows = 1000)) %>%
  write_lines("stats/220516_pmhnet-model-comparisons.txt")
