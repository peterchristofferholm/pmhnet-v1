library(tidyverse)
library(survival)
library(patchwork)
library(riskRegression)

###############################################################################
theme_set(theme_minimal(base_size = 10))
options(
  ggplot2.discrete.colour = RColorBrewer::brewer.pal(n = 8, name = "Dark2"),
  ggplot2.discrete.fill = RColorBrewer::brewer.pal(n = 8, name = "Dark2")
)

###############################################################################

con <- DBI::dbConnect(odbc::odbc(), "transdb", database = "pmhnet")
DBI::dbExecute(con, "set search_path to pmhnet")

pgs <- tbl(con, "polygenic_scores") |> 
  distinct(pid)

surv <- tbl(con, "surv_master") |> 
  semi_join(pgs, join_by(pid)) |> 
  select(pid, train, time, event) |> 
  collect() %>%
  mutate(
    event = event == "1",
    train = train == "1"
  ) %>%
  filter(!train)

###############################################################################

load_model <- function(path) {
  .data <- read_csv(path) |> 
    pivot_longer(-pid) |> 
    mutate(across(name, as.numeric)) |> 
    summarise(score = approxfun(name, 1 - value) |> list(), .by = pid)
  structure(deframe(.data), class = "pmhnet")
}

predictRisk.pmhnet <- function(object, newdata, times) {
  ms <- object[as.character(newdata$pid)]
  ps <- map(ms, exec, times) |> unlist()
  matrix(ps, nrow = NROW(newdata), byrow = TRUE)
}

###############################################################################

fit_001 <- load_model("results/study_001/test_predictions.csv")
fit_010 <- load_model("results/study_010/test_predictions.csv")
fit_009 <- load_model("results/study_009/test_predictions.csv")
fit_008 <- load_model("results/study_008/test_predictions.csv")
fit_006 <- load_model("results/study_006/test_predictions.csv")
fit_007 <- load_model("results/study_007/test_predictions.csv")
fit_011 <- load_model("results/study_011/test_predictions.csv")

results <- Score(
  object = list(
    "001" = fit_001, "010" = fit_010,
    "009" = fit_009, "008" = fit_008,
    "006" = fit_006, "007" = fit_007,
    "011" = fit_011
  ),
  metrics = c("AUC", "Brier"),
  contrasts = list(c(2, 1), c(4, 3), c(6, 5)),
  formula   = Surv(time, event) ~ 1,
  data      = surv, 
  times     = seq(0.25, 5, by = 0.125) * 365
)

p1 <- results$AUC$score |> 
  mutate(
    with_prs = model %in% c("001", "006", "009"),
    model = model |> 
      recode("010" = "001", "008" = "009", "007" = "006") |> 
      fct_relevel("006", "009", "001", "011"), 
    times = times / 365
  ) |> 
  ggplot(aes(times, AUC, color = model)) +
  geom_line(aes(linetype = with_prs)) +
  geom_text(
    data = ~ .x |> filter(row_number(times) == 1L, .by = model),
    mapping = aes(label = model), 
    nudge_x = - 0.4
  ) +
  scale_y_continuous(labels = scales::label_percent()) +
  labs(
    x = "Prediction horizon (years)",
    y = "time-dependent AUC"
  ) +
  theme(legend.position = "none")
  
p2 <- results$AUC$contrasts |> 
  mutate(
    model = model |> case_match(
        "001" ~ "001: Diagnoses", 
        "006" ~ "006: Complete", 
        "009" ~ "009: ClinicalOne + ClinicalTwo", 
      ) |> 
      fct_relevel(
        "006: Complete",
        "009: ClinicalOne + ClinicalTwo", 
        "001: Diagnoses"
      )
  ) |> 
  ggplot(aes(times/365, delta.AUC, fill = model)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = .5) +
  geom_line() + 
  geom_hline(yintercept = 0, linetype = "dashed", alpha = .5) +
  scale_y_continuous(
    labels = scales::label_percent(), breaks = scales::pretty_breaks()
  ) +
  facet_wrap(~model, ncol = 1, scales = "free_y") +
  labs(
    x = "Prediction horizon (years)",
    y = "delta time-dependent AUC\nafter addition of PRS" 
  ) +
  theme(legend.position = "none")
  
p_combined <- p1 + p2 + plot_annotation(tag_levels = "A")
ggsave(
  "plots/misc/231003_pgs-discrimination-v2.pdf", plot = p_combined,
  width = 200, height = 100, scale = .9, unit = "mm"
)