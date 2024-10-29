library(tidyverse)
library(survival)
library(patchwork)

theme_set(theme_minimal(base_size = 10, base_family = "Helvetica"))
options(
  ggplot2.discrete.colour = RColorBrewer::brewer.pal(n = 8, name = "Dark2"),
  ggplot2.discrete.fill = RColorBrewer::brewer.pal(n = 8, name = "Dark2")
)

con <- DBI::dbConnect(odbc::odbc(),
    driver = "PostgreSQL", server = "dbserver", database = "database",
    uid = Sys.getenv("DBUSER"), pwd = Sys.getenv("DBPASS")
)

yhat_train <- read_csv("results/study_006/train_predictions.csv") %>% 
  pivot_longer(
    -pid, names_to = "time", values_to = "prob", 
    names_transform = list(time = as.numeric)
  )

yhat_valid <- read_csv("results/study_006/test_predictions.csv") %>% 
  pivot_longer(
    -pid, names_to = "time", values_to = "prob", 
    names_transform = list(time = as.numeric)
  )

yhat <- bind_rows(
    training = yhat_train, validation = yhat_valid, .id = "split"
  ) %>% 
  group_by(split, pid) %>%
  summarise(
    tau = 365 * c(0.5, 1, 3, 5),
    yhat = approx(x = time, y = prob, xout = tau)$y,
    .groups = "drop"
  )

################################################################################
# plid: 210716_01
# distribution of predictions on training data
################################################################################

(p_riskgroups <- yhat %>% 
  filter(split == "training", tau == 1825L) %>% 
  ggplot(aes(yhat)) +
  geom_histogram(aes(y = after_stat(density * width)), boundary = 1, binwidth = 1/40) +
  geom_vline(
    aes(xintercept = x),
    data = tibble(x = c(0.5, 0.75, 0.875, 0.95)), size = 1/4
  ) +
  scale_x_continuous(
    labels = scales::percent_format(accuracy = 1), n.breaks = 6
  ) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(
    title = "Distribution of model predictions",
    x = "5-year survival prediction (%)", y = "Proportion of patients (%)"
  )
)

ggsave(
  "results/study_006/210712/01_yhat-histogram.pdf", plot = p_riskgroups,
  height = 6, width = 7, unit = "cm", scale = 1.3
)

################################################################################

surv <- tbl(con, dbplyr::in_schema("pmhnet", "surv_master")) %>% 
  collect() %>% 
  mutate(
    across(everything(), as.integer),
    event = if_else(time > 1825, 0L, event),
    time = if_else(time > 1825, 1825L, time)
  )

categories <- yhat %>% 
  filter(tau == 1825L) %>% 
  mutate(
    risk_group = cut(yhat, c(0, 0.5, 0.75, 0.875, 0.95, 1))
  ) %>% 
  select(pid, risk_group)

count(categories, risk_group)

################################################################################
# plid: 210716_02
# average predicted vs observed KM-estimates stratified by risk-groups
################################################################################

bootstrap_mean <- function(column, nboot = 10) {
  stat <- rlang::as_function(~mean(.x[.y]))
  bo <- boot::boot(column, statistic = stat, R = nboot)
  broom::tidy(bo, conf.int = TRUE, conf.method = "perc")
}

ave_yhat_groups_valid <- yhat_valid %>% 
  inner_join(categories, by = "pid") %>% 
  group_by(risk_group, time) %>% 
  summarise(bootstrap_mean(prob, nboot = 200)) %>% 
  rename(estimate = statistic)

pdat <- surv %>% 
  filter(train == 0L) %>% 
  inner_join(categories, by = "pid") %>% 
  group_by(risk_group) %>% 
  summarise(
    survfit(Surv(time, event) ~ 1, data = cur_data()) %>% broom::tidy()
  ) %>% 
  bind_rows(ave_yhat_groups_valid, .id = "type") %>%
  mutate(type = recode(type, `1` = "observed", `2` = "predicted")) 

(p2 <- pdat %>% 
  ggplot(aes(time/365.25, estimate)) +
  geom_ribbon(
    aes(ymin = conf.low, ymax = conf.high, fill = risk_group, linetype = type), 
    alpha = 1/5
  ) +
  geom_line(aes(linetype = type, color = risk_group)) +
  geom_text(
    aes(x = 5, y = y, label = risk_group, color = risk_group),
    data = ~ filter(.x, type == "predicted") %>% summarise(y = min(estimate)),
    hjust = 0, nudge_x = 1/20, size = 3
  ) +
  scale_y_continuous(
    labels = scales::percent_format(accuracy = 1), n.breaks = 6
  ) +
  labs(
    title = "Predictions compared to KM-estimates",
    y = "Survival (%)", x = "Time (years)", linetype = ""
  ) +
  coord_cartesian(clip = "off") +
  guides(
    fill = "none", color = "none", label = "none"
  ) +
  theme(
    plot.margin = unit(c(5, 55, 5, 5), "pt"), legend.position = "bottom"
  )
)
  
ggsave(
  "results/study_006/210712/02_yhat-vs-km.pdf", plot = p2,
  height = 6, width = 9, unit = "cm", scale = 1.3
)

################################################################################
# plid: 210719_03 
# p1 + p2 combined
################################################################################

p2 + p_riskgroups + 
  plot_layout(widths = c(9, 7)) + 
  plot_annotation(tag_levels = "A")

ggsave("results/study_006/210712/03_combined.pdf",
       height = 6, width = 17, unit = "cm", scale = 1.3)

################################################################################
# plid: 210719_001
# calibration score
################################################################################

(p3 <- surv %>%
  inner_join(yhat, by = "pid") %>% 
  filter(train == 0L, tau == 1825L) %>% 
  ggplot(aes(x = 1 - yhat, y = event)) +
  annotate("segment", x = 0, y = 0, xend = 1, yend = 1, linetype = "dashed") +
  geom_jitter(height = 1/10, size = 1/2, alpha = .5, stroke = 0) +
  geom_smooth(method = "loess", span = 0.5) +
  scale_y_continuous(breaks = 0:1, labels = c("censored", "dead")) +
  xlim(0, 1) + coord_equal() +
  labs(
    x = "Predicted risk at 5 years",
    y = "Mortality at 5 years"
  )
)

ggsave(
  "results/study_006/210712/04_calibration.pdf", plot = p3,
  height = 6, width = 6, unit = "cm", scale = 1.3
)

# calibration score is area between loess curve and the "perfect" line. Conf.
# intervals can be obtained through bootstrap resampling. 

calibration_bt <- function(data, nboot = 10) {
  data <- distinct(data, yhat, .keep_all = TRUE)  # no duplicates allowed
  calibration <- function(data, i) {
    res <- slice(data, i) %>% 
      mutate(
        loess = loess(event ~ yhat, data = cur_data()) %>% predict()
      ) %>% 
      distinct(yhat, loess) %>% 
      summarise(
        alpha = approxfun(yhat, yhat - loess) %>% 
          integrate(lower = min(yhat), upper = max(yhat)) %>% 
          `$`("val")
      )
    1 - res$alpha
  }
  bo <- boot::boot(data, statistic = calibration, R = nboot)
  broom::tidy(bo, conf.int = TRUE, conf.method = "perc")
}

surv %>%
  inner_join(yhat, by = "pid") %>% 
  filter(train == 0L, tau == 1825L) %>% 
  mutate(yhat = 1 - yhat) %>% 
  summarise(
    calibration_bt(cur_data(), nboot = 100)
  )
  
################################################################################
# concordance index and Somers' Dxy rank correlation
################################################################################

inner_join(surv, yhat, by = "pid") %>% 
  mutate(
    event = if_else(time > tau, 0L, event),
    time  = if_else(time > tau, tau, as.numeric(time))
  ) %>% 
  group_by(split, tau) %>% 
  summarise(
    Hmisc::rcorr.cens(yhat, Surv(time, event)) %>% enframe(), .groups = "drop"
  ) %>% 
  pivot_wider(names_from = name, values_from = value)

################################################################################
# time-dependent AUC with tdROC package
################################################################################

tdauc_bt <- function(data, nboot) {
  tdauc <- function(data, i) {
    ttau <- unique(data$tau)
    assertthat::assert_that(length(ttau) == 1L)
    slice(data, i) %>% 
      summarise(
        result = list(
          tdROC::tdROC(X = -yhat, Y = time, delta = event, tau = ttau)
        )
      ) %>% 
      pluck("result", 1, "AUC", "value")
  }
  bo <- boot::boot(data, statistic = tdauc, R = nboot)
  broom::tidy(bo, conf.int = TRUE, conf.method = "perc")
}

tdauc_validation_bootstrap <- surv %>% 
  inner_join(yhat, by = "pid") %>% 
  filter(split == "validation") %>% 
  group_by(tau) %>%
  summarise(
    tdauc_bt(cur_data_all(), nboot = 100)
  )

tdauc_validation_bootstrap
