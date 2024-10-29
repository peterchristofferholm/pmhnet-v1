library(tidyverse)
library(survival)
library(patchwork)
library(dbplyr)

theme_set(theme_minimal(base_size = 10, base_family = "NimbusSan"))
options(
  ggplot2.discrete.colour = RColorBrewer::brewer.pal(n = 8, name = "Dark2"),
  ggplot2.discrete.fill = RColorBrewer::brewer.pal(n = 8, name = "Dark2")
)

# get observed survival and model predictions
conn <- DBI::dbConnect(odbc::odbc(), "transdb", database = "pmhnet")
surv <- tbl(conn, in_schema("pmhnet", "surv_master"))
pred <- read_csv("results/study_006/test_predictions.csv")

pred <- pred %>% 
  pivot_longer(
    -pid, names_to = "time", values_to = "prob", 
    names_transform = list(time = as.numeric)
  ) %>% 
  group_by(pid) %>% 
  summarise(
    yhat = approx(x = time, y = prob, xout = 365*5)$y
  )

intervals <- c(0, 25, 50, 75, 90, 100) / 100

###############################################################################
# histogram of model predictions

(p1 <- pred %>% 
  group_by(
    group = cut(yhat, intervals),
    level = as.integer(group)
  ) %>% 
  ggplot(aes(x = yhat)) +
  geom_histogram(
    aes(y = after_stat(density * width)),
    boundary = 1, binwidth = 0.025
  ) +
  geom_label(
    aes(x = x, label = n, y = y, fill = group),
    data = ~ .x %>% summarise(
      x = (intervals[level[1]] + intervals[level[1] + 1])/2,
      n = scales::label_percent()(n() / 5000),
      y = .3,  .groups = "drop"
    ),
    alpha         = .5,
    label.r       = unit(1/10, "lines"),
    label.size    = unit(1/2,  "lines"),
    label.padding = unit(1/4,  "lines"),
    color         = "white"
  ) +
  geom_vline(
    aes(xintercept = x),
    data = tibble(x = intervals),
    linetype = "dashed", alpha = .5
  ) +
  scale_x_continuous(
    labels = scales::label_percent(),
    limits = c(0, 1)
  ) +
  scale_y_continuous(
    labels = scales::label_percent(), expand = expansion(add = c(0, 0.05))
  ) +
  labs(
    x = "5-year survival prediction (%)", y = "Proportion of patients (%)"
  ) +
  coord_flip(clip = "off") +
  theme(
    legend.position = "none"
  )
)
  
###############################################################################
# survival curves for each risk strata

surv <- surv %>% collect() %>% 
  inner_join(pred, by = "pid", copy = TRUE) %>% 
  mutate(
    risk_strata = cut(yhat, c(0, 25, 50, 75, 90, 100) / 100)
  )

(p2 <- surv %>% 
  group_by(risk_strata) %>% 
  summarise(
    list(survfit(Surv(time, event == "1") ~ 1, data = cur_data())) %>% 
      map(summary, times = seq(0, 365 * 5, by = 1)) %>% 
      map_dfr(`[`, c("time", "surv", "lower", "upper")) %>% 
      add_case(surv = 1, time = 0, lower = 1, upper = 1, .before = 1)
    , .groups = "drop"
  ) %>% 
  ggplot(aes(time / 365, surv, group = risk_strata)) +
  geom_ribbon(
    aes(ymin = lower, ymax = upper, fill = risk_strata),
    alpha = 1/4
  ) +
  geom_line(
    aes(color = risk_strata) 
  ) +
  scale_y_continuous(
    labels = scales::percent_format(accuracy = 1), n.breaks = 5,
    limits = c(0, 1)
  ) +
  labs(
    x = "Follow-up time (years)",
    y = "Estimated actual survival"
  ) + 
  theme(
    legend.position = "bottom"
  )
)

(p_assembled <- (p2 + p1 + plot_layout(widths = c(5, 3))) / guide_area() +
  plot_layout(
    heights = c(4, 1),
    guides = "collect"
  ) + 
  plot_annotation(tag_levels = "A")
)

ggsave(
  "plots/misc/220531_risk-groups-plot.pdf", plot = p_assembled,
  units = "mm", width = 210, height = 135, scale = .9
)
