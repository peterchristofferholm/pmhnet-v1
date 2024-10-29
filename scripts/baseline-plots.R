library(tidyverse)
library(survival)

in_schema <- dbplyr::in_schema

theme_set(theme_minimal())

con <- DBI::dbConnect(odbc::odbc(),
    driver = "PostgreSQL", server = "dbserver", database = "database",
    uid = Sys.getenv("DBUSER"), pwd = Sys.getenv("DBPASS"), bigint = "integer"
)
surv_df <- tbl(con, in_schema("pmhnet", "surv_master")) %>% collect()
sfit <- survfit(Surv(time, event == "1") ~ train, data = surv_df)


# kaplan meier curve -----------------------------------------------------------
p1 <- broom::tidy(sfit, scale = 365.25) %>%
  filter(time < 365 * 5) %>%
  mutate(
    time = time / 365.25,
    split = recode(strata, "train=0" = "test", "train=1" = "train")
  ) %>%
  ggplot(aes(time, estimate)) +
  geom_ribbon(
    aes(ymin = conf.low, ymax = conf.high, fill = split), alpha = .25
  ) +
  geom_step(aes(color = split)) +
  geom_text(
    aes(label = split, x = x, y = y, color = split),
    ~ group_by(., split) %>%
      mutate(y = if_else(split != "train", conf.high, conf.low)) %>%
      summarise(x = 4.9, y = approx(time, y, xout = 4.9)$y),
    angle = -25
  ) +
  labs(
    title = "KM-estimates of all-cause mortality",
    x = "Follow-up time (years)", y = "Survival"
  ) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_x_continuous(minor_breaks = scales::extended_breaks(n = 4 * 6)) +
  theme(legend.position = "none")

ggsave(
  snakemake@output[["p01"]], plot=p1, width=12, height=7, unit="cm", scale=1.5
)

## various statistics that might be of relevance -------------------------------

log <- file(snakemake@output[["txt"]], "wt")  # redirect stdout
sink(log)

writeLines(c("", "Restricted mean survival time:"))
print(sfit, rmean = 365*5)

writeLines(c("", "Test of difference between splits (log-rank test):"))
surv_df %>%
  mutate(
    event = if_else(time > 365.25 * 5, "0", event),
    time = if_else(time > 365.25 * 5, 365.25 * 5, as.numeric(time))
  ) %>%
  survdiff(Surv(time, event == "1") ~ train, data = .)
