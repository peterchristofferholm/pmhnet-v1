library(tidyverse)
library(dbplyr)
library(patchwork)

theme_set(theme_minimal())

con <- DBI::dbConnect(odbc::odbc(),
    driver = "PostgreSQL", server = "dbserver", database = "database",
    uid = Sys.getenv("DBUSER"), pwd = Sys.getenv("DBPASS")
)

sid <- snakemake@wildcards[["sid"]]
schema <- snakemake@config[["dbschema_optuna"]]

studies <- tbl(con, in_schema(schema, "studies"))
trials  <- tbl(con, in_schema(schema, "trials"))
values  <- tbl(con, in_schema(schema, "trial_values"))
params  <- tbl(con, in_schema(schema, "trial_params"))


## OPTIMIZATION HISTORY ########################################################

p1 <- studies %>%
  inner_join(trials, by = "study_id") %>%
  inner_join(values, by = "trial_id") %>%
  filter(value < 1.0, study_name == sid) %>%
  arrange(number) %>%
  mutate(
    best = value == cummin(value)
  ) %>%
  collect() %>%
  ggplot(aes(number, value)) +
  geom_point(alpha = .5, stroke = 0) +
  geom_line(data = ~filter(., best == "1"), color = "red") +
  geom_point(data = ~filter(., best == "1"), color = "red") +
  labs(
    x = "#Trial", y = "Negative Log-Likelihood (5-fold CV)",
    title = glue::glue("Study {sid}: Hyperparameter Optimization")
  )

ggsave(
  plot = p1, filename = snakemake@output$p1,
  width = 12, height = 6, unit = "cm", scale = 1.2
)

## SLICE PLOTS #################################################################

p2 <- studies %>%
  inner_join(trials, by = "study_id") %>%
  inner_join(values, by = "trial_id") %>%
  inner_join(params, by = "trial_id") %>%
  filter(value < 1.0, study_name == sid) %>%
  collect() %>%
  group_by(param_name) %>%
  summarise(
    plot = list(
      ggplot(cur_data(), aes(param_value, value, color = number)) +
        geom_jitter(stroke = 0, size = 1) +
        geom_vline(
          aes(xintercept = param_value),
          data = ~filter(., value == min(value)),
          linetype = "dashed", alpha = .5
        ) +
        guides(color = FALSE) +
        labs(title = cur_group()) +
        theme(
          plot.title = element_text(size = 8, hjust = 0.5),
          axis.title = element_blank()
        )
    )
  ) %>%
  pull(plot) %>%
  wrap_plots(ncol = 4)

ggsave(
  plot = p2, filename = snakemake@output$p2,
  height = (length(p2$patches$plots) %/% p2$patches$layout$ncol) * 4,
  width = 19, unit = "cm", scale = 1.4
)
