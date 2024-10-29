library(tidyverse)
source("scripts/misc/ggenearrow.R")  # library(ggenearrow)

in_schema <- dbplyr::in_schema

theme_set(theme_minimal())

conn <- DBI::dbConnect(odbc::odbc(), "transdb", database = "pmhnet")

###############################################################################

shap <- tbl(conn, in_schema("pmhnet", "shap_006")) |> 
  mutate(
    fgroup = sql("substring(feature from '^[^_]+')"),
    fname  = sql("substring(feature from '^[^_]+_([^_]+)')"),
    fvalue = sql("substring(feature from '^[^_]+_[^_]+_([^_]+)')"),
    value  = vals
  ) %>% 
  filter(timepoint == 1825L) |> 
  select(pid, fgroup, fname, fvalue, value, shap)

shap <- shap |> 
  group_by(pid, fgroup, fname) |> 
  mutate(
    fvalue = if_else(n() > 1 & all(value == 0, na.rm = TRUE), "missing", fvalue),
    shap = sum(shap)
  ) |> 
  mutate(
    keep = case_when(
      all(fvalue == "missing", na.rm = TRUE) ~ row_number() == 1L,
      n() > 1 ~ value != 0, TRUE ~ TRUE
    )
  ) |> 
  ungroup() |>  filter(keep) |>  collect() |> 
  mutate(pid = as.integer(pid))

###############################################################################

cases <- shap |>  nest(data = -pid) |>  slice_sample(n = 6)

# waterfall plot ##############################################################
aggregate_attributions_1 <- function(data, .sort = 1, .n = 10, .base = 0.821) {
  order <- function(x) { .sort |> switch (-abs(x), abs(x), x, -x ) }
  data  <- data |> 
    mutate(
      value = case_when(
        fname == "age" ~ (value * 90) + 15,
        TRUE           ~ value
      ),
      fvalue = coalesce(fvalue, as.character(round(value, 1))),
      label  = glue::glue("{fname}: {fvalue}") |> as.character(),
      label  = if_else(row_number(-abs(shap)) > .n, "rest", label)
    ) |> 
    group_by(label) |> 
    summarise(shap = sum(shap))
 
  data |>  
    arrange(order(shap)) |> 
    mutate(
      dir = if_else(shap < 0, "neg", "pos"),
      rank = row_number(-order(shap)),
      xmax = .base + cumsum(shap),
      xmin = lag(xmax, default = .base),
      yhat = .base + sum(shap)  # f(x)
    )
}

create_waterfall_plot <- function(data, .base = 0.821) {
  ggplot(data, aes(y = rank)) +
    geom_vline(xintercept = .base, linetype = "dashed", alpha = .25, size = 1/2) +
    geom_vline(
      data = ~ distinct(.x, yhat, pid), 
      mapping = aes(xintercept = yhat), alpha = .25, size = 1/2
    ) +
    geom_text(
      data = ~ filter(.x, dir == "neg"),
      mapping = aes(label = label, x = xmin), hjust = 0, nudge_x = 1/100,
      alpha = 1, size = 2.5
    ) +
    geom_text(
      data = ~ filter(.x, dir == "pos"),
      mapping = aes(label = label, x = xmin), hjust = 1, nudge_x = -1/100,
      alpha = 1, size = 2.5
    ) +
    geom_gene_arrow(
      aes(xmin = xmin, xmax = xmax, fill = dir),
      arrowhead_height = unit(5, "mm"), 
      arrow_body_height = unit(5, "mm"), 
      arrowhead_width  = unit(1, "mm"),
    ) +
    facet_wrap(~pid) +
    coord_cartesian(clip = "off") +
    scale_x_continuous(labels = scales::percent_format()) +
    labs(
      x = "Model prediction (5-year survival)", y = "",
      title = "Patient-level feature attributions (SHAP)"
    ) +
    theme(
      legend.position = "none", 
      axis.text.y = element_blank(), axis.ticks.y = element_blank(),
      plot.margin = margin(0.5, 1.5, 0.5, 1.5, "cm"),
      panel.spacing = unit(0.75, "cm"),
      strip.text = element_text(hjust = 0, face = "bold")
    )
}

cases |> 
  mutate(
    data = map(data, aggregate_attributions_1, .sort = 1, .n = 10),
    pid  = glue::glue("case {xfun::n2w(row_number())}") |> as_factor()
  ) |> 
  unnest(cols = data) |> 
  create_waterfall_plot()

ggsave("plots/shap/220622_waterfall-plot-1.pdf", unit = "mm", width = 210,
       height = 160, scale = 1.1)

cases |> 
  mutate(
    data = map(data, aggregate_attributions_1, .sort = 3, .n = 10),
    pid  = glue::glue("case {xfun::n2w(row_number())}") |> as_factor()
  ) |> 
  unnest(cols = data) |> 
  create_waterfall_plot()

ggsave("plots/shap/220622_waterfall-plot-2.pdf", unit = "mm", width = 210,
       height = 160, scale = 1.1)

# force plot ##################################################################

aggregate_attributions_2 <- function(data, .n = 10, .base = 0.81) {
  data  <- data |> 
    mutate(
      value = case_when(
        fname == "age" ~ (value * 90) + 15,
        TRUE           ~ value
      ),
      fvalue = coalesce(fvalue, as.character(round(value, 1))),
      label  = glue::glue("{fname}: {fvalue}") |> as.character(),
      label  = if_else(row_number(-abs(shap)) > .n, "rest", label)
    ) |> 
    group_by(label) |> 
    summarise(
      shap = sum(shap),
      dir = if_else(shap < 0, "neg", "pos")
    )
  
  yhat   <- .base + sum(data$shap)
  offset <- yhat - sum(filter(data, dir == "pos")$shap)
  
  data |> 
    arrange(desc(dir), shap) |> 
    mutate(
      xmax = cumsum(abs(shap)),
      xmin = lag(xmax, default = 0),
      xmax = xmax + offset, 
      xmin = xmin + offset, 
      .xmax = if_else(dir == "neg", xmin, xmax),
      .xmin = if_else(dir == "neg", xmax, xmin),
      xmax = .xmax, xmin = .xmin, yhat = yhat
    ) |> 
    select(-starts_with("."))
}

cases |> 
  mutate(
    data = map(data, aggregate_attributions_2, .n = 10),
    pid  = glue::glue("case {xfun::n2w(row_number())}") |> as_factor()
  ) |> 
  unnest(c(data)) |> 
  ggplot(aes(y = pid)) +
  geom_gene_arrow(
    aes(xmin = xmin, xmax = xmax, fill = dir),
    arrowhead_height = unit(5, "mm"), 
    arrow_body_height = unit(5, "mm"), 
    arrowhead_width  = unit(1, "mm"),
  ) +
  geom_label(
    aes(label = label, x = (xmin + xmax) / 2), hjust = .5, nudge_y = c(0.25, -0.25),
    data = ~filter(.x, abs(shap) > 0.01), size = 2.5, alpha = .5
  ) +
  coord_cartesian(clip = "off") +
  scale_x_continuous(labels = scales::percent_format()) +
  labs(
    x = "Model prediction (5-year survival)", y = "",
    title = "Patient-level feature attributions (SHAP)"
  ) +
  theme(
    legend.position = "none", 
    axis.text.y = element_blank(), axis.ticks.y = element_blank(),
    plot.margin = margin(0.5, 1.5, 0.5, 1.5, "cm"),
    panel.spacing = unit(0.75, "cm"),
    strip.text = element_text(hjust = 0, face = "bold")
  )
  
ggsave("plots/shap/220622_force-plot.pdf", unit = "mm", width = 210,
       height = 160, scale = 1.1)

