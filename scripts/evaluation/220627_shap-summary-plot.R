library(tidyverse)
library(dbplyr)
library(patchwork)
library(ggbeeswarm)

conn <- DBI::dbConnect(odbc::odbc(), "transdb", database = "pmhnet")

theme_set(theme_minimal())

standardize <- function(x) {
  ((x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE)))
}

###############################################################################
# prepare shap data

shap <- tbl(conn, in_schema("pmhnet", "shap_006")) |> 
  mutate(
    fgroup = sql("substring(feature from '^[^_]+')"),
    fname  = sql("substring(feature from '^[^_]+_([^_]+)')"),
    fvalue = sql("substring(feature from '^[^_]+_[^_]+_([^_]+)')"),
    value  = vals
  ) %>% 
  filter(timepoint == 1825L) |> 
  select(pid, fgroup, fname, fvalue, value, shap)

# aggregate one-hot encoded features
shap_onehot <- shap |> 
  group_by(pid, fgroup, fname) |>
  filter(n() > 1L) |> 
  mutate(index = row_number(fvalue)) |> 
  summarise(
    text  = str_flatten(fvalue[value == 1], collapse = ";") |> coalesce("missing"),
    value = sum(index[value == 1]) |> coalesce(0),
    shap  = sum(shap)
  ) |> 
  ungroup()

shap <- shap |> 
  anti_join(shap_onehot, by = c("pid", "fgroup", "fname")) |> 
  select(-fvalue) |> union_all(shap_onehot) |> 
  compute()

# feature labels for pretty printing
labels <- read_tsv("data/resources/code_definitions.tsv")

###############################################################################

average_impact <- shap |> 
  group_by(fgroup, fname) |> 
  summarise(ave = mean(abs(shap))) |> 
  collect() |> ungroup() |> 
  slice_max(order_by = ave, n = 25) 

assign_levels <- function (.data, ...) {
  if (n_distinct(.data$text) %in% 3:9) {
    .data$level <- .data$text |> as_factor() |> fct_relevel(sort)
    suppressWarnings(
      .data$level <- .data$level |> 
        fct_relevel("missing", "nan", after = Inf) |> 
        fct_relevel("[0,10]", after = 0L)
    ) 
    .data$level <- as.integer(.data$level)
  } else {
    .data$level <- 1L 
  }
  .data
}

local_impact <- shap |> 
  semi_join(average_impact, copy = TRUE) |> 
  collect() |> 
  left_join(labels, by = "fname") |> 
  group_by(definition) |> 
  group_map(assign_levels, .keep = TRUE) |> 
  bind_rows()


labels <- local_impact |> 
  distinct(fname, definition, text, level) |> 
  group_by(fname, definition) |> 
  arrange(level, .by_group = TRUE) |> 
  summarise(levels = str_flatten(text, collapse = "/"), .groups = "drop") |> 
  mutate(
    label = case_when(
      levels |> is.na() ~ definition,
      TRUE              ~ as.character(glue::glue("{definition}\n({levels})"))
    )
  ) |> 
  select(fname, label)
  
###############################################################################
# export the data for recreation locally

average_impact |> 
  left_join(labels, by = "fname") |> 
  mutate(label = fct_reorder(label, ave)) |> 
  write_tsv("~/export/230104_shap-top-25.tsv")

local_impact |> 
  left_join(labels, by = "fname") |> 
  mutate(
    label = fct_reorder(label, shap, .fun = function(.x) mean(abs(.x)))
  ) |> 
  group_by(label) |> 
  mutate(
    value = if (n_distinct(level) > 1) level  else  value,
    value = if_else(text %in% c("missing", "nan"), NaN, value),
    value = standardize(value),
  ) |> 
  ungroup() |> 
  select(
    fgroup, fname, ftext = text, flevel = level, normalized_value = value, 
    shap_value = shap, description = definition, label = label
  ) |> 
  write_tsv("~/export/230104_shap-top-25-values.tsv")

###############################################################################


p1 <- average_impact |> 
  left_join(labels, by = "fname") |> 
  mutate(label = fct_reorder(label, ave)) |> 
  ggplot(aes(x = ave, y = label)) +
  geom_col(fill = "#3b82e2") +
  scale_x_continuous(labels = scales::label_percent()) +
  labs(y = "", x = "mean(|SHAP value|)")

p2 <- local_impact |> 
  left_join(labels, by = "fname") |> 
  mutate(
    label = fct_reorder(label, shap, .fun = function(.x) mean(abs(.x)))
  ) |> 
  group_by(label) |> 
  mutate(
    value = if (n_distinct(level) > 1) level  else  value,
    value = if_else(text %in% c("missing", "nan"), NaN, value),
    value = standardize(value)
  ) |>
  ggplot(aes(x = shap, y = label, color = value, group = -level)) +
  geom_quasirandom(
    groupOnX = FALSE, size = .8,  
    dodge.width = .6, varwidth = TRUE, width = .3, bandwidth = 0.5
  ) +
  scale_color_gradient(
    low = "#3b82e2",  high = "#fc2a5a", na.value = "lightgrey"
  ) +
  scale_x_continuous(labels = scales::label_percent()) +
  labs(y = "", x = "SHAP value") +
  guides(
    color = guide_colorbar(
      label = FALSE, barwidth = 0.5, barheight = 20,
      title.position = "right", title = "feature\nvalue", title.hjust = 0.5
    )
  )

p_assembled <- p1 + (p2 * theme(axis.text.y = element_blank())) +
   plot_layout(widths = c(1, 3))

ggsave("plots/shap/220627_shap-summary-combined.pdf", plot = p_assembled, 
       unit = "mm", width = 210, height = 180, scale = 1.5)

###############################################################################
