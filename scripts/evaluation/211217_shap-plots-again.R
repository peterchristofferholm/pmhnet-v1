library(tidyverse)
library(patchwork)
library(vctrs)

con <- DBI::dbConnect(
    odbc::odbc(), driver = "PostgreSQL", server = "trans-db-01", 
    database = "pmhnet", uid = Sys.getenv("DBUSER"),
    pwd = Sys.getenv("DBPASS"), bigint = "integer"
)

DBI:::dbExecute(con, "SET search_path TO pmhnet;")


## global shap importance ######################################################

global_shap <- tbl(con, "shap_006") %>% 
  mutate(
    feature = sql("substring(feature from '^[^_]+_[^_]+')")
  ) %>% 
  group_by(timepoint, feature, pid) %>%
  summarise(shap = sum(shap)) %>%       # per tau, feature, patient
  summarise(shap = sum(abs(shap))) %>%  # per tau, feature
  mutate(shap = shap / sum(shap)) %>%   # per tau
  mutate(
    rank = dense_rank(desc(shap)),
    fgroup = sql("substring(feature from '^[^_]+')"),
    fname  = sql("substring(feature from '^[^_]+_([^_]+)')")
  ) %>%
  mutate(
    fgroup = if_else(fgroup == "diag-simple", "diag-1", fgroup)
  ) %>% 
  select(-feature) %>% 
  ungroup() %>% collect() %>%
  mutate(
    timepoint = timepoint %>%
      recode("183" = "6m", "365" = "1y", "1095" = "3y", "1825" = "5y") %>%
      fct_inorder(),
    fgroup_label = fgroup %>% fct_recode(
      "Biochemical" = "bioc-1",
      "Diagnoses" = "diag-1",
      "Clinical Characteristics 1" = "clnc-1",
      "Clinical Characteristics 2" = "clnc-2",
      "Procedures" = "proc-1"
    )
  )

## p1 ##########################################################################
## stacked barplot showing the relative importance of feature groups

# color mappings
colors <- c(
  "bioc-1" = "#ebb391",
  "diag-1" = "#7ba79d",
  "clnc-1" = "#e8d19d",
  "clnc-2" = "#bfd3e4",
  "proc-1" = "#cadcd8"
)

theme_set(theme_void())
theme_update(
  legend.position = "none",
  axis.title.y = element_text(
    angle = 90, margin = margin(r = 5), hjust = .06, vjust = 0
  )
)

p1 <- global_shap %>% 
  filter(timepoint == "5y") %>% 
  group_by(timepoint, fgroup_label, fgroup) %>% 
  summarise(shap = sum(shap)) %>% 
  mutate(fgroup = fct_reorder(fgroup, -shap)) %>% 
  ggplot(aes(x = timepoint, y = shap, fill = fgroup)) +
  geom_col(color = "black") +
  geom_text(
    aes(label = str_wrap(fgroup_label, 20)), 
    position = position_stack(vjust = 0.5)
  ) +
  labs(y = "Relative feature importance (SHAP)") +
  coord_fixed(ratio = 3.5) +
  scale_fill_manual(values = colors)


# we want text-labels for the top 25 features
labels <- global_shap %>% 
  filter(timepoint == "5y", rank <= 25) %>% 
  transmute(
    label = fname, fname, xpos = 2,
    ypos = 1.01 - seq(0, .95, length.out = n() + 2) %>% tail(-1) %>% head(-1)
  )

p2 <- global_shap %>% 
  filter(timepoint == "5y") %>% 
  mutate(
    xmin = 0, xmax = 1, 
    ymin = 1 - head(cumsum(c(0, shap)), n()),
    ymax = 1 - tail(cumsum(c(0, shap)), n())
  ) %>% 
  left_join(labels, by = "fname") %>% 
  ggplot() +
  geom_rect(
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = fgroup), 
    color = "black"
  ) +
  geom_text(
    data = ~filter(., !is.na(label)),
    aes(y = ypos, x = xpos, label = label), hjust = 0
  ) +
  geom_segment(
    data = ~filter(., !is.na(label)),
    aes(y = (ymin + ymax) / 2, x = 1, yend = ypos, xend = xpos)
  ) +
  expand_limits(x = 3) +
  theme(axis.title.y = element_blank()) +
  coord_fixed(ratio = 10, clip = "off") +
  scale_fill_manual(values = colors)

ggsave("plots/shap/s006-feature-importance-overview.pdf", plot = (p1 + p2),
       height = 8, width = 6, unit = "cm", scale = 2)

################################################################################
## exploring different individual features

# update theme
theme_set(theme_light())
theme_update(
  plot.margin = unit(c(.5, .5, .5, .5), "cm"),
  axis.title.y = element_text(
    angle = 90, hjust = 0, vjust = 0
  ),
  axis.title.x = element_text(
    angle = 0, hjust = 0, vjust = 0
  ),
  legend.background = element_rect(fill = "transparent"),
)

shap <- tbl(con, "shap_006") %>% 
  mutate(
    fgroup  = sql("substring(feature from '^[^_]+')"),
    fname   = sql("substring(feature from '^[^_]+_([^_]+)')")
  ) %>% 
  mutate(
    fgroup = if_else(fgroup == "diag-simple", "diag-1", fgroup)
  )

## effect of age 

shap_age <- shap %>% 
  filter(fname == "age") %>% 
  collect() %>% 
  mutate(
    timepoint = timepoint %>%
      recode("183" = "6m", "365" = "1y", "1095" = "3y", "1825" = "5y") %>%
      fct_inorder()
  )

p_age <- shap_age %>% 
  filter(timepoint == "5y") %>% 
  ggplot(aes(x = vals, y = shap, fill = shap)) +
  geom_point(size = 2, alpha = 1, shape = 21) +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_fill_gradient2(
    breaks = c(-0.2, 0, 0.2), 
    labels = c("non-survival", "unchanged", "survival")
  ) +
  labs(
    title = "Effect of age on model prediction",
    x = "Age at time zero (Z-score normalised)",
    y = "Impact on model output (SHAP value)",
    fill = ""
  ) +
  theme(
    legend.position = c(0.99, 0.99), legend.justification = c(1, 1)
  )

ggsave("plots/shap/s006-effect-of-age.pdf", plot = p_age, 
       width = 8, height = 7, unit = "cm", scale = 2)  


# effect of vessels

shap %>% 
  filter(fname == "vessels", timepoint == 1825) %>% 
  mutate(
    value  = sql("substring(feature from '[^_]+$')"),
  ) %>% 
  collect() %>% 
  group_by(pid) %>% 
  summarise(
    value = value[vals == 1],
    shap = sum(shap)
  ) -> shap_vessels

shap_vessels %>% 
  mutate(
    value = fct_relevel(value, "DIF", after = 0L)
  ) %>% 
  ggplot(aes(value, shap, fill = shap)) + 
  ggbeeswarm::geom_quasirandom(method = "pseudorandom",
    size = 2, shape = 21
  ) +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_fill_gradient2(
    breaks = ~ c(.x[1], 0, .x[2]) * 0.9,
    labels = c("non-survival", "unchanged", "survival")
  ) +
  labs(
    title = "Effect of vascular pathology on model prediction",
    x = "Number of vessels found occluded at index CAG",
    y = "Impact on model output (SHAP value)",
    fill = ""
  ) +
  theme(
    legend.position = c(0.1, 0.1), legend.justification = c(0, 0)
  ) -> p_vessels


ggsave("plots/shap/s006-effect-of-vessels.pdf", plot = p_vessels, 
       width = 8, height = 7, unit = "cm", scale = 2)  

# effect of smoking

shap %>% 
  filter(fname == "smoking", timepoint == 1825) %>% 
  mutate(
    value  = sql("substring(feature from '[^_]+$')")
  ) %>% 
  collect() %>% 
  group_by(pid) %>% 
  summarise(
    value = value[vals == 1], shap = sum(shap)
  ) -> shap_smoking

shap_smoking %>% 
  mutate(
    value = fct_relevel(value, "nan", after = 0L)
  ) %>% 
  ggplot(aes(value, shap, fill = shap)) + 
  ggbeeswarm::geom_quasirandom(method = "pseudorandom",
    size = 2, shape = 21
  ) +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_fill_gradient2(
    breaks = ~ c(.x[1], 0, .x[2]) * 0.9,
    labels = c("non-survival", "unchanged", "survival")
  ) +
  labs(
    title = "Effect of smoking on model prediction",
    x = "Smoking status reported at time zero",
    y = "Impact on model output (SHAP value)",
    fill = ""
  ) +
  theme(
    legend.position = c(0.1, 0.99), legend.justification = c(0, 1)
  ) -> p_smoking

ggsave("plots/shap/s006-effect-of-smoking.pdf", plot = p_smoking, 
       width = 8, height = 7, unit = "cm", scale = 2)


# effect of biochemical values (explore if missingness matters?)

shap %>% 
  filter(fgroup == "bioc-1", timepoint == 1825) %>% 
  mutate(
    value  = sql("substring(feature from '[^_]+$')")
  ) %>% 
  collect() %>% 
  group_by(pid, fname) %>% 
  summarise(
    shap = sum(shap),
    value = value[vals == 1] %0% "missing", 
    .groups = "drop"
  ) -> shap_bioc 

# get mappings between npu code and the component names
tbl(con, dbplyr::ident_q("core.biochem")) %>% 
  distinct(code = quantity_id, analyte = component) %>% 
  collect() -> analytes


plot_biochem_shap <- function(x) {
  shap_bioc %>% 
    filter(value == x) %>% 
    group_by(fname) %>% 
    mutate(
      median_shap = median(shap)
    ) %>% 
    ungroup() %>% 
    left_join(analytes, by = c("fname" = "code")) %>% 
    mutate(
      test = glue::glue("{analyte} ({fname})"),
      test = fct_reorder(test, shap)
    ) %>% 
    ggplot(aes(shap, test, fill = median_shap)) +
    geom_boxplot(orientation = "y", outlier.size = 1, outlier.alpha = .5) +
    labs(
      title = glue::glue("__{x}__ biochemical tests"),
      x = "Impact on model output (SHAP value)",
      y = "", fill = ""
    ) +
    scale_x_continuous(labels = scales::percent_format()) +
    scale_fill_gradient2() +
    theme(
      plot.title = ggtext::element_markdown(),
      legend.position = ""
    )
}

plots <- map(c("missing", "normal", "high"), plot_biochem_shap) 


ggsave("plots/shap/s006-effect-of-biochemical.pdf", 
       plot = wrap_plots(plots, nrow = 1, clip = FALSE), 
       width = 21, height = 12, unit = "cm", scale = 2.5)


shap_bioc %>% 
  filter(value  != "low") %>% 
  group_by(fname, value) %>% 
  summarise(shap = sum(shap)) %>% 
  ggplot(aes(fname, value, fill = shap)) +
  geom_tile() +
  scale_fill_gradient2()


set.seed(2022)
df <- shap_bioc %>% 
  filter(value  != "low") %>% 
  group_by(fname, value) %>% 
  summarise(shap = median(shap), .groups = "drop") %>% 
  pivot_wider(names_from = value, values_from = shap) %>% 
  mutate(across(-fname, ~ .x %|% 0)) 

df <- df %>% 
  select(-fname) %>% 
  kmeans(6, iter.max = 1000) %>% 
  broom::augment(df) %>%
  left_join(analytes, by = c("fname" = "code")) %>% 
  mutate(
    fname = glue::glue("{analyte} ({fname})"),
  ) %>% 
  pivot_longer(cols = high:normal) %>% 
  mutate(
    name = fct_relevel(name, "missing", after = 0L)
  ) 

p_bioc_kmeans <- df %>% 
  ggplot(aes(name, fname, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(limits = ~ c(.x[1], .x[2]) * .9, oob = scales::squish) +
  facet_grid(.cluster ~ ., scales = "free_y", space = "free_y") +
  labs(
    title = "Average model impact across different biochemical tests",
    subtitle = "Panels: k-means clusters (k=6)",
    x = "Results of biochemical tests", fill = "mean(SHAP)", y = ""
  ) +
  theme(
    plot.title.position = "panel",
    plot.title = element_text(hjust = 1),
    plot.subtitle = element_text(hjust = 1)
  )

ggsave("plots/shap/s006-shap-bioc-kmeans.pdf", 
       plot = p_bioc_kmeans,
       width = 8, height = 11, unit = "cm", scale = 2.7)


## shap summary plot ###########################################################
# top15 most impactful features as a beeswarm plot showing feature effects
################################################################################

# are features categorical or continuous?
ftypes <- shap %>% 
  distinct(fname, vals) %>% 
  group_by(fname) %>% 
  summarise(
    ftype = if_else(n() == 2L, "categorical", "continuous"), .groups = "drop"
  ) %>% 
  compute()

# top15 features in terms of global attribution
top15 <- shap %>% 
  filter(timepoint == 1825L) %>% 
  group_by(fname, pid) %>% 
  summarise(shap = sum(shap)) %>% 
  summarise(attribution = mean(abs(shap))) %>% 
  slice_max(attribution, n = 20, with_ties = FALSE) %>% 
  compute()

# collect shap values for top15 features
top_shap <- shap %>% 
  filter(timepoint == 1825L) %>% 
  inner_join(top15, by = "fname") %>% 
  inner_join(ftypes, by = "fname") %>% 
  mutate(fval = sql("substring(feature from '[^_]+$')")) %>% 
  collect() %>% 
  left_join(analytes, by = c("fname" = "code")) %>% 
  mutate(
    panel = case_when(
      !is.na(analyte) ~ paste0(analyte, " (", fname, ")"),
      fname == "DJ44" ~ "COPD DIAGNOSIS (DJ44)",
      fname == "enzymes" ~ "CARDIAC ENZYMES STATUS",
      fname == "killip" ~ "KILLIP CLASS",
      fname == "lvef" ~ "LEFT VENTRICULAR EJECTION FRACTION (LVEF)",
      fname == "nyha" ~ "NYHA CLASS",
      TRUE ~ str_to_upper(fname) 
    ) %>% str_wrap(width = 20)
  )

## plot of continuous features

p1 <- top_shap %>% 
  filter(ftype == "continuous") %>% 
  group_by(panel) %>% 
  mutate(
    value = scales::rescale(vals, c(-1, 1)),
    value = if_else(vals == 0, NaN, value)
  ) %>% 
  ungroup() %>% 
  ggplot(aes(y = "value", x = shap, color = value)) +
  ggbeeswarm::geom_quasirandom(groupOnX = FALSE, size = 1) +
  labs(x = "", y = "", title = "Continuous features", fill = "Feature value") +
  facet_grid(panel ~ ., switch = "y")  +
  scale_color_gradient(low = "blue", high = "red") +
  scale_x_continuous(limits = c(-0.4, 0.3)) +
  theme(legend.position = "bottom")

# plot of categorical features

top_shap_cate <- top_shap %>% 
  filter(ftype == "categorical") %>% 
  mutate(panel = fct_reorder(panel, attribution, .desc = TRUE)) %>% 
  group_by(panel, fname, pid) %>% 
  summarise(
    shap = sum(shap),
    value = fval[vals == 1] %0% "missing"
  ) %>%
  ungroup() %>% 
  mutate(
    value = value %>% 
      recode("nan" = "missing", "True" = "elevated", "False" = "normal") %>% 
      fct_relevel("missing", after = Inf) %>% 
      fct_relevel("[0,10]", after = 0L)
  ) %>% 
  group_by(fname) %>% 
  mutate(
    level = fct_drop(value) %>% as.integer() %>%  as.character()
  ) %>% 
  ungroup()

p2 <- top_shap_cate %>% 
  ggplot(aes(y = value, x = shap, color = level)) +
  ggbeeswarm::geom_quasirandom(groupOnX = FALSE, size = 1) +
  facet_grid(panel ~ ., scales = "free_y", space = "free_y", switch = "y") +
  scale_x_continuous(limits = c(-0.4, 0.3)) +
  labs(
    x = "Impact on model output (SHAP value)", y = "", 
    title = "Categorical features"
  ) +
  theme(legend.position = "none")

p_top15 <- (p1 + p2 +
  plot_layout(heights = c(1, 12), ncol = 1) &
  theme(
    strip.placement = "outside", 
    strip.background.y = element_rect(fill = "white"),
    strip.text.y.left = element_text(
      angle = 0, color = "black", hjust = 1, size = 7
    ),
    plot.margin = margin(),
    plot.title = element_text(size = 10),
    axis.text.y = element_text(size = 5),
    panel.spacing = unit(1/10, "lines"),
    legend.position = "none"
  )) +
  plot_annotation(
    theme = theme(plot.margin = margin(20, 20, 20, 20))
  )

ggsave("plots/shap/s006-shap-top15-summary.pdf", plot = p_top15,
       width = 12, height = 14, unit = "cm", scale = 2)