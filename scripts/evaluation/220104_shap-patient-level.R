library(tidyverse)
source("scripts/misc/ggenearrow.R")  # library(ggenearrow)

con <- DBI::dbConnect(
    odbc::odbc(), driver = "PostgreSQL", server = "trans-db-01", 
    database = "pmhnet", uid = Sys.getenv("DBUSER"),
    pwd = Sys.getenv("DBPASS"), bigint = "integer"
)

DBI::dbExecute(con, "SET search_path TO pmhnet;")

shap <- tbl(con, "shap_006") %>% 
  mutate(
    fgroup  = sql("substring(feature from '^[^_]+')"),
    fname   = sql("substring(feature from '^[^_]+_([^_]+)')"),
    fvalue  = sql("substring(feature from '^[^_]+_[^_]+_([^_]+)')"),
    value  = vals
  ) %>% 
  filter(timepoint == 1825L) %>% 
  select(pid, fgroup, fname, fvalue, value, shap)

clean <- shap %>% 
  mutate(value = round(value, 1)) %>% 
  group_by(pid, fgroup, fname) %>% 
  mutate(
    fvalue = if_else(
      n() > 1 & all(value == 0), "normal", fvalue
    )
  ) %>% 
  mutate(
    shap = sum(shap), 
    here = case_when(
      all(fvalue == "normal")  ~ row_number() == 1L,
      n() > 1                  ~ value != 0,
      TRUE                     ~ TRUE
    )
  ) %>% 
  filter(here) %>% 
  select(-here) %>% 
  compute()

# remove all but top most impactful features
top <- clean %>% 
  group_by(pid) %>% 
  filter(row_number(-abs(shap)) <= 8)

## find patients with non-defining top features (n >= 5) #######################

# feature vectors for all patients in cohort
data <- read_tsv("data/scratch/006-input.tsv") %>% 
  mutate(across(-pid, ~round(.x, 1)))

# how many patients matches the fingerprint of the top features?
lookup <- top %>% 
  transmute(
    pid,
    names  = paste(fgroup, fname, fvalue, sep = "_"),
    values = value 
  ) %>% 
  collect() %>% 
  summarise(
    n_matches = sum(map2(names, values, ~data[[.x]] == .y) %>% pmap_lgl(all))
  )

# prepare data for the selected ones
candidates <- lookup %>% 
  filter(n_matches >= 5) %>%  
  inner_join(clean, by = "pid", copy = TRUE) %>% 
  mutate(
    value = case_when(
      fname == "age" ~ (value * 90) + 15,
      TRUE           ~ value
    )
  ) %>% 
  mutate(
    feature = paste0(fname, ": ", coalesce(fvalue, as.character(value)))
  ) %>% 
  select(pid, feature, shap) %>% 
  group_by(pid) %>% 
  mutate(
    feature = if_else(row_number(-abs(shap)) > 8, "rest", feature)
  ) %>% 
  group_by(pid, feature) %>% 
  summarise(shap = sum(shap)) %>% 
  ungroup() %>% 
  collect()

################################################################################

set.seed(666)

# select one case per prediction category
cases <- candidates %>% 
  group_by(pid) %>% 
  summarise(yhat = 0.821 + sum(shap)) %>% 
  filter(yhat <= 1.00) %>% 
  mutate(category = cut(yhat, breaks = 0:4 / 4)) %>% 
  group_by(category) %>% 
  slice_sample(n = 1) %>% 
  inner_join(candidates, by = "pid")

theme_set(theme_minimal())

# construct waterfall plot
p_shap <- cases %>% 
  arrange(-abs(shap)) %>% 
  mutate(
    dir = if_else(shap < 0, "neg", "pos"),
    rank = row_number(abs(shap)),
    xmax = 0.821 + cumsum(shap),
    xmin = lag(xmax, default = 0.821),
    category = glue::glue("Patient with prediction in {category}")
  ) %>% 
  ggplot(aes(y = rank)) +
  geom_vline(xintercept = 0.821, linetype = "dashed", alpha = .25, size = 1/2) +
  geom_vline(
    data = ~ distinct(.x, yhat, category), 
    mapping = aes(xintercept = yhat), alpha = .25, size = 1/2
  ) +
  geom_text(
    data = ~ filter(.x, dir == "neg"),
    mapping = aes(label = feature, x = xmin), hjust = 0, nudge_x = 1/100,
    alpha = 1, size = 3
  ) +
  geom_text(
    data = ~ filter(.x, dir == "pos"),
    mapping = aes(label = feature, x = xmin), hjust = 1, nudge_x = -1/100,
    alpha = 1, size = 3
  ) +
  geom_gene_arrow(
    aes(xmin = xmin, xmax = xmax, fill = dir),
    arrowhead_height = unit(5, "mm"), 
    arrow_body_height = unit(5, "mm"), 
    arrowhead_width  = unit(1, "mm"),
  ) +
  facet_wrap(~category) +
  coord_cartesian(clip = "off") +
  scale_x_continuous(labels = scales::percent_format(), limits = c(0, 1)) +
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

ggsave(
  "plots/shap/s006-personalized.pdf", plot = p_shap,
  width = 15, height = 6, unit = "cm", scale = 1.5
)
