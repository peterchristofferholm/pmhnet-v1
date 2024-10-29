library(tidyverse)
library(ggbeeswarm)
library(ggrepel)
library(scales)

in_schema <- dbplyr::in_schema

theme_set(theme_minimal())

con <- DBI::dbConnect(odbc::odbc(),
    driver = "PostgreSQL", server = "dbserver", database = "database",
    uid = Sys.getenv("DBUSER"), pwd = Sys.getenv("DBPASS"), bigint = "integer"
)

shap <- tbl(con, in_schema("pmhnet", "shap_006")) %>%
  mutate(
    fgroup = sql("substring(feature from '^[^_]+')"),
    fname  = sql("substring(feature from '^[^_]+_([^_]+)')")
  )

# global feature importance ----------------------------------------------------
global_shap <- shap %>%
  group_by(timepoint, fname, pid) %>%
  summarise(shap = sum(shap)) %>%       # per tau, fname, patient
  summarise(shap = sum(abs(shap))) %>%  # per tau, fname
  mutate(shap = shap / sum(shap)) %>%   # per tau
  mutate(rank = dense_rank(desc(shap))) %>%
  ungroup() %>% collect() %>%
  mutate(
    timepoint = timepoint %>%
      recode("183" = "6m", "365" = "1y", "1095" = "3y", "1825" = "5y") %>%
      fct_inorder()
  )

(global_shap %>%
  arrange(timepoint) %>%
  mutate(
    timepoint = timepoint %>% 
      fct_recode(
        "6 months" = "6m", "1 year" = "1y", "3 years" = "3y", "5 years" = "5y"
      )
  ) %>%
  ggplot(aes(x = timepoint, y = shap)) +
  geom_quasirandom(width = 0.45, size = .5) +
  geom_text_repel(
    data = ~filter(., rank <= 10), aes(label = fname),
    direction = "both", box.padding = .2, nudge_x = .3, point.padding = .3,
    # appearance
    size = 3,            # font size
    segment.alpha = .5,  # line alpha
    segment.size = .2    # line width
    ) +
  labs(
    x = "Prediction timepoint",
    y = "Relative feature importance (SHAP)",
    title = "Top-10 highest ranking features (global)",
    fill = ""
  ) -> p1
)

ggsave(
    "results/study_006/210809/01_shap-top-global.pdf", plot=p1,
    width=14, height=9, unit="cm", scale = 1.3
)

"relative importance of **individual features**"
global_shap %>% 
  mutate(
    relative_importance = scales::percent_format(accuracy = 0.01)(shap / sum(shap))
  ) %>% 
  filter(rank <= 10) %>% 
  print(n = 40)

"relative importance of **feature groups**"
global_shap %>%
  inner_join(
    distinct(shap, fname, fgroup), by = "fname", copy = TRUE
  ) %>%
  group_by(fgroup) %>%
  summarise(shap = sum(shap)) %>%
  transmute(
    fgroup, relative_importance = scales::percent_format()(shap / sum(shap))
  )

################################################################################
## n=1 shap plots

pred <- read_csv("results/study_006/test_predictions.csv", col_types = cols()) %>%
  pivot_longer(
    -pid, names_to = "tau", names_transform = list(tau = as.numeric)
  )

tps <- shap %>% distinct(timepoint) %>% collect()

preds <- crossing(pred, tps) %>%
  group_by(pid, timepoint) %>%
  summarise(
    yhat = approx(tau, value, xout = cur_group()$timepoint)$y, .groups = "drop"
  )

set.seed(144)

"sample 4 cases with different PMHnet scores"
(cases <- preds %>%
  filter(timepoint == 1825L) %>%
  mutate(
    category = cut(yhat, breaks = 0:4 / 4)
  ) %>%
  group_by(category) %>%
  slice_sample(n = 1) %>%
  ungroup() %>%
  transmute(
    pid, category, id = row_number(desc(category)), pred = yhat
  )
)

# helper function for aggregating features into variables
create_label <- function(group, data) {
  if (nrow(data) > 1L) {
    if (all(data$x == 0)) {
      paste(group$fname, ": missing")
    } else {
      paste(group$fname, ":", data$fval[data$x])
    }
  } else {
    if (group$fname == "age") {
      paste(group$fname, ":", round(data$x * 90 + 15, 1))
    } else {
      paste(group$fname, ":", round(data$x, 1))
    }
  }
}

## create dataframe with the aggregated shap scores
df <- shap %>%
  rename(x = vals) %>%
  inner_join(cases, copy = TRUE, by = "pid") %>%
  inner_join(preds, copy = TRUE, by = c("pid", "timepoint")) %>%
  collect() %>%
  separate(
    feature, into = c("fgroup", "fname", "fval"), sep = "_", fill = "right"
  ) %>%
  group_by(category, id, timepoint, fname) %>%
  summarise(
    yhat = first(yhat), shap = sum(shap),
    dir = if_else(shap < 0, "neg", "pos"),
    labs = create_label(cur_group(), cur_data()),
    .groups = "drop_last"
  ) %>%
  group_by(dir, .add = TRUE) %>%
  arrange(-abs(shap)) %>%
  mutate(
    y1 = yhat + cumsum(shap),
    y0 = lag(y1) %>% coalesce(yhat)
  ) %>%
  ungroup() %>%
  arrange(timepoint) %>%
  mutate(
    timepoint = timepoint %>%
      recode("183" = "6m", "365" = "1y", "1095" = "3y", "1825" = "5y") %>%
      fct_inorder()
  )

(p2 <- df %>%
  mutate(id = glue::glue("case: {id}, risk score in {category}")) %>%
  ggplot(aes(x = timepoint)) +
  geom_segment(
    aes(y = y0, yend = y1, xend = timepoint, color = shap), size = 7.5
  ) +
  geom_text_repel(
    aes(label = str_wrap(labs, width = 20), y = (y0 + y1)/2, x = as.numeric(timepoint) + .1),
    ~ group_by(., id, timepoint) %>%
        filter(row_number(-abs(shap)) <= 10 | abs(shap) > 1/20) %>%
        filter(abs(shap) > 1/50),
    size = 2, direction = "y", hjust = 0, nudge_x = 1/5,
    segment.size = 1/10, box.padding = 1/10, segment.alpha = 1/2
  ) +
  scale_color_gradient2(limits = c(-1, 1)/20, oob = scales::squish) +
  scale_x_discrete(expand = expansion(add = c(1/3, 3/2))) +
  scale_y_continuous(
    labels = scales::percent, breaks = 0:5 / 5, minor_breaks = (0:4 / 5) + 1/10
  ) +
  facet_wrap(~id, nrow = 2) +
  coord_cartesian(clip = "off") +
  theme(legend.position = "none") +
  labs(
    x = "", y = "Survival prediction", 
    title = "Patient-level SHAP explanations"
  )
)

ggsave(
  "results/study_006/210809/02_shap-individuals.pdf", plot=p2,
  width = 17, height = 12, unit = "cm", scale = 1.3
)

################################################################################




# -----------------------------------------------------------------------------
# explore some of the most predictive features

## age

(shap %>%
    filter(fname == "age") %>%
    collect() %>%
    mutate(
        timepoint = fct_reorder(timepoint, lubridate::duration(timepoint))
    ) %>%
    ggplot(aes(x = vals, y = shap)) +
    geom_jitter(size = .05, alpha = .5) +
    facet_wrap(~timepoint, nrow = 1) +
    labs(
        title = "Effect of age on model output",
        x = "Feature value",
        y = "SHAP value (impact on model output)"
    )
 -> p2)

ggsave(
    "plots/shap/s008_effect-of-age.pdf", plot=p2,
    width=14, height=8, unit="cm", scale = 1.3
)


## vessels

(shap %>%
    filter(fname == "vessels") %>%
    mutate(
        value = sql("substring(feature from '[^_]+$')")
    ) %>%
    group_by(timepoint, pid) %>%
    mutate(shap = sum(shap)) %>%
    filter(vals == 1) %>%
    collect() %>% ungroup() %>%
    mutate(
        timepoint = fct_reorder(timepoint, lubridate::duration(timepoint)),
        value = fct_relevel(value, "DIF", "1VD", "2VD", "3VD")
    ) %>%
    ggplot(aes(y = value, x = shap, color = value)) +
    geom_quasirandom(
        alpha=.3, varwidth=TRUE, size=.5, width=.5,
        groupOnX=FALSE
    ) +
    facet_wrap(~timepoint, ncol=1) +
    scale_color_brewer(palette = "Set1") +
    theme(
        legend.position = "none"
    ) +
    labs(
        title = "Effect of vessel status on model output",
        y = "",
        x = "SHAP value (impact on model output)"
    ) -> p3
)

ggsave(
    "plots/shap/s008_effect-of-vessels.pdf", plot=p3,
    width=10, height=15, unit="cm", scale = 1.3
)






shap %>%
    filter(sql("feature ~ 'height|weight|sex'")) %>%
    group_by(pid) %>%
    mutate(
        shap = sum(shap)
    ) %>%
    collect() %>%
    pivot_wider(
        names_from = feature,
        values_from = vals
    ) %>%
    mutate(
        weight = base2_weight,
        height = base2_height,
        sex = if_else(base1_sex == 1, "male", "female")
    ) %>%
    ggplot(aes(weight, height, z = shap)) +
    stat_summary_2d(
        fun = sum, bins = 50
    ) +
    labs(
        fill = "SHAP", x = "Weight (Z-score)", y = "Height (z-score)"
    ) +
    facet_wrap(~sex) +
    scale_fill_gradient2(
        high = scales::muted("blue"),
        low = scales::muted("red")
    )

shap %>%
    filter(sql("not feature ~* 'time'")) %>%
    filter(sql("feature ~ 'bioc_NPU19748'")) %>%
    group_by(pid) %>%
    mutate(shap = sum(shap)) %>%
    filter(vals != 0) %>%
    ggplot(aes("feature", shap, color = feature)) +
    geom_quasirandom(cex=.5, method = "quasirandom", alpha=.5) +
    coord_flip() +
    labs(
        x = "",
        y = "SHAP value (impact on model output)"
    )


shap %>%
    filter(sql("feature ~ 'base2_smoking'")) %>%
    group_by(pid) %>%
    mutate(
        shap = sum(shap)
    ) %>%
    filter(vals != 0) %>%
    ungroup() %>%
    collect() %>%
    mutate(
        feature = str_extract(feature, "[^_]+$")
    ) %>%
    ggplot(aes(feature, shap, color = feature)) +
    geom_quasirandom(cex=.5, method = "quasirandom", alpha=.5) +
    coord_flip() +
    labs(
        x = "",
        y = "SHAP value (impact on model output)"
    ) +
    guides(
        color = guide_legend(override.aes = list(size = 5))
    ) +
    scale_color_brewer(palette = "Set1")

Sshap %>%
    filter(sql("feature ~ 'base2_nyha'")) %>%
    group_by(pid) %>%
    mutate(
        shap = sum(shap)
    ) %>%
    filter(vals != 0) %>%
    ungroup() %>%
    collect() %>%
    mutate(
        feature = str_extract(feature, "[^_]+$")
    ) %>%
    ggplot(aes(feature, shap)) +
    geom_quasirandom(cex=.5, method = "quasirandom", alpha=.5) +
    coord_flip() +
    labs(
        x = "",
        y = "SHAP value (impact on model output)"
    ) +
    guides(
        color = guide_legend(override.aes = list(size = 5))
    ) +
    scale_color_brewer(palette = "Set1")

shap <- tbl(con, in_schema("pmhnet", "explained"))

df <- shap %>% collect()


global <- shap %>%
    group_by(timepoint, feature) %>%
    summarise(
        direction = sum(shap), shap = sum(abs(shap))
    ) %>%
    mutate(shap = shap / sum(shap)) %>%
    collect() %>%
    mutate(
        direction = sign(direction) %>% as.character() %>%
            recode("-1" = "non-survival", "1" = "survival"),
        rank = row_number(desc(shap))
    ) %>%
    ungroup()


#------------------------------------------------------------------------------
# top10 most predictive (global)
# shap distributions per timepoint

shap_top <- shap %>%
    inner_join(
        select(global, -shap), by = c("feature", "timepoint"), copy = TRUE
    ) %>%
    filter(rank <= 10) %>%
    collect()

(p2 <- shap_top %>%
    mutate(
        vals = as.character(vals),
        timepoint = fct_reorder(timepoint, lubridate::duration(timepoint))
    ) %>%
    ggplot(aes(y = shap, x = as.character(rank))) +
    geom_quasirandom(aes(color = vals), size = .1) +
    stat_summary(
        aes(label = feature), fun = ~max(.) + .01, geom = "text"
    ) +
    coord_flip() +
    facet_wrap(
        ~timepoint, ncol = 1, labeller = label_both
    ) +
    scale_x_discrete(limits = rev) +
    scale_color_brewer(palette = "Set1", direction = -1) +
    labs(
        y = "Feature importance (SHAP)",
        caption = "model_001"
    ) +
    theme(
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "none"
    )
)

ggsave(
    "plots/shap/s001-top10-global.pdf", plot=p2,
    width=10, height=15, unit="cm", scale = 1.3
)
