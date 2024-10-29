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

data <- inner_join(
    tbl(con, "shap_006"),
    tbl(con, "shap_012"),
    join_by(pid, feature, timepoint, vals)
  ) |> 
  mutate(
    feature = sql("substring(feature from '^[^_]+_[^_]+')")
  ) |> 
  summarise(
    shap.x = sum(shap.x), 
    shap.y = sum(shap.y), 
    .by = c(pid, feature)
  ) |> 
  mutate(
    fgroup = sql("substring(feature from '^[^_]+')"),
    fname  = sql("substring(feature from '^[^_]+_([^_]+)')"),
    fgroup = if_else(fgroup == "diag-simple", "diag-1", fgroup)
  ) |> 
  collect()

colors <- c(
  "bioc-1" = "#ebb391",
  "diag-1" = "#7ba79d",
  "clnc-1" = "#e8d19d",
  "clnc-2" = "#bfd3e4",
  "proc-1" = "#cadcd8"
)

data |> 
  mutate(
    delta = shap.y - shap.x,
    q50 = median(delta),
    .by = fname
  ) |> 
  mutate(
    rank = dense_rank(desc(abs(q50)))
  ) |> 
  filter(rank <= 25) |> 
  mutate(
    fname = fct_reorder(fname, delta)
  ) |> 
  ggplot(aes(x  = fname, y = delta, fill = fgroup)) +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = .5) +
  geom_boxplot(outlier.shape = NA) + 
  scale_y_continuous(
    limits = c(-0.06, 0.06), 
    labels = scales::label_percent(),
    breaks = scales::pretty_breaks()
  ) +
  scale_fill_manual(values = colors) +
  labs(
    title = "top-25 most affected features",
    y = "Difference in SHAP", x = ""
  ) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1),
    legend.position = "none"
  )

ggsave(
  "plots/misc/231004_shap-diff.pdf",
  width = 200, height = 100, scale = .9, unit = "mm"
)
  