library(tidyverse)
library(survival)
library(dbplyr)

theme_set(theme_light())

study <- snakemake@wildcards$sid

# create output dir if non-existant
dir.create(snakemake@output[[1]])

con <- DBI::dbConnect(odbc::odbc(),
    driver = "PostgreSQL", server = "dbserver", database = "database",
    uid = Sys.getenv("uid"), pwd = Sys.getenv("pwd")
)

surv <- tbl(con, ident_q("pmhnet.surv_master"))
grace <- tbl(con, ident_q("pmhnet.risk_scores"))

# model predictions on test data
p_test <- read_csv(snakemake@input$p_test, col_types = cols()) %>%
    pivot_longer(-pid, names_to = "time", values_to = "risk") %>%
    mutate(time = as.numeric(time))

# observed survival on all patients
surv <- collect(surv) %>%
    mutate(
        across(everything(), as.integer),
        event = if_else(time > 1825, 0L, event),
        time = if_else(time > 1825, 1825L, time)
    )

# observed survival on test set
s_test <- semi_join(surv, p_test, by = "pid")

# predicited risk for 0.5, 1, 3, and 5 years
p_test_tau <- p_test %>%
    group_by(pid) %>%
    summarise(
        tau = 365 * c(0.5, 1, 3, 5),
        pred = approx(x = time, y = risk, xout = tau)$y,
        .groups = "drop"
    )

# calibration plot 01 ---------------------------------------------------------
p01 <- p_test_tau %>%
    inner_join(surv, by = "pid") %>%
    group_by(tau) %>%
    mutate(
        strata = cut_number(pred, n = 20) %>% as.numeric()
    ) %>%
    group_by(strata, .add = TRUE) %>%
    summarise(
        pred = mean(pred),
        kmfit = survfit(Surv(time, event) ~ 1, data = cur_data()) %>%
            summary(time = cur_group()$tau) %>% list(),
        .groups = "drop"
    ) %>%
    hoist(kmfit, "surv", "upper", "lower") %>%
    select(-kmfit) %>%
    mutate(
        tau = as_factor(tau) %>%
            lvls_revalue(c("6 months", "1 year", "3 year", "5 year"))
    ) %>%
    ggplot(aes(pred, surv)) + geom_point(size = 1) +
    geom_errorbar(
        aes(ymin = lower, ymax = upper), width = 1/50, size = .5
    ) +
    facet_wrap(~tau, nrow = 1) +
    geom_abline(slope = 1, linetype = "dashed", alpha = .5) +
    labs(
        x = "Predicted Survival", y = "Observed Survival",
        title = glue::glue("Calibration of pmhnet model_{study}"),
        caption = "points: 5% risk quantiles"
    )

ggsave(
    filename = file.path(snakemake@output[[1]], "model_calibration_01.pdf"),
    plot = p01, height = 8, width = 20, unit = "cm", scale = 0.9
)

# calibration 02 --------------------------------------------------------------

# average predicted risk
p_test_ave <- p_test %>%
    group_by(tau = time) %>%
    summarise(surv = mean(risk), .groups = "drop")

# Kaplan-Meier estimates
km_test <- survfit(Surv(time, event) ~ 1, data = s_test) %>%
    {summary(.)[c("time", "surv", "upper", "lower")]} %>%
    bind_cols() %>%
    rename(tau = time)

# survival curves: km-estimates vs average predicted
p02 <- bind_rows(`Kaplan-Meier` = km_test, Prediction = p_test_ave, .id = "type") %>%
    mutate(type = fct_relevel(type, "Prediction")) %>%
    ggplot(aes(tau/365, surv)) +
    geom_step(aes(linetype = type)) +
    geom_ribbon(
        aes(ymin = lower, ymax = upper), alpha = .25,
        data = ~filter(., type == "Kaplan-Meier")
    ) +
    labs(
        x = "Years", y = "Survival Probability (%)",
        title = glue::glue("Calibration of pmhnet model_{study}"),
        linetype = "Estimate"
    ) +
    theme(
        legend.position = c(0.8, 0.75),
        legend.background = element_blank()
    )

ggsave(
    filename = file.path(snakemake@output[[1]], "model_calibration_02.pdf"),
    plot = p02, height = 8, width = 15, unit = "cm", scale = 0.9
)

# calibration 03 --------------------------------------------------------------

p03 <- surv %>%
    semi_join(p_test, by = "pid") %>%
    summarise(
        kmfit = survfit(Surv(time, event) ~ 1, data = cur_data()) %>%
            summary(time = p_test_ave$tau) %>%
            `[`(c("time", "surv", "upper", "lower")) %>%
            bind_cols()
    ) %>%
    unpack(kmfit) %>%
    rename(tau = time) %>%
    inner_join(p_test_ave, by = "tau") %>%
    ggplot(aes(surv.y, surv.x)) +
    geom_point() +
    geom_errorbar(aes(ymin = lower, ymax = upper)) +
    geom_abline(slope = 1, linetype = "dashed") +
    ggrepel::geom_text_repel(
        aes(label = round(tau)), direction = "both",
        nudge_y = -0.01, nudge_x = 0.01, alpha = .25
    ) +
    labs(
        x = "Predicted Survival (%)", y = "Observed Survival (%)",
        title = glue::glue("Calibration of pmhnet model_{study} over time"),
        caption = "points: mean predictions at varying timepoints"
    )

ggsave(
    filename = file.path(snakemake@output[[1]], "model_calibration_03.pdf"),
    plot = p03, height = 10, width = 10, unit = "cm", scale = 1.1
)

# model discrimination 01 -----------------------------------------------------

roc_01 <- p_test_tau %>%
    inner_join(surv, by = "pid") %>%
    group_by(tau) %>%
    summarise(
        res = list(tdROC::tdROC(
                X = -pred, Y = time, delta = event, tau = cur_group()$tau
        )), .groups = "drop"
    ) %>%
    mutate(
        roc = map(res, "ROC"),
        auc = map(res, "AUC"),
        tau = as_factor(tau) %>%
            lvls_revalue(c("6 months", "1 year", "3 year", "5 year"))
    ) %>%
    select(-res)

roc_01 %>%
    unnest(auc) %>%
    distinct(tau, auc=value) %>%
    write_tsv(file.path(snakemake@output[[1]], "auc_01.tsv"))

p04 <- roc_01 %>%
    unnest(c(roc, auc)) %>%
    ggplot(aes(1-spec, sens)) +
    geom_step() +
    geom_abline(slope = 1, linetype = "dashed", alpha = .5) +
    facet_wrap(~tau, nrow = 1) +
    geom_text(
        data = ~distinct(., tau, value),  x = 0.75, y = 0.25,
        aes(label = glue::glue("AUC: {x}", x = round(value, digits = 3)))
    ) +
    coord_fixed() +
    labs(
        x = "1 - Specificity", y = "Sensitivity",
        title = glue::glue("Performance of pmhnet model_{study}")
    )

ggsave(
    filename = file.path(snakemake@output[[1]], "model_discrimination_01.pdf"),
    plot = p04, height = 7, width = 16, unit = "cm", scale = 1.2
)

# model discrimination 02 -----------------------------------------------------
# comparing with grace-scores

grace_test <- collect(grace) %>% inner_join(s_test, by = "pid")

roc_02 <- grace_test %>%
    rename(tau = time.x) %>%
    mutate(tau = if_else(tau == 183, 182.5, as.double(tau))) %>%
    inner_join(p_test_tau, by = c("pid", "tau")) %>%
    rename(grace = value, pmhnet = pred, time = time.y) %>%
    mutate(pmhnet = -pmhnet) %>%
    pivot_longer(c(grace, pmhnet)) %>%
    group_by(tau, name) %>%
    summarise(
        res = list(tdROC::tdROC(
                X = value, Y = time, delta = event, tau = cur_group()$tau
        )), .groups = "drop"
    ) %>%
    mutate(roc = map(res, "ROC"), auc = map(res, "AUC")) %>%
    select(-res)

roc_02 %>%
    unnest(auc) %>%
    distinct(tau, name, auc=value) %>%
    write_tsv(file.path(snakemake@output[[1]], "auc_02.tsv"))

p05 <- roc_02 %>%
    mutate(
        tau = as_factor(tau) %>%
            lvls_revalue(c("6 months", "1 year", "3 year"))
    ) %>%
    unnest(c(roc, auc)) %>%
    ggplot(aes(1-spec, sens)) +
    geom_step(aes(color = name)) +
    geom_abline(slope = 1, linetype = "dashed", alpha = .5) +
    facet_wrap(~tau, nrow = 1) +
    geom_text(
        data = ~filter(., name == "grace"), stat = "unique",
        x = 0.75, y = 0.1,
        aes(label = glue::glue("AUC(g): {x}", x = round(value, digits = 3)))
    ) +
    geom_text(
        data = ~filter(., name == "pmhnet"), stat = "unique",
        x = 0.75, y = 0.2,
        aes(label = glue::glue("AUC(p): {x}", x = round(value, digits = 3)))
    ) +
    coord_fixed() +
    labs(
        x = "1 - Specificity", y = "Sensitivity", color = "",
        title = glue::glue("Performance of pmhnet model_{study}"),
        subtitle = "+ comparison with GRACE2.0"
    )

ggsave(
    filename = file.path(snakemake@output[[1]], "model_discrimination_02.pdf"),
    plot = p05, height = 7, width = 15, unit = "cm", scale = 1.3
)

# model discrimination 03 -----------------------------------------------------
# comparing with grace-scores and stratifying on imputation status

roc_03 <- grace_test %>%
    rename(tau = time.x) %>%
    mutate(tau = if_else(tau == 183, 182.5, as.double(tau))) %>%
    inner_join(p_test_tau, by = c("pid", "tau")) %>%
    rename(grace = value, pmhnet = pred, time = time.y) %>%
    mutate(pmhnet = -pmhnet) %>%
    pivot_longer(c(grace, pmhnet)) %>%
    group_by(tau, name, imp = imputed) %>%
    summarise(
        res = list(tdROC::tdROC(
                X = value, Y = time, delta = event, tau = cur_group()$tau
        )), .groups = "drop"
    ) %>%
    mutate(roc = map(res, "ROC"), auc = map(res, "AUC")) %>%
    select(-res)

roc_03 %>%
    unnest(auc) %>%
    distinct(tau, name, imp, auc=value) %>%
    write_tsv(file.path(snakemake@output[[1]], "auc_03.tsv"))

p06 <- roc_03 %>%
    mutate(
        tau = as_factor(tau) %>%
            lvls_revalue(c("6 months", "1 year", "3 year")),
        imp = if_else(imp == 0, "no imputation", "imputation")
    ) %>%
    unnest(c(roc, auc)) %>%
    ggplot(aes(1-spec, sens)) +
    geom_step(aes(color = name)) +
    geom_abline(slope = 1, linetype = "dashed", alpha = .5) +
    facet_grid(imp~tau) +
    geom_text(
        data = ~filter(., name == "grace"), stat = "unique",
        x = 0.75, y = 0.1,
        aes(label = glue::glue("AUC(g): {x}", x = round(value, digits = 3)))
    ) +
    geom_text(
        data = ~filter(., name == "pmhnet"), stat = "unique",
        x = 0.75, y = 0.2,
        aes(label = glue::glue("AUC(p): {x}", x = round(value, digits = 3)))
    ) +
    coord_fixed() +
    labs(
        x = "1 - Specificity", y = "Sensitivity", color = "",
        title = glue::glue("Performance of pmhnet model_{study}"),
        subtitle = "+ comparison with GRACE2.0"
    ) +
    theme(legend.position = "bottom")

ggsave(
    filename = file.path(snakemake@output[[1]], "model_discrimination_03.pdf"),
    plot = p06, height = 12, width = 13, unit = "cm", scale = 1.3
)
