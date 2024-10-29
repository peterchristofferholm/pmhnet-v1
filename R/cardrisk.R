library(readr)
library(tibble)
library(dplyr)
library(assertthat)

m1_probs <- function(score) {

    lookup <- "R/data/grace_6m-probs.csv" %>%
        read_csv(col_types = "ddd")

    probability <- lookup %>%
        filter(from < score, score <= to) %>%
        pull("probability")

    # Following adjustments are not explained, but are present in web-tool source
    probability <- probability * 100 * (80 / 91)

    if (probability > 90) {
        probability <- 92  # capping to 92.
    }

    probability
}

m2_probs <- function(score) {
    s0 <- 0.9983577131
    (1 - s0^exp(score)) * 100
}

m3_probs <- function(score) {
    s0 <- 0.9998715509;
    (1 - s0^exp(score)) * 100
}

grace_v2 <- function(age, pulse, sbp, crea, killip, arrest, enzymes, stemi) {

    assert_that(
        age %>% between(18, 110),
        pulse %>% between(10, 300),
        sbp %>% between(10, 300),
        crea %>% between(0, 10),
        killip %in% 1:4,
        arrest %in% 0:1,
        enzymes %in% 0:1,
        stemi %in% 0:1
    )

    # Put vars in a tibble
    vars <- c(age = age, pulse = pulse, sbp = sbp, creatinine = crea,
              killip = killip, arrest = arrest, enzymes = enzymes,
              stemi = stemi) %>%
        enframe(name = "term", value = "value") %>%
        mutate(value = as.numeric(value))

    # Lookup xbetas in the .csv file
    points_lookup <- "R/data/grace_all-points.csv" %>%
        read_csv(col_types = "ddcc")

    # Return predicted probabilities
    vars %>%
        full_join(points_lookup, by = "term") %>%
        mutate(diff = abs(value.y - value.x)) %>%
        group_by(model, term) %>%
        summarise(score = score[diff == min(diff)], .groups = "drop_last") %>%
        summarise(score = sum(score), .groups = "keep") %>%
        transmute(
            prob = case_when(
                model == "6m" ~ m1_probs(score),
                model == "1y" ~ m2_probs(score),
                model == "3y" ~ m3_probs(score)
            )
        ) %>%
        ungroup()

}
