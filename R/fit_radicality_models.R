#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(lme4)
  library(lmerTest)
  library(readr)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop(
    "Usage: Rscript R/fit_radicality_models.R <input.csv> <output.csv>",
    call. = FALSE
  )
}

input_path <- args[[1]]
output_path <- args[[2]]
required_columns <- c(
  "distance", "site", "anc_aa", "gr_rad_anc", "miyata_rad"
)

data <- read_csv(input_path, show_col_types = FALSE)
missing_columns <- setdiff(required_columns, names(data))
if (length(missing_columns) > 0) {
  stop(
    paste("Missing required columns:", paste(missing_columns, collapse = ", ")),
    call. = FALSE
  )
}

data <- data %>%
  mutate(
    site = factor(site),
    anc_aa = factor(anc_aa)
  )

model_r_squared <- function(model) {
  # The thesis compares plain R^2 for linear models with the marginal and
  # conditional R^2 of Nakagawa & Schielzeth (MuMIn) for mixed models. A linear
  # model has no random part, so both of its values are the same R^2.
  if (!inherits(model, "merMod")) {
    r_squared <- summary(model)$r.squared
    return(c(marginal = r_squared, conditional = r_squared))
  }
  if (!requireNamespace("MuMIn", quietly = TRUE)) {
    warning("MuMIn is not installed; R^2 columns are reported as NA")
    return(c(marginal = NA_real_, conditional = NA_real_))
  }
  values <- MuMIn::r.squaredGLMM(model)
  c(marginal = unname(values[1, "R2m"]), conditional = unname(values[1, "R2c"]))
}

extract_fixed_effect <- function(model, response, specification, n_obs) {
  coefficients <- coef(summary(model))
  distance <- coefficients["distance", ]
  is_mixed <- inherits(model, "merMod")
  r_squared <- model_r_squared(model)

  tibble(
    response = response,
    model = specification,
    n = n_obs,
    estimate = unname(distance[["Estimate"]]),
    std_error = unname(distance[["Std. Error"]]),
    degrees_freedom = if ("df" %in% names(distance)) {
      unname(distance[["df"]])
    } else {
      stats::df.residual(model)
    },
    statistic = if ("t value" %in% names(distance)) {
      unname(distance[["t value"]])
    } else {
      unname(distance[["z value"]])
    },
    p_value = if ("Pr(>|t|)" %in% names(distance)) {
      unname(distance[["Pr(>|t|)"]])
    } else {
      unname(distance[["Pr(>|z|)"]])
    },
    r_squared_marginal = unname(r_squared[["marginal"]]),
    r_squared_conditional = unname(r_squared[["conditional"]]),
    aic = AIC(model),
    bic = BIC(model),
    reml_criterion = if (is_mixed) lme4::REMLcrit(model) else NA_real_,
    singular = if (is_mixed) lme4::isSingular(model) else FALSE
  )
}

fit_response <- function(source, response) {
  complete <- source %>%
    filter(
      is.finite(.data[[response]]),
      is.finite(distance),
      !is.na(site),
      !is.na(anc_aa)
    ) %>%
    droplevels()

  site_filtered <- complete %>%
    group_by(site) %>%
    filter(n() >= 3) %>%
    ungroup() %>%
    droplevels()

  nested_filtered <- complete %>%
    group_by(site, anc_aa) %>%
    filter(n() >= 3) %>%
    ungroup() %>%
    droplevels()

  if (nrow(complete) < 3 || n_distinct(complete$distance) < 2) {
    stop(response, ": insufficient complete observations or distance variation")
  }
  if (n_distinct(site_filtered$site) < 2) {
    stop(response, ": the site model requires at least two retained sites")
  }
  nested_groups <- interaction(
    nested_filtered$site,
    nested_filtered$anc_aa,
    drop = TRUE
  )
  if (n_distinct(nested_groups) < 2) {
    stop(response, ": the nested model requires at least two retained groups")
  }

  fixed_formula <- reformulate("distance", response = response)
  site_formula <- as.formula(paste(response, "~ distance + (1 | site)"))
  nested_formula <- as.formula(
    paste(response, "~ distance + (1 | site / anc_aa)")
  )

  models <- list(
    list(
      fit = lm(fixed_formula, data = complete),
      name = "fixed_effect_only",
      n = nrow(complete)
    ),
    list(
      fit = lmer(site_formula, data = site_filtered, REML = TRUE),
      name = "random_intercept_site",
      n = nrow(site_filtered)
    ),
    list(
      fit = lmer(nested_formula, data = nested_filtered, REML = TRUE),
      name = "random_intercept_site_parent_aa",
      n = nrow(nested_filtered)
    )
  )

  bind_rows(lapply(models, function(item) {
    extract_fixed_effect(item$fit, response, item$name, item$n)
  }))
}

results <- bind_rows(lapply(c("gr_rad_anc", "miyata_rad"), function(response) {
  fit_response(data, response)
}))

dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
write_csv(results, output_path)
message("Wrote ", nrow(results), " model summaries to ", output_path)
