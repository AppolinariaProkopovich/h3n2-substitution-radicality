#!/usr/bin/env Rscript

cran_packages <- c(
  "ape",
  "broom",
  "castor",
  "dplyr",
  "ggplot2",
  "grantham",
  "lme4",
  "lmerTest",
  "purrr",
  "readr",
  "stringr",
  "tidyverse"
)

missing <- cran_packages[!vapply(
  cran_packages,
  requireNamespace,
  logical(1),
  quietly = TRUE
)]

if (length(missing) == 0) {
  message("All R dependencies are already installed.")
} else {
  install.packages(missing, repos = "https://cloud.r-project.org")
}
