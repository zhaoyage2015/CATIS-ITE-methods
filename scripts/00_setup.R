# ============================================================
# 00_setup.R  (Methods-only / Representative script)
#
# Purpose:
#   This script documents the main R package dependencies used
#   in the CATIS individualized treatment effect (ITE) analysis.
#
# Notes:
#   - This repository provides representative scripts for method
#     transparency; individual-level trial data are not included.
#   - This setup script is illustrative and is not intended to be a
#     fully automated, one-click reproducible environment installer.
#   - No local file paths are referenced.
# ============================================================

message("▶ Loading representative analysis dependencies (illustrative)...")

## Optional: set a CRAN mirror (safe to keep, but not required)
## If you prefer to avoid region-specific settings, you can comment this out.
options(repos = c(CRAN = "https://cran.r-project.org"))

## ---- Core dependencies used in the analysis ----
pkgs <- c(
  # Data wrangling / IO
  "tidyverse", "readxl", "openxlsx", "janitor", "skimr", "naniar", "fs", "glue",

  # Tables / reporting
  "gtsummary", "tableone", "gt", "broom",

  # Visualization
  "ggplot2", "patchwork", "cowplot", "ggdist",

  # Causal ML / uplift evaluation
  "grf", "tools4uplift",

  # Metrics / utilities
  "pROC", "DescTools",

  # ML workflow helpers / interpretability
  "rsample", "rlang", "purrr", "iml"
)

## ---- Load packages if installed ----
## (Do not auto-install here to avoid implying a one-click runnable pipeline)
missing_pkgs <- character(0)

for (p in pkgs) {
  if (requireNamespace(p, quietly = TRUE)) {
    suppressPackageStartupMessages(library(p, character.only = TRUE))
  } else {
    missing_pkgs <- c(missing_pkgs, p)
  }
}

if (length(missing_pkgs) > 0) {
  message("⚠ The following packages are not installed in this R environment:\n  - ",
          paste(missing_pkgs, collapse = "\n  - "))
  message("You may install them as needed, e.g. install.packages(c(...)).")
} else {
  message("✅ All listed packages are available and loaded.")
}

## ---- Optional: renv note (documented but not enforced) ----
## The full project uses renv for dependency management. For transparency,
## we mention it here without running renv::init/restore/snapshot automatically.
if (requireNamespace("renv", quietly = TRUE)) {
  message("ℹ renv is available. (Full pipeline can be provided upon reasonable request.)")
} else {
  message("ℹ renv is not installed in this environment (optional).")
}

message("✅ 00_setup.R completed.")
