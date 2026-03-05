#!/usr/bin/env Rscript
# =============================================================================
# 09_sensitivity_predictors_stability.R
# Post-processing for sensitivity analysis (repeated random splits):
# - Read outputs/sensitivity_repeated_split_v3_noqini/vi_grouped_long.csv
# - Compute Top10 / Top15 frequency across splits
# - Harmonize display names to match main text
# - Plot stability bar charts (bright blue, no title) and save (Top10 + Top15)
#
# Outputs (under outputs/sensitivity_repeated_split_v3_noqini/):
#   - vi_top10_freq.csv
#   - vi_top15_freq.csv
#   - vi_top10_stability_mainstyle.png
#   - vi_top15_stability_mainstyle.png
#
# Project root:
#   - Option A: set env var CATIS_ITE_DIR
#   - Option B: auto-detect by walking up from script location / working dir
# =============================================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(scales)
  library(stringr)
})

# --------------------------- Project root ---------------------------
get_proj_dir <- function() {
  env <- Sys.getenv("CATIS_ITE_DIR", unset = NA_character_)
  if (!is.na(env) && nzchar(env) && dir.exists(env)) {
    return(normalizePath(env, winslash = "/", mustWork = TRUE))
  }

  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  script_dir <- if (length(file_arg) > 0) dirname(sub("^--file=", "", file_arg)) else getwd()

  markers <- c("data", "models", "outputs", "figures", ".git")
  cur <- normalizePath(script_dir, winslash = "/", mustWork = TRUE)

  for (i in 1:10) {
    if (any(file.exists(file.path(cur, markers)))) return(cur)
    parent <- dirname(cur)
    if (identical(parent, cur)) break
    cur <- parent
  }
  normalizePath(getwd(), winslash = "/", mustWork = TRUE)
}

proj_dir <- get_proj_dir()

# ------------------------------ Paths ------------------------------
OUTDIR <- file.path(proj_dir, "outputs", "sensitivity_repeated_split_v3_noqini")
IN_VI  <- file.path(OUTDIR, "vi_grouped_long.csv")

OUT_CSV_10 <- file.path(OUTDIR, "vi_top10_freq.csv")
OUT_CSV_15 <- file.path(OUTDIR, "vi_top15_freq.csv")

OUT_PNG_10 <- file.path(OUTDIR, "vi_top10_stability_mainstyle.png")
OUT_PNG_15 <- file.path(OUTDIR, "vi_top15_stability_mainstyle.png")

dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)
if (!file.exists(IN_VI)) stop("Missing input file: ", IN_VI, call. = FALSE)

# ------------------------ Display harmonization ------------------------
# Only affects plotting/tables (not modeling)
pretty_display <- c(
  "ty_g"         = "TyG",
  "hs-CRP"       = "hs-CRP",
  "Pro-BNP"      = "NT-proBNP",   # unify with main text
  "TC"           = "TC",
  "Creatinine"   = "Creatinine",
  "BMI"          = "BMI",
  "HCT"          = "HCT",
  "Glucose"      = "Glucose",
  "WBC"          = "WBC",
  "Uric acid"    = "Uric acid",
  "Age"          = "Age",
  "WC"           = "WC",
  "Hcy"          = "Hcy",
  "Heart rate"   = "Heart rate",
  "Baseline SBP" = "Baseline SBP",
  "Baseline DBP" = "Baseline DBP",
  "eGFR"         = "eGFR"
)

harmonize_display <- function(x) {
  out <- x
  hit <- out %in% names(pretty_display)
  out[hit] <- unname(pretty_display[out[hit]])
  out
}

# ------------------------------ Read ------------------------------
vi <- read_csv(IN_VI, show_col_types = FALSE)

# Basic schema check
needed_cols <- c("iter", "rank", "display", "base")
miss <- setdiff(needed_cols, names(vi))
if (length(miss) > 0) stop("Input missing required columns: ", paste(miss, collapse = ", "), call. = FALSE)

# Fill missing/blank display with base
missing_disp <- is.na(vi$display) | str_trim(vi$display) == ""
vi$display[missing_disp] <- vi$base[missing_disp]

n_splits <- length(unique(vi$iter))
if (n_splits <= 1) warning("n_splits seems small (", n_splits, "). Check 'iter' column.")

# -------------------------- Frequency tables --------------------------
calc_top_freq <- function(df, k) {
  nm_n <- paste0("n_top", k)
  df %>%
    filter(.data$rank <= k) %>%
    count(.data$display, name = nm_n) %>%
    mutate(
      n_splits = n_splits,
      freq = .data[[nm_n]] / n_splits
    ) %>%
    arrange(desc(freq))
}

vi_top10_freq <- calc_top_freq(vi, 10)
vi_top15_freq <- calc_top_freq(vi, 15)

write_csv(vi_top10_freq, OUT_CSV_10)
write_csv(vi_top15_freq, OUT_CSV_15)

message("Saved: ", OUT_CSV_10)
message("Saved: ", OUT_CSV_15)

# ------------------------------ Plotting ------------------------------
# Style per your requirement: bright blue, no title
make_stability_plot <- function(df_freq, ylab, freq_threshold = 0.5) {
  df_plot <- df_freq %>%
    mutate(display = harmonize_display(display)) %>%
    filter(freq >= freq_threshold)

  ggplot(df_plot, aes(x = reorder(display, freq), y = freq)) +
    geom_col(fill = "#1E90FF") +
    coord_flip() +
    scale_y_continuous(
      limits = c(0, 1),
      breaks = seq(0, 1, 0.1),
      labels = percent_format(accuracy = 1)
    ) +
    labs(x = NULL, y = ylab) +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank()
    )
}

p10 <- make_stability_plot(
  df_freq = vi_top10_freq,
  ylab = "Proportion of splits in Top 10",
  freq_threshold = 0.5
)

p15 <- make_stability_plot(
  df_freq = vi_top15_freq,
  ylab = "Proportion of splits in Top 15",
  freq_threshold = 0.5
)

ggsave(OUT_PNG_10, plot = p10, width = 8, height = 5, dpi = 600)
ggsave(OUT_PNG_15, plot = p15, width = 8, height = 5, dpi = 600)

message("Saved figure: ", OUT_PNG_10)
message("Saved figure: ", OUT_PNG_15)
message("Done. Project root used: ", proj_dir)
