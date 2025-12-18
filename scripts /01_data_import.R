# ============================================================
# 01_data_import.R  (Methods-only / Representative script)
#
# Purpose:
#   Illustrative data import and basic QC workflow for the CATIS ITE analysis.
#
# Notes:
#   - Individual-level trial data are not included in this repository.
#   - File paths and filenames below are placeholders for transparency.
# ============================================================

source("scripts/00_setup.R")

suppressPackageStartupMessages({
  library(fs)
  library(readxl)
  library(tidyverse)
  library(janitor)
  library(purrr)
  library(stringr)
  library(vctrs)
  library(readr)
  library(glue)
})

# ----------------------------- Paths (placeholders) --------------------------
dir_raw <- file.path("data", "raw")     # placeholder location for local data (not included)
dir_out <- file.path("outputs", "qc")   # QC summaries written here

fs::dir_create(c(dir_raw, dir_out))

# Placeholder file names (DATA NOT INCLUDED in this repository)
f_train      <- file.path(dir_raw, "train.xlsx")
f_valid      <- file.path(dir_raw, "valid.xlsx")
f_train_spec <- file.path(dir_raw, "train_with_additional_markers.xlsx")
f_valid_spec <- file.path(dir_raw, "valid_with_additional_markers.xlsx")

if (!all(file.exists(c(f_train, f_valid, f_train_spec, f_valid_spec)))) {
  missing_files <- c(f_train, f_valid, f_train_spec, f_valid_spec)[
    !file.exists(c(f_train, f_valid, f_train_spec, f_valid_spec))
  ]
  stop(
    "Data files are not available in this public repository.\n",
    "Missing placeholder files:\n  - ",
    paste(missing_files, collapse = "\n  - "),
    "\n\nPlease place the corresponding local files under: ", dir_raw,
    "\n(or adapt the file paths for your local environment)."
  )
}

# ----------------------------- Utility functions -----------------------------
norm_id <- function(df, id_col = "id") {
  if (!id_col %in% names(df)) {
    stop(glue("Missing column '{id_col}'. Please ensure an ID column exists in the input data."))
  }
  df %>%
    mutate(
      "{id_col}" := as.character(.data[[id_col]]),
      "{id_col}" := stringr::str_trim(.data[[id_col]])
    )
}

check_id <- function(df, tag, id_col = "id") {
  if (anyNA(df[[id_col]])) {
    nmis <- sum(is.na(df[[id_col]]))
    stop(glue("{tag}: '{id_col}' contains {nmis} missing values."))
  }
  if (any(duplicated(df[[id_col]]))) {
    ndups <- sum(duplicated(df[[id_col]]))
    stop(glue("{tag}: '{id_col}' contains {ndups} duplicated values."))
  }
  invisible(df)
}

col_profile <- function(df, tag) {
  tibble(
    table        = tag,
    var          = names(df),
    dtype        = map_chr(df, ~vctrs::vec_ptype_abbr(.x)),
    class        = map_chr(df, ~paste(class(.x), collapse = ",")),
    n_missing    = map_int(df, ~sum(is.na(.x))),
    missing_rate = map_dbl(df, ~mean(is.na(.x))),
    n_unique     = map_int(df, ~dplyr::n_distinct(.x, na.rm = TRUE))
  )
}

nzv_check <- function(df, tag) {
  numv <- df %>% select(where(is.numeric))
  if (ncol(numv) == 0) return(tibble())
  sds <- sapply(numv, sd, na.rm = TRUE)
  tibble(table = tag, var = names(sds), sd = as.numeric(sds)) %>%
    mutate(flag_zero_var = sd == 0 | is.na(sd))
}

miss_tbl <- function(df, tag) {
  tibble(
    table = tag,
    var = names(df),
    n_missing = colSums(is.na(df)),
    missing_rate = colMeans(is.na(df))
  ) %>% arrange(desc(missing_rate))
}

num_summary <- function(df, tag) {
  numv <- df %>% select(where(is.numeric))
  if (ncol(numv) == 0) return(tibble())
  map_dfr(names(numv), function(v) {
    x <- numv[[v]]
    tibble(
      table = tag, var = v,
      n   = sum(!is.na(x)),
      mean= mean(x, na.rm = TRUE),
      sd  = sd(x, na.rm = TRUE),
      p0  = quantile(x, 0,   na.rm = TRUE),
      p25 = quantile(x, .25, na.rm = TRUE),
      p50 = quantile(x, .50, na.rm = TRUE),
      p75 = quantile(x, .75, na.rm = TRUE),
      p100= quantile(x, 1,   na.rm = TRUE)
    )
  })
}

compare_cat <- function(df1, df2) {
  cats <- intersect(
    names(select(df1, where(~is.character(.x) || is.factor(.x)))),
    names(select(df2, where(~is.character(.x) || is.factor(.x))))
  )
  if (length(cats) == 0) return(tibble())
  map_dfr(cats, function(v) {
    a <- janitor::tabyl(df1[[v]]) %>%
      janitor::adorn_percentages("col") %>%
      rename(level = 1, train_prop = percent) %>% select(level, train_prop)
    b <- janitor::tabyl(df2[[v]]) %>%
      janitor::adorn_percentages("col") %>%
      rename(level = 1, valid_prop = percent) %>% select(level, valid_prop)
    full_join(a, b, by = "level") %>%
      mutate(var = v, diff = train_prop - valid_prop) %>%
      relocate(var)
  })
}

# ----------------------------- Read and harmonize ----------------------------
train      <- read_excel(f_train)
valid      <- read_excel(f_valid)
train_spec <- read_excel(f_train_spec)
valid_spec <- read_excel(f_valid_spec)

train      <- train      %>% clean_names() %>% norm_id("id")
valid      <- valid      %>% clean_names() %>% norm_id("id")
train_spec <- train_spec %>% clean_names() %>% norm_id("id")
valid_spec <- valid_spec %>% clean_names() %>% norm_id("id")

check_id(train, "train", "id")
check_id(valid, "valid", "id")
check_id(train_spec, "train_spec", "id")
check_id(valid_spec, "valid_spec", "id")

# Column name consistency (bidirectional)
stopifnot(length(setdiff(names(train), names(valid))) == 0)
stopifnot(length(setdiff(names(valid), names(train))) == 0)

stopifnot(length(setdiff(names(train_spec), names(valid_spec))) == 0)
stopifnot(length(setdiff(names(valid_spec), names(train_spec))) == 0)

# Align column order
valid      <- valid[, names(train)]
valid_spec <- valid_spec[, names(train_spec)]

# Type differences report (train vs valid)
type_cmp <- tibble(
  var       = names(train),
  cls_train = map_chr(train, ~paste(class(.x), collapse = ",")),
  cls_valid = map_chr(valid, ~paste(class(.x), collapse = ","))
) %>% filter(cls_train != cls_valid)

if (nrow(type_cmp) > 0) {
  warning("Type mismatch detected (train vs valid):\n",
          paste0(capture.output(print(type_cmp)), collapse = "\n"))
}

# ----------------------------- QC summaries ---------------------------------
profile_all <- bind_rows(
  col_profile(train, "train"),
  col_profile(valid, "valid"),
  col_profile(train_spec, "train_spec"),
  col_profile(valid_spec, "valid_spec")
)
write_csv(profile_all, file.path(dir_out, "column_profile.csv"))

nzv_all <- bind_rows(
  nzv_check(train, "train"),
  nzv_check(valid, "valid"),
  nzv_check(train_spec, "train_spec"),
  nzv_check(valid_spec, "valid_spec")
)
write_csv(nzv_all, file.path(dir_out, "zero_variance_check.csv"))

miss_all <- bind_rows(
  miss_tbl(train, "train"),
  miss_tbl(valid, "valid"),
  miss_tbl(train_spec, "train_spec"),
  miss_tbl(valid_spec, "valid_spec")
)
write_csv(miss_all, file.path(dir_out, "missing_report.csv"))

numsum_all <- bind_rows(
  num_summary(train, "train"),
  num_summary(valid, "valid"),
  num_summary(train_spec, "train_spec"),
  num_summary(valid_spec, "valid_spec")
)
write_csv(numsum_all, file.path(dir_out, "numeric_summary.csv"))

# Numeric compare (train vs valid)
qc_compare_num <- train %>%
  select(where(is.numeric)) %>%
  summarize(across(everything(), list(
    train_mean = ~mean(., na.rm = TRUE),
    train_p50  = ~quantile(., .5, na.rm = TRUE)
  ))) %>%
  pivot_longer(everything(), names_to = c("var", ".value"),
               names_pattern = "(.*)_(train_.*)") %>%
  left_join(
    valid %>% select(where(is.numeric)) %>%
      summarize(across(everything(), list(
        valid_mean = ~mean(., na.rm = TRUE),
        valid_p50  = ~quantile(., .5, na.rm = TRUE)
      ))) %>%
      pivot_longer(everything(), names_to = c("var", ".value"),
                   names_pattern = "(.*)_(valid_.*)"),
    by = "var"
  ) %>%
  mutate(
    mean_diff = train_mean - valid_mean,
    p50_diff  = train_p50 - valid_p50
  )

write_csv(qc_compare_num, file.path(dir_out, "train_valid_numeric_compare.csv"))

qc_compare_cat <- compare_cat(train, valid)
write_csv(qc_compare_cat, file.path(dir_out, "train_valid_categorical_compare.csv"))

# Illustrative range checks (adjust variable names as needed)
range_check <- list(
  age = c(0, 120),
  baseline_sbp = c(60, 260),
  baseline_dbp = c(30, 180),
  heart_rate   = c(30, 220),
  bmi          = c(10, 80)
)

range_flags <- purrr::map_dfr(names(range_check), function(v) {
  if (!v %in% names(train)) return(tibble(var = v, n_out = NA_integer_))
  rng <- range_check[[v]]
  tibble(
    var = v,
    n_out = sum(train[[v]] < rng[1] | train[[v]] > rng[2], na.rm = TRUE)
  )
})

write_csv(range_flags, file.path(dir_out, "range_outlier_train.csv"))

# ----------------------------- Save harmonized data --------------------------
# Saved locally only; data are not included in this public repository.
saveRDS(train,      file.path(dir_out, "train_raw.rds"))
saveRDS(valid,      file.path(dir_out, "valid_raw.rds"))
saveRDS(train_spec, file.path(dir_out, "train_spec_raw.rds"))
saveRDS(valid_spec, file.path(dir_out, "valid_spec_raw.rds"))
