# ============================================================
# 02_data_preprocess.R  (Methods-only / Representative script)
#
# Purpose:
#   Preprocessing workflow used in the CATIS ITE analysis:
#   - Align columns/types between training and validation sets
#   - Complete-case on the primary outcome only
#   - Train-only threshold learning for:
#       (a) "error -> NA" for selected variables (default 0.1%/99.9%)
#       (b) Winsorization for continuous predictors (default 0.5%/99.5%)
#   - No imputation for covariates (leave NA for GRF/MIA)
#   - Remove zero-variance predictors (recipes::step_zv)
#
# Notes:
#   - Individual-level data are not included in this repository.
#   - File paths below are placeholders; adapt for local use.
# ============================================================

source("scripts/00_setup.R")

suppressPackageStartupMessages({
  library(tidyverse)
  library(janitor)
  library(recipes)
  library(readr)
  library(purrr)
  library(stringr)
})

set.seed(42)

# ----------------------------- Paths (placeholders) --------------------------
# Prefer reading harmonized RDS produced by scripts/01_data_import.R (public repo writes to outputs/qc/)
DIR_QC    <- file.path("outputs", "qc")
DIR_CLEAN <- file.path("outputs", "clean")
DIR_CFG   <- file.path("outputs", "config")

dir.create(DIR_CLEAN, recursive = TRUE, showWarnings = FALSE)
dir.create(DIR_CFG,   recursive = TRUE, showWarnings = FALSE)

DEFAULT_FILES <- list(
  train = file.path(DIR_QC, "train_raw.rds"),
  valid = file.path(DIR_QC, "valid_raw.rds")
)

# Optional local fallback (not included in public repo)
# If you want local execution, put your data under data/raw and adapt names here.
FALLBACK_XLSX <- list(
  train = file.path("data", "raw", "train.xlsx"),
  valid = file.path("data", "raw", "valid.xlsx")
)

# ----------------------------- Variable roles -------------------------------
OUTCOME_MAIN      <- "primary_discharge"   # 1 = poor outcome, 0 = good
OUTCOME_SECONDARY <- "primary_m3"          # optional
TREATMENT         <- "treatment"           # 1 = early BP lowering, 0 = control

BINARY_VARS <- c(
  "male","hypertension","hyperlipidemia","dm","chd","af",
  "smoking","drinking","family_history_of_stroke","laa","ce","svo"
)

ID_CAND <- c("id","patient_id","subject_id","study_id")

# Exclude "special variables" (optional). Provide a text file with one variable name per line.
# Variable names should be in clean_names() format.
SPECIAL_VARS <- character(0)
SPECIAL_VARS_FILE <- file.path(DIR_CFG, "special_vars.txt")

# Variables more likely to have unit/entry issues (used for "error -> NA")
ERRORY_VARS <- c(
  "wbc","rbc","hb","hct","plt",
  "glucose","creatinine","e_gfr","uric_acid",
  "tc","tg","ldl_c","hdl_c",
  "pro_bnp","hs_crp","hcy",
  "aip","ty_g",
  "bmi","heart_rate"
)

# Default thresholds
USE_ERR2NA    <- TRUE
USE_WINSOR    <- TRUE
ERR_PROBS     <- c(0.001, 0.999)
WINSOR_PROBS  <- c(0.005, 0.995)

# ----------------------------- Helpers --------------------------------------
read_safely <- function(path_rds, path_xlsx = NULL) {
  if (!is.null(path_rds) && file.exists(path_rds)) {
    return(readRDS(path_rds))
  }
  if (!is.null(path_xlsx) && file.exists(path_xlsx)) {
    suppressPackageStartupMessages(library(readxl))
    return(suppressWarnings(readxl::read_xlsx(path_xlsx)))
  }
  stop("Input data not found. Provide local files or adapt paths in this script.")
}

read_special_vars <- function(file) {
  if (!file.exists(file)) return(character(0))
  x <- readr::read_lines(file, progress = FALSE)
  x <- x[nzchar(x)]
  janitor::make_clean_names(x)
}

is_binary_numeric <- function(x) is.numeric(x) && all(na.omit(unique(x)) %in% c(0,1))
is_binary01 <- function(x) {
  if (!is.numeric(x)) return(FALSE)
  ux <- unique(na.omit(as.numeric(x)))
  length(ux) <= 2 && all(sort(ux) %in% c(0,1))
}

# Normalize common yes/no strings to 0/1 (optional; conservative)
char_binary_dict <- list(
  yes = c("yes","y","是","有","true","t"),
  no  = c("no","n","否","无","false","f"),
  male = c("male","m","男"),
  female = c("female","f","女")
)

char_to_binary <- function(x) {
  if (!is.character(x)) return(x)
  xl <- tolower(x)
  dplyr::case_when(
    xl %in% char_binary_dict$yes ~ 1L,
    xl %in% char_binary_dict$no  ~ 0L,
    TRUE ~ NA_integer_
  )
}

gender_to_factor <- function(x) {
  if (!is.character(x)) return(x)
  xl <- tolower(x)
  dplyr::case_when(
    xl %in% char_binary_dict$male   ~ "male",
    xl %in% char_binary_dict$female ~ "female",
    TRUE ~ NA_character_
  ) |> factor(levels = c("male","female"))
}

assert_binary <- function(df, vars) {
  for (v in vars) if (v %in% names(df)) {
    df[[v]] <- if (is.factor(df[[v]])) as.integer(as.character(df[[v]])) else suppressWarnings(as.integer(df[[v]]))
    df[[v]][!df[[v]] %in% c(0,1)] <- NA_integer_
  }
  df
}

num_predictors <- function(df) {
  cand <- names(df)[vapply(df, is.numeric, TRUE)]
  setdiff(cand, c(OUTCOME_MAIN, OUTCOME_SECONDARY, TREATMENT))
}

# Error -> NA limits learned from training only
get_error_limits <- function(df, vars, probs) {
  vars <- intersect(vars, names(df))
  if (!length(vars)) return(tibble(variable=character(), lower=numeric(), upper=numeric()))
  tibble(
    variable = vars,
    lower = vapply(vars, function(v) quantile(df[[v]], probs[1], na.rm=TRUE, names=FALSE), numeric(1)),
    upper = vapply(vars, function(v) quantile(df[[v]], probs[2], na.rm=TRUE, names=FALSE), numeric(1))
  )
}

apply_error_to_na <- function(df, lims) {
  for (i in seq_len(nrow(lims))) {
    v  <- lims$variable[i]
    lo <- lims$lower[i]; hi <- lims$upper[i]
    if (v %in% names(df)) df[[v]][ df[[v]] < lo | df[[v]] > hi ] <- NA_real_
  }
  df
}

# Winsor limits learned from training only
get_winsor_limits <- function(df, probs) {
  vs <- num_predictors(df)
  vs <- vs[!vapply(df[vs], is_binary01, TRUE)]
  if (!length(vs)) return(tibble(variable=character(), lower=numeric(), upper=numeric()))
  tibble(
    variable = vs,
    lower = vapply(vs, function(v) quantile(df[[v]], probs[1], na.rm=TRUE, names=FALSE), numeric(1)),
    upper = vapply(vs, function(v) quantile(df[[v]], probs[2], na.rm=TRUE, names=FALSE), numeric(1))
  )
}

winsorize_by_limits <- function(x, lo, hi) {
  if (!is.numeric(x)) return(x)
  pmin(pmax(x, lo), hi)
}

apply_winsor <- function(df, lims) {
  for (i in seq_len(nrow(lims))) {
    v  <- lims$variable[i]
    lo <- lims$lower[i]; hi <- lims$upper[i]
    if (v %in% names(df)) df[[v]] <- winsorize_by_limits(df[[v]], lo, hi)
  }
  df
}

# Simple before/after summary written to outputs
compare_before_after <- function(df_before, df_after, tag, dir_out) {
  vars <- intersect(names(df_before), names(df_after))
  vars <- vars[vapply(df_before[vars], is.numeric, TRUE)]
  out <- purrr::map_dfr(vars, function(v) {
    b <- df_before[[v]]; a <- df_after[[v]]
    tibble(
      var = v,
      changed_prop = mean((!is.na(b) & !is.na(a)) & (b != a), na.rm = TRUE),
      set_na_prop  = mean(!is.na(b) & is.na(a)),
      mean_before  = mean(b, na.rm = TRUE),
      mean_after   = mean(a, na.rm = TRUE),
      p50_before   = quantile(b, 0.5, na.rm = TRUE, names = FALSE),
      p50_after    = quantile(a, 0.5, na.rm = TRUE, names = FALSE)
    )
  }) |> arrange(desc(changed_prop))
  write_csv(out, file.path(dir_out, paste0("qa_before_after_", tag, ".csv")))
}

# ----------------------------- Read ------------------------------------------
train0 <- read_safely(DEFAULT_FILES$train, FALLBACK_XLSX$train) %>% clean_names()
valid0 <- read_safely(DEFAULT_FILES$valid, FALLBACK_XLSX$valid) %>% clean_names()

# Identify ID column if present
id_col <- {
  tmp <- intersect(ID_CAND, names(train0))
  if (length(tmp) > 0) tmp[1] else NA_character_
}

# ----------------------------- Basic type handling ---------------------------
coerce_types <- function(df) {
  df %>%
    mutate(across(everything(), ~ if (is.character(.x)) stringr::str_squish(.x) else .x)) %>%
    mutate(across(where(is.character), ~ na_if(.x, ""))) %>%
    mutate(across(where(is_binary_numeric), as.integer))
}

train1 <- coerce_types(train0)
valid1 <- coerce_types(valid0)

# Preserve ID as character if present
if (!is.na(id_col)) {
  train1[[id_col]] <- as.character(train0[[id_col]])
  valid1[[id_col]] <- as.character(valid0[[id_col]])
}

# Conservative character conversion (optional): apply to non-ID character columns only
char_cols_train <- setdiff(names(train1)[vapply(train1, is.character, TRUE)], id_col)
char_cols_valid <- setdiff(names(valid1)[vapply(valid1, is.character, TRUE)], id_col)

train1 <- train1 %>%
  mutate(across(all_of(char_cols_train), char_to_binary)) %>%
  mutate(across(all_of(char_cols_train), gender_to_factor))

valid1 <- valid1 %>%
  mutate(across(all_of(char_cols_valid), char_to_binary)) %>%
  mutate(across(all_of(char_cols_valid), gender_to_factor))

# Enforce key 0/1 variables
force_binary <- c(OUTCOME_MAIN, OUTCOME_SECONDARY, TREATMENT, BINARY_VARS)
train1 <- assert_binary(train1, intersect(force_binary, names(train1)))
valid1 <- assert_binary(valid1, intersect(force_binary, names(valid1)))

# ----------------------------- Align columns/factors -------------------------
# (1) Keep columns present in training; add missing ones to validation as NA
all_vars <- names(train1)
valid2 <- valid1 %>%
  mutate(across(setdiff(all_vars, names(.)), ~NA)) %>%
  select(all_of(all_vars))

train2 <- train1

# (2) Align factor levels (train defines levels)
align_factors <- function(tr, va) {
  fac_vars <- names(tr)[purrr::map_lgl(tr, is.factor)]
  for (v in fac_vars) {
    lv <- levels(tr[[v]])
    va[[v]] <- factor(as.character(va[[v]]), levels = lv)
  }
  list(train = tr, valid = va)
}

aligned <- align_factors(train2, valid2)
train3 <- aligned$train
valid3 <- aligned$valid

# ----------------------------- Exclude "special variables" -------------------
special_from_file <- read_special_vars(SPECIAL_VARS_FILE)
SPECIAL_VARS_ALL <- unique(c(SPECIAL_VARS, special_from_file))
SPECIAL_VARS_ALL <- setdiff(SPECIAL_VARS_ALL, c(OUTCOME_MAIN, OUTCOME_SECONDARY, TREATMENT, id_col))

KEEP <- setdiff(names(train3), SPECIAL_VARS_ALL)
train3 <- train3 %>% select(all_of(KEEP))
valid3 <- valid3 %>% select(all_of(KEEP))

# ----------------------------- Complete-case on primary outcome --------------
train3 <- train3 %>% filter(!is.na(.data[[OUTCOME_MAIN]]))
valid3 <- valid3 %>% filter(!is.na(.data[[OUTCOME_MAIN]]))

# ----------------------------- Outlier handling (train-only learning) --------
train_before <- train3
valid_before <- valid3

if (USE_ERR2NA) {
  err_lims <- get_error_limits(train3, ERRORY_VARS, ERR_PROBS)
  write_csv(err_lims, file.path(DIR_CLEAN, "error_to_na_limits_train.csv"))
  train3 <- apply_error_to_na(train3, err_lims)
  valid3 <- apply_error_to_na(valid3, err_lims)
}

if (USE_WINSOR) {
  win_lims <- get_winsor_limits(train3, WINSOR_PROBS)
  write_csv(win_lims, file.path(DIR_CLEAN, "winsor_limits_train.csv"))
  train3 <- apply_winsor(train3, win_lims)
  valid3 <- apply_winsor(valid3, win_lims)
}

compare_before_after(train_before, train3, "train", DIR_CLEAN)
compare_before_after(valid_before, valid3, "valid", DIR_CLEAN)

# ----------------------------- Recipes: remove zero-variance only ------------
# No imputation / no scaling; leave NA for GRF/MIA
NON_PREDICTORS <- c(id_col, OUTCOME_SECONDARY)
NON_PREDICTORS <- NON_PREDICTORS[!is.na(NON_PREDICTORS)]

rec <- recipe(as.formula(paste(OUTCOME_MAIN, "~ .")),
              data = train3 %>% select(-any_of(NON_PREDICTORS))) %>%
  update_role(all_of(TREATMENT), new_role = "predictor") %>%
  step_zv(all_predictors())

prep_obj <- prep(rec, training = train3, retain = TRUE)

order_vec <- c(id_col, OUTCOME_MAIN, OUTCOME_SECONDARY, TREATMENT)
order_vec <- order_vec[!is.na(order_vec)]

train_clean <- bake(prep_obj, new_data = train3 %>% select(-any_of(NON_PREDICTORS))) %>%
  bind_cols(train3 %>% select(any_of(NON_PREDICTORS))) %>%
  relocate(all_of(order_vec))

valid_clean <- bake(prep_obj, new_data = valid3 %>% select(-any_of(NON_PREDICTORS))) %>%
  bind_cols(valid3 %>% select(any_of(NON_PREDICTORS))) %>%
  relocate(all_of(order_vec))

# ----------------------------- Save outputs ---------------------------------
saveRDS(train_clean, file.path(DIR_CLEAN, "train_clean.rds"))
saveRDS(valid_clean, file.path(DIR_CLEAN, "valid_clean.rds"))
saveRDS(prep_obj,    file.path(DIR_CLEAN, "preprocess_recipe.rds"))

# Minimal snapshot table (written to disk rather than verbose console output)
qc_tbl <- tibble(
  dataset = c("train", "valid"),
  n_rows  = c(nrow(train_clean), nrow(valid_clean)),
  n_cols  = c(ncol(train_clean), ncol(valid_clean))
)
write_csv(qc_tbl, file.path(DIR_CLEAN, "qc_snapshot.csv"))
