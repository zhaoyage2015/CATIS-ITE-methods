suppressPackageStartupMessages({
  library(tidyverse)
  library(janitor)
  library(gtsummary)
  library(readr)
  library(officer)
  library(flextable)
})

# ------------------------------------------------------------------------------
# 05_baseline_by_cate_quintiles.R
#
# Purpose:
#   Generate a baseline characteristics table stratified by predicted
#   individualized treatment-effect (CATE) quintiles in the validation cohort.
#   The table is exported in HTML, CSV, and editable Word (.docx) formats.
#
# Required inputs:
#   - data/clean/valid_clean.rds
#   - outputs/cate_valid_primary_discharge.csv
#     (first column = subject identifier; includes tau_hat)
#
# Notes:
#   - This repository provides representative scripts for method transparency.
#   - Individual-level trial data are not included.
# ------------------------------------------------------------------------------

# ---------------- Parameters ----------------
ntiles     <- 5
path_valid <- "data/clean/valid_clean.rds"
path_tau   <- "outputs/cate_valid_primary_discharge.csv"

out_html <- "outputs/baseline_by_cate_quintiles.html"
out_csv  <- "outputs/baseline_by_cate_quintiles.csv"
out_docx <- "outputs/baseline_by_cate_quintiles.docx"

dir.create(dirname(out_html), showWarnings = FALSE, recursive = TRUE)
dir.create(dirname(out_csv),  showWarnings = FALSE, recursive = TRUE)
dir.create(dirname(out_docx), showWarnings = FALSE, recursive = TRUE)

# ---------------- Read data ----------------
valid <- readr::read_rds(path_valid) %>% clean_names()
tau   <- readr::read_csv(path_tau, show_col_types = FALSE) %>% clean_names()

# ---------------- Identify / create ID column ----------------
id_col <- intersect(
  c("id", "patient_id", "subject_id", "study_id", "row_id"),
  names(valid)
)[1]

if (is.na(id_col)) {
  id_col <- "row_id"
  valid[[id_col]] <- seq_len(nrow(valid))
}

# Treat first column of tau file as ID
tau <- tau %>%
  rename(!!id_col := names(.)[1]) %>%
  mutate(!!id_col := as.character(.data[[id_col]]))

# ---------------- Merge and define quintiles ----------------
df <- valid %>%
  mutate(!!id_col := as.character(.data[[id_col]])) %>%
  left_join(tau, by = id_col) %>%
  filter(!is.na(tau_hat)) %>%
  mutate(
    quintile = ntile(tau_hat, ntiles),
    quintile = factor(
      quintile,
      levels = seq_len(ntiles),
      labels = paste0("Q", seq_len(ntiles))
    )
  )

# ---------------- Automatic variable selection ----------------
exclude_vars <- c(
  id_col, "treatment", "primary_discharge",
  "primary_m3", "tau_hat", "quintile"
)

candidates <- setdiff(names(df), exclude_vars)

is_continuous <- function(x) {
  is.numeric(x) && dplyr::n_distinct(x, na.rm = TRUE) > 10
}

cont_vars <- candidates[vapply(df[candidates], is_continuous, logical(1))]
cat_vars  <- setdiff(candidates, cont_vars)

# Convert common binary indicators (0/1) and character variables
bin_numeric <- cat_vars[vapply(
  df[cat_vars],
  function(x) is.numeric(x) && all(na.omit(x) %in% c(0, 1)),
  logical(1)
)]

char_vars <- cat_vars[vapply(df[cat_vars], is.character, logical(1))]

df <- df %>%
  mutate(across(
    all_of(bin_numeric),
    ~ factor(.x, levels = c(0, 1), labels = c("No", "Yes"))
  )) %>%
  mutate(across(all_of(char_vars), factor))

# ---------------- Build baseline table by quintile ----------------
tbl_base <-
  gtsummary::tbl_summary(
    data = df,
    by   = quintile,
    type = list(
      all_continuous() ~ "continuous2",
      all_categorical() ~ "categorical"
    ),
    statistic = list(
      all_continuous() ~ "{median} ({p25}, {p75})",
      all_categorical() ~ "{n} ({p}%)"
    ),
    digits  = list(all_continuous() ~ 2),
    missing = "no",
    include = c(all_of(cont_vars), all_of(cat_vars))
  ) %>%
  add_overall(last = FALSE) %>%
  add_p(
    test = list(
      all_continuous() ~ "kruskal.test",
      all_categorical() ~ "chisq.test"
    ),
    pvalue_fun = ~ style_pvalue(.x, digits = 3)
  ) %>%
  bold_labels()

# ---------------- Export: HTML + CSV ----------------
gtsummary::as_gt(tbl_base) %>%
  gtsummary::gtsave(out_html)

readr::write_csv(as_tibble(tbl_base), out_csv)

# ---------------- Export: Editable Word (.docx) ----------------
ft <- gtsummary::as_flex_table(tbl_base) %>%
  flextable::set_table_properties(width = 1, layout = "autofit") %>%
  flextable::fontsize(size = 10, part = "all") %>%
  flextable::align(align = "center", part = "all") %>%
  flextable::bold(part = "header") %>%
  flextable::autofit()

doc <- officer::read_docx() %>%
  officer::body_add_par(
    "Baseline characteristics stratified by predicted treatment-effect quintiles",
    style = "heading 1"
  ) %>%
  flextable::body_add_flextable(ft) %>%
  officer::body_add_par("", style = "Normal")

print(doc, target = out_docx)

message("Saved: ", out_html)
message("Saved: ", out_csv)
message("Saved: ", out_docx)
