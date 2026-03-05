#!/usr/bin/env Rscript
# =============================================================================
# 08_patient_ite_waterfall.R
# Batch ITE Waterfall (SHAP-like via fastshap) for representative patients
# Convention: treated − control (negative = treatment benefit)
#
# Outputs:
#   figures/waterfall_patient_ite/waterfall_<timestamp>/{id}.pdf|png
#
# Requirements:
#   fastshap (>= 0.1.1), grf, ggplot2, dplyr, janitor, glue, grid
#
# Project root:
#   - Option A: set env var CATIS_ITE_DIR
#   - Option B: auto-detect by walking up from script location / working dir
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(janitor)
  library(stringr)
  library(glue)
  library(ggplot2)
  library(fastshap)
  library(grf)
  library(grid)   # for mm <-> data unit conversion
})

`%||%` <- function(x, y) if (is.null(x)) y else x

# --------------------------- Project root ---------------------------
get_proj_dir <- function() {
  env <- Sys.getenv("CATIS_ITE_DIR", unset = NA_character_)
  if (!is.na(env) && nzchar(env) && dir.exists(env)) {
    return(normalizePath(env, winslash = "/", mustWork = TRUE))
  }

  # Try script directory (when run via Rscript)
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  script_dir <- if (length(file_arg) > 0) dirname(sub("^--file=", "", file_arg)) else getwd()

  markers <- c("data", "models", "outputs", "figures", ".git")
  cur <- normalizePath(script_dir, winslash = "/", mustWork = TRUE)

  for (i in 1:10) {
    ok <- any(file.exists(file.path(cur, markers)))
    if (ok) return(cur)
    parent <- dirname(cur)
    if (identical(parent, cur)) break
    cur <- parent
  }
  normalizePath(getwd(), winslash = "/", mustWork = TRUE)
}

proj_dir <- get_proj_dir()

# ----------------------------- Paths -----------------------------
path_valid <- file.path(proj_dir, "data/clean/valid_clean.rds")
path_art   <- file.path(proj_dir, "models/cf_ensemble_primary_discharge.rds")
path_cand  <- file.path(proj_dir, "outputs/representative_patients_candidates.csv")

# Output (timestamped, no overwrite)
run_tag <- format(Sys.time(), "%Y%m%d_%H%M%S")
out_dir <- file.path(proj_dir, "figures", "waterfall_patient_ite", paste0("waterfall_", run_tag))
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# --------------------------- Parameters ---------------------------
TOP_POS <- 3L
TOP_NEG <- 3L
TOP_CAP <- 5L
NSIM    <- 128

# Sign convention: treated − control; negative = benefit
TAU_SIGN <- -1

# -------------------------- Pretty labels --------------------------
pretty_labels <- c(
  wc = "WC", pro_bnp = "Pro-BNP", uric_acid = "Uric acid", hcy = "Hcy",
  hs_crp = "hs-CRP", bmi = "BMI", ldl_c = "LDL-C", aip = "AIP",
  hdl_c = "HDL-C", tyg = "TyG", tg = "TG", tc = "TC",
  creatinine = "Creatinine", hct = "HCT", rbc = "RBC", plt = "PLT",
  glucose = "Glucose", hb = "Hb", wbc = "WBC", heart_rate = "Heart rate",
  baseline_nihss = "Baseline NIHSS", ota = "OTA", age = "Age",
  baseline_sbp = "Baseline SBP", baseline_dbp = "Baseline DBP",
  e_gfr = "eGFR", baseline_m_rs = "Baseline mRS", male = "Male",
  htn = "HTN", hld = "HLD", dm = "DM", chd = "CHD", af = "AF",
  smoking = "Smoking", drinking = "Drinking", fhs = "FHS",
  laa = "LAA", ce = "CE", svo = "SVO"
)

# ----------------------------- Utils ------------------------------
assert_exists <- function(p, what = "file") {
  if (!file.exists(p)) stop("Missing ", what, ": ", p, call. = FALSE)
}

fmt1 <- function(x) {
  if (is.na(x)) return("NA")
  if (is.numeric(x)) return(format(round(as.numeric(x), 1), nsmall = 1, trim = TRUE))
  as.character(x)
}

# Build numeric model matrix aligned to reference columns (x_cols)
mk_mm <- function(df, ref_cols, exclude = character()) {
  use_cols <- setdiff(names(df), exclude)
  X <- df[, use_cols, drop = FALSE]

  build_one <- function(v, nm) {
    if (is.logical(v)) v <- as.integer(v)
    else if (inherits(v, "Date") || inherits(v, "POSIXct") || inherits(v, "difftime")) v <- as.numeric(v)
    else if (is.list(v)) v <- as.character(v)

    if (is.character(v) || is.factor(v)) {
      f <- addNA(factor(v))
      mm1 <- model.matrix(~ 0 + f, data = data.frame(f = f), na.action = na.pass)
      colnames(mm1) <- paste0(nm, "_", sub("^f", "", colnames(mm1)))
      storage.mode(mm1) <- "double"
      return(mm1)
    } else {
      m <- matrix(as.numeric(v), ncol = 1L, dimnames = list(NULL, nm))
      storage.mode(m) <- "double"
      return(m)
    }
  }

  mats <- lapply(seq_along(use_cols), function(j) build_one(X[[use_cols[j]]], use_cols[j]))
  mm   <- do.call(cbind, mats)

  add <- setdiff(ref_cols, colnames(mm))
  if (length(add) > 0) {
    mm <- cbind(mm, matrix(0, nrow(mm), length(add), dimnames = list(NULL, add)))
  }
  mm <- mm[, ref_cols, drop = FALSE]
  storage.mode(mm) <- "double"
  mm
}

# Aggregate shap at "base variable" level when it is clearly one-hot of a categorical parent
aggregate_shap <- function(shap_vec, original_df) {
  nm <- names(shap_vec)
  bases <- vapply(nm, function(n) {
    if (n %in% names(original_df)) return(n)
    prefix <- sub("_[^_]+$", "", n)
    if (prefix %in% names(original_df)) {
      v <- original_df[[prefix]]
      if (is.character(v) || is.factor(v)) return(prefix)
    }
    n
  }, character(1))

  tibble(base = bases, shap = as.numeric(shap_vec)) %>%
    group_by(base) %>%
    summarise(contrib = sum(shap, na.rm = TRUE), .groups = "drop")
}

select_waterfall_terms <- function(contrib_tbl, top_pos = 3L, top_neg = 3L, top_cap = 5L) {
  pos_tbl <- contrib_tbl %>% filter(contrib > 0) %>% arrange(desc(contrib)) %>% slice_head(n = top_pos)
  neg_tbl <- contrib_tbl %>% filter(contrib < 0) %>% arrange(contrib)       %>% slice_head(n = top_neg)
  comb    <- bind_rows(pos_tbl, neg_tbl)

  if (nrow(comb) > top_cap) {
    comb <- comb %>% arrange(desc(abs(contrib))) %>% slice_head(n = top_cap)
  }
  comb %>% arrange(desc(contrib))
}

# ------------------------ Plotting settings ------------------------
COL_POS <- "#FF3B30"  # positive contribution (harm)
COL_NEG <- "#2F80FF"  # negative contribution (benefit)

BAR_HALF        <- 0.25
Y_EXPAND        <- c(0.08, 0.52)
LABEL_SIZE      <- 4.3
LABEL_OFFSET_MM <- 6

plot_waterfall <- function(contrib_tbl, f0, pred, id_label, bucket_label, value_labels) {
  ord <- c(setdiff(contrib_tbl$base, "All other variables"), "All other variables")

  dat <- contrib_tbl %>%
    mutate(base = factor(base, levels = ord)) %>%
    arrange(base) %>%
    mutate(
      direction = ifelse(contrib >= 0, "increase", "decrease"),
      cum_start = f0 + c(0, cumsum(head(contrib, -1))),
      cum_end   = f0 + cumsum(contrib),
      ymin      = pmin(cum_start, cum_end),
      ymax      = pmax(cum_start, cum_end),
      lab       = sprintf("%+.3f", contrib),
      var       = base
    )

  ord_display <- rev(ord)
  dat <- dat %>% mutate(var_display = factor(var, levels = ord_display))

  # Convert fixed mm offset to data units (always place labels to the right end of bar)
  rng <- range(c(dat$ymin, dat$ymax), finite = TRUE)
  p_tmp <- ggplot(dat, aes(x = var_display, y = ymax)) +
    geom_blank() +
    coord_flip(clip = "off") +
    scale_y_continuous(expand = expansion(mult = Y_EXPAND)) +
    theme_bw(base_size = 12) +
    theme(plot.margin = margin(24, 120, 16, 100))

  gt <- ggplot_gtable(ggplot_build(p_tmp))
  panel_mm <- sum(convertWidth(gt$widths[gt$layout$name == "panel"], "mm", valueOnly = TRUE))
  offset_data <- (LABEL_OFFSET_MM / panel_mm) * diff(rng)
  dat <- dat %>% mutate(y_lab = pmax(ymin, ymax) + offset_data)

  label_map <- setNames(value_labels, contrib_tbl$base)
  lab_fun <- function(x) unname(label_map[as.character(x)])

  ggplot(dat, aes(x = var_display, ymin = ymin, ymax = ymax, fill = direction)) +
    geom_rect(
      aes(
        xmin = as.numeric(var_display) - BAR_HALF,
        xmax = as.numeric(var_display) + BAR_HALF
      ),
      color = NA, alpha = 0.98
    ) +
    geom_hline(yintercept = 0, linetype = "dotted", color = "grey50", linewidth = 0.6) +
    geom_text(aes(x = var_display, y = y_lab, label = lab),
              size = LABEL_SIZE, inherit.aes = FALSE) +
    coord_flip(clip = "off") +
    scale_x_discrete(limits = ord_display, labels = lab_fun) +
    scale_fill_manual(values = c(increase = COL_POS, decrease = COL_NEG)) +
    scale_y_continuous(expand = expansion(mult = Y_EXPAND)) +
    labs(
      title    = str_wrap(glue("Patient {id_label} — ITE Waterfall (bucket: {bucket_label})"), width = 60),
      subtitle = glue("predicted ITE = {sprintf('%.3f', pred)}"),
      x = NULL,
      y = expression(atop(bold("Contribution to ITE"), italic("Benefit treat   |   Benefit control")))
    ) +
    theme_bw(base_size = 12) +
    theme(
      legend.position = "none",
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank(),
      axis.text.y = element_text(size = 10, color = "black"),
      plot.margin = margin(24, 120, 16, 100),
      plot.title = element_text(hjust = 0.0, face = "bold")
    )
}

# ------------------------- Load data/model -------------------------
assert_exists(path_valid, "valid_clean.rds")
assert_exists(path_art,   "model artifact rds")
assert_exists(path_cand,  "representative candidates csv")

valid <- readRDS(path_valid) %>% clean_names()
if (!"id" %in% names(valid)) valid$id <- seq_len(nrow(valid))
valid$id <- as.character(valid$id)

cand <- read_csv(path_cand, show_col_types = FALSE) %>% clean_names()
if (!"id" %in% names(cand)) stop("Candidates file must contain an 'id' column.", call. = FALSE)
cand$id <- as.character(cand$id)
if (!"bucket" %in% names(cand)) cand$bucket <- ""

artifact <- readRDS(path_art)
models   <- artifact$models
x_cols   <- artifact$x_cols
scale_a  <- artifact$scale_a
shift_b  <- artifact$shift_b

if (is.null(models) || length(models) == 0) stop("artifact$models is empty.", call. = FALSE)
if (is.null(x_cols) || length(x_cols) == 0) stop("artifact$x_cols is missing.", call. = FALSE)
if (is.null(scale_a) || is.null(shift_b)) stop("artifact scale/shift is missing.", call. = FALSE)

# Build X matrix aligned to x_cols
exclude_cols <- c("treatment", "primary_discharge", "primary_m3", "id",
                  "patient_id", "subject_id", "study_id", "row_id")
X_valid <- mk_mm(valid, ref_cols = x_cols, exclude = exclude_cols)

# --------------------------- Prediction ---------------------------
predict_tau_ensemble <- function(newX) {
  newX <- as.matrix(newX[, x_cols, drop = FALSE])
  preds <- vapply(models, function(m) stats::predict(m, newX)$predictions, numeric(nrow(newX)))
  if (is.null(dim(preds))) preds <- matrix(preds, ncol = 1)
  rowMeans(preds)
}

calibrate <- function(tau) as.numeric(shift_b + scale_a * tau)

# Black-box wrapper for fastshap
f_blackbox <- function(object, newdata) {
  newX <- as.matrix(newdata[, x_cols, drop = FALSE])
  tau_raw <- calibrate(predict_tau_ensemble(newX))
  TAU_SIGN * tau_raw
}

# Global baseline (for waterfall accumulation start)
tau_all <- f_blackbox(NULL, as.data.frame(X_valid[, x_cols, drop = FALSE]))
f0_global <- mean(tau_all, na.rm = TRUE)

message(glue("Project root: {proj_dir}"))
message(glue("Candidates: {nrow(cand)}; output dir: {out_dir}"))
message("Start generating waterfall plots...")

# ------------------------------ Loop -----------------------------
for (i in seq_len(nrow(cand))) {
  id_i     <- cand$id[i]
  bucket_i <- cand$bucket[i] %||% ""

  idx <- which(valid$id == id_i)
  if (length(idx) != 1L) {
    warning("Skip id=", id_i, " (not found or duplicated in valid).")
    next
  }

  x_row <- X_valid[idx, , drop = FALSE]
  pred_tc <- f_blackbox(NULL, as.data.frame(x_row[, x_cols, drop = FALSE]))[1]

  # Make randomness reproducible but patient-specific
  set.seed(1000 + i)

  phi <- fastshap::explain(
    object       = NULL,
    X            = as.data.frame(X_valid[, x_cols, drop = FALSE]),
    pred_wrapper = f_blackbox,
    newdata      = as.data.frame(x_row[, x_cols, drop = FALSE]),
    nsim         = NSIM
  )

  shap_vec <- as.numeric(phi[1, ])
  names(shap_vec) <- colnames(phi)

  agg_tbl <- aggregate_shap(shap_vec, original_df = valid)

  wf_terms_core <- select_waterfall_terms(agg_tbl, TOP_POS, TOP_NEG, TOP_CAP)
  remainder <- sum(agg_tbl$contrib, na.rm = TRUE) - sum(wf_terms_core$contrib, na.rm = TRUE)

  wf_terms <- bind_rows(
    wf_terms_core,
    tibble(base = "All other variables", contrib = remainder)
  )

  # Build y-axis labels: "PrettyName value (median)" (NA explicitly shown)
  value_labels <- vapply(wf_terms$base, function(v) {
    if (v == "All other variables") return(v)

    pretty <- pretty_labels[[v]] %||% v
    if (v %in% names(valid)) {
      if (is.numeric(valid[[v]])) {
        raw_val <- suppressWarnings(as.numeric(valid[[v]][idx]))
        val_s <- ifelse(is.na(raw_val), "NA", fmt1(raw_val))
        med_s <- fmt1(stats::median(valid[[v]], na.rm = TRUE))
      } else {
        vv <- as.character(valid[[v]][idx])
        val_s <- ifelse(is.na(vv) || vv == "", "NA", vv)
        med_s <- NA_character_
      }
    } else {
      val_s <- "NA"; med_s <- NA_character_
    }

    parts <- c(val_s, if (!is.na(med_s)) paste0("(", med_s, ")") else NULL)
    lbl <- trimws(paste(pretty, paste(parts, collapse = " ")))
    if (lbl == "") pretty else lbl
  }, character(1))

  p <- plot_waterfall(
    contrib_tbl   = wf_terms,
    f0            = f0_global,
    pred          = pred_tc,
    id_label      = id_i,
    bucket_label  = bucket_i,
    value_labels  = value_labels
  )

  fp_pdf <- file.path(out_dir, glue("{id_i}.pdf"))
  fp_png <- file.path(out_dir, glue("{id_i}.png"))

  ggsave(fp_pdf, p, width = 7.0, height = 4.8)
  ggsave(fp_png, p, width = 7.0, height = 4.8, dpi = 600)

  message(glue("✓ Exported: {basename(fp_pdf)} / {basename(fp_png)} (id={id_i}, ITE={sprintf('%.3f', pred_tc)})"))
}

message(glue("Done! Output directory: {out_dir}"))
