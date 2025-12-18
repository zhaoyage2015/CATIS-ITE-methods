# ============================================================
# 03_model_causal_forest.R  (Methods-only / Representative script)
#
# Purpose:
#   Representative R script illustrating the main workflow for
#   causal forest–based individualized treatment effect (ITE) modeling:
#     - build design matrix (incl. one-hot encoding)
#     - fit an ensemble of causal forests (honest GRF)
#     - pooled ATE and best linear projection (BLP)
#     - calibration using DR pseudo-outcome (global scale/offset)
#     - validation: DR risk difference by CATE quantiles + Qini (DR gain)
#     - variable importance: pretty label mapping + aggregation across one-hot
#
# Notes:
#   - Individual-level CATIS trial data are not included in this repository.
#   - This script assumes preprocessed datasets exist under data/clean/
#     (e.g., train_clean.rds and valid_clean.rds) with harmonized columns.
#   - The code is provided for methodological transparency.
# ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(readr)
  library(stringr)
  library(glue)
  library(grf)
  library(janitor)
  library(forcats)
})

`%||%` <- function(x, y) if (is.null(x)) y else x

# --------------------------- Configuration (illustrative) ----------------------
set.seed(20251020)

CFG <- list(
  var_treat     = "treatment",
  var_outcome   = "primary_discharge",
  id_candidates = c("id","patient_id","subject_id","study_id","row_id"),
  n_runs      = 5L,
  n_trees     = 2000L,
  sample_frac = 0.632,
  min_node    = 20L,
  mtry_rule   = function(p) max(1L, min(p, round(1.5 * sqrt(p)))),
  out_prefix  = "primary_discharge",
  quantiles_q = 5L,
  min_per_arm = 30L,
  rare_min_prop = 0.01
)

# ----------------------------- Paths (relative) --------------------------------
proj_dir       <- getwd()
path_clean     <- file.path(proj_dir, "data", "clean")
out_dir_models <- file.path(proj_dir, "models")
out_dir_fig    <- file.path(proj_dir, "figures")
out_dir_tbl    <- file.path(proj_dir, "outputs")

invisible(dir.create(out_dir_models, recursive = TRUE, showWarnings = FALSE))
invisible(dir.create(out_dir_fig,    recursive = TRUE, showWarnings = FALSE))
invisible(dir.create(out_dir_tbl,    recursive = TRUE, showWarnings = FALSE))

# ----------------------------- Utilities ---------------------------------------
try_load_dataset <- function(which = c("train","valid")) {
  which <- match.arg(which)
  candidates <- c(
    file.path(path_clean, glue("{which}_clean.rds")),
    file.path(path_clean, glue("{which}_clean.csv"))
  )
  found <- candidates[file.exists(candidates)]
  if (length(found) == 0L) {
    stop(glue("No {which} dataset found under data/clean/. Expected {which}_clean.rds or {which}_clean.csv."))
  }
  f <- found[1]
  message(glue("Read file: {f}"))
  ext <- tools::file_ext(f)
  if (ext == "rds") return(readRDS(f))
  if (ext == "csv") return(readr::read_csv(f, show_col_types = FALSE))
  stop("Unsupported file extension.")
}

as_binary01 <- function(x) {
  if (is.logical(x)) return(as.integer(x))
  if (is.numeric(x)) {
    ux <- sort(unique(x[!is.na(x)]))
    if (all(ux %in% c(0,1))) return(as.numeric(x))
  }
  as.integer(factor(x)) - 1L
}

collapse_rare_levels <- function(x, min_prop = 0.01, other_label = "other") {
  if (!(is.character(x) || is.factor(x))) return(x)
  f <- factor(x, exclude = NULL)
  tab <- prop.table(table(f, useNA = "ifany"))
  rare_lvls <- names(tab)[tab < min_prop & names(tab) != "(Missing)"]
  if (length(rare_lvls) == 0) return(f)
  f2 <- f; f2[f %in% rare_lvls] <- other_label
  factor(f2, exclude = NULL)
}

mk_mm <- function(df, ref_cols = NULL, exclude = character(),
                  collapse_cats = TRUE, min_prop = 0.01, drop_nzv = TRUE) {
  use_cols <- setdiff(names(df), exclude)
  X <- df[, use_cols, drop = FALSE]

  build_one <- function(v, nm) {
    if (is.logical(v)) v <- as.integer(v)
    else if (inherits(v, "Date") || inherits(v, "POSIXct") || inherits(v, "difftime")) v <- as.numeric(v)
    else if (is.list(v)) v <- as.character(v)

    if (is.character(v) || is.factor(v)) {
      if (collapse_cats) v <- collapse_rare_levels(v, min_prop = min_prop)
      f <- addNA(factor(v))
      mm1 <- model.matrix(~ 0 + f, data = data.frame(f = f), na.action = na.pass)
      colnames(mm1) <- paste0(nm, "_", sub("^f", "", colnames(mm1)))
      storage.mode(mm1) <- "double"
      return(mm1)
    } else {
      vec <- as.numeric(v)
      m <- matrix(vec, ncol = 1L, dimnames = list(NULL, nm))
      storage.mode(m) <- "double"
      return(m)
    }
  }

  mats <- lapply(seq_along(use_cols), function(j) build_one(X[[use_cols[j]]], use_cols[j]))
  mm <- do.call(cbind, mats)

  if (!is.null(ref_cols)) {
    add <- setdiff(ref_cols, colnames(mm))
    if (length(add) > 0) mm <- cbind(mm, matrix(0, nrow(mm), length(add), dimnames = list(NULL, add)))
    mm <- mm[, ref_cols, drop = FALSE]
  }

  if (drop_nzv && is.null(ref_cols)) {
    sdv <- apply(mm, 2, stats::sd, na.rm = TRUE)
    keep <- is.finite(sdv) & sdv > 1e-8
    if (sum(!keep) > 0) {
      mm <- mm[, keep, drop = FALSE]
      message(glue("Dropped {sum(!keep)} near-zero-variance columns."))
    }
  }

  storage.mode(mm) <- "double"
  mm
}

save_table <- function(df, name) {
  fp <- file.path(out_dir_tbl, name)
  readr::write_csv(df, fp)
  message(glue("Saved table: {fp}"))
}

save_plot <- function(p, name, w = 8, h = 5) {
  fp <- file.path(out_dir_fig, name)
  ggplot2::ggsave(fp, plot = p, width = w, height = h, dpi = 300)
  message(glue("Saved figure: {fp}"))
}

eval_by_quantiles_dr <- function(df_valid, tau_hat, outcome, treat, q = 5,
                                 p_trt = NULL, min_per_arm = 30) {
  stopifnot(length(tau_hat) == nrow(df_valid))
  brks <- stats::quantile(tau_hat, probs = seq(0, 1, length.out = q + 1), na.rm = TRUE)
  bins <- if (any(duplicated(brks))) dplyr::ntile(tau_hat, q) else cut(
    tau_hat, breaks = brks, include.lowest = TRUE, labels = FALSE)
  bins <- as.integer(bins)
  dv <- dplyr::mutate(df_valid, bin = bins)

  y <- as.numeric(dv[[outcome]])
  w <- as.numeric(dv[[treat]])
  if (is.null(p_trt)) p_trt <- mean(w == 1, na.rm = TRUE)
  pc <- 1 - p_trt

  Z <- y * (w / p_trt) - y * ((1 - w) / pc)

  tab <- dv |> dplyr::group_by(bin) |>
    dplyr::summarise(n=n(), n_trt=sum(w==1,na.rm=TRUE), n_ctl=sum(w==0,na.rm=TRUE), .groups="drop")

  if (any(tab$n_trt < min_per_arm | tab$n_ctl < min_per_arm)) {
    new_bin <- bins
    for (b in sort(unique(bins))) {
      idx <- which(new_bin == b); if (length(idx) == 0) next
      n_trt_b <- sum(w[idx] == 1, na.rm = TRUE)
      n_ctl_b <- sum(w[idx] == 0, na.rm = TRUE)
      if (n_trt_b < min_per_arm || n_ctl_b < min_per_arm) {
        new_bin[idx] <- if (b > 1) b - 1L else b + 1L
      }
    }
    dv$bin <- new_bin
  }

  dv |>
    dplyr::group_by(bin) |>
    dplyr::summarise(
      n = dplyr::n(),
      n_trt = sum(w==1,na.rm = TRUE),
      n_ctl = sum(w==0,na.rm = TRUE),
      rd_dr = mean(Z, na.rm = TRUE),
      se = stats::sd(Z, na.rm = TRUE)/sqrt(n),
      ci_l = rd_dr - 1.96*se,
      ci_u = rd_dr + 1.96*se,
      tau_mean = mean(tau_hat[bin == unique(bin)], na.rm = TRUE),
      .groups="drop"
    ) |>
    dplyr::arrange(bin)
}

qini_dr_curve <- function(df_valid, tau_hat, outcome, treat, p_trt = NULL) {
  y <- as.numeric(df_valid[[outcome]])
  w <- as.numeric(df_valid[[treat]])
  if (is.null(p_trt)) p_trt <- mean(w == 1, na.rm = TRUE)
  pc <- 1 - p_trt

  Z <- y * (w / p_trt) - y * ((1 - w) / pc)
  ord <- order(tau_hat, decreasing = TRUE)
  Z_ord <- Z[ord]

  tibble(frac = seq_along(Z_ord) / length(Z_ord),
         cum_gain = cumsum(Z_ord))
}

# ----------------------------- Load datasets -----------------------------------
train <- try_load_dataset("train") |> clean_names()
valid <- try_load_dataset("valid") |> clean_names()

var_treat   <- CFG$var_treat
var_outcome <- CFG$var_outcome

if (!all(c(var_treat, var_outcome) %in% names(train))) stop("Train missing key columns.")
if (!all(c(var_treat, var_outcome) %in% names(valid))) stop("Valid missing key columns.")

id_col <- { hit <- intersect(CFG$id_candidates, names(train)); if (length(hit) > 0) hit[1] else "row_id" }
if (!id_col %in% names(train)) train[[id_col]] <- seq_len(nrow(train))
if (!id_col %in% names(valid)) valid[[id_col]] <- seq_len(nrow(valid))

train[[var_outcome]] <- as_binary01(train[[var_outcome]])
valid[[var_outcome]] <- as_binary01(valid[[var_outcome]])
train[[var_treat]]   <- as_binary01(train[[var_treat]])
valid[[var_treat]]   <- as_binary01(valid[[var_treat]])

# ----------------------------- Build X / Y / W --------------------------------
covariate_exclude <- c(var_treat, var_outcome, "primary_m3", CFG$id_candidates)

X_train <- mk_mm(train, exclude = covariate_exclude,
                 collapse_cats = TRUE, min_prop = CFG$rare_min_prop, drop_nzv = TRUE)

X_valid <- mk_mm(valid, ref_cols = colnames(X_train), exclude = covariate_exclude,
                 collapse_cats = TRUE, min_prop = CFG$rare_min_prop, drop_nzv = FALSE)

if (ncol(X_train) == 0) stop("X_train has no covariates.")

W_train <- as.numeric(train[[var_treat]])
Y_train <- as.numeric(train[[var_outcome]])

# ----------------------------- Train ensemble CF -------------------------------
num_threads <- max(1L, parallel::detectCores(logical = TRUE) - 1L)
p <- ncol(X_train)

models <- vector("list", CFG$n_runs)
vi_mat <- matrix(NA_real_, nrow = p, ncol = CFG$n_runs, dimnames = list(colnames(X_train), NULL))

pt_train <- mean(W_train == 1)

rf_y <- regression_forest(
  X_train, Y_train,
  num.trees = 1000,
  honesty = TRUE,
  num.threads = num_threads
)
Y_hat_train <- predict(rf_y)$predictions

for (i in seq_len(CFG$n_runs)) {
  set.seed(20251020 + i)

  cf_i <- causal_forest(
    X = X_train, Y = Y_train, W = W_train,
    Y.hat = Y_hat_train,
    W.hat = rep(pt_train, length(W_train)),
    num.trees        = CFG$n_trees,
    mtry             = CFG$mtry_rule(ncol(X_train)),
    min.node.size    = CFG$min_node,
    sample.fraction  = CFG$sample_frac,
    honesty          = TRUE,
    honesty.fraction = 0.5,
    stabilize.splits = TRUE,
    tune.parameters  = "all",
    tune.num.trees   = 200,
    num.threads      = num_threads
  )

  models[[i]] <- cf_i
  vi_mat[, i] <- variable_importance(cf_i)
}

predict_tau_ensemble <- function(newX) {
  rowMeans(sapply(models, function(m) predict(m, newX)$predictions))
}

# ----------------------------- Pooled ATE --------------------------------------
ate_runs <- lapply(models, function(m) average_treatment_effect(m, target.sample = "all"))
est_vec <- sapply(ate_runs, function(x) unname(x[["estimate"]]))
se_vec  <- sapply(ate_runs, function(x) unname(x[["std.err"]]))

pooled_var <- mean(se_vec^2, na.rm = TRUE) + stats::var(est_vec, na.rm = TRUE)
pooled_se  <- sqrt(pooled_var)
pooled_est <- mean(est_vec, na.rm = TRUE)

ate_tbl <- tibble(outcome  = CFG$out_prefix, estimate = pooled_est, std.err = pooled_se) |>
  mutate(ci.lower = estimate - 1.96 * std.err,
         ci.upper = estimate + 1.96 * std.err)

save_table(ate_tbl, glue("ate_{CFG$out_prefix}.csv"))

# ----------------------------- BLP (pooled) ------------------------------------
extract_blp_tbl <- function(bx, X) {
  if (is.list(bx) && !is.null(bx$coef.estimates)) {
    mat <- as.matrix(bx$coef.estimates)
    term <- rownames(mat); est <- as.numeric(mat[,1])
    se   <- if (ncol(mat) >= 2) as.numeric(mat[,2]) else rep(NA_real_, nrow(mat))
    if (is.null(term)) term <- colnames(X)[seq_len(nrow(mat))]
    return(tibble(term=term, estimate=est, std.err=se))
  }
  if (is.list(bx) && !is.null(bx$beta.hat)) {
    est <- as.numeric(bx$beta.hat)
    nm  <- names(bx$beta.hat); if (is.null(nm)) nm <- colnames(X)[seq_along(est)]
    se  <- if (!is.null(bx$beta.se)) as.numeric(bx$beta.se) else rep(NA_real_, length(est))
    if (length(se) != length(est)) se <- rep(NA_real_, length(est))
    return(tibble(term=nm, estimate=est, std.err=se))
  }
  tibble(term = colnames(X), estimate = NA_real_, std.err = NA_real_)
}

get_blp_tbl <- function(m, X) {
  bx <- tryCatch(best_linear_projection(m, X), error = function(e) e)
  if (inherits(bx, "error")) {
    warning("best_linear_projection failed: ", bx$message)
    return(tibble(term = colnames(X), estimate = NA_real_, std.err = NA_real_))
  }
  if (inherits(bx, "coeftest")) bx <- list(coef.estimates = bx)
  extract_blp_tbl(bx, X)
}

blp_runs <- lapply(models, function(m) get_blp_tbl(m, X_train))
nr <- sapply(blp_runs, nrow)
if (length(unique(nr)) != 1L) stop(glue("BLP row size differs across runs: {paste(nr, collapse=', ')}"))

blp_long <- purrr::imap_dfr(blp_runs, ~ dplyr::mutate(.x, run = .y))

blp_pooled <- blp_long |>
  dplyr::group_by(term) |>
  dplyr::summarise(
    est_mean  = mean(estimate, na.rm = TRUE),
    pooled_se = sqrt(mean(std.err^2, na.rm = TRUE) + stats::var(estimate, na.rm = TRUE)),
    .groups = "drop"
  ) |>
  dplyr::mutate(z = est_mean / pooled_se, p = 2 * pnorm(-abs(z))) |>
  dplyr::arrange(p)

save_table(
  blp_pooled |>
    dplyr::transmute(term,
                     estimate = round(est_mean, 6),
                     std.err  = round(pooled_se, 6),
                     z        = round(z, 3),
                     p        = signif(p, 3)),
  glue("blp_{CFG$out_prefix}.csv")
)

# ----------------------------- Variable importance: map + aggregate ------------
pretty_labels <- c(
  wc = "WC", pro_bnp = "Pro-BNP", uric_acid = "Uric acid", hcy = "Hcy",
  hs_crp = "hs-CRP", bmi = "BMI", ldl_c = "LDL-C", aip = "AIP", hdl_c = "HDL-C",
  tyg = "TyG", tg = "TG", tc = "TC", creatinine = "Creatinine", hct = "HCT",
  rbc = "RBC", plt = "PLT", glucose = "Glucose", hb = "Hb", wbc = "WBC",
  heart_rate = "Heart rate", baseline_nihss = "Baseline NIHSS", ota = "OTA",
  age = "Age", baseline_sbp = "Baseline SBP", baseline_dbp = "Baseline DBP",
  e_gfr = "eGFR", baseline_m_rs = "Baseline mRS", male = "Male", htn = "HTN",
  hld = "HLD", dm = "DM", chd = "CHD", af = "AF", smoking = "Smoking",
  drinking = "Drinking", fhs = "FHS", laa = "LAA", ce = "CE", svo = "SVO",
  ota_time = "OTA time"
)

vi_avg <- rowMeans(vi_mat, na.rm = TRUE)
vi_raw_tbl <- tibble(variable_raw = colnames(X_train), importance = as.numeric(vi_avg))

strip_to_base_if_categorical <- function(v) {
  base_try <- sub("_[^_]+$", "", v)
  if (base_try %in% names(pretty_labels) && !(v %in% names(pretty_labels))) base_try else v
}

vi_grouped <- vi_raw_tbl |>
  mutate(base = vapply(variable_raw, strip_to_base_if_categorical, character(1))) |>
  group_by(base) |>
  summarise(importance = sum(importance, na.rm = TRUE), .groups = "drop") |>
  mutate(display = unname(pretty_labels[base]) %||% base) |>
  arrange(desc(importance)) |>
  mutate(
    importance_norm = 100 * importance / max(importance),
    cum_importance  = cumsum(importance_norm) / sum(importance_norm) * 100,
    rank = dense_rank(desc(importance))
  )

save_table(vi_raw_tbl |> arrange(desc(importance)),
           glue("vi_full_rawcols_{CFG$out_prefix}.csv"))
save_table(vi_grouped |> select(rank, variable = display, importance,
                                importance_norm, cum_importance),
           glue("vi_full_{CFG$out_prefix}.csv"))

p_vi_top20 <- vi_grouped |>
  slice_head(n = 20) |>
  mutate(variable = fct_reorder(display, importance_norm)) |>
  ggplot(aes(x = variable, y = importance_norm)) +
  geom_col() +
  coord_flip() +
  geom_text(aes(label = sprintf("%.1f%%", importance_norm)),
            hjust = -0.1, size = 3) +
  labs(title = "Top 20 Variable Importances (Causal Forest)",
       subtitle = glue("{CFG$out_prefix} — normalized"),
       x = NULL, y = "Importance (%)") +
  theme_minimal(base_size = 12) +
  theme(panel.grid.major.y = element_blank())

save_plot(p_vi_top20, glue("vi_top20_{CFG$out_prefix}.png"))

# ----------------------------- Calibration (DR pseudo-outcome) -----------------
pt <- mean(W_train == 1); pc <- 1 - pt
tau_train <- rowMeans(sapply(models, function(m) predict(m, X_train)$predictions))
Z_train <- Y_train * (W_train/pt) - Y_train * ((1 - W_train)/pc)

cal <- lm(Z_train ~ tau_train)
scale_a <- unname(coef(cal)[2])
shift_b <- unname(coef(cal)[1])
calibrate_tau <- function(tau) as.numeric(shift_b + scale_a * tau)

# ----------------------------- Validation: CATE + DR by quantiles --------------
tau_valid_raw <- predict_tau_ensemble(X_valid)
tau_valid_cal <- calibrate_tau(tau_valid_raw)

tau_tbl <- tibble(!!id_col := valid[[id_col]], tau_hat = as.numeric(tau_valid_cal))
save_table(tau_tbl, glue("cate_valid_{CFG$out_prefix}.csv"))

qres <- eval_by_quantiles_dr(
  valid, tau_valid_cal,
  outcome = var_outcome, treat = var_treat,
  q = CFG$quantiles_q,
  p_trt = mean(W_train == 1),
  min_per_arm = CFG$min_per_arm
)
save_table(qres, glue("valid_quantiles_dr_{CFG$out_prefix}.csv"))

p_dec <- qres |>
  mutate(bin = factor(bin, levels = bin)) |>
  ggplot(aes(x = bin, y = rd_dr, group = 1)) +
  geom_line() + geom_point() +
  geom_errorbar(aes(ymin = ci_l, ymax = ci_u), width = 0.15) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(title = "Validation — DR risk difference by predicted CATE quantiles",
       x = "Quantile of predicted CATE (low → high)",
       y = "Risk difference (treated − control, DR)") +
  theme_minimal(base_size = 12)

save_plot(p_dec, glue("valid_quantiles_dr_{CFG$out_prefix}.png"))

# ----------------------------- Validation: Qini (DR gain curve) ----------------
qini_df <- qini_dr_curve(
  valid, tau_valid_cal,
  outcome = var_outcome, treat = var_treat,
  p_trt = mean(W_train == 1)
)
save_table(qini_df, glue("qini_dr_{CFG$out_prefix}.csv"))

p_qini <- ggplot(qini_df, aes(x = frac, y = cum_gain)) +
  geom_line() + geom_point(size = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(title = "Qini (Cumulative DR gain) — sorted by predicted CATE",
       x = "Targeting fraction (highest CATE first)",
       y = "Cumulative DR gain") +
  theme_minimal(base_size = 12)

save_plot(p_qini, glue("qini_dr_{CFG$out_prefix}.png"))

# ----------------------------- Save artifact (optional) ------------------------
artifact <- list(
  models   = models,
  x_cols   = colnames(X_train),
  scale_a  = scale_a,
  shift_b  = shift_b,
  cfg      = CFG
)
saveRDS(artifact, file.path(out_dir_models, glue("cf_ensemble_{CFG$out_prefix}.rds")))

# Save model column names for transparency
readr::write_csv(tibble(col = colnames(X_train)),
                 file.path(out_dir_tbl, glue("all_column_names_{CFG$out_prefix}.csv")))

message("Completed: ensemble causal forest, pooled ATE/BLP, VI (mapped & aggregated), calibration, validation DR quantiles, and Qini outputs.")
