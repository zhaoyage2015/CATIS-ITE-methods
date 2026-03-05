#!/usr/bin/env Rscript
# ============================================================
# 06_pdp_top6features_plot.R
# Strict Partial Dependence Plots (PDP) for top-6 features (ITE = treated − control)
#
# Key behaviors:
# - Mirror curve around X-axis only: y -> -y (swap lo/hi safely)
# - Fixed panel order: Creatinine, BMI, hs-CRP, NT-proBNP, HCT, Hcy
# - Asymmetric Y-axis: upper = 0.005; lower at least -0.02
# - Output to a timestamped folder under figures/pdp_strict (no overwrite)
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(purrr)
  library(tibble)
  library(grf)
})

# --------------------------- Config ---------------------------
# Project root:
# - Option A (recommended): set env var CATIS_ITE_DIR
# - Option B: auto-detect by walking up from this script location / working dir
get_proj_dir <- function() {
  env <- Sys.getenv("CATIS_ITE_DIR", unset = NA_character_)
  if (!is.na(env) && nzchar(env) && dir.exists(env)) return(normalizePath(env, winslash = "/", mustWork = TRUE))

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

# Inputs (relative to project root)
path_valid <- file.path(proj_dir, "data/clean/valid_clean.rds")
path_train <- file.path(proj_dir, "data/clean/train_clean.rds")
path_vi1   <- file.path(proj_dir, "outputs/vi_full_primary_discharge.csv")
path_vi2   <- file.path(proj_dir, "outputs/vi_primary_discharge.csv")
model_path <- file.path(proj_dir, "models/cf_ensemble_primary_discharge.rds")

# Parameters
K_top  <- 6
B_boot <- 400
N_grid <- 40

# Fixed top-6 raw feature order (must match column names in data)
desired_raw <- c("creatinine", "bmi", "hs_crp", "pro_bnp", "hct", "hcy")

# Pretty labels (raw -> display)
label_map <- c(
  creatinine = "Creatinine",
  bmi        = "BMI",
  hs_crp     = "hs-CRP",
  pro_bnp    = "NT-proBNP",
  hct        = "HCT",
  hcy        = "Hcy"
)
pretty_var <- function(v) {
  v <- as.character(v)
  hit <- v %in% names(label_map)
  v[hit] <- unname(label_map[v[hit]])
  v
}

# Y-axis settings
POS_Y_MAX <- 0.005
NEG_MIN_AT_LEAST <- 0.02
Y_BREAKS <- c(-0.02, -0.01, 0, 0.005)

# Y-axis title
y_title_main <- "Mean predicted ITE (DR = treated − control)"
y_title_sub  <- "Favors Treatment        Favors Control   "
y_title_full <- paste0(y_title_main, "\n", y_title_sub)

# Output paths (timestamped, no overwrite)
fig_dir_base <- file.path(proj_dir, "figures", "pdp_strict")
run_tag      <- format(Sys.time(), "%Y%m%d_%H%M%S")
subdir_name  <- paste0("pdp_strict_top6_reflectY_", run_tag)
fig_dir_new  <- file.path(fig_dir_base, subdir_name)
dir.create(fig_dir_new, showWarnings = FALSE, recursive = TRUE)

out_multi_pdf <- file.path(fig_dir_new, sprintf("pdp_strict_top6_multi_%s.pdf", run_tag))
out_multi_png <- file.path(fig_dir_new, sprintf("pdp_strict_top6_multi_%s.png", run_tag))

out_single_dir <- file.path(fig_dir_new, "single")
dir.create(out_single_dir, showWarnings = FALSE, recursive = TRUE)

# ------------------------ Helper utils ------------------------
assert_exists <- function(p, what = "file") {
  if (!file.exists(p)) stop("Missing ", what, ": ", p, call. = FALSE)
}

to_numeric_matrix <- function(df) {
  df2 <- df
  bad <- names(df2)[!sapply(df2, is.numeric)]
  if (length(bad) > 0) {
    for (nm in bad) {
      v <- df2[[nm]]
      if (is.logical(v)) {
        df2[[nm]] <- as.numeric(v)
      } else if (is.factor(v)) {
        if (all(levels(v) %in% c("0", "1"))) {
          df2[[nm]] <- as.numeric(as.character(v))
        } else {
          suppressWarnings(num <- as.numeric(as.character(v)))
          if (any(is.na(num) & !is.na(v))) {
            stop("Non-numeric factor column '", nm, "'. Please one-hot encode in preprocessing.", call. = FALSE)
          }
          df2[[nm]] <- num
        }
      } else if (is.character(v)) {
        suppressWarnings(num <- as.numeric(v))
        if (any(is.na(num) & !is.na(v))) stop("Non-numeric characters in '", nm, "'.", call. = FALSE)
        df2[[nm]] <- num
      } else {
        suppressWarnings(df2[[nm]] <- as.numeric(df2[[nm]]))
      }
    }
  }
  data.matrix(df2)
}

align_features <- function(df, feats) {
  out <- df
  miss <- setdiff(feats, names(out))
  if (length(miss) > 0) for (nm in miss) out[[nm]] <- 0
  out <- out[, feats, drop = FALSE]
  out
}

collect_grf <- function(obj, depth = 0, max_depth = 10) {
  if (depth > max_depth) return(list())
  if (inherits(obj, "grf")) return(list(obj))
  if (!is.list(obj)) return(list())
  out <- list()
  slots <- c("members", "models", "forests", "estimators", "ensemble", "cf_list", "objects")
  for (nm in slots) if (!is.null(obj[[nm]])) out <- c(out, collect_grf(obj[[nm]], depth + 1, max_depth))
  for (el in obj) out <- c(out, collect_grf(el, depth + 1, max_depth))
  out
}

save_pdf_safe <- function(path, plot, width, height) {
  # Try cairo_pdf; fallback to default pdf device
  ok <- FALSE
  try({
    ggsave(path, plot, width = width, height = height, device = cairo_pdf)
    ok <- TRUE
  }, silent = TRUE)
  if (!ok) ggsave(path, plot, width = width, height = height)
}

# ------------------------- Load assets -------------------------
assert_exists(path_valid, "valid_clean.rds")
assert_exists(path_train, "train_clean.rds")
assert_exists(model_path, "model rds")

valid <- readRDS(path_valid)
train <- readRDS(path_train)
model <- readRDS(model_path)

# Feature names
drop_cols <- c("treatment", "primary_discharge", "primary_m3", "id", "patient_id")
feature_names <- unique(setdiff(names(train), drop_cols))

# Align Xvalid
miss_valid <- setdiff(feature_names, names(valid))
if (length(miss_valid) > 0) message("⚠️ Valid missing columns filled with 0: ", paste(miss_valid, collapse = ", "))
Xvalid <- align_features(valid, feature_names)

# Optional read VI (kept for reference; not used for selecting top6)
vi_path <- if (file.exists(path_vi1)) path_vi1 else path_vi2
if (file.exists(vi_path)) {
  vi <- read_csv(vi_path, show_col_types = FALSE)
  invisible(vi)
}

# Fixed top6 selection (ensure numeric and exists)
top_vars <- desired_raw[desired_raw %in% names(Xvalid)]
top_vars <- top_vars[sapply(Xvalid[top_vars], is.numeric)]
if (length(top_vars) == 0) stop("No desired variables found in Xvalid. Check column names in train/valid.", call. = FALSE)
if (length(top_vars) > K_top) top_vars <- top_vars[seq_len(K_top)]
stopifnot(is.character(top_vars))

message("Top variables (fixed order): ", paste(top_vars, collapse = ", "))
message("Display labels: ", paste(pretty_var(top_vars), collapse = ", "))

# --------------------- Robust predict wrapper ------------------
.predict_cf_safe <- function(cf, newX_mat) {
  out <- tryCatch(predict(cf, newdata = newX_mat), error = function(e) e)
  if (inherits(out, "error")) out <- getS3method("predict", "causal_forest")(cf, newdata = newX_mat)
  as.numeric(out$predictions)
}

predict_fun <- function(newdata_df) {
  newX_df  <- align_features(as.data.frame(newdata_df), feature_names)
  newX_mat <- to_numeric_matrix(newX_df)

  if (inherits(model, "grf")) return(.predict_cf_safe(model, newX_mat))

  members <- Filter(function(x) inherits(x, "grf"), collect_grf(model))
  if (length(members) > 0) {
    preds_list <- lapply(members, function(m) .predict_cf_safe(m, newX_mat))
    preds_mat  <- do.call(cbind, preds_list)
    return(rowMeans(preds_mat, na.rm = TRUE))
  }

  as.numeric(stats::predict(model, newdata = newX_mat))
}

# Smoke test
invisible(predict_fun(head(Xvalid, 5)))

# ---------------------- Strict PDP compute ---------------------
strict_pdp_one <- function(var, dfX, n_grid = 40, B = 400) {
  xv <- dfX[[var]]
  grid <- quantile(xv, probs = seq(0.02, 0.98, length.out = n_grid), na.rm = TRUE) |> as.numeric()

  preds_mat <- sapply(grid, function(g) {
    newdata <- dfX
    newdata[[var]] <- g
    predict_fun(newdata)
  })

  n <- nrow(dfX)
  means_boot <- replicate(B, {
    idx <- sample.int(n, n, replace = TRUE)
    colMeans(preds_mat[idx, , drop = FALSE])
  })

  tibble(
    var = var,
    x   = grid,
    y   = rowMeans(means_boot, na.rm = TRUE),
    lo  = apply(means_boot, 1, function(z) quantile(z, 0.025, na.rm = TRUE)),
    hi  = apply(means_boot, 1, function(z) quantile(z, 0.975, na.rm = TRUE))
  )
}

pdp_df <- map_dfr(top_vars, ~ strict_pdp_one(.x, dfX = Xvalid, n_grid = N_grid, B = B_boot))

# Ensure facet order
pdp_df <- pdp_df %>%
  mutate(
    var     = factor(var, levels = top_vars),
    var_lab = factor(pretty_var(as.character(var)), levels = pretty_var(top_vars))
  )

# Mirror around X-axis only: y -> -y, swap CI bounds safely
pdp_df <- pdp_df %>%
  mutate(
    y   = -y,
    tmp = -lo,
    lo  = -hi,
    hi  = tmp
  ) %>% select(-tmp)

# Asymmetric y-limits
neg_span <- max(
  abs(pdp_df$lo[pdp_df$lo < 0]),
  abs(pdp_df$hi[pdp_df$hi < 0]),
  abs(pdp_df$y [pdp_df$y  < 0]),
  na.rm = TRUE
)
if (!is.finite(neg_span)) neg_span <- NEG_MIN_AT_LEAST
NEG_Y_MIN <- -max(NEG_MIN_AT_LEAST, neg_span)
ylim_use  <- c(NEG_Y_MIN, POS_Y_MAX)

# ------------------------ Plot: multi --------------------------
p_multi <- ggplot(pdp_df, aes(x = x, y = y)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.35, fill = "grey40", na.rm = TRUE) +
  geom_line(linewidth = 0.9, color = "black", na.rm = TRUE) +
  geom_hline(yintercept = 0, linetype = 2, color = "black") +
  facet_wrap(~ var_lab, scales = "free_x", ncol = 2) +
  labs(x = "Predictor Variable Setting", y = y_title_full) +
  scale_y_continuous(limits = ylim_use, breaks = Y_BREAKS) +
  theme_bw(base_size = 12) +
  theme(
    text        = element_text(color = "black", face = "bold"),
    axis.text   = element_text(color = "black", face = "bold", size = 11),
    axis.title  = element_text(color = "black", face = "bold", size = 12),
    strip.text  = element_text(face = "bold", size = 13),
    plot.title  = element_text(face = "bold", hjust = 0.5, size = 14),
    axis.ticks  = element_line(linewidth = 0.5, color = "black")
  )

save_pdf_safe(out_multi_pdf, p_multi, width = 8.5, height = 9.0)
ggsave(out_multi_png, p_multi, width = 8.5, height = 9.0, dpi = 320)

# ------------------------ Plot: single -------------------------
plot_one_var <- function(df, var_raw) {
  d <- df %>% filter(.data$var == !!var_raw)
  ggplot(d, aes(x = x, y = y)) +
    geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.35, fill = "grey40", na.rm = TRUE) +
    geom_line(linewidth = 0.9, color = "black", na.rm = TRUE) +
    geom_hline(yintercept = 0, linetype = 2, color = "black") +
    scale_y_continuous(limits = ylim_use, breaks = Y_BREAKS) +
    labs(
      title = pretty_var(var_raw),
      x = "Predictor Variable Setting",
      y = y_title_full
    ) +
    theme_bw(base_size = 12) +
    theme(
      text        = element_text(color = "black", face = "bold"),
      axis.text   = element_text(color = "black", face = "bold", size = 11),
      axis.title  = element_text(color = "black", face = "bold", size = 12),
      plot.title  = element_text(face = "bold", hjust = 0.5, size = 14),
      axis.ticks  = element_line(linewidth = 0.5, color = "black")
    )
}

for (v in top_vars) {
  p_one <- plot_one_var(pdp_df, v)
  pdf_path <- file.path(out_single_dir, sprintf("pdp_strict_%s_%s.pdf", v, run_tag))
  png_path <- file.path(out_single_dir, sprintf("pdp_strict_%s_%s.png", v, run_tag))
  save_pdf_safe(pdf_path, p_one, width = 5.2, height = 4.2)
  ggsave(png_path, p_one, width = 5.2, height = 4.2, dpi = 320)
}

message("✅ Multi saved: ", out_multi_png, " / ", out_multi_pdf)
message("✅ Singles saved to: ", out_single_dir)
message("📁 Output folder: ", fig_dir_new)
message("📌 Project root used: ", proj_dir)
