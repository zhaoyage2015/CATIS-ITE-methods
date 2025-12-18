suppressPackageStartupMessages({
library(dplyr)
library(readr)
library(janitor)
library(ggplot2)
library(tibble)
})

------------------------------------------------------------------------------
04_validate_q5_dr_with_rates.R (Methods-only / Representative script)
Validation cohort:
- Stratify patients into quintiles based on predicted CATE (tau_hat)
- Within each quintile, estimate observed treatment effect as DR/AIPW risk difference
- Report event counts, event rates by randomized arm, and DR with 95% CI
- Add an Overall row for the entire validation cohort
Required inputs:
- data/clean/valid_clean.rds (must contain: treatment, primary_discharge, and an ID)
- outputs/cate_valid_primary_discharge.csv (first column = ID; includes tau_hat)
------------------------------------------------------------------------------
---------------- Parameters ----------------

nq <- 5
path_valid <- "data/clean/valid_clean.rds"
path_tau <- "outputs/cate_valid_primary_discharge.csv"
out_csv <- "outputs/validation_rd_by_cate_q5_with_rates.csv"
out_fig <- "figures/validation_rd_by_cate_q5.pdf"

dir.create(dirname(out_csv), showWarnings = FALSE, recursive = TRUE)
dir.create(dirname(out_fig), showWarnings = FALSE, recursive = TRUE)

---------------- Read data ----------------

valid <- readr::read_rds(path_valid) %>% clean_names()
tau <- readr::read_csv(path_tau, show_col_types = FALSE) %>% clean_names()

---------------- Identify / create ID column ----------------

id_col <- intersect(c("id","patient_id","subject_id","study_id","row_id"), names(valid))[1]
if (is.na(id_col)) {
id_col <- "row_id"
valid[[id_col]] <- seq_len(nrow(valid))
}

First column of tau is treated as ID

tau <- tau %>%
rename(!!id_col := names(.)[1]) %>%
mutate(!!id_col := as.character(.data[[id_col]]))

valid <- valid %>% mutate(!!id_col := as.character(.data[[id_col]]))

---------------- Merge tau_hat onto validation data ----------------

df <- valid %>% left_join(tau, by = id_col)

n_missing_tau <- sum(is.na(df$tau_hat))
if (n_missing_tau > 0) {
message("Missing tau_hat after join: ", n_missing_tau, " (these rows will be excluded).")
}
df <- df %>% filter(!is.na(tau_hat))

Basic checks

stopifnot(all(c("treatment","primary_discharge","tau_hat") %in% names(df)))

---------------- Core function: DR by predicted-CATE quantiles ----------------

compute_dr_by_quantile <- function(y, w, tau_hat, q = 5, p_t = NULL) {
stopifnot(length(y) == length(w), length(y) == length(tau_hat))
y <- as.numeric(y); w <- as.numeric(w); tau_hat <- as.numeric(tau_hat)

if (is.null(p_t)) p_t <- mean(w == 1, na.rm = TRUE)
pc <- 1 - p_t

brks <- quantile(tau_hat, probs = seq(0, 1, length.out = q + 1), na.rm = TRUE)
bins <- if (anyDuplicated(brks)) {
dplyr::ntile(tau_hat, q)
} else {
cut(tau_hat, breaks = brks, include.lowest = TRUE, labels = FALSE)
}

Z <- y * (w / p_t) - y * ((1 - w) / pc)

tibble(bin = as.integer(bins), y, w, tau = tau_hat, Z) %>%
group_by(bin) %>%
summarise(
n = n(),
n_trt = sum(w == 1, na.rm = TRUE),
n_ctl = sum(w == 0, na.rm = TRUE),

  events_trt = sum(y == 1 & w == 1, na.rm = TRUE),
  events_ctl = sum(y == 1 & w == 0, na.rm = TRUE),

  rate_trt = if_else(n_trt > 0, events_trt / n_trt, NA_real_),
  rate_ctl = if_else(n_ctl > 0, events_ctl / n_ctl, NA_real_),
  rate_overall = mean(y, na.rm = TRUE),

  tau_mean = mean(tau, na.rm = TRUE),

  rd_dr = mean(Z, na.rm = TRUE),
  n_eff = sum(!is.na(Z)),
  se    = sd(Z, na.rm = TRUE) / sqrt(n_eff),
  ci_l  = rd_dr - 1.96 * se,
  ci_u  = rd_dr + 1.96 * se,
  .groups = "drop"
) %>%
arrange(bin) %>%
mutate(
  q_f = factor(paste0("Q", bin), levels = paste0("Q", 1:q)),
  rate_trt_pct = sprintf("%.1f%%", 100 * rate_trt),
  rate_ctl_pct = sprintf("%.1f%%", 100 * rate_ctl)
)


}

---------------- Quintile results ----------------

qres <- compute_dr_by_quantile(
y = df$primary_discharge,
w = df$treatment,
tau_hat = df$tau_hat,
q = nq,
p_t = mean(df$treatment == 1, na.rm = TRUE)
)

---------------- Overall row ----------------

p_t <- mean(df$treatment == 1, na.rm = TRUE)
pc <- 1 - p_t
Z_all <- with(df,
primary_discharge * (treatment / p_t) -
primary_discharge * ((1 - treatment) / pc))

overall_counts <- df %>%
summarise(
n = n(),
n_trt = sum(treatment == 1, na.rm = TRUE),
n_ctl = sum(treatment == 0, na.rm = TRUE),
events_trt = sum(primary_discharge == 1 & treatment == 1, na.rm = TRUE),
events_ctl = sum(primary_discharge == 1 & treatment == 0, na.rm = TRUE),
rate_trt = if_else(n_trt > 0, events_trt / n_trt, NA_real_),
rate_ctl = if_else(n_ctl > 0, events_ctl / n_ctl, NA_real_),
rate_overall = mean(primary_discharge, na.rm = TRUE),
.groups = "drop"
)

overall_row <- overall_counts %>%
mutate(
tau_mean = mean(df$tau_hat, na.rm = TRUE),
rd_dr = mean(Z_all, na.rm = TRUE),
n_eff = sum(!is.na(Z_all)),
se = sd(Z_all, na.rm = TRUE) / sqrt(n_eff),
ci_l = rd_dr - 1.96 * se,
ci_u = rd_dr + 1.96 * se,
bin = 0L,
q_f = factor("Overall", levels = c("Overall", paste0("Q", 1:nq))),
rate_trt_pct = sprintf("%.1f%%", 100 * rate_trt),
rate_ctl_pct = sprintf("%.1f%%", 100 * rate_ctl)
) %>%
select(
bin, n, n_trt, n_ctl, events_trt, events_ctl,
rate_trt, rate_ctl, rate_overall,
tau_mean, rd_dr, n_eff, se, ci_l, ci_u, q_f, rate_trt_pct, rate_ctl_pct
)

qres_sel <- qres %>%
select(
bin, n, n_trt, n_ctl, events_trt, events_ctl,
rate_trt, rate_ctl, rate_overall,
tau_mean, rd_dr, n_eff, se, ci_l, ci_u, q_f, rate_trt_pct, rate_ctl_pct
)

final_tbl <- bind_rows(overall_row, qres_sel) %>%
mutate(q_f = factor(as.character(q_f), levels = c("Overall", paste0("Q", 1:nq))))

write_csv(final_tbl, out_csv)
message("Saved table: ", out_csv)

---------------- Plot (vector PDF) ----------------

plot_dat <- final_tbl %>%
select(q_f, bin, rd_dr, ci_l, ci_u) %>%
arrange(bin)

g <- ggplot(plot_dat, aes(x = q_f, y = rd_dr)) +
geom_point(size = 2.6) +
geom_errorbar(aes(ymin = ci_l, ymax = ci_u), width = 0.15, linewidth = 0.8) +
geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.8) +
labs(
x = "Quintile of predicted treatment effect (τ̂)",
y = "Observed treatment effect: DR (treated − control)"
) +
theme_classic(base_size = 12) +
theme(
axis.title.x = element_text(face = "bold", size = 12, margin = margin(t = 10)),
axis.title.y = element_text(face = "bold", size = 12, margin = margin(r = 10)),
axis.text = element_text(size = 11),
axis.ticks = element_line(linewidth = 0.8),
axis.line = element_line(linewidth = 0.8)
)

pdf_save <- function(path, plot, w = 7.2, h = 4.8) {
ok <- FALSE
try({
ggsave(path, plot, width = w, height = h, dpi = 600, device = cairo_pdf)
ok <- TRUE
}, silent = TRUE)

if (!ok) {
ggsave(path, plot, width = w, height = h, dpi = 600, device = "pdf")
}
}

pdf_save(out_fig, g)
message("Saved figure: ", out_fig)_
