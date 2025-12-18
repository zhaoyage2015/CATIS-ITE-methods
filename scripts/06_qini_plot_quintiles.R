suppressPackageStartupMessages({
library(dplyr)
library(readr)
library(janitor)
library(ggplot2)
library(grid)
library(scales)
library(tidyr)
})

------------------------------------------------------------------------------
06_qini_plot_quintiles.R (Methods-only / Representative script)
Purpose:
Create a publication-style Qini/uplift curve in the validation cohort using
DR/AIPW pseudo-outcomes, with curve values reported at quintile nodes
(0%, 20%, 40%, 60%, 80%, 100% targeted).
Required inputs:
- data/clean/valid_clean.rds (must contain: treatment, primary_discharge, and an ID)
- outputs/cate_valid_primary_discharge.csv (first column = ID; includes tau_hat)
Outputs:
- outputs/qini_curve_quintiles.csv
- figures/qini_plot_quintiles.{pdf,tiff,png}
Notes:
- Outcome is coded as 1 = poor outcome (unfavorable).
- DR/AIPW pseudo-outcome Z yields a risk-difference scale (treated − control).
- Sorting direction is auto-checked so that larger scores correspond to
"more favorable for treatment" in the plotted targeting order.
------------------------------------------------------------------------------
---------------- Parameters ----------------

path_valid <- "data/clean/valid_clean.rds"
path_tau <- "outputs/cate_valid_primary_discharge.csv"

out_curve <- "outputs/qini_curve_quintiles.csv"
out_pdf <- "figures/qini_plot_quintiles.pdf"
out_tiff <- "figures/qini_plot_quintiles.tiff"
out_png <- "figures/qini_plot_quintiles.png"

dir.create(dirname(out_curve), showWarnings = FALSE, recursive = TRUE)
dir.create(dirname(out_pdf), showWarnings = FALSE, recursive = TRUE)
dir.create(dirname(out_tiff), showWarnings = FALSE, recursive = TRUE)
dir.create(dirname(out_png), showWarnings = FALSE, recursive = TRUE)

---------------- Read and merge ----------------

valid <- readr::read_rds(path_valid) %>% clean_names()
tau <- readr::read_csv(path_tau, show_col_types = FALSE) %>% clean_names()

id_col <- intersect(c("id","patient_id","subject_id","study_id","row_id"), names(valid))[1]
if (is.na(id_col)) {
id_col <- "row_id"
valid[[id_col]] <- seq_len(nrow(valid))
}

tau <- tau %>%
rename(!!id_col := names(.)[1]) %>%
mutate(!!id_col := as.character(.data[[id_col]]))

df <- valid %>%
mutate(!!id_col := as.character(.data[[id_col]])) %>%
left_join(tau, by = id_col)

stopifnot(all(c("treatment", "primary_discharge", "tau_hat") %in% names(df)))

Keep rows with tau_hat

df <- df %>% filter(!is.na(tau_hat))

---------------- DR/AIPW pseudo-outcome Z (risk-difference scale) ----------------

p_t <- mean(df$treatment == 1, na.rm = TRUE)
pc <- 1 - p_t

df <- df %>%
mutate(
Z = primary_discharge * (treatment / p_t) -
primary_discharge * ((1 - treatment) / pc)
)

---------------- Direction check: ensure "higher score = prioritize treatment" ----------------

correlation <- suppressWarnings(cor(df$tau_hat, df$Z, use = "complete.obs"))
dir_sign <- sign(correlation)

decreasing_sort <- ifelse(is.na(dir_sign) || dir_sign >= 0, TRUE, FALSE)
df_ord <- df[order(df$tau_hat, decreasing = decreasing_sort), , drop = FALSE]
n <- nrow(df_ord)

---------------- Quintile nodes (0,20,40,60,80,100%) ----------------

K <- 5
props <- c(0, (1:K) / K)

idx_m <- floor(props * n)
idx_m[1] <- 0

uplift_model <- sapply(idx_m, function(m) {
if (m <= 0) return(0)
mean(df_ord$Z[1:m], na.rm = TRUE) * 100
})
uplift_model[1] <- 0

overall_rd_pp <- mean(df$Z, na.rm = TRUE) * 100
uplift_random <- props * overall_rd_pp

"Oracle" upper bound based on observed Z (not plotted; used for adjusted Qini)

ord_oracle <- order(df$Z, decreasing = TRUE)
uplift_oracle <- sapply(idx_m, function(m) {
if (m <= 0) return(0)
mean(df$Z[ord_oracle][1:m], na.rm = TRUE) * 100
})
uplift_oracle[1] <- 0

---------------- Area-based Qini (trapezoid over the 6 nodes) ----------------

trapz <- function(x, y) sum((head(y, -1) + tail(y, -1)) / 2 * diff(x))

A_model <- trapz(props, uplift_model)
A_random <- trapz(props, uplift_random)
A_oracle <- trapz(props, uplift_oracle)

qini_raw <- A_model - A_random
qini_adj <- if (A_oracle > A_random) (A_model - A_random) / (A_oracle - A_random) else NA_real_

message(sprintf("Qini (area difference, pp): %.3f", qini_raw))
message(sprintf("Adjusted Qini (relative to oracle): %s",
ifelse(is.na(qini_adj), "NA", sprintf("%.3f", qini_adj))))

---------------- Export curve data ----------------

curve_dat <- tibble(
prop = props,
uplift_model = uplift_model,
uplift_random = uplift_random,
uplift_oracle = uplift_oracle
)
readr::write_csv(curve_dat, out_curve)

---------------- Plot settings (fixed axis range as in manuscript) ----------------

y_lim <- c(-3, 5)
y_breaks <- seq(-2.5, 5, 2.5)
y_top <- y_lim[2]
y_bot <- y_lim[1]

Long format for legend

curve_long <- curve_dat %>%
select(
prop,
Model-based targeting = uplift_model,
Random targeting = uplift_random
) %>%
pivot_longer(-prop, names_to = "curve", values_to = "uplift")

g <- ggplot(curve_long, aes(x = prop, y = uplift, linetype = curve, shape = curve)) +
geom_hline(yintercept = 0, linetype = "dotted", linewidth = 0.7, inherit.aes = FALSE) +
geom_line(linewidth = 1.0) +
geom_point(data = subset(curve_long, curve == "Model-based targeting"), size = 2.2) +
geom_point(data = subset(curve_long, curve == "Random targeting"), size = 2.2, alpha = 0) +
scale_x_continuous(
labels = scales::number_format(accuracy = 1, scale = 100),
breaks = seq(0, 1, 0.25)
) +
scale_y_continuous(limits = y_lim, breaks = y_breaks) +
scale_linetype_manual(values = c("Model-based targeting" = "solid",
"Random targeting" = "dashed")) +
scale_shape_manual(values = c("Model-based targeting" = 16,
"Random targeting" = NA)) +
labs(
x = "Proportion of population targeted (%)",
y = "Incremental uplift (risk difference, pp)"
) +
annotate(
"text", x = 0.82, y = y_top,
label = sprintf("Qini: %.2f", qini_raw),
hjust = 0, vjust = -0.4, size = 4.2
) +
theme_classic(base_size = 12) +
theme(
axis.title.x = element_text(face = "bold", margin = margin(t = 8)),
axis.title.y = element_text(face = "bold", margin = margin(r = 8)),
plot.margin = unit(c(5.5, 30, 5.5, 40), "pt"),
legend.position = c(0.30, 0.14),
legend.background = element_rect(fill = "white", color = "grey80"),
legend.title = element_blank()
) +
coord_cartesian(clip = "off") +
annotate(
"segment", x = -0.02, xend = -0.02, y = y_top, yend = 0,
arrow = arrow(type = "closed", length = unit(5, "pt"))
) +
annotate(
"segment", x = -0.02, xend = -0.02, y = y_bot, yend = 0,
arrow = arrow(type = "closed", length = unit(5, "pt"))
) +
annotate("text", x = -0.06, y = 0.6 * y_top, label = "Favors Treatment",
angle = 90, size = 3.3) +
annotate("text", x = -0.06, y = 0.6 * y_bot, label = "Favors Control",
angle = 90, size = 3.3)

---------------- Save figures ----------------

save_pdf <- function(path, plot, w = 6.8, h = 4.6) {
ok <- FALSE
try({
ggsave(path, plot, width = w, height = h, device = cairo_pdf, dpi = 600)
ok <- TRUE
}, silent = TRUE)
if (!ok) ggsave(path, plot, width = w, height = h, device = "pdf", dpi = 600)
}

save_pdf(out_pdf, g)
ggsave(out_tiff, g, width = 6.8, height = 4.6, dpi = 600, compression = "lzw")
ggsave(out_png, g, width = 6.8, height = 4.6, dpi = 600)

message("Saved: ", out_curve)
message("Saved: ", out_pdf)
message("Saved: ", out_tiff)
message("Saved: ", out_png)
