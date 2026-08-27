#!/usr/bin/env Rscript
# 10: Calibrate the expression- and splicing-response estimators against known ground truth.
#
# Run: Rscript scripts/10_layer_estimator_calibration.R
# Simulation only -- reads nothing except config. Deterministic (seed fixed in config).
#
# WHY THIS EXISTS
#
# The load-bearing claim in the UBL5 analysis is a comparison BETWEEN two layers: "the
# knockout retains X% of the expression response but only Y% of the splicing response". That
# comparison is only meaningful if both estimators respond to a change in true signal the
# same way. They do not.
#
# Expression response is estimated as |log2FC| between condition means of three replicates.
# Averaging three replicates before taking the absolute value gives a high signal-to-noise
# ratio, so the estimate tracks the truth almost linearly.
#
# Splicing response is a total variation distance -- a sum of absolute values over isoforms.
# Absolute values do not average out noise, so for signal small relative to noise,
# E|s + n| ~ E|n| + O(s^2). Any noise-subtracted TVD therefore responds roughly QUADRATICALLY
# to signal near zero, and compresses ratios: a genotype with 60% of the true signal reads
# as ~44%.
#
# Consequence: comparing an almost-unbiased expression ratio against a compressed splicing
# ratio manufactures an apparent asymmetry even when both layers are reduced by exactly the
# same factor. This script measures the compression for each estimator and inverts it, so an
# observed ratio can be converted to a true-signal ratio before the layers are compared.
#
# Two splicing estimators are calibrated, because both appear in this project's history:
#   A  mean pairwise across-condition TVD minus mean within-condition TVD (external review)
#   B  TVD of the mean per-pair dIF vector, minus a paired permutation null (script 08)

bt <- c("utils/bootstrap.R", file.path("..", "utils", "bootstrap.R"))
b_file <- if (any(f <- vapply(bt, file.exists, NA))) { bt[which(f)[1L]] } else { NA_character_ }
if (is.na(b_file)) {
  stop("Could not find utils/bootstrap.R. Run from project root or set ISOFORM_PROJECT_ROOT.")
}
source(b_file, local = FALSE, chdir = FALSE)

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(ggplot2)
})

root <- if (exists("PROJECT_ROOT", inherits = FALSE)) {
  get("PROJECT_ROOT", inherits = FALSE)
} else {
  find_project_root()
}
args_r <- commandArgs(trailingOnly = TRUE)
config_rel <- if (length(args_r) > 0L) args_r[[1L]] else "config/config.yml"
cfg <- load_yaml_config(config_rel, root = root)
paths <- cfg$paths %||% list()

results_tables <- ensure_dir(resolve_path(paths$results_tables %||% "results/tables", root = root))
results_figures <- ensure_dir(resolve_path(paths$results_figures %||% "results/figures", root = root))
fig_dir <- ensure_dir(file.path(results_figures, "layers"))

n_genes <- as.integer(analysis_param(cfg, "calibration_n_genes", 4000L))
n_iso <- as.integer(analysis_param(cfg, "calibration_n_isoforms", 6L))
n_rep <- as.integer(analysis_param(cfg, "calibration_n_replicates", 3L))
if_sigma <- as.numeric(analysis_param(cfg, "calibration_if_noise_sd", 0.05))
expr_sigma <- as.numeric(analysis_param(cfg, "calibration_expr_noise_sd", 0.35))
seed <- as.integer(analysis_param(cfg, "calibration_seed", 42L))
grid <- as.numeric(analysis_param(cfg, "calibration_ratio_grid",
                                  c(1, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2)))
theme_set(theme_bw(base_size = 11))
set.seed(seed)

# True per-gene composition shift: two isoforms move in opposite directions, magnitude drawn
# from a distribution matched to the observed |dIF| scale.
true_shift <- matrix(0, n_genes, n_iso)
for (g in seq_len(n_genes)) {
  a <- sample.int(n_iso, 2)
  s <- abs(stats::rnorm(1, 0, 0.08))
  true_shift[g, a[1]] <- s
  true_shift[g, a[2]] <- -s
}
true_lfc <- stats::rnorm(n_genes, 0, 0.6)

simulate <- function(scale) {
  base <- matrix(stats::runif(n_genes * n_iso, 0, 1), n_genes, n_iso)
  base <- base / rowSums(base)
  minus <- plus <- array(NA_real_, c(n_genes, n_iso, n_rep))
  e_minus <- e_plus <- matrix(NA_real_, n_genes, n_rep)
  for (i in seq_len(n_rep)) {
    minus[, , i] <- base + matrix(stats::rnorm(n_genes * n_iso, 0, if_sigma), n_genes, n_iso)
    plus[, , i] <- base + scale * true_shift +
      matrix(stats::rnorm(n_genes * n_iso, 0, if_sigma), n_genes, n_iso)
    e_minus[, i] <- stats::rnorm(n_genes, 0, expr_sigma)
    e_plus[, i] <- scale * true_lfc + stats::rnorm(n_genes, 0, expr_sigma)
  }
  list(minus = minus, plus = plus, e_minus = e_minus, e_plus = e_plus)
}

tvd <- function(a, b) 0.5 * rowSums(abs(a - b))

# Expression: |log2FC| between condition means -- the estimator actually used in the project.
est_expression <- function(d) abs(rowMeans(d$e_plus) - rowMeans(d$e_minus))

# A: mean pairwise across-condition TVD minus mean within-condition TVD.
est_splicing_A <- function(d) {
  acr <- rowMeans(sapply(seq_len(n_rep), function(i) {
    rowMeans(sapply(seq_len(n_rep), function(j) tvd(d$minus[, , i], d$plus[, , j])))
  }))
  prs <- utils::combn(n_rep, 2, simplify = FALSE)
  wit <- rowMeans(cbind(
    sapply(prs, function(p) tvd(d$minus[, , p[1]], d$minus[, , p[2]])),
    sapply(prs, function(p) tvd(d$plus[, , p[1]], d$plus[, , p[2]]))
  ))
  acr - wit
}

# B: TVD of the mean per-pair dIF vector, minus the paired permutation null mean.
est_splicing_B <- function(d) {
  dif <- lapply(seq_len(n_rep), function(i) d$plus[, , i] - d$minus[, , i])
  obs <- 0.5 * rowSums(abs(Reduce(`+`, dif) / n_rep))
  signs <- cbind(1, as.matrix(expand.grid(rep(list(c(1, -1)), n_rep - 1L))))
  nulls <- sapply(2:nrow(signs), function(k) {
    v <- Reduce(`+`, lapply(seq_len(n_rep), function(i) signs[k, i] * dif[[i]])) / n_rep
    0.5 * rowSums(abs(v))
  })
  obs - rowMeans(nulls)
}

message("Calibrating estimators over ", length(grid), " true-ratio values ...")
ref <- simulate(1)
rows <- lapply(grid, function(r) {
  alt <- simulate(r)
  tibble(
    true_ratio = r,
    expression = mean(est_expression(alt)) / mean(est_expression(ref)),
    splicing_A = mean(est_splicing_A(alt)) / mean(est_splicing_A(ref)),
    splicing_B = mean(est_splicing_B(alt)) / mean(est_splicing_B(ref))
  )
})
cal <- bind_rows(rows)
write_table_pair(cal, results_tables, "layer_estimator_calibration", cfg = cfg)
print(as.data.frame(cal))

# Invert: given an observed ratio, what true ratio does it imply?
invert <- function(observed, column) {
  x <- cal[[column]]; y <- cal$true_ratio
  o <- order(x)
  stats::approx(x[o], y[o], xout = observed, rule = 2)$y
}

obs_path <- file.path(results_tables, "layer_retention.csv")
ext_path <- resolve_path("docs/external_review/tables/ubl5_retention_by_layer.csv", root = root)
obs_rows <- list()
if (file.exists(obs_path)) {
  lr <- utils::read.csv(obs_path, stringsAsFactors = FALSE)
  e <- lr$retained_in_KO[lr$layer == "Expression"][1]
  s <- lr$retained_in_KO[lr$layer == "Splicing"][1]
  obs_rows[[1]] <- tibble(
    source = "this pipeline (script 08, estimator B)",
    observed_expression = e, observed_splicing = s,
    true_expression = invert(e, "expression"), true_splicing = invert(s, "splicing_B")
  )
}
if (file.exists(ext_path)) {
  er <- utils::read.csv(ext_path, stringsAsFactors = FALSE)
  er <- er[er$transcript_set == "all transcripts", ]
  e <- er$ratio[grepl("^Expression", er$layer)][1]
  s <- er$ratio[grepl("^Splicing", er$layer)][1]
  obs_rows[[length(obs_rows) + 1L]] <- tibble(
    source = "external review (estimator A)",
    observed_expression = e, observed_splicing = s,
    true_expression = invert(e, "expression"), true_splicing = invert(s, "splicing_A")
  )
}
if (length(obs_rows)) {
  corrected <- bind_rows(obs_rows) |>
    mutate(
      observed_gap = .data$observed_expression - .data$observed_splicing,
      corrected_gap = .data$true_expression - .data$true_splicing,
      note = paste(
        "Observed ratios are compressed by the estimator; the corrected columns invert the",
        "calibration. Compare layers only after correction."
      )
    )
  write_table_pair(corrected, results_tables, "layer_retention_corrected", cfg = cfg)
  print(as.data.frame(corrected |> select(-.data$note)))
}

long <- cal |>
  tidyr::pivot_longer(c("expression", "splicing_A", "splicing_B"),
                      names_to = "estimator", values_to = "observed") |>
  mutate(estimator = dplyr::recode(
    .data$estimator,
    expression = "Expression |log2FC| (3-replicate means)",
    splicing_A = "Splicing A: pairwise TVD, across minus within",
    splicing_B = "Splicing B: TVD of mean dIF, minus paired null"
  ))
p <- ggplot(long, aes(x = .data$true_ratio, y = .data$observed, colour = .data$estimator)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey40") +
  geom_line(linewidth = 0.9) +
  geom_point(size = 1.8) +
  scale_colour_manual(values = c(
    "Expression |log2FC| (3-replicate means)" = "#3182bd",
    "Splicing A: pairwise TVD, across minus within" = "#e6550d",
    "Splicing B: TVD of mean dIF, minus paired null" = "#756bb1"
  ), name = NULL) +
  coord_equal(xlim = c(0, 1.05), ylim = c(0, 1.05)) +
  labs(
    title = "Both splicing estimators compress the ratio; the expression estimator does not",
    subtitle = "Dashed line = unbiased. A layer comparison is invalid unless both are on this line.",
    x = "True ratio of signal (KO / WT)", y = "Ratio reported by the estimator"
  ) +
  theme(legend.position = "bottom", legend.direction = "vertical")
ggsave(file.path(fig_dir, "fig_estimator_calibration.png"), p, width = 8, height = 7, dpi = 200)

write_run_manifest("10_layer_estimator_calibration.R", cfg, root)
message("10_layer_estimator_calibration.R: done")
