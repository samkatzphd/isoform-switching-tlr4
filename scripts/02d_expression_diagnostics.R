#!/usr/bin/env Rscript
# 02d: Expression diagnostics, and a data-driven abundance floor for switching calls.
#
# Run (after 01): Rscript scripts/02d_expression_diagnostics.R
#
# WHY THIS EXISTS
#
# dIF is a ratio. When a gene is barely expressed the denominator is tiny, so isoform
# fraction swings between replicates for reasons that have nothing to do with the condition.
# The symptom is visible without any modelling: switching calls are ~4x MORE common in the
# lowest expression quintile than the highest, which is backwards from what statistical
# power predicts and is the signature of an unstable estimator rather than more biology.
#
# This script measures that directly. For every isoform it computes the within-condition
# standard deviation of IF across replicates -- pure noise, since condition is held fixed --
# and relates it to abundance. The floor is then chosen so that the |dIF| threshold sits a
# stated number of noise SDs above zero, rather than being picked by eye.
#
# It also compares baseline expression between the genotypes sharing a reference, since a
# systematic difference there would confound any "genotype X switches less" claim.

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

processed_dir <- resolve_path(paths$processed_dir %||% "data/processed", root = root)
results_tables <- ensure_dir(resolve_path(paths$results_tables %||% "results/tables", root = root))
results_figures <- ensure_dir(resolve_path(paths$results_figures %||% "results/figures", root = root))
fig_dir <- ensure_dir(file.path(results_figures, "expression"))

min_abs_dif <- as.numeric((cfg$significance %||% list())$min_abs_dif %||% 0.15)
# The floor is chosen so the effect threshold is at least this many replicate-noise SDs.
target_snr <- as.numeric(analysis_param(cfg, "expression_floor_snr", 3))

theme_set(theme_bw(base_size = 11))

sample_condition <- function(nms) {
  ifelse(grepl("_minus", nms, fixed = TRUE), "cond1",
    ifelse(grepl("_plus", nms, fixed = TRUE), "cond2", NA_character_)
  )
}

#' Within-condition SD of isoform fraction across replicates.
#' Condition is fixed, so this is replicate noise, not signal.
#'
#' Noise is related to GENE expression, not isoform expression. IF is isoform/gene, so it is
#' the gene's total abundance that determines how precisely the fraction can be estimated.
#' Keying this to isoform expression gives a misleading answer: a barely-expressed isoform
#' has a near-zero IF and therefore a tiny absolute SD, which looks stable but simply means
#' it cannot move.
if_noise <- function(rep_if) {
  if (is.null(rep_if)) return(NULL)
  value_cols <- setdiff(names(rep_if), "isoform_id")
  cond <- sample_condition(value_cols)
  c1 <- value_cols[cond %in% "cond1"]
  c2 <- value_cols[cond %in% "cond2"]
  if (length(c1) < 2L || length(c2) < 2L) return(NULL)
  m1 <- as.matrix(rep_if[, c1, drop = FALSE])
  m2 <- as.matrix(rep_if[, c2, drop = FALSE])
  sd1 <- apply(m1, 1, stats::sd, na.rm = TRUE)
  sd2 <- apply(m2, 1, stats::sd, na.rm = TRUE)
  tibble(
    isoform_id = as.character(rep_if$isoform_id),
    if_sd_within = sqrt((sd1^2 + sd2^2) / 2),
    n_reps_per_condition = min(length(c1), length(c2))
  )
}

datasets <- names(cfg$datasets %||% list())
stab_rows <- list()
floor_rows <- list()
expr_rows <- list()

for (k in datasets) {
  ds <- cfg$datasets[[k]]
  label <- sanitize(as.character(ds$label %||% k)[1L])
  ctx <- load_context_table(processed_dir, label, "features")
  rep_if <- load_context_table(processed_dir, label, "rep_if")
  if (is.null(ctx) || is.null(rep_if)) {
    message("[", k, "] no context / replicate IF; skipping.")
    next
  }
  noise <- if_noise(rep_if)
  if (is.null(noise)) {
    message("[", k, "] could not infer conditions from replicate column names; skipping.")
    next
  }
  gene_expr <- ctx |>
    group_by(.data$gene_id) |>
    summarize(
      gene_expr = sum(.data$iso_value_1, na.rm = TRUE) + sum(.data$iso_value_2, na.rm = TRUE),
      .groups = "drop"
    )
  d <- ctx |>
    transmute(
      isoform_id = as.character(.data$isoform_id),
      gene_id = as.character(.data$gene_id),
      IF_overall = suppressWarnings(as.numeric(.data$IF_overall)),
      dIF = suppressWarnings(as.numeric(.data$dIF)),
      iso_q = suppressWarnings(as.numeric(.data$isoform_switch_q_value)),
      gene_q = suppressWarnings(as.numeric(.data$gene_switch_q_value))
    ) |>
    left_join(gene_expr, by = "gene_id") |>
    inner_join(noise, by = "isoform_id") |>
    filter(is.finite(.data$gene_expr), is.finite(.data$if_sd_within))
  d$is_switching <- is.finite(d$iso_q) & d$iso_q < 0.05 &
    is.finite(d$gene_q) & d$gene_q < 0.05 &
    is.finite(d$dIF) & abs(d$dIF) >= min_abs_dif

  # Estimate noise only on isoforms that could plausibly show a threshold-sized shift.
  # Isoforms pinned near IF 0 or 1 contribute near-zero SD regardless of abundance and
  # would drag the estimate down without being informative.
  mid <- d |> filter(.data$gene_expr > 0, .data$IF_overall >= 0.05, .data$IF_overall <= 0.95)
  brk <- unique(stats::quantile(mid$gene_expr, probs = seq(0, 1, 0.1), na.rm = TRUE))
  mid$bin <- cut(mid$gene_expr, breaks = brk, include.lowest = TRUE)
  stab <- mid |>
    group_by(.data$bin) |>
    summarize(
      dataset = k,
      n = dplyr::n(),
      expr_min = min(.data$gene_expr),
      expr_median = stats::median(.data$gene_expr),
      if_sd_median = stats::median(.data$if_sd_within, na.rm = TRUE),
      pct_switching = 100 * mean(.data$is_switching, na.rm = TRUE),
      .groups = "drop"
    ) |>
    mutate(snr_at_threshold = min_abs_dif / .data$if_sd_median)
  stab_rows[[k]] <- stab

  # Floor: the lowest expression at which the |dIF| threshold is >= target_snr noise SDs,
  # and stays that way for all higher bins (so a single noisy bin cannot set it).
  ok <- stab$snr_at_threshold >= target_snr
  first_ok <- if (any(ok)) {
    idx <- which(ok)
    stable <- idx[vapply(idx, function(i) all(ok[i:length(ok)]), NA)]
    if (length(stable)) min(stable) else NA_integer_
  } else {
    NA_integer_
  }
  floor_val <- if (is.na(first_ok)) NA_real_ else stab$expr_min[first_ok]
  n_lost <- sum(d$is_switching & (d$gene_expr < (floor_val %||% Inf)), na.rm = TRUE)
  floor_rows[[k]] <- tibble(
    dataset = k,
    min_abs_dif = min_abs_dif,
    target_snr = target_snr,
    recommended_gene_expr_floor = floor_val,
    if_sd_at_floor = if (is.na(first_ok)) NA_real_ else stab$if_sd_median[first_ok],
    n_isoforms_total = nrow(d),
    n_isoforms_below_floor = sum(d$gene_expr < (floor_val %||% Inf), na.rm = TRUE),
    pct_isoforms_below_floor = 100 * mean(d$gene_expr < (floor_val %||% Inf), na.rm = TRUE),
    n_switching_calls = sum(d$is_switching, na.rm = TRUE),
    n_switching_calls_below_floor = n_lost,
    pct_switching_calls_below_floor = if (sum(d$is_switching, na.rm = TRUE)) {
      100 * n_lost / sum(d$is_switching, na.rm = TRUE)
    } else {
      NA_real_
    }
  )

  expr_rows[[k]] <- ctx |>
    group_by(.data$gene_id) |>
    summarize(
      gene_expr = sum(.data$iso_value_1, na.rm = TRUE) + sum(.data$iso_value_2, na.rm = TRUE),
      .groups = "drop"
    ) |>
    mutate(dataset = k, reference = as.character(ds$reference_transcriptome %||% NA))
  message(
    "[", k, "] floor = ",
    if (is.na(floor_val)) "none found" else signif(floor_val, 3),
    "  (drops ", floor_rows[[k]]$pct_switching_calls_below_floor |> round(1),
    "% of switching calls)"
  )
}

if (!length(stab_rows)) stop("No datasets had both context and replicate IF tables.")

stability <- bind_rows(stab_rows)
write_table_pair(stability, results_tables, "expression_if_stability", cfg = cfg)
floors <- bind_rows(floor_rows)

# A single floor is applied pipeline-wide rather than one per dataset: differing thresholds
# would make the datasets non-comparable, which is the class of problem this work has been
# removing. Per-dataset values scatter mostly because decile boundaries fall in different
# places, so the consensus is the median rounded to a round number.
consensus <- stats::median(floors$recommended_gene_expr_floor, na.rm = TRUE)
consensus_round <- if (is.na(consensus)) NA_real_ else round(consensus / 2) * 2
floors <- floors |> mutate(consensus_floor_applied = consensus_round)
write_table_pair(floors, results_tables, "expression_floor_recommendation", cfg = cfg)
print(as.data.frame(floors))
message(
  "\nPer-dataset floors: ",
  paste(floors$dataset, signif(floors$recommended_gene_expr_floor, 3), sep = "=", collapse = ", "),
  "\nConsensus floor to apply (config significance.min_gene_expression): ", consensus_round,
  "\nCurrently configured: ",
  (cfg$significance %||% list())$min_gene_expression %||% "not set"
)

# ---- Between-group baseline expression ------------------------------------------------
# A systematic abundance difference between genotypes sharing a reference would confound
# any "this genotype switches less" statement, so it is measured rather than assumed.
expr_all <- bind_rows(expr_rows)
cmp_rows <- list()
for (ref in unique(stats::na.omit(expr_all$reference))) {
  ds_in_ref <- unique(expr_all$dataset[expr_all$reference == ref])
  if (length(ds_in_ref) != 2L) next
  a <- expr_all |> filter(.data$dataset == ds_in_ref[1]) |> select(.data$gene_id, A = .data$gene_expr)
  b <- expr_all |> filter(.data$dataset == ds_in_ref[2]) |> select(.data$gene_id, B = .data$gene_expr)
  m <- inner_join(a, b, by = "gene_id") |> filter(is.finite(.data$A), is.finite(.data$B))
  if (!nrow(m)) next
  lr <- log2((m$B + 0.1) / (m$A + 0.1))
  w <- suppressWarnings(stats::wilcox.test(m$B, m$A, paired = TRUE))
  cmp_rows[[ref]] <- tibble(
    reference = ref,
    dataset_A = ds_in_ref[1], dataset_B = ds_in_ref[2],
    n_genes_in_both = nrow(m),
    median_expr_A = stats::median(m$A), median_expr_B = stats::median(m$B),
    median_log2_B_over_A = stats::median(lr),
    pct_genes_lower_in_B = 100 * mean(m$B < m$A),
    paired_wilcox_p = w$p.value,
    note = "Large n makes p uninformative; read the median log2 ratio."
  )
}
if (length(cmp_rows)) {
  between <- bind_rows(cmp_rows)
  write_table_pair(between, results_tables, "expression_between_groups", cfg = cfg)
  print(as.data.frame(between))
}

# ---- Figures ---------------------------------------------------------------------------
p1 <- ggplot(stability, aes(x = .data$expr_median, y = .data$if_sd_median, colour = .data$dataset)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.6) +
  geom_hline(yintercept = min_abs_dif / target_snr, linetype = 2, colour = "grey30") +
  scale_x_log10() +
  labs(
    title = "Isoform-fraction noise falls with abundance",
    subtitle = paste0(
      "Within-condition SD of IF across replicates (condition fixed, so this is noise). ",
      "Dashed line = |dIF| threshold / ", target_snr
    ),
    x = "Gene expression (median of bin, log10)", y = "Median within-condition IF SD"
  )
ggsave(file.path(fig_dir, "fig_expr_if_noise.png"), p1, width = 8, height = 5, dpi = 200)

p2 <- ggplot(stability, aes(x = .data$expr_median, y = .data$pct_switching, colour = .data$dataset)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.6) +
  scale_x_log10() +
  labs(
    title = "Switching calls concentrate at low abundance",
    subtitle = "Backwards from what statistical power predicts -- the signature of an unstable ratio",
    x = "Gene expression (median of bin, log10)", y = "% of isoforms called switching"
  )
ggsave(file.path(fig_dir, "fig_expr_switching_rate.png"), p2, width = 8, height = 5, dpi = 200)

if (exists("between") && nrow(between)) {
  p3 <- ggplot(expr_all |> filter(.data$gene_expr > 0), aes(x = .data$gene_expr, colour = .data$dataset)) +
    stat_ecdf(linewidth = 0.8) +
    scale_x_log10() +
    facet_wrap(~ .data$reference, scales = "free_x") +
    labs(
      title = "Baseline gene expression by dataset",
      subtitle = "Near-identical distributions within a reference argue against a global power difference",
      x = "Gene expression (log10)", y = "ECDF"
    )
  ggsave(file.path(fig_dir, "fig_expr_between_groups.png"), p3, width = 9, height = 4.5, dpi = 200)
}

write_run_manifest("02d_expression_diagnostics.R", cfg, root)
message("02d_expression_diagnostics.R: done")
