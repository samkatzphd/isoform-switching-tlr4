#!/usr/bin/env Rscript
# 00: Explore how switching gene/isoform counts respond to min |dIF| after q-value
#     filtering. Does not change config defaults.
#
# Run (after 01): Rscript scripts/00_threshold_scan_abs_dif.R
#
# IMPORTANT INTERPRETATION NOTE
# The input ISA objects were saved after reduceToSwitchingGenes = TRUE, so every gene
# in data/processed already carries a significant gene-level q. This scan therefore
# describes how the |dIF| floor *prunes* an already-significant set -- it is not a
# sensitivity/specificity curve against a full tested background. See
# results/tables/input_object_reduction_check.csv and docs/REVIEW_CHANGES.md.

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
fig_dir <- ensure_dir(file.path(results_figures, "threshold_scan"))

iq <- as.numeric((cfg$significance %||% list())$isoform_q %||% 0.05)
gq <- as.numeric((cfg$significance %||% list())$gene_q %||% 0.05)
thresholds <- as.numeric(analysis_param(cfg, "threshold_scan", c(0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3)))

theme_set(theme_bw(base_size = 11))

# Score once at min_abs_dif = 0 (q filters only), then vary the |dIF| floor below.
# Gene symbols are collapsed by score_isoforms(), so gene counts here are gene counts,
# not row counts -- this script previously used its own ad-hoc symbol helpers.
cfg_q_only <- cfg
cfg_q_only$significance$min_abs_dif <- 0

score <- function(path, label) {
  z <- score_isoforms(readRDS(path), cfg_q_only, dataset_key = label, dataset_label = label)
  z$sig_q <- z$is_switching
  z$dataset <- label
  z
}

files <- list.files(processed_dir, pattern = "^isoformFeatures_.*\\.rds$", full.names = TRUE)
if (!length(files)) {
  stop("No isoformFeatures_*.rds in ", processed_dir, ". Run scripts/01_load_data.R first.")
}
all <- bind_rows(lapply(files, function(f) {
  score(f, sub("^isoformFeatures_(.+)\\.rds$", "\\1", basename(f)))
}))
sig <- all |> filter(.data$sig_q)

message("Significant isoforms by dataset (q filters only):")
print(count(sig, .data$dataset))

message("|dIF| quantiles among q-significant isoforms:")
print(
  sig |>
    group_by(.data$dataset) |>
    summarize(
      n = n(),
      q10 = quantile(.data$abs_dIF, 0.1, na.rm = TRUE),
      q25 = quantile(.data$abs_dIF, 0.25, na.rm = TRUE),
      q50 = quantile(.data$abs_dIF, 0.5, na.rm = TRUE),
      q75 = quantile(.data$abs_dIF, 0.75, na.rm = TRUE),
      q90 = quantile(.data$abs_dIF, 0.9, na.rm = TRUE),
      .groups = "drop"
    )
)

rows <- list()
for (th in thresholds) {
  for (lab in unique(sig$dataset)) {
    s <- sig |> filter(.data$dataset == lab, is.finite(.data$abs_dIF), .data$abs_dIF >= th)
    rows[[length(rows) + 1L]] <- tibble(
      dataset = lab,
      min_abs_dif = th,
      n_sig_isoforms = nrow(s),
      n_sig_genes = dplyr::n_distinct(s$gene_id)
    )
  }
}
scan <- bind_rows(rows)
print(scan, n = 100)
write_table_pair(scan, results_tables, "min_abs_dif_threshold_scan", cfg = cfg)

p1 <- ggplot(sig, aes(x = .data$abs_dIF)) +
  geom_histogram(bins = 40, fill = "#6a51a3", color = "white", linewidth = 0.1) +
  geom_vline(xintercept = c(0.1, 0.2), linetype = 2, color = "grey30") +
  facet_wrap(~ .data$dataset, scales = "free_y") +
  labs(
    title = "|dIF| among q-significant isoforms",
    subtitle = "Dashed lines mark |dIF| = 0.10 and 0.20; genes shown are already switching-significant",
    x = "|dIF|", y = "Count"
  )
ggsave(file.path(fig_dir, "fig_sig_abs_dif_hist.png"), p1, width = 9, height = 6, dpi = 180)

p2 <- ggplot(scan, aes(x = .data$min_abs_dif, y = .data$n_sig_genes, color = .data$dataset)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2) +
  labs(
    title = "Switching gene count vs min |dIF| (after q filter)",
    subtitle = "Pruning within an already-significant set, not a power curve",
    x = "min |dIF| threshold",
    y = "Genes with >= 1 switching isoform"
  )
ggsave(file.path(fig_dir, "fig_gene_count_vs_min_abs_dif.png"), p2, width = 8, height = 4.8, dpi = 180)

p3 <- ggplot(scan, aes(x = .data$min_abs_dif, y = .data$n_sig_isoforms, color = .data$dataset)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2) +
  labs(
    title = "Switching isoform count vs min |dIF| (after q filter)",
    x = "min |dIF| threshold",
    y = "Significant isoforms"
  )
ggsave(file.path(fig_dir, "fig_isoform_count_vs_min_abs_dif.png"), p3, width = 8, height = 4.8, dpi = 180)

# ---- UT overlap stability ----
# Symbols come from the shared rule in utils/ (XLOC_* and ENS* are not symbols).
t <- score(file.path(processed_dir, "isoformFeatures_T_UT.rds"), "T_UT")
u <- score(file.path(processed_dir, "isoformFeatures_U_UT.rds"), "U_UT")

ov <- list()
for (th in thresholds) {
  tg <- unique(t$gene_name[t$sig_q & is.finite(t$abs_dIF) & t$abs_dIF >= th & is_real_gene_symbol(t$gene_name)])
  ug <- unique(u$gene_name[u$sig_q & is.finite(u$abs_dIF) & u$abs_dIF >= th & is_real_gene_symbol(u$gene_name)])
  shared <- length(intersect(tg, ug))
  t_only <- length(setdiff(tg, ug))
  u_only <- length(setdiff(ug, tg))
  ov[[length(ov) + 1L]] <- tibble(
    min_abs_dif = th,
    shared = shared,
    T_only = t_only,
    U_only = u_only,
    n_T = shared + t_only,
    n_U = shared + u_only,
    jaccard = if ((shared + t_only + u_only) > 0L) shared / (shared + t_only + u_only) else NA_real_
  )
}
ov <- bind_rows(ov)
print(ov)
write_table_pair(ov, results_tables, "ut_overlap_vs_min_abs_dif", cfg = cfg)

ov_long <- bind_rows(
  ov |> transmute(min_abs_dif, class = "shared", n = .data$shared),
  ov |> transmute(min_abs_dif, class = "T_only", n = .data$T_only),
  ov |> transmute(min_abs_dif, class = "U_only", n = .data$U_only)
)
p4 <- ggplot(ov_long, aes(x = .data$min_abs_dif, y = .data$n, color = .data$class)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2) +
  labs(
    title = "UT gene overlap classes vs min |dIF|",
    subtitle = "T_only/U_only are dominated by genes absent from the other saved object",
    x = "min |dIF|", y = "Genes"
  )
ggsave(file.path(fig_dir, "fig_ut_overlap_vs_min_abs_dif.png"), p4, width = 7.5, height = 4.5, dpi = 180)

write_run_manifest("00_threshold_scan_abs_dif.R", cfg, root)
message("Wrote threshold scan tables/figures under results/")
