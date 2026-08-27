#!/usr/bin/env Rscript
# 03: Novel (PacBio) isoform analysis
#     - Compare novel vs known isoforms at isoform and gene level
#     - Summarize class_code / effect sizes among novel switching events
#     - Rank top novel-involved switching genes
#     - Write tables + figures under results/
#
# Run (after 01 + 02): Rscript scripts/03_novel_isoform_analysis.R

bt <- c("utils/bootstrap.R", file.path("..", "utils", "bootstrap.R"))
b_file <- if (any(f <- vapply(bt, file.exists, NA))) { bt[which(f)[1L]] } else { NA_character_ }
if (is.na(b_file)) {
  stop("Could not find utils/bootstrap.R. Run from project root or set ISOFORM_PROJECT_ROOT.")
}
source(b_file, local = FALSE, chdir = FALSE)

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Install ggplot2: install.packages('ggplot2')")
  }
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
fig_dir <- ensure_dir(file.path(results_figures, "novel"))

sig <- cfg$significance %||% list()
iso_q_cutoff <- as.numeric(sig$isoform_q %||% 0.05)
gene_q_cutoff <- as.numeric(sig$gene_q %||% 0.05)
min_abs_dif <- as.numeric(sig$min_abs_dif %||% 0.0)
top_n <- as.integer(analysis_param(cfg, "top_n_genes", 30L))

theme_set(
  theme_bw(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold", size = 12),
      strip.text = element_text(face = "bold")
    )
)

# sanitize(), min_finite(), max_abs_finite(), score_isoforms(), is_real_gene_symbol()
# and pick_gene_symbol() now come from utils/helper_functions.R. This script used to
# carry its own copies; they had drifted from the versions in 04/06 (notably a weaker
# gene-symbol rule that let ENSG identifiers through as symbols, and no symbol
# collapsing, which split genes across rows).

iso_files <- list.files(
  processed_dir,
  pattern = "^isoformFeatures_.*\\.rds$",
  full.names = TRUE
)
if (!length(iso_files)) {
  stop("No isoformFeatures_*.rds in ", processed_dir, ". Run scripts/01_load_data.R first.")
}

scored_list <- list()
dataset_summary_rows <- list()
class_code_rows <- list()
gene_novel_rows <- list()
top_iso_rows <- list()
top_gene_rows <- list()

for (rds in iso_files) {
  label <- sub("^isoformFeatures_(.+)\\.rds$", "\\1", basename(rds), perl = TRUE, ignore.case = TRUE)
  message("[", label, "] Scoring isoforms ...")
  src <- load_scoring_table(processed_dir, label)
  if (is.null(src)) next
  iso <- score_isoforms(src, cfg, dataset_key = label, dataset_label = label)
  scored_list[[label]] <- iso

  n_iso <- nrow(iso)
  n_novel <- sum(iso$is_novel, na.rm = TRUE)
  n_switching <- sum(iso$is_switching, na.rm = TRUE)
  n_novel_switching <- sum(iso$is_novel & iso$is_switching, na.rm = TRUE)

  dataset_summary_rows[[label]] <- tibble(
    dataset_label = label,
    n_isoforms = n_iso,
    n_novel_isoforms = n_novel,
    pct_novel_isoforms = if (n_iso > 0L) 100 * n_novel / n_iso else NA_real_,
    n_switching_isoforms = n_switching,
    n_novel_switching_isoforms = n_novel_switching,
    pct_switching_that_are_novel = if (n_switching > 0L) {
      100 * n_novel_switching / n_switching
    } else {
      NA_real_
    },
    median_abs_dIF_novel = suppressWarnings(median(iso$abs_dIF[iso$is_novel], na.rm = TRUE)),
    median_abs_dIF_known = suppressWarnings(median(iso$abs_dIF[!iso$is_novel], na.rm = TRUE)),
    median_abs_dIF_novel_switching = suppressWarnings(
      median(iso$abs_dIF[iso$is_novel & iso$is_switching], na.rm = TRUE)
    ),
    median_abs_dIF_known_switching = suppressWarnings(
      median(iso$abs_dIF[!iso$is_novel & iso$is_switching], na.rm = TRUE)
    ),
    isoform_q_cutoff = iso_q_cutoff,
    gene_q_cutoff = gene_q_cutoff,
    min_abs_dif_cutoff = min_abs_dif
  )

  class_code_rows[[label]] <- iso |>
    mutate(
      class_code = ifelse(
        is.na(.data$class_code) | !nzchar(as.character(.data$class_code)),
        "NA",
        as.character(.data$class_code)
      ),
      novelty = ifelse(.data$is_novel, "novel_PB", "known")
    ) |>
    count(.data$dataset_label, .data$novelty, .data$class_code, .data$is_switching, name = "n") |>
    arrange(.data$novelty, desc(.data$n))

  gene_sum <- iso |>
    group_by(.data$gene_id, .data$gene_name) |>
    summarize(
      n_isoforms = n(),
      n_novel_isoforms = sum(.data$is_novel, na.rm = TRUE),
      n_switching_isoforms = sum(.data$is_switching, na.rm = TRUE),
      n_novel_switching_isoforms = sum(.data$is_novel & .data$is_switching, na.rm = TRUE),
      n_known_switching_isoforms = sum((!.data$is_novel) & .data$is_switching, na.rm = TRUE),
      novel_involved = any(.data$is_novel & .data$is_switching, na.rm = TRUE),
      has_novel_isoform = any(.data$is_novel, na.rm = TRUE),
      max_abs_dif_switching = max_abs_finite(.data$dIF_n[.data$is_switching]),
      max_abs_dif_novel_switching = max_abs_finite(.data$dIF_n[.data$is_novel & .data$is_switching]),
      min_isoform_switch_q = min_finite(.data$q_i[.data$is_switching]),
      min_novel_switch_q = min_finite(.data$q_i[.data$is_novel & .data$is_switching]),
      .groups = "drop"
    ) |>
    mutate(dataset_label = label)

  gene_novel_rows[[label]] <- gene_sum

  top_novel_iso <- iso |>
    filter(.data$is_novel, .data$is_switching) |>
    arrange(desc(.data$abs_dIF), .data$q_i)
  top_iso_rows[[label]] <- top_novel_iso |>
    slice_head(n = min(top_n, nrow(top_novel_iso))) |>
    transmute(
      dataset_label = .data$dataset_label,
      gene_id = .data$gene_id,
      gene_name = .data$gene_name,
      isoform_id = .data$isoform_id,
      oId = .data$oId,
      class_code = .data$class_code,
      cmp_ref = .data$cmp_ref,
      dIF = .data$dIF_n,
      abs_dIF = .data$abs_dIF,
      isoform_switch_q = .data$q_i,
      gene_switch_q = .data$q_g
    )

  top_novel_gene <- gene_sum |>
    filter(.data$novel_involved) |>
    arrange(desc(.data$max_abs_dif_novel_switching), desc(.data$n_novel_switching_isoforms), .data$min_novel_switch_q) |>
    mutate(rank = row_number())
  top_gene_rows[[label]] <- top_novel_gene |>
    slice_head(n = min(top_n, nrow(top_novel_gene)))
}

dataset_summary <- bind_rows(dataset_summary_rows)
class_code_summary <- bind_rows(class_code_rows)
gene_novel <- bind_rows(gene_novel_rows)
top_novel_isoforms <- bind_rows(top_iso_rows)
top_novel_genes <- bind_rows(top_gene_rows)
all_iso <- bind_rows(scored_list)

# ---- Write tables ----
write_out <- function(x, stem, csv = NULL) {
  write_table_pair(x, results_tables, stem, cfg = cfg, csv = csv)
}

write_out(class_code_summary, "novel_isoform_class_code_counts")
write_out(gene_novel, "novel_isoform_gene_summary")
write_out(top_novel_isoforms, "top_novel_switching_isoforms")
write_out(top_novel_genes, "top_novel_switching_genes")

# ---- Figures ----
p_pct <- ggplot(dataset_summary, aes(x = .data$dataset_label, y = .data$pct_novel_isoforms)) +
  geom_col(fill = "#8856a7", color = "grey25", linewidth = 0.2) +
  labs(
    title = "Fraction of isoforms tagged as novel PacBio (PB)",
    x = NULL, y = "% novel isoforms"
  ) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))
ggsave(file.path(fig_dir, "fig_novel_pct_isoforms_by_dataset.png"), p_pct, width = 7.2, height = 4.2, dpi = 200)

pct_sw <- dataset_summary |>
  transmute(
    dataset_label = .data$dataset_label,
    pct = .data$pct_switching_that_are_novel
  )
p_sw <- ggplot(pct_sw, aes(x = .data$dataset_label, y = .data$pct)) +
  geom_col(fill = "#6a51a3", color = "grey25", linewidth = 0.2) +
  labs(
    title = "Among switching isoforms, % that are novel PacBio",
    x = NULL, y = "% of switching isoforms"
  ) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))
ggsave(file.path(fig_dir, "fig_novel_pct_among_switching.png"), p_sw, width = 7.2, height = 4.2, dpi = 200)

iso_plot <- all_iso |>
  mutate(novelty = ifelse(.data$is_novel, "Novel (PB)", "Known"))

p_dif <- ggplot(
  iso_plot |> filter(is.finite(.data$abs_dIF)),
  aes(x = .data$novelty, y = .data$abs_dIF, fill = .data$novelty)
) +
  geom_violin(trim = TRUE, alpha = 0.55, color = "grey35") +
  geom_boxplot(width = 0.18, outlier.alpha = 0.25, alpha = 0.8) +
  scale_fill_manual(values = c("Known" = "#b3cde3", "Novel (PB)" = "#8856a7"), guide = "none") +
  facet_wrap(~ .data$dataset_label, scales = "free_y", ncol = 2) +
  labs(
    title = "Effect size (|dIF|) for novel vs known isoforms",
    x = NULL, y = "|dIF|"
  )
ggsave(file.path(fig_dir, "fig_novel_vs_known_abs_dif.png"), p_dif, width = 8, height = 6, dpi = 200)

sw_plot <- iso_plot |>
  filter(.data$is_switching, is.finite(.data$abs_dIF))
if (nrow(sw_plot) > 0L) {
  p_dif_sw <- ggplot(sw_plot, aes(x = .data$novelty, y = .data$abs_dIF, fill = .data$novelty)) +
    geom_violin(trim = TRUE, alpha = 0.55, color = "grey35") +
    geom_boxplot(width = 0.18, outlier.alpha = 0.25, alpha = 0.8) +
    scale_fill_manual(values = c("Known" = "#b3cde3", "Novel (PB)" = "#8856a7"), guide = "none") +
    facet_wrap(~ .data$dataset_label, scales = "free_y", ncol = 2) +
    labs(
      title = "Effect size among switching isoforms: novel vs known",
      x = NULL, y = "|dIF|"
    )
  ggsave(
    file.path(fig_dir, "fig_novel_vs_known_abs_dif_switching.png"),
    p_dif_sw, width = 8, height = 6, dpi = 200
  )
}

cc_plot <- class_code_summary |>
  filter(.data$novelty == "novel_PB") |>
  group_by(.data$dataset_label, .data$class_code) |>
  summarize(n = sum(.data$n), .groups = "drop")
if (nrow(cc_plot) > 0L) {
  p_cc <- ggplot(cc_plot, aes(x = reorder(.data$class_code, .data$n), y = .data$n, fill = .data$class_code)) +
    geom_col(color = "grey25", linewidth = 0.15, show.legend = FALSE) +
    coord_flip() +
    facet_wrap(~ .data$dataset_label, scales = "free_x", ncol = 2) +
    labs(
      title = "Novel PacBio isoforms by annotation class_code",
      x = "class_code", y = "Count"
    )
  ggsave(file.path(fig_dir, "fig_novel_class_code_counts.png"), p_cc, width = 8, height = 6, dpi = 200)
}

gene_plot <- gene_novel |>
  mutate(
    gene_class = case_when(
      .data$novel_involved ~ "Novel involved in switch",
      .data$has_novel_isoform ~ "Has novel isoform (not switching)",
      TRUE ~ "No novel isoform"
    )
  ) |>
  count(.data$dataset_label, .data$gene_class)

p_gene <- ggplot(gene_plot, aes(x = .data$dataset_label, y = .data$n, fill = .data$gene_class)) +
  geom_col(position = "fill", color = "grey30", linewidth = 0.2) +
  scale_y_continuous(labels = function(x) paste0(round(100 * x), "%")) +
  scale_fill_manual(
    values = c(
      "Novel involved in switch" = "#8856a7",
      "Has novel isoform (not switching)" = "#c994c7",
      "No novel isoform" = "#b3cde3"
    ),
    name = NULL
  ) +
  labs(
    title = "Gene-level novel isoform involvement",
    x = NULL, y = "Fraction of genes"
  ) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))
ggsave(file.path(fig_dir, "fig_novel_gene_involvement.png"), p_gene, width = 8, height = 4.5, dpi = 200)

# ---- HT vs UT novelty contrast ----
# Map processed labels back to config datasets / reference transcriptomes when possible.
label_to_ref <- list()
datasets_cfg <- cfg$datasets %||% list()
for (dk in names(datasets_cfg)) {
  ds <- datasets_cfg[[dk]] %||% list()
  lab <- sanitize(as.character(ds$label %||% dk)[1L])
  ref <- as.character(ds$reference_transcriptome %||% NA_character_)[1L]
  if (!is.na(ref) && nzchar(ref)) label_to_ref[[lab]] <- ref
}
# Fallback from label string if config map missed anything
for (lab in unique(dataset_summary$dataset_label)) {
  if (is.null(label_to_ref[[lab]])) {
    label_to_ref[[lab]] <- if (grepl("_HT|_ht|HT_", lab)) {
      "HT"
    } else if (grepl("_UT|_ut|UT_", lab)) {
      "UT"
    } else {
      "unknown"
    }
  }
}

dataset_summary <- dataset_summary |>
  mutate(
    reference_transcriptome = vapply(
      .data$dataset_label,
      function(x) label_to_ref[[x]] %||% "unknown",
      character(1)
    )
  )

ref_summary <- dataset_summary |>
  group_by(.data$reference_transcriptome) |>
  summarize(
    n_datasets = n(),
    n_isoforms = sum(.data$n_isoforms),
    n_novel_isoforms = sum(.data$n_novel_isoforms),
    pct_novel_isoforms = if (sum(.data$n_isoforms) > 0L) {
      100 * sum(.data$n_novel_isoforms) / sum(.data$n_isoforms)
    } else {
      NA_real_
    },
    n_switching_isoforms = sum(.data$n_switching_isoforms),
    n_novel_switching_isoforms = sum(.data$n_novel_switching_isoforms),
    pct_switching_that_are_novel = if (sum(.data$n_switching_isoforms) > 0L) {
      100 * sum(.data$n_novel_switching_isoforms) / sum(.data$n_switching_isoforms)
    } else {
      NA_real_
    },
    .groups = "drop"
  )

write_out(dataset_summary, "novel_isoform_dataset_summary")
write_out(ref_summary, "novel_isoform_HT_vs_UT_summary")

# Gene-name overlap of novel-involved genes across reference transcriptomes
# (XLOC gene_ids are transcriptome-specific; gene symbols are comparable).
gene_novel_ref <- gene_novel |>
  mutate(
    reference_transcriptome = vapply(
      .data$dataset_label,
      function(x) label_to_ref[[x]] %||% "unknown",
      character(1)
    ),
    # Shared symbol rule: this used to filter only ^XLOC_, so ENSG/ENST identifiers
    # leaked into the HT-vs-UT "gene symbol" overlap as if they were symbols.
    gene_name_clean = ifelse(
      is_real_gene_symbol(.data$gene_name),
      as.character(.data$gene_name),
      NA_character_
    )
  ) |>
  filter(!is.na(.data$gene_name_clean), .data$novel_involved)

ht_novel_genes <- unique(gene_novel_ref$gene_name_clean[gene_novel_ref$reference_transcriptome == "HT"])
ut_novel_genes <- unique(gene_novel_ref$gene_name_clean[gene_novel_ref$reference_transcriptome == "UT"])
shared_novel_genes <- sort(intersect(ht_novel_genes, ut_novel_genes))
ht_only_novel <- sort(setdiff(ht_novel_genes, ut_novel_genes))
ut_only_novel <- sort(setdiff(ut_novel_genes, ht_novel_genes))

gene_overlap_summary <- tibble(
  n_novel_involved_genes_HT = length(ht_novel_genes),
  n_novel_involved_genes_UT = length(ut_novel_genes),
  n_shared_gene_names = length(shared_novel_genes),
  n_HT_only = length(ht_only_novel),
  n_UT_only = length(ut_only_novel),
  # Both inputs were reduced to significant switching genes before saving, so
  # "HT only" largely means "absent from the UT objects", not "tested in UT and not
  # novel-involved". Symbol overlap here is bounded by object retention.
  n_symbols_present_in_both_references = length(intersect(
    unique(gene_novel_ref$gene_name_clean[gene_novel_ref$reference_transcriptome == "HT"]),
    unique(gene_novel_ref$gene_name_clean[gene_novel_ref$reference_transcriptome == "UT"])
  )),
  overlap_is_bounded_by_object_retention = TRUE
)
write_out(gene_overlap_summary, "novel_involved_genes_HT_vs_UT_overlap_summary")

shared_detail <- gene_novel_ref |>
  filter(.data$gene_name_clean %in% shared_novel_genes) |>
  group_by(.data$gene_name_clean, .data$reference_transcriptome) |>
  summarize(
    n_datasets = dplyr::n_distinct(.data$dataset_label),
    max_abs_dif_novel_switching = max_abs_finite(.data$max_abs_dif_novel_switching),
    max_n_novel_switching = max(.data$n_novel_switching_isoforms, na.rm = TRUE),
    .groups = "drop"
  )
# Wide summary without requiring tidyr
shared_ht <- shared_detail |>
  filter(.data$reference_transcriptome == "HT") |>
  transmute(
    gene_name = .data$gene_name_clean,
    n_datasets_HT = .data$n_datasets,
    max_abs_dif_novel_switching_HT = .data$max_abs_dif_novel_switching,
    max_n_novel_switching_HT = .data$max_n_novel_switching
  )
shared_ut <- shared_detail |>
  filter(.data$reference_transcriptome == "UT") |>
  transmute(
    gene_name = .data$gene_name_clean,
    n_datasets_UT = .data$n_datasets,
    max_abs_dif_novel_switching_UT = .data$max_abs_dif_novel_switching,
    max_n_novel_switching_UT = .data$max_n_novel_switching
  )
shared_detail <- full_join(shared_ht, shared_ut, by = "gene_name") |>
  arrange(
    desc(.data$max_abs_dif_novel_switching_HT),
    desc(.data$max_abs_dif_novel_switching_UT)
  )
write_out(shared_detail, "novel_involved_genes_shared_HT_UT")

# Contrast figures
p_ref <- ggplot(ref_summary, aes(x = .data$reference_transcriptome, y = .data$pct_novel_isoforms, fill = .data$reference_transcriptome)) +
  geom_col(color = "grey25", linewidth = 0.2, show.legend = FALSE) +
  scale_fill_manual(values = c(HT = "#3182bd", UT = "#e6550d", unknown = "grey60")) +
  labs(
    title = "HT vs UT: fraction of isoforms tagged novel PacBio",
    x = "Reference transcriptome", y = "% novel isoforms"
  )
ggsave(file.path(fig_dir, "fig_novel_HT_vs_UT_pct_isoforms.png"), p_ref, width = 5.5, height = 4.2, dpi = 200)

p_ref_sw <- ggplot(ref_summary, aes(x = .data$reference_transcriptome, y = .data$pct_switching_that_are_novel, fill = .data$reference_transcriptome)) +
  geom_col(color = "grey25", linewidth = 0.2, show.legend = FALSE) +
  scale_fill_manual(values = c(HT = "#3182bd", UT = "#e6550d", unknown = "grey60")) +
  labs(
    title = "HT vs UT: % of switching isoforms that are novel",
    x = "Reference transcriptome", y = "% of switching isoforms"
  )
ggsave(file.path(fig_dir, "fig_novel_HT_vs_UT_pct_among_switching.png"), p_ref_sw, width = 5.5, height = 4.2, dpi = 200)

overlap_bar <- tibble(
  category = c("Shared gene names", "HT only", "UT only"),
  n = c(length(shared_novel_genes), length(ht_only_novel), length(ut_only_novel))
)
p_ov <- ggplot(overlap_bar, aes(x = .data$category, y = .data$n, fill = .data$category)) +
  geom_col(color = "grey25", linewidth = 0.2, show.legend = FALSE) +
  scale_fill_manual(values = c(
    "Shared gene names" = "#756bb1",
    "HT only" = "#3182bd",
    "UT only" = "#e6550d"
  )) +
  labs(
    title = "Novel-involved genes by gene symbol: HT vs UT",
    x = NULL, y = "Number of genes"
  )
ggsave(file.path(fig_dir, "fig_novel_HT_vs_UT_gene_overlap.png"), p_ov, width = 6.5, height = 4.2, dpi = 200)

write_run_manifest("03_novel_isoform_analysis.R", cfg, root)
message("03_novel_isoform_analysis.R: done")
message("  Tables: ", results_tables)
message("  Figures: ", fig_dir)
print(dataset_summary)
print(ref_summary)
print(gene_overlap_summary)
