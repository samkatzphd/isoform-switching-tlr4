#!/usr/bin/env Rscript
# 07: Interferon-stimulated gene (ISG) profile of isoform switching.
#
# Run (after 01): Rscript scripts/07_isg_analysis.R
#
# WHY A CURATED SET RATHER THAN GO ENRICHMENT
#
# Script 05 asks an open question -- is anything enriched? -- and GO's nested hierarchy
# scatters interferon biology across dozens of overlapping terms, so a real ISG signal is
# both diluted and double-counted. This script asks a closed question instead: take one
# curated gene set, and test it once per dataset against the same tested background.
#
# GENE SET PROVENANCE
#
# Default sets come from GO annotations in org.Hs.eg.db (offline and reproducible):
#   type I   = GO:0034340 (response to type I interferon) + GO:0060337 (signalling)
#   type II  = GO:0034341 (response to type II interferon)
# Read with keytype "GOALL" so descendant terms are included.
#
# CAVEAT worth carrying into any writeup: "response to type I interferon" is not the same
# thing as "induced by interferon". These sets include upstream signalling components
# (CDC37, FADD, HDAC4 ...) alongside classical effectors. A functional screen-derived list
# (Schoggins) or Interferome would be a tighter definition of an ISG. Point
# `isg.custom_set_path` at a CSV with a `gene_name` column (optionally `set`) to use one --
# nothing else in this script changes.

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
fig_dir <- ensure_dir(file.path(results_figures, "isg"))

isg_cfg <- cfg$isg %||% list()
custom_path <- isg_cfg$custom_set_path %||% NULL
theme_set(theme_bw(base_size = 11))

# ---- Build the gene sets ---------------------------------------------------------------
isg_sets <- list()
set_source <- NA_character_
if (!is.null(custom_path) && nzchar(as.character(custom_path)[1L]) &&
  !identical(as.character(custom_path)[1L], "null")) {
  p <- resolve_path(as.character(custom_path)[1L], root = root)
  if (!file.exists(p)) stop("isg.custom_set_path does not exist: ", p)
  cust <- utils::read.csv(p, stringsAsFactors = FALSE, check.names = FALSE)
  if (!"gene_name" %in% names(cust)) stop("Custom ISG file needs a `gene_name` column: ", p)
  if ("set" %in% names(cust)) {
    isg_sets <- split(unique(cust$gene_name), cust$set)
  } else {
    isg_sets <- list(custom = unique(cust$gene_name))
  }
  set_source <- paste0("custom file: ", basename(p))
} else {
  if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
    stop('Need org.Hs.eg.db: BiocManager::install("org.Hs.eg.db")')
  }
  go_syms <- function(ids) {
    s <- suppressMessages(AnnotationDbi::select(
      org.Hs.eg.db::org.Hs.eg.db, keys = ids, keytype = "GOALL", columns = "SYMBOL"
    )$SYMBOL)
    sort(unique(s[!is.na(s)]))
  }
  isg_sets <- list(
    `Type I IFN` = go_syms(c("GO:0034340", "GO:0060337")),
    `Type II IFN` = go_syms("GO:0034341")
  )
  isg_sets[["Any IFN"]] <- sort(unique(unlist(isg_sets, use.names = FALSE)))
  set_source <- "GO (org.Hs.eg.db, GOALL): GO:0034340 + GO:0060337; GO:0034341"
}
message("ISG set source: ", set_source)
for (nm in names(isg_sets)) message("  ", nm, ": ", length(isg_sets[[nm]]), " genes")

# ---- Per-dataset scoring ----------------------------------------------------------------
gene_table <- function(k) {
  ds <- cfg$datasets[[k]]
  label <- sanitize(as.character(ds$label %||% k)[1L])
  src <- load_scoring_table(processed_dir, label, quiet = TRUE)
  if (is.null(src)) return(NULL)
  z <- score_isoforms(src, cfg, k, label)
  z |>
    filter(is_real_gene_symbol(.data$gene_name)) |>
    group_by(.data$gene_name) |>
    summarize(
      switching = any(.data$is_switching, na.rm = TRUE),
      max_abs_dif = max_abs_finite(.data$dIF_n[.data$is_switching]),
      max_abs_dif_any = max_abs_finite(.data$dIF_n),
      gene_expression = suppressWarnings(max(.data$gene_expression, na.rm = TRUE)),
      novel_involved = any(.data$is_switching & .data$is_novel, na.rm = TRUE),
      .groups = "drop"
    ) |>
    mutate(dataset = k, reference = as.character(ds$reference_transcriptome %||% NA))
}

datasets <- names(cfg$datasets %||% list())
gt <- lapply(datasets, gene_table)
names(gt) <- datasets
gt <- gt[!vapply(gt, is.null, NA)]
if (!length(gt)) stop("No datasets could be scored. Run 01 first.")

# ---- Enrichment: is switching over-represented among ISGs? -------------------------------
# Universe is the tested gene set for that dataset, and the ISG set is intersected with it
# so that ISGs never measured cannot inflate or deflate the result.
enrich_rows <- list()
for (k in names(gt)) {
  d <- gt[[k]]
  for (nm in names(isg_sets)) {
    in_set <- d$gene_name %in% isg_sets[[nm]]
    tab <- matrix(
      c(sum(!in_set & !d$switching), sum(!in_set & d$switching),
        sum(in_set & !d$switching), sum(in_set & d$switching)),
      nrow = 2, byrow = TRUE,
      dimnames = list(ISG = c("FALSE", "TRUE"), switching = c("FALSE", "TRUE"))
    )
    ft <- tryCatch(stats::fisher.test(tab), error = function(e) NULL)
    enrich_rows[[length(enrich_rows) + 1L]] <- tibble(
      dataset = k,
      isg_set = nm,
      n_tested = nrow(d),
      n_isg_tested = sum(in_set),
      n_switching = sum(d$switching),
      n_isg_switching = sum(in_set & d$switching),
      pct_isg_switching = if (sum(in_set)) 100 * sum(in_set & d$switching) / sum(in_set) else NA_real_,
      pct_nonisg_switching = if (sum(!in_set)) 100 * sum(!in_set & d$switching) / sum(!in_set) else NA_real_,
      odds_ratio = if (!is.null(ft)) unname(ft$estimate) else NA_real_,
      conf_low = if (!is.null(ft)) ft$conf.int[[1L]] else NA_real_,
      conf_high = if (!is.null(ft)) ft$conf.int[[2L]] else NA_real_,
      p_value = if (!is.null(ft)) ft$p.value else NA_real_
    )
  }
}
enrich <- bind_rows(enrich_rows) |>
  group_by(.data$isg_set) |>
  mutate(p_adj_across_datasets = stats::p.adjust(.data$p_value, "BH")) |>
  ungroup() |>
  mutate(isg_set_source = set_source)
write_table_pair(enrich, results_tables, "isg_enrichment", cfg = cfg)
print(as.data.frame(
  enrich |> filter(.data$isg_set == "Any IFN" | length(isg_sets) == 1L) |>
    transmute(.data$dataset, .data$isg_set, .data$n_isg_tested, .data$n_isg_switching,
              pct_isg = round(.data$pct_isg_switching, 2),
              pct_other = round(.data$pct_nonisg_switching, 2),
              OR = round(.data$odds_ratio, 2), p = signif(.data$p_value, 3))
))

# ---- Which ISGs switch, and where -------------------------------------------------------
any_isg <- isg_sets[[length(isg_sets)]]
profile <- bind_rows(gt) |>
  filter(.data$gene_name %in% any_isg) |>
  mutate(
    isg_class = dplyr::case_when(
      .data$gene_name %in% (isg_sets[["Type I IFN"]] %||% character()) &
        .data$gene_name %in% (isg_sets[["Type II IFN"]] %||% character()) ~ "Type I + II",
      .data$gene_name %in% (isg_sets[["Type I IFN"]] %||% character()) ~ "Type I",
      .data$gene_name %in% (isg_sets[["Type II IFN"]] %||% character()) ~ "Type II",
      TRUE ~ "ISG"
    )
  )
write_table_pair(profile, results_tables, "isg_gene_profile", cfg = cfg)

switching_isgs <- profile |>
  filter(.data$switching) |>
  arrange(.data$dataset, desc(.data$max_abs_dif))
write_table_pair(switching_isgs, results_tables, "isg_switching_genes", cfg = cfg)

# ---- Paired profiles: T vs H (HT) and WT vs KO (UT) ---------------------------------------
pair_profile <- function(a, b, pair_name) {
  if (!all(c(a, b) %in% names(gt))) return(NULL)
  A <- gt[[a]]; B <- gt[[b]]
  both <- intersect(A$gene_name, B$gene_name)
  isg_both <- intersect(both, any_isg)
  sa <- A$gene_name[A$switching]; sb <- B$gene_name[B$switching]
  cls <- tibble(
    gene_name = isg_both,
    in_A = isg_both %in% sa,
    in_B = isg_both %in% sb
  ) |>
    mutate(
      status = dplyr::case_when(
        .data$in_A & .data$in_B ~ "Switching in both",
        .data$in_A ~ paste0("Switching in ", a, " only"),
        .data$in_B ~ paste0("Switching in ", b, " only"),
        TRUE ~ "Switching in neither"
      ),
      pair = pair_name, dataset_A = a, dataset_B = b
    ) |>
    left_join(A |> select(.data$gene_name, A_max_abs_dif = .data$max_abs_dif_any,
                          A_expr = .data$gene_expression), by = "gene_name") |>
    left_join(B |> select(.data$gene_name, B_max_abs_dif = .data$max_abs_dif_any,
                          B_expr = .data$gene_expression), by = "gene_name")
  cls
}
pairs <- bind_rows(
  pair_profile("T_HT", "H_HT", "HT: T vs H"),
  pair_profile("T_UT", "U_UT", "UT: WT vs UBL5 KO")
)
if (nrow(pairs)) {
  write_table_pair(pairs, results_tables, "isg_pair_profile", cfg = cfg)
  pair_summary <- pairs |>
    count(.data$pair, .data$dataset_A, .data$dataset_B, .data$status, name = "n_isgs") |>
    group_by(.data$pair) |>
    mutate(pct = 100 * .data$n_isgs / sum(.data$n_isgs)) |>
    ungroup()
  write_table_pair(pair_summary, results_tables, "isg_pair_summary", cfg = cfg)
  print(as.data.frame(pair_summary))
}

# ---- Effect size: do ISG switches differ in magnitude? ----------------------------------
eff <- bind_rows(gt) |>
  filter(.data$switching, is.finite(.data$max_abs_dif)) |>
  mutate(group = ifelse(.data$gene_name %in% any_isg, "ISG", "Other"))
eff_rows <- list()
for (k in unique(eff$dataset)) {
  d <- eff |> filter(.data$dataset == k)
  a <- d$max_abs_dif[d$group == "ISG"]; b <- d$max_abs_dif[d$group == "Other"]
  w <- if (length(a) >= 3L && length(b) >= 3L) {
    suppressWarnings(stats::wilcox.test(a, b))
  } else {
    NULL
  }
  eff_rows[[k]] <- tibble(
    dataset = k, n_isg = length(a), n_other = length(b),
    median_abs_dif_isg = suppressWarnings(stats::median(a, na.rm = TRUE)),
    median_abs_dif_other = suppressWarnings(stats::median(b, na.rm = TRUE)),
    wilcox_p = if (!is.null(w)) w$p.value else NA_real_
  )
}
eff_tbl <- bind_rows(eff_rows)
write_table_pair(eff_tbl, results_tables, "isg_effect_size", cfg = cfg)
print(as.data.frame(eff_tbl))

# ---- Figures -----------------------------------------------------------------------------
p1 <- ggplot(
  enrich |> filter(is.finite(.data$odds_ratio)),
  aes(x = .data$dataset, y = .data$odds_ratio, colour = .data$isg_set)
) +
  geom_hline(yintercept = 1, linetype = 2, colour = "grey40") +
  geom_pointrange(
    aes(ymin = .data$conf_low, ymax = pmin(.data$conf_high, 50)),
    position = position_dodge(width = 0.5), fatten = 2.5
  ) +
  scale_y_log10() +
  labs(
    title = "Is isoform switching over-represented among interferon-response genes?",
    subtitle = "Odds ratio vs the tested background of each dataset; CI capped at 50 for display",
    x = NULL, y = "Odds ratio (log scale)", colour = "Gene set"
  )
ggsave(file.path(fig_dir, "fig_isg_enrichment.png"), p1, width = 8.5, height = 5, dpi = 200)

pct_long <- enrich |>
  select(.data$dataset, .data$isg_set, .data$pct_isg_switching, .data$pct_nonisg_switching) |>
  tidyr::pivot_longer(c("pct_isg_switching", "pct_nonisg_switching"),
                      names_to = "group", values_to = "pct") |>
  mutate(group = ifelse(.data$group == "pct_isg_switching", "ISG", "Other genes"))
p2 <- ggplot(pct_long, aes(x = .data$dataset, y = .data$pct, fill = .data$group)) +
  geom_col(position = "dodge", colour = "grey25", linewidth = 0.2) +
  facet_wrap(~ .data$isg_set) +
  scale_fill_manual(values = c("ISG" = "#e6550d", "Other genes" = "#b3cde3"), name = NULL) +
  labs(
    title = "Percentage of tested genes with an isoform switch",
    x = NULL, y = "% of tested genes switching"
  ) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
ggsave(file.path(fig_dir, "fig_isg_switching_rate.png"), p2, width = 9, height = 4.5, dpi = 200)

if (nrow(pairs)) {
  p3 <- ggplot(
    pair_summary |> filter(.data$status != "Switching in neither"),
    aes(x = .data$status, y = .data$n_isgs, fill = .data$pair)
  ) +
    geom_col(colour = "grey25", linewidth = 0.2, show.legend = FALSE) +
    geom_text(aes(label = .data$n_isgs), vjust = -0.3, size = 3.2) +
    facet_wrap(~ .data$pair, scales = "free_x") +
    labs(
      title = "Interferon-response genes with isoform switches, by pair",
      subtitle = "Restricted to ISGs tested in both members of each pair",
      x = NULL, y = "ISGs"
    ) +
    theme(axis.text.x = element_text(angle = 20, hjust = 1))
  ggsave(file.path(fig_dir, "fig_isg_pair_overlap.png"), p3, width = 10, height = 5, dpi = 200)
}

if (nrow(eff)) {
  p4 <- ggplot(eff, aes(x = .data$group, y = .data$max_abs_dif, fill = .data$group)) +
    geom_violin(trim = TRUE, alpha = 0.55, colour = "grey35") +
    geom_boxplot(width = 0.18, outlier.alpha = 0.3, alpha = 0.85) +
    scale_fill_manual(values = c("ISG" = "#e6550d", "Other" = "#b3cde3"), guide = "none") +
    facet_wrap(~ .data$dataset, nrow = 1) +
    labs(
      title = "Effect size of switching genes: ISG versus the rest",
      x = NULL, y = "max |dIF| among switching isoforms"
    )
  ggsave(file.path(fig_dir, "fig_isg_effect_size.png"), p4, width = 9, height = 4.5, dpi = 200)
}

write_run_manifest("07_isg_analysis.R", cfg, root, extra = list(isg_set_source = set_source))
message("07_isg_analysis.R: done")
