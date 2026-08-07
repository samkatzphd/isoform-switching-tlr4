#!/usr/bin/env Rscript
# 09: What does the splicing layer see that differential expression is blind to?
#
# Run (after 01): Rscript scripts/09_hidden_layer.R
# Works from committed data/processed/; the external drive is not needed.
#
# TWO CASES, both posed as "gene-level expression is misleading here":
#
#   CRYPTIC SWITCHERS -- the gene's total output barely moves, but the isoform composition
#   shifts substantially. Invisible to differential expression by construction.
#
#   COMPOSITIONAL RESCUE -- the gene's total output falls, but the fall is carried by minor
#   isoforms while the baseline-dominant one holds or rises. Read as "down" by differential
#   expression, when the functional transcript is being concentrated rather than lost.
#
# TWO DESIGN POINTS THAT SHAPE THE ANALYSIS
#
# 1. Dominance is defined at BASELINE (unstimulated), never after stimulation. Defining the
#    dominant isoform post-LPS selects on the outcome and manufactures the second class.
#
# 2. Reproduction is checked against the OTHER reference transcriptome for the same
#    libraries (T_HT vs T_UT are the same six libraries). That is a genuine re-annotation
#    replication -- the transcript spaces differ, since HT = T+H transcripts and
#    UT = T+U transcripts -- but NOT a biological replication: same RNA, same cells. It
#    tests robustness to annotation, nothing more.

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
fig_dir <- ensure_dir(file.path(results_figures, "hidden_layer"))

expr_floor <- as.numeric((cfg$significance %||% list())$min_gene_expression %||% 0)
cryptic_lfc <- as.numeric(analysis_param(cfg, "cryptic_max_abs_lfc", 0.25))
drop_lfc <- as.numeric(analysis_param(cfg, "rescue_min_gene_drop", 0.3))
dominant_min_if <- as.numeric(analysis_param(cfg, "dominant_min_if", 0.3))
theme_set(theme_bw(base_size = 11))

lab_of <- function(k) sanitize(as.character((cfg$datasets[[k]] %||% list())$label %||% k)[1L])

gene_iso <- function(k) {
  src <- load_scoring_table(processed_dir, lab_of(k), quiet = TRUE)
  if (is.null(src)) return(NULL)
  score_isoforms(src, cfg, k, lab_of(k)) |>
    mutate(
      gene_lfc = log2((as.numeric(.data$gene_value_2) + 1) / (as.numeric(.data$gene_value_1) + 1)),
      iso_lfc = log2((as.numeric(.data$iso_value_2) + 1) / (as.numeric(.data$iso_value_1) + 1))
    ) |>
    filter(.data$gene_expression >= expr_floor)
}

classify <- function(z, k) {
  # --- cryptic switchers -------------------------------------------------------------
  genes <- z |>
    group_by(.data$gene_id, .data$gene_name) |>
    summarize(
      switching = any(.data$is_switching, na.rm = TRUE),
      max_abs_dif = max_abs_finite(.data$dIF_n[.data$is_switching]),
      gene_lfc = dplyr::first(.data$gene_lfc),
      gene_expression = dplyr::first(.data$gene_expression),
      .groups = "drop"
    )
  cryptic <- genes |>
    filter(.data$switching, abs(.data$gene_lfc) < cryptic_lfc) |>
    mutate(dataset = k, class = "cryptic switcher")

  # --- compositional rescue ----------------------------------------------------------
  # Dominant isoform chosen on the UNSTIMULATED profile only.
  dom <- z |>
    group_by(.data$gene_id) |>
    filter(is.finite(.data$IF1)) |>
    slice_max(.data$IF1, n = 1, with_ties = FALSE) |>
    ungroup() |>
    select(.data$gene_id, dom_iso = .data$isoform_id, dom_IF1 = .data$IF1,
           dom_IF2 = .data$IF2, dom_lfc = .data$iso_lfc)
  minor <- z |>
    left_join(dom |> select(.data$gene_id, .data$dom_iso), by = "gene_id") |>
    filter(.data$isoform_id != .data$dom_iso) |>
    group_by(.data$gene_id) |>
    summarize(
      minor_expr_1 = sum(as.numeric(.data$iso_value_1), na.rm = TRUE),
      minor_expr_2 = sum(as.numeric(.data$iso_value_2), na.rm = TRUE),
      .groups = "drop"
    ) |>
    mutate(minor_lfc = log2((.data$minor_expr_2 + 1) / (.data$minor_expr_1 + 1)))

  rescue <- genes |>
    inner_join(dom, by = "gene_id") |>
    inner_join(minor, by = "gene_id") |>
    filter(
      .data$gene_lfc < -drop_lfc,          # gene output falls
      .data$dom_IF1 >= dominant_min_if,    # there is a real dominant isoform at baseline
      .data$dom_lfc > -0.1,                # dominant holds or rises
      .data$minor_lfc < .data$gene_lfc     # minor isoforms fall harder than the gene
    ) |>
    mutate(dataset = k, class = "compositional rescue")

  list(genes = genes |> mutate(dataset = k), cryptic = cryptic, rescue = rescue)
}

res <- lapply(c("T_UT", "U_UT", "T_HT", "H_HT"), function(k) {
  z <- gene_iso(k)
  if (is.null(z)) return(NULL)
  classify(z, k)
})
res <- res[!vapply(res, is.null, NA)]
names(res) <- vapply(res, function(r) r$genes$dataset[1], character(1))

all_genes <- bind_rows(lapply(res, `[[`, "genes"))
cryptic <- bind_rows(lapply(res, `[[`, "cryptic"))
rescue <- bind_rows(lapply(res, `[[`, "rescue"))

summary_tbl <- bind_rows(lapply(names(res), function(k) {
  g <- res[[k]]$genes
  tibble(
    dataset = k,
    n_genes_tested = nrow(g),
    n_switching = sum(g$switching, na.rm = TRUE),
    n_cryptic = nrow(res[[k]]$cryptic),
    pct_of_switchers_cryptic = if (sum(g$switching, na.rm = TRUE)) {
      100 * nrow(res[[k]]$cryptic) / sum(g$switching, na.rm = TRUE)
    } else NA_real_,
    n_genes_dropping = sum(g$gene_lfc < -drop_lfc, na.rm = TRUE),
    n_rescue = nrow(res[[k]]$rescue),
    pct_of_dropping_rescue = if (sum(g$gene_lfc < -drop_lfc, na.rm = TRUE)) {
      100 * nrow(res[[k]]$rescue) / sum(g$gene_lfc < -drop_lfc, na.rm = TRUE)
    } else NA_real_
  )
}))
write_table_pair(summary_tbl, results_tables, "hidden_layer_summary", cfg = cfg)
write_table_pair(
  cryptic |> arrange(.data$dataset, desc(.data$max_abs_dif)) |>
    select(.data$dataset, .data$gene_name, .data$gene_id, .data$gene_lfc,
           .data$max_abs_dif, .data$gene_expression),
  results_tables, "hidden_layer_cryptic_switchers", cfg = cfg
)
write_table_pair(
  rescue |> arrange(.data$dataset, .data$gene_lfc) |>
    select(.data$dataset, .data$gene_name, .data$gene_id, .data$gene_lfc,
           .data$dom_lfc, .data$minor_lfc, .data$dom_IF1, .data$dom_IF2, .data$gene_expression),
  results_tables, "hidden_layer_compositional_rescue", cfg = cfg
)
print(as.data.frame(summary_tbl))

# ---- Re-annotation reproduction ----------------------------------------------------------
# T_HT and T_UT are the same libraries under different transcript spaces. Reproduction here
# means robust to annotation, NOT biologically replicated.
reproduce <- function(cls, a = "T_UT", b = "T_HT") {
  if (!all(c(a, b) %in% names(res))) return(NULL)
  set_a <- cls |> filter(.data$dataset == a) |> pull(.data$gene_name) |> unique()
  set_b <- cls |> filter(.data$dataset == b) |> pull(.data$gene_name) |> unique()
  testable_b <- res[[b]]$genes$gene_name
  eligible <- intersect(set_a[is_real_gene_symbol(set_a)], testable_b)
  base_rate <- length(intersect(set_b, testable_b)) / length(unique(testable_b))
  hit <- sum(eligible %in% set_b)
  tibble(
    class_from = a, checked_in = b,
    n_eligible = length(eligible),
    n_reproduced = hit,
    pct_reproduced = if (length(eligible)) 100 * hit / length(eligible) else NA_real_,
    base_rate_pct = 100 * base_rate,
    enrichment_over_base = if (base_rate > 0 && length(eligible)) {
      (hit / length(eligible)) / base_rate
    } else NA_real_
  )
}
rep_tbl <- bind_rows(
  reproduce(cryptic) |> mutate(class = "cryptic switcher"),
  reproduce(rescue) |> mutate(class = "compositional rescue")
)
if (nrow(rep_tbl)) {
  write_table_pair(rep_tbl, results_tables, "hidden_layer_reproduction", cfg = cfg)
  print(as.data.frame(rep_tbl))
}

# ---- Figures -------------------------------------------------------------------------------
gl <- all_genes |> filter(.data$dataset %in% c("T_UT", "U_UT"))
p1 <- ggplot(gl, aes(x = .data$gene_lfc, y = .data$max_abs_dif)) +
  geom_point(data = gl |> filter(!.data$switching), colour = "grey80", alpha = 0.25, size = 0.5) +
  geom_point(data = gl |> filter(.data$switching), aes(colour = abs(.data$gene_lfc) < cryptic_lfc),
             alpha = 0.75, size = 1.3) +
  geom_vline(xintercept = c(-cryptic_lfc, cryptic_lfc), linetype = 2, colour = "grey30") +
  scale_colour_manual(values = c("TRUE" = "#e6550d", "FALSE" = "#3182bd"),
                      labels = c("TRUE" = "cryptic (gene barely moves)", "FALSE" = "also changes expression"),
                      name = NULL) +
  facet_wrap(~ .data$dataset) +
  coord_cartesian(xlim = c(-4, 4)) +
  labs(
    title = "Isoform switching is partly invisible to differential expression",
    subtitle = paste0("Orange: switching genes whose total output moves less than ",
                      cryptic_lfc, " log2 units"),
    x = "Gene-level LPS response (log2 fold-change)", y = "max |dIF| among switching isoforms"
  )
ggsave(file.path(fig_dir, "fig_cryptic_switchers.png"), p1, width = 9.5, height = 5, dpi = 200)

if (nrow(rescue)) {
  rl <- rescue |> filter(.data$dataset %in% c("T_UT", "U_UT"))
  if (nrow(rl)) {
    p2 <- ggplot(rl, aes(x = .data$gene_lfc, y = .data$dom_lfc)) +
      geom_hline(yintercept = 0, linetype = 2, colour = "grey40") +
      geom_abline(slope = 1, intercept = 0, linetype = 3, colour = "grey60") +
      geom_point(aes(size = .data$dom_IF2 - .data$dom_IF1), colour = "#1b7837", alpha = 0.75) +
      scale_size_continuous(name = "gain in dominant IF") +
      facet_wrap(~ .data$dataset) +
      labs(
        title = "Genes called 'down' whose dominant isoform is not going down",
        subtitle = "Dotted line = dominant isoform tracking the gene; points above it are concentrating, not lost",
        x = "Gene-level log2 fold-change", y = "Dominant isoform log2 fold-change"
      )
    ggsave(file.path(fig_dir, "fig_compositional_rescue.png"), p2, width = 9.5, height = 5, dpi = 200)
  }
}

p3 <- ggplot(summary_tbl, aes(x = .data$dataset)) +
  geom_col(aes(y = .data$pct_of_switchers_cryptic), fill = "#e6550d",
           colour = "grey25", width = 0.6) +
  geom_text(aes(y = .data$pct_of_switchers_cryptic,
                label = paste0(.data$n_cryptic, "/", .data$n_switching)),
            vjust = -0.35, size = 3.2) +
  labs(
    title = "Share of switching genes that differential expression would miss",
    subtitle = paste0("|gene log2FC| < ", cryptic_lfc),
    x = NULL, y = "% of switching genes"
  )
ggsave(file.path(fig_dir, "fig_cryptic_share.png"), p3, width = 7.5, height = 4.5, dpi = 200)

write_run_manifest("09_hidden_layer.R", cfg, root)
message("09_hidden_layer.R: done")
