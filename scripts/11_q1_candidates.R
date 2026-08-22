#!/usr/bin/env Rscript
# 11: Candidate genes for Q1 -- the wildtype LPS switching response, on its own terms.
#
# Run (after 01): Rscript scripts/11_q1_candidates.R
# Works from committed data/processed/; the external drive is not needed.
#
# WHY THIS SCRIPT EXISTS
#
# The candidate set carried through the external review (IRAK3, SOCS4, RAB7B, SPRING1) was
# built for Q2: it starts from the 23 CRYPTIC switchers in T_UT and then hand-picks four
# genes for TLR-pathway relevance. Two consequences that this script exists to fix:
#
#   1. It was drawn from T_UT -- the wildtype arm of the Q2 experiment -- and never from
#      T_HT, which is the designated Q1 ANCHOR. Both are wildtype +/- LPS on the same six
#      libraries, so the biology is Q1 either way, but the anchor was never mined.
#
#   2. It is restricted to the cryptic window (|gene log2FC| < 0.25). That is the right
#      filter for "what would differential expression miss", and the wrong one for "what
#      does LPS do to isoform usage" -- it discards the largest switches in the dataset
#      (NCOA7, RSAD2, DDX58, IFIT1 ...) purely because the gene also changes in abundance.
#
# So: no ranking of Q1 candidates exists. This produces one.
#
# WHAT COUNTS AS SUPPORT HERE, AND WHAT DOES NOT
#
#   T_HT vs T_UT agreement is ANNOTATION ROBUSTNESS, not replication. They are the same six
#   libraries -- same RNA, same cells -- quantified against two merged transcriptomes
#   (HT = T+H PacBio, UT = T+U). A gene switching in both survived a change of transcript
#   space, which is a real and useful filter for assay design. It is NOT independent
#   evidence, and must never be described as reproducing in a second experiment.
#
#   H_HT agreement IS cross-genotype support -- different donors, primary cells, so a gene
#   switching there too is genuinely more interesting. But H_HT alone was batch-corrected
#   (REVIEW_CHANGES.md 0g), so its calls sit on a different scale. It is reported as a
#   supporting flag ONLY, never as a filter, and no H_HT effect size is compared to a
#   T effect size anywhere in this script.
#
#   Ranking by effect size selects on the outcome. The ranked table is descriptive; it
#   carries no test and no p-value, by design.

bt <- c("utils/bootstrap.R", file.path("..", "utils", "bootstrap.R"))
b_file <- if (any(f <- vapply(bt, file.exists, NA))) { bt[which(f)[1L]] } else { NA_character_ }
if (is.na(b_file)) {
  stop("Could not find utils/bootstrap.R. Run from project root or set ISOFORM_PROJECT_ROOT.")
}
source(b_file, local = FALSE, chdir = FALSE)

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
})
has_ggplot <- requireNamespace("ggplot2", quietly = TRUE)
if (has_ggplot) suppressPackageStartupMessages(library(ggplot2))

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
fig_dir <- ensure_dir(file.path(results_figures, "q1_candidates"))

expr_floor <- as.numeric((cfg$significance %||% list())$min_gene_expression %||% 0)
cryptic_lfc <- as.numeric(analysis_param(cfg, "cryptic_max_abs_lfc", 0.25))
shortlist_n <- as.integer(analysis_param(cfg, "q1_shortlist_n", 20))
if (has_ggplot) theme_set(theme_bw(base_size = 11))

lab_of <- function(k) sanitize(as.character((cfg$datasets[[k]] %||% list())$label %||% k)[1L])

# The Q1 arms: the anchor, its re-annotation, and the cross-genotype check.
ANCHOR <- "T_HT"
REANNOT <- "T_UT"
CROSS <- "H_HT"

gene_table <- function(k) {
  src <- load_scoring_table(processed_dir, lab_of(k), quiet = TRUE)
  if (is.null(src)) {
    message("  ", k, ": no scoring table; skipped")
    return(NULL)
  }
  score_isoforms(src, cfg, k, lab_of(k)) |>
    mutate(
      gene_lfc = log2((as.numeric(.data$gene_value_2) + 1) / (as.numeric(.data$gene_value_1) + 1))
    ) |>
    filter(.data$gene_expression >= expr_floor) |>
    group_by(.data$gene_id, .data$gene_name) |>
    summarize(
      switching = any(.data$is_switching, na.rm = TRUE),
      max_abs_dif = max_abs_finite(.data$dIF_n[.data$is_switching]),
      n_switching_isoforms = sum(.data$is_switching, na.rm = TRUE),
      gene_lfc = dplyr::first(.data$gene_lfc),
      gene_expression = dplyr::first(.data$gene_expression),
      .groups = "drop"
    ) |>
    filter(is_real_gene_symbol(.data$gene_name))
}

message("11_q1_candidates.R: scoring the wildtype LPS arms")
tabs <- setNames(lapply(c(ANCHOR, REANNOT, CROSS), gene_table), c(ANCHOR, REANNOT, CROSS))
if (is.null(tabs[[ANCHOR]])) {
  stop("The Q1 anchor (", ANCHOR, ") has no scoring table. Run 01_load_data.R first.")
}

# Gene symbols are the only cross-reference key: TCONS_* ids are assigned per reference.
sym_of <- function(k, switching_only = TRUE) {
  if (is.null(tabs[[k]])) return(character(0))
  d <- tabs[[k]]
  if (switching_only) d <- d[d$switching %in% TRUE, , drop = FALSE]
  unique(d$gene_name)
}

tested_anchor <- unique(tabs[[ANCHOR]]$gene_name)
tested_reannot <- if (is.null(tabs[[REANNOT]])) character(0) else unique(tabs[[REANNOT]]$gene_name)
tested_cross <- if (is.null(tabs[[CROSS]])) character(0) else unique(tabs[[CROSS]]$gene_name)
sw_reannot <- sym_of(REANNOT)
sw_cross <- sym_of(CROSS)

# One row per gene switching in EITHER wildtype annotation -- the Q1 response as a whole,
# not the cryptic subset.
anchor_sw <- tabs[[ANCHOR]] |> filter(.data$switching)
reannot_sw <- if (is.null(tabs[[REANNOT]])) NULL else tabs[[REANNOT]] |> filter(.data$switching)

cand <- tibble(
  gene_name = union(
    anchor_sw$gene_name,
    if (is.null(reannot_sw)) character(0) else reannot_sw$gene_name
  )
)

# Attach each dataset's effect columns by symbol. A symbol can carry more than one gene_id
# (the references disagree on locus boundaries), so collapse to the strongest row per
# symbol before joining.
#
# NB: coerce non-finite to NA over NUMERIC columns only. An earlier version used
# across(everything()), which ran is.finite() on gene_name, silently turned every symbol
# into NA, and left the joined effect columns empty while the shortlist -- computed from
# the symbol sets, not this join -- still looked right.
join_side <- function(d, k, prefix) {
  dif_col <- paste0("dif_", prefix)
  lfc_col <- paste0("lfc_", prefix)
  if (is.null(tabs[[k]])) {
    d[[dif_col]] <- NA_real_
    d[[lfc_col]] <- NA_real_
    return(d)
  }
  s <- tabs[[k]] |>
    transmute(
      gene_name = .data$gene_name,
      !!dif_col := as.numeric(.data$max_abs_dif),
      !!lfc_col := as.numeric(.data$gene_lfc)
    ) |>
    group_by(.data$gene_name) |>
    summarize(across(where(is.numeric), ~ {
      v <- .x[is.finite(.x)]
      if (length(v)) v[which.max(abs(v))] else NA_real_
    }), .groups = "drop")
  left_join(d, s, by = "gene_name")
}
cand <- cand |>
  join_side(ANCHOR, "anchor") |>
  join_side(REANNOT, "reannot") |>
  join_side(CROSS, "cross") |>
  left_join(
    anchor_sw |> transmute(
      gene_name = .data$gene_name,
      expr_anchor = .data$gene_expression,
      n_iso_anchor = .data$n_switching_isoforms
    ) |>
      group_by(.data$gene_name) |>
      summarize(across(where(is.numeric), ~ suppressWarnings(max(.x, na.rm = TRUE))),
                .groups = "drop"),
    by = "gene_name"
  )

cand <- cand |>
  mutate(
    switches_anchor = .data$gene_name %in% sym_of(ANCHOR),
    testable_anchor = .data$gene_name %in% tested_anchor,
    switches_reannot = .data$gene_name %in% sw_reannot,
    testable_reannot = .data$gene_name %in% tested_reannot,
    switches_cross = .data$gene_name %in% sw_cross,
    testable_cross = .data$gene_name %in% tested_cross,
    # Annotation-robust: called in BOTH transcript spaces for the same libraries.
    annotation_robust = .data$switches_anchor & .data$switches_reannot,
    # Failing robustness has two very different causes and the project's standing rule is
    # never to conflate them: a gene TESTED in the other annotation and not called is a
    # real failure, while a gene below the expression floor there was never given the
    # chance. IRAK3 is the worked example -- gene expression 10.3 in T_HT against a floor
    # of 12, so its absence from the anchor is missing data, not disagreement.
    annotation_status = case_when(
      .data$switches_anchor & .data$switches_reannot ~ "robust (called in both)",
      .data$switches_anchor & !.data$testable_reannot ~ paste0("not testable in ", REANNOT),
      !.data$switches_anchor & !.data$testable_anchor ~ paste0("not testable in ", ANCHOR),
      TRUE ~ "tested in both, called in one"
    ),
    called_in = case_when(
      .data$switches_anchor & .data$switches_reannot ~ "both annotations",
      .data$switches_anchor ~ paste0(ANCHOR, " only"),
      TRUE ~ paste0(REANNOT, " only")
    ),
    # The effect used for ranking: the anchor where available, else the re-annotation.
    # lfc_rank follows dif_rank's source, so the visibility call below never mixes the
    # dIF from one annotation with the fold-change from the other.
    dif_rank = dplyr::coalesce(.data$dif_anchor, .data$dif_reannot),
    lfc_rank = ifelse(is.finite(.data$dif_anchor), .data$lfc_anchor, .data$lfc_reannot),
    visibility = ifelse(abs(.data$lfc_rank) < cryptic_lfc,
                        "cryptic (DE would miss it)", "also changes expression"),
    cross_genotype_support = case_when(
      !.data$testable_cross ~ "not testable in H_HT",
      .data$switches_cross ~ "also switches in H_HT",
      TRUE ~ "tested in H_HT, not switching"
    )
  ) |>
  arrange(desc(.data$dif_rank)) |>
  mutate(rank = row_number()) |>
  # Plain names here: select() is a tidyselect context, where .data$ is deprecated.
  # `.data$` stays in the data-masking verbs above, per the project style.
  select("rank", "gene_name", "dif_rank", "lfc_rank", "visibility",
         "called_in", "annotation_robust", "annotation_status", "cross_genotype_support",
         "dif_anchor", "dif_reannot", "dif_cross",
         "lfc_anchor", "lfc_reannot",
         "expr_anchor", "n_iso_anchor",
         "switches_anchor", "testable_anchor", "switches_reannot", "testable_reannot",
         "switches_cross", "testable_cross")

write_table_pair(cand, results_tables, "q1_candidates", cfg = cfg)

# ---- Two lists, two questions ------------------------------------------------------------
# These are deliberately NOT one ranked table with a flag. They answer different questions
# and a gene can be a strong answer to one and irrelevant to the other:
#
#   LIST A -- "what does LPS do to isoform usage in wildtype?"  The main Q1 question.
#   Ranked purely by effect size. Whether the gene also changes in abundance is beside the
#   point: NCOA7 shifts dIF 0.65 and that is the biggest switch in the data whether or not
#   DE would also have found it.
#
#   LIST B -- "which LPS switches would differential expression miss?"  A sub-question of
#   Q1, and the one the hidden-layer analysis is about. Membership is the finding here;
#   within the set, effect size is secondary. Ranking A by effect and then filtering to
#   cryptic would bury the point of B, which is why they are separate outputs.
#
# Both use the same cryptic definition as script 09 (config analysis.cryptic_max_abs_lfc),
# so there is one threshold in the project, not two.

list_a <- cand |>
  filter(.data$annotation_robust) |>
  arrange(desc(.data$dif_rank)) |>
  head(shortlist_n) |>
  mutate(list_rank = row_number())

write_table_pair(list_a, results_tables, "q1_response_candidates", cfg = cfg)

# List B keeps every cryptic candidate in the union, not just the annotation-robust ones:
# at this set size the robust subset is small, and `annotation_status` already says which
# non-robust entries are genuine disagreements and which were never testable.
list_b <- cand |>
  filter(grepl("^cryptic", .data$visibility)) |>
  arrange(desc(.data$annotation_robust), desc(.data$dif_rank)) |>
  mutate(list_rank = row_number())

write_table_pair(list_b, results_tables, "q1_cryptic_candidates", cfg = cfg)

# ---- Summary --------------------------------------------------------------------------
cryptic_all <- grepl("^cryptic", cand$visibility)
summary_tbl <- tibble(
  metric = c(
    "genes switching in the anchor (T_HT)",
    "genes switching in the re-annotation (T_UT)",
    "union across the two wildtype annotations",
    "annotation-robust (switching in both)",
    "  never testable in one annotation",
    "  tested in both, called in one",
    "LIST A -- response candidates (robust, ranked)",
    "LIST B -- cryptic candidates (invisible to DE)",
    "  of those, annotation-robust",
    "  of those, also switching in H_HT"
  ),
  n = c(
    sum(cand$switches_anchor),
    sum(cand$switches_reannot),
    nrow(cand),
    sum(cand$annotation_robust),
    sum(grepl("^not testable", cand$annotation_status)),
    sum(cand$annotation_status == "tested in both, called in one"),
    nrow(list_a),
    nrow(list_b),
    sum(cryptic_all & cand$annotation_robust),
    sum(cryptic_all & cand$switches_cross)
  )
)
write_table_pair(summary_tbl, results_tables, "q1_candidates_summary", cfg = cfg)
print(as.data.frame(summary_tbl))

message("\nLIST A -- the LPS response, ranked by effect:")
print(as.data.frame(list_a |> select("list_rank", "gene_name", "dif_rank", "lfc_rank",
                                     "cross_genotype_support") |> head(12)))
message("\nLIST B -- switches differential expression would miss:")
print(as.data.frame(list_b |> select("list_rank", "gene_name", "dif_rank", "lfc_rank",
                                     "annotation_status") |> head(12)))

# ---- Figures --------------------------------------------------------------------------
if (has_ggplot) {
  both <- cand |> filter(.data$testable_reannot, is.finite(.data$dif_anchor),
                         is.finite(.data$dif_reannot))
  if (nrow(both)) {
    p1 <- ggplot(both, aes(x = .data$dif_anchor, y = .data$dif_reannot)) +
      geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey60") +
      geom_hline(yintercept = 0.15, linetype = 3, colour = "grey40") +
      geom_vline(xintercept = 0.15, linetype = 3, colour = "grey40") +
      geom_point(aes(colour = .data$annotation_robust), alpha = 0.8, size = 2) +
      scale_colour_manual(
        values = c("TRUE" = "#1f78b4", "FALSE" = "grey65"),
        labels = c("TRUE" = "called in both", "FALSE" = "called in one"),
        name = NULL
      ) +
      labs(
        title = "Q1 candidates: the same six libraries under two annotations",
        subtitle = paste0("Agreement here is robustness to transcript space, NOT replication.",
                          " Dotted lines = the |dIF| ", 0.15, " threshold."),
        x = paste0("max |dIF| in ", ANCHOR, " (anchor)"),
        y = paste0("max |dIF| in ", REANNOT, " (re-annotation)")
      )
    ggsave(file.path(fig_dir, "fig_q1_annotation_agreement.png"), p1,
           width = 8, height = 5.5, dpi = 200)
  }

  lollipop <- function(d, title, sub, colour_by, palette) {
    d <- d |> mutate(gene_name = factor(.data$gene_name, levels = rev(.data$gene_name)))
    ggplot(d, aes(x = .data$dif_rank, y = .data$gene_name)) +
      geom_segment(aes(x = 0, xend = .data$dif_rank, yend = .data$gene_name),
                   colour = "grey75") +
      geom_point(aes(colour = .data[[colour_by]]), size = 3) +
      scale_colour_manual(values = palette, name = NULL) +
      labs(title = title, subtitle = sub,
           x = "max |dIF| among switching isoforms", y = NULL) +
      theme(legend.position = "top")
  }

  if (nrow(list_a)) {
    p2 <- lollipop(
      list_a,
      paste0("LIST A -- the wildtype LPS response, top ", nrow(list_a)),
      "Annotation-robust: switching in both transcript spaces for these libraries",
      "cross_genotype_support",
      c("also switches in H_HT" = "#1f78b4",
        "tested in H_HT, not switching" = "grey60",
        "not testable in H_HT" = "grey82")
    )
    ggsave(file.path(fig_dir, "fig_q1_response_candidates.png"), p2,
           width = 8.5, height = 0.32 * nrow(list_a) + 2.2, dpi = 200)
  }

  if (nrow(list_b)) {
    p3 <- lollipop(
      list_b,
      paste0("LIST B -- switches differential expression would miss (", nrow(list_b), ")"),
      paste0("Wildtype LPS switchers with |gene log2FC| < ", cryptic_lfc,
             "; membership is the finding, rank is secondary"),
      "annotation_status",
      c("robust (called in both)" = "#e6550d",
        "tested in both, called in one" = "grey60",
        "not testable in T_UT" = "grey82",
        "not testable in T_HT" = "grey82")
    )
    ggsave(file.path(fig_dir, "fig_q1_cryptic_candidates.png"), p3,
           width = 8.5, height = 0.28 * nrow(list_b) + 2.2, dpi = 200)
  }
  message("Figures: ", fig_dir)
} else {
  message("ggplot2 not available; tables written, figures skipped.")
}

write_run_manifest("11_q1_candidates.R", cfg, root)
message("11_q1_candidates.R: done")
