#!/usr/bin/env Rscript
# 19: Which findings need PacBio, and which would short-read Illumina have found?
#
# Run (after 03, 13, 15, 18): Rscript scripts/19_novelty_partition.R
# Pure recode + join over committed tables -- no network, no drive, no GFF3.
#
# WHY THIS EXISTS
#
# `is_novel_pacbio` is assigned by `tag_novel_isoforms()` from the TRANSCRIPT ID PREFIX
# alone -- `oId` starting "PB" is novel, "ENST" is known. That is a statement about which
# pipeline emitted the transcript model, NOT about whether its structure is in the
# reference. gffcompare's `class_code` is the structural statement, and the two disagree
# for thousands of isoforms: a PacBio transcript with class_code "=" has an intron chain
# that EXACTLY matches a named Ensembl transcript. Short-read Illumina quantified against
# Ensembl would have that transcript in its index and would have found that switch.
#
# So `pct_switching_that_are_novel` in novel_isoform_dataset_summary.csv (03) counts
# pipeline provenance, not discovery. This script recounts it structurally.
#
# THE FOUR TIERS (isoform level)
#
#   reference     oId ENST*          carried reference model; always class "="
#   pb_redundant  oId PB*, class "="  PacBio model, intron chain already in Ensembl
#   pb_contained  oId PB*, class "c"  contained in a reference transcript -- a FRAGMENT.
#                                     Novel only if its 5' end is real; this is exactly
#                                     the population the cap-selection question governs.
#   pb_novel      oId PB*, class j/k/m/n/o/u/x/i/y/e/p/s
#                                     structure absent from the reference
#
# `illumina_findable` = tier in (reference, pb_redundant).
#
# THE SWITCH IS THE UNIT, NOT THE ISOFORM
#
# A switch is resolvable by short reads only if BOTH of its isoforms are separable in the
# reference annotation. Two ways that fails:
#   (a) either side is a novel structure (class != "="), or
#   (b) both sides match the SAME reference transcript (cmp_ref_candidate ==
#       cmp_ref_partner) -- then Ensembl holds one transcript where PacBio resolves two,
#       and no annotation-based short-read quantifier can see the switch at all.
# (b) is the case that the isoform-level view misses entirely, and it is not rare here.
#
# CAVEAT ON "=" THAT THIS SCRIPT DOES NOT RESOLVE
#
# gffcompare "=" requires a matching INTRON CHAIN; the 5' and 3' ends may differ freely.
# A pb_redundant transcript can therefore still carry a PacBio-discovered TSS. Where both
# sides are "=" against DISTINCT Ensembl transcripts the switch is findable regardless,
# because both intron chains are annotated. Where they are "=" against the SAME transcript
# case (b) already catches it. What remains unresolved is the magnitude of a TSS shift
# within a matched intron chain, which needs the GFF3 and belongs to 15, not here.

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
fig_dir <- ensure_dir(file.path(results_figures, "novelty"))

theme_set(
  theme_bw(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold", size = 12),
      strip.text = element_text(face = "bold")
    )
)

# Class codes that denote a structure absent from the reference. "c" is deliberately NOT
# here -- it is tracked separately because a contained transcript may be a 5'-incomplete
# assembly rather than a discovery.
NOVEL_CLASS_CODES <- c("j", "k", "m", "n", "o", "u", "x", "i", "y", "e", "p", "s")

#' Assign the four-tier detectability class from provenance + structure.
#' @param oId original transcript id (ENST* = carried reference, PB* = PacBio model)
#' @param class_code gffcompare class code against the reference
assign_detect_tier <- function(oId, class_code) {
  o <- as.character(oId)
  cc <- as.character(class_code)
  cc[is.na(cc) | !nzchar(cc)] <- NA_character_
  is_pb <- !is.na(o) & startsWith(o, "PB")
  out <- rep(NA_character_, length(o))
  out[!is_pb & !is.na(o)] <- "reference"
  out[is_pb & cc %in% "="] <- "pb_redundant"
  out[is_pb & cc %in% "c"] <- "pb_contained"
  out[is_pb & cc %in% NOVEL_CLASS_CODES] <- "pb_novel"
  out
}

TIER_LEVELS <- c("reference", "pb_redundant", "pb_contained", "pb_novel")
TIER_LABELS <- c(
  reference    = "Reference model (ENST)",
  pb_redundant = "PacBio, chain already in Ensembl",
  pb_contained = "PacBio, contained fragment (c)",
  pb_novel     = "PacBio, novel structure"
)
TIER_COLS <- c(
  reference    = "#4F7CAC",
  pb_redundant = "#9DBEDC",
  pb_contained = "#F2B134",
  pb_novel     = "#C1443C"
)

# ---------------------------------------------------------------------------------------
# Part A -- recount switching-isoform novelty structurally, per dataset
# ---------------------------------------------------------------------------------------

message("\n==== Part A: dataset-level recount ====")

ds_keys <- names(cfg$datasets %||% list())
lab_of <- function(k) sanitize(as.character((cfg$datasets[[k]] %||% list())$label %||% k)[1L])

iso_tiers <- list()
for (k in ds_keys) {
  lab <- lab_of(k)
  src <- load_scoring_table(processed_dir, lab, quiet = TRUE)
  if (is.null(src)) {
    message("  [", lab, "] no scoring table; skipped.")
    next
  }
  iso <- score_isoforms(src, cfg, dataset_key = k, dataset_label = lab)
  iso$detect_tier <- assign_detect_tier(iso$oId, iso$class_code)
  iso_tiers[[lab]] <- iso |>
    transmute(
      dataset_label = lab,
      isoform_id = .data$isoform_id,
      gene_name = .data$gene_name,
      oId = .data$oId,
      class_code = .data$class_code,
      cmp_ref = .data$cmp_ref,
      is_novel_id_prefix = .data$is_novel,
      detect_tier = .data$detect_tier,
      is_switching = .data$is_switching,
      abs_dIF = .data$abs_dIF
    )
  message(sprintf("  [%-20s] %6d isoforms, %4d switching", lab, nrow(iso), sum(iso$is_switching)))
}
if (!length(iso_tiers)) stop("No dataset produced a scoring table. Run 01_load_data.R first.")
iso_all <- bind_rows(iso_tiers)

unassigned <- sum(is.na(iso_all$detect_tier))
if (unassigned > 0L) {
  message("  NOTE: ", unassigned, " isoforms could not be tiered (missing oId or class_code).")
}

# Per-dataset tier counts among SWITCHING isoforms, plus the three nested novelty rates.
tier_counts <- iso_all |>
  filter(.data$is_switching) |>
  count(.data$dataset_label, .data$detect_tier, name = "n") |>
  mutate(detect_tier = factor(.data$detect_tier, levels = TIER_LEVELS))

detect_summary <- iso_all |>
  filter(.data$is_switching) |>
  group_by(.data$dataset_label) |>
  summarize(
    n_switching = dplyr::n(),
    n_reference = sum(.data$detect_tier %in% "reference"),
    n_pb_redundant = sum(.data$detect_tier %in% "pb_redundant"),
    n_pb_contained = sum(.data$detect_tier %in% "pb_contained"),
    n_pb_novel = sum(.data$detect_tier %in% "pb_novel"),
    .groups = "drop"
  ) |>
  mutate(
    # as published by 03: provenance only
    pct_novel_id_prefix = 100 * (.data$n_pb_redundant + .data$n_pb_contained + .data$n_pb_novel) /
      .data$n_switching,
    # structural, counting contained fragments as novel
    pct_novel_structural = 100 * (.data$n_pb_contained + .data$n_pb_novel) / .data$n_switching,
    # structural and strict: contained fragments excluded
    pct_novel_strict = 100 * .data$n_pb_novel / .data$n_switching,
    pct_illumina_findable = 100 * (.data$n_reference + .data$n_pb_redundant) / .data$n_switching,
    overcount_pp = .data$pct_novel_id_prefix - .data$pct_novel_strict
  ) |>
  arrange(.data$dataset_label)

write_table_pair(detect_summary, results_tables, "novelty_detectability_summary", cfg = cfg)
write_table_pair(tier_counts, results_tables, "novelty_detectability_tier_counts", cfg = cfg)
write_table_pair(
  iso_all |> filter(.data$is_switching),
  results_tables, "novelty_switching_isoform_tiers", cfg = cfg,
  csv = isTRUE((cfg$output %||% list())$csv_twin_isoform_level %||% FALSE)
)

message("\n---- switching isoforms: published vs structural novelty ----")
print(as.data.frame(
  detect_summary |>
    transmute(
      dataset = .data$dataset_label,
      n_sw = .data$n_switching,
      published_pct_novel = round(.data$pct_novel_id_prefix, 1),
      structural_pct_novel = round(.data$pct_novel_structural, 1),
      strict_pct_novel = round(.data$pct_novel_strict, 1),
      overcount_pp = round(.data$overcount_pp, 1)
    )
), row.names = FALSE)

# ---------------------------------------------------------------------------------------
# Part B -- switch-level detectability for the candidate shortlist
# ---------------------------------------------------------------------------------------

message("\n==== Part B: candidate switch detectability ====")

pairs <- readRDS(file.path(results_tables, "candidate_switch_pairs.rds"))
prof_path <- file.path(results_tables, "candidate_profiles.rds")
prof <- if (file.exists(prof_path)) readRDS(prof_path) else NULL

# The candidate PAIRS were resolved in the UT annotation space (13/15 read U_T.gff3 and
# fa/transcripts_U-T.fa), so oId must be looked up there, not in the HT space. T_UT is the
# wildtype +/- LPS arm, which is the same biology as the T_HT anchor the symbols came from.
PAIR_SPACE <- "T_UT"
ctx_ids <- iso_all |>
  filter(.data$dataset_label == PAIR_SPACE) |>
  select("isoform_id", "oId", "class_code", "cmp_ref", "detect_tier")
if (!nrow(ctx_ids)) {
  stop("No scoring table for ", PAIR_SPACE, "; cannot resolve candidate oIds.")
}
need_ids <- unique(c(pairs$isoform_id, pairs$partner_id))
miss <- setdiff(need_ids, ctx_ids$isoform_id)
if (length(miss)) {
  stop(
    "Candidate pair ids absent from ", PAIR_SPACE, ": ", length(miss), " of ",
    length(need_ids), ". The pair table and the context layer are in different ",
    "annotation spaces -- do not join them."
  )
}
message("  resolved all ", length(need_ids), " candidate isoform ids in ", PAIR_SPACE)

cand <- pairs |>
  left_join(
    ctx_ids |> select(isoform_id = "isoform_id", oId_candidate = "oId",
                      tier_candidate = "detect_tier"),
    by = "isoform_id"
  ) |>
  left_join(
    ctx_ids |> select(partner_id = "isoform_id", oId_partner = "oId",
                      tier_partner = "detect_tier"),
    by = "partner_id"
  ) |>
  mutate(
    cc_cand = as.character(.data$class_code_candidate),
    cc_part = as.character(.data$class_code_partner),
    ref_cand = as.character(.data$cmp_ref_candidate),
    ref_part = as.character(.data$cmp_ref_partner),
    # (a) either side is a structure the reference does not contain
    side_novel_structure = !(.data$cc_cand %in% "=") | !(.data$cc_part %in% "="),
    # (b) both sides collapse onto one reference transcript
    same_reference_transcript = !is.na(.data$ref_cand) & !is.na(.data$ref_part) &
      .data$ref_cand == .data$ref_part,
    illumina_resolvable = !.data$side_novel_structure & !.data$same_reference_transcript,
    detect_reason = dplyr::case_when(
      .data$same_reference_transcript & .data$side_novel_structure ~
        "PacBio only: novel structure, and both sides map to one Ensembl transcript",
      .data$same_reference_transcript ~
        "PacBio only: both sides map to the same Ensembl transcript",
      .data$side_novel_structure ~
        "PacBio only: one or both sides are a novel structure",
      TRUE ~
        "Illumina-findable: both sides are distinct annotated Ensembl transcripts"
    ),
    rising_is_novel_structure = !(.data$cc_cand %in% "="),
    detect_class = ifelse(.data$illumina_resolvable, "illumina_findable", "pacbio_required")
  )

# Carry the decision columns from 18 so the partition can be read against priority, and
# keep the locus-fusion disqualification visible -- a fused locus is not a switch at all,
# so its detectability call is meaningless and must not be read as a candidate.
if (!is.null(prof)) {
  cand <- cand |>
    left_join(
      # `coding_change` (and 17 other columns) exist in BOTH tables -- 18 carried them up
      # from 15. Selecting it here too would silently produce coding_change.x/.y and the
      # any_of() below would then match NEITHER, dropping the column without an error.
      # Take the structural columns from `pairs` and only the DECISION columns from `prof`.
      prof |> select(
        "gene_name", "priority_score", "assayable", "protein_changing",
        "changes_annotated_site", "structurally_ambiguous", "locus_fusion",
        "longest_unique_bp", "n_isoform_papers"
      ),
      by = "gene_name"
    )
  stopifnot(!any(grepl("\\.(x|y)$", names(cand))))
} else {
  message("  (candidate_profiles.rds absent; priority columns omitted)")
}

cand_out <- cand |>
  select(
    "gene_name",
    isoform_id = "isoform_id", partner_id = "partner_id",
    "oId_candidate", "oId_partner",
    class_code_candidate = "cc_cand", class_code_partner = "cc_part",
    cmp_ref_candidate = "ref_cand", cmp_ref_partner = "ref_part",
    "tier_candidate", "tier_partner",
    "side_novel_structure", "same_reference_transcript",
    "rising_is_novel_structure",
    "illumina_resolvable", "detect_class", "detect_reason",
    "alt_first_exon", "as_events", "contained_caveat",
    "coding_change", "frame_relationship",
    dplyr::any_of(c(
      "priority_score", "assayable", "protein_changing", "changes_annotated_site",
      "structurally_ambiguous", "locus_fusion", "longest_unique_bp",
      "n_isoform_papers", "coding_change"
    ))
  ) |>
  arrange(dplyr::desc(.data$detect_class == "pacbio_required"),
          dplyr::desc(dplyr::coalesce(.data$priority_score, -1L)),
          .data$gene_name)

write_table_pair(cand_out, results_tables, "candidate_detectability", cfg = cfg)

# ONE ROW PER SWITCH. `candidate_switch_pairs` is keyed on the CANDIDATE isoform, so a gene
# whose two isoforms are both significant ("both sides") contributes two rows describing the
# SAME switch with the roles swapped (A->B and B->A). Counting those twice would inflate
# every total. Collapse on the unordered {isoform, partner} pair; the detectability call is
# symmetric, so either row's verdict stands for the switch.
cand_switch <- cand_out |>
  mutate(
    pair_key = paste(
      pmin(.data$isoform_id, .data$partner_id),
      pmax(.data$isoform_id, .data$partner_id),
      sep = "|"
    ),
    both_sides_significant = FALSE
  ) |>
  group_by(.data$pair_key) |>
  mutate(both_sides_significant = dplyr::n() > 1L) |>
  # keep the row whose candidate side is the novel structure where there is one, so the
  # retained row reads in the informative direction
  arrange(dplyr::desc(.data$rising_is_novel_structure), .by_group = TRUE) |>
  dplyr::slice(1L) |>
  ungroup() |>
  arrange(dplyr::desc(.data$detect_class == "pacbio_required"),
          dplyr::desc(dplyr::coalesce(.data$priority_score, -1L)),
          .data$gene_name)

write_table_pair(cand_switch, results_tables, "candidate_detectability_by_switch", cfg = cfg)

# Cross-tab: detectability against the locus-fusion disqualification, per switch.
cand_xtab <- cand_switch |>
  mutate(fused = dplyr::coalesce(.data$locus_fusion, FALSE)) |>
  count(.data$detect_class, .data$detect_reason, .data$fused, name = "n") |>
  arrange(.data$detect_class, dplyr::desc(.data$n))
write_table_pair(cand_xtab, results_tables, "candidate_detectability_xtab", cfg = cfg)

# ---------------------------------------------------------------------------------------
# Figures
# ---------------------------------------------------------------------------------------

p_tiers <- ggplot(
  tier_counts,
  aes(x = .data$dataset_label, y = .data$n, fill = .data$detect_tier)
) +
  geom_col(width = 0.68) +
  scale_fill_manual(
    values = TIER_COLS, labels = TIER_LABELS, name = NULL, drop = FALSE
  ) +
  labs(
    title = "Switching isoforms: pipeline provenance vs structural novelty",
    subtitle = paste(
      "Blue = short-read Illumina against Ensembl would have found it.",
      "Red = structure absent from the reference.",
      "Amber = contained fragment, novel only if its 5' end is real.",
      sep = "\n"
    ),
    x = NULL, y = "switching isoforms"
  ) +
  theme(legend.position = "right")
ggsave(file.path(fig_dir, "switching_isoform_detectability.png"),
       plot = p_tiers, width = 9.5, height = 4.6, dpi = 200)

cand_plot <- cand_switch |>
  filter(!dplyr::coalesce(.data$locus_fusion, FALSE)) |>
  filter(!is.na(.data$priority_score)) |>
  mutate(
    gene_name = stats::reorder(.data$gene_name, .data$priority_score),
    detect_class = factor(
      .data$detect_class,
      levels = c("pacbio_required", "illumina_findable"),
      labels = c("PacBio required", "Illumina-findable")
    )
  )

if (nrow(cand_plot)) {
  # A lollipop, not geom_col: one gene (CD86) carries TWO distinct switch pairs, and
  # geom_col would STACK them into a single bar of twice the priority score. Dodging also
  # keeps genes with priority_score 0 (SOCS4, R3HCC1L) visible, which a zero-length bar
  # would not.
  dodge <- position_dodge2(width = 0.75, padding = 0.1)
  p_cand <- ggplot(
    cand_plot,
    aes(x = .data$priority_score, y = .data$gene_name, colour = .data$detect_class)
  ) +
    geom_segment(
      aes(x = 0, xend = .data$priority_score,
          y = .data$gene_name, yend = .data$gene_name),
      position = dodge, linewidth = 1.6
    ) +
    geom_point(position = dodge, size = 2.6) +
    scale_colour_manual(
      values = c("PacBio required" = "#C1443C", "Illumina-findable" = "#4F7CAC"),
      name = NULL
    ) +
    labs(
      title = "Candidate switches by priority and by what platform can see them",
      subtitle = paste0(
        "Locus fusions excluded (not switches). A switch needs PacBio when either side is a\n",
        "novel structure, or when both sides map to the same Ensembl transcript."
      ),
      x = "priority_score (sort key from 18, not a ranking)", y = NULL
    ) +
    theme(legend.position = "top")
  ggsave(file.path(fig_dir, "candidate_detectability.png"),
         plot = p_cand, width = 8.2, height = 7.4, dpi = 200)
}

# ---------------------------------------------------------------------------------------

message("\n---- candidate switches ----")
message("  candidate-isoform rows in 15:            ", nrow(cand_out))
message("  distinct switches after collapsing:      ", nrow(cand_switch),
        "  (", sum(cand_switch$both_sides_significant), " with both sides significant)")
message("  locus fusions (not switches, excluded):  ",
        sum(cand_switch$locus_fusion %in% TRUE))
real <- cand_switch |> filter(!dplyr::coalesce(.data$locus_fusion, FALSE))
message("  real switches:                           ", nrow(real),
        " across ", dplyr::n_distinct(real$gene_name), " genes")
message("    Illumina-findable:                     ", sum(real$detect_class == "illumina_findable"))
message("    PacBio required:                       ", sum(real$detect_class == "pacbio_required"))
message("      ... a side is a novel structure:     ", sum(real$side_novel_structure))
message("      ... both map to one ENST:            ", sum(real$same_reference_transcript))
message("  rising (LPS-induced) side is novel:      ", sum(real$rising_is_novel_structure))
message("  PacBio-required AND assayable:           ",
        sum(real$detect_class == "pacbio_required" & real$assayable %in% TRUE))

message("\n---- PacBio-dependent switches, by priority ----")
print(as.data.frame(
  real |>
    filter(.data$detect_class == "pacbio_required") |>
    transmute(
      gene = .data$gene_name,
      prio = .data$priority_score,
      cc = paste0(.data$class_code_candidate, "/", .data$class_code_partner),
      same_ENST = .data$same_reference_transcript,
      assayable = .data$assayable,
      unique_bp = .data$longest_unique_bp
    )
), row.names = FALSE)

write_run_manifest(script = "19_novelty_partition.R", cfg = cfg, root = root)
message("\n19_novelty_partition.R: done")
