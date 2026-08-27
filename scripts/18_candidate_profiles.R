#!/usr/bin/env Rscript
# 18: One row per candidate gene, joining everything 13/15/16/17 know about it.
#
# Run (after 13, 15, 16, 17): Rscript scripts/18_candidate_profiles.R
# Pure join over committed tables -- no network, no drive, no GFF3.
#
# WHY THIS EXISTS
#
# The facts needed to decide whether a candidate is worth a primer order are spread across
# five tables with three different keys: 13 is per isoform, 15 is per switch pair, 16 is per
# gene/feature/pair, 17 is per gene. Anyone answering "should we do RBKS or SAMD14 first"
# was joining those by hand. This does it once.
#
# THE PRIORITY SCORE IS A SORT KEY, NOT A RANKING
#
# `priority_score` adds up four things the project already treats as evidence: the switch is
# measurable at the bench (assayable, and better if both sides are), it changes the protein
# rather than only the UTR, it removes or adds an ANNOTATED functional site, and it is not
# resting on a contained isoform that could be a 5'-incomplete assembly. Those are not
# commensurable quantities and the weights are a judgement, not a measurement. Use it to
# order the table; read the columns to decide.

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

need <- c("candidate_switch_pairs", "q1_assayability")
opt <- c("candidate_functional_sites", "candidate_functional_site_summary",
         "candidate_uniprot_annotation", "candidate_literature_counts",
         "candidate_splice_sites")
rd <- function(stem, required = FALSE) {
  f <- file.path(results_tables, paste0(stem, ".rds"))
  if (!file.exists(f)) {
    if (required) stop("Missing ", stem, ".rds -- run scripts 13, 15, 16, 17 first.")
    message("  (optional input ", stem, " absent; those columns will be NA)")
    return(NULL)
  }
  readRDS(f)
}
pairs <- rd("candidate_switch_pairs", TRUE)
asy <- rd("q1_assayability", TRUE)
sites <- rd("candidate_functional_sites")
site_sum <- rd("candidate_functional_site_summary")
uni <- rd("candidate_uniprot_annotation")
lit <- rd("candidate_literature_counts")
ss <- rd("candidate_splice_sites")

message("18_candidate_profiles.R: assembling profiles for ",
        dplyr::n_distinct(pairs$gene_name), " candidate genes")

# ---- One row per gene, describing its LARGEST switch ----------------------------------------
# A gene can contribute several candidate isoforms. The profile describes the pair with the
# biggest |dIF|, which is the one an assay would target, and `n_candidate_isoforms` records
# that there were others.
top_pair <- pairs |>
  group_by(.data$gene_name) |>
  slice_max(abs(.data$dIF_candidate), n = 1, with_ties = FALSE) |>
  ungroup()

# Direction: name the rising and falling side, as 16 does, so the profile reads as biology.
top_pair <- top_pair |>
  mutate(
    rising_isoform = ifelse(.data$dIF_candidate > 0, .data$isoform_id, .data$partner_id),
    falling_isoform = ifelse(.data$dIF_candidate > 0, .data$partner_id, .data$isoform_id),
    rising_orf_aa = ifelse(.data$dIF_candidate > 0, .data$orf_aa_candidate, .data$orf_aa_partner),
    falling_orf_aa = ifelse(.data$dIF_candidate > 0, .data$orf_aa_partner, .data$orf_aa_candidate),
    rising_PTC = ifelse(.data$dIF_candidate > 0, .data$PTC_candidate, .data$PTC_partner),
    max_abs_dIF = abs(.data$dIF_candidate)
  )

per_gene_counts <- pairs |>
  group_by(.data$gene_name) |>
  summarize(
    n_candidate_isoforms = dplyr::n(),
    both_sides_candidates = dplyr::n() > 1L,
    .groups = "drop"
  )

asy_g <- asy |>
  group_by(.data$gene_name) |>
  summarize(
    n_assayable_isoforms = sum(.data$assayable %in% TRUE),
    assayable = any(.data$assayable %in% TRUE),
    longest_unique_bp = suppressWarnings(max(.data$longest_unique_bp, na.rm = TRUE)),
    .groups = "drop"
  ) |>
  mutate(longest_unique_bp = ifelse(is.finite(.data$longest_unique_bp),
                                    .data$longest_unique_bp, NA_integer_))

# Splice sites: the only thing worth surfacing per gene is whether ANY differential junction
# is non-canonical, because that is a reason to distrust the transcript model.
ss_g <- if (is.null(ss)) NULL else {
  ss |>
    group_by(.data$gene_name) |>
    summarize(
      n_differential_junctions = dplyr::n(),
      all_junctions_canonical = all(grepl("^GT-AG|^GC-AG", .data$intron_class)),
      junction_classes = paste(sort(unique(.data$intron_class)), collapse = "; "),
      .groups = "drop"
    )
}

lost_gained <- if (is.null(sites)) NULL else {
  sites |>
    filter(.data$status %in% c("LOST on switching", "partially lost on switching")) |>
    group_by(.data$gene_name) |>
    summarize(
      sites_lost = paste(unique(paste0(
        .data$feature_type,
        ifelse(nzchar(.data$feature_description), paste0(": ", .data$feature_description), "")
      )), collapse = " | "),
      .groups = "drop"
    )
}
gained_tab <- if (is.null(sites)) NULL else {
  sites |>
    filter(.data$status %in% c("GAINED on switching", "partially gained on switching")) |>
    group_by(.data$gene_name) |>
    summarize(
      sites_gained = paste(unique(paste0(
        .data$feature_type,
        ifelse(nzchar(.data$feature_description), paste0(": ", .data$feature_description), "")
      )), collapse = " | "),
      .groups = "drop"
    )
}

prof <- top_pair |>
  select("gene_name", "rising_isoform", "falling_isoform", "max_abs_dIF",
         "dIF_candidate", "IF1_candidate", "IF2_candidate",
         "rising_orf_aa", "falling_orf_aa", "rising_PTC",
         "as_events", "tss_distance_bp", "alt_first_exon", "exon_skipping",
         "intron_retention", "alt_donor", "alt_acceptor",
         "frame_relationship", "coding_change", "coding_change_detail",
         "class_code_candidate", "class_code_partner", "contained_caveat",
         "cmp_ref_candidate", "cmp_ref_partner") |>
  left_join(per_gene_counts, by = "gene_name") |>
  left_join(asy_g, by = "gene_name")
if (!is.null(ss_g)) prof <- prof |> left_join(ss_g, by = "gene_name")
if (!is.null(uni)) {
  prof <- prof |> left_join(
    uni |> select("gene_name", "uniprot_accession", "protein_name", "canonical_length_aa",
                  "uniprot_function", "subcellular_location"),
    by = "gene_name"
  )
}
if (!is.null(site_sum)) {
  prof <- prof |> left_join(
    site_sum |> select("gene_name", "n_features_changed", "lost", "gained"), by = "gene_name"
  )
}
if (!is.null(lost_gained)) prof <- prof |> left_join(lost_gained, by = "gene_name")
if (!is.null(gained_tab)) prof <- prof |> left_join(gained_tab, by = "gene_name")
if (!is.null(lit)) {
  prof <- prof |> left_join(
    lit |> select("gene_name", "n_papers_total", "n_isoform_papers", "n_immunity_papers",
                  "n_isoform_and_immunity"),
    by = "gene_name"
  )
}

# ---- Priority score ---------------------------------------------------------------------------
prof <- prof |>
  mutate(
    protein_changing = !.data$coding_change %in% c("UTR only", "unknown (ORF not called)"),
    changes_annotated_site = !is.na(.data$n_features_changed) & .data$n_features_changed > 0,
    structurally_ambiguous = grepl("^c-class AND differs only", .data$contained_caveat),
    priority_score =
      2L * as.integer(.data$assayable %in% TRUE) +
      1L * as.integer(.data$both_sides_candidates %in% TRUE) +
      2L * as.integer(.data$protein_changing %in% TRUE) +
      2L * as.integer(.data$changes_annotated_site %in% TRUE) -
      2L * as.integer(.data$structurally_ambiguous %in% TRUE)
  ) |>
  arrange(desc(.data$priority_score), desc(.data$max_abs_dIF))

write_table_pair(prof, results_tables, "candidate_profiles", cfg = cfg)

message("\n---- candidate profiles, by priority score ----")
print(as.data.frame(
  prof |>
    transmute(
      gene = .data$gene_name,
      score = .data$priority_score,
      dIF = round(.data$max_abs_dIF, 3),
      assay = .data$assayable,
      sides = .data$n_candidate_isoforms,
      consequence = .data$coding_change,
      sites = ifelse(is.na(.data$n_features_changed), 0L, .data$n_features_changed),
      iso_lit = .data$n_isoform_papers,
      ambiguous = .data$structurally_ambiguous
    )
), row.names = FALSE)

message("\n---- headline counts ----")
message("  genes profiled:                          ", nrow(prof))
message("  assayable at the bench:                  ", sum(prof$assayable %in% TRUE))
message("  switch changes the protein:              ", sum(prof$protein_changing %in% TRUE))
message("  switch changes an annotated site:        ", sum(prof$changes_annotated_site %in% TRUE))
message("  UTR-only (protein unchanged):            ",
        sum(prof$coding_change %in% "UTR only"))
message("  structurally ambiguous (c-class, 5'):    ",
        sum(prof$structurally_ambiguous %in% TRUE))
message("  no isoform-level prior art at all:       ",
        sum(prof$n_isoform_papers %in% 0L, na.rm = TRUE))

write_run_manifest(script = "18_candidate_profiles.R", cfg = cfg, root = root)
message("\n18_candidate_profiles.R: done")
