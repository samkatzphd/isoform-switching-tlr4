#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------------------
# STANDALONE explorer for the novelty / platform-detectability partition.
#
#   Rscript explore/19_novelty_partition_explore.R
#   # or paste into an R session with the project root as the working directory
#
# Needs only dplyr + ggplot2 and the CSVs under results/tables/. It does NOT source
# utils/bootstrap.R, does not read config, does not touch /Volumes/Expansion, and writes
# nothing. Every block is independently runnable -- change a filter and re-run it.
#
# WHAT THIS ANSWERS
#   Q: of the LPS isoform switches this project found, which ones required PacBio, and
#      which would short-read Illumina against Ensembl have found on its own?
#
# THE ONE THING TO UNDERSTAND BEFORE READING ANY NUMBER
#   The pipeline's `is_novel_pacbio` flag is set from the transcript ID prefix (oId starting
#   "PB"), which says which assembler emitted the model -- not whether its structure is new.
#   gffcompare's class_code is the structural statement. "=" means the intron chain exactly
#   matches a named Ensembl transcript, so Illumina would have had it in its index.
#   The two definitions disagree by 34-43 percentage points. See below.
# ---------------------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})

TBL <- "results/tables"
if (!dir.exists(TBL)) {
  stop("Run from the project root (expected ", TBL, "/). Current wd: ", getwd())
}
rd <- function(stem) {
  f <- file.path(TBL, paste0(stem, ".csv"))
  if (!file.exists(f)) stop("Missing ", f, " -- run: Rscript scripts/19_novelty_partition.R")
  utils::read.csv(f, stringsAsFactors = FALSE)
}

summary_tbl <- rd("novelty_detectability_summary")
tiers       <- rd("novelty_detectability_tier_counts")
switches    <- rd("candidate_detectability_by_switch")
by_isoform  <- rd("candidate_detectability")     # one row per candidate isoform (pre-collapse)

# =========================================================================================
# CONCLUSION 1 -- the published novelty rate counts provenance, not discovery
# =========================================================================================
# `pct_switching_that_are_novel` in novel_isoform_dataset_summary.csv is reproduced here as
# `pct_novel_id_prefix`. The structural rates are strictly smaller:
#
#   pct_novel_id_prefix    any PacBio-emitted model                  (as published by 03)
#   pct_novel_structural   PacBio AND structure not in the reference (c-class counted novel)
#   pct_novel_strict       PacBio AND a genuinely novel junction set (c-class excluded)
#
# c-class = "contained in a reference transcript" = a FRAGMENT. It is novel only if its 5'
# end is real, which is precisely the open cap-selection question. Hence two rates.

cat("\n== Conclusion 1: switching-isoform novelty, three definitions ==\n")
summary_tbl |>
  transmute(
    dataset = dataset_label,
    n_switching,
    published = round(pct_novel_id_prefix, 1),
    structural = round(pct_novel_structural, 1),
    strict = round(pct_novel_strict, 1),
    overcount_pp = round(overcount_pp, 1),
    illumina_findable_pct = round(pct_illumina_findable, 1)
  ) |>
  as.data.frame() |>
  print(row.names = FALSE)

# The counts behind it, if you want to re-derive any percentage by hand:
cat("\n-- tier counts among switching isoforms --\n")
summary_tbl |>
  select(dataset_label, n_switching, n_reference, n_pb_redundant,
         n_pb_contained, n_pb_novel) |>
  as.data.frame() |>
  print(row.names = FALSE)

# FIGURE 1 -- provenance vs structure, per dataset
TIER_LEVELS <- c("reference", "pb_redundant", "pb_contained", "pb_novel")
p1 <- tiers |>
  mutate(detect_tier = factor(detect_tier, levels = TIER_LEVELS)) |>
  ggplot(aes(x = dataset_label, y = n, fill = detect_tier)) +
  geom_col(width = 0.68) +
  scale_fill_manual(
    values = c(reference = "#4F7CAC", pb_redundant = "#9DBEDC",
               pb_contained = "#F2B134", pb_novel = "#C1443C"),
    labels = c("Reference model (ENST)", "PacBio, chain already in Ensembl",
               "PacBio, contained fragment (c)", "PacBio, novel structure"),
    name = NULL, drop = FALSE
  ) +
  labs(title = "Switching isoforms: provenance vs structural novelty",
       x = NULL, y = "switching isoforms") +
  theme_bw(base_size = 11)
print(p1)

# =========================================================================================
# CONCLUSION 2 -- the candidate shortlist splits roughly in half by platform
# =========================================================================================
# A switch is Illumina-findable only if BOTH isoforms are separable in Ensembl. It needs
# PacBio when either side is a novel structure (class != "="), OR when both sides map to the
# SAME reference transcript -- in that case Ensembl holds one transcript where PacBio
# resolves two, and no annotation-based short-read quantifier can see the switch at all.
#
# Locus fusions are excluded throughout: those are two neighbouring genes merged into one
# XLOC by the assembler, so they are not switches and their dIF is a ratio of two genes.

real <- switches |> filter(!locus_fusion)

cat("\n== Conclusion 2: candidate switches by platform ==\n")
cat("  distinct switches (fusions excluded):", nrow(real),
    "across", n_distinct(real$gene_name), "genes\n")
real |>
  count(detect_class, detect_reason, name = "n") |>
  arrange(detect_class, desc(n)) |>
  as.data.frame() |>
  print(row.names = FALSE)

cat("\n-- the PacBio-dependent switches that are also assayable at the bench --\n")
real |>
  filter(detect_class == "pacbio_required", assayable) |>
  arrange(desc(priority_score)) |>
  transmute(gene_name, priority_score,
            class_codes = paste0(class_code_candidate, "/", class_code_partner),
            same_reference_transcript, longest_unique_bp,
            coding_change, n_isoform_papers) |>
  as.data.frame() |>
  print(row.names = FALSE)

cat("\n-- Illumina-findable, assayable, highest priority --\n")
real |>
  filter(detect_class == "illumina_findable", assayable) |>
  arrange(desc(priority_score)) |>
  transmute(gene_name, priority_score,
            cmp_ref_candidate, cmp_ref_partner,
            longest_unique_bp, coding_change, n_isoform_papers) |>
  as.data.frame() |>
  print(row.names = FALSE)

# FIGURE 2 -- priority against platform. Dodged, not stacked: CD86 carries two distinct
# switch pairs, and geom_col would add their priority scores together.
p2 <- real |>
  filter(!is.na(priority_score)) |>
  mutate(
    gene_name = stats::reorder(gene_name, priority_score),
    detect_class = factor(detect_class,
                          levels = c("pacbio_required", "illumina_findable"),
                          labels = c("PacBio required", "Illumina-findable"))
  ) |>
  ggplot(aes(x = priority_score, y = gene_name, colour = detect_class)) +
  geom_segment(aes(x = 0, xend = priority_score, y = gene_name, yend = gene_name),
               position = position_dodge2(width = 0.75, padding = 0.1), linewidth = 1.6) +
  geom_point(position = position_dodge2(width = 0.75, padding = 0.1), size = 2.6) +
  scale_colour_manual(
    values = c("PacBio required" = "#C1443C", "Illumina-findable" = "#4F7CAC"), name = NULL
  ) +
  labs(title = "Candidate switches by priority and by platform",
       x = "priority_score (sort key, not a ranking)", y = NULL) +
  theme_bw(base_size = 11) + theme(legend.position = "top")
print(p2)

# =========================================================================================
# THINGS TO TRY -- edit and re-run any of these
# =========================================================================================

# (a) Switches where the LPS-INDUCED (rising) isoform is itself the novel structure. These
#     are the ones where PacBio did not just resolve the switch better -- it found the
#     induced transcript. Strongest claim available for "a novel isoform responds to LPS".
if (FALSE) {
  real |>
    filter(rising_is_novel_structure) |>
    transmute(gene_name, priority_score, class_code_candidate, tier_candidate,
              assayable, longest_unique_bp, coding_change) |>
    arrange(desc(priority_score)) |>
    as.data.frame() |> print(row.names = FALSE)
}

# (b) The hardest category: both sides collapse onto ONE Ensembl transcript. No short-read
#     reanalysis of a public dataset can corroborate these, so they need either PacBio or an
#     isoform-specific assay as their only evidence.
if (FALSE) {
  real |> filter(same_reference_transcript) |>
    transmute(gene_name, cmp_ref_candidate, class_code_candidate, class_code_partner,
              assayable, longest_unique_bp) |>
    as.data.frame() |> print(row.names = FALSE)
}

# (c) Promoter switches only (alt_first_exon). 54 of 63 candidate pairs differ by a distinct
#     first exon, so most of this shortlist turns on where transcription STARTS, not on
#     spliceosome catalysis. Cross that against platform.
if (FALSE) {
  real |> count(alt_first_exon, detect_class) |> as.data.frame() |> print(row.names = FALSE)
}

# (d) Pre-collapse view: one row per candidate isoform rather than per switch. Use this when
#     you want the A->B and B->A rows separately (19 switches have both sides significant).
if (FALSE) {
  by_isoform |> filter(gene_name == "IRAK3") |> t() |> print()
}

# (e) UTR-only switches: the protein is unchanged, so a western blot cannot read them out
#     and the assay has to be RNA-level.
if (FALSE) {
  real |> filter(coding_change == "UTR only") |>
    transmute(gene_name, detect_class, assayable, longest_unique_bp) |>
    as.data.frame() |> print(row.names = FALSE)
}

cat("\nDone. Tables read from ", TBL, "/ -- nothing was written.\n", sep = "")
