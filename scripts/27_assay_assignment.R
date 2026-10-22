#!/usr/bin/env Rscript
# 27: Which assay can actually resolve each switch, and what it would read out.
#
# Run (after 26): Rscript scripts/27_assay_assignment.R
# Pure join + routing over committed tables. No GFF3, no FASTA, no drive.
#
# WHY AN ASSAY MATRIX RATHER THAN A PRIMER LIST
#
# 26 found that 29 of the 44 surviving switches are alternative-first-exon events and only 13
# are qPCR-addressable: 10 candidates are nested inside a sibling isoform and others have no
# candidate-only stretch long enough to prime into. An isoform that differs only by where it
# STARTS cannot be distinguished by a unique-region primer, because it has no unique region.
# Designing qPCR for everything would have covered 30% of the list and silently dropped the
# dominant biology -- including KDM6A and NCOA7, the two largest effects.
#
# The assay has to match the structural feature. The routing below is therefore driven by
# `primary_event` and `qpcr_addressable` from 26, plus PTC and ORF coordinates from the
# context layer.
#
# WHAT EACH ASSAY ACTUALLY MEASURES -- and what it cannot
#
#   qPCR / ddPCR over a unique region   Relative (qPCR) or absolute (ddPCR) abundance of one
#       isoform. Cheap, quantitative, and the only one of these that scales to many samples.
#       Reads abundance, NOT structure: it confirms the isoform is there and moves, and says
#       nothing about where it starts. ddPCR is preferable for a dIF of 0.15-0.25, where a
#       2^-ddCt ratio on a shifting denominator is fragile.
#   Junction-spanning qPCR              Same, but the discriminator is a junction rather than
#       a stretch of sequence. Appropriate for exon skipping and alternative splice sites.
#       An intron-spanning pair is the standard readout for intron retention.
#   5' RACE                             Directly sequences the 5' end. THE assay for an
#       alternative first exon: it returns the actual TSS rather than inferring it. Low
#       throughput and semi-quantitative; template-switching RACE can produce artefactual
#       ends, so a known-TSS positive control matters.
#   CAGE / nanoCAGE                     Genome-wide, cap-dependent TSS quantification. More
#       expensive and not per-candidate, but it answers the whole shortlist at once AND the
#       standing cap-selection question -- whether these 5' ends are real caps or assembly
#       truncation. For a list that is 2/3 promoter switches this is arguably the single
#       highest-value experiment available.
#   Targeted long-read cDNA             Amplify the full transcript with flanking primers and
#       sequence it (PacBio or ONT). The universal fallback: it reads the ENTIRE structure in
#       one molecule, so it resolves nested isoforms, multi-difference pairs and phasing that
#       no short amplicon can. Quantitative only crudely, and PCR over long templates biases
#       toward short products -- which matters here, because the shorter isoform is often the
#       one being measured.
#   3' RACE / 3'-end seq (QuantSeq)     The 3' equivalent, for alternative last exons.
#   NMD inhibition (UPF1 KD or CHX)     For PTC-positive isoforms: if the isoform is an NMD
#       target, blocking NMD should raise it specifically. Tests the predicted FATE rather
#       than the structure. CHX is blunt and has pleiotropic effects; UPF1 knockdown is
#       cleaner but slower.
#   Isoform-resolving western / MS      For switches that change the protein. Requires an
#       antibody whose epitope lies in the differing region -- which must be checked per
#       candidate. I CANNOT validate antibodies or confirm availability; these are specified
#       with the epitope requirement stated, for a human to check.
#   Ribosome / polysome profiling       For UTR-only switches, where no protein changes and
#       the plausible consequence is translational efficiency, not composition.
#   Minigene splicing reporter          Asks whether the event is cis-encoded and
#       LPS-responsive, rather than merely correlated with LPS.
#   Alternative-promoter luciferase     For promoter switches: is the alternative promoter
#       itself LPS-inducible?
#   CRISPRi / promoter deletion         The functional test: silence the alternative promoter
#       and ask whether the LPS response at that locus changes.

bt <- c("utils/bootstrap.R", file.path("..", "utils", "bootstrap.R"))
b_file <- if (any(f <- vapply(bt, file.exists, NA))) { bt[which(f)[1L]] } else { NA_character_ }
if (is.na(b_file)) stop("Could not find utils/bootstrap.R. Run from project root.")
source(b_file, local = FALSE, chdir = FALSE)
suppressPackageStartupMessages({ library(dplyr); library(tibble) })

root <- if (exists("PROJECT_ROOT", inherits = FALSE)) get("PROJECT_ROOT", inherits = FALSE) else find_project_root()
args_r <- commandArgs(trailingOnly = TRUE)
cfg <- load_yaml_config(if (length(args_r)) args_r[[1L]] else "config/config.yml", root = root)
paths <- cfg$paths %||% list()
results_tables <- ensure_dir(resolve_path(paths$results_tables %||% "results/tables", root = root))
processed_dir <- resolve_path(paths$processed_dir %||% "data/processed", root = root)

anat <- readRDS(file.path(results_tables, "ht_switch_anatomy.rds"))

# PTC and ORF coordinates for candidate and partner, from the unfiltered context layer.
ctx <- bind_rows(lapply(c("T_HT_anchor", "H_HT_supplementary"), function(lab) {
  f <- file.path(processed_dir, paste0("isoformContext_", lab, ".rds"))
  if (!file.exists(f)) return(NULL)
  x <- readRDS(f)
  tibble(isoform_id = x$isoform_id, PTC = x$PTC,
         orf_start = suppressWarnings(as.numeric(x$orfTransciptStart)),
         orf_len = suppressWarnings(as.numeric(x$orfTransciptLength)))
})) |> distinct(.data$isoform_id, .keep_all = TRUE)

a <- anat |>
  left_join(ctx |> rename(PTC_candidate = "PTC", orf_start_candidate = "orf_start",
                          orf_len_candidate = "orf_len"), by = "isoform_id") |>
  left_join(ctx |> select(partner = "isoform_id", PTC_partner = "PTC",
                          orf_start_partner = "orf_start", orf_len_partner = "orf_len"),
            by = "partner") |>
  mutate(
    # Compare ORF LENGTH only. orfTransciptStart is TRANSCRIPT-RELATIVE, so two isoforms with
    # different 5' ends have different ORF starts even when the CDS is identical -- requiring
    # equal starts as well made `orf_identical` zero for all 44, which is impossible when 29
    # of them are alternative-first-exon events. Length is origin-independent.
    # Caveat kept explicit: equal length implies equal protein LENGTH, not equal sequence;
    # confirming identity needs the translated ORF, which the context layer does not carry.
    orf_len_delta = .data$orf_len_candidate - .data$orf_len_partner,
    same_orf_length = !is.na(.data$orf_len_delta) & .data$orf_len_delta == 0,
    protein_changing = !is.na(.data$orf_len_delta) & .data$orf_len_delta != 0
  )

# ---- routing ------------------------------------------------------------------------------
# primary_assay is the ONE assay that resolves the structural difference. Everything else is
# orthogonal: it either confirms by a different principle or tests a consequence.
a <- a |>
  mutate(
    primary_assay = dplyr::case_when(
      .data$primary_event == "intron retention" ~ "intron-spanning qPCR / ddPCR",
      .data$primary_event %in% c("exon skipping", "alternative splice site") ~
        "junction-spanning qPCR / ddPCR",
      .data$primary_event == "alternative last exon (3' end)" & .data$qpcr_addressable ~
        "3'-region qPCR / ddPCR",
      .data$primary_event == "alternative last exon (3' end)" ~ "3' RACE",
      .data$primary_event %in% c("alternative first exon (TSS)", "both ends differ") &
        .data$qpcr_addressable ~ "unique-region qPCR / ddPCR (+ 5' RACE to fix the TSS)",
      .data$primary_event %in% c("alternative first exon (TSS)", "both ends differ") ~
        "5' RACE",
      TRUE ~ "targeted long-read cDNA"
    ),
    # a nested candidate with no unique stretch cannot be amplified apart from its host at
    # all, so the full-length read is the only way to see it as a distinct molecule
    fallback_assay = ifelse(.data$is_nested %in% TRUE | !.data$qpcr_addressable,
                            "targeted long-read cDNA (full structure in one molecule)",
                            "targeted long-read cDNA (confirmatory)"),
    orthogonal_assays = vapply(seq_len(dplyr::n()), function(i) {
      o <- character(0)
      if (grepl("first exon|both ends", a$primary_event[i])) {
        o <- c(o, "CAGE/nanoCAGE (global TSS; also settles cap-selection)",
               "alternative-promoter luciferase (is that promoter LPS-inducible?)",
               "CRISPRi of the alternative promoter (functional)")
      }
      if (a$primary_event[i] %in% c("exon skipping", "alternative splice site",
                                    "intron retention")) {
        o <- c(o, "minigene splicing reporter (cis-encoded and LPS-responsive?)")
      }
      if (a$primary_event[i] == "alternative last exon (3' end)") {
        o <- c(o, "3'-end sequencing (QuantSeq/PAS-seq) for polyA site usage")
      }
      if (isTRUE(a$PTC_candidate[i])) {
        o <- c(o, "NMD inhibition (UPF1 knockdown preferred, CHX blunt) - rising isoform is PTC+")
      }
      if (isTRUE(a$protein_changing[i])) {
        o <- c(o, "isoform-resolving western or targeted MS - EPITOPE MUST LIE IN THE DIFFERING REGION; antibody availability unverified")
      }
      if (isTRUE(a$same_orf_length[i])) {
        o <- c(o, "ribosome/polysome profiling - ORF LENGTH unchanged, so the consequence may be translational (confirm sequence identity first)")
      }
      if (!length(o)) o <- "none beyond the primary and fallback"
      paste(o, collapse = " | ")
    }, character(1))
  )

out <- a |>
  select("gene_name", "gene_id", "isoform_id", "partner", "discovered_in", "set", "tier",
         "own_dIF", "other_dIF", "primary_event", "is_nested", "nested_in",
         "longest_unique_bp", "unique_region", "qpcr_addressable",
         "tss_distance_bp", "tts_distance_bp",
         "PTC_candidate", "PTC_partner", "orf_len_candidate", "orf_len_partner",
         "orf_len_delta", "same_orf_length", "protein_changing",
         "primary_assay", "fallback_assay", "orthogonal_assays") |>
  arrange(dplyr::desc(abs(.data$own_dIF)))
write_table_pair(out, results_tables, "assay_assignment", cfg = cfg)

message("\n---- primary assay assignment for the ", nrow(out), " survivors ----")
print(as.data.frame(out |> count(.data$primary_assay, name = "candidates") |>
                      arrange(dplyr::desc(.data$candidates))), row.names = FALSE)
message("\n  PTC-positive rising isoforms (NMD route):   ", sum(out$PTC_candidate %in% TRUE))
message("  protein-changing (western/MS route):        ", sum(out$protein_changing %in% TRUE))
message("  same ORF length (translational route):     ", sum(out$same_orf_length %in% TRUE))
message("  nested -> long-read is the only full view:  ", sum(out$is_nested %in% TRUE))

message("\n---- top 12 by effect, with the assay that can read them ----")
print(as.data.frame(out |> head(12) |>
  transmute(gene = .data$gene_name, dIF = round(.data$own_dIF, 2),
            event = substr(.data$primary_event, 1, 26),
            assay = substr(.data$primary_assay, 1, 40))), row.names = FALSE)

write_run_manifest(script = "27_assay_assignment.R", cfg = cfg, root = root)
message("\n27_assay_assignment.R: done")
