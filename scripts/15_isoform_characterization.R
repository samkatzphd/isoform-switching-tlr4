#!/usr/bin/env Rscript
# 15: What ARE the candidate isoforms, and how do the two sides of each switch differ?
#
# Run (after 01, 11, 13): Rscript scripts/15_isoform_characterization.R
# Needs gff3/U_T.gff3 and fa/transcripts_U-T.fa (LOCAL but gitignored, 1.1 GB) and
# BSgenome.Hsapiens.UCSC.hg38. The external drive is NOT needed. Guarded: warns and skips
# rather than failing, the way 13 guards the GFF3.
#
# WHAT THIS ADDS THAT 13 DOES NOT
#
# 13 asked where a discriminating PRIMER can sit -- a question about unique sequence. This
# asks what the isoforms ARE: how the two sides of each switch differ structurally, whether
# that difference reaches the protein, and what splice sites carry it.
#
# THE COMPARISON IS PAIRWISE, NOT ISA-STYLE
#
# IsoformSwitchAnalyzeR's analyzeAlternativeSplicing types each isoform against a
# "hypothetical pre-RNA" built from the union of the gene's exons, which answers "what does
# this isoform do relative to everything the gene can make". A bench scientist validating a
# switch needs the other question -- "how does THIS isoform differ from the one it is
# replacing" -- so every event here is called between an ordered PAIR of isoforms and reads
# directionally: `exon_skipping` means the candidate skips an exon the partner has. The two
# schemes will not agree isoform-for-isoform and are not meant to.
#
# THE PARTNER RULE
#
# The partner is the isoform of the same gene with the largest OPPOSITE-SIGN dIF, because a
# switch is only measurable as a switch when both sides are named. Genes whose candidate has
# no opposite-moving partner fall back to the baseline-dominant isoform (max IF1, the same
# rule 09 uses), and the fallback is recorded in `partner_rule` -- a fallback pair is a
# description of the isoform, not of a switch.
#
# THE READING-FRAME TEST, which is the part that turns structure into consequence
#
# Two isoforms can share coding sequence and still make unrelated proteins if they read it in
# different frames. For every genomic position coding in BOTH isoforms this computes each
# one's codon offset; equal offsets across the shared region mean the shared exons make the
# same protein and the difference is confined to what is added or removed. Unequal means a
# frameshift, and the shared exons downstream encode something different in each. That
# distinction cannot be recovered from exon overlap alone.
#
# STANDING CAVEAT CARRIED INTO EVERY ROW
#
# `class_code` c = the isoform is CONTAINED within another. Those are either genuine
# alternative-internal-promoter products or 5'-incomplete assemblies, and without cap
# selection the second class is produced in quantity (docs/NEXT_SESSION.md, 2026-08-23). A
# `c`-class candidate whose only difference from its partner is where it STARTS is exactly
# the ambiguous case; `contained_caveat` flags those rather than leaving it to the reader.

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
has_bs <- requireNamespace("BSgenome.Hsapiens.UCSC.hg38", quietly = TRUE) &&
  requireNamespace("Biostrings", quietly = TRUE)
has_bios <- requireNamespace("Biostrings", quietly = TRUE)

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
results_seq <- ensure_dir(resolve_path("results/sequences", root = root))

asy <- cfg$assayability %||% list()
gff_rel <- asy$gff3_path %||% "gff3/U_T.gff3"
gff_path <- resolve_path(gff_rel, root = root)
fa_rel <- (cfg$characterization %||% list())$fasta_path %||% "fa/transcripts_U-T.fa"
fa_path <- resolve_path(fa_rel, root = root)

# The candidate isoform set is defined by 13, so the two stay in step by construction: one
# row per isoform 13 screened, candidate or partner.
asy_f <- file.path(results_tables, "q1_assayability.rds")
if (!file.exists(asy_f)) {
  stop("Missing q1_assayability.rds. Run: Rscript scripts/13_assayability_screen.R")
}
screened <- readRDS(asy_f)

message("15_isoform_characterization.R: characterizing ", nrow(screened),
        " candidate isoforms across ", dplyr::n_distinct(screened$gene_name), " genes")

# ---- Scoring table: usage, expression and the ORF annotation ------------------------------
# T_UT, because the candidate ids are TCONS_* under the UT reference (13 screens against
# gff3/U_T.gff3) and TCONS ids are NOT comparable across references.
ctx <- load_scoring_table(processed_dir, "T_UT", quiet = TRUE)
if (is.null(ctx)) stop("No T_UT scoring table in ", processed_dir, "; run 01_load_data.R.")
scored <- score_isoforms(ctx, cfg, "T_UT", "T_UT")
if (!"orfTransciptStart" %in% names(scored)) {
  stop("The context layer carries no ORF columns. Re-run 01_load_data.R with the drive ",
       "mounted so extract_context_tables() can pull orfAnalysis across.")
}

gene_ids <- unique(scored$gene_id[scored$gene_name %in% unique(screened$gene_name)])
all_iso <- scored |> filter(.data$gene_id %in% gene_ids)
message("  ", nrow(all_iso), " isoforms belong to the ", length(gene_ids), " candidate genes")

# ---- Exon structures ----------------------------------------------------------------------
if (!file.exists(gff_path)) {
  message("GFF3 not found at ", gff_rel, "; cannot characterize structure. Nothing written.")
  message("15_isoform_characterization.R: done (skipped)")
  quit(save = "no", status = 0)
}
ex <- parse_gff3_exons(gff_path, gene_ids)
message("  ", nrow(ex), " exon records for ", dplyr::n_distinct(ex$isoform_id), " isoforms")

# ---- Per-isoform structure ----------------------------------------------------------------
iso_struct <- ex |>
  group_by(.data$gene_id, .data$isoform_id) |>
  summarize(
    chr = dplyr::first(.data$chr), strand = dplyr::first(.data$strand),
    n_exons = dplyr::n(),
    locus_start = min(.data$start), locus_end = max(.data$end),
    transcript_length = sum(.data$end - .data$start + 1L),
    .groups = "drop"
  )

struct <- all_iso |>
  select("isoform_id", "gene_id", "gene_name", "IF1", "IF2", "dIF_n", "abs_dIF",
         "is_switching", "iso_value_1", "iso_value_2", "gene_expression",
         "class_code", "cmp_ref", "oId", "is_novel_pacbio", "PTC",
         "orfTransciptStart", "orfTransciptEnd", "orfTransciptLength",
         "stopDistanceToLastJunction", "orf_origin") |>
  left_join(iso_struct, by = c("isoform_id", "gene_id")) |>
  mutate(
    orf_length_aa = ifelse(is.na(.data$orfTransciptLength), NA_integer_,
                           as.integer(.data$orfTransciptLength / 3L) - 1L),
    # ISA's PTC flag restated as the rule it encodes, so a reader can see the threshold:
    # a stop codon more than 50 nt upstream of the final exon-exon junction is read by the
    # EJC as premature. `NA` means no ORF was called, which is NOT "not an NMD target".
    nmd_call = dplyr::case_when(
      is.na(.data$PTC) ~ "unknown (no ORF called)",
      .data$PTC ~ "predicted NMD target",
      TRUE ~ "not a predicted NMD target"
    )
  )

# ---- Pairwise switch partners --------------------------------------------------------------
pair_for <- function(iso_id) {
  me <- struct |> filter(.data$isoform_id == iso_id)
  if (!nrow(me)) return(NULL)
  sibs <- struct |> filter(.data$gene_id == me$gene_id[1], .data$isoform_id != iso_id,
                           is.finite(.data$dIF_n))
  if (!nrow(sibs)) return(NULL)
  opp <- sibs |> filter(sign(.data$dIF_n) == -sign(me$dIF_n[1]), .data$abs_dIF > 0)
  if (nrow(opp)) {
    p <- opp |> slice_max(.data$abs_dIF, n = 1, with_ties = FALSE)
    rule <- "largest opposite-sign dIF"
  } else {
    p <- sibs |> filter(is.finite(.data$IF1)) |> slice_max(.data$IF1, n = 1, with_ties = FALSE)
    rule <- "fallback: baseline-dominant (no opposite-moving isoform)"
  }
  if (!nrow(p)) return(NULL)
  tibble(isoform_id = iso_id, partner_id = p$isoform_id[1], partner_rule = rule)
}
pairs <- bind_rows(lapply(screened$isoform_id, pair_for))
message("  ", nrow(pairs), " candidate/partner pairs (",
        sum(grepl("^fallback", pairs$partner_rule)), " by fallback)")

# ---- Sequences -----------------------------------------------------------------------------
tx_seq <- NULL
if (has_bios && file.exists(fa_path)) {
  message("  reading transcript FASTA ", fa_rel, " ...")
  need <- unique(c(pairs$isoform_id, pairs$partner_id))
  # readDNAStringSet loads the whole 527 MB file; index by name and drop the rest at once.
  fa <- Biostrings::readDNAStringSet(fa_path)
  names(fa) <- sub("\\s.*$", "", names(fa))
  tx_seq <- fa[intersect(need, names(fa))]
  rm(fa); invisible(gc())
  message("    ", length(tx_seq), " of ", length(need), " transcript sequences found")
} else {
  message("  transcript FASTA unavailable (", fa_rel, "); protein output skipped")
}

genome <- NULL
if (has_bs) {
  genome <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
} else {
  message("  BSgenome.Hsapiens.UCSC.hg38 unavailable; splice-site sequences skipped")
}

# ---- Build the per-pair characterization ----------------------------------------------------
ex_by_iso <- split(ex[, c("start", "end")], ex$isoform_id)
strand_by_iso <- setNames(iso_struct$strand, iso_struct$isoform_id)
chr_by_iso <- setNames(iso_struct$chr, iso_struct$isoform_id)

orf_row <- function(id) {
  r <- struct |> filter(.data$isoform_id == id)
  if (!nrow(r)) return(NULL)
  list(start = r$orfTransciptStart[1], end = r$orfTransciptEnd[1])
}

pair_rows <- list()
splice_rows <- list()
prot <- list()

for (i in seq_len(nrow(pairs))) {
  a_id <- pairs$isoform_id[i]
  b_id <- pairs$partner_id[i]
  if (is.null(ex_by_iso[[a_id]]) || is.null(ex_by_iso[[b_id]])) next
  ea <- as.matrix(ex_by_iso[[a_id]][order(ex_by_iso[[a_id]]$start), ])
  eb <- as.matrix(ex_by_iso[[b_id]][order(ex_by_iso[[b_id]]$start), ])
  strand <- strand_by_iso[[a_id]]
  chrom <- chr_by_iso[[a_id]]

  ev <- pairwise_as_events(ea, eb, strand)

  # --- protein-level consequence ---
  oa <- orf_row(a_id); ob <- orf_row(b_id)
  fr <- frame_relationship(ea, eb, strand, oa, ob)

  a_meta <- struct |> filter(.data$isoform_id == a_id)
  b_meta <- struct |> filter(.data$isoform_id == b_id)

  pair_rows[[length(pair_rows) + 1L]] <- tibble(
    gene_name = a_meta$gene_name[1],
    isoform_id = a_id, partner_id = b_id,
    partner_rule = pairs$partner_rule[i],
    dIF_candidate = a_meta$dIF_n[1], dIF_partner = b_meta$dIF_n[1],
    IF1_candidate = a_meta$IF1[1], IF2_candidate = a_meta$IF2[1],
    IF1_partner = b_meta$IF1[1], IF2_partner = b_meta$IF2[1],
    class_code_candidate = a_meta$class_code[1], class_code_partner = b_meta$class_code[1],
    cmp_ref_candidate = a_meta$cmp_ref[1], cmp_ref_partner = b_meta$cmp_ref[1],
    n_exons_candidate = a_meta$n_exons[1], n_exons_partner = b_meta$n_exons[1],
    tx_len_candidate = a_meta$transcript_length[1],
    tx_len_partner = b_meta$transcript_length[1],
    as_events = ev$summary,
    n_events = ev$n_events,
    alt_TSS = ev$alt_TSS, alt_TTS = ev$alt_TTS,
    tss_distance_bp = ev$tss_distance_bp, tts_distance_bp = ev$tts_distance_bp,
    alt_first_exon = ev$alt_first_exon, alt_last_exon = ev$alt_last_exon,
    exon_skipping = ev$exon_skipping, intron_retention = ev$intron_retention,
    alt_donor = ev$alt_donor, alt_acceptor = ev$alt_acceptor,
    mutually_exclusive = ev$mutually_exclusive,
    bp_unique_to_candidate = ev$bp_unique_a,
    bp_unique_to_partner = ev$bp_unique_b,
    shared_exonic_bp = a_meta$transcript_length[1] - ev$bp_unique_a,
    # TWO TRANSCRIPTS OF ONE GENE ALWAYS SHARE SOME EXONIC SEQUENCE. When they share NONE,
    # the more likely explanation is that the assembler merged two neighbouring genes into one
    # XLOC locus, and the "switch" is really the relative expression of two different genes.
    # Confirmed downstream in 18 by looking up whether the two cmp_ref transcripts belong to
    # the same Ensembl gene -- this flag is the cheap local screen that decides what to look up.
    locus_fusion_suspect = (a_meta$transcript_length[1] - ev$bp_unique_a) <= 0,
    orf_aa_candidate = a_meta$orf_length_aa[1], orf_aa_partner = b_meta$orf_length_aa[1],
    PTC_candidate = a_meta$PTC[1], PTC_partner = b_meta$PTC[1],
    shared_coding_bp = fr$shared_coding_bp,
    frame_relationship = fr$relationship,
    coding_change = fr$coding_change,
    coding_change_detail = fr$coding_change_detail,
    coding_bp_unique_candidate = fr$coding_bp_unique_candidate,
    coding_bp_unique_partner = fr$coding_bp_unique_partner,
    # The cap-selection exposure, made explicit rather than left to the reader. A contained
    # isoform differing ONLY in where it starts is the case a 5'-incomplete assembly would
    # also produce.
    contained_caveat = if (identical(a_meta$class_code[1], "c") &&
                           !ev$exon_skipping && !ev$intron_retention &&
                           !ev$alt_donor && !ev$alt_acceptor) {
      "c-class AND differs only at the 5' end: indistinguishable from a 5'-incomplete assembly"
    } else if (identical(a_meta$class_code[1], "c")) {
      "c-class (contained), but carries an internal structural difference too"
    } else {
      NA_character_
    }
  )

  # --- splice sites at the differential junctions ---
  if (!is.null(genome) && nrow(ev$diff_introns)) {
    for (j in seq_len(nrow(ev$diff_introns))) {
      d <- ev$diff_introns[j, ]
      ss <- splice_site_seqs(genome, chrom, d[["start"]], d[["end"]], strand)
      splice_rows[[length(splice_rows) + 1L]] <- tibble(
        gene_name = a_meta$gene_name[1],
        isoform_id = a_id, partner_id = b_id,
        present_in = if (d[["in_a"]] == 1) "candidate" else "partner",
        chr = chrom, intron_start = d[["start"]], intron_end = d[["end"]], strand = strand,
        intron_length = d[["end"]] - d[["start"]] + 1L,
        donor_dinucleotide = ss$donor_di, acceptor_dinucleotide = ss$acceptor_di,
        intron_class = ss$intron_class,
        donor_site = ss$donor, acceptor_site = ss$acceptor
      )
    }
  }

  # --- translated ORFs ---
  if (!is.null(tx_seq)) {
    for (id in c(a_id, b_id)) {
      if (!is.null(prot[[id]])) next
      o <- orf_row(id)
      if (is.null(o) || is.na(o$start) || is.na(o$end)) next
      if (!id %in% names(tx_seq)) next
      p <- translate_orf(tx_seq[[id]], o$start, o$end)
      if (!is.null(p)) prot[[id]] <- p
    }
  }
}

pair_tab <- bind_rows(pair_rows)
splice_tab <- bind_rows(splice_rows)

# ---- Write ----------------------------------------------------------------------------------
struct_out <- struct |>
  select("gene_name", "isoform_id", "gene_id", "chr", "strand", "locus_start", "locus_end",
         "n_exons", "transcript_length", "class_code", "cmp_ref", "is_novel_pacbio",
         "IF1", "IF2", "dIF_n", "is_switching", "iso_value_1", "iso_value_2",
         "orfTransciptStart", "orfTransciptEnd", "orfTransciptLength", "orf_length_aa",
         "PTC", "nmd_call", "stopDistanceToLastJunction", "orf_origin") |>
  arrange(.data$gene_name, desc(abs(.data$dIF_n)))
write_table_pair(struct_out, results_tables, "candidate_isoform_structure", cfg = cfg)

if (nrow(pair_tab)) {
  pair_tab <- pair_tab |> arrange(.data$gene_name, desc(abs(.data$dIF_candidate)))
  write_table_pair(pair_tab, results_tables, "candidate_switch_pairs", cfg = cfg)
  message("\n---- structural difference between the two sides of each switch ----")
  print(as.data.frame(
    pair_tab |>
      select("gene_name", "isoform_id", "as_events", "frame_relationship", "coding_change") |>
      head(25)
  ), row.names = FALSE)
}
if (nrow(splice_tab)) {
  write_table_pair(splice_tab, results_tables, "candidate_splice_sites", cfg = cfg)
  message("\n---- splice sites at the differential junctions ----")
  message("  ", nrow(splice_tab), " junctions; intron classes: ",
          paste(names(table(splice_tab$intron_class)), table(splice_tab$intron_class),
                sep = "=", collapse = ", "))
}
if (length(prot) && has_bios) {
  aa <- Biostrings::AAStringSet(unlist(prot))
  Biostrings::writeXStringSet(aa, file.path(results_seq, "candidate_proteins.faa"))
  message("  wrote ", length(aa), " translated ORFs to results/sequences/candidate_proteins.faa")
}

# ---- Summary of what the switches DO --------------------------------------------------------
if (nrow(pair_tab)) {
  summ <- pair_tab |>
    summarize(
      n_pairs = dplyr::n(),
      alt_TSS = sum(.data$alt_TSS), alt_TTS = sum(.data$alt_TTS),
      alt_first_exon = sum(.data$alt_first_exon), alt_last_exon = sum(.data$alt_last_exon),
      exon_skipping = sum(.data$exon_skipping), intron_retention = sum(.data$intron_retention),
      alt_donor = sum(.data$alt_donor), alt_acceptor = sum(.data$alt_acceptor),
      mutually_exclusive = sum(.data$mutually_exclusive)
    ) |>
    tidyr::pivot_longer(everything(), names_to = "event", values_to = "n_pairs")
  write_table_pair(summ, results_tables, "candidate_as_event_summary", cfg = cfg)
  message("\n---- event types across the candidate switches ----")
  print(as.data.frame(summ), row.names = FALSE)

  fr_summ <- pair_tab |> count(.data$frame_relationship, .data$coding_change, name = "n_pairs")
  write_table_pair(fr_summ, results_tables, "candidate_coding_consequence_summary", cfg = cfg)
  message("\n---- coding consequence ----")
  print(as.data.frame(fr_summ), row.names = FALSE)
}

write_run_manifest(
  script = "15_isoform_characterization.R", cfg = cfg, root = root,
  extra = list(inputs = list(gff3 = gff_rel, fasta = fa_rel,
                             scoring = "T_UT context layer"))
)
message("\n15_isoform_characterization.R: done")
