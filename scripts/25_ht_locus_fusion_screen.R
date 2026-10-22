#!/usr/bin/env Rscript
# 25: Is each shortlist locus ONE gene, or two neighbouring genes merged by the assembler?
#     The HT-space equivalent of the screen script 18 ran in UT space.
#
# Run (after 23): Rscript scripts/25_ht_locus_fusion_screen.R
# Needs gff3/H_T.gff3 (603 MB, LOCAL but gitignored). No network, no external drive.
#
# WHY THIS EXISTS
#
# Script 18 found that seven UT-space "switches" were not switches at all: the assembler had
# merged two neighbouring genes into one XLOC, so the two "isoforms" shared ZERO exonic bases
# and the dIF between them was the relative expression of two different genes. PLD6(+FLCN),
# DDX58(=TOPORS+RIGI), IPCEF1(+CNKSR3), THAP2(+TMEM19), MREG(=PECR+MREG), IFITM2(+IFITM1),
# CCL3(+CCL3L3). That screen was never re-run on the HT assembly, so the HT follow-up
# shortlist from 23 carries THAP2 and IFITM2 unresolved -- they can be neither cleared nor
# condemned, and must not be ordered as primers until they are.
#
# THE TEST THAT WORKS, AND THE ONE THAT DID NOT
#
# PRIMARY (what 18 did, reproduced here): take the candidate isoform and the isoform in the
# same locus that moves hardest in the OPPOSITE direction -- the pair whose ratio the dIF
# actually describes -- and count the exonic bases they share. Two isoforms of one gene always
# share some exonic sequence. **Zero shared exonic bases means the "switch" is comparing two
# transcripts that have nothing in common**, i.e. the relative expression of two different
# genes the assembler merged into one XLOC. That is a disqualification, not a caveat.
#
# REJECTED -- connected components. The first version asked whether a locus splits into >=2
# groups of transcripts with no exonic overlap between groups. It returned ONE component for
# all 177 loci, including CCL3 and IPCEF1 which 18 had already confirmed fused. The reason is
# structural: a fused XLOC usually exists BECAUSE a read-through transcript spans both genes,
# and that single transcript bridges the two groups into one component. The test can only ever
# find fusions that the assembler had no reason to make. Do not reinstate it.
#
# SUPPORTING, AND IT OVER-CALLS -- distinct gene symbols per locus. The GFF3 carries a
# gene_name per transcript, and ~30 of the 177 loci carry two different symbols; every
# UT-confirmed fusion is among them (THAP2+TMEM19, FLCN+PLD6, MREG+PECR, DDX58+TOPORS,
# CNKSR3+IPCEF1, CCL3+CCL3L1). But two symbols is not two genes in the damaging sense: it also
# catches hosted miRNAs (IRAK3+MIR6502, MIR922+RUBCN), near-identical paralogs that multi-map
# (APOBEC3A+APOBEC3B, CCL4+CCL4L2, IFITM1+IFITM2) and sense/antisense pairs (IRF1-AS1+RAD50).
# It is reported as a flag to read alongside the primary verdict, never as the verdict.
#
# A NOTE ON NAMES
#
# The shortlist's gene_name comes from the RSEM GeneName column; the GFF3 carries its own
# gene_name per transcript. They disagree for 10 of the 63 shortlist loci, and the
# disagreement is itself diagnostic: the locus the shortlist calls THAP2 is TMEM19 in the
# GFF3, and the one it calls IFITM2 is IFITM1 -- exactly the pairs 18 flagged. Both names are
# kept side by side rather than reconciled.

bt <- c("utils/bootstrap.R", file.path("..", "utils", "bootstrap.R"))
b_file <- if (any(f <- vapply(bt, file.exists, NA))) { bt[which(f)[1L]] } else { NA_character_ }
if (is.na(b_file)) {
  stop("Could not find utils/bootstrap.R. Run from project root or set ISOFORM_PROJECT_ROOT.")
}
source(b_file, local = FALSE, chdir = FALSE)

suppressPackageStartupMessages({ library(dplyr); library(tibble) })

root <- if (exists("PROJECT_ROOT", inherits = FALSE)) get("PROJECT_ROOT", inherits = FALSE) else find_project_root()
args_r <- commandArgs(trailingOnly = TRUE)
cfg <- load_yaml_config(if (length(args_r)) args_r[[1L]] else "config/config.yml", root = root)
paths <- cfg$paths %||% list()
results_tables <- ensure_dir(resolve_path(paths$results_tables %||% "results/tables", root = root))
results_figures <- ensure_dir(resolve_path(paths$results_figures %||% "results/figures", root = root))
fig_dir <- ensure_dir(file.path(results_figures, "fusion"))

GFF <- file.path(root, (cfg$assayability %||% list())$gff3_ht_path %||% "gff3/H_T.gff3")
if (!file.exists(GFF)) {
  message("GFF3 not found at ", GFF, " -- skipping (warn and skip, as 13 does).")
  quit(save = "no", status = 0)
}

# The seven loci script 18 confirmed fused in UT space, for cross-checking.
UT_FUSED <- c("PLD6", "DDX58", "IPCEF1", "THAP2", "MREG", "IFITM2", "CCL3")

dual <- readRDS(file.path(results_tables, "refit_dual_anchor.rds"))
targets <- sort(unique(dual$gene_id))
message("Screening ", length(targets), " loci from the dual-anchor table.")

# ---- parse the GFF3 once, for the target loci only ---------------------------------------
# 2.4M lines / 603 MB: pulling only the needed records with a streaming filter is far cheaper
# than reading the whole file into R.
tf_exon <- tempfile(fileext = ".tsv"); tf_tx <- tempfile(fileext = ".tsv")
tf_ids <- tempfile(fileext = ".txt"); writeLines(targets, tf_ids)

awk_exon <- sprintf(
  "awk -F'\\t' 'NR==FNR{w[$1]=1;next} $3==\"exon\"{gid=\"\";tid=\"\";n=split($9,a,\";\");for(i=1;i<=n;i++){if(a[i]~/^gene_id=/)gid=substr(a[i],9);else if(a[i]~/^Parent=/)tid=substr(a[i],8)} if(gid in w)print gid\"\\t\"tid\"\\t\"$1\"\\t\"$4\"\\t\"$5\"\\t\"$7}' %s %s > %s",
  shQuote(tf_ids), shQuote(GFF), shQuote(tf_exon))
awk_tx <- sprintf(
  "awk -F'\\t' 'NR==FNR{w[$1]=1;next} $3==\"transcript\"{gid=\"\";tid=\"\";cc=\"\";cr=\"\";gn=\"\";n=split($9,a,\";\");for(i=1;i<=n;i++){if(a[i]~/^gene_id=/)gid=substr(a[i],9);else if(a[i]~/^ID=/)tid=substr(a[i],4);else if(a[i]~/^class_code=/)cc=substr(a[i],12);else if(a[i]~/^cmp_ref=/)cr=substr(a[i],9);else if(a[i]~/^gene_name=/)gn=substr(a[i],11)} if(gid in w){gsub(/%%3D/,\"=\",cc);print gid\"\\t\"tid\"\\t\"cc\"\\t\"cr\"\\t\"gn}}' %s %s > %s",
  shQuote(tf_ids), shQuote(GFF), shQuote(tf_tx))
message("  extracting exons ..."); stopifnot(system(awk_exon) == 0L)
message("  extracting transcripts ..."); stopifnot(system(awk_tx) == 0L)

ex <- utils::read.delim(tf_exon, header = FALSE, stringsAsFactors = FALSE,
                        col.names = c("gene_id", "transcript_id", "seqid", "start", "end", "strand"))
tx <- utils::read.delim(tf_tx, header = FALSE, stringsAsFactors = FALSE,
                        col.names = c("gene_id", "transcript_id", "class_code", "cmp_ref", "gene_name"))
message("  ", nrow(ex), " exons, ", nrow(tx), " transcripts")
stopifnot(nrow(ex) > 0, nrow(tx) > 0)

# The GFF3's own name for each locus: the symbol carried by most of its transcripts. Kept
# alongside the shortlist's RSEM-derived name rather than reconciled -- where they disagree,
# the disagreement is itself a fusion signal (THAP2 vs TMEM19, IFITM2 vs IFITM1).
gene_label <- tx |>
  filter(nzchar(.data$gene_name)) |>
  count(.data$gene_id, .data$gene_name, name = "n") |>
  group_by(.data$gene_id) |> slice_max(.data$n, n = 1, with_ties = FALSE) |> ungroup() |>
  select("gene_id", "gene_name")

# ---- distinct gene symbols per locus (supporting flag, over-calls) -----------------------
sym <- tx |>
  filter(nzchar(.data$gene_name), !grepl("^(XLOC|ENSG)", .data$gene_name)) |>
  group_by(.data$gene_id) |>
  summarize(n_symbols = dplyr::n_distinct(.data$gene_name),
            symbols = paste(sort(unique(.data$gene_name)), collapse = "+"),
            .groups = "drop")
message("  loci carrying >1 gene symbol: ", sum(sym$n_symbols > 1), " of ", nrow(sym),
        "  (a FLAG, not a verdict -- see header)")

# ---- PRIMARY TEST: exonic overlap of the pair the dIF actually compares -------------------
exl <- split(ex[, c("start", "end")], ex$transcript_id)

#' Total exonic bases shared by two transcripts.
shared_exonic_bp <- function(a, b) {
  A <- exl[[a]]; B <- exl[[b]]
  if (is.null(A) || is.null(B)) return(NA_integer_)
  tot <- 0L
  for (i in seq_len(nrow(A))) {
    ov <- pmin(A$end[i], B$end) - pmax(A$start[i], B$start) + 1L
    tot <- tot + sum(ov[ov > 0])
  }
  as.integer(tot)
}

iso_path <- file.path(results_tables, "refit_isoform_level.rds")
sl_path <- file.path(results_tables, "followup_shortlist.csv")
if (!file.exists(iso_path) || !file.exists(sl_path)) {
  message("Need refit_isoform_level.rds and followup_shortlist.csv -- run 23 and the ",
          "convergence explorer first. Skipping the pairwise test.")
} else {
  r <- readRDS(iso_path) |> filter(.data$model == "blocked (nested replicate)")
  sl <- utils::read.csv(sl_path, stringsAsFactors = FALSE)
  lab_of <- function(x) ifelse(x == "THP-1 (T)", "T_HT_anchor", "H_HT_supplementary")

  pair <- bind_rows(lapply(seq_len(nrow(sl)), function(k) {
    g <- sl$gene_id[k]; cand <- sl$isoform_id[k]; dir <- sign(sl$own_dIF[k])
    pool <- r |> filter(.data$gene_id == g, .data$isoform_id != cand,
                        .data$dataset_label == lab_of(sl$discovered_in[k]))
    opp <- pool |> filter(sign(.data$mean_dIF) == -dir) |>
      arrange(.data$mean_dIF * dir) |> dplyr::slice(1)
    if (!nrow(opp)) {
      return(tibble(gene_id = g, candidate = cand, partner = NA_character_,
                    partner_dIF = NA_real_, shared_exonic_bp = NA_integer_))
    }
    tibble(gene_id = g, candidate = cand, partner = opp$isoform_id[1],
           partner_dIF = opp$mean_dIF[1],
           shared_exonic_bp = shared_exonic_bp(cand, opp$isoform_id[1]))
  }))

  screened <- sl |>
    left_join(pair, by = c("gene_id", "isoform_id" = "candidate")) |>
    left_join(sym, by = "gene_id") |>
    left_join(gene_label |> rename(gff3_gene_name = "gene_name"), by = "gene_id") |>
    mutate(
      n_symbols = dplyr::coalesce(.data$n_symbols, 1L),
      name_disagrees = !is.na(.data$gff3_gene_name) &
        .data$gff3_gene_name != .data$gene_name,
      verdict = dplyr::case_when(
        is.na(.data$shared_exonic_bp) ~ "no opposing isoform found",
        .data$shared_exonic_bp == 0 ~ "DISQUALIFIED: no shared exonic sequence",
        TRUE ~ "genuine isoform switch"
      ),
      passes_screen = .data$verdict == "genuine isoform switch"
    ) |>
    arrange(.data$passes_screen, dplyr::desc(abs(.data$own_dIF)))
  write_table_pair(screened, results_tables, "followup_shortlist_fusion_screened", cfg = cfg)

  message("\n---- PRIMARY TEST on the ", nrow(screened), "-gene follow-up shortlist ----")
  print(as.data.frame(screened |> count(.data$verdict, name = "genes")), row.names = FALSE)

  message("\n---- DISQUALIFIED ----")
  print(as.data.frame(
    screened |> filter(!.data$passes_screen) |>
      transmute(gene = .data$gene_name, gff3_name = .data$gff3_gene_name,
                symbols = .data$symbols, set = .data$set,
                dIF = round(.data$own_dIF, 3), partner_dIF = round(.data$partner_dIF, 3),
                shared_bp = .data$shared_exonic_bp)
  ), row.names = FALSE)

  message("\n---- survivors, top 15 by |dIF| ----")
  print(as.data.frame(
    screened |> filter(.data$passes_screen) |> arrange(dplyr::desc(abs(.data$own_dIF))) |>
      head(15) |>
      transmute(gene = .data$gene_name, set = .data$set, dIF = round(.data$own_dIF, 3),
                partner_dIF = round(.data$partner_dIF, 3),
                shared_bp = .data$shared_exonic_bp, tier = substr(.data$tier, 1, 28))
  ), row.names = FALSE)

  # figure: shared exonic sequence per shortlist gene
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    library(ggplot2)
    pdf_df <- screened |>
      filter(!is.na(.data$shared_exonic_bp), grepl("^[A-Z]", .data$gene_name)) |>
      arrange(dplyr::desc(abs(.data$own_dIF))) |> head(34) |>
      mutate(gene_name = stats::reorder(.data$gene_name, abs(.data$own_dIF)),
             status = ifelse(.data$passes_screen, "genuine switch",
                             "no shared exonic sequence"))
    # Zero cannot be drawn on a log axis, and a length-1 bar renders as a blank row -- the
    # reader sees absence rather than disqualification, and the red legend key never appears.
    # Draw the passes as bars and the failures as an explicit marker at the axis floor.
    pass_df <- pdf_df |> filter(.data$passes_screen)
    fail_df <- pdf_df |> filter(!.data$passes_screen)
    p <- ggplot(pdf_df, aes(y = .data$gene_name)) +
      geom_col(data = pass_df, aes(x = .data$shared_exonic_bp, fill = .data$status),
               width = 0.7) +
      geom_point(data = fail_df, aes(x = 1.15, colour = .data$status), shape = 4,
                 size = 2.6, stroke = 1.3) +
      geom_text(data = fail_df, aes(x = 1.45, label = "0 shared bp - two genes merged"),
                hjust = 0, size = 2.7, colour = "#C1443C") +
      scale_fill_manual(values = c("genuine switch" = "#1F6F54"), name = NULL) +
      scale_colour_manual(values = c("no shared exonic sequence" = "#C1443C"), name = NULL) +
      scale_x_log10(labels = scales::label_comma(),
                    expand = expansion(mult = c(0.02, 0.08))) +
      labs(title = "HT locus screen: does the compared pair share any exonic sequence?",
           subtitle = paste0("Two isoforms of one gene always share exonic sequence. Zero ",
                             "means the dIF compares two transcripts with\nnothing in common ",
                             "-- two genes merged into one XLOC by the assembler, not an ",
                             "isoform switch."),
           x = "exonic bases shared by candidate and opposing isoform (log)", y = NULL) +
      theme_bw(base_size = 11) +
      theme(plot.title = element_text(face = "bold", size = 12),
            plot.subtitle = element_text(size = 9, colour = "grey30"),
            legend.position = "top", panel.grid.minor = element_blank())
    ggsave(file.path(fig_dir, "ht_fusion_screen.png"), plot = p,
           width = 9.5, height = 8.2, dpi = 200, bg = "white")
    message("\nFigure: ", file.path(fig_dir, "ht_fusion_screen.png"))
  }
}

write_run_manifest(script = "25_ht_locus_fusion_screen.R", cfg = cfg, root = root)
message("\n25_ht_locus_fusion_screen.R: done")
