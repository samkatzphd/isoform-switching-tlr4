#!/usr/bin/env Rscript
# 23: Refit H and T from raw counts with the replicate as an explicit blocking factor,
#     then compare the two systems as PRIVATE SETS rather than anchor-and-check.
#
# Run: Rscript scripts/23_donor_blocked_refit.R
# Reads only committed data (counts/H_T_pipe/). No external drive. Needs limma.
#
# WHY THIS EXISTS
#
# `importRdata` found 2 surrogate variables for H_HT and applied limma::removeBatchEffect to
# its expression, deriving all of H's IF/dIF from the corrected matrix (AGENT_CONTEXT fact 4).
# That did not just reduce noise -- it inverted the biology. Measured identically on raw
# counts, three donors collected on three days are MORE variable at baseline than one THP-1
# line split three ways (median IF SD 0.0451 vs 0.0389, ratio 1.16), as they must be. In the
# pipeline's stored replicate IF the ratio is 0.71: H appears LESS variable than T, which is
# impossible. H's stored baseline spread sits 49% below its raw value; T's only 16% below.
# See explore/donor_variance_diagnostic.R.
#
# Consequence: H's 2.07%-vs-0.87% switching rate, its lower recommended abundance floor
# (6.02 vs 10.07), and every count derived from them are artefacts of an over-correction.
# They cannot be fixed by relabelling.
#
# THE FIX, AND WHY PAIRING IS THE BLOCKING FACTOR
#
# Donor identity is largely an INTERCEPT: people differ in baseline isoform ratio more than
# in response. Taking the within-donor paired change, IF(+LPS) - IF(-LPS), cancels that
# intercept exactly, which is what an explicit donor block buys you. It keeps donor out of
# the treatment estimate WITHOUT deleting it from the residual variance used for inference --
# precisely what removeBatchEffect on the whole matrix destroyed. On the paired quantity H is
# only 1.10x noisier than T, so H is not badly underpowered; it was mis-analysed.
#
# T1-T3 were matured, treated and extracted separately, so they are genuine paired batches
# too. The SAME estimator is therefore applied to both systems -- otherwise a donor-blocked H
# would be compared against an ISA-processed T, and any difference would be the method.
#
# MODEL, AND ONE APPROACH THAT WAS TRIED AND DISCARDED
#
# Significance comes from limma-voom + `diffSplice` on scaledTPM. diffSplice tests, for each
# transcript, whether its log-fold-change differs from its own gene's overall log-fold-change
# -- which is exactly differential transcript USAGE rather than differential expression. The
# scale and the DTU filters are script 12's (§0h).
#
# NOT satuRn, although script 12 uses it: satuRn fits a quasi-binomial GLM per transcript, and
# serially (BiocParallel cannot bind its cluster ports in this sandbox) that did not finish on
# ~90,000 transcripts. Worse, the conda Bioconductor stack here has an ABI mismatch --
# Biobase.so fails to load with "Symbol not found: _NO_REFERENCES" -- which takes
# SummarizedExperiment and satuRn with it. limma-voom + diffSplice needs neither, is
# vectorised rather than a per-transcript loop, and is a standard DTU routine. Anyone
# re-running this on a machine with a working satuRn should expect similar but not identical
# calls; the two tests are not the same statistic.
#
# Effect size is reported as mean dIF -- the average of the three paired IF differences --
# because satuRn's estimate is on a log-odds usage scale and the project's thresholds, floors
# and prior results are all on the dIF scale. So: model-based q, interpretable effect size.
#
# REJECTED: fitting the three paired dIF values directly with an intercept-only limma model
# (one-sample moderated t). The dIF estimates it produced were right -- they correlate with
# the pipeline's stored values at r = 0.985 for T -- but it called ZERO switching isoforms in
# either system against 202/276 stored. Reducing each isoform to three ratios discards the
# read-depth information that makes DTU detectable: an isoform with 10,000 reads is weighted
# the same as one with 10, and 2 residual df across ~143,000 BH-corrected tests has no power
# left. A count-based model with gene offsets is the field-standard routine for exactly this
# reason. Do not reinstate the ratio-only fit.

bt <- c("utils/bootstrap.R", file.path("..", "utils", "bootstrap.R"))
b_file <- if (any(f <- vapply(bt, file.exists, NA))) { bt[which(f)[1L]] } else { NA_character_ }
if (is.na(b_file)) {
  stop("Could not find utils/bootstrap.R. Run from project root or set ISOFORM_PROJECT_ROOT.")
}
source(b_file, local = FALSE, chdir = FALSE)

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  for (p in c("limma", "edgeR")) {
    if (!requireNamespace(p, quietly = TRUE)) {
      stop("Install ", p, ": BiocManager::install('", p, "')")
    }
  }
  library(limma)
  library(edgeR)
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
results_figures <- ensure_dir(resolve_path(paths$results_figures %||% "results/figures", root = root))
fig_dir <- ensure_dir(file.path(results_figures, "refit"))

sig <- cfg$significance %||% list()
I_Q    <- as.numeric(sig$isoform_q %||% 0.05)
M_DIF  <- as.numeric(sig$min_abs_dif %||% 0.15)
MIN_GE <- as.numeric(sig$min_gene_expression %||% 12)

LIBS <- list(
  H = list(minus = c("H1_minus_S21_001", "H2_minus_S23_001", "H3_minus_S25_001"),
           plus  = c("H1_plus_S22_001",  "H2_plus_S24_001",  "H3_plus_S26_001"),
           label = "H_HT_supplementary", block = "donor"),
  T = list(minus = c("T1_minus_S15_001", "T2_minus_S17_001", "T3_minus_S19_001"),
           plus  = c("T1_plus_S16_001",  "T2_plus_S18_001",  "T3_plus_S20_001"),
           label = "T_HT_anchor", block = "batch")
)

dtu <- cfg$dtu %||% list()
MIN_COUNT <- as.numeric(dtu$min_count %||% 10)
MIN_SAMPLES <- as.numeric(dtu$min_samples %||% 3)
MIN_GENE_COUNT <- as.numeric(dtu$min_gene_count %||% 20)
FDR <- as.numeric(dtu$fdr %||% 0.05)
N_WORKERS <- max(1L, as.integer(dtu$workers %||% 2))

cnt_dir <- file.path(root, "counts", "H_T_pipe")
iso_f <- file.path(cnt_dir, "RSEM.isoforms.TPM.all_samples.txt")
cnt_f <- file.path(cnt_dir, "RSEM.isoforms.expected_count.all_samples.txt")
for (f in c(iso_f, cnt_f)) if (!file.exists(f)) stop("Missing ", f)
message("Reading TPM and expected counts ...")
tpm <- utils::read.delim(iso_f, check.names = FALSE, stringsAsFactors = FALSE)
cnt <- utils::read.delim(cnt_f, check.names = FALSE, stringsAsFactors = FALSE)
stopifnot(identical(tpm$transcript_id, cnt$transcript_id))
message("  ", nrow(tpm), " isoforms x ", ncol(tpm) - 3L, " libraries")

ANN <- tibble(isoform_id = tpm$transcript_id, gene_id = tpm$gene_id,
              gene_name = tpm$GeneName)

# ---- one combined fit over all 12 HT libraries -------------------------------------------
#
# WHY COMBINED RATHER THAN TWO SEPARATE FITS
#
# Fitting H and T separately leaves 6 samples and 2 residual df each, and satuRn then calls
# essentially nothing (1 and 0 switching isoforms against 272/200 stored). Fitting all twelve
# together estimates dispersion across 12 libraries while still extracting a WITHIN-SYSTEM
# LPS contrast for each, which is what the dual-anchor comparison needs.
#
# WHY THE BLOCKING FACTOR IS NESTED, NOT SHARED
#
# Script 12 hit this and documented it: a 6-level pair factor nests entirely inside the
# system, so `~ 0 + group + pair` is rank deficient. Its workaround was a SHARED 3-level
# `rep3`, which assumes T1 and U1 are a matched pair. Here that would assume donor H1 is
# matched to culture batch T1, which is meaningless -- H's replicates are people and T's are
# flasks. So the block is nested WITHIN system: two dummies per system, giving
# `~ 0 + group + rH2 + rH3 + rT2 + rT3` (8 parameters, 12 libraries, 4 residual df).
# This blocks donor within H and batch within T while assuming nothing across them.
#
# Both the blocked model and an unblocked `~ 0 + group` are fitted, so the cost of honest
# blocking is visible rather than asserted.

build_inputs <- function() {
  libs <- unlist(lapply(LIBS, function(p) c(rbind(p$minus, p$plus))), use.names = FALSE)
  tpm_m <- as.matrix(tpm[, libs, drop = FALSE])
  cnt_m <- as.matrix(cnt[, libs, drop = FALSE])
  libsize <- colSums(cnt_m, na.rm = TRUE)
  mat <- sweep(tpm_m, 2, libsize / 1e6, `*`)          # scaledTPM, as in script 12

  system <- ifelse(startsWith(libs, "H"), "H", "T")
  rep_id <- sub("^[HT]([0-9]).*$", "\\1", libs)
  treat <- ifelse(grepl("_plus_", libs), "plus", "minus")
  meta <- data.frame(
    library = libs, system = system, rep_id = rep_id, treat = treat,
    group = factor(paste(system, treat, sep = "_"),
                   levels = c("H_minus", "H_plus", "T_minus", "T_plus")),
    # Nested replicate dummies, numeric so they cannot alias against `group`. A shared
    # replicate factor would be rank deficient AND would assert that donor H1 pairs with
    # culture batch T1, which is meaningless.
    rH2 = as.integer(system == "H" & rep_id == "2"),
    rH3 = as.integer(system == "H" & rep_id == "3"),
    rT2 = as.integer(system == "T" & rep_id == "2"),
    rT3 = as.integer(system == "T" & rep_id == "3"),
    row.names = libs
  )
  stopifnot(nrow(meta) == 12L, all(table(meta$group) == 3L))

  keep <- rowSums(mat >= MIN_COUNT) >= MIN_SAMPLES
  gene_tot <- rowsum(rowSums(mat), ANN$gene_id)
  keep <- keep & ANN$gene_id %in% rownames(gene_tot)[gene_tot[, 1] >= MIN_GENE_COUNT]
  npg <- table(ANN$gene_id[keep])
  keep <- keep & ANN$gene_id %in% names(npg)[npg >= 2]
  a <- ANN[keep, ]
  m <- round(mat[keep, , drop = FALSE]); rownames(m) <- a$isoform_id
  message("  filter: ", nrow(m), " transcripts across ",
          dplyr::n_distinct(a$gene_id), " multi-isoform genes (of ", nrow(mat), ")")
  list(ann = a, mat = m, meta = meta)
}

#' Per-isoform mean dIF within a system, paired across its three replicates.
dif_for <- function(m, ann, p) {
  gs <- rowsum(m, group = ann$gene_id, reorder = FALSE)
  gs <- gs[match(ann$gene_id, rownames(gs)), , drop = FALSE]
  IFm <- ifelse(gs > 0, m / gs, NA_real_)
  d <- IFm[, p$plus, drop = FALSE] - IFm[, p$minus, drop = FALSE]
  list(mean_dIF = rowMeans(d, na.rm = TRUE),
       sd_dIF = apply(d, 1, sd, na.rm = TRUE),
       n_pos = apply(d, 1, function(z) sum(z > 0, na.rm = TRUE)),
       n_neg = apply(d, 1, function(z) sum(z < 0, na.rm = TRUE)))
}

fits <- (function() {
  message("\nBuilding the combined 12-library inputs ...")
  B <- build_inputs()
  a <- B$ann; m <- B$mat; meta <- B$meta

  dif <- lapply(names(LIBS), function(k) dif_for(m, a, LIBS[[k]]))
  names(dif) <- names(LIBS)

  dge <- edgeR::DGEList(counts = m, genes = as.data.frame(a))
  dge <- edgeR::calcNormFactors(dge)

  fit_formula <- function(form, tag) {
    design <- stats::model.matrix(form, meta)
    r <- qr(design)$rank
    message("  [", tag, "] design ", paste(dim(design), collapse = "x"),
            ", rank ", r, ", residual df ", nrow(design) - r)
    if (r < ncol(design)) {
      stop("[", tag, "] design is rank deficient -- the block is not properly nested.")
    }
    v <- limma::voom(dge, design, plot = FALSE)
    fit <- limma::lmFit(v, design)
    cn <- colnames(design)
    L <- matrix(0, nrow = ncol(design), ncol = 2,
                dimnames = list(cn, c("H_LPS", "T_LPS")))
    L["groupH_plus", "H_LPS"] <- 1; L["groupH_minus", "H_LPS"] <- -1
    L["groupT_plus", "T_LPS"] <- 1; L["groupT_minus", "T_LPS"] <- -1
    cf <- limma::contrasts.fit(fit, L)
    # diffSplice contrasts each transcript's logFC against its own gene's logFC: DTU, not DE.
    ds <- limma::diffSplice(cf, geneid = a$gene_id, verbose = FALSE)

    bind_rows(lapply(names(LIBS), function(k) {
      cn2 <- paste0(k, "_LPS")
      tt <- limma::topSplice(ds, coef = cn2, test = "t", number = Inf, sort.by = "none")
      d <- dif[[k]]
      # topSplice may reorder/rename; align on the transcript id column it carries.
      idcol <- intersect(c("isoform_id", "ExonID", "GeneID"), names(tt))
      key <- if (length(idcol)) as.character(tt[[idcol[1]]]) else rownames(tt)
      i <- match(a$isoform_id, key)
      stopifnot(!anyNA(i))
      tibble(
        model = tag, system = k, dataset_label = LIBS[[k]]$label,
        blocked_on = if (grepl("blocked", tag)) LIBS[[k]]$block else "none",
        isoform_id = a$isoform_id, gene_id = a$gene_id, gene_name = a$gene_name,
        usage_logFC = as.numeric(tt$logFC[i]),
        t = as.numeric(tt$t[i]),
        p_value = as.numeric(tt$P.Value[i]),
        q_value = as.numeric(tt$FDR[i]),
        mean_dIF = as.numeric(d$mean_dIF), sd_dIF = as.numeric(d$sd_dIF),
        n_pos = as.integer(d$n_pos), n_neg = as.integer(d$n_neg)
      )
    }))
  }

  bind_rows(
    fit_formula(~ 0 + group + rH2 + rH3 + rT2 + rT3, "blocked (nested replicate)"),
    fit_formula(~ 0 + group, "unblocked (~0+group)")
  )
})()

# The DTU filter above (min_count in min_samples libraries, gene total, >=2 isoforms) already
# enforces abundance. The pipeline's separate gene-expression floor of 12 is on ISA's scale,
# not scaledTPM, so re-applying the bare number here would be a scale error; it is recorded
# for reference and deliberately NOT used as a second filter.
refit <- fits |>
  mutate(
    n_concordant = pmax(.data$n_pos, .data$n_neg),
    all_pairs_agree = .data$n_concordant == 3L,
    is_switching = !is.na(.data$q_value) & .data$q_value < FDR &
      abs(.data$mean_dIF) >= M_DIF
  )
PRIMARY_MODEL <- "blocked (nested replicate)"
# NOT COMMITTED -- .gitignore excludes this one. It is ~13 MB (two models x ~112,000
# transcripts) and nothing reads it: the figures and explore/23_refit_dual_anchor_explore.R
# work from the summary CSVs below, and re-running this script regenerates it. Committing it
# once took the push pack to 33 MB and the push failed. Everything downstream is derived here,
# so the file is a convenience for interactive digging, not an input.
write_table_pair(refit, results_tables, "refit_isoform_level", cfg = cfg, csv = FALSE)

# ---- per-dataset summary ----------------------------------------------------------------

refit_summary <- refit |>
  group_by(.data$model, .data$dataset_label, .data$blocked_on) |>
  summarize(
    n_isoforms_fitted = dplyr::n(),
    n_switching_isoforms = sum(.data$is_switching),
    n_switching_genes = dplyr::n_distinct(.data$gene_id[.data$is_switching]),
    n_genes_fitted = dplyr::n_distinct(.data$gene_id),
    median_sd_dIF = stats::median(.data$sd_dIF, na.rm = TRUE),
    pct_all_pairs_agree = 100 * mean(.data$all_pairs_agree),
    .groups = "drop"
  ) |>
  mutate(pct_switching_genes = 100 * .data$n_switching_genes / .data$n_genes_fitted)
write_table_pair(refit_summary, results_tables, "refit_summary", cfg = cfg)

message("\n---- refit summary ----")
print(as.data.frame(refit_summary |>
  transmute(model = .data$model, dataset = .data$dataset_label,
            blocked_on = .data$blocked_on, n_fitted = .data$n_isoforms_fitted,
            sw_iso = .data$n_switching_isoforms, sw_genes = .data$n_switching_genes,
            pct_sw_genes = round(.data$pct_switching_genes, 2),
            median_sd_dIF = round(.data$median_sd_dIF, 4),
            pct_3of3_agree = round(.data$pct_all_pairs_agree, 1))), row.names = FALSE)

# ---- before / after for H ----------------------------------------------------------------
# What the over-correction cost, per isoform, on the same isoforms.

processed_dir <- resolve_path(paths$processed_dir %||% "data/processed", root = root)
ba_rows <- list()
for (k in names(LIBS)) {
  lab <- LIBS[[k]]$label
  old <- load_scoring_table(processed_dir, lab, quiet = TRUE)
  if (is.null(old)) next
  o <- tibble(isoform_id = old$isoform_id,
              old_dIF = as.numeric(old$dIF),
              old_q = as.numeric(old$isoform_switch_q_value))
  ba_rows[[lab]] <- refit |>
    filter(.data$dataset_label == lab, .data$model == PRIMARY_MODEL) |>
    select("dataset_label", "isoform_id", "gene_name", "mean_dIF", "q_value",
           "is_switching", "n_concordant") |>
    inner_join(o, by = "isoform_id")
}
before_after <- bind_rows(ba_rows)
write_table_pair(before_after, results_tables, "refit_before_after", cfg = cfg, csv = FALSE)

ba_summary <- before_after |>
  group_by(.data$dataset_label) |>
  summarize(
    n_matched = dplyr::n(),
    dIF_correlation = stats::cor(.data$mean_dIF, .data$old_dIF, use = "complete.obs"),
    median_abs_dIF_refit = stats::median(abs(.data$mean_dIF), na.rm = TRUE),
    median_abs_dIF_stored = stats::median(abs(.data$old_dIF), na.rm = TRUE),
    n_switching_refit = sum(.data$is_switching),
    n_switching_stored = sum(.data$old_q < I_Q & abs(.data$old_dIF) >= M_DIF, na.rm = TRUE),
    .groups = "drop"
  )
write_table_pair(ba_summary, results_tables, "refit_before_after_summary", cfg = cfg)

message("\n---- stored vs refit, matched isoforms ----")
print(as.data.frame(ba_summary |>
  transmute(dataset = .data$dataset_label, n = .data$n_matched,
            dIF_r = round(.data$dIF_correlation, 3),
            sw_stored = .data$n_switching_stored, sw_refit = .data$n_switching_refit)),
  row.names = FALSE)

# ---- dual anchor: private and shared sets --------------------------------------------------
# Neither system is privileged. A gene is reported as shared, H-private or T-private, and for
# the private ones the OTHER system's evidence is graded rather than treated as a pass/fail.

# All of this uses the PRIMARY (blocked) model only.
PRI <- refit |> filter(.data$model == PRIMARY_MODEL)

sw_genes <- function(lab) {
  PRI |> filter(.data$dataset_label == lab, .data$is_switching) |>
    pull(.data$gene_id) |> unique()
}
gH <- sw_genes("H_HT_supplementary"); gT <- sw_genes("T_HT_anchor")
fitted_in <- function(lab) PRI |> filter(.data$dataset_label == lab) |>
  pull(.data$gene_id) |> unique()
fH <- fitted_in("H_HT_supplementary"); fT <- fitted_in("T_HT_anchor")

# For every switching gene in either system, grade the other system's support on the SAME
# isoform. The point of the ladder is that "not called" is not one thing: an isoform that is
# absent, one that moves the other way, and one that moves the same way in all three
# replicates but misses the threshold are completely different kinds of evidence.
grade_other <- function(genes, own_lab, other_lab) {
  own <- PRI |> filter(.data$dataset_label == own_lab, .data$gene_id %in% genes,
                       .data$is_switching) |>
    group_by(.data$gene_id) |>
    slice_max(abs(.data$mean_dIF), n = 1, with_ties = FALSE) |>
    ungroup() |>
    select("gene_id", "gene_name", "isoform_id",
           own_dIF = "mean_dIF", own_q = "q_value", own_concordant = "n_concordant")
  oth <- PRI |> filter(.data$dataset_label == other_lab) |>
    select("isoform_id", other_dIF = "mean_dIF", other_q = "q_value",
           other_concordant = "n_concordant", other_switching = "is_switching")
  own |>
    left_join(oth, by = "isoform_id") |>
    mutate(
      discovery = own_lab,
      tier = dplyr::case_when(
        is.na(.data$other_dIF)                                ~ "0 not quantified",
        sign(.data$other_dIF) != sign(.data$own_dIF)          ~ "1 opposite direction",
        .data$other_concordant < 3L                           ~ "2 same direction, replicates split",
        !.data$other_switching & abs(.data$other_dIF) < M_DIF ~ "3 all replicates agree, below |dIF|",
        !.data$other_switching                                ~ "4 all agree, |dIF| ok, not significant",
        TRUE                                                  ~ "5 switching in both"
      )
    )
}
dual <- bind_rows(
  grade_other(gT, "T_HT_anchor", "H_HT_supplementary"),
  grade_other(gH, "H_HT_supplementary", "T_HT_anchor")
) |>
  mutate(set = dplyr::case_when(
    .data$gene_id %in% gH & .data$gene_id %in% gT ~ "shared",
    .data$gene_id %in% gT ~ "T-private",
    TRUE ~ "H-private"
  )) |>
  # satuRn returns some columns as lists inside rowData; make sure nothing list-valued
  # reaches count()/write.csv, which fail opaquely on them.
  mutate(across(where(is.list), ~ vapply(.x, function(z) as.character(z)[1L], character(1))))
write_table_pair(dual, results_tables, "refit_dual_anchor", cfg = cfg)

dual_counts <- dual |> count(.data$discovery, .data$set, .data$tier, name = "n") |>
  arrange(.data$discovery, .data$tier)
write_table_pair(dual_counts, results_tables, "refit_dual_anchor_tiers", cfg = cfg)

message("\n---- dual anchor ----")
message("  switching genes, H: ", length(gH), " | T: ", length(gT),
        " | shared: ", length(intersect(gH, gT)))
message("  genes fitted in both systems: ", length(intersect(fH, fT)))
message("\n  evidence tier in the OTHER system, per discovery set:")
print(as.data.frame(dual |> count(.data$discovery, .data$tier, name = "n")), row.names = FALSE)

# ---- is the direction concordance more than chance? ---------------------------------------
# With three replicates, "all three agree" happens by chance 2*(1/2)^3 = 25% of the time.
# The genome-wide rate is close to that, so the concordance of a SELECTED set has to be
# tested against the dataset's own background, not against zero. This conditions on
# significance in the discovery system, so it is descriptive -- hence the selconf_ prefix.
bg_rate <- PRI |>
  group_by(.data$dataset_label) |>
  summarize(background_3of3 = mean(.data$all_pairs_agree), n_background = dplyr::n(),
            .groups = "drop")

selconf_concordance <- dual |>
  mutate(other_label = ifelse(.data$discovery == "T_HT_anchor",
                              "H_HT_supplementary", "T_HT_anchor")) |>
  filter(!is.na(.data$other_dIF)) |>
  group_by(discovery = .data$discovery, other_system = .data$other_label) |>
  summarize(
    n_tested = dplyr::n(),
    n_same_direction = sum(sign(.data$other_dIF) == sign(.data$own_dIF)),
    n_all_three_agree = sum(.data$other_concordant == 3L &
                              sign(.data$other_dIF) == sign(.data$own_dIF)),
    .groups = "drop"
  ) |>
  left_join(bg_rate, by = c("other_system" = "dataset_label")) |>
  mutate(
    pct_same_direction = 100 * .data$n_same_direction / .data$n_tested,
    pct_all_three_agree = 100 * .data$n_all_three_agree / .data$n_tested,
    background_pct = 100 * .data$background_3of3,
    # one-sided binomial: is 3/3 concordance commoner in the selected set than at large?
    p_vs_background = mapply(function(x, n, p) {
      stats::binom.test(x, n, p = p, alternative = "greater")$p.value
    }, .data$n_all_three_agree, .data$n_tested, .data$background_3of3),
    # direction agreement alone, against the 50% coin-flip null
    p_direction_vs_half = mapply(function(x, n) {
      stats::binom.test(x, n, p = 0.5, alternative = "greater")$p.value
    }, .data$n_same_direction, .data$n_tested)
  )
write_table_pair(selconf_concordance, results_tables, "refit_selconf_concordance", cfg = cfg)

message("\n---- direction concordance in the OTHER system (selection-conditioned) ----")
print(as.data.frame(selconf_concordance |>
  transmute(discovery = .data$discovery, n = .data$n_tested,
            pct_same_dir = round(.data$pct_same_direction, 1),
            p_dir = signif(.data$p_direction_vs_half, 2),
            pct_3of3 = round(.data$pct_all_three_agree, 1),
            background = round(.data$background_pct, 1),
            p_vs_bg = signif(.data$p_vs_background, 2))), row.names = FALSE)

# ---- figure --------------------------------------------------------------------------------
if (requireNamespace("ggplot2", quietly = TRUE)) {
  library(ggplot2)
  th <- theme_bw(base_size = 11) +
    theme(plot.subtitle = element_text(face = "bold"),
          strip.text = element_text(face = "bold"),
          panel.grid.minor = element_blank())
  pretty <- c(H_HT_supplementary = "Primary MDM (H)", T_HT_anchor = "THP-1 (T)")

  # (a) switching rate: stored vs blocked refit
  rate_df <- bind_rows(
    refit_summary |> filter(.data$model == PRIMARY_MODEL) |>
      transmute(dataset = .data$dataset_label, source = "donor-blocked refit",
                pct = .data$pct_switching_genes),
    tibble(dataset = c("T_HT_anchor", "H_HT_supplementary"),
           source = "pipeline (stored)", pct = c(0.87, 2.07))
  ) |> mutate(dataset = pretty[.data$dataset],
              source = factor(.data$source,
                              levels = c("pipeline (stored)", "donor-blocked refit")))
  pa <- ggplot(rate_df, aes(x = .data$source, y = .data$pct, fill = .data$source)) +
    geom_col(width = 0.6) +
    geom_text(aes(label = sprintf("%.2f%%", .data$pct)), vjust = -0.35, size = 3) +
    facet_wrap(~dataset, nrow = 1) +
    scale_fill_manual(values = c("pipeline (stored)" = "#9E9E9E",
                                 "donor-blocked refit" = "#1F6F54"), guide = "none") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +
    labs(subtitle = "Switching-gene rate: the H/T gap closes under honest blocking",
         x = NULL, y = "% of genes switching") + th +
    theme(axis.text.x = element_text(angle = 15, hjust = 1))

  # (b) the evidence ladder
  pb <- ggplot(dual |> mutate(discovery = pretty[.data$discovery]),
               aes(y = .data$tier, fill = .data$discovery)) +
    geom_bar(width = 0.68, position = position_dodge2(padding = 0.12)) +
    scale_fill_manual(values = c("Primary MDM (H)" = "#C1443C",
                                 "THP-1 (T)" = "#4F7CAC"), name = "discovered in") +
    scale_y_discrete(limits = rev) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.08))) +
    labs(subtitle = "Evidence for each discovery in the other system",
         x = "genes", y = NULL) + th + theme(legend.position = "top")

  # (c) private vs shared
  pc <- ggplot(dual |> distinct(.data$gene_id, .data$set) |> count(.data$set),
               aes(x = stats::reorder(.data$set, -.data$n), y = .data$n, fill = .data$set)) +
    geom_col(width = 0.6) +
    geom_text(aes(label = .data$n), vjust = -0.35, size = 3.2) +
    scale_fill_manual(values = c("shared" = "#1F6F54", "T-private" = "#4F7CAC",
                                 "H-private" = "#C1443C"), guide = "none") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
    labs(subtitle = "Private and shared switching genes", x = NULL, y = "genes") + th

  if (requireNamespace("patchwork", quietly = TRUE)) {
    g <- patchwork::wrap_plots(list(pa, pb, pc), design = "AABB\nCCBB",
                               heights = c(1, 1)) +
      patchwork::plot_annotation(
        title = "Donor-blocked refit of H and T, and the two systems as private sets",
        theme = theme(plot.title = element_text(face = "bold", size = 13)))
    ggsave(file.path(fig_dir, "refit_dual_anchor.png"), plot = g,
           width = 12.5, height = 8.6, dpi = 200, bg = "white")
  } else {
    ggsave(file.path(fig_dir, "refit_rate.png"), plot = pa, width = 7, height = 4, dpi = 200)
    ggsave(file.path(fig_dir, "refit_ladder.png"), plot = pb, width = 8, height = 5, dpi = 200)
  }
  message("Figures written to ", fig_dir)
}

write_run_manifest(script = "23_donor_blocked_refit.R", cfg = cfg, root = root)
message("\n23_donor_blocked_refit.R: done")
