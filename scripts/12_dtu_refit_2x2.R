#!/usr/bin/env Rscript
# 12: Model-based DTU refit -- genotype x treatment as ONE 2x2 interaction.
#
# Run: Rscript scripts/12_dtu_refit_2x2.R
# REQUIRES THE EXTERNAL DRIVE (RSEM count matrices). Warns and exits 0 without it.
#
# WHY THIS EXISTS -- it is the blocking step for the whole UBL5 question.
#
# Everything upstream of this script derives from IsoformSwitchAnalyzeR's condition means
# and replicate isoform fractions, not a fitted model. That is what produced the retracted
# layer asymmetry: the expression estimator averages replicates BEFORE taking an absolute
# value and so tracks truth almost linearly, while the splicing estimator is a sum of
# absolute values and is roughly quadratic near zero. Comparing an upward-biased number to
# a downward-biased one manufactured an asymmetry that inverted under calibration
# (scripts/10, REVIEW_CHANGES.md; corrected gap -0.04 vs -0.27 depending on estimator).
#
# The fix is not a better distance measure. It is to stop comparing two ad-hoc estimators
# and fit one model that returns a per-gene effect on a comparable scale. Because the U-T
# RSEM matrix carries ALL TWELVE libraries -- T1-T3 and U1-U3, each +/- LPS -- the layer
# question becomes a single genotype x treatment interaction rather than a comparison of
# two separately-fit ISA runs. That removes the cross-run normalisation problem at the
# same time.
#
# WHAT THE INTERACTION MEANS HERE
#
#   WT LPS effect   : does isoform usage shift on LPS in wildtype?
#   KO LPS effect   : does it shift in the UBL5 knockout?
#   INTERACTION     : does the LPS-driven shift DIFFER between genotypes?
#
# The interaction is the actual Q2 question and it has never been tested directly. A gene
# significant for the interaction is one where losing UBL5 changes how LPS remodels isoform
# usage -- which is what "UBL5 shapes isoform selection" would have to mean.
#
# TWO DESIGN DECISIONS, BOTH RECORDED ON PURPOSE
#
# 1. COUNT SCALE. ISA's isoformCountMatrix is not RSEM expected counts: column totals match
#    but per-isoform values differ by 0.18-2.3x, because ISA derives counts from abundance
#    (scaledTPM = TPM * libsize / 1e6) rather than carrying RSEM's effective-length-weighted
#    expected counts. For DTU, scaledTPM is the RECOMMENDED scale -- raw expected counts are
#    biased by effective-length differences between isoforms of the SAME gene, which is
#    exactly the comparison DTU makes. This script uses scaledTPM and writes the raw-count
#    result alongside it as a sensitivity check (config: dtu.count_scale). Open question #1
#    asked for this choice to be deliberate and recorded; this is that record.
#
# 2. BLOCKING. The design is paired (REVIEW_CHANGES.md 0e, confirmed with the
#    experimentalist), but a 6-level pair factor is RANK DEFICIENT in the 2x2: pair nests
#    entirely inside genotype, so pairU1+pairU2+pairU3 already equals the U indicator that
#    `group` spans (verified: ncol 9, rank 8). Two models are therefore fitted:
#      primary     ~ 0 + group          -- no blocking, no assumption
#      sensitivity ~ 0 + group + rep3   -- shared 3-level replicate, which ASSUMES T1 and
#                                          U1 were processed as matched batches
#    The assumption behind the sensitivity model has not been confirmed, so the primary
#    model is the one reported. Agreement between them is the check.

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
for (p in c("satuRn", "SummarizedExperiment", "edgeR")) {
  if (!requireNamespace(p, quietly = TRUE)) {
    message("12_dtu_refit_2x2.R: ", p, " not installed; nothing written.")
    quit(save = "no", status = 0L)
  }
}
suppressPackageStartupMessages({
  library(satuRn)
  library(SummarizedExperiment)
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

results_tables <- ensure_dir(resolve_path(paths$results_tables %||% "results/tables", root = root))
results_figures <- ensure_dir(resolve_path(paths$results_figures %||% "results/figures", root = root))
fig_dir <- ensure_dir(file.path(results_figures, "dtu_refit"))

dtu <- cfg$dtu %||% list()
counts_dir <- dtu$counts_dir %||% "/Volumes/Expansion/IsoformSwitchAnalyzer/Counts"
min_count <- as.numeric(dtu$min_count %||% 10)
min_samples <- as.integer(dtu$min_samples %||% 3)
min_gene_count <- as.numeric(dtu$min_gene_count %||% 20)
fdr <- as.numeric(dtu$fdr %||% 0.05)
n_workers <- as.integer(dtu$workers %||% 2)
if (has_ggplot) theme_set(theme_bw(base_size = 11))

iso_dir <- file.path(counts_dir, "U-T", "Isoforms")
f_tpm <- file.path(iso_dir, "RSEM.isoforms.TPM.all_samples.txt")
f_cnt <- file.path(iso_dir, "RSEM.isoforms.expected_count.all_samples.txt")

if (!all(file.exists(f_tpm, f_cnt))) {
  message("12_dtu_refit_2x2.R: RSEM matrices not found under ", iso_dir)
  message("  Mount the external drive and re-run. Nothing written; existing tables left in place.")
  quit(save = "no", status = 0L)
}

# ---- Load ---------------------------------------------------------------------------------
message("Reading RSEM isoform matrices (all twelve U-T libraries)")
tpm <- utils::read.delim(f_tpm, check.names = FALSE, stringsAsFactors = FALSE)
cnt <- utils::read.delim(f_cnt, check.names = FALSE, stringsAsFactors = FALSE)

key <- c("gene_id", "GeneName", "transcript_id")
stopifnot(all(key %in% names(tpm)), all(key %in% names(cnt)))
stopifnot(identical(tpm$transcript_id, cnt$transcript_id))

libs <- setdiff(names(tpm), key)
message("  ", nrow(tpm), " transcripts x ", length(libs), " libraries: ",
        paste(libs, collapse = ", "))

meta <- tibble(
  sample = libs,
  genotype = factor(sub("^([TU]).*$", "\\1", libs), levels = c("T", "U")),
  treat = factor(ifelse(grepl("_plus", libs), "plus", "minus"), levels = c("minus", "plus")),
  rep3 = factor(sub("^[TU]([0-9]).*$", "\\1", libs))
)
meta$group <- factor(paste(meta$genotype, meta$treat, sep = "_"),
                     levels = c("T_minus", "T_plus", "U_minus", "U_plus"))
stopifnot(nrow(meta) == 12L, all(table(meta$group) == 3L))
print(as.data.frame(meta))

tpm_m <- as.matrix(tpm[, libs, drop = FALSE])
cnt_m <- as.matrix(cnt[, libs, drop = FALSE])

# scaledTPM: TPM rescaled so each column sums to that library's mapped-read total. This is
# the DTU-appropriate scale -- see design decision 1 in the header.
libsize <- colSums(cnt_m, na.rm = TRUE)
scaled_m <- sweep(tpm_m, 2, libsize / 1e6, `*`)
message("  library sizes (M reads): ",
        paste(sprintf("%.1f", libsize / 1e6), collapse = ", "))

count_scale <- dtu$count_scale %||% "scaledTPM"
mats <- list(scaledTPM = scaled_m, expected_count = cnt_m)
if (!count_scale %in% names(mats)) stop("dtu.count_scale must be scaledTPM or expected_count")

# ---- Fit one scale ------------------------------------------------------------------------
run_scale <- function(mat, scale_name) {
  message("\n=== ", scale_name, " ===")
  ann <- tibble(isoform_id = tpm$transcript_id, gene_id = tpm$gene_id,
                gene_name = tpm$GeneName)

  # DTU filter: a transcript needs real support, and a gene needs >=2 surviving transcripts
  # or it carries no usage information at all.
  keep_tx <- rowSums(mat >= min_count) >= min_samples
  gene_tot <- rowsum(rowSums(mat), ann$gene_id)
  ok_gene <- rownames(gene_tot)[gene_tot[, 1] >= min_gene_count]
  keep_tx <- keep_tx & ann$gene_id %in% ok_gene
  n_per_gene <- table(ann$gene_id[keep_tx])
  multi <- names(n_per_gene)[n_per_gene >= 2]
  keep_tx <- keep_tx & ann$gene_id %in% multi

  message("  filter: ", sum(keep_tx), " transcripts across ",
          length(unique(ann$gene_id[keep_tx])), " multi-isoform genes",
          " (from ", nrow(mat), " transcripts)")

  m <- round(mat[keep_tx, , drop = FALSE])
  a <- ann[keep_tx, ]
  rownames(m) <- a$isoform_id

  se <- SummarizedExperiment(
    assays = list(counts = m),
    colData = as.data.frame(meta),
    rowData = as.data.frame(a)
  )

  fit_one <- function(formula, tag) {
    design <- model.matrix(formula, colData(se))
    if (qr(design)$rank < ncol(design)) {
      message("  [", tag, "] design is rank deficient; skipped")
      return(NULL)
    }
    message("  [", tag, "] fitting ", nrow(se), " transcripts, design ",
            paste(dim(design), collapse = "x"))
    obj <- satuRn::fitDTU(
      object = se, formula = formula,
      parallel = n_workers > 1,
      BPPARAM = BiocParallel::MulticoreParam(n_workers),
      verbose = FALSE
    )
    cn <- colnames(design)
    gcol <- function(g) { i <- match(paste0("group", g), cn); stopifnot(!is.na(i)); i }
    L <- matrix(0, nrow = ncol(design), ncol = 3,
                dimnames = list(cn, c("WT_LPS", "KO_LPS", "interaction")))
    L[gcol("T_plus"), "WT_LPS"] <- 1;  L[gcol("T_minus"), "WT_LPS"] <- -1
    L[gcol("U_plus"), "KO_LPS"] <- 1;  L[gcol("U_minus"), "KO_LPS"] <- -1
    # (KO_plus - KO_minus) - (WT_plus - WT_minus)
    L[gcol("U_plus"), "interaction"] <- 1;  L[gcol("U_minus"), "interaction"] <- -1
    L[gcol("T_plus"), "interaction"] <- -1; L[gcol("T_minus"), "interaction"] <- 1

    obj <- satuRn::testDTU(object = obj, contrasts = L, diagplot1 = FALSE,
                           diagplot2 = FALSE, sort = FALSE)
    out <- lapply(colnames(L), function(cn2) {
      r <- rowData(obj)[[paste0("fitDTUResult_", cn2)]]
      tibble(
        isoform_id = rownames(obj), gene_id = a$gene_id, gene_name = a$gene_name,
        contrast = cn2,
        estimate = r$estimates, se = r$se, t = r$t,
        p = r$pval,
        # TWO multiple-testing corrections, and they disagree completely here -- see the
        # empirical-null note below. q_bh is the reported one.
        q_bh = r$regular_FDR,
        q_emp = r$empirical_FDR
      )
    })
    bind_rows(out) |> mutate(model = tag, scale = scale_name)
  }

  bind_rows(
    fit_one(~ 0 + group, "primary (~0+group)"),
    fit_one(~ 0 + group + rep3, "sensitivity (+rep3)")
  )
}

res_primary <- run_scale(mats[[count_scale]], count_scale)
alt_scale <- setdiff(names(mats), count_scale)[1L]
res_alt <- run_scale(mats[[alt_scale]], alt_scale)
all_res <- bind_rows(res_primary, res_alt)

if (!nrow(all_res)) stop("No DTU results produced.")

write_table_pair(
  all_res |> arrange(.data$scale, .data$model, .data$contrast, .data$q_bh),
  results_tables, "dtu_refit_isoform_results", cfg = cfg,
  csv = isFALSE((cfg$output %||% list())$csv_twin_isoform_level)
)

# ---- satuRn's empirical null did not fit; record that rather than believing it ------------
# satuRn offers two corrections. `regular_FDR` is BH on the theoretical null. `empirical_FDR`
# re-estimates the null from the z-statistics with locfdr, which is normally the better
# choice because the theoretical null is often too narrow for these tests.
#
# Here it fails. locfdr emits "f(z) misfit" warnings on every contrast and the resulting
# empirical FDR never drops below ~0.30 -- it returns ZERO genes even for the wildtype LPS
# effect, which ISA independently calls for 118 genes and which is the single least
# controversial signal in this dataset. A method that finds nothing where the positive
# control is strongest is not being conservative, it is broken at this sample size: with
# n = 3 per group there are too few informative z-statistics to estimate an empirical null,
# so locfdr inflates the null variance and absorbs the real signal.
#
# So BH is reported and the empirical FDR is carried alongside as a recorded failure. This
# is a limitation of the correction at n = 3, NOT evidence of no effect -- do not read the
# empirical column as a negative result.
emp_check <- all_res |>
  group_by(.data$scale, .data$model, .data$contrast) |>
  summarize(
    min_q_bh = min(.data$q_bh, na.rm = TRUE),
    min_q_emp = min(.data$q_emp, na.rm = TRUE),
    n_iso_bh = sum(.data$q_bh < fdr, na.rm = TRUE),
    n_iso_emp = sum(.data$q_emp < fdr, na.rm = TRUE),
    .groups = "drop"
  ) |>
  mutate(empirical_null_usable = .data$min_q_emp < fdr)
write_table_pair(emp_check, results_tables, "dtu_refit_fdr_comparison", cfg = cfg)
message("\n---- BH vs satuRn empirical FDR (empirical null fails at n=3) ----")
print(as.data.frame(emp_check))

# ---- Gene-level roll-up -------------------------------------------------------------------
# A gene is called for a contrast if any of its isoforms is. satuRn returns per-isoform
# tests; the gene-level question is whether usage was remodelled anywhere in the gene.
gene_res <- all_res |>
  filter(is.finite(.data$q_bh)) |>
  group_by(.data$scale, .data$model, .data$contrast, .data$gene_id, .data$gene_name) |>
  summarize(
    min_q = min(.data$q_bh, na.rm = TRUE),
    min_q_emp = min(.data$q_emp, na.rm = TRUE),
    max_abs_estimate = max(abs(.data$estimate), na.rm = TRUE),
    n_iso = dplyr::n(),
    n_sig_iso = sum(.data$q_bh < fdr, na.rm = TRUE),
    .groups = "drop"
  ) |>
  mutate(significant = .data$min_q < fdr)

write_table_pair(gene_res |> arrange(.data$scale, .data$model, .data$contrast, .data$min_q),
                 results_tables, "dtu_refit_gene_results", cfg = cfg)

summary_tbl <- gene_res |>
  group_by(.data$scale, .data$model, .data$contrast) |>
  summarize(
    genes_tested = dplyr::n(),
    genes_significant = sum(.data$significant),
    pct_significant = round(100 * sum(.data$significant) / dplyr::n(), 2),
    .groups = "drop"
  ) |>
  arrange(.data$scale, .data$model, match(.data$contrast, c("WT_LPS", "KO_LPS", "interaction")))

write_table_pair(summary_tbl, results_tables, "dtu_refit_summary", cfg = cfg)
message("\n---- DTU refit summary ----")
print(as.data.frame(summary_tbl))

# ---- Retention, on ONE scale -- and the honest limits of estimating it --------------------
# This is what the refit was for. The retracted 64%/38% compared an expression estimator to a
# splicing estimator biased the other way. Here both arms come from the SAME model and the
# SAME estimator, so a KO/WT ratio is finally like-for-like. Slope through the origin,
# KO = beta * WT, so beta is the fraction of the WT LPS effect retained in the knockout.
#
# It is fitted two ways because NEITHER is unbiased, and they fail in opposite directions:
#
#   ALL isoforms, no selection. Ruined by REGRESSION DILUTION. The reliability of the WT
#   estimate -- the share of its observed variance that is signal rather than measurement
#   error, (var(est) - mean(se^2)) / var(est) -- is about 0.007 here. That is, ~99% of the
#   spread in the transcriptome-wide WT estimates is noise, because most isoforms have no
#   real LPS effect at n = 3. A slope with a near-zero-reliability predictor is driven to
#   zero regardless of the truth, and the standard correction (divide by reliability) is
#   unstable at this magnitude -- it returns beta > 10, which is meaningless. This column is
#   reported as a DIAGNOSTIC, never as a retention estimate.
#
#   WT-SIGNIFICANT isoforms only. Real signal, but selected on the WT outcome, so the WT
#   estimates carry a winner's curse and are inflated in magnitude. The KO estimates for
#   those same isoforms are not selected on and stay unbiased. An inflated denominator
#   biases the slope DOWNWARD, which makes this figure a conservative LOWER BOUND on
#   retention, not an overestimate.
#
# So: the selected slope is the reportable number and it bounds retention from below.
attenuation <- function(scale_name, model_tag) {
  w <- all_res |>
    filter(.data$scale == scale_name, .data$model == model_tag,
           .data$contrast %in% c("WT_LPS", "KO_LPS")) |>
    select("isoform_id", "contrast", "estimate", "se", "q_bh") |>
    tidyr::pivot_wider(names_from = "contrast",
                       values_from = c("estimate", "se", "q_bh")) |>
    filter(is.finite(.data$estimate_WT_LPS), is.finite(.data$estimate_KO_LPS),
           is.finite(.data$se_WT_LPS))
  if (!nrow(w)) return(NULL)

  v <- stats::var(w$estimate_WT_LPS)
  reliability <- (v - mean(w$se_WT_LPS^2)) / v

  fit_all <- stats::lm(estimate_KO_LPS ~ 0 + estimate_WT_LPS, data = w)
  sel <- w |> filter(.data$q_bh_WT_LPS < fdr)
  if (nrow(sel) < 3) return(NULL)
  fit_sel <- stats::lm(estimate_KO_LPS ~ 0 + estimate_WT_LPS, data = sel)
  ci_sel <- stats::confint(fit_sel)

  tibble(
    scale = scale_name, model = model_tag,
    n_isoforms_wt_sig = nrow(sel),
    retention_lower_bound = unname(coef(fit_sel)[1]),
    ci_low = ci_sel[1, 1], ci_high = ci_sel[1, 2],
    n_isoforms_all = nrow(w),
    wt_estimate_reliability = reliability,
    diagnostic_slope_all = unname(coef(fit_all)[1]),
    note = paste0("retention_lower_bound: WT-significant isoforms; winner's curse on the ",
                  "denominator biases it DOWN, so true retention is at least this. ",
                  "diagnostic_slope_all is regression dilution at reliability ",
                  sprintf("%.3f", reliability), " -- not a retention estimate.")
  )
}
att <- bind_rows(lapply(
  unique(all_res$scale),
  function(s) bind_rows(lapply(unique(all_res$model), function(m) attenuation(s, m)))
))
if (nrow(att)) {
  write_table_pair(att, results_tables, "dtu_refit_attenuation", cfg = cfg)
  message("\n---- LPS-response retention in the KO (single estimator, lower bound) ----")
  print(as.data.frame(att |> select("scale", "model", "n_isoforms_wt_sig",
                                    "retention_lower_bound", "ci_low", "ci_high",
                                    "wt_estimate_reliability", "diagnostic_slope_all")))
}

# ---- The headline: how many genes show a genotype x treatment interaction? ----------------
prim <- gene_res |> filter(.data$scale == count_scale, grepl("^primary", .data$model))
inter <- prim |> filter(.data$contrast == "interaction", .data$significant) |>
  arrange(.data$min_q)
write_table_pair(inter, results_tables, "dtu_refit_interaction_genes", cfg = cfg)

message("\nInteraction genes (genotype x treatment, ", count_scale, ", primary model): ",
        nrow(inter))
if (nrow(inter)) print(as.data.frame(head(inter |> select("gene_name", "min_q",
                                                          "max_abs_estimate", "n_sig_iso"), 15)))

# ---- Figures -------------------------------------------------------------------------------
if (has_ggplot) {
  p1 <- ggplot(summary_tbl, aes(x = .data$contrast, y = .data$genes_significant,
                                fill = .data$model)) +
    geom_col(position = position_dodge(width = .8), width = .7) +
    geom_text(aes(label = .data$genes_significant),
              position = position_dodge(width = .8), vjust = -0.3, size = 3) +
    facet_wrap(~ .data$scale) +
    scale_fill_brewer(palette = "Set2", name = NULL) +
    labs(
      title = "Model-based DTU: genes called per contrast",
      subtitle = paste0("satuRn, BH FDR < ", fdr,
                        "; the interaction is the direct test of a UBL5 splicing effect"),
      x = NULL, y = "significant genes"
    ) +
    theme(legend.position = "top")
  ggsave(file.path(fig_dir, "fig_dtu_refit_counts.png"), p1, width = 9, height = 5, dpi = 200)

  wide <- prim |>
    select("contrast", "gene_id", "gene_name", "min_q") |>
    tidyr::pivot_wider(names_from = "contrast", values_from = "min_q")
  if (all(c("WT_LPS", "KO_LPS") %in% names(wide))) {
    p2 <- ggplot(wide, aes(x = -log10(.data$WT_LPS), y = -log10(.data$KO_LPS))) +
      geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey60") +
      geom_point(alpha = .35, size = 1.4, colour = "#1f78b4") +
      geom_hline(yintercept = -log10(fdr), linetype = 3) +
      geom_vline(xintercept = -log10(fdr), linetype = 3) +
      labs(
        title = "LPS remodelling of isoform usage, wildtype vs UBL5 knockout",
        subtitle = "Gene-level minimum FDR per arm; dotted lines mark the cutoff",
        x = "-log10 FDR, WT LPS effect", y = "-log10 FDR, KO LPS effect"
      )
    ggsave(file.path(fig_dir, "fig_dtu_refit_arms.png"), p2, width = 7.5, height = 6, dpi = 200)
  }
  message("Figures: ", fig_dir)
}

write_run_manifest("12_dtu_refit_2x2.R", cfg, root,
                   extra = list(dtu = list(count_scale = count_scale, fdr = fdr,
                                           min_count = min_count, min_samples = min_samples)))
message("12_dtu_refit_2x2.R: done")
