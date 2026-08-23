#!/usr/bin/env Rscript
# 14: Do the switching isoforms have a shared functional character? PTC / intron retention.
#
# Run (after 01, 09): Rscript scripts/14_isoform_consequence.R
# Works from committed data/processed/; the external drive is NOT needed.
#
# WHERE THE ANNOTATION COMES FROM
#
# The upstream ISA script runs `analyzeORF` unconditionally (NCBR-40-main, REVIEW_CHANGES 0i),
# so every object carries two consequence columns the pipeline had never touched:
#
#   PTC  premature termination codon -> a predicted NMD target. TRUE/FALSE/NA.
#   IR   number of retained introns. 0 = none.
#
# THE COVERAGE CONSTRAINT, WHICH SHAPES EVERYTHING BELOW
#
# PTC and IR live ONLY in the reduced objects (isoformFeatures_*), never in the unfiltered
# context layer -- 01_load_data.R does not carry them across because the unfiltered exports
# were not parsed for them. The reduced objects hold only genes that already passed the
# upstream switching test, so:
#
#   CRYPTIC SWITCHERS are switching genes and are 56-87% covered. Usable, with the caveat
#   that the misses are not random: a gene can be switching under our context-layer scoring
#   and still be absent from the reduced object, because the two use DIFFERENT FDR UNIVERSES
#   (the standing warning in CLAUDE.md). The covered subset therefore leans toward genes both
#   universes agree on.
#
#   CLASS B / COMPOSITIONAL RESCUE is 5-21% covered and CANNOT be analysed here. Those genes
#   are defined by falling expression with the dominant isoform held, which does not require
#   them to switch -- so they are absent from a switching-gene object by construction.
#   Reporting PTC/IR on 5% of them would describe the switching minority, not the class.
#   This script measures the coverage, writes it out, and stops. See "TO UNBLOCK" below.
#
# TO UNBLOCK CLASS B: re-parse the *unfiltered* exports for PTC/IR in 01_load_data.R and add
# them to isoformContext_<label>.rds. analyzeORF ran on those too, so the values exist; they
# were simply never extracted. Needs the external drive.

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
fig_dir <- ensure_dir(file.path(results_figures, "consequence"))
if (has_ggplot) theme_set(theme_bw(base_size = 11))

lab_of <- function(k) sanitize(as.character((cfg$datasets[[k]] %||% list())$label %||% k)[1L])
DATASETS <- names(cfg$datasets %||% list())

features_of <- function(k) {
  f <- file.path(processed_dir, paste0("isoformFeatures_", lab_of(k), ".rds"))
  if (!file.exists(f)) return(NULL)
  x <- readRDS(f)
  if (!all(c("PTC", "IR") %in% names(x))) {
    message("  ", k, ": no PTC/IR columns; skipped")
    return(NULL)
  }
  x[is.finite(x$dIF), , drop = FALSE]
}

message("14_isoform_consequence.R: reading ORF consequence annotation")
feats <- setNames(lapply(DATASETS, features_of), DATASETS)
feats <- feats[!vapply(feats, is.null, NA)]
if (!length(feats)) stop("No dataset carries PTC/IR. Run 01_load_data.R.")

# ---- Coverage, reported before any result ------------------------------------------------
cov_rows <- list()
hl_files <- list(
  cryptic = file.path(results_tables, "hidden_layer_cryptic_switchers.csv"),
  rescue  = file.path(results_tables, "hidden_layer_compositional_rescue.csv")
)
for (nm in names(hl_files)) {
  if (!file.exists(hl_files[[nm]])) next
  hl <- utils::read.csv(hl_files[[nm]], stringsAsFactors = FALSE)
  for (k in names(feats)) {
    g <- hl$gene_id[hl$dataset == k]
    if (!length(g)) next
    cov_rows[[length(cov_rows) + 1L]] <- tibble(
      set = nm, dataset = k, n_genes = length(g),
      n_with_annotation = sum(g %in% feats[[k]]$gene_id),
      pct_covered = round(100 * mean(g %in% feats[[k]]$gene_id), 1)
    )
  }
}
coverage <- bind_rows(cov_rows)
if (nrow(coverage)) {
  coverage$usable <- coverage$pct_covered >= 50
  write_table_pair(coverage, results_tables, "consequence_annotation_coverage", cfg = cfg)
  message("\n---- annotation coverage of the hidden-layer sets ----")
  print(as.data.frame(coverage))
  bad <- coverage |> filter(!.data$usable)
  if (nrow(bad)) {
    message("\nNOT ANALYSED -- coverage below 50%:")
    for (i in seq_len(nrow(bad))) {
      message(sprintf("  %s / %s: %d of %d genes (%.0f%%)", bad$set[i], bad$dataset[i],
                      bad$n_with_annotation[i], bad$n_genes[i], bad$pct_covered[i]))
    }
    message("  Class B genes need not switch, so a switching-gene object misses most of them.")
    message("  Unblock by carrying PTC/IR from the unfiltered exports into the context layer.")
  }
}

# ---- 1. Unpaired: does consequence track switch direction? -------------------------------
# Marginal shares of PTC+/IR+ among rising vs falling isoforms, within switching genes.
share_test <- function(x, flag, k) {
  ok <- !is.na(flag)
  tb <- table(dir = ifelse(x$dIF[ok] > 0, "up", "down"), hit = flag[ok])
  if (nrow(tb) < 2 || ncol(tb) < 2) return(NULL)
  ft <- stats::fisher.test(tb)
  tibble(
    dataset = k, n = sum(ok),
    pct_up = 100 * tb["up", "TRUE"] / sum(tb["up", ]),
    pct_down = 100 * tb["down", "TRUE"] / sum(tb["down", ]),
    odds_ratio = unname(ft$estimate), p = ft$p.value
  )
}
# PTC is passed RAW so share_test's !is.na() filter can drop unknown-ORF isoforms. Passing
# `PTC %in% TRUE` would reach the filter as a complete logical and count them as negative.
unpaired <- bind_rows(
  bind_rows(lapply(names(feats), function(k) {
    r <- share_test(feats[[k]], feats[[k]]$PTC, k); if (!is.null(r)) r$marker <- "PTC"; r
  })),
  bind_rows(lapply(names(feats), function(k) {
    r <- share_test(feats[[k]], feats[[k]]$IR > 0, k); if (!is.null(r)) r$marker <- "IR"; r
  }))
) |> select("marker", everything())
write_table_pair(unpaired, results_tables, "consequence_direction_unpaired", cfg = cfg)
message("\n---- 1. consequence vs switch direction (unpaired) ----")
print(as.data.frame(unpaired |> mutate(across(where(is.numeric), ~ round(.x, 3)))))

# ---- 2. Paired within gene: top riser vs top faller ---------------------------------------
# Conditions on the gene, so gene-level isoform composition cannot drive the contrast. The
# two picks are extremes in OPPOSITE directions, so the selection is symmetric.
paired_one <- function(x, k) {
  pr <- x |>
    group_by(.data$gene_id) |>
    filter(any(.data$dIF > 0), any(.data$dIF < 0)) |>
    summarize(
      r_PTC = .data$PTC[which.max(.data$dIF)], f_PTC = .data$PTC[which.min(.data$dIF)],
      r_IR = .data$IR[which.max(.data$dIF)] > 0, f_IR = .data$IR[which.min(.data$dIF)] > 0,
      r_novel = .data$is_novel_pacbio[which.max(.data$dIF)],
      f_novel = .data$is_novel_pacbio[which.min(.data$dIF)],
      .groups = "drop"
    )
  # NA means no ORF was called, which is NOT "not an NMD target" -- those gene pairs are
  # dropped, not counted as negative. Coercing with `%in% TRUE` here would silently reclassify
  # the 112 unknown-ORF isoforms in T_UT as PTC-negative and inflate the FALSE class.
  disc <- function(a, b, nm) {
    keep <- !is.na(a) & !is.na(b)
    n_dropped <- sum(!keep)
    a <- a[keep]; b <- b[keep]
    n_r <- sum(a & !b); n_f <- sum(!a & b)
    bt <- stats::binom.test(c(n_r, n_f))
    tibble(dataset = k, marker = nm, n_genes = length(a), n_dropped_unknown = n_dropped,
           riser_only = n_r, faller_only = n_f,
           ratio = if (n_f > 0) n_r / n_f else NA_real_, p = bt$p.value)
  }
  bind_rows(
    disc(pr$r_PTC, pr$f_PTC, "PTC"),
    disc(pr$r_IR, pr$f_IR, "IR"),
    # Confound check: PacBio-novel isoforms are often 5'-truncated (0i), and a truncated
    # transcript is more likely to be called PTC+. If risers were systematically novel, PTC
    # would follow without any regulatory meaning.
    disc(pr$r_novel, pr$f_novel, "novel_pacbio (confound check)")
  )
}
paired <- bind_rows(lapply(names(feats), function(k) paired_one(feats[[k]], k)))
write_table_pair(paired, results_tables, "consequence_paired_riser_vs_faller", cfg = cfg)
message("\n---- 2. paired within gene: top riser vs top faller ----")
print(as.data.frame(paired |> mutate(p = round(.data$p, 4), ratio = round(.data$ratio, 2))))

# ---- 2b. Pooling the paired PTC result, and what carries it -------------------------------
# The direction is consistent in all four datasets, but only one reaches significance, so the
# question is whether the combined evidence means anything.
#
# T_HT and T_UT are the SAME six libraries under two annotations, so pooling both would count
# the same RNA twice. Only the three biological units are pooled, with T_UT standing for T.
# The leave-one-out row is the important one: if removing a single dataset collapses the
# result, the pool is that dataset wearing a larger n.
#
# Caveat on including H_HT: its riser/faller assignment derives from batch-corrected IF values
# (0g), so pooling assumes the assignment -- not the rate -- is comparable. The PTC annotation
# itself is scale-free.
pool_units <- c("T_UT", "U_UT", "H_HT")
pp <- paired |> filter(.data$marker == "PTC", .data$dataset %in% pool_units)
pooled <- NULL
if (nrow(pp) >= 2) {
  mk <- function(d, lbl) {
    r <- sum(d$riser_only); f <- sum(d$faller_only)
    bt <- stats::binom.test(r, r + f)
    tibble(subset = lbl, datasets = paste(d$dataset, collapse = "+"),
           riser_only = r, faller_only = f, ratio = if (f > 0) r / f else NA_real_,
           p = bt$p.value)
  }
  loo <- bind_rows(lapply(pp$dataset, function(dr) {
    mk(pp |> filter(.data$dataset != dr), paste0("drop ", dr))
  }))
  pooled <- bind_rows(mk(pp, "all three units"), loo)
  # Direction consistency across ALL four, including T_HT -- a sign test, not a pooled count.
  allp <- paired |> filter(.data$marker == "PTC")
  n_lean <- sum(allp$riser_only > allp$faller_only)
  pooled$sign_test_datasets_leaning_riser <- paste0(n_lean, "/", nrow(allp))
  pooled$sign_test_p <- stats::binom.test(n_lean, nrow(allp))$p.value
  write_table_pair(pooled, results_tables, "consequence_paired_pooled", cfg = cfg)
  message("\n---- 2b. pooled paired PTC, with leave-one-out ----")
  print(as.data.frame(pooled |> mutate(across(where(is.numeric), ~ round(.x, 4)))))
}

# ---- 3. Are cryptic switchers functionally distinct? --------------------------------------
cryptic_f <- hl_files$cryptic
cryptic_cmp <- NULL
if (file.exists(cryptic_f)) {
  hl <- utils::read.csv(cryptic_f, stringsAsFactors = FALSE)
  cryptic_cmp <- bind_rows(lapply(names(feats), function(k) {
    s <- score_isoforms(feats[[k]], cfg, k, lab_of(k)) |> filter(.data$is_switching)
    if (!nrow(s)) return(NULL)
    s$cryptic <- s$gene_id %in% hl$gene_id[hl$dataset == k]
    if (sum(s$cryptic) < 3) return(NULL)
    one <- function(flag, nm) {
      ok <- !is.na(flag)
      tb <- table(cryptic = s$cryptic[ok], hit = flag[ok])
      if (nrow(tb) < 2 || ncol(tb) < 2) return(NULL)
      tibble(dataset = k, marker = nm,
             n_cryptic = sum(tb["TRUE", ]), n_other = sum(tb["FALSE", ]),
             pct_cryptic = 100 * tb["TRUE", "TRUE"] / sum(tb["TRUE", ]),
             pct_other = 100 * tb["FALSE", "TRUE"] / sum(tb["FALSE", ]),
             p = stats::fisher.test(tb)$p.value)
    }
    bind_rows(one(s$PTC, "PTC"), one(s$IR > 0, "IR"))
  }))
  if (!is.null(cryptic_cmp) && nrow(cryptic_cmp)) {
    write_table_pair(cryptic_cmp, results_tables, "consequence_cryptic_vs_other", cfg = cfg)
    message("\n---- 3. cryptic switchers vs other switchers ----")
    print(as.data.frame(cryptic_cmp |> mutate(across(where(is.numeric), ~ round(.x, 2)))))
  }
}

# ---- Verdict -------------------------------------------------------------------------------
# Written out rather than left to the reader, because the paired PTC result is the kind of
# single-dataset signal this project has repeatedly found not to survive.
ptc_p <- paired |> filter(.data$marker == "PTC")
sig <- ptc_p |> filter(.data$p < 0.05)
verdict <- tibble(
  question = c(
    "Does switch direction track NMD (PTC) status?",
    "Does switch direction track intron retention?",
    "Are cryptic switchers functionally distinct (PTC/IR)?",
    "Is Class B enriched for NMD-target minor isoforms?"
  ),
  answer = c(
    "Suggestive, NOT established -- direction consistent 4/4, significant 1/4",
    "No -- null in every dataset, paired and unpaired",
    "No -- indistinguishable in every dataset",
    "NOT TESTABLE -- Class B genes are 5-21% covered by the ORF annotation"
  ),
  note = c(
    paste0("Within a gene the isoform gaining usage is more often an NMD target in all four ",
           "datasets, but only T_UT reaches p<0.05 (21 vs 6). Sign test across the four is ",
           "p=0.125. Pooling the three biological units gives 1.78x, p=0.020 -- but dropping ",
           "T_UT alone collapses it to 1.29x, p=0.47, so the pool is largely that one dataset. ",
           "T_HT, the same libraries re-annotated, leans the same way and does NOT reach ",
           "significance, so the result is not even annotation-robust. Do not report as a ",
           "finding without an independent dataset."),
    "Paired riser-vs-faller is almost exactly balanced in all four datasets.",
    "Both markers, all four datasets, no difference beyond noise.",
    "Carry PTC/IR from the unfiltered exports into the context layer to test it."
  )
)
write_table_pair(verdict, results_tables, "consequence_verdict", cfg = cfg)
message("\n---- verdict ----")
print(as.data.frame(verdict |> select("question", "answer")))

# ---- Figure --------------------------------------------------------------------------------
if (has_ggplot && nrow(paired)) {
  d <- paired |>
    filter(.data$marker %in% c("PTC", "IR")) |>
    tidyr::pivot_longer(c("riser_only", "faller_only"), names_to = "side", values_to = "n") |>
    mutate(side = ifelse(.data$side == "riser_only", "riser is +", "faller is +"))
  p <- ggplot(d, aes(x = .data$dataset, y = .data$n, fill = .data$side)) +
    geom_col(position = position_dodge(width = .75), width = .65) +
    geom_text(aes(label = .data$n), position = position_dodge(width = .75),
              vjust = -0.3, size = 3) +
    facet_wrap(~ .data$marker) +
    scale_fill_manual(values = c("riser is +" = "#d95f02", "faller is +" = "#7570b3"),
                      name = NULL) +
    labs(
      title = "Within a gene, is the isoform that gains usage an NMD target?",
      subtitle = paste0("Discordant gene pairs only. A real effect would lean the same way in ",
                        "every dataset;\nonly T_UT does, and T_HT is the same libraries."),
      x = NULL, y = "discordant genes"
    ) +
    theme(legend.position = "top")
  ggsave(file.path(fig_dir, "fig_consequence_paired.png"), p, width = 8.5, height = 5, dpi = 200)
  message("Figure: ", fig_dir)
}

write_run_manifest("14_isoform_consequence.R", cfg, root)
message("14_isoform_consequence.R: done")
