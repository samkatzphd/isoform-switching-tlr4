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
#   CLASS B / COMPOSITIONAL RESCUE used to be 5-21% covered and was NOT analysable here.
#   Those genes are defined by falling expression with the dominant isoform held, which does
#   not require them to switch -- so they were absent from a switching-gene object by
#   construction. UNBLOCKED 2026-08-27: `extract_context_tables()` now carries PTC and the
#   whole orfAnalysis slot across from the unfiltered exports, so section 4 below reads the
#   CONTEXT layer and covers ~92% of every gene set. See §"Class B" in REVIEW_CHANGES.md.
#
# WHICH SOURCE EACH SECTION USES, AND WHY IT DIFFERS
#
#   Sections 1-3 ask about SWITCHING isoforms, so the reduced object is the right source:
#   those isoforms are in it by definition, and it is the only place `IR` exists.
#
#   Section 4 (Class B) asks about NON-switching isoforms -- the minor transcripts a gene
#   sheds while holding its dominant one. Those exist only in the context layer. It therefore
#   uses PTC alone: `IR` is a column of the reduced objects only, and was not re-derived for
#   the unfiltered universe. IR is null on every switching test above, so its absence here
#   costs the weaker of the two markers.

bt <- c("utils/bootstrap.R", file.path("..", "utils", "bootstrap.R"))
b_file <- if (any(f <- vapply(bt, file.exists, NA))) { bt[which(f)[1L]] } else { NA_character_ }
if (is.na(b_file)) {
  stop("Could not find utils/bootstrap.R. Run from project root or set ISOFORM_PROJECT_ROOT.")
}
source(b_file, local = FALSE, chdir = FALSE)

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(tidyr)
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
# Section 4 scores the context layer itself, so it must apply the same abundance floor every
# other context-scored analysis applies. Sections 1-3 read the reduced objects, which were
# already filtered upstream.
expr_floor <- as.numeric((cfg$significance %||% list())$min_gene_expression %||% 0)

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

# ---- 4. Class B / compositional rescue: is the gene shedding NMD targets? -----------------
# THE CLAIM UNDER TEST. A Class B gene's total output falls, but its baseline-dominant isoform
# holds while the minor isoforms fall harder. One reading is "purification": the transcripts
# being lost are non-productive -- NMD targets -- so the gene is concentrating its functional
# output rather than being switched off. The alternative reading is arithmetic: minor isoforms
# are low-abundance, low-abundance transcripts are noisier and fall more often by chance, and
# PTC status has nothing to do with which ones go.
#
# SOURCE. The context layer, not the reduced objects. Class B genes need not switch, so a
# switching-gene object covers 5-21% of them; the context covers ~92% and is the only source
# in which a gene's NON-switching minor isoforms exist at all.
#
# THREE TESTS, because "enriched" is ambiguous until you name the comparison:
#   4a  WITHIN gene, dominant vs minor. Do minor isoforms carry PTC more often than the
#       dominant one does? Paired, so gene-level composition cancels. But every gene has this
#       contrast -- dominant isoforms are full-length more often everywhere -- so 4a alone
#       cannot distinguish Class B from any gene.
#   4b  BETWEEN gene sets. Is the minor-isoform PTC rate HIGHER in Class B genes than in
#       matched background genes? This is the one that isolates the class.
#   4c  WITHIN Class B genes, does PTC status predict WHICH minor isoforms fall? The sharpest
#       form of the claim: purification requires the NMD targets to be the ones lost, not
#       merely present. Paired within gene, so gene-level decline cancels.
message("\n---- 4. Class B (compositional rescue): PTC status of the shed minor isoforms ----")

rescue_f <- hl_files$rescue
classb <- classb_bg <- classb_fall <- NULL
if (!file.exists(rescue_f)) {
  message("  no hidden_layer_compositional_rescue.csv; run 09_hidden_layer.R first")
} else {
  hl_r <- utils::read.csv(rescue_f, stringsAsFactors = FALSE)

  # Dominance is defined on the UNSTIMULATED profile, exactly as 09 defines it. Re-deriving it
  # here from the same rule (max IF1) rather than trusting a column keeps the two in step; if
  # 09's rule ever changes, this follows it. Baseline dominance is the whole design point --
  # picking the dominant isoform post-LPS would select on the outcome.
  classb_one <- function(k) {
    ctx <- load_scoring_table(processed_dir, lab_of(k), quiet = TRUE)
    if (is.null(ctx) || !"PTC" %in% names(ctx)) return(NULL)
    z <- score_isoforms(ctx, cfg, k, lab_of(k)) |>
      filter(.data$gene_expression >= expr_floor, is.finite(.data$IF1)) |>
      mutate(iso_lfc = log2((as.numeric(.data$iso_value_2) + 1) /
                            (as.numeric(.data$iso_value_1) + 1)))
    if (!nrow(z)) return(NULL)
    z |>
      group_by(.data$gene_id) |>
      # A gene needs at least one minor isoform for "dominant vs minor" to mean anything.
      filter(dplyr::n() > 1L) |>
      # which.max returns ONE index, so a tie at max IF1 cannot make two isoforms dominant
      # -- the same with_ties = FALSE behaviour 09 uses.
      mutate(is_dominant = seq_len(dplyr::n()) == which.max(.data$IF1)) |>
      ungroup() |>
      mutate(
        dataset = k,
        class_b = .data$gene_id %in% hl_r$gene_id[hl_r$dataset == k]
      )
  }
  cb <- bind_rows(lapply(names(feats), classb_one))

  if (is.null(cb) || !nrow(cb)) {
    message("  no context layer carries PTC; run 01_load_data.R with the drive mounted")
  } else {
    # Coverage of the class in its new source, reported before any test.
    cov_b <- cb |>
      filter(.data$class_b) |>
      group_by(.data$dataset) |>
      summarize(n_genes = dplyr::n_distinct(.data$gene_id),
                n_isoforms = dplyr::n(),
                pct_isoforms_with_PTC_call = round(100 * mean(!is.na(.data$PTC)), 1),
                .groups = "drop")
    message("  context coverage of Class B genes:")
    print(as.data.frame(cov_b))

    # ---- 4a. within gene: dominant vs minor ------------------------------------------------
    pair_dom_minor <- function(d, lbl) {
      g <- d |>
        group_by(.data$dataset, .data$gene_id) |>
        summarize(
          dom_ptc = suppressWarnings(.data$PTC[.data$is_dominant][1]),
          any_minor_ptc = if (any(!.data$is_dominant & !is.na(.data$PTC))) {
            any(.data$PTC[!.data$is_dominant], na.rm = TRUE)
          } else {
            NA
          },
          .groups = "drop"
        ) |>
        filter(!is.na(.data$dom_ptc), !is.na(.data$any_minor_ptc))
      if (!nrow(g)) return(NULL)
      g |>
        group_by(.data$dataset) |>
        summarize(
          set = lbl, n_genes = dplyr::n(),
          # Discordant pairs only -- concordant genes carry no information about direction.
          minor_only = sum(!.data$dom_ptc & .data$any_minor_ptc),
          dom_only = sum(.data$dom_ptc & !.data$any_minor_ptc),
          p = stats::binom.test(c(sum(!.data$dom_ptc & .data$any_minor_ptc),
                                  sum(.data$dom_ptc & !.data$any_minor_ptc)))$p.value,
          .groups = "drop"
        )
    }
    classb <- bind_rows(
      pair_dom_minor(cb |> filter(.data$class_b), "Class B"),
      pair_dom_minor(cb |> filter(!.data$class_b), "background")
    )
    if (!is.null(classb) && nrow(classb)) {
      classb <- classb |>
        mutate(ratio_minor_over_dom = ifelse(.data$dom_only > 0,
                                             .data$minor_only / .data$dom_only, NA_real_)) |>
        arrange(.data$dataset, .data$set)
      write_table_pair(classb, results_tables, "consequence_class_b_dominant_vs_minor", cfg = cfg)
      message("\n  4a. within gene -- is the minor side more often PTC+ than the dominant one?")
      print(as.data.frame(classb |> mutate(across(where(is.numeric), ~ round(.x, 4)))))
    }

    # ---- 4b. between gene sets: minor-isoform PTC rate --------------------------------------
    # The comparison that actually isolates the class. Restricted to MINOR isoforms in both
    # sets, so the dominant/minor asymmetry of 4a cannot leak in.
    classb_bg <- cb |>
      filter(!.data$is_dominant, !is.na(.data$PTC)) |>
      group_by(.data$dataset) |>
      group_modify(function(d, key) {
        tb <- table(class_b = d$class_b, ptc = d$PTC)
        if (nrow(tb) < 2 || ncol(tb) < 2) return(tibble())
        ft <- stats::fisher.test(tb)
        tibble(
          n_minor_class_b = sum(tb["TRUE", ]), n_minor_background = sum(tb["FALSE", ]),
          pct_ptc_class_b = 100 * tb["TRUE", "TRUE"] / sum(tb["TRUE", ]),
          pct_ptc_background = 100 * tb["FALSE", "TRUE"] / sum(tb["FALSE", ]),
          odds_ratio = unname(ft$estimate), p = ft$p.value
        )
      }) |>
      ungroup()
    if (nrow(classb_bg)) {
      write_table_pair(classb_bg, results_tables, "consequence_class_b_vs_background", cfg = cfg)
      message("\n  4b. between sets -- minor-isoform PTC rate, Class B vs background")
      print(as.data.frame(classb_bg |> mutate(across(where(is.numeric), ~ round(.x, 3)))))
    }

    # ---- 4c. within Class B: does PTC predict which minor isoforms are shed? ----------------
    # THE TRAP THIS TEST FELL INTO FIRST, recorded because the naive version looked like a
    # clean four-of-four result in the WRONG direction and would have been reported as one.
    #
    # The obvious response variable is iso_lfc = log2((v2+1)/(v1+1)). Class B minor isoforms
    # sit at 0.07-0.40 expression units, so the +1 pseudocount DOMINATES the ratio: an isoform
    # at 0.07 cannot produce a fold-change of any size, while one at 0.40 can. Section 4d
    # below shows PTC+ minors are 2-6x LESS abundant than PTC- minors in every dataset. The
    # naive test therefore reported "PTC+ isoforms fall less" (median delta +0.08 to +0.20,
    # p = 0.002-0.07) when what it had measured was that they start lower. Do not restore it.
    #
    # THE FIX, two parts:
    #   response  dIF = IF2 - IF1, a difference of two ratios ISA computes without any
    #             pseudocount, expressed relative to IF1 so it reads as proportional change in
    #             usage. An isoform at any abundance can lose all of its usage share.
    #   adjust    baseline abundance enters the model explicitly, and gene enters as a fixed
    #             effect so only the WITHIN-gene contrast is estimated. Both the unadjusted
    #             and adjusted coefficients are reported side by side: if they disagree, the
    #             disagreement is the result.
    minors <- cb |>
      filter(.data$class_b, !.data$is_dominant, !is.na(.data$PTC),
             is.finite(.data$IF1), .data$IF1 > 0, is.finite(.data$dIF)) |>
      mutate(rel_dIF = .data$dIF / .data$IF1)

    # ---- 4d. the confound itself, reported whether or not 4c is run ------------------------
    classb_conf <- minors |>
      group_by(.data$dataset, .data$PTC) |>
      summarize(n = dplyr::n(),
                median_expr_baseline = median(as.numeric(.data$iso_value_1), na.rm = TRUE),
                median_IF1 = median(.data$IF1, na.rm = TRUE), .groups = "drop") |>
      tidyr::pivot_wider(names_from = "PTC", values_from = c("n", "median_expr_baseline",
                                                             "median_IF1")) |>
      mutate(abundance_ratio_nonPTC_over_PTC =
               .data$median_expr_baseline_FALSE / .data$median_expr_baseline_TRUE)
    write_table_pair(classb_conf, results_tables, "consequence_class_b_abundance_confound",
                     cfg = cfg)
    message("\n  4d. CONFOUND -- baseline abundance of PTC+ vs PTC- minor isoforms")
    print(as.data.frame(classb_conf |> mutate(across(where(is.numeric), ~ round(.x, 3)))))
    message("     PTC+ minors are the LESS abundant class everywhere. Any pseudocount-based")
    message("     fold-change comparison between the two groups measures this, not NMD.")

    classb_fall <- bind_rows(lapply(sort(unique(minors$dataset)), function(k) {
      d <- minors |> filter(.data$dataset == k)
      # Only genes carrying both a PTC+ and a PTC- minor contribute to a within-gene contrast.
      keep <- d |> group_by(.data$gene_id) |>
        filter(any(.data$PTC), any(!.data$PTC)) |> ungroup()
      if (dplyr::n_distinct(keep$gene_id) < 5) return(NULL)
      unadj <- stats::lm(rel_dIF ~ PTC + factor(gene_id), data = keep)
      adj <- stats::lm(rel_dIF ~ PTC + log10(IF1) + factor(gene_id), data = keep)
      co <- function(m) {
        r <- summary(m)$coefficients
        rn <- grep("^PTCTRUE$", rownames(r), value = TRUE)
        if (!length(rn)) return(c(NA_real_, NA_real_))
        c(r[rn, "Estimate"], r[rn, "Pr(>|t|)"])
      }
      u <- co(unadj); a <- co(adj)
      tibble(
        dataset = k,
        n_genes_informative = dplyr::n_distinct(keep$gene_id),
        n_isoforms = nrow(keep),
        median_rel_dIF_ptc = median(keep$rel_dIF[keep$PTC]),
        median_rel_dIF_nonptc = median(keep$rel_dIF[!keep$PTC]),
        ptc_coef_unadjusted = u[1], p_unadjusted = u[2],
        ptc_coef_abundance_adjusted = a[1], p_adjusted = a[2]
      )
    }))
    if (!is.null(classb_fall) && nrow(classb_fall)) {
      write_table_pair(classb_fall, results_tables, "consequence_class_b_ptc_decline", cfg = cfg)
      message("\n  4c. within Class B -- does PTC predict proportional LOSS of usage?")
      message("      (negative coefficient = PTC+ minors shed MORE usage = purification)")
      print(as.data.frame(classb_fall |> mutate(across(where(is.numeric), ~ round(.x, 4)))))
      message("      Gene enters both models as a fixed effect, so a gene's overall decline")
      message("      cancels and only the within-gene PTC contrast is estimated.")
    }
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
    "No -- not enrichment, and not the transcripts being shed (answered 2026-08-27)"
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
    paste0("Testable since the context layer gained PTC (86-93% coverage vs 5-21%). Three ",
           "readings, none supporting purification. (a) Minor isoforms ARE more often PTC+ ",
           "than the dominant one -- but so are they in background genes (OR 75-319), so it ",
           "is a property of being minor, not of Class B. (b) Minor-isoform PTC RATE in ",
           "Class B vs background: 14.4/14.9% vs ~20% in T_HT/T_UT (depleted, p=0.012/0.042) ",
           "and null in H_HT/U_UT -- and T_HT/T_UT are the SAME six libraries re-annotated, ",
           "so that is one observation, not two. No enrichment anywhere. (c) Whether PTC ",
           "predicts WHICH minors are shed: null in all four (p=0.27-0.98) once measured on ",
           "dIF/IF1 and adjusted for baseline abundance. An earlier version of (c) used ",
           "iso_lfc and reported a consistent four-of-four effect in the OPPOSITE direction; ",
           "that was the +1 pseudocount acting on a class that is 1.5-4.1x less abundant ",
           "(consequence_class_b_abundance_confound.csv), not biology.")
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
