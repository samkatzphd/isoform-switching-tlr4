#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------------------
# STANDALONE explorer for the donor-blocked refit and the dual-anchor comparison.
#
#   Rscript explore/23_refit_dual_anchor_explore.R
#
# Needs only dplyr + ggplot2 and the CSVs under results/tables/. No bootstrap, no config, no
# raw counts, no Bioconductor. Writes nothing. Run from the project root.
#
# WHAT THIS IS FOR
#   H_HT is three human donors collected on three different days; T_HT is one THP-1 line split
#   and processed simultaneously. H must be more variable -- and in raw counts it is (median
#   baseline IF SD 0.0451 vs 0.0389). But `importRdata` applied limma::removeBatchEffect to H
#   using 2 surrogate variables, and in the pipeline's STORED values H looks LESS variable than
#   T (0.0232 vs 0.0327), which is impossible. Script 23 refits both systems from raw counts
#   with the replicate as an explicit nested blocking factor, so the two are measured the same
#   way, and then reports them as private sets rather than anchor-and-check.
#
#   Run explore/donor_variance_diagnostic.R for the variance evidence itself.
# ---------------------------------------------------------------------------------------

suppressPackageStartupMessages({ library(dplyr); library(ggplot2) })

TBL <- "results/tables"
rd <- function(stem) {
  f <- file.path(TBL, paste0(stem, ".csv"))
  if (!file.exists(f)) stop("Missing ", f, " -- run: Rscript scripts/23_donor_blocked_refit.R")
  utils::read.csv(f, stringsAsFactors = FALSE)
}
summ   <- rd("refit_summary")
dual   <- rd("refit_dual_anchor")
tiers  <- rd("refit_dual_anchor_tiers")
conc   <- rd("refit_selconf_concordance")
ba     <- rd("refit_before_after_summary")

PRIMARY <- "blocked (nested replicate)"
pretty <- c(H_HT_supplementary = "Primary MDM (H)", T_HT_anchor = "THP-1 (T)")

# =========================================================================================
# CONCLUSION 1 -- blocking closes the H/T gap, and it also INCREASES power
# =========================================================================================
# Removing a real nuisance component (donor, batch) from the residual buys power; it does not
# cost it. Both systems gain calls relative to the unblocked fit, and the switching rates end
# up essentially equal -- so the 2.07%-vs-0.87% gap was the correction, not the biology.

cat("\n== Conclusion 1: switching rate by model ==\n")
summ |>
  transmute(model, dataset = pretty[dataset_label], blocked_on,
            switching_isoforms = n_switching_isoforms,
            switching_genes = n_switching_genes,
            pct_switching_genes = round(pct_switching_genes, 2)) |>
  arrange(desc(model), dataset) |>
  as.data.frame() |> print(row.names = FALSE)

cat("\n   for reference, the pipeline's stored rates were H 2.07% and T 0.87%\n")

cat("\n-- stored vs refit on matched isoforms --\n")
ba |>
  transmute(dataset = pretty[dataset_label], n_matched,
            dIF_correlation = round(dIF_correlation, 3),
            switching_stored = n_switching_stored, switching_refit = n_switching_refit) |>
  as.data.frame() |> print(row.names = FALSE)
cat("   T correlates at ~0.97 (it was never corrected); H at ~0.75 (it was).\n")

p1 <- bind_rows(
  summ |> filter(model == PRIMARY) |>
    transmute(dataset = pretty[dataset_label], source = "donor-blocked refit",
              pct = pct_switching_genes),
  tibble(dataset = unname(pretty[c("H_HT_supplementary", "T_HT_anchor")]),
         source = "pipeline (stored)", pct = c(2.07, 0.87))
) |>
  mutate(source = factor(source, levels = c("pipeline (stored)", "donor-blocked refit"))) |>
  ggplot(aes(x = source, y = pct, fill = source)) +
  geom_col(width = 0.6) +
  geom_text(aes(label = sprintf("%.2f%%", pct)), vjust = -0.35, size = 3) +
  facet_wrap(~dataset, nrow = 1) +
  scale_fill_manual(values = c("pipeline (stored)" = "#9E9E9E",
                               "donor-blocked refit" = "#1F6F54"), guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +
  labs(title = "The H/T switching-rate gap closes under honest blocking",
       x = NULL, y = "% of genes switching") +
  theme_bw(base_size = 11) + theme(axis.text.x = element_text(angle = 15, hjust = 1))
print(p1)

# =========================================================================================
# CONCLUSION 2 -- "not replicated" is not one thing
# =========================================================================================
# Only 8 genes are switching in both systems. Read as a binary that looks like failure. Graded
# by what the other system actually shows, most candidates have real supporting evidence that
# simply does not clear a significance threshold at n = 3.

cat("\n== Conclusion 2: evidence ladder in the other system ==\n")
dual |>
  count(discovery = pretty[discovery], tier, name = "genes") |>
  tidyr::pivot_wider(names_from = discovery, values_from = genes, values_fill = 0L) |>
  as.data.frame() |> print(row.names = FALSE)

cat("\n-- private vs shared --\n")
dual |> distinct(gene_id, set) |> count(set, name = "genes") |>
  as.data.frame() |> print(row.names = FALSE)

p2 <- dual |>
  mutate(discovery = pretty[discovery]) |>
  ggplot(aes(y = tier, fill = discovery)) +
  geom_bar(width = 0.68, position = position_dodge2(padding = 0.12)) +
  scale_fill_manual(values = c("Primary MDM (H)" = "#C1443C", "THP-1 (T)" = "#4F7CAC"),
                    name = "discovered in") +
  scale_y_discrete(limits = rev) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.08))) +
  labs(title = "What the other system shows for each discovery",
       x = "genes", y = NULL) +
  theme_bw(base_size = 11) + theme(legend.position = "top")
print(p2)

# =========================================================================================
# CONCLUSION 3 -- direction replicates; significance does not
# =========================================================================================
# This is the empirical case for grading on DIRECTION rather than on "called in the other
# system". Two nulls matter and both are handled here:
#   - direction agreement alone has a 50% null (a coin flip)
#   - "all three replicates agree" has a 25% null, and the genome-wide rate is close to it,
#     so a selected set must be tested against the dataset's OWN background, not against zero
# Everything here conditions on significance in the discovery system, so it is descriptive.

cat("\n== Conclusion 3: is the direction concordance more than chance? ==\n")
conc |>
  transmute(discovery = pretty[discovery], checked_in = pretty[other_system],
            n_tested,
            pct_same_direction = round(pct_same_direction, 1),
            p_vs_coinflip = signif(p_direction_vs_half, 2),
            pct_all_three_agree = round(pct_all_three_agree, 1),
            background_pct = round(background_pct, 1),
            p_vs_background = signif(p_vs_background, 2)) |>
  as.data.frame() |> print(row.names = FALSE)

cat("\n  Direction agreement runs ~69% against a 50% null in BOTH directions (p < 1e-3).\n")
cat("  For T-discovered candidates, all-three-donor agreement in H is 42.7% against a\n")
cat("  25.7% background (p = 3.6e-4). The reverse is 34.8% vs 29.7% and is not significant.\n")

# =========================================================================================
# THINGS TO TRY
# =========================================================================================

# (a) The candidates worth most: discovered in T, and all three H DONORS agree in direction.
#     These are the ones primary cells actively support rather than merely fail to contradict.
if (FALSE) {
  dual |>
    filter(discovery == "T_HT_anchor", !is.na(other_dIF),
           sign(other_dIF) == sign(own_dIF), other_concordant == 3) |>
    transmute(gene_name, own_dIF = round(own_dIF, 3), other_dIF = round(other_dIF, 3),
              tier, set) |>
    arrange(desc(abs(own_dIF))) |> as.data.frame() |> print(row.names = FALSE)
}

# (b) Outright contradictions: opposite direction in the other system. Treat these as a
#     cell-system difference to explain, not as noise to ignore.
if (FALSE) {
  dual |> filter(grepl("^1 ", tier)) |>
    transmute(discovery, gene_name, own_dIF = round(own_dIF, 3),
              other_dIF = round(other_dIF, 3)) |>
    as.data.frame() |> print(row.names = FALSE)
}

# (c) What blocking changed, gene by gene: switching under the blocked model but not the
#     unblocked one (or vice versa). Needs refit_isoform_level.rds, which is .rds-only
#     because of its size.
if (FALSE) {
  r <- readRDS(file.path(TBL, "refit_isoform_level.rds"))
  r |> group_by(gene_id, gene_name, dataset_label) |>
    summarize(blocked = any(is_switching[model == PRIMARY]),
              unblocked = any(is_switching[model != PRIMARY]), .groups = "drop") |>
    filter(blocked != unblocked) |> count(dataset_label, blocked, unblocked) |>
    as.data.frame() |> print(row.names = FALSE)
}

cat("\nDone. Nothing was written.\n")
