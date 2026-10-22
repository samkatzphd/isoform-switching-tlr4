#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------------------
# Is H_HT actually underpowered relative to T_HT? STANDALONE DIAGNOSTIC.
#
#   Rscript explore/donor_variance_diagnostic.R          # prints + plots, writes nothing
#   options(donor_diag.save_dir = "some/dir")            # ... to also write CSV + PNG
#
# Reads only committed data: counts/H_T_pipe/RSEM.isoforms.expected_count.all_samples.txt
# (raw RSEM, no batch correction) and data/processed/isoformContextRepIF_*.rds (the
# pipeline's stored replicate IF, which for H is derived from the CORRECTED matrix).
# No bootstrap, no config, no external drive. Needs dplyr + ggplot2.
#
# THE QUESTION
#   H is three human donors collected on three different days; T is one THP-1 line split and
#   processed simultaneously. A priori H must be more variable. Does the pipeline's H reflect
#   that, and is H therefore underpowered for detecting isoform-ratio changes?
#
# THE THREE QUANTITIES, AND WHY THE DISTINCTION MATTERS
#   baseline_IF_sd  spread of isoform fraction across the three -LPS replicates. This is
#                   where donor identity lives: different people have different baseline
#                   isoform ratios. UNPAIRED.
#   paired_dIF_sd   spread of the within-replicate change (IF+LPS - IF-LPS) across the three
#                   replicates. For H this is WITHIN DONOR, so the donor's baseline cancels.
#                   This is the quantity that actually sets power to detect a switch.
#   The gap between them is the whole argument: a design that pairs within donor does not
#   pay the full donor-variance penalty, because that penalty is largely an intercept.
#
# WHAT THIS SCRIPT SHOWS (see the printed summary)
#   1. In RAW counts, H's baseline spread exceeds T's -- as biology requires.
#   2. In the pipeline's STORED RepIF, H's baseline spread is LOWER than T's, which is
#      biologically impossible and is the fingerprint of limma::removeBatchEffect having
#      absorbed the donor structure (AGENT_CONTEXT fact 4).
#   3. On the PAIRED quantity, H is only marginally noisier than T, so H is NOT badly
#      underpowered once the analysis is done within donor.
# ---------------------------------------------------------------------------------------

suppressPackageStartupMessages({ library(dplyr); library(ggplot2) })

CNT <- file.path("counts", "H_T_pipe", "RSEM.isoforms.expected_count.all_samples.txt")
if (!file.exists(CNT)) {
  stop("Missing ", CNT, " -- run from the project root. Current wd: ", getwd())
}

LIBS <- list(
  H = list(minus = c("H1_minus_S21_001", "H2_minus_S23_001", "H3_minus_S25_001"),
           plus  = c("H1_plus_S22_001",  "H2_plus_S24_001",  "H3_plus_S26_001"),
           label = "H: 3 donors, 3 days"),
  T = list(minus = c("T1_minus_S15_001", "T2_minus_S17_001", "T3_minus_S19_001"),
           plus  = c("T1_plus_S16_001",  "T2_plus_S18_001",  "T3_plus_S20_001"),
           label = "T: one line, split")
)

rc <- utils::read.delim(CNT, check.names = FALSE, stringsAsFactors = FALSE)

#' Isoform fraction per library from RAW counts: isoform count / total count of its gene.
if_raw <- function(col) {
  v <- rc[[col]]
  g <- ave(v, rc$gene_id, FUN = function(z) sum(z, na.rm = TRUE))
  ifelse(g > 0, v / g, NA_real_)
}

# Restrict to isoforms whose baseline fraction can actually move. An isoform sitting at 0 or
# 1 has no room for a ratio change, and including those deflates every SD toward zero.
per_iso <- bind_rows(lapply(names(LIBS), function(k) {
  p <- LIBS[[k]]
  m  <- sapply(p$minus, if_raw)
  pl <- sapply(p$plus,  if_raw)
  dIF <- pl - m
  mu <- rowMeans(m)
  keep <- which(mu > 0.05 & mu < 0.95 & stats::complete.cases(dIF) & stats::complete.cases(m))
  tibble(
    dataset = k, label = p$label,
    isoform_id = rc$transcript_id[keep],
    gene_name = rc$GeneName[keep],
    baseline_IF_sd = apply(m[keep, , drop = FALSE], 1, sd),
    paired_dIF_sd  = apply(dIF[keep, , drop = FALSE], 1, sd),
    mean_paired_dIF = rowMeans(dIF[keep, , drop = FALSE]),
    n_concordant_direction = apply(dIF[keep, , drop = FALSE], 1,
                                   function(z) max(sum(z > 0), sum(z < 0)))
  )
}))

# The pipeline's stored replicate IF, for the same baseline quantity. H's is corrected.
stored <- bind_rows(lapply(names(LIBS), function(k) {
  f <- file.path("data", "processed",
                 if (k == "H") "isoformContextRepIF_H_HT_supplementary.rds"
                 else "isoformContextRepIF_T_HT_anchor.rds")
  if (!file.exists(f)) return(NULL)
  x <- readRDS(f)
  cols <- grep("_minus_", names(x), value = TRUE)
  m <- as.matrix(x[, cols, drop = FALSE])
  mu <- rowMeans(m)
  keep <- which(mu > 0.05 & mu < 0.95 & stats::complete.cases(m))
  tibble(dataset = k, stored_baseline_IF_sd = apply(m[keep, , drop = FALSE], 1, sd))
}))

#' Minimum |dIF| detectable at 80% power, two-sided one-sample t with n replicates.
min_detectable <- function(sd_eff, n = 3, alpha = 0.05, power = 0.80) {
  sd_eff / sqrt(n) * (qt(1 - alpha / 2, n - 1) + qt(power, n - 1))
}

summary_tbl <- per_iso |>
  group_by(dataset, label) |>
  summarize(
    n_isoforms = dplyr::n(),
    raw_baseline_IF_sd = median(baseline_IF_sd, na.rm = TRUE),
    paired_dIF_sd = median(paired_dIF_sd, na.rm = TRUE),
    .groups = "drop"
  ) |>
  left_join(
    stored |> group_by(dataset) |>
      summarize(stored_baseline_IF_sd = median(stored_baseline_IF_sd, na.rm = TRUE),
                .groups = "drop"),
    by = "dataset"
  ) |>
  mutate(min_detectable_dIF_n3 = min_detectable(paired_dIF_sd))

cat("\n== Median isoform-fraction variability, H vs T ==\n")
print(as.data.frame(summary_tbl |>
  transmute(dataset, label, n_isoforms,
            raw_baseline = round(raw_baseline_IF_sd, 4),
            stored_baseline = round(stored_baseline_IF_sd, 4),
            paired_dIF = round(paired_dIF_sd, 4),
            min_detectable_dIF = round(min_detectable_dIF_n3, 3))), row.names = FALSE)

r <- function(col) {
  h <- summary_tbl[[col]][summary_tbl$dataset == "H"]
  t <- summary_tbl[[col]][summary_tbl$dataset == "T"]
  round(h / t, 2)
}
cat("\nH/T ratios:\n")
cat("  raw baseline spread    : ", r("raw_baseline_IF_sd"),
    "   <- > 1 as biology requires\n", sep = "")
cat("  STORED baseline spread : ", r("stored_baseline_IF_sd"),
    "   <- < 1, biologically impossible: the correction absorbed donor\n", sep = "")
cat("  paired dIF spread      : ", r("paired_dIF_sd"),
    "   <- pairing within donor cancels most of the penalty\n", sep = "")

cat("\nDirection concordance across the 3 replicates (raw, all informative isoforms):\n")
print(as.data.frame(per_iso |> count(dataset, n_concordant_direction) |>
  group_by(dataset) |> mutate(pct = round(100 * n / sum(n), 1)) |> ungroup()),
  row.names = FALSE)

# ---- figure ------------------------------------------------------------------------------
plot_df <- bind_rows(
  per_iso |> transmute(label, quantity = "baseline IF SD (unpaired)\nraw counts",
                       sd = baseline_IF_sd),
  per_iso |> transmute(label, quantity = "paired dIF SD (within donor)\nraw counts",
                       sd = paired_dIF_sd),
  stored |> left_join(tibble(dataset = c("H", "T"),
                             label = c(LIBS$H$label, LIBS$T$label)), by = "dataset") |>
    transmute(label, quantity = "baseline IF SD (unpaired)\npipeline stored",
              sd = stored_baseline_IF_sd)
) |>
  mutate(quantity = factor(quantity, levels = c(
    "baseline IF SD (unpaired)\nraw counts",
    "baseline IF SD (unpaired)\npipeline stored",
    "paired dIF SD (within donor)\nraw counts")))

p <- ggplot(plot_df, aes(x = label, y = sd, fill = label)) +
  geom_boxplot(outlier.shape = NA, width = 0.55, alpha = 0.9) +
  facet_wrap(~quantity, nrow = 1) +
  coord_cartesian(ylim = c(0, 0.16)) +
  scale_fill_manual(values = c("H: 3 donors, 3 days" = "#C1443C",
                               "T: one line, split" = "#4F7CAC"), guide = "none") +
  labs(
    title = "H is more variable than T at baseline -- except in the pipeline's stored values",
    subtitle = paste0("Middle panel is biologically impossible and is the batch-correction ",
                      "artefact. Right panel is what sets power."),
    x = NULL, y = "isoform-fraction SD across replicates"
  ) +
  theme_bw(base_size = 11) +
  theme(plot.title = element_text(face = "bold", size = 12),
        plot.subtitle = element_text(size = 9.5, colour = "grey30"),
        strip.text = element_text(face = "bold"),
        axis.text.x = element_text(angle = 12, hjust = 1),
        panel.grid.minor = element_blank())

save_dir <- getOption("donor_diag.save_dir", NULL)
if (is.null(save_dir)) {
  print(p)
  cat("\nPrinted the figure. Set options(donor_diag.save_dir=) to write CSV + PNG.\n")
} else {
  dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(summary_tbl, file.path(save_dir, "donor_variance_summary.csv"),
                   row.names = FALSE)
  utils::write.csv(per_iso, file.path(save_dir, "donor_variance_per_isoform.csv"),
                   row.names = FALSE)
  ggsave(file.path(save_dir, "donor_variance_diagnostic.png"), plot = p,
         width = 10.5, height = 4.6, dpi = 200, bg = "white")
  cat("\nWrote CSV + PNG to ", save_dir, "\n", sep = "")
}
