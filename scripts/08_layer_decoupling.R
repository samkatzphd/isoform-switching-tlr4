#!/usr/bin/env Rscript
# 08: Are the expression and splicing layers of the LPS response decoupled in the KO?
#
# Run (after 01): Rscript scripts/08_layer_decoupling.R
# Works from committed data/processed/; the external drive is not needed.
#
# THE QUESTION
#
# Two models make the same qualitative prediction and must be separated quantitatively:
#   (A) downstream  -- UBL5 loss blunts the transcriptional response; less isoform switching
#                      follows mechanically, because switches need a condition-driven shift.
#   (B) upstream    -- UBL5 is a spliceosome-associated regulator and splicing is part of how
#                      the LPS response is executed; the expression deficit is itself partly
#                      a consequence of impaired splicing.
#
# Under (A) both layers should degrade together. Under (B) the splicing layer should degrade
# disproportionately. So the discriminating quantity is the RATIO OF RETENTION between
# layers -- which only means something if both layers are measured on the same footing.
#
# WHY THE OBVIOUS COMPARISON IS INVALID
#
# Expression response (|log2FC|) is unbounded and its noise is roughly symmetric. Splicing
# response (total variation distance over isoform fractions, 0.5 * sum|dIF|) is bounded in
# [0,1] and built from absolute values, so noise can only ever inflate it. Comparing raw
# "percent retained" across the two therefore compares quantities with different noise
# geometry, and any noise correction is a subtraction of two similar numbers whose ratio is
# violently sensitive to how the null was built.
#
# Concretely, a null built from within-condition replicate PAIRS gives a mean TVD of ~0.146
# against an observed condition TVD of ~0.104 -- the null exceeds the signal, because the
# condition means average three replicates while the replicate pairs do not. Corrections
# built on that null are meaningless and can even reverse the sign of the difference.
#
# WHAT THIS SCRIPT DOES INSTEAD
#
# The design is PAIRED: T1_minus and T1_plus are the same biological sample before and after
# LPS, and likewise for T2/T3, U1-U3, H1-H3. The correct null therefore permutes the
# condition label WITHIN each pair, not across all six libraries. That removes between-sample
# variance from the comparison, which is the whole point of having paired samples.
#
# Consequences, both of which shape the implementation:
#
# 1. There are 2^3 = 8 sign patterns, collapsing to 4 distinct once a global flip is quotiented
#    out (a magnitude statistic is unchanged by negating every pair). One is the truth, so
#    only 3 null patterns exist. That is far too few to estimate a per-gene SD, so the
#    z-score standardisation used for the unpaired null is not available. This script uses
#    EXCESS OVER THE NULL MEAN instead -- observed minus the mean of the 3 nulls -- which is
#    identically constructed for both layers and needs no variance estimate.
#
# 2. Per-pair TVD is invariant to swapping that pair, because it is already an absolute
#    value; using it would produce a null identical to the observed statistic. The splicing
#    statistic must therefore be the TVD of the MEAN per-pair dIF vector, i.e.
#    0.5 * sum_j |mean_i (s_i * dIF_ij)|, which does respond to the sign pattern.
#
# The unpaired result is reported alongside so the cost of having ignored the pairing is
# visible rather than asserted.

bt <- c("utils/bootstrap.R", file.path("..", "utils", "bootstrap.R"))
b_file <- if (any(f <- vapply(bt, file.exists, NA))) { bt[which(f)[1L]] } else { NA_character_ }
if (is.na(b_file)) {
  stop("Could not find utils/bootstrap.R. Run from project root or set ISOFORM_PROJECT_ROOT.")
}
source(b_file, local = FALSE, chdir = FALSE)

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(ggplot2)
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

processed_dir <- resolve_path(paths$processed_dir %||% "data/processed", root = root)
results_tables <- ensure_dir(resolve_path(paths$results_tables %||% "results/tables", root = root))
results_figures <- ensure_dir(resolve_path(paths$results_figures %||% "results/figures", root = root))
fig_dir <- ensure_dir(file.path(results_figures, "layers"))

expr_floor <- as.numeric((cfg$significance %||% list())$min_gene_expression %||% 0)
theme_set(theme_bw(base_size = 11))

lab_of <- function(k) sanitize(as.character((cfg$datasets[[k]] %||% list())$label %||% k)[1L])

#' Per-sample gene expression and isoform fraction, from the replicate expression matrix.
#' IF is recomputed from expression rather than taken from repIF so that both layers derive
#' from one source and the permutation applies identically to each.
sample_level <- function(k) {
  lab <- lab_of(k)
  ctx <- load_context_table(processed_dir, lab, "features")
  re <- load_context_table(processed_dir, lab, "rep_expr")
  if (is.null(ctx) || is.null(re)) return(NULL)
  map <- ctx |> select(.data$isoform_id, .data$gene_id)
  m <- inner_join(map, re, by = "isoform_id")
  samples <- setdiff(names(m), c("isoform_id", "gene_id"))
  cond <- ifelse(grepl("_minus", samples, fixed = TRUE), "A",
                 ifelse(grepl("_plus", samples, fixed = TRUE), "B", NA))
  if (anyNA(cond)) return(NULL)
  list(iso = m, samples = samples, cond = cond, gene_ids = unique(m$gene_id))
}

#' Observed and paired-permutation-null statistics for both layers.
#'
#' @param paired TRUE for sign-flips within pairs (correct for this design); FALSE for the
#'   unpaired 3-vs-3 split null, retained only for comparison.
layer_stats <- function(sl, k, paired = TRUE) {
  m <- sl$iso
  samples <- sl$samples
  gexp <- m |>
    group_by(.data$gene_id) |>
    summarize(across(all_of(samples), ~ sum(.x, na.rm = TRUE)), .groups = "drop")
  gmat <- as.matrix(gexp[, samples, drop = FALSE])
  rownames(gmat) <- gexp$gene_id
  gidx <- match(m$gene_id, gexp$gene_id)
  ifmat <- ifelse(gmat[gidx, , drop = FALSE] > 0,
                  as.matrix(m[, samples, drop = FALSE]) / gmat[gidx, , drop = FALSE], 0)

  if (paired) {
    pair <- sub("_(minus|plus)_S[0-9]+$", "", samples)
    pair <- sub("_(minus|plus)$", "", pair)
    is_plus <- grepl("_plus", samples, fixed = TRUE)
    upair <- unique(pair)
    if (length(upair) < 2L || !all(table(pair) == 2L)) {
      stop("[", k, "] could not resolve a paired design from sample names: ",
           paste(samples, collapse = ", "))
    }
    # per-pair difference (plus - minus)
    lfc_p <- vapply(upair, function(u) {
      a <- samples[pair == u & !is_plus]; b <- samples[pair == u & is_plus]
      log2((gmat[, b] + 1) / (gmat[, a] + 1))
    }, numeric(nrow(gmat)))
    dif_p <- vapply(upair, function(u) {
      a <- samples[pair == u & !is_plus]; b <- samples[pair == u & is_plus]
      ifmat[, b] - ifmat[, a]
    }, numeric(nrow(ifmat)))

    # sign patterns, quotienting out the global flip by fixing the first pair to +1
    n_p <- length(upair)
    grid <- as.matrix(expand.grid(rep(list(c(1, -1)), n_p - 1L)))
    signs <- cbind(1, grid)
    stat_for <- function(sv) {
      e <- abs(as.vector(lfc_p %*% sv) / n_p)
      d <- as.vector(dif_p %*% sv) / n_p
      list(lfc = e, tvd = tapply(abs(d), m$gene_id, function(x) 0.5 * sum(x, na.rm = TRUE))[rownames(gmat)])
    }
    is_true <- apply(signs, 1, function(sv) all(sv == 1))
    obs <- stat_for(signs[which(is_true), ])
    nulls <- lapply(which(!is_true), function(i) stat_for(signs[i, ]))
  } else {
    a_true <- samples[sl$cond == "A"]
    splits <- utils::combn(samples, 3, simplify = FALSE)
    seen <- character(); keep <- list()
    for (sp in splits) {
      ckey <- paste(sort(setdiff(samples, sp)), collapse = "|")
      if (!(ckey %in% seen)) {
        keep[[length(keep) + 1L]] <- sp
        seen <- c(seen, paste(sort(sp), collapse = "|"))
      }
    }
    is_true <- vapply(keep, function(sp) setequal(sp, a_true), NA)
    stat_for <- function(grpA) {
      grpB <- setdiff(samples, grpA)
      lfc <- abs(log2((rowMeans(gmat[, grpB, drop = FALSE]) + 1) /
                        (rowMeans(gmat[, grpA, drop = FALSE]) + 1)))
      d <- rowMeans(ifmat[, grpB, drop = FALSE]) - rowMeans(ifmat[, grpA, drop = FALSE])
      list(lfc = lfc, tvd = tapply(abs(d), m$gene_id, function(x) 0.5 * sum(x, na.rm = TRUE))[rownames(gmat)])
    }
    obs <- stat_for(keep[[which(is_true)]])
    nulls <- lapply(keep[!is_true], stat_for)
  }

  null_lfc <- do.call(cbind, lapply(nulls, `[[`, "lfc"))
  null_tvd <- do.call(cbind, lapply(nulls, `[[`, "tvd"))
  tibble(
    dataset = k,
    gene_id = rownames(gmat),
    gene_expression = rowMeans(gmat) * 2,
    obs_lfc = obs$lfc, null_lfc_mean = rowMeans(null_lfc, na.rm = TRUE),
    null_lfc_sd = apply(null_lfc, 1, stats::sd, na.rm = TRUE),
    obs_tvd = obs$tvd, null_tvd_mean = rowMeans(null_tvd, na.rm = TRUE),
    null_tvd_sd = apply(null_tvd, 1, stats::sd, na.rm = TRUE),
    # Excess over the null mean: identically built for both layers, no variance estimate
    # needed, which matters because the paired null has only 3 patterns.
    zE = obs$lfc - rowMeans(null_lfc, na.rm = TRUE),
    zS = obs$tvd - rowMeans(null_tvd, na.rm = TRUE),
    n_null_patterns = ncol(null_lfc),
    null_type = if (paired) "paired (sign-flip within pairs)" else "unpaired (3v3 split)"
  ) |>
    filter(.data$gene_expression >= expr_floor)
}

message("Computing paired permutation nulls (sign-flips within pairs) ...")
sl_t <- sample_level("T_UT"); sl_u <- sample_level("U_UT")
if (is.null(sl_t) || is.null(sl_u)) stop("Need context + replicate expression for both UT arms. Run 01.")
st <- layer_stats(sl_t, "T_UT", paired = TRUE)
su <- layer_stats(sl_u, "U_UT", paired = TRUE)
st_up <- layer_stats(sl_t, "T_UT", paired = FALSE)
su_up <- layer_stats(sl_u, "U_UT", paired = FALSE)
message("  paired null patterns: ", st$n_null_patterns[1], " | unpaired: ", st_up$n_null_patterns[1])

j <- inner_join(
  st |> select(.data$gene_id, zE_T = .data$zE, zS_T = .data$zS,
               oE_T = .data$obs_lfc, oS_T = .data$obs_tvd, nS_T = .data$null_tvd_mean),
  su |> select(.data$gene_id, zE_U = .data$zE, zS_U = .data$zS,
               oE_U = .data$obs_lfc, oS_U = .data$obs_tvd, nS_U = .data$null_tvd_mean),
  by = "gene_id"
) |>
  filter(is.finite(.data$zE_T), is.finite(.data$zE_U), is.finite(.data$zS_T), is.finite(.data$zS_U))
message("  genes with finite statistics in both arms: ", nrow(j))

# ---- Retention per layer, on identical footing ------------------------------------------
boot_ratio <- function(a, b, B = 2000) {
  # ratio of means, percentile CI; deterministic given the data order (no RNG seed needed
  # because we resample indices with a fixed generator state set here)
  set.seed(1)
  n <- length(a)
  r <- vapply(seq_len(B), function(i) {
    ix <- sample.int(n, n, replace = TRUE)
    mean(b[ix]) / mean(a[ix])
  }, numeric(1))
  c(ratio = mean(b) / mean(a), lo = unname(stats::quantile(r, 0.025)),
    hi = unname(stats::quantile(r, 0.975)))
}
pos <- function(x) pmax(x, 0)  # retention is about signal present, not sign
ret_e <- boot_ratio(pos(j$zE_T), pos(j$zE_U))
ret_s <- boot_ratio(pos(j$zS_T), pos(j$zS_U))
ret_raw_e <- boot_ratio(j$oE_T, j$oE_U)
ret_raw_s <- boot_ratio(j$oS_T, j$oS_U)

ju <- inner_join(
  st_up |> select(.data$gene_id, zE_T = .data$zE, zS_T = .data$zS),
  su_up |> select(.data$gene_id, zE_U = .data$zE, zS_U = .data$zS),
  by = "gene_id"
) |> filter(is.finite(.data$zE_T), is.finite(.data$zE_U), is.finite(.data$zS_T), is.finite(.data$zS_U))
ret_e_up <- boot_ratio(pos(ju$zE_T), pos(ju$zE_U))
ret_s_up <- boot_ratio(pos(ju$zS_T), pos(ju$zS_U))

retention <- tibble(
  layer = c("Expression", "Splicing", "Expression (raw, uncorrected)", "Splicing (raw, uncorrected)",
            "Expression (UNPAIRED null)", "Splicing (UNPAIRED null)"),
  statistic = c("excess over paired null", "excess over paired null",
                "mean |log2FC|", "mean TVD",
                "excess over unpaired null", "excess over unpaired null"),
  retained_in_KO = c(ret_e["ratio"], ret_s["ratio"], ret_raw_e["ratio"], ret_raw_s["ratio"],
                     ret_e_up["ratio"], ret_s_up["ratio"]),
  ci_low = c(ret_e["lo"], ret_s["lo"], ret_raw_e["lo"], ret_raw_s["lo"],
             ret_e_up["lo"], ret_s_up["lo"]),
  ci_high = c(ret_e["hi"], ret_s["hi"], ret_raw_e["hi"], ret_raw_s["hi"],
              ret_e_up["hi"], ret_s_up["hi"]),
  n_genes = c(rep(nrow(j), 4), rep(nrow(ju), 2))
)
write_table_pair(retention, results_tables, "layer_retention", cfg = cfg)
print(as.data.frame(retention))

# ---- Coupling between layers -------------------------------------------------------------
# If splicing is purely downstream of expression, the two genotypes should lie on the SAME
# line: the KO simply sits lower along it because its expression response is smaller. A
# genotype difference in SLOPE means the splicing response per unit of expression response
# has changed, which downstream-only cannot produce.
long <- bind_rows(
  j |> transmute(gene_id = .data$gene_id, genotype = "WT (T_UT)", zE = .data$zE_T, zS = .data$zS_T),
  j |> transmute(gene_id = .data$gene_id, genotype = "KO (U_UT)", zE = .data$zE_U, zS = .data$zS_U)
)
fit <- stats::lm(zS ~ zE * genotype, data = long)
co <- summary(fit)$coefficients
inter_row <- grep(":", rownames(co), value = TRUE)[1]
coupling <- tibble(
  slope_WT = unname(co["zE", "Estimate"]),
  slope_difference_KO_minus_WT = unname(co[inter_row, "Estimate"]),
  slope_KO = unname(co["zE", "Estimate"] + co[inter_row, "Estimate"]),
  interaction_p = unname(co[inter_row, "Pr(>|t|)"]),
  n_genes = nrow(j),
  interpretation = paste(
    "A non-zero interaction means the splicing response per unit expression response differs",
    "between genotypes, which a purely downstream model cannot produce."
  )
)
write_table_pair(coupling, results_tables, "layer_coupling", cfg = cfg)
print(as.data.frame(coupling |> select(-.data$interpretation)))

# ---- Matched-response subset -------------------------------------------------------------
# The strictest test: genes whose EXPRESSION response is closely matched between genotypes.
# If splicing is downstream, matching expression should abolish the splicing difference.
# Thresholds are QUANTILE-based: the excess-over-null statistic is in natural units
# (log2 units for expression, TVD units for splicing), so a fixed cutoff is not portable
# between layers or between the paired and unpaired nulls.
qsel <- function(x, y, q, tol_mult = 0.5) {
  tx <- stats::quantile(x, q, na.rm = TRUE); ty <- stats::quantile(y, q, na.rm = TRUE)
  tol <- tol_mult * stats::sd(c(x, y), na.rm = TRUE)
  x > tx & y > ty & abs(x - y) < tol
}
matched <- j |> filter(qsel(.data$zE_T, .data$zE_U, 0.90))
mt <- if (nrow(matched) >= 10L) {
  w <- suppressWarnings(stats::wilcox.test(matched$zS_T, matched$zS_U, paired = TRUE))
  tibble(
    n_matched_genes = nrow(matched),
    median_z_splicing_WT = stats::median(matched$zS_T),
    median_z_splicing_KO = stats::median(matched$zS_U),
    ratio = stats::median(matched$zS_U) / stats::median(matched$zS_T),
    paired_wilcox_p = w$p.value
  )
} else {
  tibble(n_matched_genes = nrow(matched))
}
# RECIPROCAL CONTROL. Matching on one layer and comparing the other will show a deficit in
# BOTH directions whenever the KO is generally lower and the layers are correlated -- that is
# regression to the mean, not evidence of specificity. What discriminates is whether the
# deficit is ASYMMETRIC: splicing should be hit harder when matching on expression than
# expression is when matching on splicing.
recip <- bind_rows(lapply(c(0.80, 0.90, 0.95), function(thr) {
  mE <- j |> filter(qsel(.data$zE_T, .data$zE_U, thr))
  mS <- j |> filter(qsel(.data$zS_T, .data$zS_U, thr))
  bind_rows(
    tibble(threshold = thr, matched_on = "expression", compared = "splicing",
           n = nrow(mE),
           median_WT = stats::median(mE$zS_T), median_KO = stats::median(mE$zS_U),
           ratio_KO_over_WT = stats::median(mE$zS_U) / stats::median(mE$zS_T),
           paired_p = if (nrow(mE) > 10) suppressWarnings(stats::wilcox.test(mE$zS_T, mE$zS_U, paired = TRUE))$p.value else NA_real_),
    tibble(threshold = thr, matched_on = "splicing", compared = "expression",
           n = nrow(mS),
           median_WT = stats::median(mS$zE_T), median_KO = stats::median(mS$zE_U),
           ratio_KO_over_WT = stats::median(mS$zE_U) / stats::median(mS$zE_T),
           paired_p = if (nrow(mS) > 10) suppressWarnings(stats::wilcox.test(mS$zE_T, mS$zE_U, paired = TRUE))$p.value else NA_real_)
  )
})) |>
  mutate(note = "Deficits in both directions are expected from regression to the mean; the ASYMMETRY is the signal.")
write_table_pair(recip, results_tables, "layer_reciprocal_matching", cfg = cfg)
write_table_pair(mt, results_tables, "layer_matched_response", cfg = cfg)
print(as.data.frame(mt))
print(as.data.frame(recip |> select(-.data$note)))

p4 <- ggplot(recip, aes(x = factor(.data$threshold), y = .data$ratio_KO_over_WT,
                        fill = paste0("matched on ", .data$matched_on))) +
  geom_col(position = "dodge", colour = "grey25", linewidth = 0.2) +
  geom_hline(yintercept = 1, linetype = 2, colour = "grey40") +
  scale_fill_manual(values = c("matched on expression" = "#e6550d",
                               "matched on splicing" = "#3182bd"), name = NULL) +
  labs(
    title = "Reciprocal matching: which layer is hit harder?",
    subtitle = "Deficits appear both ways (regression to the mean); the asymmetry is the evidence",
    x = "Matching threshold (quantile of response, in both genotypes)",
    y = "Other layer retained in KO"
  )
ggsave(file.path(fig_dir, "fig_layer_reciprocal.png"), p4, width = 8, height = 5, dpi = 200)

# ---- Baseline composition ----------------------------------------------------------------
# Is there a constitutive splicing defect at rest? Compare unstimulated samples only,
# calibrated against within-genotype replicate variability so the scale is meaningful.
baseline_gap <- function(sl_a, sl_b, name_a, name_b) {
  # Align on isoform_id: the two arms' context tables do not contain identical isoform sets,
  # so gene-level filtering alone leaves vectors of different length.
  prep <- function(sl) {
    m <- sl$iso
    s <- sl$samples[sl$cond == "A"]           # unstimulated only
    g <- m |> group_by(.data$gene_id) |>
      summarize(across(all_of(s), ~ sum(.x, na.rm = TRUE)), .groups = "drop")
    gm <- as.matrix(g[, s, drop = FALSE]); rownames(gm) <- g$gene_id
    d <- gm[match(m$gene_id, g$gene_id), , drop = FALSE]
    f <- ifelse(d > 0, as.matrix(m[, s, drop = FALSE]) / d, 0)
    list(key = m$isoform_id, gene = m$gene_id, f = f,
         gene_expr = stats::setNames(rowMeans(gm) * 2, rownames(gm)))
  }
  A <- prep(sl_a); B <- prep(sl_b)
  shared_iso <- intersect(A$key, B$key)
  ia <- match(shared_iso, A$key); ib <- match(shared_iso, B$key)
  gid <- A$gene[ia]
  keep_g <- intersect(names(A$gene_expr)[A$gene_expr >= expr_floor],
                      names(B$gene_expr)[B$gene_expr >= expr_floor])
  tv <- function(x, y) {
    v <- tapply(abs(x - y), gid, function(z) 0.5 * sum(z, na.rm = TRUE))
    v[intersect(names(v), keep_g)]
  }
  across_g <- tv(rowMeans(A$f[ia, , drop = FALSE]), rowMeans(B$f[ib, , drop = FALSE]))
  pair_mean <- function(M, idx) {
    prs <- utils::combn(ncol(M), 2, simplify = FALSE)
    mean(vapply(prs, function(p) mean(tv(M[idx, p[1]], M[idx, p[2]]), na.rm = TRUE), numeric(1)))
  }
  within_a <- pair_mean(A$f, ia); within_b <- pair_mean(B$f, ib)
  tibble(
    contrast = paste0(name_b, " vs ", name_a, " (unstimulated only)"),
    n_genes = sum(is.finite(across_g)),
    n_shared_isoforms = length(shared_iso),
    mean_tvd_across_genotype = mean(across_g, na.rm = TRUE),
    mean_tvd_within_replicates = mean(c(within_a, within_b)),
    excess_over_replicate_noise = mean(across_g, na.rm = TRUE) - mean(c(within_a, within_b)),
    ratio_across_over_within = mean(across_g, na.rm = TRUE) / mean(c(within_a, within_b))
  )
}
base_rows <- list(baseline_gap(sl_t, sl_u, "WT", "KO"))
sl_th <- sample_level("T_HT"); sl_h <- sample_level("H_HT")
if (!is.null(sl_th) && !is.null(sl_h)) {
  base_rows[[2]] <- baseline_gap(sl_th, sl_h, "T_HT", "H_HT")
}
baseline <- bind_rows(base_rows)
write_table_pair(baseline, results_tables, "layer_baseline_composition", cfg = cfg)
print(as.data.frame(baseline))

# ---- Figures ------------------------------------------------------------------------------
p1 <- ggplot(
  retention |> filter(!grepl("raw|UNPAIRED", .data$layer)),
  aes(x = .data$layer, y = 100 * .data$retained_in_KO)
) +
  geom_col(fill = "#756bb1", colour = "grey25", width = 0.55) +
  geom_errorbar(aes(ymin = 100 * .data$ci_low, ymax = 100 * .data$ci_high), width = 0.15) +
  geom_hline(yintercept = 100, linetype = 2, colour = "grey40") +
  labs(
    title = "How much of the wildtype LPS response does the knockout retain?",
    subtitle = "Excess over a paired permutation null (sign-flips within pairs), identically built for both layers",
    x = NULL, y = "% of wildtype response retained in KO"
  )
ggsave(file.path(fig_dir, "fig_layer_retention.png"), p1, width = 7.5, height = 5, dpi = 200)

p2 <- ggplot(long |> filter(is.finite(.data$zE), is.finite(.data$zS)),
             aes(x = .data$zE, y = .data$zS, colour = .data$genotype)) +
  geom_point(alpha = 0.08, size = 0.5) +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE, linewidth = 0.9) +
  coord_cartesian(xlim = c(-2, 8), ylim = c(-2, 6)) +
  scale_colour_manual(values = c("WT (T_UT)" = "#3182bd", "KO (U_UT)" = "#e6550d"), name = NULL) +
  labs(
    title = "Coupling between the expression and splicing layers",
    subtitle = "Same line = splicing is downstream of expression; different slope = it is not",
    x = "Expression response (standardised)", y = "Splicing response (standardised)"
  )
ggsave(file.path(fig_dir, "fig_layer_coupling.png"), p2, width = 8, height = 5.5, dpi = 200)

p3 <- ggplot(baseline, aes(x = .data$contrast)) +
  geom_col(aes(y = .data$mean_tvd_across_genotype, fill = "across genotypes"),
           width = 0.4, position = position_nudge(x = -0.11), colour = "grey25") +
  geom_col(aes(y = .data$mean_tvd_within_replicates, fill = "within replicates"),
           width = 0.4, position = position_nudge(x = 0.11), colour = "grey25") +
  scale_fill_manual(values = c("across genotypes" = "#e6550d", "within replicates" = "#b3cde3"),
                    name = NULL) +
  labs(
    title = "Baseline isoform composition, unstimulated cells only",
    subtitle = "A constitutive splicing defect would show as across-genotype distance far above replicate noise",
    x = NULL, y = "Mean total variation distance"
  ) +
  theme(axis.text.x = element_text(angle = 12, hjust = 1))
ggsave(file.path(fig_dir, "fig_layer_baseline.png"), p3, width = 8, height = 5, dpi = 200)

# ---- Specification sensitivity -----------------------------------------------------------
# The matched-response result turns out to depend on how the effect is normalised, not on
# pairing. Standardising by the null SD (a z-score) shows a splicing deficit; subtracting the
# null mean (excess) does not. With only 3 paired null patterns the SD cannot be estimated
# reliably, so the z-based version is the less trustworthy of the two -- but the honest
# summary is that this design cannot resolve the question, and reporting either number alone
# would overstate what the data support.
spec <- expand.grid(
  paired = c(TRUE, FALSE), statistic = c("z (divide by null SD)", "excess (subtract null mean)"),
  stringsAsFactors = FALSE
)
spec_rows <- list()
for (i in seq_len(nrow(spec))) {
  a <- layer_stats(sl_t, "T_UT", paired = spec$paired[i])
  b <- layer_stats(sl_u, "U_UT", paired = spec$paired[i])
  jj <- inner_join(
    a |> select(.data$gene_id, oE_T = .data$obs_lfc, nE_T = .data$null_lfc_mean,
                sdE_T = .data$null_lfc_sd, oS_T = .data$obs_tvd,
                nS_T = .data$null_tvd_mean, sdS_T = .data$null_tvd_sd),
    b |> select(.data$gene_id, oE_U = .data$obs_lfc, nE_U = .data$null_lfc_mean,
                sdE_U = .data$null_lfc_sd, oS_U = .data$obs_tvd,
                nS_U = .data$null_tvd_mean, sdS_U = .data$null_tvd_sd),
    by = "gene_id"
  )
  use_z <- grepl("^z ", spec$statistic[i])
  dv <- function(x, sdv) if (use_z) { sdv[!is.finite(sdv) | sdv <= 0] <- NA_real_; x / sdv } else x
  eT <- dv(jj$oE_T - jj$nE_T, jj$sdE_T); eU <- dv(jj$oE_U - jj$nE_U, jj$sdE_U)
  sT <- dv(jj$oS_T - jj$nS_T, jj$sdS_T); sU <- dv(jj$oS_U - jj$nS_U, jj$sdS_U)
  ok <- is.finite(eT) & is.finite(eU) & is.finite(sT) & is.finite(sU)
  k <- ok & qsel(eT, eU, 0.90)
  k[is.na(k)] <- FALSE
  spec_rows[[i]] <- tibble(
    paired = spec$paired[i], statistic = spec$statistic[i],
    n_matched = sum(k),
    splicing_ratio_KO_over_WT = stats::median(sU[k], na.rm = TRUE) / stats::median(sT[k], na.rm = TRUE)
  )
}
spec_tbl <- bind_rows(spec_rows) |>
  mutate(note = "Result depends on normalisation, not on pairing; see script header.")
write_table_pair(spec_tbl, results_tables, "layer_specification_sensitivity", cfg = cfg)
print(as.data.frame(spec_tbl |> select(-.data$note)))

write_run_manifest("08_layer_decoupling.R", cfg, root)
message("08_layer_decoupling.R: done")
