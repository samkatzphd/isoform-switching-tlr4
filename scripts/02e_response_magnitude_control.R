#!/usr/bin/env Rscript
# 02e: Is reduced isoform switching in one arm explained by a smaller LPS response overall?
#
# Run (after 01): Rscript scripts/02e_response_magnitude_control.R
# Works from committed data/processed/; the external drive is not needed.
#
# WHY
#
# 02d rules out a *baseline abundance* confounder -- the two UT arms have near-identical
# expression distributions. That is the wrong confounder for a post-LPS switching claim.
# The relevant one is the size of the LPS response itself: isoform switching needs a
# condition-driven shift in isoform fraction to clear a fixed |dIF| threshold, so if the
# whole perturbation runs at reduced amplitude fewer shifts clear it, for reasons that have
# nothing to do with splicing.
#
# This script measures response amplitude per gene, stratifies on it, and asks whether the
# switching deficit survives. It does so for BOTH pairs, so the arms are treated alike.
#
# INTERPRETATION -- read before quoting the adjusted odds ratio.
#
# Response magnitude is plausibly a MEDIATOR, not a confounder: if losing UBL5 blunts the
# LPS response, and a blunted response yields fewer switches, then adjusting for response
# magnitude removes part of the very effect being measured. A non-significant adjusted OR
# therefore does NOT establish that there is no splicing-specific role -- it establishes
# that this design cannot separate the two. The parsimonious reading is that the global
# blunting accounts for the deficit; a direct splicing effect on top of it remains possible
# and would need a different experiment to detect.

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
fig_dir <- ensure_dir(file.path(results_figures, "response"))

expr_floor <- as.numeric((cfg$significance %||% list())$min_gene_expression %||% 0)
strata <- as.numeric(analysis_param(cfg, "response_strata", c(0.1, 0.25, 0.5, 1)))
brk <- c(-Inf, sort(unique(strata)), Inf)
lab <- c(
  paste0("<", brk[2]),
  vapply(seq_len(length(brk) - 3L), function(i) paste0(brk[i + 1], "-", brk[i + 2]), character(1)),
  paste0(">", brk[length(brk) - 1])
)

#' Per-gene LPS response amplitude and switching status.
#'
#' gene_value_1 / gene_value_2 are the condition means and are constant within gene_id;
#' that is asserted rather than assumed. Grouping is by gene_id ALONE -- grouping by
#' (gene_id, gene_name) splits genes whose isoforms disagree on the symbol and produces a
#' many-to-many join downstream, which silently inflates every count.
gene_response <- function(label) {
  src <- load_scoring_table(processed_dir, label, quiet = TRUE)
  if (is.null(src)) return(NULL)
  if (!all(c("gene_value_1", "gene_value_2") %in% names(src))) {
    warning("[", label, "] no gene_value_1/2 columns; skipping.")
    return(NULL)
  }
  chk <- src |>
    group_by(.data$gene_id) |>
    summarize(k = dplyr::n_distinct(.data$gene_value_1) + dplyr::n_distinct(.data$gene_value_2),
              .groups = "drop")
  if (any(chk$k != 2L)) {
    stop("[", label, "] gene_value_1/2 are not constant within gene_id; the log2FC below would be ill-defined.")
  }
  z <- score_isoforms(src, cfg, label, label)
  z |>
    group_by(.data$gene_id) |>
    summarize(
      v1 = dplyr::first(.data$gene_value_1),
      v2 = dplyr::first(.data$gene_value_2),
      switching = any(.data$is_switching, na.rm = TRUE),
      .groups = "drop"
    ) |>
    mutate(
      gene_expression = .data$v1 + .data$v2,
      lfc = log2((.data$v2 + 1) / (.data$v1 + 1)),
      abs_lfc = abs(.data$lfc)
    ) |>
    filter(.data$gene_expression >= expr_floor) |>
    mutate(stratum = cut(.data$abs_lfc, breaks = brk, labels = lab), dataset = label)
}

analyse_pair <- function(key_ref, key_alt, pair_name) {
  a <- gene_response(sanitize(as.character((cfg$datasets[[key_ref]] %||% list())$label %||% key_ref)[1L]))
  b <- gene_response(sanitize(as.character((cfg$datasets[[key_alt]] %||% list())$label %||% key_alt)[1L]))
  if (is.null(a) || is.null(b)) return(NULL)

  both <- inner_join(
    a |> select(.data$gene_id, lfc_ref = .data$lfc),
    b |> select(.data$gene_id, lfc_alt = .data$lfc),
    by = "gene_id"
  )
  fit <- stats::lm(lfc_alt ~ lfc_ref, data = both)
  top <- both |> arrange(desc(abs(.data$lfc_ref))) |> utils::head(200)

  amp <- tibble(
    pair = pair_name, reference_arm = key_ref, comparison_arm = key_alt,
    n_genes_tested_in_both = nrow(both),
    sd_lfc_ref = stats::sd(both$lfc_ref), sd_lfc_alt = stats::sd(both$lfc_alt),
    sd_ratio = stats::sd(both$lfc_alt) / stats::sd(both$lfc_ref),
    iqr_lfc_ref = stats::IQR(both$lfc_ref), iqr_lfc_alt = stats::IQR(both$lfc_alt),
    iqr_ratio = stats::IQR(both$lfc_alt) / stats::IQR(both$lfc_ref),
    mean_abs_lfc_ref = mean(abs(both$lfc_ref)), mean_abs_lfc_alt = mean(abs(both$lfc_alt)),
    mean_abs_ratio = mean(abs(both$lfc_alt)) / mean(abs(both$lfc_ref)),
    regression_slope_alt_on_ref = unname(stats::coef(fit)[2]),
    spearman_ref_vs_alt = stats::cor(both$lfc_ref, both$lfc_alt, method = "spearman"),
    pct_magnitude_retained_top200 = 100 * mean(abs(top$lfc_alt)) / mean(abs(top$lfc_ref))
  )

  st <- full_join(
    a |> group_by(.data$stratum) |>
      summarize(n_ref = dplyr::n(), sw_ref = sum(.data$switching, na.rm = TRUE), .groups = "drop"),
    b |> group_by(.data$stratum) |>
      summarize(n_alt = dplyr::n(), sw_alt = sum(.data$switching, na.rm = TRUE), .groups = "drop"),
    by = "stratum"
  ) |>
    mutate(
      pair = pair_name,
      pct_ref = 100 * .data$sw_ref / .data$n_ref,
      pct_alt = 100 * .data$sw_alt / .data$n_alt
    ) |>
    filter(!is.na(.data$stratum))

  crude <- tryCatch(
    stats::fisher.test(matrix(
      c(sum(b$switching), nrow(b) - sum(b$switching),
        sum(a$switching), nrow(a) - sum(a$switching)),
      nrow = 2, byrow = TRUE
    )),
    error = function(e) NULL
  )
  arr <- array(0L, c(2, 2, nrow(st)))
  for (i in seq_len(nrow(st))) {
    arr[, , i] <- matrix(
      c(st$sw_alt[i], st$n_alt[i] - st$sw_alt[i], st$sw_ref[i], st$n_ref[i] - st$sw_ref[i]),
      nrow = 2, byrow = TRUE
    )
  }
  mh <- tryCatch(stats::mantelhaen.test(arr), error = function(e) NULL)

  ctrl <- tibble(
    pair = pair_name, reference_arm = key_ref, comparison_arm = key_alt,
    n_ref = nrow(a), n_alt = nrow(b),
    n_switching_ref = sum(a$switching), n_switching_alt = sum(b$switching),
    crude_odds_ratio = if (!is.null(crude)) unname(crude$estimate) else NA_real_,
    crude_conf_low = if (!is.null(crude)) crude$conf.int[[1L]] else NA_real_,
    crude_conf_high = if (!is.null(crude)) crude$conf.int[[2L]] else NA_real_,
    crude_p = if (!is.null(crude)) crude$p.value else NA_real_,
    mh_odds_ratio = if (!is.null(mh)) unname(mh$estimate) else NA_real_,
    mh_conf_low = if (!is.null(mh)) mh$conf.int[[1L]] else NA_real_,
    mh_conf_high = if (!is.null(mh)) mh$conf.int[[2L]] else NA_real_,
    mh_p = if (!is.null(mh)) mh$p.value else NA_real_,
    note = paste(
      "Response magnitude may be a mediator rather than a confounder;",
      "a non-significant adjusted OR shows the design cannot separate a splicing-specific",
      "effect from a globally blunted response, not that none exists."
    )
  ) |>
    bind_cols(amp |> select(-.data$pair, -.data$reference_arm, -.data$comparison_arm))

  message(sprintf(
    "[%s] amplitude ratio (slope) %.3f | crude OR %.3f | MH OR %.3f (p = %.3f)",
    pair_name, amp$regression_slope_alt_on_ref,
    ctrl$crude_odds_ratio, ctrl$mh_odds_ratio, ctrl$mh_p
  ))
  list(strata = st, control = ctrl, genes = bind_rows(a, b))
}

res <- list(
  analyse_pair("T_UT", "U_UT", "UT: WT vs UBL5 KO"),
  analyse_pair("T_HT", "H_HT", "HT: T vs H")
)
res <- res[!vapply(res, is.null, NA)]
if (!length(res)) stop("No pairs could be analysed. Run 01 first.")

strata_tbl <- bind_rows(lapply(res, `[[`, "strata"))
control_tbl <- bind_rows(lapply(res, `[[`, "control"))
write_table_pair(strata_tbl, results_tables, "ut_response_magnitude_strata", cfg = cfg)
write_table_pair(control_tbl, results_tables, "ut_response_magnitude_control", cfg = cfg)
print(as.data.frame(strata_tbl))
print(as.data.frame(
  control_tbl |> select(.data$pair, .data$crude_odds_ratio, .data$mh_odds_ratio,
                        .data$mh_conf_low, .data$mh_conf_high, .data$mh_p,
                        .data$regression_slope_alt_on_ref)
))

genes <- bind_rows(lapply(res, `[[`, "genes"))

p1 <- ggplot(
  strata_tbl |>
    tidyr::pivot_longer(c("pct_ref", "pct_alt"), names_to = "arm", values_to = "pct") |>
    mutate(arm = ifelse(.data$arm == "pct_ref", "reference arm", "comparison arm")),
  aes(x = .data$stratum, y = .data$pct, fill = .data$arm)
) +
  geom_col(position = "dodge", colour = "grey25", linewidth = 0.2) +
  facet_wrap(~ .data$pair, scales = "free_y") +
  scale_fill_manual(values = c("reference arm" = "#3182bd", "comparison arm" = "#e6550d"), name = NULL) +
  labs(
    title = "Switching rate within strata of LPS response magnitude",
    subtitle = "Each gene binned by its own |log2FC|, measured in its own dataset",
    x = "|log2 fold-change| stratum", y = "% of genes with a switch"
  )
ggsave(file.path(fig_dir, "fig_response_strata.png"), p1, width = 9.5, height = 4.5, dpi = 200)

p2 <- ggplot(genes, aes(x = .data$abs_lfc, colour = .data$dataset)) +
  stat_ecdf(linewidth = 0.8) +
  coord_cartesian(xlim = c(0, 3)) +
  labs(
    title = "Magnitude of the LPS response by dataset",
    subtitle = "A shifted curve means the whole perturbation runs at different amplitude",
    x = "|log2 fold-change| (LPS+ vs LPS-)", y = "ECDF"
  )
ggsave(file.path(fig_dir, "fig_response_amplitude.png"), p2, width = 8, height = 5, dpi = 200)

write_run_manifest("02e_response_magnitude_control.R", cfg, root)
message("02e_response_magnitude_control.R: done")
