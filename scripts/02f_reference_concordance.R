#!/usr/bin/env Rscript
# 02f: Same-sample concordance across reference transcriptomes -- the technical ceiling.
#
# Run (after 01): Rscript scripts/02f_reference_concordance.R
# Works from committed data/processed/; the external drive is not needed.
#
# WHY
#
# T_HT and T_UT are the SAME six sequencing libraries quantified against two different
# reference transcriptomes. Replicate columns carry identical sample names including the
# S-numbers (T1_minus_S15 ... T3_plus_S20), and gene-level expression correlates at r ~ 0.99
# between the two for matched samples. Biological variability between them is therefore
# zero: every disagreement is annotation plus thresholding.
#
# That makes the pair the best technical control in the project, and it sets a ceiling on
# how much cross-dataset replication can ever be expected. A replication rate should be
# judged against this ceiling, not against 100%.
#
# NOTE ON IDENTIFIERS: TCONS_* isoform ids are assigned per reference and are NOT comparable
# across them -- the same id denotes different transcripts in HT and UT. All comparison here
# is on gene symbols. (Correlating the raw expression matrices by isoform_id gives r ~ 0,
# which looks like different libraries but is purely an id-collision artefact.)

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
fig_dir <- ensure_dir(file.path(results_figures, "reference_concordance"))

min_abs_dif <- as.numeric((cfg$significance %||% list())$min_abs_dif %||% 0.15)
expr_floor <- as.numeric((cfg$significance %||% list())$min_gene_expression %||% 0)
near_miss_lo <- as.numeric(analysis_param(cfg, "near_miss_floor", 0.10))

lab_of <- function(k) sanitize(as.character((cfg$datasets[[k]] %||% list())$label %||% k)[1L])

#' Same-library check, driven off the replicate column names rather than hardcoded.
same_libraries <- function(la, lb) {
  ra <- load_context_table(processed_dir, la, "rep_if")
  rb <- load_context_table(processed_dir, lb, "rep_if")
  if (is.null(ra) || is.null(rb)) return(list(ok = NA, shared = character(), a = character(), b = character()))
  sa <- setdiff(names(ra), "isoform_id")
  sb <- setdiff(names(rb), "isoform_id")
  list(ok = identical(sort(sa), sort(sb)), shared = intersect(sa, sb), a = sa, b = sb)
}

gene_level <- function(k) {
  src <- load_scoring_table(processed_dir, lab_of(k), quiet = TRUE)
  if (is.null(src)) return(NULL)
  z <- score_isoforms(src, cfg, k, lab_of(k))
  z |>
    filter(is_real_gene_symbol(.data$gene_name)) |>
    group_by(.data$gene_name) |>
    summarize(
      switching = any(.data$is_switching, na.rm = TRUE),
      max_abs_dif = max_abs_finite(.data$dIF_n),
      gene_expression = suppressWarnings(max(.data$gene_expression, na.rm = TRUE)),
      above_floor = any(!(.data$low_expression %in% TRUE)),
      .groups = "drop"
    ) |>
    mutate(dataset = k)
}

pair_a <- "T_HT"
pair_b <- "T_UT"
libs <- same_libraries(lab_of(pair_a), lab_of(pair_b))
if (isTRUE(libs$ok)) {
  message(
    "Library check: ", pair_a, " and ", pair_b,
    " share all ", length(libs$shared), " replicate sample names -- same libraries."
  )
} else {
  warning(
    "Library check FAILED: ", pair_a, " and ", pair_b, " do not share replicate sample names (",
    paste(utils::head(libs$a, 3), collapse = ","), " vs ",
    paste(utils::head(libs$b, 3), collapse = ","),
    "). The concordance below is NOT a same-sample technical control -- ",
    "it mixes biological and technical variation. Interpret accordingly."
  )
}

A <- gene_level(pair_a)
B <- gene_level(pair_b)
if (is.null(A) || is.null(B)) stop("Could not score both datasets. Run 01 first.")

m <- inner_join(
  A |> select(.data$gene_name, sw_A = .data$switching, dif_A = .data$max_abs_dif,
              floor_A = .data$above_floor),
  B |> select(.data$gene_name, sw_B = .data$switching, dif_B = .data$max_abs_dif,
              floor_B = .data$above_floor),
  by = "gene_name"
)

n_both <- sum(m$sw_A & m$sw_B)
n_a <- sum(m$sw_A)
n_b <- sum(m$sw_B)
jacc <- if ((n_a + n_b - n_both) > 0) n_both / (n_a + n_b - n_both) else NA_real_
expected <- nrow(m) * (n_a / nrow(m)) * (n_b / nrow(m))
ft <- tryCatch(
  stats::fisher.test(matrix(
    c(sum(!m$sw_A & !m$sw_B), sum(!m$sw_A & m$sw_B),
      sum(m$sw_A & !m$sw_B), n_both), nrow = 2, byrow = TRUE
  )),
  error = function(e) NULL
)
rho <- suppressWarnings(stats::cor(m$dif_A, m$dif_B, method = "spearman", use = "complete.obs"))

# Why do the discordant genes disagree? Classified by their status in the reference that
# did NOT call them, which is what determines whether the disagreement is a threshold
# artefact or a genuine contradiction.
disc <- m |>
  filter(.data$sw_A != .data$sw_B) |>
  mutate(
    other_dif = ifelse(.data$sw_A, .data$dif_B, .data$dif_A),
    other_floor = ifelse(.data$sw_A, .data$floor_B, .data$floor_A),
    reason = dplyr::case_when(
      !(.data$other_floor %in% TRUE) ~ "below the expression floor",
      .data$other_dif >= near_miss_lo & .data$other_dif < min_abs_dif ~
        paste0("near-miss (|dIF| ", near_miss_lo, "-", min_abs_dif, ", above floor)"),
      TRUE ~ paste0("clearly absent (|dIF| < ", near_miss_lo, ", above floor)")
    )
  )

summary_tbl <- tibble(
  dataset_A = pair_a, dataset_B = pair_b,
  same_libraries = libs$ok,
  n_shared_replicate_samples = length(libs$shared),
  n_symbols_testable_in_both = nrow(m),
  n_switching_A = n_a, n_switching_B = n_b, n_switching_both = n_both,
  jaccard = jacc,
  expected_overlap_if_independent = expected,
  odds_ratio = if (!is.null(ft)) unname(ft$estimate) else NA_real_,
  p_value = if (!is.null(ft)) ft$p.value else NA_real_,
  spearman_max_abs_dif = rho,
  n_discordant = nrow(disc),
  median_abs_dif_in_non_calling_reference = suppressWarnings(stats::median(disc$other_dif, na.rm = TRUE)),
  pct_discordant_reaching_near_miss = 100 * mean(disc$other_dif >= near_miss_lo, na.rm = TRUE),
  note = paste(
    "Same libraries, two references: all disagreement is annotation + thresholding.",
    "Judge cross-dataset replication rates against this ceiling, not against 100%."
  )
)
reason_tbl <- disc |> count(.data$reason, name = "n") |> arrange(desc(.data$n))

write_table_pair(summary_tbl, results_tables,
                 paste0("reference_concordance_", pair_a, "_vs_", pair_b), cfg = cfg)
write_table_pair(reason_tbl, results_tables,
                 paste0("reference_discordance_reasons_", pair_a, "_vs_", pair_b), cfg = cfg)
write_table_pair(disc |> arrange(desc(.data$other_dif)), results_tables,
                 paste0("reference_discordant_genes_", pair_a, "_vs_", pair_b), cfg = cfg)
print(as.data.frame(summary_tbl |> select(-.data$note)))
print(as.data.frame(reason_tbl))

p1 <- ggplot(reason_tbl, aes(x = stats::reorder(.data$reason, .data$n), y = .data$n)) +
  geom_col(fill = "#756bb1", colour = "grey25", linewidth = 0.2) +
  geom_text(aes(label = .data$n), hjust = -0.2, size = 3.4) +
  coord_flip() +
  expand_limits(y = max(reason_tbl$n) * 1.15) +
  labs(
    title = paste0("Why the same libraries disagree across references (", pair_a, " vs ", pair_b, ")"),
    subtitle = "Zero biological variability: all of this is annotation plus thresholding",
    x = NULL, y = "Genes"
  )
ggsave(file.path(fig_dir, "fig_reference_discordance_reasons.png"), p1, width = 8.5, height = 4, dpi = 200)

p2 <- ggplot(m |> filter(is.finite(.data$dif_A), is.finite(.data$dif_B)),
             aes(x = .data$dif_A, y = .data$dif_B)) +
  geom_point(alpha = 0.15, size = 0.6) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey40") +
  geom_hline(yintercept = min_abs_dif, linetype = 3, colour = "#e6550d") +
  geom_vline(xintercept = min_abs_dif, linetype = 3, colour = "#e6550d") +
  coord_cartesian(xlim = c(0, 0.6), ylim = c(0, 0.6)) +
  labs(
    title = "Gene-level max |dIF|, same libraries under two references",
    subtitle = paste0("Spearman ", sprintf("%.3f", rho),
                      "; dotted lines mark the ", min_abs_dif, " calling threshold"),
    x = paste0("max |dIF| in ", pair_a), y = paste0("max |dIF| in ", pair_b)
  )
ggsave(file.path(fig_dir, "fig_reference_dif_scatter.png"), p2, width = 6.5, height = 6, dpi = 200)

write_run_manifest("02f_reference_concordance.R", cfg, root)
message("02f_reference_concordance.R: done")
