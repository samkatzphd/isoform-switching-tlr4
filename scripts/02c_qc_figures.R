#!/usr/bin/env Rscript
# 02c: QC and overview figures (gene- and isoform-level)
# Run: Rscript scripts/02c_qc_figures.R  (from project root; after 01, 02)

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("Install ggplot2: install.packages('ggplot2')")
}
has_scales <- requireNamespace("scales", quietly = TRUE)

bt <- c("utils/bootstrap.R", file.path("..", "utils", "bootstrap.R"))
b_file <- if (any(f <- vapply(bt, file.exists, NA))) { bt[which(f)[1L]] } else { NA_character_ }
if (is.na(b_file)) {
  stop("Run from project root, or set ISOFORM_PROJECT_ROOT.")
}
source(b_file, local = FALSE, chdir = FALSE)
root <- if (exists("PROJECT_ROOT", inherits = FALSE)) {
  get("PROJECT_ROOT", inherits = FALSE)
} else {
  find_project_root()
}

cfg <- load_yaml_config("config/config.yml", root = root)
paths <- cfg$paths %||% list()
out_dir <- ensure_dir(resolve_path(paths$results_figures %||% "results/figures", root = root))
tab_dir <- resolve_path(paths$results_tables %||% "results/tables", root = root)
proc_dir <- resolve_path(paths$processed_dir %||% "data/processed", root = root)

ggplot2::theme_set(
  ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 12),
      strip.text = ggplot2::element_text(face = "bold")
    )
)

# ---- Gene-level: novel involvement, max |dIF|, n switching isoforms ----
gs <- list.files(
  tab_dir, pattern = "^gene_level_summary_.*\\.csv$", full.names = TRUE
)
if (length(gs) < 1L) {
  stop("No gene_level_summary_*.csv in ", tab_dir, ". Run 02 first.")
}
gene_long <- lapply(gs, function(f) {
  d <- utils::read.csv(f, check.names = FALSE, stringsAsFactors = FALSE)
  d$._dataset <- sub("^gene_level_summary_(.+)\\.csv$", "\\1", basename(f), perl = TRUE)
  d
})
gdf <- do.call(rbind, gene_long)
gdf$novel_flag <- gdf$novel_involved %in% TRUE
gdf$n_switching_isoforms <- as.numeric(gdf$n_switching_isoforms)
gdf$max_abs_dif <- as.numeric(gdf$max_abs_dif)
if (!"min_isoform_switch_q" %in% names(gdf)) {
  gdf$min_isoform_switch_q <- NA_real_
} else {
  gdf$min_isoform_switch_q <- as.numeric(gdf$min_isoform_switch_q)
}

novel_tally <- gdf |>
  dplyr::group_by(.data$._dataset, .data$novel_flag) |>
  dplyr::summarise(n = dplyr::n(), .groups = "drop")

p_novel <- ggplot2::ggplot(
  novel_tally,
  ggplot2::aes(x = .data$._dataset, y = .data$n, fill = .data$novel_flag)
) +
  ggplot2::geom_col(position = "fill", color = "grey30", linewidth = 0.2) +
  ggplot2::scale_y_continuous(
    labels = if (has_scales) scales::percent else function(x) paste0(round(100 * x, 1), "%")
  ) +
  ggplot2::scale_fill_manual(
    values = c("FALSE" = "#b3cde3", "TRUE" = "#8856a7"),
    labels = c("FALSE" = "Novel not involved", "TRUE" = "Novel involved"),
    name = NULL
  ) +
  ggplot2::labs(
    title = "Gene-level: fraction with PacBio novel isoform involvement",
    x = NULL, y = "Fraction of genes"
  ) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 35, hjust = 1))

ggplot2::ggsave(
  file.path(out_dir, "fig_qc_novel_fraction_stacked.png"),
  p_novel, width = 7.2, height = 4.2, dpi = 200
)

p_violin_dif <- ggplot2::ggplot(
  gdf,
  ggplot2::aes(x = .data$._dataset, y = .data$max_abs_dif)
) +
  ggplot2::geom_violin(trim = TRUE, alpha = 0.5, fill = "grey80") +
  ggplot2::geom_boxplot(width = 0.15, outlier.alpha = 0.3) +
  ggplot2::labs(
    title = "Distribution of max |dIF| per gene",
    x = NULL, y = "max |dIF| (gene level)"
  ) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 35, hjust = 1))

ggplot2::ggsave(
  file.path(out_dir, "fig_qc_gene_max_abs_dif_violin.png"),
  p_violin_dif, width = 7.2, height = 4.5, dpi = 200
)

p_nsw <- ggplot2::ggplot(
  gdf,
  ggplot2::aes(x = .data$._dataset, y = .data$n_switching_isoforms)
) +
  ggplot2::geom_violin(trim = TRUE, alpha = 0.5, fill = "#e5f5f9") +
  ggplot2::geom_boxplot(width = 0.15, outlier.alpha = 0.3) +
  ggplot2::labs(
    title = "Switching isoforms per gene (counts meeting pipeline thresholds)",
    x = NULL, y = "n switching isoforms"
  ) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 35, hjust = 1))

ggplot2::ggsave(
  file.path(out_dir, "fig_qc_gene_n_switching_isoforms_violin.png"),
  p_nsw, width = 7.2, height = 4.5, dpi = 200
)

# ---- Isoform-level: |dIF| density + volcano (facets) ----
iso_files <- list.files(
  proc_dir, pattern = "^isoformFeatures_.*\\.rds$", full.names = TRUE
)
if (length(iso_files) < 1L) {
  message("No isoformFeatures_*.rds — skip isoform figures.")
} else {
  idf <- lapply(iso_files, function(f) {
    d <- readRDS(f)
    id <- sub("^isoformFeatures_(.+)\\.rds$", "\\1", basename(f), perl = TRUE, ignore.case = TRUE)
    if (!"dIF" %in% names(d) || !"isoform_switch_q_value" %in% names(d)) {
      return(NULL)
    }
    d2 <- tibble::as_tibble(d)
    d2$.dataset <- id
    d2$dif <- as.numeric(d2$dIF)
    d2$nlp <- -log10(pmax(as.numeric(d2$isoform_switch_q_value), 1e-300, na.rm = FALSE))
    if ("is_novel_pacbio" %in% names(d2)) {
      d2$is_novel <- d2$is_novel_pacbio %in% TRUE
    } else {
      d2$is_novel <- NA
    }
    d2[, c("gene_id", "isoform_id", "dif", "nlp", "is_novel", ".dataset")]
  })
  idf <- idf[!vapply(idf, is.null, NA)]
  if (length(idf) < 1L) {
    message("Isoform tables missing dIF or q — skip isoform figures.")
  } else {
    idf <- do.call(dplyr::bind_rows, idf)
    idf$novel_lab <- dplyr::case_when(
      is.na(idf$is_novel) ~ "Unknown",
      idf$is_novel ~ "Novel (PB)",
      !idf$is_novel ~ "Other"
    )
    p_hist <- ggplot2::ggplot(
      idf, ggplot2::aes(x = abs(.data$dif))
    ) +
      ggplot2::geom_histogram(bins = 50, fill = "steelblue", color = "white", linewidth = 0.1) +
      ggplot2::facet_wrap(~ .dataset, scales = "free_y", ncol = 2) +
      ggplot2::labs(
        title = "Isoform-level: distribution of |dIF|",
        x = "|dIF|", y = "Count"
      )

    ggplot2::ggsave(
      file.path(out_dir, "fig_qc_isoform_abs_dif_histogram.png"),
      p_hist, width = 7.5, height = 6, dpi = 200
    )

    p_vol <- ggplot2::ggplot(
      idf, ggplot2::aes(x = .data$dif, y = .data$nlp, color = .data$novel_lab)
    ) +
      ggplot2::geom_point(alpha = 0.2, size = 0.6) +
      ggplot2::scale_color_manual(
        values = c("Other" = "grey50", "Novel (PB)" = "darkorange", "Unknown" = "grey45"),
        name = NULL
      ) +
      ggplot2::facet_wrap(~ .dataset, scales = "free", ncol = 2) +
      ggplot2::labs(
        title = "Isoform-level volcano: dIF vs -log10(isoform q)",
        x = "dIF", y = expression(-log[10](q))
      ) +
      ggplot2::guides(
        color = ggplot2::guide_legend(override.aes = list(alpha = 1, size = 2))
      )

    ggplot2::ggsave(
      file.path(out_dir, "fig_qc_isoform_volcano_by_dataset.png"),
      p_vol, width = 7.5, height = 6, dpi = 200
    )
  }
}

message("02c_qc_figures.R: done — figures in: ", out_dir)
