#!/usr/bin/env Rscript
# 06: HT-focused ranking + IsoformSwitchAnalyzeR switch plots
#     - Targets datasets T_HT and H_HT
#     - Ranks genes with significant isoform-fraction changes (q-cutoff)
#     - Writes top-gene tables
#     - Draws transcript-structure switch plots for top genes
#     - Builds T_HT vs H_HT isoform overlap tables (including non-significant H)
#
# Run: Rscript scripts/06_visualization.R

bt <- c("utils/bootstrap.R", file.path("..", "utils", "bootstrap.R"))
b_file <- if (any(f <- vapply(bt, file.exists, NA))) { bt[which(f)[1L]] } else { NA_character_ }
if (is.na(b_file)) {
  stop("Could not find utils/bootstrap.R. See README (ISOFORM_PROJECT_ROOT or `cd` into project).")
}
source(b_file, local = FALSE, chdir = FALSE)

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(tibble)
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
out_fig_dir <- ensure_dir(file.path(results_figures, "ht_top_switch_plots"))

sig <- cfg$significance %||% list()
iso_q_cutoff <- as.numeric(sig$isoform_q %||% 0.05)
gene_q_cutoff <- as.numeric(sig$gene_q %||% 0.05)
min_abs_dif <- as.numeric(sig$min_abs_dif %||% 0.0)

target_datasets <- c("T_HT", "H_HT")
top_n_genes <- 30L
plot_top_n <- 30L

sanitize <- function(x) {
  gsub("[^A-Za-z0-9_]+", "_", x, perl = TRUE)
}
min_finite <- function(v) {
  if (!length(v)) return(NA_real_)
  m <- suppressWarnings(min(v, na.rm = TRUE))
  if (is.infinite(m)) NA_real_ else m
}

compute_switching_gene_rank <- function(iso_tbl, dataset_key, dataset_label) {
  z <- as_tibble(iso_tbl)
  if (!"gene_id" %in% names(z) || !"isoform_id" %in% names(z)) {
    stop("Isoform table must include gene_id and isoform_id.")
  }
  if (!"gene_name" %in% names(z)) z$gene_name <- NA_character_
  if (!"dIF" %in% names(z)) z$dIF <- NA_real_
  if (!"isoform_switch_q_value" %in% names(z) && !"gene_switch_q_value" %in% names(z)) {
    stop("Need isoform_switch_q_value or gene_switch_q_value.")
  }
  z$q_i <- if ("isoform_switch_q_value" %in% names(z)) as.numeric(z$isoform_switch_q_value) else as.numeric(z$gene_switch_q_value)
  z$q_g <- if ("gene_switch_q_value" %in% names(z)) as.numeric(z$gene_switch_q_value) else NA_real_
  z$dIF_n <- as.numeric(z$dIF)
  z$abs_dIF <- abs(z$dIF_n)
  z$is_switching <- is.finite(z$q_i) & z$q_i < iso_q_cutoff &
    (is.na(z$q_g) | (is.finite(z$q_g) & z$q_g < gene_q_cutoff)) &
    (is.na(z$dIF_n) | abs(z$dIF_n) >= min_abs_dif)
  z$is_switching <- replace(z$is_switching, is.na(z$is_switching), FALSE)

  ranks <- z |>
    group_by(.data$gene_id, .data$gene_name) |>
    summarize(
      n_isoforms = n(),
      n_switching_isoforms = sum(.data$is_switching, na.rm = TRUE),
      n_switching_up = sum(.data$is_switching & is.finite(.data$dIF_n) & .data$dIF_n > 0, na.rm = TRUE),
      n_switching_down = sum(.data$is_switching & is.finite(.data$dIF_n) & .data$dIF_n < 0, na.rm = TRUE),
      max_abs_dif_switching = {
        v <- .data$abs_dIF[.data$is_switching]
        if (!length(v) || all(!is.finite(v))) NA_real_ else max(v, na.rm = TRUE)
      },
      min_isoform_switch_q = min_finite(.data$q_i[.data$is_switching]),
      min_gene_switch_q = min_finite(.data$q_g[.data$is_switching]),
      .groups = "drop"
    ) |>
    mutate(
      has_opposing_direction_switches = (.data$n_switching_up > 0L) & (.data$n_switching_down > 0L),
      has_multi_isoform_opposing_switch = (.data$n_switching_isoforms >= 2L) &
        (.data$n_switching_up > 0L) & (.data$n_switching_down > 0L)
    ) |>
    filter(.data$n_switching_isoforms > 0L) |>
    arrange(
      desc(.data$has_multi_isoform_opposing_switch),
      desc(.data$max_abs_dif_switching),
      desc(.data$n_switching_isoforms),
      .data$min_isoform_switch_q
    ) |>
    mutate(
      rank = row_number(),
      dataset_key = dataset_key,
      dataset_label = dataset_label,
      isoform_q_cutoff = iso_q_cutoff,
      gene_q_cutoff = gene_q_cutoff,
      min_abs_dif_cutoff = min_abs_dif
    )
  list(iso = z, ranks = ranks)
}

plot_gene_fallback <- function(gene_row, iso_tbl, dataset_label) {
  gid <- as.character(gene_row$gene_id[[1]])
  gname <- as.character(gene_row$gene_name[[1]])
  if (is.na(gname) || !nzchar(gname)) gname <- gid

  iso_tbl |>
    filter(.data$gene_id == gid) |>
    arrange(desc(abs(.data$dIF_n)), .data$q_i) |>
    slice_head(n = 12L) |>
    mutate(isoform_id = as.character(.data$isoform_id)) |>
    ggplot(aes(x = reorder(.data$isoform_id, .data$dIF_n), y = .data$dIF_n, fill = .data$is_switching)) +
    geom_col() +
    coord_flip() +
    scale_fill_manual(values = c(`TRUE` = "#1b7837", `FALSE` = "#bdbdbd")) +
    labs(
      title = paste0(dataset_label, " | ", gname, " (", gid, ")"),
      subtitle = "Fallback dIF view (switchPlot failed for this gene)",
      x = "Isoform",
      y = "dIF",
      fill = "Significant"
    ) +
    theme_bw(base_size = 10)
}

make_switch_plot <- function(isa_obj, gene_id, gene_name, out_png) {
  if (!requireNamespace("IsoformSwitchAnalyzeR", quietly = TRUE)) {
    return(list(ok = FALSE, message = "IsoformSwitchAnalyzeR is not installed."))
  }
  switch_fn <- get("switchPlot", envir = asNamespace("IsoformSwitchAnalyzeR"))
  candidates <- unique(na.omit(c(as.character(gene_id), as.character(gene_name))))
  if (!length(candidates)) {
    return(list(ok = FALSE, message = "No gene identifier available."))
  }

  # switchPlot draws directly, so render inside png device.
  for (g in candidates) {
    ok <- tryCatch({
      png(filename = out_png, width = 2400, height = 1500, res = 180)
      on.exit(dev.off(), add = TRUE)
      # Try the most common function signatures across package versions.
      tryCatch(
        do.call(switch_fn, list(switchAnalyzeRlist = isa_obj, gene = g)),
        error = function(e1) do.call(switch_fn, list(isa_obj, gene = g))
      )
      TRUE
    }, error = function(e) {
      FALSE
    })
    if (ok && file.exists(out_png) && file.info(out_png)$size > 0) {
      return(list(ok = TRUE, message = paste0("switchPlot succeeded for ", g)))
    }
  }
  list(ok = FALSE, message = "switchPlot failed for both gene_id and gene_name.")
}

load_isa_for_dataset <- function(ds) {
  path_raw <- ds$isa_path %||% NULL
  if (is.null(path_raw) || !nzchar(as.character(path_raw)[1L])) return(NULL)
  abs_in <- resolve_path(as.character(path_raw)[1L], root = root)
  if (!file.exists(abs_in)) return(NULL)
  obj_name <- ds$object_name %||% NULL
  load_isa_input(abs_in, object_name = obj_name)
}

extract_t_h_overlap <- function(t_iso, h_iso, iso_q_cutoff, gene_q_cutoff, min_abs_dif) {
  select_cols <- c("isoform_id", "gene_id", "gene_name", "dIF_n", "q_i", "q_g", "is_switching")
  t_keep <- t_iso[, intersect(select_cols, names(t_iso)), drop = FALSE] |>
    mutate(
      isoform_id = as.character(.data$isoform_id),
      gene_id = as.character(.data$gene_id),
      gene_name = as.character(.data$gene_name)
    ) |>
    rename_with(~ paste0("T_", .x), -c("isoform_id", "gene_id", "gene_name"))
  h_keep <- h_iso[, intersect(select_cols, names(h_iso)), drop = FALSE] |>
    mutate(
      isoform_id = as.character(.data$isoform_id),
      gene_id = as.character(.data$gene_id),
      gene_name = as.character(.data$gene_name)
    ) |>
    rename_with(~ paste0("H_", .x), -c("isoform_id", "gene_id", "gene_name"))

  merged <- full_join(t_keep, h_keep, by = c("isoform_id", "gene_id", "gene_name")) |>
    mutate(
      in_T = !is.na(.data$T_dIF_n) | !is.na(.data$T_q_i),
      in_H = !is.na(.data$H_dIF_n) | !is.na(.data$H_q_i),
      T_significant = .data$T_is_switching %in% TRUE,
      H_significant = .data$H_is_switching %in% TRUE,
      present_in_both = .data$in_T & .data$in_H,
      same_direction = dplyr::case_when(
        is.finite(.data$T_dIF_n) & is.finite(.data$H_dIF_n) ~ sign(.data$T_dIF_n) == sign(.data$H_dIF_n),
        TRUE ~ NA
      ),
      abs_dIF_delta = dplyr::case_when(
        is.finite(.data$T_dIF_n) & is.finite(.data$H_dIF_n) ~ abs(.data$T_dIF_n - .data$H_dIF_n),
        TRUE ~ NA_real_
      )
    )

  t_sig_in_h <- merged |>
    filter(.data$T_significant, .data$present_in_both) |>
    mutate(
      H_passes_significance = .data$H_significant,
      overlap_note = case_when(
        .data$H_significant ~ "Present in H and significant in H",
        .data$in_H ~ "Present in H but not significant in H",
        TRUE ~ "Not found in H"
      )
    ) |>
    arrange(desc(abs(.data$T_dIF_n)), .data$T_q_i)

  summary <- tibble(
    isoform_q_cutoff = iso_q_cutoff,
    gene_q_cutoff = gene_q_cutoff,
    min_abs_dif_cutoff = min_abs_dif,
    n_isoforms_T = sum(!is.na(t_keep$isoform_id)),
    n_isoforms_H = sum(!is.na(h_keep$isoform_id)),
    n_isoforms_shared = sum(merged$present_in_both, na.rm = TRUE),
    n_T_significant = sum(merged$T_significant, na.rm = TRUE),
    n_H_significant = sum(merged$H_significant, na.rm = TRUE),
    n_T_significant_present_in_H = sum(merged$T_significant & merged$present_in_both, na.rm = TRUE),
    n_T_significant_and_H_significant = sum(merged$T_significant & merged$H_significant, na.rm = TRUE),
    n_T_significant_present_in_H_same_direction = sum(
      merged$T_significant & merged$present_in_both & (merged$same_direction %in% TRUE),
      na.rm = TRUE
    )
  )
  list(all = merged, t_sig_in_h = t_sig_in_h, summary = summary)
}

all_rankings <- list()
all_iso <- list()

for (ds_key in target_datasets) {
  ds <- cfg$datasets[[ds_key]] %||% NULL
  if (is.null(ds)) {
    warning("Dataset not found in config: ", ds_key, " (skipping)")
    next
  }
  ds_label <- as.character(ds$label %||% ds_key)[1L]
  label_clean <- sanitize(ds_label)
  iso_rds <- file.path(processed_dir, paste0("isoformFeatures_", label_clean, ".rds"))
  if (!file.exists(iso_rds)) {
    warning("Processed isoform table missing for ", ds_key, ": ", iso_rds, " (run scripts/01_load_data.R)")
    next
  }

  message("[", ds_key, "] Reading ", iso_rds)
  iso_tbl <- readRDS(iso_rds)
  out <- compute_switching_gene_rank(iso_tbl, ds_key, ds_label)
  ranks <- out$ranks
  iso_scored <- out$iso
  all_rankings[[ds_key]] <- ranks
  all_iso[[ds_key]] <- iso_scored

  top_tbl <- ranks |>
    slice_head(n = min(top_n_genes, nrow(ranks)))
  out_csv <- file.path(results_tables, paste0("top_switching_genes_", label_clean, ".csv"))
  out_rds <- file.path(results_tables, paste0("top_switching_genes_", label_clean, ".rds"))
  utils::write.csv(top_tbl, out_csv, row.names = FALSE, fileEncoding = "UTF-8", na = "")
  saveRDS(top_tbl, out_rds, compress = "xz")
  message("[", ds_key, "] Wrote top genes: ", out_csv)

  isa_obj <- load_isa_for_dataset(ds)
  if (is.null(isa_obj)) {
    warning("[", ds_key, "] Could not load ISA object; switchPlot outputs will be skipped.")
  }

  gene_rows <- top_tbl |>
    slice_head(n = min(plot_top_n, nrow(top_tbl)))
  if (!nrow(gene_rows)) {
    message("[", ds_key, "] No significant switching genes to plot.")
    next
  }
  for (i in seq_len(nrow(gene_rows))) {
    gr <- gene_rows[i, , drop = FALSE]
    gene_stub <- sanitize(as.character(gr$gene_name[[1]] %||% gr$gene_id[[1]]))
    out_png <- file.path(out_fig_dir, paste0("switch_plot_", label_clean, "_rank", sprintf("%02d", i), "_", gene_stub, ".png"))

    made <- FALSE
    if (!is.null(isa_obj)) {
      sp <- make_switch_plot(
        isa_obj = isa_obj,
        gene_id = as.character(gr$gene_id[[1]]),
        gene_name = as.character(gr$gene_name[[1]]),
        out_png = out_png
      )
      made <- isTRUE(sp$ok)
      if (!made) {
        message("[", ds_key, "] switchPlot unavailable for ", gr$gene_id[[1]], "; writing fallback plot.")
      }
    }
    if (!made) {
      p <- plot_gene_fallback(gr, iso_scored, ds_label)
      ggsave(out_png, p, width = 10, height = 5, units = "in", dpi = 150)
    }
  }
}

if (length(all_rankings)) {
  combined <- bind_rows(all_rankings) |>
    arrange(.data$dataset_key, .data$rank)
  utils::write.csv(
    combined,
    file.path(results_tables, "top_switching_genes_HT_combined.csv"),
    row.names = FALSE,
    fileEncoding = "UTF-8",
    na = ""
  )
}

if (all(c("T_HT", "H_HT") %in% names(all_iso))) {
  ov <- extract_t_h_overlap(
    t_iso = all_iso[["T_HT"]],
    h_iso = all_iso[["H_HT"]],
    iso_q_cutoff = iso_q_cutoff,
    gene_q_cutoff = gene_q_cutoff,
    min_abs_dif = min_abs_dif
  )
  utils::write.csv(
    ov$all,
    file.path(results_tables, "isoform_overlap_T_HT_vs_H_HT_all.csv"),
    row.names = FALSE,
    fileEncoding = "UTF-8",
    na = ""
  )
  utils::write.csv(
    ov$t_sig_in_h,
    file.path(results_tables, "isoform_overlap_T_significant_in_H_context.csv"),
    row.names = FALSE,
    fileEncoding = "UTF-8",
    na = ""
  )
  utils::write.csv(
    ov$summary,
    file.path(results_tables, "isoform_overlap_T_HT_vs_H_HT_summary.csv"),
    row.names = FALSE,
    fileEncoding = "UTF-8",
    na = ""
  )
}

message("06_visualization.R: done")
