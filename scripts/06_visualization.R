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
top_n_genes <- as.integer(analysis_param(cfg, "top_n_genes", 30L))
plot_top_n <- as.integer(analysis_param(cfg, "plot_top_n", 30L))

# sanitize(), min_finite(), is_real_gene_symbol(), pick_gene_symbol() and
# score_isoforms() come from utils/helper_functions.R -- this script used to keep its
# own copies alongside near-identical ones in 03 and 04.

compute_switching_gene_rank <- function(iso_tbl, dataset_key, dataset_label) {
  z <- score_isoforms(iso_tbl, cfg, dataset_key, dataset_label)

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

# ISA objects often store gene_name as XLOC placeholders. Real symbols and PB
# novelty live in our processed isoform table (after GFF3 join). Patch the
# in-memory ISA list so switchPlot titles/labels use symbols and mark novels.
annotate_isa_for_plotting <- function(isa_obj, processed_iso) {
  if (is.null(isa_obj) || is.null(processed_iso)) return(isa_obj)
  proc <- as_tibble(processed_iso)
  if (!all(c("gene_id", "isoform_id") %in% names(proc))) return(isa_obj)

  # Prefer symbols from processed table over XLOC/ENSG placeholders in ISA.
  if ("gene_name" %in% names(proc) && "isoformFeatures" %in% names(isa_obj)) {
    gmap <- proc |>
      group_by(.data$gene_id) |>
      summarize(
        gene_name = pick_gene_symbol(.data$gene_name),
        .groups = "drop"
      ) |>
      filter(is_real_gene_symbol(.data$gene_name)) |>
      transmute(
        gene_id = as.character(.data$gene_id),
        gene_name = as.character(.data$gene_name)
      )
    feat <- isa_obj$isoformFeatures
    feat$gene_id <- as.character(feat$gene_id)
    m <- match(feat$gene_id, gmap$gene_id)
    hit <- !is.na(m)
    if (any(hit)) {
      feat$gene_name[hit] <- gmap$gene_name[m[hit]]
      isa_obj$isoformFeatures <- feat
    }
  }

  # Relabel novel PacBio isoforms as "<TCONS> (PB)" across ID-bearing slots.
  if (!"is_novel_pacbio" %in% names(proc)) return(isa_obj)
  novel_ids <- unique(as.character(proc$isoform_id[proc$is_novel_pacbio %in% TRUE]))
  novel_ids <- novel_ids[!is.na(novel_ids) & nzchar(novel_ids)]
  if (!length(novel_ids)) return(isa_obj)
  # Avoid double-tagging on re-runs
  novel_ids <- novel_ids[!grepl(" \\(PB\\)$", novel_ids)]
  if (!length(novel_ids)) return(isa_obj)
  new_ids <- paste0(novel_ids, " (PB)")
  names(new_ids) <- novel_ids

  remap_ids <- function(ids) {
    ids <- as.character(ids)
    hit <- ids %in% names(new_ids)
    ids[hit] <- unname(new_ids[ids[hit]])
    ids
  }

  for (nm in names(isa_obj)) {
    x <- isa_obj[[nm]]
    if (is.data.frame(x) && "isoform_id" %in% names(x)) {
      x$isoform_id <- remap_ids(x$isoform_id)
      isa_obj[[nm]] <- x
    } else if (inherits(x, "GRanges")) {
      mc <- S4Vectors::mcols(x)
      if ("isoform_id" %in% names(mc)) {
        mc$isoform_id <- remap_ids(mc$isoform_id)
        S4Vectors::mcols(x) <- mc
        isa_obj[[nm]] <- x
      }
    } else if (inherits(x, "DNAStringSet") || inherits(x, "XStringSet")) {
      nms <- names(x)
      if (!is.null(nms)) {
        names(x) <- remap_ids(nms)
        isa_obj[[nm]] <- x
      }
    }
  }
  isa_obj
}

plot_gene_fallback <- function(gene_row, iso_tbl, dataset_label) {
  gid <- as.character(gene_row$gene_id[[1]])
  gname <- as.character(gene_row$gene_name[[1]])
  if (!is_real_gene_symbol(gname)) gname <- gid

  iso_tbl |>
    filter(.data$gene_id == gid) |>
    arrange(desc(abs(.data$dIF_n)), .data$q_i) |>
    slice_head(n = 12L) |>
    mutate(
      isoform_label = paste0(
        as.character(.data$isoform_id),
        ifelse(.data$is_novel_pacbio %in% TRUE, " (PB)", "")
      )
    ) |>
    ggplot(aes(x = reorder(.data$isoform_label, .data$dIF_n), y = .data$dIF_n, fill = .data$is_switching)) +
    geom_col() +
    coord_flip() +
    scale_fill_manual(values = c(`TRUE` = "#1b7837", `FALSE` = "#bdbdbd")) +
    labs(
      title = paste0(dataset_label, " | ", gname),
      subtitle = paste0("Fallback dIF view (switchPlot failed); gene_id=", gid),
      x = "Isoform",
      y = "dIF",
      fill = "Significant"
    ) +
    theme_bw(base_size = 10)
}

#' Draw an ISA switchPlot, trying the gene symbol first and the gene_id as a fallback.
#'
#' Device handling detail that used to be wrong here: `on.exit()` called inside a
#' `tryCatch({...})` registers on the *enclosing function* frame, not per iteration.
#' Devices therefore stayed open until the function returned and then closed LIFO, so
#' (a) the `file.info()$size` check ran before the PNG had been flushed, and (b) when
#' the symbol attempt failed and the gene_id attempt succeeded, the failed device --
#' opened first, closed last, on the same path -- overwrote the good plot. Each attempt
#' now renders to its own temp file, closes its device immediately, and is promoted to
#' the destination only after it is verified non-empty.
make_switch_plot <- function(isa_obj, gene_id, gene_name, out_png) {
  if (!requireNamespace("IsoformSwitchAnalyzeR", quietly = TRUE)) {
    return(list(ok = FALSE, message = "IsoformSwitchAnalyzeR is not installed."))
  }
  switch_fn <- get("switchPlot", envir = asNamespace("IsoformSwitchAnalyzeR"))
  # Prefer human-readable gene_name so the plot title is a symbol, not XLOC.
  candidates <- character()
  if (is_real_gene_symbol(gene_name)) candidates <- c(candidates, as.character(gene_name)[1L])
  if (!is.na(gene_id) && nzchar(as.character(gene_id)[1L])) {
    candidates <- c(candidates, as.character(gene_id)[1L])
  }
  candidates <- unique(candidates)
  if (!length(candidates)) {
    return(list(ok = FALSE, message = "No gene identifier available."))
  }

  for (g in candidates) {
    tmp <- tempfile(fileext = ".png")
    ok <- tryCatch({
      grDevices::png(filename = tmp, width = 2400, height = 1500, res = 180)
      on.exit(if (grDevices::dev.cur() > 1L) grDevices::dev.off(), add = TRUE)
      do.call(
        switch_fn,
        list(
          switchAnalyzeRlist = isa_obj,
          gene = g,
          plotTopology = FALSE
        )
      )
      grDevices::dev.off()
      TRUE
    }, error = function(e) {
      FALSE
    })
    if (grDevices::dev.cur() > 1L) grDevices::dev.off()
    if (ok && file.exists(tmp) && isTRUE(file.info(tmp)$size > 0)) {
      file.copy(tmp, out_png, overwrite = TRUE)
      unlink(tmp)
      return(list(ok = TRUE, message = paste0("switchPlot succeeded for ", g)))
    }
    unlink(tmp)
  }
  list(ok = FALSE, message = "switchPlot failed for gene_name and gene_id candidates.")
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
  write_table_pair(
    top_tbl, results_tables, paste0("top_switching_genes_", label_clean), cfg = cfg
  )

  isa_obj <- load_isa_for_dataset(ds)
  if (is.null(isa_obj)) {
    # Skip plotting entirely rather than replacing existing switchPlots with fallback
    # figures. The ISA objects live on an external drive; running this script without
    # it mounted used to delete the real plots and write dIF bar charts over them.
    warning(
      "[", ds_key, "] ISA object unavailable (is the external drive mounted?). ",
      "Keeping existing figures and skipping all plotting for this dataset."
    )
    next
  }
  message("[", ds_key, "] Annotating ISA object with gene symbols and PB isoform labels ...")
  isa_obj <- annotate_isa_for_plotting(isa_obj, iso_scored)

  # Drop stale rank plots from earlier ranking schemas before rewriting.
  old_plots <- list.files(
    out_fig_dir,
    pattern = paste0("^switch_plot_", label_clean, "_rank[0-9]{2}_.*\\.png$"),
    full.names = TRUE
  )
  if (length(old_plots)) {
    unlink(old_plots)
  }
  # Remove any leftover test plots
  unlink(file.path(out_fig_dir, "_test_NCOA7_labels.png"))

  gene_rows <- top_tbl |>
    slice_head(n = min(plot_top_n, nrow(top_tbl)))
  if (!nrow(gene_rows)) {
    message("[", ds_key, "] No significant switching genes to plot.")
    next
  }
  for (i in seq_len(nrow(gene_rows))) {
    gr <- gene_rows[i, , drop = FALSE]
    gene_stub <- sanitize(
      if (is_real_gene_symbol(gr$gene_name[[1]])) gr$gene_name[[1]] else gr$gene_id[[1]]
    )
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
        message("[", ds_key, "] switchPlot unavailable for ", gene_stub, "; writing fallback plot.")
      } else {
        message("[", ds_key, "] Wrote ", basename(out_png), " (", sp$message, ")")
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
  write_table_pair(combined, results_tables, "top_switching_genes_HT_combined", cfg = cfg)
}

if (all(c("T_HT", "H_HT") %in% names(all_iso))) {
  ov <- extract_t_h_overlap(
    t_iso = all_iso[["T_HT"]],
    h_iso = all_iso[["H_HT"]],
    iso_q_cutoff = iso_q_cutoff,
    gene_q_cutoff = gene_q_cutoff,
    min_abs_dif = min_abs_dif
  )
  write_table_pair(
    ov$all, results_tables, "isoform_overlap_T_HT_vs_H_HT_all", cfg = cfg,
    csv = isTRUE((cfg$output %||% list())$csv_twin_isoform_level %||% FALSE)
  )
  write_table_pair(
    ov$t_sig_in_h, results_tables, "isoform_overlap_T_significant_in_H_context", cfg = cfg
  )
  write_table_pair(
    ov$summary, results_tables, "isoform_overlap_T_HT_vs_H_HT_summary", cfg = cfg
  )
}

write_run_manifest("06_visualization.R", cfg, root)
message("06_visualization.R: done")
