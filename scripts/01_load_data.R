#!/usr/bin/env Rscript
# 01: Load SwitchAnalyzeR list objects, extract isoform features, tag novel isoforms, save
#     processed isoform tibbles for downstream scripts.
#
# Run: Rscript scripts/01_load_data.R
#      Rscript scripts/01_load_data.R  path/to/optional/config.yml
#  From project root, or: cd scripts && Rscript 01_load_data.R  (see README)

bootstrap_b <- c(
  "utils/bootstrap.R",
  file.path("..", "utils", "bootstrap.R")
)
b_file <- {
  w <- vapply(bootstrap_b, file.exists, NA)
  if (any(w)) { bootstrap_b[which(w)[1L]] } else { NA_character_ }
}
if (is.na(b_file)) {
  stop(
    "Run from project root or from scripts/ so that utils/bootstrap.R is found, ",
    "or set ISOFORM_PROJECT_ROOT to the project path."
  )
}
source(b_file, local = FALSE, chdir = FALSE)
root <- if (exists("PROJECT_ROOT", inherits = FALSE)) {
  get("PROJECT_ROOT", inherits = FALSE)
} else {
  find_project_root()
}
args_r <- commandArgs(trailingOnly = TRUE)
config_rel <- if (length(args_r) > 0L) { args_r[[1L]] } else { "config/config.yml" }
cfg <- load_yaml_config(config_rel, root = root)
paths <- cfg$paths %||% list()
if (is.null(paths$processed_dir) || (is.character(paths$processed_dir) && !nzchar(paths$processed_dir))) {
  paths$processed_dir <- "data/processed"
}
proc <- ensure_dir(resolve_path(paths$processed_dir, root = root))
nov <- cfg$novel %||% list()
novel_id_col <- nov$id_column %||% "oId"
novel_pre <- nov$pb_prefix %||% "PB"
datasets <- cfg$datasets %||% NULL
if (is.null(datasets) || !length(datasets)) {
  # Backward compatibility with original T/U/H style config
  labels <- as.list(cfg$labels %||% list(T = "T", U = "U", H = "H"))
  object_names <- as.list(cfg$object_names %||% list(T = NULL, U = NULL, H = NULL))
  annotation_paths <- as.list(cfg$annotation_paths %||% list(T = NULL, U = NULL, H = NULL))
  datasets <- list(
    T = list(
      isa_path = paths$T_isa_list %||% NULL,
      label = labels$T %||% "T",
      object_name = object_names$T %||% NULL,
      annotation_path = annotation_paths$T %||% NULL
    ),
    U = list(
      isa_path = paths$U_isa_list %||% NULL,
      label = labels$U %||% "U",
      object_name = object_names$U %||% NULL,
      annotation_path = annotation_paths$U %||% NULL
    ),
    H = list(
      isa_path = paths$H_isa_list %||% NULL,
      label = labels$H %||% "H",
      object_name = object_names$H %||% NULL,
      annotation_path = annotation_paths$H %||% NULL
    )
  )
}
reduction_rows <- list()

for (k in names(datasets)) {
  ds <- datasets[[k]] %||% list()
  raw <- ds$isa_path %||% ds$path %||% NULL
  is_empty_path <- is.null(raw) || identical(as.character(raw)[1L], "null") ||
    (is.character(raw) && (length(raw) < 1L || !any(nzchar(na.omit(as.character(raw))))))
  if (is_empty_path) {
    message("Skipping ", k, " (isa_path missing).")
    next
  }
  abs_in <- resolve_path(as.character(raw)[1L], root = root)
  if (is.null(abs_in) || is.na(abs_in) || !file.exists(abs_in)) {
    stop("Input for dataset ", k, " not found: ", as.character(raw)[1L], " (resolved: ", abs_in, ")")
    next
  }
  key_label <- as.character(ds$label %||% k)[1L]
  message("Loading [", k, " / ", key_label, "] from ", abs_in, " ...")
  obj_name <- ds$object_name %||% NULL
  isa <- load_isa_input(abs_in, object_name = obj_name)
  iso <- extract_isoform_features(isa, .verbose = TRUE)

  # Map transcript metadata from the reference GFF3 (transcript_id -> oId / cmp_ref /
  # class_code / gene symbol). Previously this whole block was skipped whenever the ISA
  # table already carried `oId`, which would silently drop class_code, cmp_ref and the
  # XLOC -> symbol repair too. Now each column is joined only if it is actually missing.
  ann_raw <- ds$annotation_path %||% NULL
  ann_path <- if (!is.null(ann_raw) && nzchar(as.character(ann_raw)[1L])) {
    resolve_path(as.character(ann_raw)[1L], root = root)
  } else {
    NULL
  }
  if ("isoform_id" %in% names(iso) && !is.null(ann_path) && file.exists(ann_path)) {
    tx_map <- build_transcript_map_from_gtf(ann_path)
    want <- c(novel_id_col, "cmp_ref", "class_code", "gene_name_from_gtf")
    want <- intersect(unique(want), names(tx_map))
    missing_cols <- setdiff(want, names(iso))
    if (length(missing_cols)) {
      message(
        "  Joining annotation map from ", ann_path,
        " (adding: ", paste(missing_cols, collapse = ", "), ")"
      )
      iso <- dplyr::left_join(
        iso,
        tx_map[, c("transcript_id", missing_cols), drop = FALSE],
        by = c("isoform_id" = "transcript_id")
      )
    } else {
      message("  Annotation columns already present; no join needed.")
    }
    # Prefer gene symbols from annotation over XLOC/ENS placeholders in ISA.
    if ("gene_name_from_gtf" %in% names(iso)) {
      if (!"gene_name_original" %in% names(iso) && "gene_name" %in% names(iso)) {
        iso$gene_name_original <- iso$gene_name
      }
      isa_name <- if ("gene_name" %in% names(iso)) as.character(iso$gene_name) else rep(NA_character_, nrow(iso))
      gtf_name <- as.character(iso$gene_name_from_gtf)
      isa_placeholder <- !is_real_gene_symbol(isa_name)
      gtf_ok <- is_real_gene_symbol(gtf_name)
      if (!"gene_name" %in% names(iso)) {
        iso$gene_name <- gtf_name
      } else {
        iso$gene_name <- ifelse(isa_placeholder & gtf_ok, gtf_name, isa_name)
      }
    }
  } else if (is.null(ann_path) || !file.exists(ann_path)) {
    message("  No usable annotation_path for ", k, "; skipping transcript map join.")
  }
  iso <- tag_novel_isoforms(
    iso, col = novel_id_col, prefix = novel_pre, out_name = "is_novel_pacbio"
  )
  # Record whether this object was already reduced to significant switching genes
  # before it was saved (ISA's reduceToSwitchingGenes = TRUE default). Everything
  # downstream that counts genes or defines a background depends on this.
  red <- detect_isa_reduction(iso, cfg)
  reduction_rows[[k]] <- data.frame(
    dataset_key = k,
    dataset_label = key_label,
    n_isoforms = nrow(iso),
    n_genes = red$n_genes,
    n_genes_gene_q_significant = red$n_genes_significant,
    pct_genes_gene_q_significant = round(red$pct_genes_significant, 2),
    looks_reduced_to_switching_genes = red$looks_reduced,
    stringsAsFactors = FALSE
  )
  if (isTRUE(red$looks_reduced)) {
    message(
      "  NOTE: all ", red$n_genes, " genes already pass gene q < ",
      cfg$significance$gene_q, " -- this object was saved after ",
      "reduceToSwitchingGenes = TRUE. 'Genes present' therefore means 'genes already ",
      "called significant'; see docs/REVIEW_CHANGES.md."
    )
  }

  write_table_pair(
    iso, proc, paste0("isoformFeatures_", sanitize(key_label)),
    cfg = cfg,
    csv = isTRUE((cfg$output %||% list())$csv_twin_isoform_level %||% FALSE)
  )
}

reduction_summary <- do.call(rbind, reduction_rows)
if (!is.null(reduction_summary)) {
  tab_dir <- ensure_dir(resolve_path(paths$results_tables %||% "results/tables", root = root))
  write_table_pair(reduction_summary, tab_dir, "input_object_reduction_check", cfg = cfg)
  print(reduction_summary)
}
write_run_manifest(
  "01_load_data.R", cfg, root,
  extra = list(inputs = reduction_rows)
)
message("01_load_data.R: done")
