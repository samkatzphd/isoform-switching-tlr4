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
  # If ISA table lacks oId (or configured novel id column), map transcript metadata
  # from GTF/GFF-derived transcript annotations (transcript_id -> oId/cmp_ref/class_code).
  if (!novel_id_col %in% names(iso) && "isoform_id" %in% names(iso)) {
    ann_raw <- ds$annotation_path %||% NULL
    ann_path <- if (!is.null(ann_raw) && nzchar(as.character(ann_raw)[1L])) {
      resolve_path(as.character(ann_raw)[1L], root = root)
    } else {
      NULL
    }
    if (!is.null(ann_path) && file.exists(ann_path)) {
      message("  Joining annotation map from ", ann_path)
      tx_map <- build_transcript_map_from_gtf(ann_path)
      iso <- dplyr::left_join(
        iso,
        tx_map,
        by = c("isoform_id" = "transcript_id")
      )
    } else {
      message("  No annotation_path for ", k, "; skipping transcript map join.")
    }
  }
  iso <- tag_novel_isoforms(
    iso, col = novel_id_col, prefix = novel_pre, out_name = "is_novel_pacbio"
  )
  base <- file.path(
    proc,
    paste0(
      "isoformFeatures_",
      gsub("[^A-Za-z0-9_]+", "_", key_label, perl = TRUE)
    )
  )
  rds_name <- paste0(base, ".rds")
  csv_name <- paste0(base, ".csv")
  saveRDS(iso, rds_name, compress = "xz")
  utils::write.csv(iso, file = csv_name, row.names = FALSE, fileEncoding = "UTF-8", na = "")
  message("  Wrote: ", rds_name)
  message("  Wrote: ", csv_name)
}
message("01_load_data.R: done")
