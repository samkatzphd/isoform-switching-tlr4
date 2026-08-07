#!/usr/bin/env Rscript
# 02: Collapse isoform-level tables to gene-level summaries (n switching isoforms,
#     min q, novel involvement) using thresholds from config.
#
# Run: Rscript scripts/02_gene_level_summary.R
#  After 01_load_data.R. Reads data/processed/isoformFeatures_*.rds and writes
#  results/tables/gene_level_summary_*.[rds|csv]

bt <- c("utils/bootstrap.R", file.path("..", "utils", "bootstrap.R"))
b_file <- if (any(f <- vapply(bt, file.exists, NA))) { bt[which(f)[1L]] } else { NA_character_ }
if (is.na(b_file)) {
  stop("Could not find utils/bootstrap.R. See README (ISOFORM_PROJECT_ROOT or `cd` into project).")
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
in_dir <- resolve_path(
  paths$processed_dir %||% "data/processed", root = root
)
if (!dir.exists(in_dir)) {
  stop("Processed directory not found: ", in_dir, ". Run 01_load_data.R first.")
}
out_tab <- ensure_dir(
  resolve_path(
    paths$results_tables %||% "results/tables", root = root
  )
)
fl <- list.files(
  in_dir, pattern = "^isoformFeatures_.*\\.rds$", full.names = TRUE, ignore.case = FALSE
)
if (length(fl) < 1L) {
  stop("No isoformFeatures_*.rds under ", in_dir, ". Run 01_load_data.R first.")
}
for (rds in fl) {
  bn <- sub("\\.rds$", "", basename(rds), ignore.case = TRUE)
  # bn like isoformFeatures_T  ->  suffix = T
  key <- if (grepl("^isoformFeatures_(.+)$", bn, perl = TRUE)) {
    sub("^isoformFeatures_(.+)$", "\\1", bn, perl = TRUE)
  } else {
    "unknown"
  }
  message("Summarizing genes: [", key, "] ...")
  # Prefer the unfiltered context so significance and presence share one FDR universe.
  iso <- load_scoring_table(in_dir, key)
  if (is.null(iso)) next
  gtab <- summarize_genes_from_isoform_table(iso, cfg)
  # summarize_genes_from_isoform_table() collapses gene symbols first, so this must
  # hold. It did not before: disagreeing gene_name values split a gene across rows.
  stopifnot(nrow(gtab) == length(unique(gtab$gene_id)))
  write_table_pair(gtab, out_tab, paste0("gene_level_summary_", sanitize(key)), cfg = cfg)
}
write_run_manifest("02_gene_level_summary.R", cfg, root)
message("02_gene_level_summary.R: done")
