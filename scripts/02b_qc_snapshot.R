#!/usr/bin/env Rscript
# 02b: Quick QA snapshot of gene-level outputs.
# - Percent genes with novel isoform involvement
# - Top 10 genes per dataset (ranked by switching strength + effect size)

bt <- c("utils/bootstrap.R", file.path("..", "utils", "bootstrap.R"))
b_file <- if (any(f <- vapply(bt, file.exists, NA))) { bt[which(f)[1L]] } else { NA_character_ }
if (is.na(b_file)) {
  stop("Could not find utils/bootstrap.R. Run from project root or set ISOFORM_PROJECT_ROOT.")
}
source(b_file, local = FALSE, chdir = FALSE)
root <- if (exists("PROJECT_ROOT", inherits = FALSE)) {
  get("PROJECT_ROOT", inherits = FALSE)
} else {
  find_project_root()
}

cfg <- load_yaml_config("config/config.yml", root = root)
paths <- cfg$paths %||% list()
in_dir <- resolve_path(paths$results_tables %||% "results/tables", root = root)
if (!dir.exists(in_dir)) {
  stop("Results table directory not found: ", in_dir)
}

fl <- list.files(
  in_dir,
  pattern = "^gene_level_summary_.*\\.csv$",
  full.names = TRUE,
  ignore.case = FALSE
)
if (!length(fl)) {
  stop("No gene_level_summary_*.csv found in ", in_dir, ". Run scripts/02_gene_level_summary.R first.")
}

qa_rows <- list()
top_rows <- list()

for (csv in fl) {
  d <- utils::read.csv(csv, stringsAsFactors = FALSE, check.names = FALSE)
  key <- sub("^gene_level_summary_(.+)\\.csv$", "\\1", basename(csv), perl = TRUE)
  req <- c("gene_id", "gene_name", "n_switching_isoforms", "max_abs_dif", "novel_involved")
  miss <- setdiff(req, names(d))
  if (length(miss)) {
    warning("Skipping ", basename(csv), " (missing columns: ", paste(miss, collapse = ", "), ")")
    next
  }

  d$novel_involved <- d$novel_involved %in% TRUE
  d$n_switching_isoforms <- as.numeric(d$n_switching_isoforms)
  d$max_abs_dif <- as.numeric(d$max_abs_dif)
  d$min_isoform_switch_q <- if ("min_isoform_switch_q" %in% names(d)) as.numeric(d$min_isoform_switch_q) else NA_real_

  n_genes <- nrow(d)
  n_novel <- sum(d$novel_involved, na.rm = TRUE)
  pct_novel <- if (n_genes > 0L) 100 * n_novel / n_genes else NA_real_
  n_switching <- sum(d$n_switching_isoforms > 0, na.rm = TRUE)

  qa_rows[[length(qa_rows) + 1L]] <- data.frame(
    dataset = key,
    n_genes = n_genes,
    n_switching_genes = n_switching,
    n_novel_involved = n_novel,
    pct_novel_involved = round(pct_novel, 2),
    stringsAsFactors = FALSE
  )

  ord <- order(
    -d$n_switching_isoforms,
    -d$max_abs_dif,
    d$min_isoform_switch_q,
    na.last = TRUE
  )
  top <- d[ord, c("gene_id", "gene_name", "n_switching_isoforms", "max_abs_dif", "min_isoform_switch_q", "novel_involved")]
  top <- utils::head(top, 10L)
  top$dataset <- key
  top <- top[, c("dataset", "gene_id", "gene_name", "n_switching_isoforms", "max_abs_dif", "min_isoform_switch_q", "novel_involved")]
  top_rows[[length(top_rows) + 1L]] <- top
}

qa <- do.call(rbind, qa_rows)
top10 <- do.call(rbind, top_rows)

qa_out <- file.path(in_dir, "qc_snapshot_gene_level.csv")
top_out <- file.path(in_dir, "qc_top10_genes_per_dataset.csv")
utils::write.csv(qa, qa_out, row.names = FALSE, fileEncoding = "UTF-8", na = "")
utils::write.csv(top10, top_out, row.names = FALSE, fileEncoding = "UTF-8", na = "")

message("Wrote: ", qa_out)
message("Wrote: ", top_out)
print(qa)
