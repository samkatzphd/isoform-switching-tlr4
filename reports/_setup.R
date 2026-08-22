# Shared setup for every report in reports/.
#
# Each .qmd used to carry its own copy of the project-root discovery, and the four
# copies had drifted into three different strategies: absolute paths via proj_path()
# (overview), knitr root.dir (ht / novel_isoform) and normalizePath("..") (ut).
# Source this instead:
#
#   ```{r}
#   #| label: setup
#   source("_setup.R")
#   fig <- make_fig("novel")
#   ```
#
# Paths are absolute. Do NOT set knitr's root.dir here: it breaks the combination of
# include_graphics() and embed-resources, which is why the reports diverged originally.

.infer_project_root <- function(start = getwd()) {
  if (nzchar(Sys.getenv("ISOFORM_PROJECT_ROOT", ""))) {
    return(normalizePath(Sys.getenv("ISOFORM_PROJECT_ROOT"), mustWork = TRUE))
  }
  p <- normalizePath(start, mustWork = FALSE)
  for (i in seq_len(20L)) {
    if (file.exists(file.path(p, "config", "config.yml"))) {
      return(normalizePath(p, mustWork = TRUE))
    }
    parent <- dirname(p)
    if (identical(p, parent)) break
    p <- parent
  }
  stop(
    "Could not find config/config.yml. Render from inside the repo, ",
    "or set ISOFORM_PROJECT_ROOT."
  )
}

PROJECT_ROOT <- .infer_project_root()

proj_path <- function(...) normalizePath(file.path(PROJECT_ROOT, ...), mustWork = FALSE)
tab <- function(...) proj_path("results", "tables", ...)
figs <- function(...) proj_path("results", "figures", ...)

#' Build a figure-path helper bound to one subdirectory of results/figures
make_fig <- function(subdir = NULL) {
  if (is.null(subdir)) {
    function(...) figs(...)
  } else {
    function(...) figs(subdir, ...)
  }
}

#' Stop with an actionable message when a report's inputs are missing
require_outputs <- function(files, hint) {
  abs <- vapply(files, tab, character(1))
  miss <- abs[!file.exists(abs)]
  if (length(miss)) {
    stop(hint, "\nMissing:\n", paste(miss, collapse = "\n"), call. = FALSE)
  }
  abs
}

read_tab <- function(name) {
  utils::read.csv(tab(name), stringsAsFactors = FALSE, check.names = FALSE)
}

# --- External review artefacts -------------------------------------------------
# The 2026-08 review wrote its own tables and figures under docs/external_review/
# rather than results/, because they came from analyses that are not (yet) pipeline
# scripts. A report that tells the whole story has to cite both trees, so give the
# review side the same three helpers rather than letting each .qmd invent paths.

xr_tab  <- function(...) proj_path("docs", "external_review", "tables", ...)
xr_figs <- function(...) proj_path("docs", "external_review", "figures", ...)

xr_read <- function(name) {
  utils::read.csv(xr_tab(name), stringsAsFactors = FALSE, check.names = FALSE)
}

#' Pull one value out of a two-column statistic/value table by its key.
#' Returns NA rather than erroring, so a renamed row degrades to a visible gap
#' instead of killing the render.
xr_stat <- function(d, key, col = "value", key_col = "statistic") {
  i <- match(key, d[[key_col]])
  if (is.na(i)) NA_real_ else as.numeric(d[[col]][i])
}

#' Keep only figure paths that exist, so a report never renders a broken image
existing_figs <- function(paths) paths[file.exists(paths)]

#' Standard caveat about the input objects, rendered into any report that counts
#' genes or compares datasets. Reads the flags recorded by 01_load_data.R rather
#' than hard-coding the claim.
reduction_caveat <- function() {
  f <- tab("input_object_reduction_check.csv")
  if (!file.exists(f)) {
    return("")
  }
  d <- utils::read.csv(f, stringsAsFactors = FALSE, check.names = FALSE)
  reduced <- d$looks_reduced_to_switching_genes %in% c(TRUE, "TRUE", "true")
  if (!any(reduced)) {
    return("")
  }
  paste0(
    "::: {.callout-important title=\"Inputs are already reduced to switching genes\"}\n",
    "All ", sum(reduced), " of ", nrow(d), " input ISA objects contain **only genes that ",
    "already pass the gene-level q cutoff** (saved after ",
    "`isoformSwitchTestDEXSeq(reduceToSwitchingGenes = TRUE)`).\n\n",
    "Consequences for everything below: gene counts are *genes retained*, not *genes ",
    "tested*; \"present in the other dataset\" means *retained in that saved object*, not ",
    "*tested and non-significant*; and any percentage of genes has a denominator that was ",
    "already selected on the outcome. Overlap between datasets is bounded by how many ",
    "symbols survive in both objects. No enrichment test against this background is ",
    "reported. See `docs/REVIEW_CHANGES.md`.\n",
    ":::\n"
  )
}
