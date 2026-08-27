#!/usr/bin/env Rscript
# 17: What does the literature already say about each candidate -- and specifically about
#     its ISOFORMS?
#
# Run (after 15): Rscript scripts/17_candidate_literature.R
# Needs network on FIRST run only -- PubMed responses are cached under
# data/external/pubmed/ and committed, so re-runs are offline.
#
# WHY A SCRIPT AND NOT A READING LIST
#
# The question a candidate has to survive before anyone orders primers is "has someone
# already characterized this isoform?". That is a different question from "is this gene
# studied in macrophages", and the two come apart hard: a gene can carry 4,000 immunology
# papers and none about its transcript structure. Splitting the two counts makes the
# distinction visible per gene instead of leaving it to whoever reads the table.
#
# THREE QUERIES PER GENE
#
#   isoform    <gene> AND (isoform OR splicing OR "splice variant" OR "alternative promoter")
#              -- prior art on the transcript structure itself. This is the one that decides
#              whether a candidate is novel or already described.
#   immunity   <gene> AND (macrophage OR LPS OR "innate immunity" OR inflammation OR TLR)
#              -- is the gene known in this system at all?
#   both       the intersection: an isoform characterized IN this system. A non-zero count
#              here is the strongest possible prior, and the rarest.
#
# WHAT THE COUNTS ARE AND ARE NOT
#
# These are PubMed hit counts on a title/abstract query, so they are a crude proxy for
# attention, sensitive to gene-symbol ambiguity (short symbols collide with English words and
# with other species' genes) and blind to anything phrased differently. `query_note` flags
# symbols short enough to be collision-prone. Use the counts to TRIAGE and the returned titles
# to judge; never report a count as a measure of how well understood a gene is.

bt <- c("utils/bootstrap.R", file.path("..", "utils", "bootstrap.R"))
b_file <- if (any(f <- vapply(bt, file.exists, NA))) { bt[which(f)[1L]] } else { NA_character_ }
if (is.na(b_file)) {
  stop("Could not find utils/bootstrap.R. Run from project root or set ISOFORM_PROJECT_ROOT.")
}
source(b_file, local = FALSE, chdir = FALSE)

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
})
if (!requireNamespace("jsonlite", quietly = TRUE)) {
  stop("jsonlite is required for 17_candidate_literature.R")
}

root <- if (exists("PROJECT_ROOT", inherits = FALSE)) {
  get("PROJECT_ROOT", inherits = FALSE)
} else {
  find_project_root()
}
args_r <- commandArgs(trailingOnly = TRUE)
config_rel <- if (length(args_r) > 0L) args_r[[1L]] else "config/config.yml"
cfg <- load_yaml_config(config_rel, root = root)
paths <- cfg$paths %||% list()

results_tables <- ensure_dir(resolve_path(paths$results_tables %||% "results/tables", root = root))
cache_dir <- ensure_dir(resolve_path("data/external/pubmed", root = root))

pairs_f <- file.path(results_tables, "candidate_switch_pairs.rds")
if (!file.exists(pairs_f)) {
  stop("Missing candidate_switch_pairs.rds. Run: Rscript scripts/15_isoform_characterization.R")
}
genes <- sort(unique(readRDS(pairs_f)$gene_name))
message("17_candidate_literature.R: querying PubMed for ", length(genes), " candidate genes")

EUTILS <- "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
ISOFORM_TERMS <- '("isoform" OR "splicing" OR "splice variant" OR "alternative promoter" OR "transcript variant")'
IMMUNE_TERMS <- '("macrophage" OR "LPS" OR "lipopolysaccharide" OR "innate immunity" OR "inflammation" OR "toll-like receptor")'

# NCBI asks for <=3 requests/second without a key. One gene needs three esearch calls plus one
# esummary, so a fixed pause keeps this well inside the limit for a 40-gene run.
polite <- function() Sys.sleep(0.4)

cached_json <- function(f, url) {
  if (file.exists(f)) {
    return(tryCatch(jsonlite::fromJSON(f, simplifyVector = FALSE), error = function(e) NULL))
  }
  ok <- tryCatch({
    utils::download.file(url, f, quiet = TRUE, method = "libcurl")
    TRUE
  }, error = function(e) FALSE, warning = function(w) FALSE)
  polite()
  if (!ok || !file.exists(f)) return(NULL)
  tryCatch(jsonlite::fromJSON(f, simplifyVector = FALSE), error = function(e) NULL)
}

esearch <- function(gene, tag, term) {
  f <- file.path(cache_dir, paste0(gene, "__", tag, ".json"))
  url <- paste0(EUTILS, "/esearch.fcgi?db=pubmed&retmode=json&retmax=8&sort=relevance&term=",
                utils::URLencode(term, reserved = TRUE))
  j <- cached_json(f, url)
  if (is.null(j) || is.null(j$esearchresult)) return(list(n = NA_integer_, ids = character(0)))
  list(
    n = suppressWarnings(as.integer(j$esearchresult$count %||% NA)),
    ids = unlist(j$esearchresult$idlist %||% list())
  )
}

esummary <- function(gene, tag, ids) {
  if (!length(ids)) return(tibble())
  f <- file.path(cache_dir, paste0(gene, "__", tag, "__summary.json"))
  url <- paste0(EUTILS, "/esummary.fcgi?db=pubmed&retmode=json&id=", paste(ids, collapse = ","))
  j <- cached_json(f, url)
  if (is.null(j) || is.null(j$result)) return(tibble())
  bind_rows(lapply(ids, function(id) {
    r <- j$result[[id]]
    if (is.null(r)) return(NULL)
    tibble(
      pmid = id,
      year = suppressWarnings(substr(as.character(r$pubdate %||% ""), 1, 4)),
      journal = as.character(r$source %||% NA_character_),
      title = as.character(r$title %||% NA_character_)
    )
  }))
}

counts <- list(); refs <- list()
for (g in genes) {
  ta <- function(x) paste0(x, "[Title/Abstract]")
  q_iso <- paste0(ta(g), " AND ", ISOFORM_TERMS)
  q_imm <- paste0(ta(g), " AND ", IMMUNE_TERMS)
  q_both <- paste0(ta(g), " AND ", ISOFORM_TERMS, " AND ", IMMUNE_TERMS)

  r_iso <- esearch(g, "isoform", q_iso)
  r_imm <- esearch(g, "immune", q_imm)
  r_both <- esearch(g, "both", q_both)
  r_all <- esearch(g, "all", ta(g))

  counts[[length(counts) + 1L]] <- tibble(
    gene_name = g,
    n_papers_total = r_all$n,
    n_isoform_papers = r_iso$n,
    n_immunity_papers = r_imm$n,
    n_isoform_and_immunity = r_both$n,
    # A three-letter symbol matches a great deal of text that is not this gene.
    query_note = if (nchar(g) <= 4L) {
      "short symbol -- counts likely inflated by ambiguous matches"
    } else {
      NA_character_
    }
  )
  for (tag in c("both", "isoform")) {
    ids <- if (tag == "both") r_both$ids else r_iso$ids
    s <- esummary(g, tag, ids)
    if (nrow(s)) {
      refs[[length(refs) + 1L]] <- s |> mutate(gene_name = g, query = tag)
    }
  }
  message("  ", g, ": ", r_all$n %||% NA, " total / ", r_iso$n %||% NA, " isoform / ",
          r_imm$n %||% NA, " immunity / ", r_both$n %||% NA, " both")
}

lit <- bind_rows(counts)
ref_tab <- bind_rows(refs)
if (nrow(ref_tab)) {
  ref_tab <- ref_tab |> select("gene_name", "query", "pmid", "year", "journal", "title")
}

write_table_pair(lit, results_tables, "candidate_literature_counts", cfg = cfg)
if (nrow(ref_tab)) {
  write_table_pair(ref_tab, results_tables, "candidate_literature_references", cfg = cfg)
}

message("\n---- literature attention, ordered by isoform-specific prior art ----")
print(as.data.frame(
  lit |> arrange(desc(.data$n_isoform_and_immunity), desc(.data$n_isoform_papers)) |>
    select("gene_name", "n_papers_total", "n_isoform_papers", "n_immunity_papers",
           "n_isoform_and_immunity")
), row.names = FALSE)

gap <- lit |> filter(.data$n_isoform_papers %in% c(0L, 1L, 2L), .data$n_immunity_papers >= 20L)
if (nrow(gap)) {
  message("\n---- studied in this system, but essentially NO isoform-level prior art ----")
  message("     (the openings: a well-known gene whose transcript structure nobody has looked at)")
  print(as.data.frame(gap |> select("gene_name", "n_immunity_papers", "n_isoform_papers")),
        row.names = FALSE)
}

write_run_manifest(
  script = "17_candidate_literature.R", cfg = cfg, root = root,
  extra = list(pubmed_cache = "data/external/pubmed", n_genes = length(genes))
)
message("\n17_candidate_literature.R: done")
