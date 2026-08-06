#!/usr/bin/env Rscript
# 05: Pathway enrichment (GO, Reactome, optionally KEGG) over switching gene sets.
#
# Run (after 01, 02, 04): Rscript scripts/05_pathway_enrichment.R
#
# THE UNIVERSE IS THE POINT OF THIS SCRIPT.
#
# This step was deliberately left unimplemented while the only gene lists available came
# from ISA objects saved with reduceToSwitchingGenes = TRUE. Every gene in those objects is
# already significant, so an enrichment test against them compares a significant set to a
# significant set and reports whatever the annotation happens to be enriched for. Testing
# against "all human genes" instead is just as wrong in the other direction: it treats genes
# that were never expressed in these cells as candidates that failed.
#
# The unfiltered context tables written by 01 give the correct denominator -- the genes
# actually quantified and tested in each experiment (8.2k-12.0k per dataset). That is the
# universe used throughout. If a dataset has no context table the script refuses to run for
# it rather than silently substituting a whole-genome background.
#
# Gene sets tested:
#   - per dataset: genes with >=1 switching isoform, against that dataset's tested genes
#   - UT overlap classes (shared / WT-only / KO-only), against genes tested in BOTH UT arms

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

need_pkgs <- c("clusterProfiler", "org.Hs.eg.db", "AnnotationDbi")
missing <- need_pkgs[!vapply(need_pkgs, requireNamespace, NA, quietly = TRUE)]
if (length(missing)) {
  stop(
    "Missing required packages: ", paste(missing, collapse = ", "), "\n",
    'Install with: BiocManager::install(c("clusterProfiler", "org.Hs.eg.db"))'
  )
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

processed_dir <- resolve_path(paths$processed_dir %||% "data/processed", root = root)
results_tables <- ensure_dir(resolve_path(paths$results_tables %||% "results/tables", root = root))
results_figures <- ensure_dir(resolve_path(paths$results_figures %||% "results/figures", root = root))
fig_dir <- ensure_dir(file.path(results_figures, "pathway"))

pw <- cfg$pathway %||% list()
ontologies <- as.character(pw$ontologies %||% c("BP", "MF"))
p_adj_method <- as.character(pw$p_adjust_method %||% "BH")[1L]
q_cut <- as.numeric(pw$q_value_cutoff %||% 0.05)
p_cut <- as.numeric(pw$p_value_cutoff %||% 0.05)
min_gs <- as.integer(pw$min_gs_size %||% 10L)
max_gs <- as.integer(pw$max_gs_size %||% 500L)
top_n_plot <- as.integer(pw$top_n_plot %||% 20L)
run_kegg <- isTRUE(pw$run_kegg %||% TRUE)
run_reactome <- isTRUE(pw$run_reactome %||% TRUE) &&
  requireNamespace("ReactomePA", quietly = TRUE)

theme_set(theme_bw(base_size = 11))

# ---- Symbol -> ENTREZ mapping -------------------------------------------------------
# Enrichment is done on ENTREZ ids. Mapping is lossy, so the rate is recorded per analysis:
# a low rate means the result describes only the mappable fraction and should be discounted.
map_symbols <- function(symbols) {
  symbols <- unique(as.character(symbols))
  symbols <- symbols[is_real_gene_symbol(symbols)]
  if (!length(symbols)) {
    return(tibble(SYMBOL = character(), ENTREZID = character()))
  }
  suppressMessages(suppressWarnings(
    tryCatch(
      AnnotationDbi::select(
        org.Hs.eg.db::org.Hs.eg.db,
        keys = symbols, keytype = "SYMBOL", columns = "ENTREZID"
      ) |>
        as_tibble() |>
        filter(!is.na(.data$ENTREZID)) |>
        distinct(.data$SYMBOL, .keep_all = TRUE),
      error = function(e) tibble(SYMBOL = character(), ENTREZID = character())
    )
  ))
}

enrich_one <- function(gene_entrez, universe_entrez, analysis, source_label) {
  out <- list()
  for (ont in ontologies) {
    r <- tryCatch(
      clusterProfiler::enrichGO(
        gene = gene_entrez, universe = universe_entrez,
        OrgDb = org.Hs.eg.db::org.Hs.eg.db, keyType = "ENTREZID", ont = ont,
        pAdjustMethod = p_adj_method, pvalueCutoff = p_cut, qvalueCutoff = q_cut,
        minGSSize = min_gs, maxGSSize = max_gs, readable = TRUE
      ),
      error = function(e) {
        message("    GO:", ont, " failed: ", conditionMessage(e))
        NULL
      }
    )
    if (!is.null(r) && nrow(as.data.frame(r))) {
      out[[paste0("GO_", ont)]] <- as.data.frame(r) |>
        mutate(database = paste0("GO:", ont), analysis = analysis, source = source_label)
    }
  }
  if (run_reactome) {
    r <- tryCatch(
      ReactomePA::enrichPathway(
        gene = gene_entrez, universe = universe_entrez, organism = "human",
        pAdjustMethod = p_adj_method, pvalueCutoff = p_cut, qvalueCutoff = q_cut,
        minGSSize = min_gs, maxGSSize = max_gs, readable = TRUE
      ),
      error = function(e) {
        message("    Reactome failed: ", conditionMessage(e))
        NULL
      }
    )
    if (!is.null(r) && nrow(as.data.frame(r))) {
      out[["Reactome"]] <- as.data.frame(r) |>
        mutate(database = "Reactome", analysis = analysis, source = source_label)
    }
  }
  if (run_kegg) {
    # KEGG needs network access; skipped rather than failing the run when offline.
    r <- tryCatch(
      clusterProfiler::enrichKEGG(
        gene = gene_entrez, universe = universe_entrez, organism = "hsa",
        pAdjustMethod = p_adj_method, pvalueCutoff = p_cut, qvalueCutoff = q_cut,
        minGSSize = min_gs, maxGSSize = max_gs
      ),
      error = function(e) {
        message("    KEGG unavailable (", conditionMessage(e), ") -- skipped.")
        NULL
      }
    )
    if (!is.null(r) && nrow(as.data.frame(r))) {
      # enrichKEGG has no `readable` argument, so its geneID column stays as ENTREZ ids.
      # Left alone, the results table mixes symbols and numeric ids for the same gene
      # (6348 and CCL3 both appear), which breaks any per-gene tallying downstream.
      r <- tryCatch(
        clusterProfiler::setReadable(r, OrgDb = org.Hs.eg.db::org.Hs.eg.db, keyType = "ENTREZID"),
        error = function(e) r
      )
      out[["KEGG"]] <- as.data.frame(r) |>
        mutate(database = "KEGG", analysis = analysis, source = source_label)
    }
  }
  if (!length(out)) {
    return(NULL)
  }
  common <- Reduce(intersect, lapply(out, names))
  bind_rows(lapply(out, function(d) d[, common, drop = FALSE]))
}

#' Run one enrichment analysis and record how much of it survived ID mapping.
run_analysis <- function(gene_symbols, universe_symbols, analysis, source_label) {
  gs <- intersect(unique(gene_symbols), unique(universe_symbols))
  if (length(gs) < length(unique(gene_symbols))) {
    message(
      "    ", length(unique(gene_symbols)) - length(gs),
      " gene(s) not in the tested universe were dropped (they must be a subset)."
    )
  }
  u_map <- map_symbols(universe_symbols)
  g_map <- map_symbols(gs)
  stats <- tibble(
    analysis = analysis,
    source = source_label,
    n_genes_in = length(unique(gs)),
    n_genes_mapped = nrow(g_map),
    pct_genes_mapped = if (length(gs)) 100 * nrow(g_map) / length(unique(gs)) else NA_real_,
    n_universe_in = length(unique(universe_symbols)),
    n_universe_mapped = nrow(u_map),
    pct_universe_mapped = if (length(universe_symbols)) {
      100 * nrow(u_map) / length(unique(universe_symbols))
    } else {
      NA_real_
    }
  )
  message(sprintf(
    "  [%s] %d genes (%.0f%% mapped) vs universe of %d (%.0f%% mapped)",
    analysis, stats$n_genes_in, stats$pct_genes_mapped,
    stats$n_universe_in, stats$pct_universe_mapped
  ))
  if (nrow(g_map) < min_gs) {
    message("    too few mapped genes (<", min_gs, "); skipping.")
    return(list(results = NULL, stats = stats |> mutate(n_terms = 0L)))
  }
  res <- enrich_one(g_map$ENTREZID, u_map$ENTREZID, analysis, source_label)
  n_terms <- if (is.null(res)) 0L else nrow(res)
  message("    ", n_terms, " enriched term(s) at q < ", q_cut)
  list(results = res, stats = stats |> mutate(n_terms = n_terms))
}

# ---- Assemble gene sets and universes ------------------------------------------------
context_universe <- function(label) {
  ctx <- load_context_table(processed_dir, label, "features")
  if (is.null(ctx)) {
    return(NULL)
  }
  # Collapse to one symbol per tested gene_id with the same rule used for the gene sets.
  # Taking unique gene_name over raw isoform rows would count symbol *variants* rather than
  # genes and inflate the universe (12,235 "genes" from 10,907 gene_ids in T_HT), biasing
  # every enrichment p-value against the gene sets, which are collapsed.
  g <- ctx |>
    group_by(.data$gene_id) |>
    summarize(gene_name = pick_gene_symbol(.data$gene_name), .groups = "drop")
  unique(as.character(g$gene_name[is_real_gene_symbol(g$gene_name)]))
}

analyses <- list()
for (k in names(cfg$datasets %||% list())) {
  ds <- cfg$datasets[[k]]
  label <- sanitize(as.character(ds$label %||% k)[1L])
  iso_rds <- file.path(processed_dir, paste0("isoformFeatures_", label, ".rds"))
  if (!file.exists(iso_rds)) {
    warning("[", k, "] processed table missing; run 01. Skipping.")
    next
  }
  universe <- context_universe(label)
  if (is.null(universe)) {
    warning(
      "[", k, "] no unfiltered context table, so no defensible universe exists. ",
      "Skipping rather than substituting a whole-genome background. ",
      "Set isa_unfiltered_path in config and rerun 01."
    )
    next
  }
  z <- score_isoforms(readRDS(iso_rds), cfg, k, label)
  switching <- z |>
    filter(.data$is_switching, is_real_gene_symbol(.data$gene_name)) |>
    pull(.data$gene_name) |>
    unique()
  analyses[[length(analyses) + 1L]] <- list(
    analysis = paste0(k, "_switching"),
    source = paste0("switching genes in ", k),
    genes = switching,
    universe = universe
  )
}

# UT overlap classes against the genes tested in BOTH UT arms. This is the contrast the
# project actually cares about: are knockout-specific switches enriched for anything the
# shared ones are not?
ut_ok <- file.exists(file.path(results_tables, "ut_gene_overlap_T_vs_U.rds"))
u_t <- context_universe(sanitize(as.character((cfg$datasets$T_UT %||% list())$label %||% "T_UT")))
u_u <- context_universe(sanitize(as.character((cfg$datasets$U_UT %||% list())$label %||% "U_UT")))
if (ut_ok && !is.null(u_t) && !is.null(u_u)) {
  shared_universe <- intersect(u_t, u_u)
  gover <- readRDS(file.path(results_tables, "ut_gene_overlap_T_vs_U.rds"))
  for (cls in c("Shared (T and U)", "T-only", "U-only")) {
    g <- gover |>
      filter(.data$overlap_class == cls, is_real_gene_symbol(.data$gene_name)) |>
      pull(.data$gene_name) |>
      unique()
    analyses[[length(analyses) + 1L]] <- list(
      analysis = paste0("UT_", sanitize(cls)),
      source = paste0(cls, " (universe = genes tested in both UT arms)"),
      genes = g,
      universe = shared_universe
    )
  }
} else {
  message("Skipping UT overlap-class enrichment (need 04 outputs and both UT contexts).")
}

if (!length(analyses)) {
  stop("No analyses could be assembled. Run 01 (with isa_unfiltered_path set) and 04 first.")
}

# ---- Run ------------------------------------------------------------------------------
message("Running enrichment for ", length(analyses), " gene set(s) ...")
all_res <- list()
all_stats <- list()
for (a in analyses) {
  out <- run_analysis(a$genes, a$universe, a$analysis, a$source)
  if (!is.null(out$results)) all_res[[a$analysis]] <- out$results
  all_stats[[a$analysis]] <- out$stats
}

stats_tbl <- bind_rows(all_stats)
write_table_pair(stats_tbl, results_tables, "pathway_enrichment_summary", cfg = cfg)
print(as.data.frame(stats_tbl))

if (!length(all_res)) {
  message(
    "\nNo enriched terms anywhere at q < ", q_cut, ".\n",
    "That is a result, not a failure: against a properly matched tested background, ",
    "switching gene sets are not enriched for the tested annotation categories."
  )
  write_run_manifest("05_pathway_enrichment.R", cfg, root)
  quit(status = 0L)
}

res_tbl <- bind_rows(all_res) |>
  arrange(.data$analysis, .data$database, .data$p.adjust)
write_table_pair(res_tbl, results_tables, "pathway_enrichment_results", cfg = cfg)

top_tbl <- res_tbl |>
  group_by(.data$analysis, .data$database) |>
  slice_min(.data$p.adjust, n = 10, with_ties = FALSE) |>
  ungroup()
write_table_pair(top_tbl, results_tables, "pathway_enrichment_top_terms", cfg = cfg)

# How many DISTINCT genes actually drive the enriched terms, and which ones.
#
# A term count is a bad summary on its own: GO and Reactome are nested hierarchies, so a
# handful of genes in one tight functional module produces dozens of overlapping terms. The
# shared UT set returns ~150 terms from 18 genes, but five chemokines/cytokines account for
# the large majority of them. Reporting the drivers alongside the count keeps that visible.
# NOTE: clusterProfiler attaches AnnotationDbi, whose select() masks dplyr's. Namespace
# the ambiguous verbs explicitly here.
terms_per_analysis <- res_tbl |>
  dplyr::count(.data$analysis, name = "n_terms_total")

driver_tbl <- res_tbl |>
  dplyr::select(.data$analysis, .data$database, .data$Description, .data$geneID) |>
  dplyr::mutate(gene = strsplit(as.character(.data$geneID), "/", fixed = TRUE)) |>
  tidyr::unnest(.data$gene) |>
  dplyr::count(.data$analysis, .data$gene, name = "n_terms") |>
  dplyr::group_by(.data$analysis) |>
  dplyr::mutate(n_distinct_driver_genes = dplyr::n()) |>
  dplyr::arrange(.data$analysis, dplyr::desc(.data$n_terms), .by_group = TRUE) |>
  dplyr::ungroup() |>
  dplyr::left_join(terms_per_analysis, by = "analysis") |>
  dplyr::mutate(pct_of_terms = 100 * .data$n_terms / .data$n_terms_total)
write_table_pair(driver_tbl, results_tables, "pathway_enrichment_driver_genes", cfg = cfg)

driver_summary <- driver_tbl |>
  dplyr::group_by(.data$analysis) |>
  dplyr::summarize(
    n_terms_total = dplyr::first(.data$n_terms_total),
    n_distinct_driver_genes = dplyr::first(.data$n_distinct_driver_genes),
    top_drivers = paste(utils::head(.data$gene, 5), collapse = ", "),
    top5_share_of_gene_hits = 100 * sum(utils::head(.data$n_terms, 5)) / sum(.data$n_terms),
    .groups = "drop"
  )
write_table_pair(driver_summary, results_tables, "pathway_enrichment_driver_summary", cfg = cfg)
print(as.data.frame(driver_summary))

# ---- Figures --------------------------------------------------------------------------
# Drop figures from a previous run first: an analysis that loses all its enriched terms
# (as T_HT did when the universe was corrected) would otherwise keep a stale plot on disk
# that no longer corresponds to any row in the results table.
stale <- list.files(fig_dir, pattern = "^fig_pathway_.*\\.png$", full.names = TRUE)
if (length(stale)) unlink(stale)

parse_ratio <- function(x) {
  p <- strsplit(as.character(x), "/", fixed = TRUE)
  vapply(p, function(v) if (length(v) == 2L) as.numeric(v[1]) / as.numeric(v[2]) else NA_real_, numeric(1))
}
for (a in unique(res_tbl$analysis)) {
  d <- res_tbl |>
    filter(.data$analysis == a) |>
    group_by(.data$database) |>
    slice_min(.data$p.adjust, n = top_n_plot, with_ties = FALSE) |>
    ungroup() |>
    mutate(
      gene_ratio = parse_ratio(.data$GeneRatio),
      Description = substr(as.character(.data$Description), 1, 60)
    )
  if (!nrow(d)) next
  p <- ggplot(d, aes(x = .data$gene_ratio, y = stats::reorder(.data$Description, .data$gene_ratio))) +
    geom_point(aes(size = .data$Count, colour = .data$p.adjust)) +
    scale_colour_gradient(low = "#6a51a3", high = "#c994c7", name = "q") +
    scale_size_continuous(name = "genes") +
    facet_wrap(~ .data$database, scales = "free_y") +
    labs(
      title = paste0("Enriched terms: ", a),
      subtitle = "Universe = genes quantified and tested in the same experiment",
      x = "Gene ratio", y = NULL
    ) +
    theme(axis.text.y = element_text(size = 7))
  n_rows <- max(4, min(nrow(d), top_n_plot * dplyr::n_distinct(d$database)))
  ggsave(
    file.path(fig_dir, paste0("fig_pathway_", sanitize(a), ".png")),
    p, width = 11, height = 2 + 0.22 * n_rows, dpi = 200, limitsize = FALSE
  )
}

write_run_manifest("05_pathway_enrichment.R", cfg, root)
message("05_pathway_enrichment.R: done")
message("  Tables: ", results_tables)
message("  Figures: ", fig_dir)
