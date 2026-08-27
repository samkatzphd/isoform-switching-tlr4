#!/usr/bin/env Rscript
# 16: Which annotated functional sites does each switch actually remove or add?
#
# Run (after 15): Rscript scripts/16_functional_sites.R
# Needs network on FIRST run only -- UniProt responses are cached under
# data/external/uniprot/ and committed, so re-runs are offline. The external drive is
# not needed.
#
# WHY THIS IS A SEPARATE STEP FROM 15
#
# 15 says the switch removes 390 codons from the N-terminus. That is a fact about the
# transcript. Whether it matters depends on what those codons encoded, and the pipeline has
# no protein-feature annotation of its own: analyzeORF ran upstream but Pfam, SignalP and
# IUPred2A never did, so `domainAnalysis`, `signalPeptideAnalysis` and `idrAnalysis` are
# absent from all eight objects. UniProt is the substitute, and a good one for this purpose --
# these are all well-studied human genes with curated domain boundaries.
#
# HOW A FEATURE IS ASSIGNED TO AN ISOFORM
#
# Both sides of a switch are translated (15) and aligned to the UniProt canonical sequence by
# Smith-Waterman. Walking the alignment gives, for each canonical residue, whether that
# residue is actually present in that isoform. A feature is then scored by the FRACTION of its
# residues present, not by whether its start and end fall inside some interval -- an internal
# exon skip can gut the middle of a domain while leaving both boundaries intact, and an
# endpoint test would call that "retained".
#
# THE FRAMING IS DIRECTIONAL, AND THAT IS THE POINT
#
# Results are expressed as what happens ON SWITCHING: the isoform GAINING usage (dIF > 0) and
# the one LOSING it are identified per pair, so a row reads "the Death domain is lost when
# this switch runs", not "isoform A lacks the Death domain". Which isoform is "candidate" is
# an artefact of how 13 built its list; which one is rising is biology.
#
# WHAT THIS IS NOT
#
# UniProt features describe the CANONICAL protein. A region absent from the canonical -- novel
# exons, extended termini, anything a PacBio-only transcript contributes -- has no annotation
# to lose or gain and will simply not appear. Absence of a listed feature change is therefore
# not evidence that a switch is functionally silent, and `aligned_pct_canonical` is reported
# per isoform so a poorly-matching isoform is visible rather than silently under-annotated.

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
has_json <- requireNamespace("jsonlite", quietly = TRUE)
has_aln <- requireNamespace("pwalign", quietly = TRUE) ||
  requireNamespace("Biostrings", quietly = TRUE)
if (!has_json) stop("jsonlite is required for 16_functional_sites.R")
if (!requireNamespace("Biostrings", quietly = TRUE)) {
  stop("Biostrings is required for 16_functional_sites.R")
}
# pairwiseAlignment moved from Biostrings to pwalign in Bioconductor 3.19; accept either.
align_fn <- if (requireNamespace("pwalign", quietly = TRUE)) {
  list(pa = pwalign::pairwiseAlignment, ap = pwalign::alignedPattern,
       as_ = pwalign::alignedSubject, sub = pwalign::subject, mat = "BLOSUM62")
} else {
  list(pa = Biostrings::pairwiseAlignment, ap = Biostrings::alignedPattern,
       as_ = Biostrings::alignedSubject, sub = Biostrings::subject, mat = "BLOSUM62")
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
results_seq <- resolve_path("results/sequences", root = root)
cache_dir <- ensure_dir(resolve_path("data/external/uniprot", root = root))

pairs_f <- file.path(results_tables, "candidate_switch_pairs.rds")
prot_f <- file.path(results_seq, "candidate_proteins.faa")
if (!file.exists(pairs_f)) {
  stop("Missing candidate_switch_pairs.rds. Run: Rscript scripts/15_isoform_characterization.R")
}
pairs <- readRDS(pairs_f)
if (!file.exists(prot_f)) {
  message("No translated proteins at results/sequences/candidate_proteins.faa; ",
          "re-run 15 with the transcript FASTA available. Nothing written.")
  message("16_functional_sites.R: done (skipped)")
  quit(save = "no", status = 0)
}
prot <- Biostrings::readAAStringSet(prot_f)
names(prot) <- sub("\\s.*$", "", names(prot))

genes <- sort(unique(pairs$gene_name))
message("16_functional_sites.R: mapping switch consequences onto UniProt features for ",
        length(genes), " genes")

# ---- UniProt, cached ------------------------------------------------------------------------
# Cached per gene and committed, so the script is reproducible and offline after the first
# run. A stale cache is a real risk with a live database, so the fetch date is recorded in
# each file and reported in the output table.
UNIPROT_FIELDS <- paste(
  c("accession", "id", "protein_name", "sequence", "gene_primary", "gene_names", "cc_function", "cc_subcellular_location",
    "ft_domain", "ft_region", "ft_motif", "ft_binding", "ft_act_site", "ft_site",
    "ft_transmem", "ft_signal", "ft_coiled", "ft_zn_fing", "ft_dna_bind", "ft_mod_res",
    "ft_lipid", "ft_carbohyd", "ft_disulfid", "ft_compbias"),
  collapse = ","
)
fetch_uniprot <- function(gene) {
  f <- file.path(cache_dir, paste0(gene, ".json"))
  if (file.exists(f)) {
    return(tryCatch(jsonlite::fromJSON(f, simplifyVector = FALSE), error = function(e) NULL))
  }
  url <- paste0(
    "https://rest.uniprot.org/uniprotkb/search?query=",
    utils::URLencode(paste0("gene_exact:", gene, " AND organism_id:9606 AND reviewed:true"),
                     reserved = TRUE),
    "&format=json&size=10&fields=", UNIPROT_FIELDS
  )
  ok <- tryCatch({
    utils::download.file(url, f, quiet = TRUE, method = "libcurl")
    TRUE
  }, error = function(e) FALSE, warning = function(w) FALSE)
  if (!ok || !file.exists(f)) {
    message("  ", gene, ": UniProt fetch failed; skipped")
    return(NULL)
  }
  Sys.sleep(0.2)
  tryCatch(jsonlite::fromJSON(f, simplifyVector = FALSE), error = function(e) NULL)
}

# ---- Which canonical residues does an isoform actually contain? -----------------------------
# Smith-Waterman, then walk the alignment. Only columns where BOTH sides carry a residue count
# as present: a gap in the query opposite a canonical residue means that residue is missing
# from the isoform, which is exactly the case a domain-boundary test would miss.
covered_residues <- function(iso_aa, canon_aa) {
  if (is.null(iso_aa) || !nchar(iso_aa) || !nchar(canon_aa)) return(NULL)
  al <- tryCatch(
    align_fn$pa(Biostrings::AAString(iso_aa), Biostrings::AAString(canon_aa),
                type = "local", substitutionMatrix = align_fn$mat,
                gapOpening = 10, gapExtension = 0.5),
    error = function(e) NULL
  )
  if (is.null(al)) return(NULL)
  ps <- strsplit(as.character(align_fn$ap(al)), "")[[1]]
  ss <- strsplit(as.character(align_fn$as_(al)), "")[[1]]
  spos <- as.integer(BiocGenerics::start(align_fn$sub(al)))
  cov <- logical(nchar(canon_aa))
  p <- spos
  for (i in seq_along(ss)) {
    if (ss[i] == "-") next          # insertion in the isoform; no canonical residue consumed
    if (ps[i] != "-") cov[p] <- TRUE
    p <- p + 1L
  }
  cov
}

feature_rows <- function(entry) {
  fts <- entry$features
  if (is.null(fts) || !length(fts)) return(tibble())
  bind_rows(lapply(fts, function(f) {
    tibble(
      feature_type = f$type %||% NA_character_,
      feature_description = (f$description %||% "") ,
      feature_start = as.integer(f$location$start$value %||% NA),
      feature_end = as.integer(f$location$end$value %||% NA)
    )
  })) |>
    filter(!is.na(.data$feature_start), !is.na(.data$feature_end))
}

first_text <- function(x) {
  if (is.null(x) || !length(x)) return(NA_character_)
  as.character(x)[1L]
}
uniprot_comment <- function(entry, type) {
  cs <- entry$comments
  if (is.null(cs)) return(NA_character_)
  for (cm in cs) {
    if (!identical(cm$commentType, type)) next
    if (!is.null(cm$texts) && length(cm$texts)) return(first_text(cm$texts[[1]]$value))
    if (!is.null(cm$subcellularLocations) && length(cm$subcellularLocations)) {
      return(paste(vapply(cm$subcellularLocations, function(s) {
        first_text(s$location$value)
      }, NA_character_), collapse = "; "))
    }
  }
  NA_character_
}

# ---- Walk the pairs -------------------------------------------------------------------------
site_rows <- list()
gene_rows <- list()

for (g in genes) {
  entry_all <- fetch_uniprot(g)
  res <- if (is.null(entry_all)) NULL else entry_all$results
  if (is.null(res) || !length(res)) {
    message("  ", g, ": no reviewed human UniProt entry; skipped")
    next
  }
  # PICK THE ENTRY WHOSE PRIMARY GENE NAME IS ACTUALLY THIS SYMBOL. `gene_exact` also matches
  # SYNONYMS, and a symbol that is another gene's legacy name silently returns the wrong
  # protein: querying LIG1 returns Q96JA1 (LRIG1, 1093 aa), which lists LIG1 as an old
  # alternative name, instead of DNA ligase 1 (919 aa). Every feature assignment downstream of
  # that is meaningless, and nothing about the output looks wrong.
  primary_name <- function(e) {
    g <- e$genes
    if (is.null(g) || !length(g)) return(NA_character_)
    as.character(g[[1]]$geneName$value %||% NA_character_)
  }
  hit <- Filter(function(e) identical(primary_name(e), g), res)
  if (length(hit)) {
    entry <- hit[[1]]
    symbol_match <- "primary"
  } else if (length(res) == 1L) {
    # No entry carries this symbol as its primary name, and exactly one came back. That is the
    # signature of a RENAME (DDX58 -> RIGI, CBWD3 -> ZNG1C) rather than a synonym collision:
    # a collision like LIG1 -> LRIG1 is caught above, because a correctly-primary LIG1 entry
    # also exists and the filter finds it. Accepted, but the rename is recorded in the output
    # rather than silently absorbed.
    entry <- res[[1]]
    symbol_match <- paste0("renamed: UniProt primary symbol is ", primary_name(entry))
    message("  ", g, ": UniProt primary symbol is ", primary_name(entry), " (rename; accepted)")
  } else {
    message("  ", g, ": ", length(res), " UniProt entries, none primary for this symbol ",
            "and no unambiguous rename; skipped")
    next
  }
  canon <- entry$sequence$value
  acc <- entry$primaryAccession %||% NA_character_
  fts <- feature_rows(entry)

  gene_rows[[length(gene_rows) + 1L]] <- tibble(
    gene_name = g, uniprot_accession = acc,
    uniprot_primary_symbol = primary_name(entry), symbol_match = symbol_match,
    protein_name = first_text(entry$proteinDescription$recommendedName$fullName$value),
    canonical_length_aa = nchar(canon),
    n_annotated_features = nrow(fts),
    uniprot_function = uniprot_comment(entry, "FUNCTION"),
    subcellular_location = uniprot_comment(entry, "SUBCELLULAR LOCATION")
  )

  gp <- pairs |> filter(.data$gene_name == g)
  for (i in seq_len(nrow(gp))) {
    a_id <- gp$isoform_id[i]; b_id <- gp$partner_id[i]
    # Direction: name the rising and falling side rather than candidate/partner.
    if (is.na(gp$dIF_candidate[i])) next
    if (gp$dIF_candidate[i] > 0) {
      gain_id <- a_id; lose_id <- b_id
    } else {
      gain_id <- b_id; lose_id <- a_id
    }
    aa_gain <- if (gain_id %in% names(prot)) as.character(prot[[gain_id]]) else NULL
    aa_lose <- if (lose_id %in% names(prot)) as.character(prot[[lose_id]]) else NULL
    cov_gain <- covered_residues(aa_gain, canon)
    cov_lose <- covered_residues(aa_lose, canon)
    # A side with no ORF called has no protein to align, and treating that as "no residues
    # present" makes every feature read as LOST -- which is how KDM6B and CD86 first came
    # back with 36 and 17 "changed" sites while their coding consequence was
    # "unknown (ORF not called)". Missing is not absent. Such a pair is recorded once, with
    # the reason, and contributes no feature rows.
    if (is.null(cov_gain) || is.null(cov_lose)) {
      site_rows[[length(site_rows) + 1L]] <- tibble(
        gene_name = g, uniprot_accession = acc,
        rising_isoform = gain_id, falling_isoform = lose_id,
        feature_type = "(not assessable)",
        feature_description = paste0(
          "no ORF called for the ",
          if (is.null(cov_gain)) "rising" else "falling",
          " isoform, so no protein could be aligned"),
        feature_start = NA_integer_, feature_end = NA_integer_,
        pct_present_rising = NA_real_, pct_present_falling = NA_real_,
        status = "not assessable (no ORF called)"
      )
      next
    }

    if (!nrow(fts)) next
    for (j in seq_len(nrow(fts))) {
      s <- max(1L, fts$feature_start[j]); e <- min(nchar(canon), fts$feature_end[j])
      if (s > e) next
      fg <- mean(cov_gain[s:e]); fl <- mean(cov_lose[s:e])
      status <- dplyr::case_when(
        fg >= 0.9 & fl >= 0.9 ~ "retained through the switch",
        fg < 0.1 & fl >= 0.9 ~ "LOST on switching",
        fg >= 0.9 & fl < 0.1 ~ "GAINED on switching",
        fg < 0.1 & fl < 0.1 ~ "absent from both isoforms",
        fg < fl ~ "partially lost on switching",
        fg > fl ~ "partially gained on switching",
        TRUE ~ "equally partial in both"
      )
      site_rows[[length(site_rows) + 1L]] <- tibble(
        gene_name = g, uniprot_accession = acc,
        rising_isoform = gain_id, falling_isoform = lose_id,
        feature_type = fts$feature_type[j],
        feature_description = fts$feature_description[j],
        feature_start = fts$feature_start[j], feature_end = fts$feature_end[j],
        pct_present_rising = round(100 * fg, 1),
        pct_present_falling = round(100 * fl, 1),
        status = status
      )
    }
    # The extent each isoform covers, not just how much of it -- two isoforms can cover 35% of
    # the canonical and mean completely different things depending on WHICH 35%. The first and
    # last covered residue are what a domain track has to be drawn from.
    rng <- function(v) if (any(v)) range(which(v)) else c(NA_integer_, NA_integer_)
    rg <- rng(cov_gain); rl <- rng(cov_lose)
    site_rows[[length(site_rows) + 1L]] <- tibble(
      gene_name = g, uniprot_accession = acc,
      rising_isoform = gain_id, falling_isoform = lose_id,
      feature_type = "(whole protein)", feature_description = "canonical coverage",
      feature_start = 1L, feature_end = nchar(canon),
      pct_present_rising = round(100 * mean(cov_gain), 1),
      pct_present_falling = round(100 * mean(cov_lose), 1),
      status = "coverage reference",
      rising_first_residue = rg[1], rising_last_residue = rg[2],
      falling_first_residue = rl[1], falling_last_residue = rl[2]
    )
  }
}

sites <- bind_rows(site_rows)
gene_tab <- bind_rows(gene_rows)

if (!nrow(sites)) {
  message("No functional-site rows produced. Nothing written.")
  message("16_functional_sites.R: done")
  quit(save = "no", status = 0)
}

# 13's candidate list contains BOTH sides of a switch for the five both-sides genes, so the
# pairs table carries that switch twice (a vs b and b vs a). Once each row is expressed as
# rising-vs-falling those two are the same statement, so they are collapsed here -- otherwise
# every feature of those genes is counted twice in the tallies.
sites <- sites |>
  distinct(.data$gene_name, .data$rising_isoform, .data$falling_isoform, .data$feature_type,
           .data$feature_description, .data$feature_start, .data$feature_end,
           .keep_all = TRUE) |>
  arrange(.data$gene_name, .data$feature_start)
write_table_pair(sites, results_tables, "candidate_functional_sites", cfg = cfg)
write_table_pair(gene_tab, results_tables, "candidate_uniprot_annotation", cfg = cfg)

# ---- What changes, per gene ------------------------------------------------------------------
changed <- sites |>
  filter(.data$status %in% c("LOST on switching", "GAINED on switching",
                             "partially lost on switching", "partially gained on switching"))
gene_summary <- changed |>
  group_by(.data$gene_name) |>
  summarize(
    n_features_changed = dplyr::n(),
    lost = sum(grepl("^LOST|^partially lost", .data$status)),
    gained = sum(grepl("^GAINED|^partially gained", .data$status)),
    features = paste(unique(paste0(.data$feature_type,
                                   ifelse(nzchar(.data$feature_description),
                                          paste0(" (", .data$feature_description, ")"), ""))),
                     collapse = "; "),
    .groups = "drop"
  ) |>
  arrange(desc(.data$n_features_changed))
write_table_pair(gene_summary, results_tables, "candidate_functional_site_summary", cfg = cfg)

message("\n---- feature status across ", nrow(sites), " gene/feature/pair rows ----")
print(as.data.frame(sites |> count(.data$status, sort = TRUE)), row.names = FALSE)
message("\n---- genes whose switch changes an annotated site ----")
print(as.data.frame(gene_summary |> mutate(features = substr(.data$features, 1, 90))),
      row.names = FALSE)

write_run_manifest(
  script = "16_functional_sites.R", cfg = cfg, root = root,
  extra = list(uniprot_cache = "data/external/uniprot", n_genes = length(genes))
)
message("\n16_functional_sites.R: done")
