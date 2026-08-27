#!/usr/bin/env Rscript
# 13: Can these candidate isoforms actually be measured? Nesting + unique-sequence screen.
#
# Run (after 11): Rscript scripts/13_assayability_screen.R
# Needs gff3/U_T.gff3 and fa/transcripts_U-T.fa, which are LOCAL but gitignored (3.9 GB).
# Warns and exits 0 without them. The external drive is NOT required.
#
# WHY THIS RUNS BEFORE ANY PRIMER IS ORDERED
#
# The external review applied this screen to the 7 Q2 candidate isoforms and eliminated 6 of
# them: 5 were fully NESTED inside a longer isoform of the same gene, and a 6th had no unique
# stretch long enough to design against. Only IRAK3 TCONS_00065141 survived. The Q1 lists from
# script 11 have never been screened, so they are biological shortlists, not assay-ready ones.
#
# The structural reason this keeps happening is in the event-structure result: ~22% of
# switching events are PURE promoter switches and much of the rest involves 5' truncation
# (event_structure_summary_gff3.csv). An isoform that differs only by where it starts is,
# by construction, contained inside its longer partner -- so nesting is the expected case
# here, not an unlucky one.
#
# THE TWO TESTS
#
# 1. NESTING. Candidate isoform X is nested if some other isoform Y of the same gene contains
#    every base X has AND every junction X has. If so, no primer pair can distinguish X from
#    Y: any amplicon inside X is also amplifiable from Y. Junctions matter as well as bases --
#    an isoform whose exons all fall inside Y's but which splices differently is NOT nested,
#    because a junction-spanning primer can discriminate it.
#
# 2. UNIQUE SEQUENCE. Even when not nested, an isoform can be tiled over by the union of the
#    others. Subtract every other isoform's exonic footprint from X's and ask for a remaining
#    contiguous run of at least `min_unique_bp` (60 by default, matching the review). Below
#    that there is nowhere to put a primer pair.
#
# WHAT THIS SCRIPT DOES NOT DO. It does not design primers or run the empirical
# specificity match against all 315,349 transcripts. That is the step AFTER this one, and it
# only makes sense for isoforms that survive here. Survivors are written with their unique
# interval so primer design can start from a coordinate.

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
has_ggplot <- requireNamespace("ggplot2", quietly = TRUE)
if (has_ggplot) suppressPackageStartupMessages(library(ggplot2))

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
fig_dir <- ensure_dir(file.path(results_figures, "assayability"))

asy <- cfg$assayability %||% list()
gff_rel <- asy$gff3_path %||% "gff3/U_T.gff3"
min_unique_bp <- as.integer(asy$min_unique_bp %||% 60)
cand_min_abs_dif <- as.numeric(asy$candidate_min_abs_dif %||% 0.15)
gff_path <- resolve_path(gff_rel, root = root)

if (!file.exists(gff_path)) {
  message("13_assayability_screen.R: ", gff_path, " not found.")
  message("  Copy the merged GFF3 in from the external drive and re-run. Nothing written.")
  quit(save = "no", status = 0L)
}
if (has_ggplot) theme_set(theme_bw(base_size = 11))

# ---- Candidate genes -> candidate isoforms ------------------------------------------------
lists <- list(
  response = file.path(results_tables, "q1_response_candidates.csv"),
  cryptic  = file.path(results_tables, "q1_cryptic_candidates.csv")
)
missing <- lists[!file.exists(unlist(lists))]
if (length(missing)) {
  stop("Missing candidate lists. Run: Rscript scripts/11_q1_candidates.R\n  ",
       paste(unlist(missing), collapse = "\n  "))
}
cand_genes <- bind_rows(lapply(names(lists), function(nm) {
  utils::read.csv(lists[[nm]], stringsAsFactors = FALSE) |>
    transmute(gene_name = .data$gene_name, list = nm, dif_rank = .data$dif_rank)
})) |>
  group_by(.data$gene_name) |>
  summarize(lists = paste(sort(unique(.data$list)), collapse = "+"),
            dif_rank = max(.data$dif_rank, na.rm = TRUE), .groups = "drop")

message("Candidate genes across both Q1 lists: ", nrow(cand_genes))

# The isoform set per gene: the switching isoforms plus their partners. You cannot validate a
# switch by measuring only one side of it, which is why the rising partner is admitted on
# effect size even when it misses the isoform-level q cutoff (the review's rule).
ctx <- load_scoring_table(processed_dir, "T_UT", quiet = TRUE)
if (is.null(ctx)) stop("No T_UT scoring table in ", processed_dir, "; run 01_load_data.R.")
scored <- score_isoforms(ctx, cfg, "T_UT", "T_UT")

cand_iso <- scored |>
  filter(.data$gene_name %in% cand_genes$gene_name) |>
  mutate(abs_dIF = abs(.data$dIF_n)) |>
  filter(.data$is_switching %in% TRUE | .data$abs_dIF > cand_min_abs_dif) |>
  filter(is.finite(.data$abs_dIF)) |>
  select("isoform_id", "gene_id", "gene_name", "dIF_n", "abs_dIF", "is_switching",
         "class_code", "oId", "is_novel_pacbio") |>
  left_join(cand_genes, by = "gene_name")

message("Candidate isoforms to screen: ", nrow(cand_iso),
        " across ", length(unique(cand_iso$gene_id)), " genes")

# ---- Exon structures for those genes only -------------------------------------------------
# 2.1M exon records in the GFF3; pull just the genes of interest with grep before parsing.
want_genes <- unique(cand_iso$gene_id)
message("Extracting exon records from ", basename(gff_path), " ...")
pat <- paste0("gene_id=(", paste(want_genes, collapse = "|"), ");")
tmp <- tempfile(fileext = ".tsv")
on.exit(unlink(tmp), add = TRUE)
ok <- system2("grep", c("-E", shQuote(pat), shQuote(gff_path)), stdout = tmp, stderr = FALSE)
if (!file.exists(tmp) || file.info(tmp)$size == 0) {
  stop("No GFF3 records matched the candidate gene ids.")
}
gff <- utils::read.delim(tmp, header = FALSE, comment.char = "#", quote = "",
                         stringsAsFactors = FALSE)
names(gff)[c(1, 3, 4, 5, 7, 9)] <- c("chr", "feature", "start", "end", "strand", "attr")
ex <- gff |>
  filter(.data$feature == "exon") |>
  mutate(
    isoform_id = sub('^.*Parent=([^;]+).*$', "\\1", .data$attr),
    gene_id = sub('^.*gene_id=([^;]+).*$', "\\1", .data$attr)
  ) |>
  filter(.data$gene_id %in% want_genes) |>
  select("gene_id", "isoform_id", "chr", "start", "end", "strand")
message("  ", nrow(ex), " exon records for ", length(unique(ex$isoform_id)), " isoforms")

# ---- Interval helpers ---------------------------------------------------------------------
# Base-R interval algebra on small per-gene sets; no GenomicRanges dependency.
norm_iv <- function(s, e) {
  if (!length(s)) return(matrix(numeric(0), ncol = 2))
  o <- order(s); s <- s[o]; e <- e[o]
  keep_s <- s[1]; keep_e <- e[1]; out <- list()
  for (i in seq_along(s)[-1]) {
    if (s[i] <= keep_e + 1) {
      keep_e <- max(keep_e, e[i])
    } else {
      out[[length(out) + 1L]] <- c(keep_s, keep_e); keep_s <- s[i]; keep_e <- e[i]
    }
  }
  out[[length(out) + 1L]] <- c(keep_s, keep_e)
  do.call(rbind, out)
}
iv_subtract <- function(a, b) {
  # a, b: 2-col matrices of closed intervals. Returns a minus b.
  if (!nrow(a)) return(a)
  if (!nrow(b)) return(a)
  out <- list()
  for (i in seq_len(nrow(a))) {
    seg <- list(c(a[i, 1], a[i, 2]))
    for (j in seq_len(nrow(b))) {
      nxt <- list()
      for (s in seg) {
        if (b[j, 2] < s[1] || b[j, 1] > s[2]) { nxt[[length(nxt) + 1L]] <- s; next }
        if (b[j, 1] > s[1]) nxt[[length(nxt) + 1L]] <- c(s[1], b[j, 1] - 1)
        if (b[j, 2] < s[2]) nxt[[length(nxt) + 1L]] <- c(b[j, 2] + 1, s[2])
      }
      seg <- nxt
      if (!length(seg)) break
    }
    if (length(seg)) out <- c(out, seg)
  }
  if (!length(out)) return(matrix(numeric(0), ncol = 2))
  do.call(rbind, out)
}
iv_covers <- function(outer, inner) {
  # does `outer` contain every base of `inner`?
  nrow(iv_subtract(inner, outer)) == 0
}
junctions_of <- function(d) {
  if (nrow(d) < 2) return(character(0))
  d <- d[order(d$start), ]
  paste0(d$end[-nrow(d)], "-", d$start[-1])
}

# ---- Screen -------------------------------------------------------------------------------
by_gene <- split(ex, ex$gene_id)
screen_one <- function(iso_id, gene_id) {
  g <- by_gene[[gene_id]]
  if (is.null(g)) {
    return(tibble(nested_in = NA_character_, is_nested = NA,
                  longest_unique_bp = NA_real_, unique_region = NA_character_,
                  n_other_isoforms = NA_integer_))
  }
  me <- g[g$isoform_id == iso_id, , drop = FALSE]
  if (!nrow(me)) {
    return(tibble(nested_in = NA_character_, is_nested = NA,
                  longest_unique_bp = NA_real_, unique_region = NA_character_,
                  n_other_isoforms = NA_integer_))
  }
  others <- unique(g$isoform_id[g$isoform_id != iso_id])
  my_iv <- norm_iv(me$start, me$end)
  my_j <- junctions_of(me)

  # 1. nesting
  host <- NA_character_
  for (o in others) {
    od <- g[g$isoform_id == o, , drop = FALSE]
    if (iv_covers(norm_iv(od$start, od$end), my_iv) && all(my_j %in% junctions_of(od))) {
      host <- o; break
    }
  }

  # 2. unique sequence against the UNION of all other isoforms
  if (length(others)) {
    od <- g[g$isoform_id %in% others, , drop = FALSE]
    uniq <- iv_subtract(my_iv, norm_iv(od$start, od$end))
  } else {
    uniq <- my_iv
  }
  if (nrow(uniq)) {
    w <- uniq[, 2] - uniq[, 1] + 1
    k <- which.max(w)
    longest <- w[k]
    region <- paste0(me$chr[1], ":", uniq[k, 1], "-", uniq[k, 2])
  } else {
    longest <- 0; region <- NA_character_
  }

  tibble(nested_in = host, is_nested = !is.na(host),
         longest_unique_bp = longest, unique_region = region,
         n_other_isoforms = length(others))
}

message("Screening ...")
res <- cand_iso |>
  rowwise() |>
  mutate(scr = list(screen_one(.data$isoform_id, .data$gene_id))) |>
  ungroup() |>
  tidyr::unnest("scr") |>
  mutate(
    assayable = !is.na(.data$is_nested) & !.data$is_nested &
      .data$longest_unique_bp >= min_unique_bp,
    verdict = case_when(
      is.na(.data$is_nested) ~ "not in GFF3",
      .data$is_nested ~ "nested -- no primer pair can discriminate",
      .data$longest_unique_bp < min_unique_bp ~
        paste0("tiled over -- longest unique run ", .data$longest_unique_bp, " bp"),
      TRUE ~ "designable"
    )
  ) |>
  arrange(desc(.data$assayable), desc(.data$dif_rank), .data$gene_name)

write_table_pair(res, results_tables, "q1_assayability", cfg = cfg)

# ---- Class code vs assayability -----------------------------------------------------------
# gffcompare class codes, and why they belong in this table:
#   =  matches a reference transcript end to end
#   c  CONTAINED in a reference transcript -- i.e. nested, structurally
#   j  novel junction combination, shares at least one junction
#
# `c` is the one to watch. It is the expected code for a genuine alternative internal
# promoter, but it is EQUALLY the code a 5'-incomplete assembly lands in -- a PacBio read
# that never reached the 5' end, or degraded input. Without cap selection those are produced
# in quantity. So a `c`-class candidate is either a real truncated isoform or an artefact of
# assembly, and this data cannot tell which.
#
# The correspondence with the screen is total and not a coincidence: "contained in another
# transcript" is nearly the definition of nested, so `c` predicts undesignable without
# running any of the interval algebra above.
cc <- res |>
  mutate(cls = ifelse(.data$assayable, "designable", "not designable")) |>
  count(.data$class_code, .data$cls) |>
  tidyr::pivot_wider(names_from = "cls", values_from = "n", values_fill = 0)
write_table_pair(cc, results_tables, "q1_assayability_by_class_code", cfg = cfg)
message("\n---- class code vs assayability ----")
print(as.data.frame(cc))
n_c <- sum(res$class_code == "c", na.rm = TRUE)
if (n_c) {
  message("NOTE: ", n_c, " candidate isoforms are class `c` (contained). All are nested by ",
          "construction.\n  These are either real internal-promoter isoforms or ",
          "5'-incomplete assemblies;\n  distinguishing them needs cap-selection status for ",
          "the PacBio library.")
}

# ---- Gene-level verdict -------------------------------------------------------------------
# A gene is assayable if at least one side of its switch can be measured. Being able to
# measure both sides is better, so it is counted separately.
gene_res <- res |>
  group_by(.data$gene_name, .data$lists) |>
  summarize(
    n_candidate_iso = dplyr::n(),
    n_designable = sum(.data$assayable, na.rm = TRUE),
    n_nested = sum(.data$is_nested %in% TRUE),
    max_unique_bp = suppressWarnings(max(.data$longest_unique_bp, na.rm = TRUE)),
    class_codes = paste(sort(unique(.data$class_code)), collapse = ","),
    n_contained = sum(.data$class_code == "c", na.rm = TRUE),
    .groups = "drop"
  ) |>
  mutate(
    max_unique_bp = ifelse(is.finite(.data$max_unique_bp), .data$max_unique_bp, NA_real_),
    gene_verdict = case_when(
      .data$n_designable >= 2 ~ "both sides designable",
      .data$n_designable == 1 ~ "one side designable",
      TRUE ~ "no designable isoform"
    )
  ) |>
  arrange(desc(.data$n_designable), desc(.data$max_unique_bp))

write_table_pair(gene_res, results_tables, "q1_assayability_by_gene", cfg = cfg)

summary_tbl <- tibble(
  metric = c("candidate isoforms screened", "  designable",
             "  nested (undiscriminable)", "  tiled over (<min unique bp)",
             "  absent from GFF3",
             "candidate genes", "  both sides designable", "  one side designable",
             "  no designable isoform"),
  n = c(nrow(res), sum(res$assayable, na.rm = TRUE),
        sum(res$is_nested %in% TRUE),
        sum(!res$is_nested %in% TRUE & res$longest_unique_bp < min_unique_bp, na.rm = TRUE),
        sum(is.na(res$is_nested)),
        nrow(gene_res), sum(gene_res$gene_verdict == "both sides designable"),
        sum(gene_res$gene_verdict == "one side designable"),
        sum(gene_res$gene_verdict == "no designable isoform"))
)
write_table_pair(summary_tbl, results_tables, "q1_assayability_summary", cfg = cfg)
print(as.data.frame(summary_tbl))

message("\nGenes with at least one designable isoform:")
print(as.data.frame(gene_res |> filter(.data$n_designable > 0) |>
                      select("gene_name", "lists", "n_designable", "max_unique_bp") |>
                      head(25)))

# ---- Figure --------------------------------------------------------------------------------
if (has_ggplot && nrow(res)) {
  d <- res |>
    mutate(cls = case_when(
      is.na(.data$is_nested) ~ "not in GFF3",
      .data$is_nested ~ "nested",
      .data$longest_unique_bp < min_unique_bp ~ "tiled over",
      TRUE ~ "designable"
    )) |>
    count(.data$lists, .data$cls)
  p <- ggplot(d, aes(x = .data$lists, y = .data$n, fill = .data$cls)) +
    geom_col(width = .7) +
    geom_text(aes(label = .data$n), position = position_stack(vjust = .5), size = 3) +
    scale_fill_manual(values = c("designable" = "#1b9e77", "nested" = "#d95f02",
                                 "tiled over" = "#e6ab02", "not in GFF3" = "grey70"),
                      name = NULL) +
    labs(
      title = "Assayability of the Q1 candidate isoforms",
      subtitle = paste0("Nested isoforms cannot be discriminated by any primer pair; ",
                        "designable = not nested and >= ", min_unique_bp, " bp unique"),
      x = "candidate list", y = "candidate isoforms"
    )
  ggsave(file.path(fig_dir, "fig_q1_assayability.png"), p, width = 8, height = 5, dpi = 200)
  message("Figure: ", fig_dir)
}

write_run_manifest("13_assayability_screen.R", cfg, root,
                   extra = list(assayability = list(min_unique_bp = min_unique_bp,
                                                    gff3 = gff_rel)))
message("13_assayability_screen.R: done")
