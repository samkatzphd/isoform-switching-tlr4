#!/usr/bin/env Rscript
# 26: What structurally distinguishes each surviving switch -- and therefore which assay can
#     resolve it at all.
#
# Run (after 23, 25): Rscript scripts/26_ht_switch_anatomy.R
# Needs gff3/H_T.gff3 (LOCAL, gitignored). Warns and exits 0 without it.
#
# WHY THIS DECIDES THE ASSAY, NOT JUST THE DESCRIPTION
#
# Script 13's header records the failure mode that keeps recurring: the external review
# screened 7 Q2 candidates and eliminated 6, five of them because the candidate isoform was
# fully NESTED inside a longer isoform of the same gene. That is not bad luck. ~22% of
# switching events here are pure promoter switches and much of the rest is 5' truncation, and
# an isoform that differs only by where it STARTS is contained inside its longer partner by
# construction. There is no candidate-only stretch to put a primer in.
#
# The previous framing treated nesting as a disqualification. That was the wrong conclusion:
# it disqualifies ONE assay. A 5'-end difference is exactly what 5' RACE and CAGE measure
# directly, and a fully nested isoform is still readable by a targeted long-read amplicon
# that sequences the whole molecule. So nesting is recorded here as a ROUTING signal.
#
# WHAT IS CLASSIFIED
#
# Per candidate/partner pair, from exon coordinates and the junction sets they imply:
#   alt_first_exon      terminal 5' exons do not overlap -> promoter / TSS choice
#   alt_last_exon       terminal 3' exons do not overlap -> polyA / 3'UTR choice
#   exon_skipping       an exon of one sits inside an intron of the other
#   intron_retention    an exon of one spans a whole intron of the other
#   alt_donor           two junctions share an acceptor but differ at the donor
#   alt_acceptor        two junctions share a donor but differ at the acceptor
#   mutually_exclusive  both have an exon in the same window, non-overlapping
# plus nesting against every sibling isoform, and the longest candidate-only exonic run
# (what a discriminating qPCR primer would need), against assayability.min_unique_bp.
#
# Strand matters: "first exon" is the lowest coordinate on + and the highest on -.

bt <- c("utils/bootstrap.R", file.path("..", "utils", "bootstrap.R"))
b_file <- if (any(f <- vapply(bt, file.exists, NA))) { bt[which(f)[1L]] } else { NA_character_ }
if (is.na(b_file)) stop("Could not find utils/bootstrap.R. Run from project root.")
source(b_file, local = FALSE, chdir = FALSE)

suppressPackageStartupMessages({ library(dplyr); library(tibble) })

root <- if (exists("PROJECT_ROOT", inherits = FALSE)) get("PROJECT_ROOT", inherits = FALSE) else find_project_root()
args_r <- commandArgs(trailingOnly = TRUE)
cfg <- load_yaml_config(if (length(args_r)) args_r[[1L]] else "config/config.yml", root = root)
paths <- cfg$paths %||% list()
results_tables <- ensure_dir(resolve_path(paths$results_tables %||% "results/tables", root = root))
results_figures <- ensure_dir(resolve_path(paths$results_figures %||% "results/figures", root = root))
fig_dir <- ensure_dir(file.path(results_figures, "assays"))
asy <- cfg$assayability %||% list()
MIN_UNIQUE_BP <- as.integer(asy$min_unique_bp %||% 60)

GFF <- file.path(root, "gff3", "H_T.gff3")
if (!file.exists(GFF)) {
  message("gff3/H_T.gff3 not found -- warn and skip, as 13 does."); quit(save = "no", status = 0)
}

scr_f <- file.path(results_tables, "followup_shortlist_fusion_screened.rds")
if (!file.exists(scr_f)) stop("Run scripts/25_ht_locus_fusion_screen.R first.")
scr <- readRDS(scr_f) |> filter(.data$passes_screen)
message("Surviving candidates to classify: ", nrow(scr))

# ---- exons for every isoform of the relevant loci -----------------------------------------
tf_ids <- tempfile(); writeLines(sort(unique(scr$gene_id)), tf_ids)
tf_ex <- tempfile(fileext = ".tsv")
cmd <- sprintf(
  "awk -F'\\t' 'NR==FNR{w[$1]=1;next} $3==\"exon\"{gid=\"\";tid=\"\";n=split($9,a,\";\");for(i=1;i<=n;i++){if(a[i]~/^gene_id=/)gid=substr(a[i],9);else if(a[i]~/^Parent=/)tid=substr(a[i],8)} if(gid in w)print gid\"\\t\"tid\"\\t\"$1\"\\t\"$4\"\\t\"$5\"\\t\"$7}' %s %s > %s",
  shQuote(tf_ids), shQuote(GFF), shQuote(tf_ex))
message("  extracting exons ..."); stopifnot(system(cmd) == 0L)
ex <- utils::read.delim(tf_ex, header = FALSE, stringsAsFactors = FALSE,
  col.names = c("gene_id", "transcript_id", "seqid", "start", "end", "strand"))
message("  ", nrow(ex), " exons across ", dplyr::n_distinct(ex$transcript_id), " transcripts")

exon_tbl <- function(tid) {
  d <- ex[ex$transcript_id == tid, , drop = FALSE]
  d[order(d$start), c("start", "end", "strand", "seqid")]
}
#' Introns implied by an exon table (sorted by coordinate).
introns_of <- function(d) {
  if (nrow(d) < 2L) return(data.frame(start = integer(0), end = integer(0)))
  data.frame(start = utils::head(d$end, -1L) + 1L, end = d$start[-1L] - 1L)
}
#' Junctions as donor/acceptor coordinate pairs.
junctions_of <- function(d) {
  i <- introns_of(d)
  if (!nrow(i)) return(character(0))
  paste0(i$start, "-", i$end)
}
overlaps <- function(a1, a2, b1, b2) a1 <= b2 & a2 >= b1

#' Classify one candidate/partner pair.
classify_pair <- function(cand, part) {
  A <- exon_tbl(cand); B <- exon_tbl(part)
  if (!nrow(A) || !nrow(B)) return(NULL)
  strand <- A$strand[1]
  # 5'/3' terminal exons, strand aware
  a5 <- if (strand == "+") A[1, ] else A[nrow(A), ]
  b5 <- if (strand == "+") B[1, ] else B[nrow(B), ]
  a3 <- if (strand == "+") A[nrow(A), ] else A[1, ]
  b3 <- if (strand == "+") B[nrow(B), ] else B[1, ]

  jA <- junctions_of(A); jB <- junctions_of(B)
  iA <- introns_of(A);   iB <- introns_of(B)

  # exon of one wholly inside an intron of the other
  skip <- any(vapply(seq_len(nrow(B)), function(j)
    any(iA$start <= B$start[j] & iA$end >= B$end[j]), logical(1))) ||
          any(vapply(seq_len(nrow(A)), function(j)
    any(iB$start <= A$start[j] & iB$end >= A$end[j]), logical(1)))
  # exon of one spanning a whole intron of the other
  ir <- any(vapply(seq_len(nrow(A)), function(j)
    any(A$start[j] <= iB$start & A$end[j] >= iB$end), logical(1))) ||
        any(vapply(seq_len(nrow(B)), function(j)
    any(B$start[j] <= iA$start & B$end[j] >= iA$end), logical(1)))

  # junctions sharing one end but not the other
  da <- function(x) as.integer(sub("-.*$", "", x))
  ac <- function(x) as.integer(sub("^.*-", "", x))
  onlyA <- setdiff(jA, jB); onlyB <- setdiff(jB, jA)
  alt_donor <- FALSE; alt_acceptor <- FALSE
  if (length(onlyA) && length(onlyB)) {
    alt_donor    <- any(outer(ac(onlyA), ac(onlyB), "==") & outer(da(onlyA), da(onlyB), "!="))
    alt_acceptor <- any(outer(da(onlyA), da(onlyB), "==") & outer(ac(onlyA), ac(onlyB), "!="))
  }
  tibble(
    isoform_id = cand, partner = part, strand = strand,
    n_exons_candidate = nrow(A), n_exons_partner = nrow(B),
    alt_first_exon = !overlaps(a5$start, a5$end, b5$start, b5$end),
    alt_last_exon  = !overlaps(a3$start, a3$end, b3$start, b3$end),
    exon_skipping = skip, intron_retention = ir,
    alt_donor = alt_donor, alt_acceptor = alt_acceptor,
    n_junctions_candidate_only = length(onlyA),
    n_junctions_partner_only = length(onlyB),
    tss_distance_bp = abs(if (strand == "+") a5$start - b5$start else a5$end - b5$end),
    tts_distance_bp = abs(if (strand == "+") a3$end - b3$end else a3$start - b3$start)
  )
}

#' Is the candidate wholly contained in a sibling, and what is its longest unique run?
nesting_and_unique <- function(cand, gene) {
  A <- exon_tbl(cand)
  sibs <- setdiff(unique(ex$transcript_id[ex$gene_id == gene]), cand)
  host <- NA_character_
  for (s in sibs) {
    S <- exon_tbl(s)
    covered <- all(vapply(seq_len(nrow(A)), function(j)
      any(S$start <= A$start[j] & S$end >= A$end[j]), logical(1)))
    if (covered) { host <- s; break }
  }
  # candidate-only exonic sequence: subtract the union of all siblings' exons
  keep <- data.frame(start = A$start, end = A$end)
  for (s in sibs) {
    S <- exon_tbl(s); out <- list()
    for (j in seq_len(nrow(keep))) {
      seg <- c(keep$start[j], keep$end[j])
      pieces <- list(seg)
      for (k in seq_len(nrow(S))) {
        np <- list()
        for (pz in pieces) {
          if (S$end[k] < pz[1] || S$start[k] > pz[2]) { np[[length(np) + 1L]] <- pz; next }
          if (S$start[k] > pz[1]) np[[length(np) + 1L]] <- c(pz[1], S$start[k] - 1L)
          if (S$end[k] < pz[2])   np[[length(np) + 1L]] <- c(S$end[k] + 1L, pz[2])
        }
        pieces <- np
      }
      for (pz in pieces) out[[length(out) + 1L]] <- pz
    }
    keep <- if (length(out)) {
      data.frame(start = vapply(out, `[`, numeric(1), 1), end = vapply(out, `[`, numeric(1), 2))
    } else data.frame(start = integer(0), end = integer(0))
    if (!nrow(keep)) break
  }
  longest <- if (nrow(keep)) max(keep$end - keep$start + 1L) else 0L
  region <- if (nrow(keep)) {
    i <- which.max(keep$end - keep$start + 1L)
    paste0(A$seqid[1], ":", keep$start[i], "-", keep$end[i])
  } else NA_character_
  tibble(isoform_id = cand, nested_in = host, is_nested = !is.na(host),
         longest_unique_bp = as.integer(longest), unique_region = region,
         n_sibling_isoforms = length(sibs))
}

message("  classifying ", nrow(scr), " pairs ...")
anat <- bind_rows(lapply(seq_len(nrow(scr)), function(k)
  classify_pair(scr$isoform_id[k], scr$partner[k])))
nest <- bind_rows(lapply(seq_len(nrow(scr)), function(k)
  nesting_and_unique(scr$isoform_id[k], scr$gene_id[k])))

out <- scr |>
  select("gene_name", "gene_id", "isoform_id", "discovered_in", "set", "tier",
         "own_dIF", "other_dIF", "shared_exonic_bp") |>
  left_join(anat, by = "isoform_id") |>
  left_join(nest, by = "isoform_id") |>
  mutate(
    # the single structural headline, in the order that determines assay choice
    primary_event = dplyr::case_when(
      .data$alt_first_exon & .data$alt_last_exon ~ "both ends differ",
      .data$alt_first_exon ~ "alternative first exon (TSS)",
      .data$alt_last_exon ~ "alternative last exon (3' end)",
      .data$intron_retention ~ "intron retention",
      .data$exon_skipping ~ "exon skipping",
      .data$alt_donor | .data$alt_acceptor ~ "alternative splice site",
      TRUE ~ "other / unresolved"
    ),
    qpcr_addressable = !.data$is_nested & .data$longest_unique_bp >= MIN_UNIQUE_BP
  ) |>
  arrange(dplyr::desc(abs(.data$own_dIF)))
write_table_pair(out, results_tables, "ht_switch_anatomy", cfg = cfg)

message("\n---- structural composition of the ", nrow(out), " survivors ----")
print(as.data.frame(out |> count(.data$primary_event, name = "candidates") |>
                      arrange(dplyr::desc(.data$candidates))), row.names = FALSE)
message("\n  nested inside a sibling:        ", sum(out$is_nested %in% TRUE),
        "  <- cannot be addressed by a unique-region qPCR primer")
message("  longest unique run >= ", MIN_UNIQUE_BP, " bp:    ",
        sum(out$longest_unique_bp >= MIN_UNIQUE_BP, na.rm = TRUE))
message("  qPCR-addressable:              ", sum(out$qpcr_addressable %in% TRUE),
        " of ", nrow(out))
message("\n  median TSS distance between the pair: ",
        stats::median(out$tss_distance_bp, na.rm = TRUE), " bp")

message("\n---- qPCR-addressable candidates ----")
print(as.data.frame(out |> filter(.data$qpcr_addressable) |>
  transmute(gene = .data$gene_name, event = .data$primary_event,
            dIF = round(.data$own_dIF, 3), unique_bp = .data$longest_unique_bp,
            region = .data$unique_region)), row.names = FALSE)

if (requireNamespace("ggplot2", quietly = TRUE)) {
  library(ggplot2)
  p <- out |>
    count(.data$primary_event, .data$qpcr_addressable, name = "n") |>
    mutate(primary_event = stats::reorder(.data$primary_event, .data$n, FUN = sum),
           route = ifelse(.data$qpcr_addressable, "qPCR-addressable",
                          "needs RACE / long-read")) |>
    ggplot(aes(x = .data$n, y = .data$primary_event, fill = .data$route)) +
    geom_col(width = 0.68) +
    scale_fill_manual(values = c("qPCR-addressable" = "#1F6F54",
                                 "needs RACE / long-read" = "#F2B134"), name = NULL) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.08))) +
    labs(title = "What distinguishes each surviving switch, and what can read it",
         subtitle = paste0("A nested or 5'-only difference has no candidate-only stretch to ",
                           "prime into, so qPCR cannot resolve it;\n5' RACE and targeted ",
                           "long-read cDNA can."),
         x = "candidates", y = NULL) +
    theme_bw(base_size = 11) +
    theme(plot.title = element_text(face = "bold", size = 12),
          plot.subtitle = element_text(size = 9, colour = "grey30"),
          legend.position = "top", panel.grid.minor = element_blank())
  ggsave(file.path(fig_dir, "switch_anatomy.png"), plot = p,
         width = 9.5, height = 5.2, dpi = 200, bg = "white")
  message("\nFigure: ", file.path(fig_dir, "switch_anatomy.png"))
}

write_run_manifest(script = "26_ht_switch_anatomy.R", cfg = cfg, root = root)
message("\n26_ht_switch_anatomy.R: done")
