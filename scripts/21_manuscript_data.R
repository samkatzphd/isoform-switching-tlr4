#!/usr/bin/env Rscript
# 21: The HT-only manuscript data spine -- one tidy table per figure panel.
#
# Run (after 01, 02d, 03, 19): Rscript scripts/21_manuscript_data.R
#
# WHY THIS EXISTS, AND WHAT IT DELIBERATELY EXCLUDES
#
# This builds the data for a manuscript that stands WITHOUT the UBL5 knockout. Scope is the
# HT reference transcriptome and the H and T libraries only. `T_UT` and `U_UT` are excluded
# even though `T_UT` is the same six wildtype libraries as `T_HT` -- it is quantified against
# the UT annotation, and mixing annotation spaces is what §0b of REVIEW_CHANGES warns about.
# The exclusion is ASSERTED below, not just intended: the script stops if any UT-space label
# reaches a panel table.
#
# THE H_HT FRAMING CONSTRAINT IS LOAD-BEARING (AGENT_CONTEXT fact 4)
#
# `H_HT` is human PRIMARY macrophages, three donors, three days, and it is the only dataset
# `importRdata` batch-corrected (2 surrogate variables, limma::removeBatchEffect). Its
# replicate IF noise is halved relative to T_HT as a result, which is why it shows 2.07% vs
# 0.87% switching. **Its rates and dIF distributions are not comparable to T_HT's.** Every
# panel table that carries both datasets also carries `comparable_across_datasets`, and the
# figure scripts must facet rather than overlay wherever that flag is FALSE. H enters the
# argument as gene-level corroboration in an independent cell system -- not as a second
# measurement of the same quantity.
#
# OUTPUT
#
# results/tables/manuscript/panel_<id>.csv, one per panel, plus manuscript_panel_index.csv
# mapping panel -> figure -> source -> status (real | placeholder). The figure scripts read
# ONLY those CSVs: no .rds, no context objects, no external drive.

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
man_dir <- ensure_dir(file.path(results_tables, "manuscript"))

# ---- scope: HT reference only -----------------------------------------------------------

HT_KEYS <- names(Filter(
  function(d) identical(as.character(d$reference_transcriptome %||% ""), "HT"),
  cfg$datasets %||% list()
))
if (!setequal(HT_KEYS, c("T_HT", "H_HT"))) {
  stop("Expected exactly T_HT and H_HT on the HT reference; got: ", paste(HT_KEYS, collapse = ", "))
}
lab_of <- function(k) sanitize(as.character((cfg$datasets[[k]] %||% list())$label %||% k)[1L])
HT_LABELS <- vapply(HT_KEYS, lab_of, character(1))
message("HT scope: ", paste(HT_LABELS, collapse = ", "))

# A panel table must never carry UT-space information. Scan EVERY character/factor column,
# not just the ones named `dataset*`: the first version of this check looked only at dataset
# columns and missed `called_in` in q1_candidates, whose values ("T_UT only", "both
# annotations") encode the T_HT-vs-T_UT annotation comparison and are therefore UT-dependent
# even though the column name says nothing about a dataset.
assert_ht_only <- function(x, panel) {
  for (col in names(x)) {
    v <- x[[col]]
    if (!is.character(v) && !is.factor(v)) next
    bad <- grep("T_UT|U_UT", unique(as.character(v)), value = TRUE)
    if (length(bad)) {
      stop("Panel ", panel, ": UT-space data leaked into column '", col, "': ",
           paste(utils::head(bad, 5), collapse = ", "))
    }
  }
  invisible(TRUE)
}

index_rows <- list()

#' Panel tables are read back by plain read.csv, so their content must survive the round
#' trip. Non-ASCII characters do not: write.csv(fileEncoding = "UTF-8") truncated a field at
#' a section sign and left an unterminated quote, which silently corrupted the whole index.
#' Fail loudly instead.
assert_ascii <- function(x, panel) {
  for (col in names(x)) {
    v <- x[[col]]
    if (!is.character(v) && !is.factor(v)) next
    v <- as.character(v)
    bad <- v[!is.na(v) & grepl("[^\x01-\x7F]", v)]
    if (length(bad)) {
      stop("Panel ", panel, ": non-ASCII content in column '", col,
           "' will not survive write.csv/read.csv: ", bad[1])
    }
  }
  invisible(TRUE)
}

emit <- function(x, panel, figure, title, source, status = "real", notes = "") {
  assert_ht_only(x, panel)
  x <- as.data.frame(x)
  assert_ascii(x, panel)
  assert_ascii(data.frame(title = title, source = source, notes = notes), panel)
  utils::write.csv(x, file.path(man_dir, paste0("panel_", panel, ".csv")),
                   row.names = FALSE, fileEncoding = "UTF-8", na = "")
  index_rows[[panel]] <<- tibble(
    panel = panel, figure = figure, title = title, source = source,
    status = status, n_rows = nrow(x), notes = notes
  )
  message(sprintf("  %-5s %-10s %-4s n=%-6d %s", panel, figure, status, nrow(x), title))
}

# Keep a single named vector for display labels so every figure agrees.
DS_PRETTY <- c(T_HT_anchor = "THP-1 (T), HT ref", H_HT_supplementary = "Primary MDM (H), HT ref")

# ---- load the two HT context layers once ------------------------------------------------

message("\n==== loading HT context layers ====")
scored <- list()
for (k in HT_KEYS) {
  lab <- lab_of(k)
  src <- load_scoring_table(processed_dir, lab, quiet = TRUE)
  if (is.null(src)) stop("No scoring table for ", lab, " -- run scripts/01_load_data.R first.")
  z <- score_isoforms(src, cfg, dataset_key = k, dataset_label = lab)
  z$detect_tier <- {
    o <- as.character(z$oId); cc <- as.character(z$class_code)
    is_pb <- !is.na(o) & startsWith(o, "PB")
    t <- rep(NA_character_, length(o))
    t[!is_pb & !is.na(o)] <- "reference"
    t[is_pb & cc %in% "="] <- "pb_redundant"
    t[is_pb & cc %in% "c"] <- "pb_contained"
    t[is_pb & cc %in% c("j","k","m","n","o","u","x","i","y","e","p","s")] <- "pb_novel"
    t
  }
  scored[[lab]] <- z
  message(sprintf("  [%-20s] %6d isoforms, %4d switching", lab, nrow(z), sum(z$is_switching)))
}

# =========================================================================================
# FIGURE 1 -- design and sequencing
# =========================================================================================
message("\n==== Figure 1 ====")

# --- 1A: design schematic. PLACEHOLDER. ---
# What IS known is the library structure; what is NOT known is the protocol, and those gaps
# are the standing wet-lab questions (lab meeting update §6).
f1a <- tibble(
  element = c("cell system", "cell system", "LPS", "replication", "replication",
              "PMA differentiation", "LPS dose", "LPS timepoint", "library chemistry",
              "cap verification"),
  value = c("THP-1 monocytic line (T), one clone", "Human primary MDM (H), 3 donors",
            "paired minus / plus per replicate", "T1-T3 matured, treated, extracted separately",
            "H1-H3 three donors on three days",
            "UNKNOWN", "UNKNOWN", "UNKNOWN", "UNKNOWN (polyA vs ribo-depleted)", "UNKNOWN"),
  known = c(TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE)
)
emit(f1a, "f1a", "fig1", "Experimental design schematic",
     # Keep panel-table CONTENT pure ASCII: a literal section sign here was truncated
     # mid-field by write.csv(fileEncoding = "UTF-8"), producing an unterminated quote
     # that made the whole index unreadable by read.csv.
     "docs/AGENT_CONTEXT.md fact 4 + lab meeting update sec 6", status = "placeholder",
     notes = "Needs PMA protocol, LPS dose, timepoint, library chemistry, cap verification.")

# --- 1B: per-library sequencing depth, H and T only ---
cnt_f <- file.path(root, "counts", "H_T_pipe", "RSEM.isoforms.expected_count.all_samples.txt")
if (file.exists(cnt_f)) {
  rc <- utils::read.delim(cnt_f, check.names = FALSE, stringsAsFactors = FALSE)
  samp <- setdiff(names(rc), c("gene_id", "GeneName", "transcript_id"))
  # The HT-reference quantification covers the U libraries too; drop them -- out of scope.
  keep <- samp[grepl("^[HT][0-9]_", samp)]
  message("  library columns: ", length(samp), " total, ", length(keep), " in HT scope")
  f1b <- tibble(
    library = keep,
    cell_system = ifelse(startsWith(keep, "H"), "Primary MDM (H)", "THP-1 (T)"),
    replicate = sub("^([HT][0-9]).*$", "\\1", keep),
    lps = ifelse(grepl("_plus_", keep), "plus", "minus"),
    total_expected_count = vapply(keep, function(s) sum(rc[[s]], na.rm = TRUE), numeric(1)),
    n_isoforms_detected = vapply(keep, function(s) sum(rc[[s]] > 0, na.rm = TRUE), numeric(1)),
    n_genes_detected = vapply(keep, function(s) {
      dplyr::n_distinct(rc$gene_id[rc[[s]] > 0]) }, numeric(1))
  )
  emit(f1b, "f1b", "fig1", "Per-library depth and detection",
       "counts/H_T_pipe/RSEM.isoforms.expected_count.all_samples.txt")
} else {
  emit(tibble(note = "counts/H_T_pipe/RSEM.isoforms.expected_count.all_samples.txt absent"),
       "f1b", "fig1", "Per-library depth and detection", "RSEM expected counts",
       status = "placeholder", notes = "RSEM count matrix not found at counts/H_T_pipe/.")
}

# --- 1C: composition of the merged HT PacBio annotation (ALL isoforms, not just switching) ---
f1c <- bind_rows(lapply(names(scored), function(lab) {
  scored[[lab]] |>
    mutate(class_code = ifelse(is.na(.data$class_code) | !nzchar(as.character(.data$class_code)),
                               "unassigned", as.character(.data$class_code))) |>
    count(.data$class_code, .data$detect_tier, name = "n") |>
    mutate(dataset_label = lab)
})) |>
  select("dataset_label", "class_code", "detect_tier", "n") |>
  mutate(comparable_across_datasets = FALSE)
emit(f1c, "f1c", "fig1", "Merged HT annotation composition by class code",
     "data/processed/isoformContext_{T_HT_anchor,H_HT_supplementary}.rds",
     notes = "All quantified isoforms. Annotation is shared; per-dataset counts differ by what was quantified.")

# --- 1D: the pre-reduction caveat ---
f1d <- readRDS(file.path(results_tables, "input_object_reduction_check.rds")) |>
  as_tibble() |>
  filter(.data$dataset_label %in% HT_LABELS)
emit(f1d, "f1d", "fig1", "Primary objects are pre-reduced to switching genes",
     "results/tables/input_object_reduction_check.rds",
     notes = "Why significance is taken from the unfiltered context layer, not the primary object.")

# =========================================================================================
# FIGURE 2 -- the isoform-usage response to LPS
# =========================================================================================
message("\n==== Figure 2 ====")

# --- 2A: switching genes and isoforms. NOT comparable across datasets. ---
f2a <- readRDS(file.path(results_tables, "qc_snapshot_gene_level.rds")) |>
  as_tibble() |>
  filter(.data$dataset %in% HT_LABELS) |>
  left_join(
    bind_rows(lapply(names(scored), function(lab) tibble(
      dataset = lab,
      n_isoforms_quantified = nrow(scored[[lab]]),
      n_switching_isoforms = sum(scored[[lab]]$is_switching)
    ))),
    by = "dataset"
  ) |>
  mutate(comparable_across_datasets = FALSE)
emit(f2a, "f2a", "fig2", "Switching genes and isoforms per dataset",
     "results/tables/qc_snapshot_gene_level.rds + context layers",
     notes = "H_HT is batch-corrected primary cells; rates are NOT comparable to T_HT. Facet, never overlay.")

# --- 2B: dIF distribution of switching isoforms ---
f2b <- bind_rows(lapply(names(scored), function(lab) {
  scored[[lab]] |>
    filter(.data$is_switching) |>
    transmute(dataset = lab, gene_name = .data$gene_name, isoform_id = .data$isoform_id,
              dIF = .data$dIF_n, abs_dIF = .data$abs_dIF, detect_tier = .data$detect_tier)
})) |>
  mutate(comparable_across_datasets = FALSE)
emit(f2b, "f2b", "fig2", "dIF distribution of switching isoforms",
     "context layers, score_isoforms()",
     notes = "Signed dIF per switching isoform. Facet by dataset.")

# --- 2C: the abundance floor, derived from replicate IF noise ---
f2c <- readRDS(file.path(results_tables, "expression_if_stability.rds")) |>
  as_tibble() |>
  filter(.data$dataset %in% c("T_HT", "H_HT", HT_LABELS))
floor_rec <- readRDS(file.path(results_tables, "expression_floor_recommendation.rds")) |>
  as_tibble() |>
  filter(.data$dataset %in% c("T_HT", "H_HT", HT_LABELS))
emit(f2c, "f2c", "fig2", "Replicate IF noise against gene abundance",
     "results/tables/expression_if_stability.rds (script 02d)",
     notes = paste0("Consensus floor applied = ",
                    paste(unique(floor_rec$consensus_floor_applied), collapse = "/"),
                    ". snr_at_threshold = min_abs_dif / if_sd_median."))
emit(floor_rec, "f2c_floor", "fig2", "Abundance floor recommendation",
     "results/tables/expression_floor_recommendation.rds (script 02d)")

# --- 2D: where the 0.15 dIF threshold sits ---
f2d <- readRDS(file.path(results_tables, "min_abs_dif_threshold_scan.rds")) |>
  as_tibble() |>
  filter(.data$dataset %in% HT_LABELS) |>
  mutate(is_chosen = abs(.data$min_abs_dif - as.numeric((cfg$significance %||% list())$min_abs_dif %||% 0.15)) < 1e-9)
emit(f2d, "f2d", "fig2", "Sensitivity to the |dIF| threshold",
     "results/tables/min_abs_dif_threshold_scan.rds (script 00)")

# =========================================================================================
# FIGURE 3 -- what PacBio actually added
# =========================================================================================
message("\n==== Figure 3 ====")

# --- 3A: switching isoforms by detectability tier ---
f3a <- bind_rows(lapply(names(scored), function(lab) {
  scored[[lab]] |> filter(.data$is_switching) |> count(.data$detect_tier, name = "n") |>
    mutate(dataset_label = lab)
})) |>
  select("dataset_label", "detect_tier", "n") |>
  mutate(comparable_across_datasets = FALSE)
emit(f3a, "f3a", "fig3", "Switching isoforms by detectability tier",
     "scripts/19_novelty_partition.R tier logic, recomputed HT-only")

# --- 3B: the three novelty definitions side by side ---
f3b <- readRDS(file.path(results_tables, "novelty_detectability_summary.rds")) |>
  as_tibble() |>
  filter(.data$dataset_label %in% HT_LABELS) |>
  select("dataset_label", "n_switching", "pct_novel_id_prefix", "pct_novel_structural",
         "pct_novel_strict", "pct_illumina_findable", "overcount_pp") |>
  tidyr::pivot_longer(cols = dplyr::starts_with("pct_"), names_to = "definition",
                      values_to = "pct") |>
  mutate(definition = factor(
    .data$definition,
    levels = c("pct_novel_id_prefix", "pct_novel_structural", "pct_novel_strict",
               "pct_illumina_findable"),
    labels = c("published (ID prefix)", "structural", "strict (c excluded)",
               "Illumina-findable")
  ))
emit(f3b, "f3b", "fig3", "Novelty under three definitions",
     "results/tables/novelty_detectability_summary.rds (script 19)",
     notes = "published = provenance only; structural = class_code != '='; strict also excludes class c.")

# --- 3C: class code composition, all isoforms vs switching ---
f3c <- bind_rows(lapply(names(scored), function(lab) {
  z <- scored[[lab]] |>
    mutate(class_code = ifelse(is.na(.data$class_code) | !nzchar(as.character(.data$class_code)),
                               "unassigned", as.character(.data$class_code)))
  bind_rows(
    z |> count(.data$class_code, name = "n") |> mutate(set = "all quantified"),
    z |> filter(.data$is_switching) |> count(.data$class_code, name = "n") |> mutate(set = "switching")
  ) |> mutate(dataset_label = lab)
})) |>
  group_by(.data$dataset_label, .data$set) |>
  mutate(pct = 100 * .data$n / sum(.data$n)) |>
  ungroup() |>
  select("dataset_label", "set", "class_code", "n", "pct")
emit(f3c, "f3c", "fig3", "Class code composition, all vs switching",
     "context layers")

# --- 3D: effect size novel vs known, WITHIN dataset only ---
f3d <- bind_rows(lapply(names(scored), function(lab) {
  scored[[lab]] |>
    filter(.data$is_switching, !is.na(.data$detect_tier)) |>
    mutate(structurally_novel = .data$detect_tier %in% c("pb_contained", "pb_novel")) |>
    group_by(.data$structurally_novel) |>
    summarize(n = dplyr::n(),
              median_abs_dIF = stats::median(.data$abs_dIF, na.rm = TRUE),
              q25 = stats::quantile(.data$abs_dIF, 0.25, na.rm = TRUE),
              q75 = stats::quantile(.data$abs_dIF, 0.75, na.rm = TRUE),
              .groups = "drop") |>
    mutate(dataset_label = lab)
})) |>
  select("dataset_label", "structurally_novel", "n", "median_abs_dIF", "q25", "q75") |>
  mutate(comparable_across_datasets = FALSE)
emit(f3d, "f3d", "fig3", "Effect size, structurally novel vs reference-matched",
     "context layers",
     notes = "Compare WITHIN a dataset only. Lead with the median difference, not a p-value.")

# =========================================================================================
# FIGURE 4 -- candidate selection
# =========================================================================================
message("\n==== Figure 4 ====")

sig <- cfg$significance %||% list()
g_q <- as.numeric(sig$gene_q %||% 0.05)
i_q <- as.numeric(sig$isoform_q %||% 0.05)
m_dif <- as.numeric(sig$min_abs_dif %||% 0.15)
min_ge <- as.numeric(sig$min_gene_expression %||% 12)

# --- 4A: the selection funnel, stage by stage ---
# The funnel must use the SAME predicates score_isoforms() uses, or the stages will not
# nest. In particular the abundance floor is applied to `gene_expression` -- which
# score_isoforms derives as the per-gene sum of iso_value_1 + iso_value_2 -- and NOT to
# `gene_overall_mean`, which the context layer does not carry at all. Using the wrong column
# silently yields an all-NA stage and a funnel that widens instead of narrowing.
f4a <- bind_rows(lapply(names(scored), function(lab) {
  z <- scored[[lab]]
  stopifnot("gene_expression" %in% names(z), "low_expression" %in% names(z))
  s1 <- rep(TRUE, nrow(z))
  s2 <- s1 & is.finite(z$q_g) & z$q_g < g_q
  s3 <- s2 & !(z$low_expression %in% TRUE)
  s4 <- s3 & is.finite(z$q_i) & z$q_i < i_q & is.finite(z$dIF_n) & abs(z$dIF_n) >= m_dif
  # s4 must reproduce is_switching exactly; if it does not, the predicates have drifted.
  stopifnot(identical(
    sort(unique(z$gene_id[s4])), sort(unique(z$gene_id[z$is_switching]))
  ))
  tibble(
    dataset_label = lab,
    stage = c("genes quantified", "gene q < cutoff", "above abundance floor",
              "|dIF| >= cutoff (switching)", "assayable at the bench"),
    stage_order = 1:5,
    n_genes = c(
      dplyr::n_distinct(z$gene_id[s1]),
      dplyr::n_distinct(z$gene_id[s2]),
      dplyr::n_distinct(z$gene_id[s3]),
      dplyr::n_distinct(z$gene_id[s4]),
      NA_integer_
    ),
    status = c(rep("real", 4L), "placeholder")
  )
}))
emit(f4a, "f4a", "fig4", "Candidate selection funnel",
     "context layers + config/config.yml thresholds", status = "placeholder",
     notes = paste0("Stages 1-4 real (gene_q<", g_q, ", floor=", min_ge, ", |dIF|>=", m_dif,
                    "). Stage 5 needs an HT-space assayability screen: H_T.gff3 + ",
                    "transcripts_H-T.fa, both on /Volumes/Expansion."))

# --- 4B: anchor (T_HT) against cross-system corroboration (H_HT) ---
# `q1_candidates` carries three arms: anchor = T_HT, reannot = T_UT, cross = H_HT. Only the
# anchor and cross arms are in scope; every *_reannot column is UT space and is dropped here.
q1 <- readRDS(file.path(results_tables, "q1_candidates.rds")) |> as_tibble()
dropped <- grep("reannot", names(q1), value = TRUE)
message("  dropping UT-space columns from q1_candidates: ", paste(dropped, collapse = ", "))
f4b <- q1 |>
  select(-dplyr::all_of(dropped), -dplyr::any_of(c("annotation_robust", "annotation_status"))) |>
  mutate(
    corroboration = dplyr::case_when(
      .data$switches_anchor & .data$switches_cross ~ "both T_HT and H_HT",
      .data$switches_anchor & !.data$testable_cross ~ "T_HT only (not testable in H_HT)",
      .data$switches_anchor ~ "T_HT only (testable in H_HT, not called)",
      .data$switches_cross ~ "H_HT only",
      TRUE ~ "neither"
    )
  ) |>
  # `called_in` is dropped with the *_reannot columns: its values ("T_UT only", "both
  # annotations") describe agreement between the T_HT and T_UT ANNOTATIONS, so the column is
  # UT-dependent regardless of its name. `visibility` is anchor-derived and stays.
  select("rank", "gene_name", "dif_anchor", "lfc_anchor", "dif_cross", "expr_anchor",
         "n_iso_anchor", "class_codes", "n_switching_novel_pb", "n_switching_contained",
         "switches_anchor", "testable_anchor", "switches_cross", "testable_cross",
         "cross_genotype_support", "corroboration", "visibility")
emit(f4b, "f4b", "fig4", "T_HT anchor against H_HT cross-system corroboration",
     "results/tables/q1_candidates.rds (script 11), *_reannot columns dropped",
     notes = "Presence/absence corroboration only. Do NOT read dif_anchor vs dif_cross as comparable effect sizes.")

# --- 4C: switch anatomy / AS event composition. PLACEHOLDER. ---
f4c <- tibble(
  as_event = c("alternative first exon", "exon skipping", "alternative donor",
               "alternative acceptor", "intron retention", "alternative last exon"),
  n_switches = NA_integer_,
  required_input = "H_T.gff3 (exon coordinates) on /Volumes/Expansion"
)
emit(f4c, "f4c", "fig4", "Switch anatomy by AS event class", "not yet derivable in HT space",
     status = "placeholder",
     notes = paste0("Script 15 computed this from U_T.gff3 (UT space): 54 of 63 pairs were ",
                    "alternative-first-exon. The HT equivalent needs H_T.gff3, which is on ",
                    "/Volumes/Expansion and outside the current host grant."))

# --- 4D: the shortlist with platform class. PLACEHOLDER. ---
f4d <- tibble(
  gene_name = NA_character_, priority_score = NA_integer_, detect_class = NA_character_,
  assayable = NA, longest_unique_bp = NA_integer_, coding_change = NA_character_,
  required_input = "HT-space re-run of scripts 13, 15, 16 (needs H_T.gff3 + transcripts_H-T.fa)"
)[0, ]
emit(f4d, "f4d", "fig4", "Candidate shortlist with platform class",
     "not yet derivable in HT space", status = "placeholder",
     notes = paste0("The existing shortlist (candidate_detectability_by_switch) is UT space: ",
                    "all 84 candidate isoform ids resolve in T_UT/U_UT, only 42 of 80 in T_HT. ",
                    "It cannot be carried into an HT-only manuscript."))

# ---- index --------------------------------------------------------------------------

idx <- bind_rows(index_rows)
utils::write.csv(idx, file.path(man_dir, "manuscript_panel_index.csv"),
                 row.names = FALSE, fileEncoding = "UTF-8", na = "")

message("\n---- panel index ----")
print(as.data.frame(idx |> select("panel", "figure", "status", "n_rows", "title")),
      row.names = FALSE)
message("\n  real panels:        ", sum(idx$status == "real"))
message("  placeholder panels: ", sum(idx$status == "placeholder"))
message("  written to:         ", man_dir)

write_run_manifest(script = "21_manuscript_data.R", cfg = cfg, root = root)
message("\n21_manuscript_data.R: done")
