#!/usr/bin/env Rscript
# 04: UT reference comparison — T_UT (WT) vs U_UT (UBL5 knockout)
#     Question: which isoform-switch events differ between knockout (U) and wildtype (T)?
#
# Outputs under results/tables/ and results/figures/ut_t_vs_u/:
#   - per-dataset ranked switching genes
#   - isoform- and gene-level overlap (shared / T-only / U-only)
#   - concordance and delta-dIF statistics
#   - figures highlighting overlap + effect-size contrasts
#   - per-gene explorer panels for top shared / T-only / U-only genes
#   - IsoformSwitchAnalyzeR switchPlots for showcase genes
#
# Run (after 01): Rscript scripts/04_comparison_T_vs_U.R

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
fig_dir <- ensure_dir(file.path(results_figures, "ut_t_vs_u"))
gene_fig_dir <- ensure_dir(file.path(fig_dir, "gene_explorer"))
switch_fig_dir <- ensure_dir(file.path(fig_dir, "switch_plots"))

sig <- cfg$significance %||% list()
iso_q_cutoff <- as.numeric(sig$isoform_q %||% 0.05)
gene_q_cutoff <- as.numeric(sig$gene_q %||% 0.05)
min_abs_dif <- as.numeric(sig$min_abs_dif %||% 0.0)

top_n_genes <- as.integer(analysis_param(cfg, "top_n_genes", 30L))
n_gene_panels_per_class <- as.integer(analysis_param(cfg, "gene_panels_per_class", 12L))
n_switch_plots_per_class <- as.integer(analysis_param(cfg, "switch_plots_per_class", 6L))
explorer_max_isoforms <- as.integer(analysis_param(cfg, "explorer_max_isoforms", 12L))

# sanitize(), min_finite(), max_abs_finite(), is_real_gene_symbol(), pick_gene_symbol()
# and score_isoforms() come from utils/helper_functions.R -- this script used to keep
# its own copies alongside near-identical ones in 03 and 06.

theme_set(
  theme_bw(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold", size = 12),
      strip.text = element_text(face = "bold")
    )
)

gene_rank_from_scored <- function(z, dataset_key, dataset_label) {
  z |>
    group_by(.data$gene_id, .data$gene_name) |>
    summarize(
      n_isoforms = n(),
      n_switching_isoforms = sum(.data$is_switching, na.rm = TRUE),
      n_novel_switching_isoforms = sum(.data$is_switching & .data$is_novel, na.rm = TRUE),
      n_switching_up = sum(.data$is_switching & is.finite(.data$dIF_n) & .data$dIF_n > 0, na.rm = TRUE),
      n_switching_down = sum(.data$is_switching & is.finite(.data$dIF_n) & .data$dIF_n < 0, na.rm = TRUE),
      max_abs_dif_switching = max_abs_finite(.data$dIF_n[.data$is_switching]),
      min_isoform_switch_q = min_finite(.data$q_i[.data$is_switching]),
      min_gene_switch_q = min_finite(.data$q_g[.data$is_switching]),
      novel_involved = any(.data$is_switching & .data$is_novel, na.rm = TRUE),
      .groups = "drop"
    ) |>
    mutate(
      has_multi_isoform_opposing_switch = (.data$n_switching_isoforms >= 2L) &
        (.data$n_switching_up > 0L) & (.data$n_switching_down > 0L),
      dataset_key = dataset_key,
      dataset_label = dataset_label
    ) |>
    filter(.data$n_switching_isoforms > 0L) |>
    arrange(
      desc(.data$has_multi_isoform_opposing_switch),
      desc(.data$max_abs_dif_switching),
      desc(.data$n_switching_isoforms),
      .data$min_isoform_switch_q
    ) |>
    mutate(rank = row_number())
}

write_out <- function(x, stem, csv = NULL) {
  write_table_pair(x, results_tables, stem, cfg = cfg, csv = csv)
}
isoform_level_csv <- isTRUE((cfg$output %||% list())$csv_twin_isoform_level %||% FALSE)

load_processed <- function(dataset_key) {
  ds <- cfg$datasets[[dataset_key]] %||% NULL
  if (is.null(ds)) stop("Dataset missing from config: ", dataset_key)
  label <- sanitize(as.character(ds$label %||% dataset_key)[1L])
  rds <- file.path(processed_dir, paste0("isoformFeatures_", label, ".rds"))
  if (!file.exists(rds)) {
    stop("Missing processed isoform table: ", rds, " (run scripts/01_load_data.R)")
  }
  list(ds = ds, label = label, iso = readRDS(rds))
}

annotate_isa_for_plotting <- function(isa_obj, processed_iso) {
  if (is.null(isa_obj) || is.null(processed_iso)) return(isa_obj)
  proc <- as_tibble(processed_iso)
  if (!all(c("gene_id", "isoform_id") %in% names(proc))) return(isa_obj)

  if ("gene_name" %in% names(proc) && "isoformFeatures" %in% names(isa_obj)) {
    gmap <- proc |>
      group_by(.data$gene_id) |>
      summarize(gene_name = pick_gene_symbol(.data$gene_name), .groups = "drop") |>
      filter(is_real_gene_symbol(.data$gene_name)) |>
      transmute(gene_id = as.character(.data$gene_id), gene_name = as.character(.data$gene_name))
    feat <- isa_obj$isoformFeatures
    feat$gene_id <- as.character(feat$gene_id)
    m <- match(feat$gene_id, gmap$gene_id)
    hit <- !is.na(m)
    if (any(hit)) {
      feat$gene_name[hit] <- gmap$gene_name[m[hit]]
      isa_obj$isoformFeatures <- feat
    }
  }

  if (!"is_novel_pacbio" %in% names(proc)) return(isa_obj)
  novel_ids <- unique(as.character(proc$isoform_id[proc$is_novel_pacbio %in% TRUE]))
  novel_ids <- novel_ids[!is.na(novel_ids) & nzchar(novel_ids) & !grepl(" \\(PB\\)$", novel_ids)]
  if (!length(novel_ids)) return(isa_obj)
  new_ids <- paste0(novel_ids, " (PB)")
  names(new_ids) <- novel_ids
  remap_ids <- function(ids) {
    ids <- as.character(ids)
    hit <- ids %in% names(new_ids)
    ids[hit] <- unname(new_ids[ids[hit]])
    ids
  }
  for (nm in names(isa_obj)) {
    x <- isa_obj[[nm]]
    if (is.data.frame(x) && "isoform_id" %in% names(x)) {
      x$isoform_id <- remap_ids(x$isoform_id)
      isa_obj[[nm]] <- x
    } else if (inherits(x, "GRanges")) {
      mc <- S4Vectors::mcols(x)
      if ("isoform_id" %in% names(mc)) {
        mc$isoform_id <- remap_ids(mc$isoform_id)
        S4Vectors::mcols(x) <- mc
        isa_obj[[nm]] <- x
      }
    } else if (inherits(x, "DNAStringSet") || inherits(x, "XStringSet")) {
      nms <- names(x)
      if (!is.null(nms)) {
        names(x) <- remap_ids(nms)
        isa_obj[[nm]] <- x
      }
    }
  }
  isa_obj
}

#' Draw an ISA switchPlot, trying the gene symbol first and the gene_id as a fallback.
#'
#' Device handling detail that used to be wrong here: `on.exit()` called inside a
#' `tryCatch({...})` registers on the *enclosing function* frame, not per iteration.
#' Devices therefore stayed open until the function returned and then closed LIFO, so
#' (a) the `file.info()$size` check ran before the PNG had been flushed, and (b) when
#' the symbol attempt failed and the gene_id attempt succeeded, the failed device --
#' opened first, closed last, on the same path -- overwrote the good plot. Each attempt
#' now renders to its own temp file, closes its device immediately, and is promoted to
#' the destination only after it is verified non-empty.
make_switch_plot <- function(isa_obj, gene_id, gene_name, out_png) {
  if (is.null(isa_obj) || !requireNamespace("IsoformSwitchAnalyzeR", quietly = TRUE)) {
    return(FALSE)
  }
  switch_fn <- get("switchPlot", envir = asNamespace("IsoformSwitchAnalyzeR"))
  candidates <- character()
  if (is_real_gene_symbol(gene_name)) candidates <- c(candidates, as.character(gene_name)[1L])
  if (!is.na(gene_id) && nzchar(as.character(gene_id)[1L])) {
    candidates <- c(candidates, as.character(gene_id)[1L])
  }
  candidates <- unique(candidates)
  for (g in candidates) {
    tmp <- tempfile(fileext = ".png")
    ok <- tryCatch({
      grDevices::png(filename = tmp, width = 2400, height = 1500, res = 180)
      on.exit(if (grDevices::dev.cur() > 1L) grDevices::dev.off(), add = TRUE)
      do.call(switch_fn, list(switchAnalyzeRlist = isa_obj, gene = g, plotTopology = FALSE))
      grDevices::dev.off()
      TRUE
    }, error = function(e) FALSE)
    if (grDevices::dev.cur() > 1L) grDevices::dev.off()
    if (ok && file.exists(tmp) && isTRUE(file.info(tmp)$size > 0)) {
      file.copy(tmp, out_png, overwrite = TRUE)
      unlink(tmp)
      return(TRUE)
    }
    unlink(tmp)
  }
  FALSE
}

# ---- Load and score UT datasets ----
message("Loading T_UT and U_UT ...")
t_pack <- load_processed("T_UT")
u_pack <- load_processed("U_UT")
t_iso <- score_isoforms(t_pack$iso, cfg, "T_UT", t_pack$label)
u_iso <- score_isoforms(u_pack$iso, cfg, "U_UT", u_pack$label)

# Unfiltered context (written by 01 when isa_unfiltered_path is configured). Used for
# plotting and for telling "tested and not switching" apart from "not detected"; the
# switching calls above are unaffected.
ctx_t <- load_context_table(processed_dir, t_pack$label, "features")
ctx_u <- load_context_table(processed_dir, u_pack$label, "features")
rep_t <- load_context_table(processed_dir, t_pack$label, "rep_if")
rep_u <- load_context_table(processed_dir, u_pack$label, "rep_if")
has_context <- !is.null(ctx_t) && !is.null(ctx_u)
if (has_context) {
  message(
    "Context layer: T ", nrow(ctx_t), " isoforms / ", dplyr::n_distinct(ctx_t$gene_id),
    " genes; U ", nrow(ctx_u), " / ", dplyr::n_distinct(ctx_u$gene_id), " genes."
  )
} else {
  message(
    "No unfiltered context tables found. Overlap classes fall back to object retention, ",
    "and gene explorer panels can only show the dataset a gene was significant in. ",
    "Set isa_unfiltered_path in config and rerun 01 to enable it."
  )
}

t_genes <- gene_rank_from_scored(t_iso, "T_UT", t_pack$label)
u_genes <- gene_rank_from_scored(u_iso, "U_UT", u_pack$label)
write_out(t_genes |> slice_head(n = min(top_n_genes, nrow(t_genes))), "top_switching_genes_T_UT")
write_out(u_genes |> slice_head(n = min(top_n_genes, nrow(u_genes))), "top_switching_genes_U_UT")
write_out(t_genes, "switching_genes_T_UT_all")
write_out(u_genes, "switching_genes_U_UT_all")

# ---- Gene-level overlap by gene symbol (transcriptomes share symbols more than XLOC IDs) ----
t_gene_key <- t_genes |>
  filter(is_real_gene_symbol(.data$gene_name)) |>
  transmute(
    gene_name = as.character(.data$gene_name),
    T_gene_id = as.character(.data$gene_id),
    T_n_switching = .data$n_switching_isoforms,
    T_max_abs_dif = .data$max_abs_dif_switching,
    T_min_q = .data$min_isoform_switch_q,
    T_novel_involved = .data$novel_involved,
    T_opposing = .data$has_multi_isoform_opposing_switch,
    T_rank = .data$rank
  )
u_gene_key <- u_genes |>
  filter(is_real_gene_symbol(.data$gene_name)) |>
  transmute(
    gene_name = as.character(.data$gene_name),
    U_gene_id = as.character(.data$gene_id),
    U_n_switching = .data$n_switching_isoforms,
    U_max_abs_dif = .data$max_abs_dif_switching,
    U_min_q = .data$min_isoform_switch_q,
    U_novel_involved = .data$novel_involved,
    U_opposing = .data$has_multi_isoform_opposing_switch,
    U_rank = .data$rank
  )

gene_overlap <- full_join(t_gene_key, u_gene_key, by = "gene_name") |>
  mutate(
    significant_in_T = !is.na(.data$T_n_switching) & .data$T_n_switching > 0,
    significant_in_U = !is.na(.data$U_n_switching) & .data$U_n_switching > 0,
    overlap_class = case_when(
      .data$significant_in_T & .data$significant_in_U ~ "Shared (T and U)",
      .data$significant_in_T & !.data$significant_in_U ~ "T-only",
      !.data$significant_in_T & .data$significant_in_U ~ "U-only",
      TRUE ~ "Neither"
    ),
    max_abs_dif_either = pmax(.data$T_max_abs_dif, .data$U_max_abs_dif, na.rm = TRUE),
    delta_max_abs_dif = dplyr::case_when(
      is.finite(.data$T_max_abs_dif) & is.finite(.data$U_max_abs_dif) ~
        .data$U_max_abs_dif - .data$T_max_abs_dif,
      TRUE ~ NA_real_
    )
  ) |>
  filter(.data$overlap_class != "Neither") |>
  arrange(.data$overlap_class, desc(.data$max_abs_dif_either))

write_out(gene_overlap, "ut_gene_overlap_T_vs_U")

gene_overlap_summary <- gene_overlap |>
  count(.data$overlap_class, name = "n_genes") |>
  mutate(pct = 100 * .data$n_genes / sum(.data$n_genes))
write_out(gene_overlap_summary, "ut_gene_overlap_T_vs_U_summary")

# ---- Isoform-level join on isoform_id (shared TCONS within UT reference) ----
t_iso_key <- t_iso |>
  transmute(
    isoform_id = as.character(.data$isoform_id),
    gene_id = as.character(.data$gene_id),
    gene_name = as.character(.data$gene_name),
    oId = as.character(.data$oId),
    class_code = as.character(.data$class_code),
    T_dIF = .data$dIF_n,
    T_abs_dIF = .data$abs_dIF,
    T_q = .data$q_i,
    T_gene_q = .data$q_g,
    T_switching = .data$is_switching,
    T_novel = .data$is_novel
  )
u_iso_key <- u_iso |>
  transmute(
    isoform_id = as.character(.data$isoform_id),
    U_gene_id = as.character(.data$gene_id),
    U_gene_name = as.character(.data$gene_name),
    # class_code / oId are properties of the transcript in the shared UT annotation, not
    # of the experiment. They were previously carried only from the T side, so isoforms
    # significant in U but absent from T's reduced object lost them (78 of 85). Verified
    # identical for all 172 isoforms present in both tables.
    U_class_code = as.character(.data$class_code),
    U_oId = as.character(.data$oId),
    U_dIF = .data$dIF_n,
    U_abs_dIF = .data$abs_dIF,
    U_q = .data$q_i,
    U_gene_q = .data$q_g,
    U_switching = .data$is_switching,
    U_novel = .data$is_novel
  )

iso_overlap <- full_join(t_iso_key, u_iso_key, by = "isoform_id") |>
  mutate(
    present_in_T = !is.na(.data$T_dIF) | !is.na(.data$T_q),
    present_in_U = !is.na(.data$U_dIF) | !is.na(.data$U_q),
    present_in_both = .data$present_in_T & .data$present_in_U,
    gene_name = dplyr::coalesce(
      ifelse(is_real_gene_symbol(.data$gene_name), .data$gene_name, NA_character_),
      ifelse(is_real_gene_symbol(.data$U_gene_name), .data$U_gene_name, NA_character_),
      .data$gene_name,
      .data$U_gene_name
    ),
    gene_id = dplyr::coalesce(.data$gene_id, .data$U_gene_id),
    class_code = dplyr::coalesce(.data$class_code, .data$U_class_code),
    oId = dplyr::coalesce(.data$oId, .data$U_oId),
    is_novel_pb = (.data$T_novel %in% TRUE) | (.data$U_novel %in% TRUE),
    T_switching = .data$T_switching %in% TRUE,
    U_switching = .data$U_switching %in% TRUE,
    # With the unfiltered context available, "tested in the other dataset" is a real
    # statement: the isoform was quantified there and simply did not switch. Without it,
    # membership in the saved (reduced) object is all we know, so the labels say
    # "retained" instead and must not be read as tested-and-non-significant.
    tested_in_T = if (has_context) {
      .data$isoform_id %in% ctx_t$isoform_id
    } else {
      .data$present_in_T
    },
    tested_in_U = if (has_context) {
      .data$isoform_id %in% ctx_u$isoform_id
    } else {
      .data$present_in_U
    },
    overlap_class = if (has_context) {
      case_when(
        .data$T_switching & .data$U_switching ~ "Shared significant",
        .data$T_switching & !.data$U_switching & .data$tested_in_U ~ "T-significant / tested in U, not switching",
        .data$T_switching & !.data$tested_in_U ~ "T-significant / not detected in U",
        .data$U_switching & !.data$T_switching & .data$tested_in_T ~ "U-significant / tested in T, not switching",
        .data$U_switching & !.data$tested_in_T ~ "U-significant / not detected in T",
        TRUE ~ "Non-significant / other"
      )
    } else {
      case_when(
        .data$T_switching & .data$U_switching ~ "Shared significant",
        .data$T_switching & !.data$U_switching & .data$present_in_U ~ "T-significant / retained in U, not significant",
        .data$T_switching & !.data$present_in_U ~ "T-significant / not retained in U",
        .data$U_switching & !.data$T_switching & .data$present_in_T ~ "U-significant / retained in T, not significant",
        .data$U_switching & !.data$present_in_T ~ "U-significant / not retained in T",
        TRUE ~ "Non-significant / other"
      )
    },
    same_direction = dplyr::case_when(
      is.finite(.data$T_dIF) & is.finite(.data$U_dIF) ~ sign(.data$T_dIF) == sign(.data$U_dIF),
      TRUE ~ NA
    ),
    delta_dIF = dplyr::case_when(
      is.finite(.data$T_dIF) & is.finite(.data$U_dIF) ~ .data$U_dIF - .data$T_dIF,
      TRUE ~ NA_real_
    ),
    abs_delta_dIF = abs(.data$delta_dIF)
  )

write_out(iso_overlap, "ut_isoform_overlap_T_vs_U_all", csv = isoform_level_csv)

iso_focus <- iso_overlap |>
  filter(.data$T_switching | .data$U_switching) |>
  arrange(desc(pmax(.data$T_abs_dIF, .data$U_abs_dIF, na.rm = TRUE)))
write_out(iso_focus, "ut_isoform_overlap_T_vs_U_significant")

# ---- Annotation class-code composition by overlap class ----
# Does a switch that happens in only one genotype involve a structurally different kind of
# transcript than one that happens in both? gffcompare class codes describe how each
# transcript model relates to the reference annotation:
#   "="  intron chain identical to a reference transcript
#   "j"  multi-exon with at least one novel splice-junction combination
#   "c"  contained within a reference transcript (shorter/fragmentary model)
# The rest are individually rare here and are pooled as "other".
cc_label <- c(
  "=" = "= (matches reference)",
  "j" = "j (novel junction combination)",
  "c" = "c (contained in reference)"
)
group_class_code <- function(x) {
  x <- as.character(x)
  x[is.na(x) | !nzchar(x)] <- "other"
  out <- unname(cc_label[x])
  out[is.na(out)] <- "other"
  factor(out, levels = c(cc_label, "other"))
}
# Collapse the five outcome classes into the three the question is about.
overlap_group <- function(x) {
  dplyr::case_when(
    grepl("^Shared", x) ~ "Shared (both genotypes)",
    grepl("^T-significant", x) ~ "WT-only",
    grepl("^U-significant", x) ~ "KO-only",
    TRUE ~ NA_character_
  )
}

cc_dat <- iso_focus |>
  mutate(
    cc = group_class_code(.data$class_code),
    grp = overlap_group(.data$overlap_class)
  ) |>
  filter(!is.na(.data$grp))

# Fine-grained: percentage within each of the five outcome classes.
class_code_by_overlap <- cc_dat |>
  count(.data$overlap_class, .data$cc, name = "n") |>
  group_by(.data$overlap_class) |>
  mutate(pct_within_class = 100 * .data$n / sum(.data$n), n_in_class = sum(.data$n)) |>
  ungroup() |>
  arrange(.data$overlap_class, desc(.data$n))
write_out(class_code_by_overlap, "ut_class_code_by_overlap_class")

# Coarse: the three-group comparison.
class_code_by_group <- cc_dat |>
  count(.data$grp, .data$cc, name = "n") |>
  group_by(.data$grp) |>
  mutate(pct_within_group = 100 * .data$n / sum(.data$n), n_in_group = sum(.data$n)) |>
  ungroup() |>
  arrange(.data$grp, desc(.data$n))
write_out(class_code_by_group, "ut_class_code_by_overlap_group")

# Tests. Counts are small (the shared set especially), so Fisher throughout; the omnibus
# uses simulation because the table is larger than 2x2. These are exploratory and are NOT
# corrected for multiple testing -- four related comparisons on the same data.
fisher_safe <- function(tab, simulate = FALSE) {
  tryCatch(
    if (simulate) {
      stats::fisher.test(tab, simulate.p.value = TRUE, B = 20000)
    } else {
      stats::fisher.test(tab)
    },
    error = function(e) NULL
  )
}
tab_of <- function(groups) {
  d <- cc_dat[cc_dat$grp %in% groups, , drop = FALSE]
  t <- table(droplevels(d$cc), d$grp)
  t[, colSums(t) > 0, drop = FALSE]
}

tests <- list()
omni <- tab_of(c("Shared (both genotypes)", "WT-only", "KO-only"))
ft <- fisher_safe(omni, simulate = TRUE)
tests[[length(tests) + 1L]] <- tibble(
  comparison = "All three groups x class code",
  test = "Fisher (simulated p, B=20000)",
  n = sum(omni),
  odds_ratio = NA_real_, conf_low = NA_real_, conf_high = NA_real_,
  p_value = if (!is.null(ft)) ft$p.value else NA_real_
)
for (other in c("WT-only", "KO-only")) {
  tt <- tab_of(c("Shared (both genotypes)", other))
  ft <- fisher_safe(tt, simulate = ncol(tt) > 2 || nrow(tt) > 2)
  tests[[length(tests) + 1L]] <- tibble(
    comparison = paste0("Shared vs ", other, " x class code"),
    test = if (nrow(tt) > 2) "Fisher (simulated p)" else "Fisher exact",
    n = sum(tt),
    odds_ratio = NA_real_, conf_low = NA_real_, conf_high = NA_real_,
    p_value = if (!is.null(ft)) ft$p.value else NA_real_
  )
}
# Focused 2x2: novel-junction transcripts ("j") versus everything else.
for (other in c("WT-only", "KO-only")) {
  d <- cc_dat[cc_dat$grp %in% c("Shared (both genotypes)", other), , drop = FALSE]
  j2 <- table(
    factor(ifelse(grepl("^j ", as.character(d$cc)), "j", "not j"), levels = c("not j", "j")),
    factor(d$grp, levels = c("Shared (both genotypes)", other))
  )
  ft <- fisher_safe(j2)
  tests[[length(tests) + 1L]] <- tibble(
    comparison = paste0("Shared vs ", other, ": novel junction (j) vs rest"),
    test = "Fisher exact (2x2)",
    n = sum(j2),
    odds_ratio = if (!is.null(ft)) unname(ft$estimate) else NA_real_,
    conf_low = if (!is.null(ft)) ft$conf.int[[1L]] else NA_real_,
    conf_high = if (!is.null(ft)) ft$conf.int[[2L]] else NA_real_,
    p_value = if (!is.null(ft)) ft$p.value else NA_real_
  )
}
class_code_tests <- bind_rows(tests) |>
  mutate(note = "Exploratory; not corrected for multiple testing. Shared group is small.")
write_out(class_code_tests, "ut_class_code_composition_tests")
print(class_code_by_group)
print(class_code_tests)

p_cc <- ggplot(
  cc_dat |> count(.data$grp, .data$cc, name = "n") |>
    group_by(.data$grp) |>
    mutate(pct = 100 * .data$n / sum(.data$n), lab = paste0("n=", sum(.data$n))) |>
    ungroup(),
  aes(x = .data$grp, y = .data$pct, fill = .data$cc)
) +
  geom_col(colour = "grey25", linewidth = 0.2, width = 0.7) +
  geom_text(aes(x = .data$grp, y = 103, label = .data$lab), inherit.aes = FALSE,
            data = ~ distinct(.x, .data$grp, .data$lab), size = 3.2) +
  scale_fill_manual(
    values = c(
      "= (matches reference)" = "#b3cde3",
      "j (novel junction combination)" = "#8856a7",
      "c (contained in reference)" = "#9ebcda",
      "other" = "grey75"
    ),
    name = "Annotation class code"
  ) +
  labs(
    title = "Structural class of switching isoforms by outcome group",
    subtitle = "Does a genotype-specific switch involve a different kind of transcript than a shared one?",
    x = NULL, y = "% of switching isoforms in group"
  )
ggsave(file.path(fig_dir, "fig_ut_class_code_by_overlap_group.png"), p_cc, width = 8.5, height = 5, dpi = 200)

# ---- Statistics ----
#
# NO ENRICHMENT TEST IS REPORTED HERE, deliberately.
#
# This section used to run fisher.test() on a 2x2 of "significant in T" x "significant
# in U" over gene symbols present in both datasets, and report the odds ratio and p in
# the headline summary. That test is not interpretable for these inputs: both ISA
# objects were saved after reduceToSwitchingGenes = TRUE, so a gene is only "present"
# if it had already been called significant. The background was 25 genes and the table
# was 3/1/3/18 -- an association measured inside a set already selected on the outcome.
# What remains is the descriptive retention/overlap accounting, which is honest about
# what it counts. See docs/REVIEW_CHANGES.md.
bg_t <- unique(t_iso$gene_name[is_real_gene_symbol(t_iso$gene_name)])
bg_u <- unique(u_iso$gene_name[is_real_gene_symbol(u_iso$gene_name)])
bg_both <- intersect(bg_t, bg_u)
in_t_sig_both <- bg_both %in% t_gene_key$gene_name
in_u_sig_both <- bg_both %in% u_gene_key$gene_name
retention_tbl <- tibble(
  symbols_retained_in_T_object = length(bg_t),
  symbols_retained_in_U_object = length(bg_u),
  symbols_retained_in_both = length(bg_both),
  of_those_switching_in_T_only = sum(in_t_sig_both & !in_u_sig_both),
  of_those_switching_in_U_only = sum(!in_t_sig_both & in_u_sig_both),
  of_those_switching_in_both = sum(in_t_sig_both & in_u_sig_both),
  of_those_switching_in_neither = sum(!in_t_sig_both & !in_u_sig_both),
  note = paste(
    "Both objects were reduced to significant switching genes before saving;",
    "'retained' is not 'tested'. Gene overlap cannot exceed symbols_retained_in_both."
  )
)
write_out(retention_tbl, "ut_T_vs_U_retention_accounting")

# With the unfiltered context there IS a genuine tested background -- genes quantified in
# both datasets, most of which are not switching -- so a co-occurrence test is meaningful
# here in a way it was not against the reduced objects. It is computed only when the
# context is present, and reported separately from the retention accounting above so the
# two can never be confused.
background_enrichment <- NULL
if (has_context) {
  bg_genes <- intersect(unique(as.character(ctx_t$gene_id)), unique(as.character(ctx_u$gene_id)))
  t_sw <- bg_genes %in% as.character(t_genes$gene_id)
  u_sw <- bg_genes %in% as.character(u_genes$gene_id)
  bg_tab <- matrix(
    c(sum(!t_sw & !u_sw), sum(!t_sw & u_sw), sum(t_sw & !u_sw), sum(t_sw & u_sw)),
    nrow = 2, byrow = TRUE,
    dimnames = list(T_switching = c("FALSE", "TRUE"), U_switching = c("FALSE", "TRUE"))
  )
  ft <- tryCatch(stats::fisher.test(bg_tab), error = function(e) NULL)
  background_enrichment <- tibble(
    n_genes_tested_in_both = length(bg_genes),
    n_switching_in_T = sum(t_sw),
    n_switching_in_U = sum(u_sw),
    n_switching_in_both = sum(t_sw & u_sw),
    expected_in_both_if_independent = length(bg_genes) * mean(t_sw) * mean(u_sw),
    odds_ratio = if (!is.null(ft)) unname(ft$estimate) else NA_real_,
    conf_low = if (!is.null(ft)) ft$conf.int[[1L]] else NA_real_,
    conf_high = if (!is.null(ft)) ft$conf.int[[2L]] else NA_real_,
    p_value = if (!is.null(ft)) ft$p.value else NA_real_,
    background = "genes quantified in both unfiltered UT objects"
  )
  write_out(background_enrichment, "ut_T_vs_U_background_enrichment")
  write_out(
    data.frame(
      T_switching = c("FALSE", "TRUE"),
      U_switching_FALSE = bg_tab[, "FALSE"],
      U_switching_TRUE = bg_tab[, "TRUE"],
      stringsAsFactors = FALSE, row.names = NULL
    ),
    "ut_T_vs_U_background_contingency"
  )
  message(
    "Background enrichment over ", length(bg_genes), " genes tested in both: OR = ",
    signif(background_enrichment$odds_ratio, 4), ", p = ",
    signif(background_enrichment$p_value, 4)
  )
}

# Condition on significance in ONE dataset, then describe the effect in the OTHER.
# This is the defensible concordance measure: the shared-significant statistics below
# select on the outcome in both datasets, which manufactures agreement.
one_way_concordance <- function(df, sig_col, this_dIF, other_dIF, label) {
  d <- df[df[[sig_col]] %in% TRUE & is.finite(df[[this_dIF]]) & is.finite(df[[other_dIF]]), ]
  if (!nrow(d)) {
    return(tibble(
      direction = label, n = 0L, median_abs_dIF_conditioned = NA_real_,
      median_abs_dIF_other = NA_real_, pct_same_direction_in_other = NA_real_,
      pct_passing_min_abs_dif_in_other = NA_real_
    ))
  }
  tibble(
    direction = label,
    n = nrow(d),
    median_abs_dIF_conditioned = stats::median(abs(d[[this_dIF]]), na.rm = TRUE),
    median_abs_dIF_other = stats::median(abs(d[[other_dIF]]), na.rm = TRUE),
    pct_same_direction_in_other = 100 * mean(sign(d[[this_dIF]]) == sign(d[[other_dIF]]), na.rm = TRUE),
    pct_passing_min_abs_dif_in_other = 100 * mean(abs(d[[other_dIF]]) >= min_abs_dif, na.rm = TRUE)
  )
}
concordance_one_way <- bind_rows(
  one_way_concordance(iso_overlap, "T_switching", "T_dIF", "U_dIF", "significant in T -> measured in U"),
  one_way_concordance(iso_overlap, "U_switching", "U_dIF", "T_dIF", "significant in U -> measured in T")
)
write_out(concordance_one_way, "ut_T_vs_U_one_way_concordance")

jaccard_genes <- {
  a <- sum(gene_overlap$overlap_class == "Shared (T and U)")
  b <- sum(gene_overlap$overlap_class == "T-only")
  c <- sum(gene_overlap$overlap_class == "U-only")
  if ((a + b + c) > 0L) a / (a + b + c) else NA_real_
}

shared_iso_both_sig <- iso_overlap |>
  filter(.data$T_switching, .data$U_switching, .data$present_in_both)

cor_dIF <- if (nrow(shared_iso_both_sig) >= 3L) {
  suppressWarnings(stats::cor.test(shared_iso_both_sig$T_dIF, shared_iso_both_sig$U_dIF, method = "spearman"))
} else {
  NULL
}

same_dir_n <- sum(shared_iso_both_sig$same_direction %in% TRUE, na.rm = TRUE)
opp_dir_n <- sum(shared_iso_both_sig$same_direction %in% FALSE, na.rm = TRUE)
binom_dir <- if ((same_dir_n + opp_dir_n) > 0L) {
  stats::binom.test(same_dir_n, same_dir_n + opp_dir_n, p = 0.5)
} else {
  NULL
}

# Compare |dIF| distributions: T-only switching isoforms vs U-only vs shared.
#
# The Shared group used to be summarised with pmax(T_abs_dIF, U_abs_dIF) while T-only
# and U-only used a single dataset's value. A maximum of two draws is upward-biased
# relative to one draw, so "Shared isoforms have larger |dIF|" was partly an artefact
# of the summary rather than a finding. Each isoform is now measured in exactly one
# dataset: the one it was called significant in (Shared uses T, and the paired T-vs-U
# comparison for those isoforms is covered by delta_dIF below).
eff_cmp <- iso_focus |>
  mutate(
    effect_group = case_when(
      .data$T_switching & .data$U_switching ~ "Shared",
      .data$T_switching & !.data$U_switching ~ "T-only",
      .data$U_switching & !.data$T_switching ~ "U-only",
      TRUE ~ "Other"
    ),
    effect_abs = dplyr::case_when(
      .data$effect_group == "Shared" ~ .data$T_abs_dIF,
      .data$effect_group == "T-only" ~ .data$T_abs_dIF,
      .data$effect_group == "U-only" ~ .data$U_abs_dIF,
      TRUE ~ NA_real_
    ),
    effect_measured_in = dplyr::case_when(
      .data$effect_group %in% c("Shared", "T-only") ~ "T_UT",
      .data$effect_group == "U-only" ~ "U_UT",
      TRUE ~ NA_character_
    )
  ) |>
  filter(.data$effect_group != "Other", is.finite(.data$effect_abs))

kw <- if (n_distinct(eff_cmp$effect_group) >= 2L && nrow(eff_cmp) > 3L) {
  stats::kruskal.test(effect_abs ~ effect_group, data = eff_cmp)
} else {
  NULL
}

# Pairwise Wilcoxon for shared-present isoforms: |U_dIF| vs |T_dIF|
paired_present <- iso_overlap |>
  filter(.data$present_in_both, is.finite(.data$T_dIF), is.finite(.data$U_dIF))
wilcox_abs <- if (nrow(paired_present) > 3L) {
  stats::wilcox.test(abs(paired_present$U_dIF), abs(paired_present$T_dIF), paired = TRUE)
} else {
  NULL
}
wilcox_delta <- if (nrow(paired_present) > 3L) {
  stats::wilcox.test(paired_present$delta_dIF, mu = 0)
} else {
  NULL
}

stats_summary <- tibble(
  isoform_q_cutoff = iso_q_cutoff,
  gene_q_cutoff = gene_q_cutoff,
  min_abs_dif_cutoff = min_abs_dif,
  n_isoforms_T = nrow(t_iso),
  n_isoforms_U = nrow(u_iso),
  n_switching_isoforms_T = sum(t_iso$is_switching),
  n_switching_isoforms_U = sum(u_iso$is_switching),
  n_switching_genes_T = nrow(t_genes),
  n_switching_genes_U = nrow(u_genes),
  n_shared_isoforms_any = sum(iso_overlap$present_in_both, na.rm = TRUE),
  n_shared_significant_isoforms = nrow(shared_iso_both_sig),
  n_genes_shared = sum(gene_overlap$overlap_class == "Shared (T and U)"),
  n_genes_T_only = sum(gene_overlap$overlap_class == "T-only"),
  n_genes_U_only = sum(gene_overlap$overlap_class == "U-only"),
  jaccard_gene_overlap = jaccard_genes,
  # Retention ceiling: gene overlap is bounded by how many symbols survive in both
  # saved objects, not by biology.
  n_symbols_retained_in_both_objects = length(bg_both),
  # Effect sizes first -- with thousands of matched isoforms the p-values below are
  # driven by n, and the medians are the interpretable quantities.
  median_delta_dIF_matched = suppressWarnings(median(paired_present$delta_dIF, na.rm = TRUE)),
  median_abs_delta_dIF_matched = suppressWarnings(median(abs(paired_present$delta_dIF), na.rm = TRUE)),
  n_matched_isoforms_both_finite = nrow(paired_present),
  # One-way concordance: condition on significance in one dataset, measure in the other.
  pct_T_sig_same_direction_in_U = concordance_one_way$pct_same_direction_in_other[[1L]],
  pct_T_sig_passing_threshold_in_U = concordance_one_way$pct_passing_min_abs_dif_in_other[[1L]],
  pct_U_sig_same_direction_in_T = concordance_one_way$pct_same_direction_in_other[[2L]],
  pct_U_sig_passing_threshold_in_T = concordance_one_way$pct_passing_min_abs_dif_in_other[[2L]],
  # The next four are conditioned on significance in BOTH datasets. Selecting on the
  # outcome in both and then measuring agreement manufactures agreement -- these are
  # descriptive of the shared set only, not evidence of concordance.
  selconf_spearman_rho_shared_sig_dIF = if (!is.null(cor_dIF)) unname(cor_dIF$estimate) else NA_real_,
  selconf_n_shared_sig_same_direction = same_dir_n,
  selconf_n_shared_sig_opposite_direction = opp_dir_n,
  selconf_binom_p_same_direction = if (!is.null(binom_dir)) binom_dir$p.value else NA_real_,
  kruskal_p_abs_dIF_by_class = if (!is.null(kw)) kw$p.value else NA_real_,
  wilcox_paired_abs_dIF_U_vs_T_p = if (!is.null(wilcox_abs)) wilcox_abs$p.value else NA_real_,
  wilcox_delta_dIF_vs_0_p = if (!is.null(wilcox_delta)) wilcox_delta$p.value else NA_real_,
  inputs_pre_reduced_to_switching_genes = TRUE,
  # Background statistics are valid only via the unfiltered context layer.
  unfiltered_context_available = has_context,
  n_genes_tested_in_both = if (!is.null(background_enrichment)) {
    background_enrichment$n_genes_tested_in_both
  } else {
    NA_integer_
  },
  background_odds_ratio = if (!is.null(background_enrichment)) {
    background_enrichment$odds_ratio
  } else {
    NA_real_
  },
  background_p_value = if (!is.null(background_enrichment)) {
    background_enrichment$p_value
  } else {
    NA_real_
  }
)
write_out(stats_summary, "ut_T_vs_U_stats_summary")

# Descriptive 2x2 over symbols retained in both objects. Reported as counts only --
# no test is run on it (see the note above the retention accounting).
retention_contingency <- data.frame(
  switching_in_T = c("FALSE", "TRUE"),
  switching_in_U_FALSE = c(
    sum(!in_t_sig_both & !in_u_sig_both),
    sum(in_t_sig_both & !in_u_sig_both)
  ),
  switching_in_U_TRUE = c(
    sum(!in_t_sig_both & in_u_sig_both),
    sum(in_t_sig_both & in_u_sig_both)
  ),
  stringsAsFactors = FALSE,
  row.names = NULL
)
write_out(retention_contingency, "ut_T_vs_U_retention_contingency")

# ---- Figures: overlap highlight ----
class_colors <- c(
  "Shared (T and U)" = "#756bb1",
  "T-only" = "#3182bd",
  "U-only" = "#e6550d",
  "Shared significant" = "#756bb1",
  "T-significant / tested in U, not switching" = "#9ecae1",
  "T-significant / not detected in U" = "#08519c",
  "U-significant / tested in T, not switching" = "#fdd0a2",
  "U-significant / not detected in T" = "#a63603",
  "T-significant / retained in U, not significant" = "#9ecae1",
  "T-significant / not retained in U" = "#08519c",
  "U-significant / retained in T, not significant" = "#fdd0a2",
  "U-significant / not retained in T" = "#a63603",
  "Non-significant / other" = "grey70",
  "Shared" = "#756bb1",
  "Other" = "grey70"
)

p_gene_bar <- ggplot(gene_overlap_summary, aes(x = .data$overlap_class, y = .data$n_genes, fill = .data$overlap_class)) +
  geom_col(color = "grey25", linewidth = 0.2, show.legend = FALSE) +
  geom_text(aes(label = .data$n_genes), vjust = -0.3, size = 3.5) +
  scale_fill_manual(values = class_colors) +
  labs(
    title = "UT gene-level switching overlap: T (WT) vs U (knockout)",
    subtitle = paste(
      "Genes with >=1 significant isoform switch (cutoffs from config).",
      "T-only/U-only are dominated by genes absent from the other saved object."
    ),
    x = NULL, y = "Number of genes"
  ) +
  ylim(0, max(gene_overlap_summary$n_genes) * 1.15)
ggsave(file.path(fig_dir, "fig_ut_gene_overlap_counts.png"), p_gene_bar, width = 7.5, height = 4.5, dpi = 200)

iso_class_counts <- iso_focus |>
  count(.data$overlap_class, name = "n") |>
  arrange(desc(.data$n))
p_iso_bar <- ggplot(iso_class_counts, aes(x = reorder(.data$overlap_class, .data$n), y = .data$n, fill = .data$overlap_class)) +
  geom_col(color = "grey25", linewidth = 0.2, show.legend = FALSE) +
  coord_flip() +
  scale_fill_manual(values = class_colors, guide = "none") +
  labs(
    title = "UT isoform-level switching categories (T vs U)",
    x = NULL, y = "Number of isoforms"
  )
ggsave(file.path(fig_dir, "fig_ut_isoform_overlap_counts.png"), p_iso_bar, width = 8.5, height = 4.8, dpi = 200)

# Scatter of matched isoforms
scatter_df <- iso_overlap |>
  filter(.data$present_in_both, is.finite(.data$T_dIF), is.finite(.data$U_dIF)) |>
  mutate(
    point_class = .data$overlap_class
  )
p_scatter <- ggplot(scatter_df, aes(x = .data$T_dIF, y = .data$U_dIF, color = .data$point_class)) +
  geom_hline(yintercept = 0, linewidth = 0.3, color = "grey50") +
  geom_vline(xintercept = 0, linewidth = 0.3, color = "grey50") +
  geom_abline(slope = 1, intercept = 0, linetype = 2, color = "grey40") +
  geom_point(alpha = 0.65, size = 1.8) +
  scale_color_manual(values = class_colors, name = NULL) +
  labs(
    title = "Matched UT isoforms: dIF in T (WT) vs U (knockout)",
    subtitle = "Shared isoform_id within UT reference; dashed line = equal dIF",
    x = "T_UT dIF", y = "U_UT dIF"
  ) +
  coord_equal()
ggsave(file.path(fig_dir, "fig_ut_dIF_scatter_matched_isoforms.png"), p_scatter, width = 8, height = 6.5, dpi = 200)

p_delta <- ggplot(
  paired_present |> mutate(abs_delta = abs(.data$delta_dIF)),
  aes(x = .data$delta_dIF)
) +
  geom_histogram(bins = 40, fill = "#6a51a3", color = "white", linewidth = 0.15) +
  geom_vline(xintercept = 0, linetype = 2) +
  labs(
    title = "Delta dIF (U - T) for matched UT isoforms",
    x = "U_dIF - T_dIF", y = "Count"
  )
ggsave(file.path(fig_dir, "fig_ut_delta_dIF_histogram.png"), p_delta, width = 7.2, height = 4.2, dpi = 200)

p_eff <- ggplot(eff_cmp, aes(x = .data$effect_group, y = .data$effect_abs, fill = .data$effect_group)) +
  geom_violin(trim = TRUE, alpha = 0.55, color = "grey35") +
  geom_boxplot(width = 0.18, outlier.alpha = 0.25, alpha = 0.85) +
  scale_fill_manual(values = class_colors, guide = "none") +
  labs(
    title = "|dIF| among significant isoform classes",
    subtitle = "Each isoform measured in one dataset only (Shared and T-only in T, U-only in U)",
    x = NULL, y = "|dIF|"
  )
ggsave(file.path(fig_dir, "fig_ut_abs_dIF_by_overlap_class.png"), p_eff, width = 7, height = 4.5, dpi = 200)

# Top genes by class for a rank-style overview
top_by_class <- gene_overlap |>
  group_by(.data$overlap_class) |>
  slice_head(n = 15) |>
  ungroup() |>
  mutate(
    gene_label = .data$gene_name,
    sort_val = dplyr::coalesce(.data$max_abs_dif_either, 0)
  )
p_top <- ggplot(
  top_by_class,
  aes(x = reorder(.data$gene_label, .data$sort_val), y = .data$sort_val, fill = .data$overlap_class)
) +
  geom_col(color = "grey25", linewidth = 0.15, show.legend = FALSE) +
  coord_flip() +
  facet_wrap(~ .data$overlap_class, scales = "free_y", ncol = 3) +
  scale_fill_manual(values = class_colors) +
  labs(
    title = "Top genes by overlap class (max |dIF| in either dataset)",
    x = NULL, y = "max |dIF|"
  )
ggsave(file.path(fig_dir, "fig_ut_top_genes_by_overlap_class.png"), p_top, width = 11, height = 7, dpi = 200)

# ---- Gene explorer panels: isoform usage in T and U, before and after LPS ----
#
# The primary objects are reduced to significant switching genes, so a gene switching in
# T but not U has no U rows at all -- the old version of this panel could only draw the
# side that was significant, and could not distinguish "the knockout has this transcript
# but its ratio does not change" from "the transcript is not there". The unfiltered
# context tables written by 01 carry every tested isoform, so both genotypes can be drawn
# for any gene present in either.
#
# Isoform fraction is shown per condition (IF1 = LPS-, IF2 = LPS+) with an arrow between
# them, so a ratio shift is visible directly rather than collapsed into a single dIF bar.
# Per-replicate IF values are overlaid when available. Significance still comes from the
# primary scored tables -- the context layer is for display and classification only.

sample_condition <- function(nms) {
  # Replicate columns are named like T1_minus_S15 / T1_plus_S16.
  ifelse(grepl("_minus", nms, fixed = TRUE), "LPS-",
    ifelse(grepl("_plus", nms, fixed = TRUE), "LPS+", NA_character_)
  )
}

rep_points_for <- function(rep_tbl, isoform_ids) {
  if (is.null(rep_tbl) || !length(isoform_ids)) return(NULL)
  r <- rep_tbl[rep_tbl$isoform_id %in% isoform_ids, , drop = FALSE]
  if (!nrow(r)) return(NULL)
  value_cols <- setdiff(names(r), "isoform_id")
  cond <- sample_condition(value_cols)
  keep <- !is.na(cond)
  if (!any(keep)) return(NULL)
  value_cols <- value_cols[keep]
  cond <- cond[keep]
  bind_rows(lapply(seq_along(value_cols), function(i) {
    tibble(
      isoform_id = as.character(r$isoform_id),
      condition = cond[[i]],
      IF = suppressWarnings(as.numeric(r[[value_cols[[i]]]]))
    )
  }))
}

#' Pull one gene's isoforms from a context table, falling back to the scored primary
#' table when no context is configured for that dataset.
gene_rows_for <- function(gene_id, ctx, scored, ds_label) {
  if (!is.null(ctx) && !is.na(gene_id)) {
    rows <- ctx[as.character(ctx$gene_id) == as.character(gene_id), , drop = FALSE]
    if (nrow(rows)) {
      return(tibble(
        isoform_id = as.character(rows$isoform_id),
        dataset = ds_label,
        IF1 = suppressWarnings(as.numeric(rows$IF1)),
        IF2 = suppressWarnings(as.numeric(rows$IF2)),
        dIF = suppressWarnings(as.numeric(rows$dIF)),
        source = "context"
      ))
    }
  }
  rows <- scored[as.character(scored$gene_id) == as.character(gene_id), , drop = FALSE]
  if (!nrow(rows)) return(NULL)
  tibble(
    isoform_id = as.character(rows$isoform_id),
    dataset = ds_label,
    IF1 = if ("IF1" %in% names(rows)) suppressWarnings(as.numeric(rows$IF1)) else NA_real_,
    IF2 = if ("IF2" %in% names(rows)) suppressWarnings(as.numeric(rows$IF2)) else NA_real_,
    dIF = suppressWarnings(as.numeric(rows$dIF_n)),
    source = "primary"
  )
}

plot_gene_explorer <- function(gname, gene_id_t, gene_id_u, out_png) {
  gid <- dplyr::coalesce(as.character(gene_id_t), as.character(gene_id_u))
  if (is.na(gid)) return(FALSE)
  # Both datasets use the UT reference, so a gene carries the same gene_id in each.
  df <- bind_rows(
    gene_rows_for(gid, ctx_t, t_iso, "T_UT (WT)"),
    gene_rows_for(gid, ctx_u, u_iso, "U_UT (UBL5 KO)")
  )
  if (is.null(df) || !nrow(df)) return(FALSE)

  sig_map <- bind_rows(
    t_iso |> transmute(
      isoform_id = as.character(.data$isoform_id), dataset = "T_UT (WT)",
      switching = .data$is_switching, novel = .data$is_novel
    ),
    u_iso |> transmute(
      isoform_id = as.character(.data$isoform_id), dataset = "U_UT (UBL5 KO)",
      switching = .data$is_switching, novel = .data$is_novel
    )
  )
  df <- df |>
    left_join(sig_map, by = c("isoform_id", "dataset")) |>
    mutate(
      switching = .data$switching %in% TRUE,
      status = ifelse(.data$switching, "Switching (significant)", "Tested, not switching")
    )

  novel_ids <- unique(c(
    as.character(t_iso$isoform_id[t_iso$is_novel %in% TRUE]),
    as.character(u_iso$isoform_id[u_iso$is_novel %in% TRUE])
  ))
  # Rank isoforms by the largest usage shift seen in either genotype.
  ord <- df |>
    group_by(.data$isoform_id) |>
    summarize(m = max(abs(.data$dIF), na.rm = TRUE), .groups = "drop") |>
    arrange(desc(.data$m))
  keep <- utils::head(ord$isoform_id[is.finite(ord$m)], explorer_max_isoforms)
  if (!length(keep)) return(FALSE)
  df <- df |>
    filter(.data$isoform_id %in% keep) |>
    mutate(
      isoform_label = paste0(
        .data$isoform_id, ifelse(.data$isoform_id %in% novel_ids, " (PB)", "")
      )
    )
  lab_levels <- unique(df$isoform_label[order(match(df$isoform_id, keep))])
  df$isoform_label <- factor(df$isoform_label, levels = rev(lab_levels))

  long <- bind_rows(
    df |> transmute(.data$isoform_label, .data$dataset, .data$status, condition = "LPS-", IF = .data$IF1),
    df |> transmute(.data$isoform_label, .data$dataset, .data$status, condition = "LPS+", IF = .data$IF2)
  ) |>
    filter(is.finite(.data$IF))
  if (!nrow(long)) return(FALSE)

  reps <- bind_rows(
    {
      r <- rep_points_for(rep_t, keep)
      if (is.null(r)) NULL else mutate(r, dataset = "T_UT (WT)")
    },
    {
      r <- rep_points_for(rep_u, keep)
      if (is.null(r)) NULL else mutate(r, dataset = "U_UT (UBL5 KO)")
    }
  )
  if (!is.null(reps) && nrow(reps)) {
    reps <- reps |>
      left_join(
        df |> distinct(.data$isoform_id, .data$isoform_label),
        by = "isoform_id"
      ) |>
      filter(!is.na(.data$isoform_label), is.finite(.data$IF))
  }

  p <- ggplot(df, aes(y = .data$isoform_label)) +
    geom_segment(
      aes(x = .data$IF1, xend = .data$IF2, yend = .data$isoform_label, colour = .data$status),
      arrow = arrow(length = unit(0.10, "in"), type = "closed"),
      linewidth = 0.9, na.rm = TRUE
    ) +
    geom_point(
      data = long, aes(x = .data$IF, shape = .data$condition),
      size = 2.1, colour = "grey20", na.rm = TRUE
    ) +
    scale_shape_manual(values = c("LPS-" = 1, "LPS+" = 16), name = "Condition mean") +
    scale_colour_manual(
      values = c("Switching (significant)" = "#1b7837", "Tested, not switching" = "#9e9ac8"),
      name = NULL
    ) +
    facet_wrap(~ .data$dataset, ncol = 2) +
    coord_cartesian(xlim = c(0, 1)) +
    labs(
      title = paste0("Gene explorer: ", gname),
      subtitle = paste0(
        "Isoform fraction before (open) and after (filled) LPS; arrow = shift. ",
        "Both genotypes shown regardless of significance."
      ),
      x = "Isoform fraction (IF)", y = NULL,
      caption = "PB = novel PacBio isoform. Small points = individual replicates where available."
    )
  if (!is.null(reps) && nrow(reps)) {
    p <- p + geom_point(
      data = reps, aes(x = .data$IF, y = .data$isoform_label),
      size = 0.8, alpha = 0.55, colour = "grey35",
      position = position_nudge(y = 0.22), na.rm = TRUE
    )
  }
  ggsave(out_png, p, width = 11, height = 6, units = "in", dpi = 160)
  TRUE
}

explorer_targets <- bind_rows(
  gene_overlap |> filter(.data$overlap_class == "Shared (T and U)") |> slice_head(n = n_gene_panels_per_class) |> mutate(pick_from = "shared"),
  gene_overlap |> filter(.data$overlap_class == "T-only") |> slice_head(n = n_gene_panels_per_class) |> mutate(pick_from = "T_only"),
  gene_overlap |> filter(.data$overlap_class == "U-only") |> slice_head(n = n_gene_panels_per_class) |> mutate(pick_from = "U_only")
) |>
  distinct(.data$gene_name, .keep_all = TRUE)

explorer_index <- list()
for (i in seq_len(nrow(explorer_targets))) {
  gname <- explorer_targets$gene_name[[i]]
  cls <- sanitize(explorer_targets$overlap_class[[i]])
  out_png <- file.path(
    gene_fig_dir,
    paste0("gene_explorer_", cls, "_", sprintf("%02d", i), "_", sanitize(gname), ".png")
  )
  ok <- tryCatch(
    plot_gene_explorer(
      gname,
      explorer_targets$T_gene_id[[i]],
      explorer_targets$U_gene_id[[i]],
      out_png
    ),
    error = function(e) {
      message("  gene explorer failed for ", gname, ": ", conditionMessage(e))
      FALSE
    }
  )
  explorer_index[[length(explorer_index) + 1L]] <- tibble(
    gene_name = gname,
    overlap_class = explorer_targets$overlap_class[[i]],
    figure = if (ok) basename(out_png) else NA_character_,
    T_max_abs_dif = explorer_targets$T_max_abs_dif[[i]],
    U_max_abs_dif = explorer_targets$U_max_abs_dif[[i]]
  )
}
explorer_index_df <- bind_rows(explorer_index)
write_out(explorer_index_df, "ut_gene_explorer_index")

# Full explorer table (all overlapping genes) for report browsing
write_out(
  gene_overlap |>
    select(
      gene_name, overlap_class, T_gene_id, U_gene_id,
      T_n_switching, U_n_switching, T_max_abs_dif, U_max_abs_dif,
      delta_max_abs_dif, T_min_q, U_min_q, T_novel_involved, U_novel_involved,
      T_opposing, U_opposing, T_rank, U_rank
    ),
  "ut_gene_explorer_table"
)

# ---- Switch plots for showcase genes ----
message("Loading ISA objects for UT switchPlots ...")
isa_t <- tryCatch(
  load_isa_input(resolve_path(t_pack$ds$isa_path, root = root), object_name = t_pack$ds$object_name),
  error = function(e) {
    warning("Could not load T_UT ISA: ", conditionMessage(e))
    NULL
  }
)
isa_u <- tryCatch(
  load_isa_input(resolve_path(u_pack$ds$isa_path, root = root), object_name = u_pack$ds$object_name),
  error = function(e) {
    warning("Could not load U_UT ISA: ", conditionMessage(e))
    NULL
  }
)
if (!is.null(isa_t)) isa_t <- annotate_isa_for_plotting(isa_t, t_iso)
if (!is.null(isa_u)) isa_u <- annotate_isa_for_plotting(isa_u, u_iso)

# Without the ISA objects every switchPlot attempt fails, and the loop below unlinks
# the destination on failure -- which used to silently delete good plots from a
# previous run whenever the external drive was not mounted.
skip_switch_plots <- is.null(isa_t) && is.null(isa_u)
if (skip_switch_plots) {
  warning(
    "ISA objects unavailable (is the external drive mounted?). ",
    "Keeping existing switch plots and skipping switchPlot generation."
  )
}

showcase <- bind_rows(
  gene_overlap |> filter(.data$overlap_class == "Shared (T and U)") |> slice_head(n = n_switch_plots_per_class),
  gene_overlap |> filter(.data$overlap_class == "T-only") |> slice_head(n = n_switch_plots_per_class),
  gene_overlap |> filter(.data$overlap_class == "U-only") |> slice_head(n = n_switch_plots_per_class)
) |>
  distinct(.data$gene_name, .keep_all = TRUE)

switch_index <- list()
for (i in seq_len(if (skip_switch_plots) 0L else nrow(showcase))) {
  gname <- showcase$gene_name[[i]]
  cls <- sanitize(showcase$overlap_class[[i]])
  # Prefer plotting in the dataset(s) where the gene is significant
  plot_sets <- character()
  if (isTRUE(showcase$significant_in_T[[i]])) plot_sets <- c(plot_sets, "T")
  if (isTRUE(showcase$significant_in_U[[i]])) plot_sets <- c(plot_sets, "U")
  for (ps in plot_sets) {
    isa_obj <- if (ps == "T") isa_t else isa_u
    gid <- if (ps == "T") showcase$T_gene_id[[i]] else showcase$U_gene_id[[i]]
    out_png <- file.path(
      switch_fig_dir,
      paste0("switch_", ps, "_", cls, "_", sprintf("%02d", i), "_", sanitize(gname), ".png")
    )
    ok <- make_switch_plot(isa_obj, gid, gname, out_png)
    if (!ok && file.exists(out_png)) unlink(out_png)
    switch_index[[length(switch_index) + 1L]] <- tibble(
      gene_name = gname,
      overlap_class = showcase$overlap_class[[i]],
      dataset = paste0(ps, "_UT"),
      figure = if (ok) basename(out_png) else NA_character_,
      ok = ok
    )
    if (ok) message("  switchPlot: ", basename(out_png))
  }
}
if (!skip_switch_plots) {
  write_out(bind_rows(switch_index), "ut_switch_plot_index")
} else {
  message("Skipped ut_switch_plot_index (no ISA objects; existing index left in place).")
}

write_run_manifest("04_comparison_T_vs_U.R", cfg, root)
message("04_comparison_T_vs_U.R: done")
print(stats_summary)
print(gene_overlap_summary)
