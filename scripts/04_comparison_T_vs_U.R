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

top_n_genes <- 30L
n_gene_panels_per_class <- 12L
n_switch_plots_per_class <- 6L

sanitize <- function(x) gsub("[^A-Za-z0-9_]+", "_", as.character(x), perl = TRUE)
min_finite <- function(v) {
  if (!length(v)) return(NA_real_)
  m <- suppressWarnings(min(v, na.rm = TRUE))
  if (is.infinite(m)) NA_real_ else m
}
max_abs_finite <- function(v) {
  w <- abs(as.numeric(v))
  w <- w[is.finite(w)]
  if (!length(w)) NA_real_ else max(w)
}
is_real_gene_symbol <- function(x) {
  x <- as.character(x)
  !is.na(x) & nzchar(x) & !grepl("^(XLOC_|ENS[GTFP]\\d)", x, perl = TRUE)
}
pick_gene_symbol <- function(names) {
  names <- unique(as.character(names))
  names <- names[!is.na(names) & nzchar(names)]
  if (!length(names)) return(NA_character_)
  good <- names[is_real_gene_symbol(names)]
  if (length(good)) good[[1L]] else names[[1L]]
}

theme_set(
  theme_bw(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold", size = 12),
      strip.text = element_text(face = "bold")
    )
)

score_isoforms <- function(iso_tbl, dataset_key, dataset_label) {
  z <- as_tibble(iso_tbl)
  if (!all(c("gene_id", "isoform_id") %in% names(z))) {
    stop("Need gene_id and isoform_id for ", dataset_key)
  }
  if (!"gene_name" %in% names(z)) z$gene_name <- NA_character_
  if (!"dIF" %in% names(z)) z$dIF <- NA_real_
  if (!"is_novel_pacbio" %in% names(z)) z$is_novel_pacbio <- NA
  if (!"oId" %in% names(z)) z$oId <- NA_character_
  if (!"class_code" %in% names(z)) z$class_code <- NA_character_
  if (!"isoform_switch_q_value" %in% names(z) && !"gene_switch_q_value" %in% names(z)) {
    stop("Need isoform/gene switch q for ", dataset_key)
  }

  z <- z |>
    group_by(.data$gene_id) |>
    mutate(gene_name = pick_gene_symbol(.data$gene_name)) |>
    ungroup()

  z$q_i <- if ("isoform_switch_q_value" %in% names(z)) {
    as.numeric(z$isoform_switch_q_value)
  } else {
    as.numeric(z$gene_switch_q_value)
  }
  z$q_g <- if ("gene_switch_q_value" %in% names(z)) as.numeric(z$gene_switch_q_value) else NA_real_
  z$dIF_n <- as.numeric(z$dIF)
  z$abs_dIF <- abs(z$dIF_n)
  z$is_novel <- z$is_novel_pacbio %in% TRUE
  z$is_switching <- is.finite(z$q_i) & z$q_i < iso_q_cutoff &
    (is.na(z$q_g) | (is.finite(z$q_g) & z$q_g < gene_q_cutoff)) &
    (is.na(z$dIF_n) | abs(z$dIF_n) >= min_abs_dif)
  z$is_switching <- replace(z$is_switching, is.na(z$is_switching), FALSE)
  z$dataset_key <- dataset_key
  z$dataset_label <- dataset_label
  z
}

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

write_out <- function(x, stem) {
  csv <- file.path(results_tables, paste0(stem, ".csv"))
  rds <- file.path(results_tables, paste0(stem, ".rds"))
  utils::write.csv(x, csv, row.names = FALSE, fileEncoding = "UTF-8", na = "")
  saveRDS(x, rds, compress = "xz")
  message("Wrote: ", csv)
}

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
    ok <- tryCatch({
      grDevices::png(filename = out_png, width = 2400, height = 1500, res = 180)
      on.exit({
        if (grDevices::dev.cur() > 1L) grDevices::dev.off()
      }, add = TRUE)
      do.call(switch_fn, list(switchAnalyzeRlist = isa_obj, gene = g, plotTopology = FALSE))
      TRUE
    }, error = function(e) FALSE)
    if (ok && file.exists(out_png) && isTRUE(file.info(out_png)$size > 0)) return(TRUE)
  }
  FALSE
}

# ---- Load and score UT datasets ----
message("Loading T_UT and U_UT ...")
t_pack <- load_processed("T_UT")
u_pack <- load_processed("U_UT")
t_iso <- score_isoforms(t_pack$iso, "T_UT", t_pack$label)
u_iso <- score_isoforms(u_pack$iso, "U_UT", u_pack$label)

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
    T_switching = .data$T_switching %in% TRUE,
    U_switching = .data$U_switching %in% TRUE,
    overlap_class = case_when(
      .data$T_switching & .data$U_switching ~ "Shared significant",
      .data$T_switching & !.data$U_switching & .data$present_in_U ~ "T-significant / U-present non-sig",
      .data$T_switching & !.data$present_in_U ~ "T-significant / missing in U",
      .data$U_switching & !.data$T_switching & .data$present_in_T ~ "U-significant / T-present non-sig",
      .data$U_switching & !.data$present_in_T ~ "U-significant / missing in T",
      TRUE ~ "Non-significant / other"
    ),
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

write_out(iso_overlap, "ut_isoform_overlap_T_vs_U_all")

iso_focus <- iso_overlap |>
  filter(.data$T_switching | .data$U_switching) |>
  arrange(desc(pmax(.data$T_abs_dIF, .data$U_abs_dIF, na.rm = TRUE)))
write_out(iso_focus, "ut_isoform_overlap_T_vs_U_significant")

# ---- Statistics ----
# Background for co-occurrence: gene symbols observed in BOTH UT datasets.
# Keep a full 2x2 even when some cells are zero.
bg_t <- unique(t_iso$gene_name[is_real_gene_symbol(t_iso$gene_name)])
bg_u <- unique(u_iso$gene_name[is_real_gene_symbol(u_iso$gene_name)])
bg_both <- intersect(bg_t, bg_u)
in_t_sig_both <- bg_both %in% t_gene_key$gene_name
in_u_sig_both <- bg_both %in% u_gene_key$gene_name
tab_both <- matrix(
  c(
    sum(!in_t_sig_both & !in_u_sig_both),
    sum(!in_t_sig_both & in_u_sig_both),
    sum(in_t_sig_both & !in_u_sig_both),
    sum(in_t_sig_both & in_u_sig_both)
  ),
  nrow = 2,
  byrow = TRUE,
  dimnames = list(
    T_sig = c("FALSE", "TRUE"),
    U_sig = c("FALSE", "TRUE")
  )
)
fisher_gene <- tryCatch({
  ft <- fisher.test(tab_both)
  list(
    table = tab_both,
    odds_ratio = unname(ft$estimate),
    p_value = ft$p.value,
    conf_low = ft$conf.int[[1]],
    conf_high = ft$conf.int[[2]],
    n_background_both_present = length(bg_both)
  )
}, error = function(e) {
  list(
    table = tab_both,
    odds_ratio = NA_real_,
    p_value = NA_real_,
    conf_low = NA_real_,
    conf_high = NA_real_,
    n_background_both_present = length(bg_both),
    error = conditionMessage(e)
  )
})

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

# Compare |dIF| distributions: T-only switching isoforms vs U-only vs shared
eff_cmp <- iso_focus |>
  mutate(
    effect_group = case_when(
      .data$T_switching & .data$U_switching ~ "Shared",
      .data$T_switching & !.data$U_switching ~ "T-only",
      .data$U_switching & !.data$T_switching ~ "U-only",
      TRUE ~ "Other"
    ),
    effect_abs = dplyr::case_when(
      .data$effect_group == "Shared" ~ pmax(.data$T_abs_dIF, .data$U_abs_dIF, na.rm = TRUE),
      .data$effect_group == "T-only" ~ .data$T_abs_dIF,
      .data$effect_group == "U-only" ~ .data$U_abs_dIF,
      TRUE ~ NA_real_
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
  n_background_genes_present_in_both = fisher_gene$n_background_both_present %||% length(bg_both),
  fisher_odds_ratio_gene_sig_both_present = fisher_gene$odds_ratio %||% NA_real_,
  fisher_p_gene_sig_both_present = fisher_gene$p_value %||% NA_real_,
  spearman_rho_shared_sig_dIF = if (!is.null(cor_dIF)) unname(cor_dIF$estimate) else NA_real_,
  spearman_p_shared_sig_dIF = if (!is.null(cor_dIF)) cor_dIF$p.value else NA_real_,
  n_shared_sig_same_direction = same_dir_n,
  n_shared_sig_opposite_direction = opp_dir_n,
  binom_p_same_direction = if (!is.null(binom_dir)) binom_dir$p.value else NA_real_,
  kruskal_p_abs_dIF_by_class = if (!is.null(kw)) kw$p.value else NA_real_,
  wilcox_paired_abs_dIF_U_vs_T_p = if (!is.null(wilcox_abs)) wilcox_abs$p.value else NA_real_,
  wilcox_delta_dIF_vs_0_p = if (!is.null(wilcox_delta)) wilcox_delta$p.value else NA_real_,
  median_delta_dIF_shared_present = suppressWarnings(median(paired_present$delta_dIF, na.rm = TRUE)),
  median_abs_delta_dIF_shared_present = suppressWarnings(median(abs(paired_present$delta_dIF), na.rm = TRUE))
)
write_out(stats_summary, "ut_T_vs_U_stats_summary")

# Fisher table for report (explicit 2x2)
fisher_tbl <- data.frame(
  T_sig = c("FALSE", "TRUE"),
  U_sig_FALSE = tab_both[, "FALSE"],
  U_sig_TRUE = tab_both[, "TRUE"],
  stringsAsFactors = FALSE,
  row.names = NULL
)
write_out(fisher_tbl, "ut_T_vs_U_fisher_gene_contingency")

# ---- Figures: overlap highlight ----
class_colors <- c(
  "Shared (T and U)" = "#756bb1",
  "T-only" = "#3182bd",
  "U-only" = "#e6550d",
  "Shared significant" = "#756bb1",
  "T-significant / U-present non-sig" = "#9ecae1",
  "T-significant / missing in U" = "#08519c",
  "U-significant / T-present non-sig" = "#fdd0a2",
  "U-significant / missing in T" = "#a63603",
  "Shared" = "#756bb1",
  "Other" = "grey70"
)

p_gene_bar <- ggplot(gene_overlap_summary, aes(x = .data$overlap_class, y = .data$n_genes, fill = .data$overlap_class)) +
  geom_col(color = "grey25", linewidth = 0.2, show.legend = FALSE) +
  geom_text(aes(label = .data$n_genes), vjust = -0.3, size = 3.5) +
  scale_fill_manual(values = class_colors) +
  labs(
    title = "UT gene-level switching overlap: T (WT) vs U (knockout)",
    subtitle = "Genes with >=1 significant isoform switch (q cutoff from config)",
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
    point_class = case_when(
      .data$T_switching & .data$U_switching ~ "Shared significant",
      .data$T_switching & !.data$U_switching ~ "T-significant / U-present non-sig",
      .data$U_switching & !.data$T_switching ~ "U-significant / T-present non-sig",
      TRUE ~ "Non-significant / other"
    )
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

# ---- Gene explorer panels: side-by-side isoform dIF for T vs U ----
plot_gene_explorer <- function(gname, out_png) {
  t_rows <- t_iso |>
    filter(.data$gene_name == gname) |>
    transmute(
      isoform_id = as.character(.data$isoform_id),
      dataset = "T_UT",
      dIF = .data$dIF_n,
      q = .data$q_i,
      switching = .data$is_switching,
      novel = .data$is_novel
    )
  u_rows <- u_iso |>
    filter(.data$gene_name == gname) |>
    transmute(
      isoform_id = as.character(.data$isoform_id),
      dataset = "U_UT",
      dIF = .data$dIF_n,
      q = .data$q_i,
      switching = .data$is_switching,
      novel = .data$is_novel
    )
  df <- bind_rows(t_rows, u_rows)
  if (!nrow(df)) return(FALSE)
  df <- df |>
    mutate(
      isoform_label = paste0(
        .data$isoform_id,
        ifelse(.data$novel %in% TRUE, " (PB)", "")
      ),
      status = ifelse(.data$switching %in% TRUE, "Significant", "Not significant")
    )
  # Keep isoforms that exist in either and order by max |dIF|
  ord <- df |>
    group_by(.data$isoform_label) |>
    summarize(m = max(abs(.data$dIF), na.rm = TRUE), .groups = "drop") |>
    arrange(desc(.data$m))
  keep <- utils::head(ord$isoform_label, 12L)
  df <- df |>
    filter(.data$isoform_label %in% keep) |>
    mutate(isoform_label = factor(.data$isoform_label, levels = rev(keep)))

  p <- ggplot(df, aes(x = .data$isoform_label, y = .data$dIF, fill = .data$status)) +
    geom_col(width = 0.7, color = "grey25", linewidth = 0.15) +
    geom_hline(yintercept = 0, linewidth = 0.3) +
    coord_flip() +
    facet_wrap(~ .data$dataset, ncol = 2) +
    scale_fill_manual(values = c("Significant" = "#1b7837", "Not significant" = "#bdbdbd")) +
    labs(
      title = paste0("Gene explorer: ", gname),
      subtitle = "Isoform dIF in T_UT (WT) vs U_UT (knockout); PB = novel PacBio",
      x = NULL, y = "dIF", fill = NULL
    )
  ggsave(out_png, p, width = 10, height = 5.5, units = "in", dpi = 160)
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
  ok <- tryCatch(plot_gene_explorer(gname, out_png), error = function(e) FALSE)
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

showcase <- bind_rows(
  gene_overlap |> filter(.data$overlap_class == "Shared (T and U)") |> slice_head(n = n_switch_plots_per_class),
  gene_overlap |> filter(.data$overlap_class == "T-only") |> slice_head(n = n_switch_plots_per_class),
  gene_overlap |> filter(.data$overlap_class == "U-only") |> slice_head(n = n_switch_plots_per_class)
) |>
  distinct(.data$gene_name, .keep_all = TRUE)

switch_index <- list()
for (i in seq_len(nrow(showcase))) {
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
switch_index_df <- bind_rows(switch_index)
write_out(switch_index_df, "ut_switch_plot_index")

message("04_comparison_T_vs_U.R: done")
print(stats_summary)
print(gene_overlap_summary)
