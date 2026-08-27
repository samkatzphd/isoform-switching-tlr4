# Shared helpers for the isoform switching pipeline
# (SwitchAnalyzeR / IsoformSwitchAnalyzeR outputs, isoform -> gene collapse, novel flags)

# ---- Project / config (no hidden state: explicit paths) ---------------------------------

#' Resolve project root: env ISOFORM_PROJECT_ROOT, else walk up for config/config.yml
find_project_root <- function(start = getwd()) {
  if (nzchar(Sys.getenv("ISOFORM_PROJECT_ROOT", ""))) {
    return(normalizePath(Sys.getenv("ISOFORM_PROJECT_ROOT"), mustWork = TRUE))
  }
  p <- normalizePath(start, mustWork = FALSE)
  for (i in seq_len(20L)) {
    yml_file <- file.path(p, "config", "config.yml")
    if (file.exists(yml_file)) {
      return(normalizePath(p, mustWork = TRUE))
    }
    parent <- dirname(p)
    if (identical(p, parent)) {
      break
    }
    p <- parent
  }
  stop(
    "Could not locate config/config.yml. Run from the project root, or set ISOFORM_PROJECT_ROOT."
  )
}

is_absolute_path <- function(p) {
  nchar(p) > 0L && grepl("^(/|[A-Za-z]:)", p, perl = TRUE)
}

#' Load YAML project config; path is relative to `root` unless absolute
load_yaml_config <- function(
    path = "config/config.yml",
    root = find_project_root()
) {
  yml_path <- if (is_absolute_path(path)) path else file.path(root, path)
  if (!file.exists(yml_path)) {
    stop("Config not found: ", yml_path)
  }
  yaml::read_yaml(yml_path)
}

#' Resolve a path in config: relative to project `root` unless already absolute
resolve_path <- function(p, root = find_project_root()) {
  if (is.null(p) || (is.character(p) && !any(nzchar(p), na.rm = TRUE))) {
    return(NULL)
  }
  if (is_absolute_path(p)) {
    return(normalizePath(p, mustWork = FALSE))
  }
  normalizePath(file.path(root, p), mustWork = FALSE)
}

#' Create directory (recursive) and return the path
ensure_dir <- function(dir_path) {
  if (!dir.exists(dir_path)) {
    dir.create(dir_path, recursive = TRUE, showWarnings = FALSE)
  }
  dir_path
}

#' Load a SwitchAnalyzeR object from .rds or .RData/.Rdata
#'
#' For .RData files, if `object_name` is NULL and exactly one object is found,
#' it is returned. If multiple objects exist, the first object inheriting from
#' list/S4/data.frame is selected, otherwise the first loaded object is used.
load_isa_input <- function(path, object_name = NULL) {
  if (!file.exists(path)) {
    stop("Input file not found: ", path)
  }
  low <- tolower(path)
  if (grepl("\\.rds$", low, perl = TRUE)) {
    return(readRDS(path))
  }
  if (grepl("\\.rdata$|\\.rda$", low, perl = TRUE)) {
    env <- new.env(parent = emptyenv())
    loaded <- load(path, envir = env)
    if (!length(loaded)) {
      stop("No objects found in RData file: ", path)
    }
    if (!is.null(object_name)) {
      if (!object_name %in% loaded) {
        stop(
          "Requested object_name '", object_name, "' not found in ", path,
          ". Loaded objects: ", paste(loaded, collapse = ", ")
        )
      }
      return(get(object_name, envir = env, inherits = FALSE))
    }
    if (length(loaded) == 1L) {
      return(get(loaded[[1L]], envir = env, inherits = FALSE))
    }
    vals <- lapply(loaded, function(nm) get(nm, envir = env, inherits = FALSE))
    good <- vapply(vals, function(v) is.list(v) || isS4(v) || is.data.frame(v), logical(1))
    pick <- if (any(good)) which(good)[1L] else 1L
    return(vals[[pick]])
  }
  stop("Unsupported input extension for ", path, ". Use .rds or .RData/.Rdata")
}

# Parse one key from GTF/GFF attribute string, e.g. key "oId"
extract_gtf_attr <- function(attr, key) {
  # Supports both:
  # - GTF style: key "value";
  # - GFF3 style: key=value;
  m_gtf <- regexec(
    paste0("(^|;[[:space:]]*)", key, "[[:space:]]+\"([^\"]+)\""),
    attr,
    perl = TRUE
  )
  r_gtf <- regmatches(attr, m_gtf)
  out <- vapply(r_gtf, function(x) {
    if (length(x) >= 3L) x[[3L]] else NA_character_
  }, character(1))
  miss <- is.na(out) | !nzchar(out)
  if (any(miss)) {
    m_gff3 <- regexec(
      paste0("(^|;)", key, "=([^;]+)"),
      attr[miss],
      perl = TRUE
    )
    r_gff3 <- regmatches(attr[miss], m_gff3)
    out2 <- vapply(r_gff3, function(x) {
      if (length(x) >= 3L) utils::URLdecode(x[[3L]]) else NA_character_
    }, character(1))
    out[miss] <- out2
  }
  out
}

# Build transcript-level mapping table from a GTF-like file.
# Keeps one row per transcript_id with columns useful for isoform annotation.
build_transcript_map_from_gtf <- function(gtf_path) {
  if (!file.exists(gtf_path)) {
    stop("GTF/GFF file not found: ", gtf_path)
  }
  x <- utils::read.delim(
    gtf_path,
    sep = "\t",
    header = FALSE,
    comment.char = "#",
    quote = "",
    stringsAsFactors = FALSE
  )
  if (ncol(x) < 9L) {
    stop("Expected 9-column GTF-like file: ", gtf_path)
  }
  names(x)[1:9] <- c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute")
  tx <- x[x$feature %in% c("transcript", "mRNA"), c("attribute"), drop = FALSE]
  if (nrow(tx) < 1L) {
    stop("No transcript records found in: ", gtf_path)
  }
  a <- tx$attribute
  map <- data.frame(
    transcript_id = extract_gtf_attr(a, "transcript_id"),
    gene_id_from_gtf = extract_gtf_attr(a, "gene_id"),
    gene_name_from_gtf = extract_gtf_attr(a, "gene_name"),
    oId = extract_gtf_attr(a, "oId"),
    cmp_ref = extract_gtf_attr(a, "cmp_ref"),
    class_code = extract_gtf_attr(a, "class_code"),
    stringsAsFactors = FALSE
  )
  map <- map[!is.na(map$transcript_id) & nzchar(map$transcript_id), , drop = FALSE]
  map <- map[!duplicated(map$transcript_id), , drop = FALSE]
  tibble::as_tibble(map, .name_repair = "unique")
}

# ---- SwitchAnalyzeR: isoform-level table extraction -----------------------------------

#' @description
#' Extracts the per-isoform table from a SwitchAnalyzeR/IsoformSwitchAnalyzeR
#' object. Tries, in order:
#' 1) Named `isoformFeatures` if the object is a plain list
#' 2) S4 slot `isoformFeatures` if available
#' 3) `x@.Data[[1]]` when the first list element is the main isoform tibble
#' 4) `x[[1]]` on a plain list when the first element is a data frame
#'
#' @param x A SwitchList / switch list / isoform-level tibble
#' @param .verbose If TRUE, print which branch was used
#' @return A tibble with one row per isoform
extract_isoform_features <- function(x, .verbose = FALSE) {
  if (is.null(x)) {
    stop("extract_isoform_features: input is NULL.")
  }
  if (is.data.frame(x) || tibble::is_tibble(x)) {
    if (isTRUE(.verbose)) {
      message("Input is already a data frame; coercing to tibble.")
    }
    return(tibble::as_tibble(x, .name_repair = "unique"))
  }
  if (is.list(x) && !isS4(x) && "isoformFeatures" %in% names(x) &&
    !is.null(x[["isoformFeatures"]])) {
    if (isTRUE(.verbose)) {
      message("Using x$isoformFeatures.")
    }
    return(tibble::as_tibble(x[["isoformFeatures"]], .name_repair = "unique"))
  }
  if (isS4(x) && "isoformFeatures" %in% methods::slotNames(x)) {
    y <- methods::slot(x, "isoformFeatures")
    if (isTRUE(.verbose)) {
      message("Using S4 slot 'isoformFeatures'.")
    }
    return(tibble::as_tibble(y, .name_repair = "unique"))
  }
  if (isS4(x)) {
    d <- tryCatch(x@.Data, error = function(e) NULL)
    if (is.list(d) && length(d) >= 1L) {
      first <- d[[1L]]
      if (is.data.frame(first) || tibble::is_tibble(first)) {
        if (isTRUE(.verbose)) {
          message("Using S4 @.Data[[1]] (data frame).")
        }
        return(tibble::as_tibble(first, .name_repair = "unique"))
      }
    }
  }
  if (is.list(x) && !isS4(x) && length(x) >= 1L) {
    first <- x[[1L]]
    if (is.data.frame(first) || tibble::is_tibble(first)) {
      if (isTRUE(.verbose)) {
        message("Using list x[[1]] (data frame).")
      }
      return(tibble::as_tibble(first, .name_repair = "unique"))
    }
  }
  stop(
    "Could not find isoform-level data: expected isoformFeatures, @.Data[[1]], or list [[1]] as tibble/df"
  )
}

# ---- Novel (PacBio) flags --------------------------------------------------------------

#' @description
#' Mark isoform rows with novel (PacBio-style) ids when the id column starts with
#' a prefix, usually \"PB\" (IsoformSwitchAnalyzeR convention for ENST vs. PacBio).
#'
#' @param col Name of the transcript id column (often `oId`)
#' @param out_name Name of the new logical column
tag_novel_isoforms <- function(
    df,
    col = "oId",
    prefix = "PB",
    out_name = "is_novel_pacbio"
) {
  if (!is.data.frame(df)) {
    stop("tag_novel_isoforms: `df` must be a data frame.")
  }
  if (!col %in% names(df)) {
    warning(
      "No column '", col, "': cannot tag novel isoforms; '", out_name, "' set to NA."
    )
    df[[out_name]] <- NA
    return(df)
  }
  s <- as.character(df[[col]])
  present <- !is.na(s) & s != "NA" & s != "" & s != "NaN"
  if (nchar(prefix) < 1L) {
    warning("Empty `prefix`: setting '", out_name, "' to NA for all rows.")
    df[[out_name]] <- NA
    return(df)
  }
  df[[out_name]] <- present & startsWith(s, prefix)
  df
}

`%||%` <- function(x, y) {
  if (is.null(x) || (length(x) == 0L)) y else x
}

# ---- Small shared utilities -----------------------------------------------------------

#' Filesystem-safe token (used for dataset labels and gene stubs in filenames)
sanitize <- function(x) gsub("[^A-Za-z0-9_]+", "_", as.character(x), perl = TRUE)

#' min() that returns NA rather than Inf on an all-NA/empty vector
min_finite <- function(v) {
  if (!length(v)) {
    return(NA_real_)
  }
  m <- suppressWarnings(min(v, na.rm = TRUE))
  if (is.infinite(m)) NA_real_ else m
}

#' max(|v|) that returns NA rather than -Inf on an all-NA/empty vector
max_abs_finite <- function(v) {
  w <- abs(as.numeric(v))
  w <- w[is.finite(w)]
  if (!length(w)) NA_real_ else max(w)
}

#' Read an analysis knob from config$analysis, falling back to `default`
analysis_param <- function(cfg, name, default) {
  v <- (cfg$analysis %||% list())[[name]]
  if (is.null(v) || (is.character(v) && !nzchar(v))) default else v
}

# ---- Gene symbol handling -------------------------------------------------------------
# Reference transcriptomes assign their own XLOC_* gene ids, and some GFF rows carry
# ENSG/ENST placeholders in gene_name. A "real" symbol is neither. Cross-dataset joins
# use symbols (ids are transcriptome-specific), so this rule must be identical
# everywhere -- it previously differed between scripts 03 and 04/06.

is_real_gene_symbol <- function(x) {
  x <- as.character(x)
  !is.na(x) & nzchar(x) & !grepl("^(XLOC_|ENS[GTFP]\\d)", x, perl = TRUE)
}

#' Pick one representative symbol when isoforms of a gene_id disagree
pick_gene_symbol <- function(names) {
  names <- unique(as.character(names))
  names <- names[!is.na(names) & nzchar(names)]
  if (!length(names)) {
    return(NA_character_)
  }
  good <- names[is_real_gene_symbol(names)]
  if (length(good)) good[[1L]] else names[[1L]]
}

#' Force one gene_name per gene_id.
#'
#' Without this, `group_by(gene_id, gene_name)` splits a gene across several rows
#' whenever its isoforms disagree on the symbol (11-18% of gene_ids in this project),
#' which inflates gene counts and makes any later join on gene_id fan out.
collapse_gene_symbols <- function(iso) {
  iso <- tibble::as_tibble(iso, .name_repair = "unique")
  if (!"gene_id" %in% names(iso)) {
    stop("collapse_gene_symbols: 'gene_id' column required.")
  }
  if (!"gene_name" %in% names(iso)) {
    iso$gene_name <- NA_character_
    return(iso)
  }
  iso |>
    dplyr::group_by(.data$gene_id) |>
    dplyr::mutate(gene_name = pick_gene_symbol(.data$gene_name)) |>
    dplyr::ungroup()
}

# ---- Canonical switching rule ---------------------------------------------------------

#' @description
#' Score isoforms against the config significance rule. This is the single definition
#' of "switching" for the whole pipeline; scripts 02/03/04/06 all call it rather than
#' re-implementing the comparison (they used to, and had drifted apart).
#'
#' An isoform is *switching* when isoform q < `isoform_q`, gene q < `gene_q` (only
#' applied when the column exists), and `|dIF| >= min_abs_dif`.
#'
#' NA handling: with `significance.require_finite_dif: true` (the default) an isoform
#' with a missing or non-finite dIF or gene q *fails* the filter. The original code
#' let those rows through, so an isoform with no measured effect size could be called
#' switching on the q-value alone. Set the flag to false to restore the old behaviour.
#'
#' @return the input tibble plus q_i, q_g, dIF_n, abs_dIF, is_novel, is_switching
score_isoforms <- function(iso, cfg, dataset_key = NA_character_,
                           dataset_label = NA_character_, collapse_symbols = TRUE) {
  z <- tibble::as_tibble(iso, .name_repair = "unique")
  if (!all(c("gene_id", "isoform_id") %in% names(z))) {
    stop("score_isoforms: need gene_id and isoform_id (dataset: ", dataset_key, ")")
  }
  has_iso_q <- "isoform_switch_q_value" %in% names(z)
  has_gene_q <- "gene_switch_q_value" %in% names(z)
  if (!has_iso_q && !has_gene_q) {
    stop(
      "score_isoforms: need isoform_switch_q_value or gene_switch_q_value (dataset: ",
      dataset_key, ")"
    )
  }
  for (col in c("gene_name", "oId", "class_code", "cmp_ref")) {
    if (!col %in% names(z)) z[[col]] <- NA_character_
  }
  if (!"dIF" %in% names(z)) z$dIF <- NA_real_
  if (!"is_novel_pacbio" %in% names(z)) z$is_novel_pacbio <- NA

  if (isTRUE(collapse_symbols)) {
    z <- collapse_gene_symbols(z)
  }

  sig <- cfg$significance
  if (is.null(sig)) {
    stop("config must define `significance` (isoform_q, gene_q, min_abs_dif).")
  }
  i_q <- as.numeric(sig$isoform_q %||% 0.05)
  g_q <- as.numeric(sig$gene_q %||% 0.05)
  m_dif <- as.numeric(sig$min_abs_dif %||% 0.0)
  strict_na <- isTRUE(sig$require_finite_dif %||% TRUE)
  min_gexpr <- suppressWarnings(as.numeric(sig$min_gene_expression %||% 0))
  if (is.na(min_gexpr)) min_gexpr <- 0

  z$q_i <- if (has_iso_q) as.numeric(z$isoform_switch_q_value) else as.numeric(z$gene_switch_q_value)
  z$q_g <- if (has_gene_q) as.numeric(z$gene_switch_q_value) else NA_real_
  z$dIF_n <- as.numeric(z$dIF)
  z$abs_dIF <- abs(z$dIF_n)
  z$is_novel <- z$is_novel_pacbio %in% TRUE

  ok_i <- is.finite(z$q_i) & z$q_i < i_q
  ok_g <- if (!has_gene_q) {
    rep(TRUE, nrow(z))
  } else if (strict_na) {
    is.finite(z$q_g) & z$q_g < g_q
  } else {
    is.na(z$q_g) | (is.finite(z$q_g) & z$q_g < g_q)
  }
  ok_d <- if (strict_na) {
    is.finite(z$dIF_n) & abs(z$dIF_n) >= m_dif
  } else {
    !is.finite(z$dIF_n) | abs(z$dIF_n) >= m_dif
  }

  # Gene-level abundance floor. dIF is isoform/gene, so it is the GENE's total expression
  # that sets how precisely the fraction can be estimated -- not the isoform's. Isoforms of
  # a gene below the floor are excluded from switching calls and flagged, so they stay
  # inspectable without entering downstream gene sets. See 02d for how the floor is derived.
  if (all(c("iso_value_1", "iso_value_2") %in% names(z))) {
    ge <- z |>
      dplyr::group_by(.data$gene_id) |>
      dplyr::mutate(
        gene_expression = sum(as.numeric(.data$iso_value_1), na.rm = TRUE) +
          sum(as.numeric(.data$iso_value_2), na.rm = TRUE)
      ) |>
      dplyr::ungroup()
    z$gene_expression <- ge$gene_expression
  } else if (all(c("gene_value_1", "gene_value_2") %in% names(z))) {
    z$gene_expression <- as.numeric(z$gene_value_1) + as.numeric(z$gene_value_2)
  } else {
    z$gene_expression <- NA_real_
  }
  z$low_expression <- if (min_gexpr > 0 && any(is.finite(z$gene_expression))) {
    is.finite(z$gene_expression) & z$gene_expression < min_gexpr
  } else {
    FALSE
  }
  ok_e <- if (min_gexpr > 0 && any(is.finite(z$gene_expression))) {
    !(z$low_expression %in% TRUE)
  } else {
    rep(TRUE, nrow(z))
  }

  z$is_switching <- ok_i & ok_g & ok_d & ok_e
  z$is_switching <- replace(z$is_switching, is.na(z$is_switching), FALSE)
  # Retained so the cost of the floor is always visible, not silently absorbed.
  z$is_switching_before_expr_floor <- replace(ok_i & ok_g & ok_d, is.na(ok_i & ok_g & ok_d), FALSE)
  z$dataset_key <- dataset_key
  z$dataset_label <- dataset_label
  z
}

#' @description
#' Return the table a script should score, preferring the unfiltered context.
#'
#' Significance and presence must come from the SAME object. Taking q-values from the
#' reduced object while taking presence/effect from the unfiltered one mixes two
#' multiple-testing universes: across 966 shared U_UT isoforms the dIF values are identical
#' but only 132 isoform q-values match, and 29 isoforms flip across q < 0.05 in one
#' direction or the other. That mismatch misclassified genes (CD86 was called T-only despite
#' passing every criterion in U). The context is preferred because its FDR correction spans
#' the genes actually tested.
#'
#' @return the context table when available, else the primary table; the source is recorded
#'   in the "scoring_source" attribute and should be reported by callers.
load_scoring_table <- function(processed_dir, label, quiet = FALSE) {
  ctx <- load_context_table(processed_dir, label, "features")
  if (!is.null(ctx)) {
    attr(ctx, "scoring_source") <- "unfiltered context"
    if (!quiet) {
      message("  [", label, "] scoring from unfiltered context (", nrow(ctx), " isoforms)")
    }
    return(ctx)
  }
  p <- file.path(processed_dir, paste0("isoformFeatures_", sanitize(label), ".rds"))
  if (!file.exists(p)) {
    return(NULL)
  }
  prim <- readRDS(p)
  attr(prim, "scoring_source") <- "reduced primary object"
  if (!quiet) {
    message(
      "  [", label, "] no context; scoring from the REDUCED object (", nrow(prim),
      " isoforms). Counts are not comparable with context-scored datasets."
    )
  }
  prim
}

# ---- Gene-level summary from isoform tibble ------------------------------------------------

#' @description
#' One row per gene_id (guaranteed: symbols are collapsed first). *novel_involved* is
#' TRUE if any switching isoform is flagged PacBio-novel.
#'
#' Two distinct effect-size columns are reported because conflating them is easy:
#'   - `max_abs_dif_all_isoforms` -- max |dIF| over every isoform of the gene
#'   - `max_abs_dif_switching`    -- max |dIF| over switching isoforms only (NA if none)
#' The old single `max_abs_dif` column was the former while scripts 03/04/06 ranked on
#' the latter under a near-identical name.
summarize_genes_from_isoform_table <- function(iso, cfg) {
  if (!tibble::is_tibble(iso) && !is.data.frame(iso)) {
    stop("summarize_genes_from_isoform_table: `iso` must be a tibble or data frame.")
  }
  z <- score_isoforms(iso, cfg)
  z |>
    dplyr::group_by(.data$gene_id, .data$gene_name) |>
    dplyr::summarize(
      n_isoforms = dplyr::n(),
      n_switching_isoforms = sum(.data$is_switching, na.rm = TRUE),
      n_novel_isoforms = sum(.data$is_novel, na.rm = TRUE),
      min_isoform_switch_q = min_finite(.data$q_i),
      min_gene_switch_q = min_finite(.data$q_g),
      max_abs_dif_all_isoforms = max_abs_finite(.data$dIF_n),
      max_abs_dif_switching = max_abs_finite(.data$dIF_n[.data$is_switching]),
      novel_involved = any(.data$is_switching & .data$is_novel, na.rm = TRUE),
      .groups = "drop"
    )
}

# ---- Output helpers -------------------------------------------------------------------

#' Write a table as .rds (always) and .csv (unless suppressed by config).
#'
#' `output.csv_twin: false` drops every CSV twin; `output.csv_twin_max_rows` drops it
#' above a row count. Pass `csv = FALSE` at a call site for wide isoform-level tables,
#' whose CSV twins were adding megabytes to git on every regeneration.
write_table_pair <- function(x, dir, stem, cfg = NULL, csv = NULL, quiet = FALSE) {
  rds <- file.path(dir, paste0(stem, ".rds"))
  saveRDS(x, rds, compress = "xz")
  out <- cfg$output %||% list()
  want_csv <- isTRUE(out$csv_twin %||% TRUE)
  if (!is.null(csv)) want_csv <- isTRUE(csv) && isTRUE(out$csv_twin %||% TRUE)
  max_rows <- suppressWarnings(as.numeric(out$csv_twin_max_rows %||% Inf))
  if (is.na(max_rows)) max_rows <- Inf
  csv <- file.path(dir, paste0(stem, ".csv"))
  if (want_csv && nrow(x) <= max_rows) {
    utils::write.csv(x, csv, row.names = FALSE, fileEncoding = "UTF-8", na = "")
    if (!quiet) message("Wrote: ", csv, " | ", basename(rds))
  } else {
    if (file.exists(csv)) unlink(csv)
    if (!quiet) message("Wrote: ", rds, " (csv suppressed by config)")
  }
  invisible(rds)
}

# ---- Unfiltered context layer ---------------------------------------------------------

#' @description
#' Extract a slim "context" view from an UNREDUCED switchAnalyzeRlist.
#'
#' The saved analysis objects were reduced to significant switching genes, so a gene
#' that switches in one dataset but not the other is simply absent from the other -- you
#' cannot plot its isoform usage there, and you cannot tell "tested and not switching"
#' apart from "not present". The unfiltered objects carry every tested isoform, so this
#' pulls the columns needed for plotting and classification without disturbing the
#' switching calls, which continue to come from the primary (reduced) object.
#'
#' @return list(features = per-isoform context, rep_if = per-replicate IF or NULL)
extract_context_tables <- function(isa_obj) {
  f <- extract_isoform_features(isa_obj)
  f <- tibble::as_tibble(f, .name_repair = "unique")
  want <- c(
    "isoform_id", "gene_id", "gene_name", "condition_1", "condition_2",
    "IF1", "IF2", "dIF", "IF_overall",
    "iso_value_1", "iso_value_2", "iso_overall_mean",
    "gene_value_1", "gene_value_2",
    "isoform_switch_q_value", "gene_switch_q_value", "iso_q_value", "gene_q_value",
    # ORF / NMD annotation. `analyzeORF` runs unconditionally upstream, so PTC is present
    # on the unfiltered table too -- it was simply never extracted. Without it the context
    # cannot answer any consequence question about NON-switching isoforms, which is the
    # whole of the compositional-rescue (Class B) test: those genes need not switch, so the
    # reduced objects cover only 5-21% of them.
    "PTC", "IR"
  )
  keep <- intersect(want, names(f))
  if (!all(c("isoform_id", "gene_id") %in% keep)) {
    stop("extract_context_tables: unfiltered object lacks isoform_id / gene_id.")
  }
  features <- f[, keep, drop = FALSE]

  # The orfAnalysis slot carries the quantitative ORF description behind the PTC flag:
  # ORF length, and the distance from the stop codon to the last exon-exon junction, which
  # is the actual NMD rule (>50-55 nt upstream of the last junction triggers decay). Keeping
  # it lets downstream code re-derive NMD status with an explicit rule rather than trusting
  # a single boolean, and distinguishes "no ORF found" from "ORF found, not a PTC".
  orf <- if (is.list(isa_obj) && !is.null(isa_obj[["orfAnalysis"]])) {
    tibble::as_tibble(isa_obj[["orfAnalysis"]], .name_repair = "unique")
  } else {
    NULL
  }
  if (!is.null(orf) && "isoform_id" %in% names(orf)) {
    orf_want <- intersect(
      c("isoform_id", "orfTransciptStart", "orfTransciptEnd", "orfTransciptLength",
        "stopDistanceToLastJunction", "stopIndex", "orf_origin"),
      names(orf)
    )
    orf <- orf[, orf_want, drop = FALSE]
    orf <- orf[!duplicated(orf$isoform_id), , drop = FALSE]
    features <- dplyr::left_join(features, orf, by = "isoform_id")
  }

  grab_rep <- function(slot) {
    if (!is.list(isa_obj) || is.null(isa_obj[[slot]])) {
      return(NULL)
    }
    r <- tibble::as_tibble(isa_obj[[slot]], .name_repair = "unique")
    if ("isoform_id" %in% names(r)) r else NULL
  }
  # Per-replicate IF and expression. Expression is needed to relate isoform-fraction
  # instability to abundance -- IF is a ratio, so at low expression it swings wildly and a
  # large dIF can be pure noise. Without the replicate expression there is no principled
  # way to set an abundance floor.
  list(
    features = features,
    rep_if = grab_rep("isoformRepIF"),
    rep_expr = grab_rep("isoformRepExpression")
  )
}

#' Load a context table written by 01, or NULL when the dataset has no unfiltered object
load_context_table <- function(processed_dir, label, what = c("features", "rep_if", "rep_expr")) {
  what <- match.arg(what)
  stem <- switch(what,
    features = "isoformContext_",
    rep_if = "isoformContextRepIF_",
    rep_expr = "isoformContextRepExpr_"
  )
  p <- file.path(processed_dir, paste0(stem, sanitize(label), ".rds"))
  if (!file.exists(p)) NULL else readRDS(p)
}

# ---- Provenance -----------------------------------------------------------------------

#' TRUE when every gene in the table already carries a significant gene-level q.
#'
#' IsoformSwitchAnalyzeR's `isoformSwitchTestDEXSeq(reduceToSwitchingGenes = TRUE)` (the
#' default) drops non-switching genes before the object is saved. Downstream that makes
#' "genes present" mean "genes already called significant", so overlap counts, novel
#' fractions and any background-dependent test describe the subsetting rather than the
#' biology. Detect it explicitly instead of assuming a full tested background.
detect_isa_reduction <- function(iso, cfg) {
  z <- tibble::as_tibble(iso, .name_repair = "unique")
  g_q <- as.numeric((cfg$significance %||% list())$gene_q %||% 0.05)
  n_genes <- length(unique(z$gene_id))
  if (!"gene_switch_q_value" %in% names(z) || !n_genes) {
    return(list(
      n_genes = n_genes, n_genes_significant = NA_integer_,
      pct_genes_significant = NA_real_, looks_reduced = NA
    ))
  }
  sig_by_gene <- tapply(
    as.numeric(z$gene_switch_q_value), z$gene_id,
    function(v) any(is.finite(v) & v < g_q)
  )
  n_sig <- sum(sig_by_gene, na.rm = TRUE)
  list(
    n_genes = n_genes,
    n_genes_significant = as.integer(n_sig),
    pct_genes_significant = 100 * n_sig / n_genes,
    looks_reduced = isTRUE(n_sig == n_genes)
  )
}

git_sha <- function(root = find_project_root()) {
  out <- tryCatch(
    system2("git", c("-C", shQuote(root), "rev-parse", "HEAD"), stdout = TRUE, stderr = FALSE),
    error = function(e) NA_character_, warning = function(w) NA_character_
  )
  if (!length(out) || !nzchar(out[[1L]])) NA_character_ else out[[1L]]
}

#' Record what produced the current contents of results/.
#'
#' The pipeline advertises reproducibility but previously recorded nothing about the
#' code version, package versions or thresholds behind a given set of tables.
write_run_manifest <- function(script, cfg, root, extra = list()) {
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    message("jsonlite not installed; skipping run manifest.")
    return(invisible(NULL))
  }
  dir <- ensure_dir(resolve_path(
    (cfg$paths %||% list())$results_tables %||% "results/tables",
    root = root
  ))
  path <- file.path(dirname(dir), "run_manifest.json")
  prev <- if (file.exists(path)) {
    tryCatch(jsonlite::read_json(path, simplifyVector = TRUE), error = function(e) list())
  } else {
    list()
  }
  si <- utils::sessionInfo()
  sha <- git_sha(root)
  key_pkgs <- vapply(
    c("dplyr", "tibble", "ggplot2", "yaml", "IsoformSwitchAnalyzeR"),
    function(p) {
      v <- tryCatch(as.character(utils::packageVersion(p)), error = function(e) NA_character_)
      if (is.na(v)) "not installed" else v
    },
    character(1)
  )
  # A manifest entry that silently records `git_sha: null` or a missing package looks like
  # provenance but is not: an Aug-2026 rerun under a transient R 4.5.3 wrote both and nobody
  # noticed until the tables were already committed. Say so at run time instead.
  if (is.na(sha)) {
    warning(
      "write_run_manifest(): git SHA unavailable, recording null. Is `git` on PATH ",
      "and is ", root, " inside the work tree? This entry cannot be traced to a commit.",
      call. = FALSE, immediate. = TRUE
    )
  }
  if (any(key_pkgs == "not installed")) {
    warning(
      "write_run_manifest(): key package(s) not installed: ",
      paste(names(key_pkgs)[key_pkgs == "not installed"], collapse = ", "),
      ". Results may not be comparable to earlier runs.",
      call. = FALSE, immediate. = TRUE
    )
  }
  entry <- list(
    script = script,
    run_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
    git_sha = sha,
    r_version = si$R.version$version.string,
    platform = si$platform,
    key_packages = key_pkgs,
    significance = cfg$significance,
    analysis = cfg$analysis
  )
  for (nm in names(extra)) entry[[nm]] <- extra[[nm]]
  runs <- prev$runs %||% list()
  if (is.data.frame(runs)) runs <- split(runs, seq_len(nrow(runs)))
  runs <- c(Filter(function(r) !identical(r$script, script), runs), list(entry))
  jsonlite::write_json(
    list(runs = unname(runs)), path,
    auto_unbox = TRUE, pretty = TRUE, null = "null"
  )
  message("Run manifest: ", path)
  invisible(path)
}

# ---- Isoform structure: exons, pairwise splicing events, reading frame ------------------
# Used by 15_isoform_characterization.R. Everything here is PAIRWISE and directional -- an
# event describes how isoform A differs FROM isoform B -- which is deliberately not the same
# question IsoformSwitchAnalyzeR's analyzeAlternativeSplicing answers (each isoform against
# the gene's hypothetical pre-RNA). See the header of 15 for why the bench question is the
# pairwise one.

#' @description
#' Exon records for a set of gene ids, pulled out of a GTF/GFF3 with grep before parsing.
#'
#' The merged annotations here are ~600 MB with ~2.1M exon records, so reading the whole file
#' to keep 40 genes costs minutes and gigabytes. Matching `gene_id=<id>;` first turns that
#' into seconds. 13_assayability_screen.R carries an inline copy of this from before it was
#' factored out; do not add a third.
#'
#' @return tibble(gene_id, isoform_id, chr, start, end, strand), sorted by isoform and start
parse_gff3_exons <- function(gff_path, gene_ids) {
  if (!file.exists(gff_path)) {
    stop("GFF3 not found: ", gff_path)
  }
  gene_ids <- unique(as.character(gene_ids))
  if (!length(gene_ids)) {
    stop("parse_gff3_exons: no gene ids requested.")
  }
  pat <- paste0("gene_id=(", paste(gene_ids, collapse = "|"), ");")
  tmp <- tempfile(fileext = ".tsv")
  on.exit(unlink(tmp), add = TRUE)
  system2("grep", c("-E", shQuote(pat), shQuote(gff_path)), stdout = tmp, stderr = FALSE)
  if (!file.exists(tmp) || file.info(tmp)$size == 0) {
    stop("No GFF3 records matched the requested gene ids.")
  }
  g <- utils::read.delim(tmp, header = FALSE, comment.char = "#", quote = "",
                         stringsAsFactors = FALSE)
  names(g)[c(1, 3, 4, 5, 7, 9)] <- c("chr", "feature", "start", "end", "strand", "attr")
  tibble::as_tibble(g) |>
    dplyr::filter(.data$feature == "exon") |>
    dplyr::mutate(
      isoform_id = sub("^.*Parent=([^;]+).*$", "\\1", .data$attr),
      gene_id = sub("^.*gene_id=([^;]+).*$", "\\1", .data$attr)
    ) |>
    dplyr::filter(.data$gene_id %in% gene_ids) |>
    dplyr::select("gene_id", "isoform_id", "chr", "start", "end", "strand") |>
    dplyr::arrange(.data$isoform_id, .data$start)
}

#' Introns implied by a sorted 2-column exon matrix (start, end). Zero rows for single exons.
exons_to_introns <- function(e) {
  if (!is.matrix(e)) e <- as.matrix(e)
  if (nrow(e) < 2L) return(matrix(numeric(0), ncol = 2L, dimnames = list(NULL, c("start", "end"))))
  o <- order(e[, 1L])
  e <- e[o, , drop = FALSE]
  cbind(start = e[-nrow(e), 2L] + 1L, end = e[-1L, 1L] - 1L)
}

#' @description
#' Classify how isoform A differs from isoform B, both given as sorted (start, end) exon
#' matrices on the same strand.
#'
#' EVERY INTERNAL EVENT IS RESTRICTED TO THE SPAN THE TWO TRANSCRIPTS SHARE. Without that,
#' an isoform that simply starts 40 kb downstream is scored as "skipping" every upstream exon
#' of its partner, and the skipping count becomes a restatement of the TSS difference. Only
#' the overlapping interval can distinguish a real internal event from a boundary effect;
#' differences outside it are reported as alt_TSS / alt_TTS and nothing else.
#'
#' @return list of logical event flags, a printable `summary`, bp unique to each side, and
#'   `diff_introns` -- the introns present in exactly one of the two, for splice-site lookup
pairwise_as_events <- function(ea, eb, strand, end_tolerance = 50L) {
  ea <- as.matrix(ea); eb <- as.matrix(eb)
  ea <- ea[order(ea[, 1L]), , drop = FALSE]
  eb <- eb[order(eb[, 1L]), , drop = FALSE]
  plus <- !identical(as.character(strand)[1L], "-")

  ia <- exons_to_introns(ea)
  ib <- exons_to_introns(eb)

  # Transcript ends, strand-aware: the 5' end is the low coordinate on +, the high one on -.
  a_five <- if (plus) min(ea[, 1L]) else max(ea[, 2L])
  b_five <- if (plus) min(eb[, 1L]) else max(eb[, 2L])
  a_three <- if (plus) max(ea[, 2L]) else min(ea[, 1L])
  b_three <- if (plus) max(eb[, 2L]) else min(eb[, 1L])

  ov_start <- max(min(ea[, 1L]), min(eb[, 1L]))
  ov_end <- min(max(ea[, 2L]), max(eb[, 2L]))
  inside <- function(m) {
    if (!nrow(m)) return(m)
    m[m[, 1L] >= ov_start & m[, 2L] <= ov_end, , drop = FALSE]
  }
  ia_in <- inside(ia); ib_in <- inside(ib)
  ea_in <- inside(ea); eb_in <- inside(eb)

  contained_in <- function(x, box) {
    # rows of x fully inside any row of box
    if (!nrow(x) || !nrow(box)) return(logical(nrow(x)))
    vapply(seq_len(nrow(x)), function(i) {
      any(box[, 1L] <= x[i, 1L] & box[, 2L] >= x[i, 2L])
    }, NA)
  }
  # A retains an intron that B splices out: B's intron sits inside one of A's exons.
  ir <- any(contained_in(ib_in, ea))
  # A skips an exon B has: B's exon sits inside one of A's introns.
  es <- any(contained_in(eb_in, ia))
  # ...and the mirror image, needed for mutually-exclusive detection.
  es_rev <- any(contained_in(ea_in, ib))
  mxe <- es && es_rev

  # Alternative donor / acceptor: two introns that overlap and agree on exactly one end.
  # Which end is the donor depends on strand -- on the minus strand the intron's high
  # coordinate is the 5' splice site.
  alt_d <- alt_a <- FALSE
  if (nrow(ia_in) && nrow(ib_in)) {
    for (p in seq_len(nrow(ia_in))) {
      for (q in seq_len(nrow(ib_in))) {
        if (ia_in[p, 2L] < ib_in[q, 1L] || ia_in[p, 1L] > ib_in[q, 2L]) next
        same_lo <- ia_in[p, 1L] == ib_in[q, 1L]
        same_hi <- ia_in[p, 2L] == ib_in[q, 2L]
        if (same_lo && !same_hi) {
          if (plus) alt_a <- TRUE else alt_d <- TRUE
        } else if (same_hi && !same_lo) {
          if (plus) alt_d <- TRUE else alt_a <- TRUE
        }
      }
    }
  }

  key <- function(m) if (nrow(m)) paste(m[, 1L], m[, 2L], sep = "_") else character(0)
  ka <- key(ia); kb <- key(ib)
  only_a <- ia[!(ka %in% kb), , drop = FALSE]
  only_b <- ib[!(kb %in% ka), , drop = FALSE]
  diff_introns <- rbind(
    if (nrow(only_a)) cbind(only_a, in_a = 1L) else NULL,
    if (nrow(only_b)) cbind(only_b, in_a = 0L) else NULL
  )
  if (is.null(diff_introns)) {
    diff_introns <- matrix(numeric(0), ncol = 3L,
                           dimnames = list(NULL, c("start", "end", "in_a")))
  }

  exonic_bp <- function(m) {
    if (!nrow(m)) return(0L)
    n <- norm_intervals(m)
    sum(n[, 2L] - n[, 1L] + 1L)
  }
  bp_unique_a <- exonic_bp(interval_subtract(norm_intervals(ea), norm_intervals(eb)))
  bp_unique_b <- exonic_bp(interval_subtract(norm_intervals(eb), norm_intervals(ea)))

  # TSS/TTS DIFFERENCES NEED A TOLERANCE, or the flag means nothing. Assembled transcript ends
  # are almost never identical to the base, so "start coordinate differs" fires on essentially
  # every pair and stops discriminating a different PROMOTER from a ragged end. Three
  # quantities are reported instead of one flag:
  #   *_distance_bp    how far apart the two ends are
  #   alt_TSS/alt_TTS  that distance exceeds `end_tolerance` (default 50 bp)
  #   alt_first_exon   the two terminal exons do not overlap AT ALL -- the strong form, and
  #                    the only one that implies a distinct promoter or polyA site
  tss_dist <- abs(a_five - b_five)
  tts_dist <- abs(a_three - b_three)
  term_exon <- function(e, five) {
    if (plus == five) e[1L, , drop = FALSE] else e[nrow(e), , drop = FALSE]
  }
  no_overlap <- function(x, y) x[1L, 2L] < y[1L, 1L] || x[1L, 1L] > y[1L, 2L]
  alt_first <- no_overlap(term_exon(ea, TRUE), term_exon(eb, TRUE))
  alt_last <- no_overlap(term_exon(ea, FALSE), term_exon(eb, FALSE))

  flags <- c(
    alt_TSS = tss_dist > end_tolerance, alt_TTS = tts_dist > end_tolerance,
    exon_skipping = es, intron_retention = ir,
    alt_donor = alt_d, alt_acceptor = alt_a, mutually_exclusive = mxe
  )
  list(
    alt_TSS = unname(flags[["alt_TSS"]]), alt_TTS = unname(flags[["alt_TTS"]]),
    tss_distance_bp = tss_dist, tts_distance_bp = tts_dist,
    alt_first_exon = alt_first, alt_last_exon = alt_last,
    exon_skipping = es, intron_retention = ir,
    alt_donor = alt_d, alt_acceptor = alt_a, mutually_exclusive = mxe,
    n_events = sum(flags),
    summary = if (any(flags)) paste(names(flags)[flags], collapse = "+") else "no difference",
    bp_unique_a = bp_unique_a, bp_unique_b = bp_unique_b,
    diff_introns = diff_introns
  )
}

#' Merge overlapping/abutting closed intervals in a 2-column matrix.
norm_intervals <- function(m) {
  if (!is.matrix(m)) m <- as.matrix(m)
  if (!nrow(m)) return(matrix(numeric(0), ncol = 2L))
  m <- m[order(m[, 1L]), , drop = FALSE]
  s <- m[1L, 1L]; e <- m[1L, 2L]; out <- list()
  for (i in seq_len(nrow(m))[-1L]) {
    if (m[i, 1L] <= e + 1L) {
      e <- max(e, m[i, 2L])
    } else {
      out[[length(out) + 1L]] <- c(s, e); s <- m[i, 1L]; e <- m[i, 2L]
    }
  }
  out[[length(out) + 1L]] <- c(s, e)
  do.call(rbind, out)
}

#' a minus b, for 2-column closed-interval matrices.
interval_subtract <- function(a, b) {
  if (!nrow(a) || !nrow(b)) return(a)
  out <- list()
  for (i in seq_len(nrow(a))) {
    seg <- list(c(a[i, 1L], a[i, 2L]))
    for (j in seq_len(nrow(b))) {
      nxt <- list()
      for (s in seg) {
        if (b[j, 2L] < s[1L] || b[j, 1L] > s[2L]) { nxt[[length(nxt) + 1L]] <- s; next }
        if (b[j, 1L] > s[1L]) nxt[[length(nxt) + 1L]] <- c(s[1L], b[j, 1L] - 1L)
        if (b[j, 2L] < s[2L]) nxt[[length(nxt) + 1L]] <- c(b[j, 2L] + 1L, s[2L])
      }
      seg <- nxt
      if (!length(seg)) break
    }
    if (length(seg)) out <- c(out, seg)
  }
  if (!length(out)) return(matrix(numeric(0), ncol = 2L))
  do.call(rbind, out)
}

#' @description
#' Coding blocks of one isoform in GENOMIC coordinates, each tagged with how many coding
#' bases precede it. The tag is what makes a reading-frame comparison possible: two isoforms
#' share a genomic position but not necessarily the same codon boundary at it.
#'
#' @return matrix(gstart, gend, off) or NULL when no ORF was called
orf_genomic_blocks <- function(e, strand, orf_start, orf_end) {
  if (is.na(orf_start) || is.na(orf_end) || orf_end < orf_start) return(NULL)
  e <- as.matrix(e)
  plus <- !identical(as.character(strand)[1L], "-")
  # Exons in TRANSCRIPT order: 5' to 3', which is descending genomic start on the minus strand.
  e <- e[order(e[, 1L], decreasing = !plus), , drop = FALSE]
  out <- list()
  tpos <- 1L
  for (i in seq_len(nrow(e))) {
    len <- e[i, 2L] - e[i, 1L] + 1L
    t1 <- tpos; t2 <- tpos + len - 1L
    tpos <- t2 + 1L
    s <- max(t1, orf_start); en <- min(t2, orf_end)
    if (s > en) next
    if (plus) {
      gstart <- e[i, 1L] + (s - t1); gend <- e[i, 1L] + (en - t1)
    } else {
      gend <- e[i, 2L] - (s - t1); gstart <- e[i, 2L] - (en - t1)
    }
    # unname(): indexing a matrix that carries dimnames leaves the column name attached to
    # the scalar, and c(gstart = <named>) then builds "gstart.start" instead of "gstart".
    out[[length(out) + 1L]] <- c(gstart = unname(gstart), gend = unname(gend),
                                 off = unname(s - orf_start))
  }
  if (!length(out)) return(NULL)
  do.call(rbind, out)
}

#' @description
#' Do two isoforms read their SHARED coding sequence in the same frame, and where does their
#' coding sequence differ?
#'
#' This is the step that turns a structural difference into a protein claim. Two isoforms can
#' share most of their exons and still encode unrelated C-termini if an indel upstream shifted
#' the frame; conversely a large structural difference confined to UTR changes no protein at
#' all. Neither case is visible from exon overlap.
#'
#' @return list(shared_coding_bp, relationship, coding_change)
frame_relationship <- function(ea, eb, strand, oa, ob) {
  unknown <- list(shared_coding_bp = NA_integer_,
                  relationship = "unknown (ORF not called)",
                  coding_change = "unknown (ORF not called)",
                  coding_change_detail = "unknown (ORF not called)",
                  coding_bp_unique_candidate = NA_integer_,
                  coding_bp_unique_partner = NA_integer_)
  if (is.null(oa) || is.null(ob)) return(unknown)
  ba <- orf_genomic_blocks(ea, strand, oa$start, oa$end)
  bb <- orf_genomic_blocks(eb, strand, ob$start, ob$end)
  if (is.null(ba) || is.null(bb)) return(unknown)
  plus <- !identical(as.character(strand)[1L], "-")

  shared <- 0L; shifted <- FALSE; agreed <- FALSE
  for (p in seq_len(nrow(ba))) {
    for (q in seq_len(nrow(bb))) {
      x1 <- max(ba[p, "gstart"], bb[q, "gstart"])
      x2 <- min(ba[p, "gend"], bb[q, "gend"])
      if (x1 > x2) next
      shared <- shared + (x2 - x1 + 1L)
      # Codon offset at the 5'-most base of the shared block, in each isoform's own ORF.
      oa_off <- if (plus) ba[p, "off"] + (x1 - ba[p, "gstart"]) else
        ba[p, "off"] + (ba[p, "gend"] - x2)
      ob_off <- if (plus) bb[q, "off"] + (x1 - bb[q, "gstart"]) else
        bb[q, "off"] + (bb[q, "gend"] - x2)
      if ((oa_off - ob_off) %% 3 == 0) agreed <- TRUE else shifted <- TRUE
    }
  }
  if (shared == 0L) {
    return(list(shared_coding_bp = 0L, relationship = "no shared coding sequence",
                coding_change = "no shared coding sequence",
                coding_change_detail = "no shared coding sequence",
                coding_bp_unique_candidate = NA_integer_,
                coding_bp_unique_partner = NA_integer_))
  }
  rel <- if (shifted && agreed) {
    "partly frameshifted"
  } else if (shifted) {
    "frameshifted"
  } else {
    "same frame throughout"
  }

  cds_a <- norm_intervals(ba[, c("gstart", "gend"), drop = FALSE])
  cds_b <- norm_intervals(bb[, c("gstart", "gend"), drop = FALSE])
  only_a <- interval_subtract(cds_a, cds_b)
  only_b <- interval_subtract(cds_b, cds_a)
  bp <- function(m) if (nrow(m)) sum(m[, 2L] - m[, 1L] + 1L) else 0L

  # `categ` and `change` are built together and returned as a pair: the category is what you
  # tally, the detail is what you read. An earlier version set the category with `<<-` inside
  # the if/else, which assigns in the ENCLOSING environment -- the local stayed NA and every
  # tallied category came back missing.
  res <- if (shifted) {
    list(categ = "frameshift",
         change = paste0("frameshift -- the shared exons are read in different frames (",
                         shared, " bp shared)"))
  } else if (bp(only_a) == 0L && bp(only_b) == 0L) {
    list(categ = "UTR only",
         change = "coding sequence identical -- the difference is UTR only")
  } else {
    # WHERE the difference sits, judged against the span the two ORFs share rather than by
    # counting blocks. A pure N-terminal truncation spanning ten exons produces ten unique
    # blocks, so "more blocks than boundaries" -- the obvious test -- calls it internal every
    # time. Position relative to the shared coding span is what actually distinguishes them.
    shared_iv <- interval_subtract(cds_a, only_a)
    if (!nrow(shared_iv)) shared_iv <- cds_a
    sh_lo <- min(shared_iv[, 1L]); sh_hi <- max(shared_iv[, 2L])
    where <- function(m) {
      if (!nrow(m)) return(character(0))
      lab <- vapply(seq_len(nrow(m)), function(i) {
        if (m[i, 2L] < sh_lo) if (plus) "N-terminal" else "C-terminal"
        else if (m[i, 1L] > sh_hi) if (plus) "C-terminal" else "N-terminal"
        else "internal"
      }, NA_character_)
      unique(lab)
    }
    parts <- unique(c(where(only_a), where(only_b)))
    parts <- parts[order(match(parts, c("N-terminal", "internal", "C-terminal")))]
    if (!length(parts)) parts <- "internal"
    list(categ = paste(parts, collapse = " + "),
         change = paste0(paste(parts, collapse = " + "), " difference (",
                         bp(only_a), " bp coding unique to candidate, ",
                         bp(only_b), " bp unique to partner)"))
  }
  list(shared_coding_bp = as.integer(shared), relationship = rel,
       coding_change = res$categ, coding_change_detail = res$change,
       coding_bp_unique_candidate = bp(only_a), coding_bp_unique_partner = bp(only_b))
}

#' @description
#' Splice-site sequence at one intron, in transcript orientation.
#'
#' Windows follow the usual maximum-entropy convention: 3 exonic + 6 intronic bases at the
#' donor, 20 intronic + 3 exonic at the acceptor. The dinucleotides are reported separately
#' because they carry the spliceosome assignment: GT-AG and GC-AG are U2 introns, AT-AC is
#' U12 (minor spliceosome), and anything else is either a genuine non-canonical site or --
#' far more often -- an assembly artefact worth distrusting.
splice_site_seqs <- function(genome, chrom, istart, iend, strand) {
  plus <- !identical(as.character(strand)[1L], "-")
  grab <- function(s, e) {
    if (s < 1) return(NA_character_)
    sq <- Biostrings::getSeq(genome, names = chrom, start = s, end = e)
    if (!plus) sq <- Biostrings::reverseComplement(sq)
    as.character(sq)
  }
  if (plus) {
    donor <- grab(istart - 3L, istart + 5L)
    acceptor <- grab(iend - 19L, iend + 3L)
  } else {
    donor <- grab(iend - 5L, iend + 3L)
    acceptor <- grab(istart - 3L, istart + 19L)
  }
  # Window arithmetic, spelled out because an off-by-one here silently reports every intron
  # as non-canonical. Donor window is [istart-3, istart+5], 9 nt: positions 1-3 are exonic, so
  # the intron's first base is position 4 and the donor dinucleotide is 4-5.
  # Acceptor window is [iend-19, iend+3], 23 nt: positions 1-20 are intronic and 21-23 exonic,
  # so the intron's LAST base is position 20 and the acceptor dinucleotide is 19-20.
  d2 <- if (is.na(donor)) NA_character_ else substr(donor, 4L, 5L)
  a2 <- if (is.na(acceptor)) NA_character_ else substr(acceptor, 19L, 20L)
  cls <- if (is.na(d2) || is.na(a2)) {
    NA_character_
  } else if (d2 == "GT" && a2 == "AG") {
    "GT-AG (U2, canonical)"
  } else if (d2 == "GC" && a2 == "AG") {
    "GC-AG (U2, minor variant)"
  } else if (d2 == "AT" && a2 == "AC") {
    "AT-AC (U12, minor spliceosome)"
  } else {
    paste0("non-canonical (", d2, "-", a2, ")")
  }
  list(donor = donor, acceptor = acceptor, donor_di = d2, acceptor_di = a2, intron_class = cls)
}

#' Translate an ORF given transcript coordinates. Returns a character scalar, or NULL.
translate_orf <- function(dna, orf_start, orf_end) {
  if (is.na(orf_start) || is.na(orf_end)) return(NULL)
  n <- length(dna)
  if (orf_start < 1L || orf_end > n || orf_end <= orf_start) return(NULL)
  sub <- Biostrings::subseq(dna, start = orf_start, end = orf_end)
  # Trim to a whole number of codons; ISA's ORF end includes the stop codon.
  keep <- (length(sub) %/% 3L) * 3L
  if (keep < 3L) return(NULL)
  sub <- Biostrings::subseq(sub, start = 1L, end = keep)
  aa <- suppressWarnings(
    as.character(Biostrings::translate(sub, if.fuzzy.codon = "solve"))
  )
  sub("\\*$", "", aa)
}
