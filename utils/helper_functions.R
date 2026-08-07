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
    "isoform_switch_q_value", "gene_switch_q_value", "iso_q_value", "gene_q_value"
  )
  keep <- intersect(want, names(f))
  if (!all(c("isoform_id", "gene_id") %in% keep)) {
    stop("extract_context_tables: unfiltered object lacks isoform_id / gene_id.")
  }
  features <- f[, keep, drop = FALSE]

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
  entry <- list(
    script = script,
    run_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
    git_sha = git_sha(root),
    r_version = si$R.version$version.string,
    platform = si$platform,
    key_packages = vapply(
      c("dplyr", "tibble", "ggplot2", "yaml", "IsoformSwitchAnalyzeR"),
      function(p) {
        v <- tryCatch(as.character(utils::packageVersion(p)), error = function(e) NA_character_)
        if (is.na(v)) "not installed" else v
      },
      character(1)
    ),
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
