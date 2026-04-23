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

# ---- Gene-level summary from isoform tibble ------------------------------------------------

#' @description
#' One row per gene. A *switching* isoform passes isoform q, gene q (if present), and
#' optional |dIF| cutoffs. *novel_involved* is TRUE if any switching isoform is
#' flagged PacBio-novel.
summarize_genes_from_isoform_table <- function(iso, cfg) {
  if (!tibble::is_tibble(iso) && !is.data.frame(iso)) {
    stop("summarize_genes_from_isoform_table: `iso` must be a tibble or data frame.")
  }
  iso <- tibble::as_tibble(iso, .name_repair = "unique")
  if (!"gene_id" %in% names(iso)) {
    stop("Column 'gene_id' not found in isoform table.")
  }
  sig <- cfg$significance
  if (is.null(sig)) {
    stop("config must define `significance` (e.g. isoform_q, gene_q).")
  }
  i_q <- as.numeric(sig$isoform_q %||% 0.05)
  g_q <- as.numeric(sig$gene_q %||% 0.05)
  m_dif <- as.numeric(sig$min_abs_dif %||% 0.0)
  if (!"isoform_switch_q_value" %in% names(iso) && !"gene_switch_q_value" %in% names(iso)) {
    stop("Need at least one of: isoform_switch_q_value, gene_switch_q_value.")
  }
  if (!"gene_name" %in% names(iso)) {
    iso$gene_name <- NA_character_
  }
  if (!"dIF" %in% names(iso)) {
    iso$dIF <- NA_real_
  }
  if (!"is_novel_pacbio" %in% names(iso)) {
    iso$is_novel_pacbio <- NA
  }
  iso$q_i <- if ("isoform_switch_q_value" %in% names(iso)) {
    as.numeric(iso$isoform_switch_q_value)
  } else {
    as.numeric(iso$gene_switch_q_value)
  }
  iso$q_g <- if ("gene_switch_q_value" %in% names(iso)) {
    as.numeric(iso$gene_switch_q_value)
  } else {
    as.numeric(NA)
  }
  iso$dIF_n <- as.numeric(iso$dIF)
  ok_i <- is.finite(iso$q_i) & (iso$q_i < i_q)
  ok_g <- is.na(iso$q_g) | (is.finite(iso$q_g) & (iso$q_g < g_q))
  ok_d <- is.na(iso$dIF_n) | is.nan(iso$dIF_n) | (abs(iso$dIF_n) >= m_dif)
  iso$is_switching <- ok_i & ok_g & ok_d
  iso$is_switching <- replace(iso$is_switching, is.na(iso$is_switching), FALSE)
  min_finite <- function(v) {
    if (!length(v)) {
      return(NA_real_)
    }
    m <- suppressWarnings(min(v, na.rm = TRUE))
    if (is.infinite(m)) {
      as.numeric(NA)
    } else {
      m
    }
  }
  max_abs_finite <- function(v) {
    w <- abs(v)
    w <- w[is.finite(w)]
    if (!length(w)) {
      as.numeric(NA)
    } else {
      max(w, na.rm = TRUE)
    }
  }
  iso |>
    dplyr::group_by(.data$gene_id, .data$gene_name) |>
    dplyr::summarize(
      n_isoforms = dplyr::n(),
      n_switching_isoforms = sum(.data$is_switching, na.rm = TRUE),
      min_isoform_switch_q = min_finite(.data$q_i),
      min_gene_switch_q = min_finite(.data$q_g),
      max_abs_dif = max_abs_finite(.data$dIF_n),
      novel_involved = any(
        .data$is_switching & (.data$is_novel_pacbio %in% TRUE),
        na.rm = TRUE
      ),
      .groups = "drop"
    )
}
