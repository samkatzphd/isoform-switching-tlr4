#!/usr/bin/env Rscript
# 99: Verify every configured input file before the pipeline touches it.
#
# Run after any download, unzip, or copy, and before 01:
#   Rscript scripts/99_verify_inputs.R
#   Rscript scripts/99_verify_inputs.R --load    # also try loading each R object (slow)
#
# Why this exists: two ISA exports arrived damaged and neither was obvious from the file
# listing. One was truncated (130 MB of an expected 615 MB); the other was full size but
# corrupt mid-stream, and only failed when R tried to deserialise it 20 minutes into a run.
# Both are caught here in seconds.
#
# Exits non-zero if any configured file is missing or fails its integrity check, so it can
# gate a pipeline run:  Rscript scripts/99_verify_inputs.R && Rscript scripts/01_load_data.R

bt <- c("utils/bootstrap.R", file.path("..", "utils", "bootstrap.R"))
b_file <- if (any(f <- vapply(bt, file.exists, NA))) { bt[which(f)[1L]] } else { NA_character_ }
if (is.na(b_file)) {
  stop("Could not find utils/bootstrap.R. Run from project root or set ISOFORM_PROJECT_ROOT.")
}
source(b_file, local = FALSE, chdir = FALSE)

root <- if (exists("PROJECT_ROOT", inherits = FALSE)) {
  get("PROJECT_ROOT", inherits = FALSE)
} else {
  find_project_root()
}
args_r <- commandArgs(trailingOnly = TRUE)
do_load <- "--load" %in% args_r
config_rel <- {
  a <- setdiff(args_r, "--load")
  if (length(a)) a[[1L]] else "config/config.yml"
}
cfg <- load_yaml_config(config_rel, root = root)

# Compression is detected from magic bytes rather than the extension: .Rdata and .rds are
# gzip by default but may be xz or bzip2, and the wrong test would report a false failure.
detect_compression <- function(path) {
  con <- file(path, "rb")
  on.exit(close(con), add = TRUE)
  m <- readBin(con, "raw", n = 6L)
  if (length(m) >= 2L && m[1] == as.raw(0x1f) && m[2] == as.raw(0x8b)) {
    return("gzip")
  }
  if (length(m) >= 6L && identical(as.integer(m[1:6]), c(0xfd, 0x37, 0x7a, 0x58, 0x5a, 0x00))) {
    return("xz")
  }
  if (length(m) >= 3L && rawToChar(m[1:3]) == "BZh") {
    return("bzip2")
  }
  if (length(m) >= 4L && rawToChar(m[1:4]) %in% c("RDX2", "RDX3")) {
    return("none (uncompressed R data)")
  }
  "unknown"
}

test_integrity <- function(path, kind) {
  tool <- switch(kind, gzip = "gzip", xz = "xz", bzip2 = "bzip2", NULL)
  if (is.null(tool)) {
    return(list(ok = NA, detail = "no compression test applicable"))
  }
  if (!nzchar(Sys.which(tool))) {
    return(list(ok = NA, detail = paste0(tool, " not installed")))
  }
  res <- suppressWarnings(system2(tool, c("-t", shQuote(path)), stdout = TRUE, stderr = TRUE))
  status <- attr(res, "status")
  if (is.null(status) || identical(status, 0L)) {
    list(ok = TRUE, detail = "integrity OK")
  } else {
    list(ok = FALSE, detail = paste(utils::head(res, 2), collapse = "; "))
  }
}

targets <- list()
add_target <- function(dataset, role, path) {
  if (is.null(path)) return(invisible(NULL))
  p <- as.character(path)[1L]
  if (!nzchar(p) || identical(p, "null") || is.na(p)) return(invisible(NULL))
  targets[[length(targets) + 1L]] <<- list(
    dataset = dataset, role = role, path = resolve_path(p, root = root)
  )
}
for (k in names(cfg$datasets %||% list())) {
  ds <- cfg$datasets[[k]]
  add_target(k, "isa_path", ds$isa_path)
  add_target(k, "isa_unfiltered_path", ds$isa_unfiltered_path)
  add_target(k, "annotation_path", ds$annotation_path)
}
if (!length(targets)) stop("No dataset paths found in config.")

rows <- list()
for (t in targets) {
  exists_ok <- file.exists(t$path)
  size <- if (exists_ok) file.info(t$path)$size else NA_real_
  kind <- if (exists_ok) tryCatch(detect_compression(t$path), error = function(e) "unreadable") else NA_character_
  chk <- if (exists_ok) test_integrity(t$path, kind) else list(ok = FALSE, detail = "FILE MISSING")
  load_ok <- NA
  load_detail <- ""
  if (do_load && isTRUE(chk$ok) && t$role != "annotation_path") {
    obj_name <- (cfg$datasets[[t$dataset]] %||% list())$object_name %||% NULL
    load_ok <- tryCatch({
      x <- load_isa_input(t$path, object_name = obj_name)
      f <- extract_isoform_features(x)
      load_detail <- paste0(nrow(f), " isoforms")
      rm(x, f)
      invisible(gc())
      TRUE
    }, error = function(e) {
      load_detail <<- conditionMessage(e)
      FALSE
    })
  }
  status <- if (!exists_ok) {
    "MISSING"
  } else if (isFALSE(chk$ok)) {
    "CORRUPT"
  } else if (isFALSE(load_ok)) {
    "UNREADABLE"
  } else {
    "OK"
  }
  message(sprintf(
    "[%-10s] %-20s %-8s %9s  %s%s",
    t$dataset, t$role, status,
    if (is.na(size)) "-" else format(structure(size, class = "object_size"), units = "auto"),
    chk$detail,
    if (nzchar(load_detail)) paste0(" | ", load_detail) else ""
  ))
  rows[[length(rows) + 1L]] <- data.frame(
    dataset = t$dataset, role = t$role, path = t$path,
    exists = exists_ok, size_bytes = size, compression = kind,
    integrity_ok = chk$ok, integrity_detail = chk$detail,
    load_ok = load_ok, load_detail = load_detail,
    status = status, stringsAsFactors = FALSE
  )
}

out <- do.call(rbind, rows)
tab_dir <- ensure_dir(resolve_path((cfg$paths %||% list())$results_tables %||% "results/tables", root = root))
write_table_pair(out, tab_dir, "input_file_verification", cfg = cfg)

bad <- out[out$status != "OK", , drop = FALSE]
if (nrow(bad)) {
  message("\n", nrow(bad), " of ", nrow(out), " configured input(s) FAILED:")
  for (i in seq_len(nrow(bad))) {
    message("  - ", bad$dataset[i], " / ", bad$role[i], ": ", bad$status[i])
    message("      ", bad$path[i])
    message("      ", bad$integrity_detail[i])
  }
  message(
    "\nRe-copy the affected file(s). See the download/verify checklist in README.md.\n",
    "Do not run the pipeline against these inputs."
  )
  quit(status = 1L)
}
message("\nAll ", nrow(out), " configured inputs present and intact.")
