# Sourced at the start of each pipeline script. Adds helpers to the global R session.
# Usage from project root:  source("utils/bootstrap.R")
#            from scripts/:   source(file.path("..", "utils", "bootstrap.R"))
# Or: set environment variable ISOFORM_PROJECT_ROOT to the project path.

root <- (function() {
  if (nzchar(Sys.getenv("ISOFORM_PROJECT_ROOT", ""))) {
    return(normalizePath(Sys.getenv("ISOFORM_PROJECT_ROOT"), mustWork = TRUE))
  }
  p <- normalizePath(getwd(), mustWork = FALSE)
  for (i in 1:20) {
    if (file.exists(file.path(p, "config", "config.yml"))) {
      return(p)
    }
    p2 <- dirname(p)
    if (identical(p, p2)) {
      break
    }
    p <- p2
  }
  NULL
})()
if (is.null(root)) {
  stop("Could not find project root (config/config.yml). `cd` into the project, or set ISOFORM_PROJECT_ROOT")
}
h <- file.path(root, "utils", "helper_functions.R")
if (!file.exists(h)) {
  stop("Not found: ", h)
}
source(h, chdir = FALSE, local = FALSE)
assign("PROJECT_ROOT", root, envir = .GlobalEnv)
invisible(root)
