#!/usr/bin/env Rscript
# 22: Render the HT-only manuscript figures.
#
# Run (after 21): Rscript scripts/22_manuscript_figures.R
#
# This is a THIN entry point. The plotting code lives in explore/manuscript_fig[1-4]_explore.R
# so that there is exactly one copy of it: those scripts are standalone (they read only
# results/tables/manuscript/*.csv and print rather than write), and setting
# `manuscript.save_dir` is what turns them into figure producers. Editing a panel means
# editing the explore script, and the version a colleague runs by hand is the same code that
# made the figure in the paper.
#
# Scope: HT reference, H and T libraries. No UBL5 knockout data. Script 21 asserts that.

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
config_rel <- if (length(args_r) > 0L) args_r[[1L]] else "config/config.yml"
cfg <- load_yaml_config(config_rel, root = root)
paths <- cfg$paths %||% list()
results_tables <- resolve_path(paths$results_tables %||% "results/tables", root = root)
results_figures <- resolve_path(paths$results_figures %||% "results/figures", root = root)

man_tbl <- file.path(results_tables, "manuscript")
if (!dir.exists(man_tbl) ||
    !file.exists(file.path(man_tbl, "manuscript_panel_index.csv"))) {
  stop("No panel tables at ", man_tbl, " -- run: Rscript scripts/21_manuscript_data.R")
}
out_dir <- ensure_dir(file.path(results_figures, "manuscript"))

idx <- utils::read.csv(file.path(man_tbl, "manuscript_panel_index.csv"),
                       stringsAsFactors = FALSE)
message("Panels available: ", nrow(idx),
        " (", sum(idx$status == "real"), " real, ",
        sum(idx$status == "placeholder"), " placeholder)")

# The explore scripts resolve their inputs relative to the project root.
owd <- setwd(root)
on.exit(setwd(owd), add = TRUE)
old_opt <- options(manuscript.save_dir = file.path(results_figures, "manuscript"))
on.exit(options(old_opt), add = TRUE)

figs <- sprintf("explore/manuscript_fig%d_explore.R", 1:4)
missing <- figs[!file.exists(figs)]
if (length(missing)) stop("Missing figure script(s): ", paste(missing, collapse = ", "))

for (f in figs) {
  message("\n---- ", f, " ----")
  # local = new.env() so each figure's panel objects cannot leak into the next one
  source(f, local = new.env(), chdir = FALSE)
}

written <- list.files(out_dir, pattern = "\\.png$", full.names = FALSE)
message("\n---- figures in ", out_dir, " ----")
for (w in sort(written)) message("  ", w)

write_run_manifest(script = "22_manuscript_figures.R", cfg = cfg, root = root)
message("\n22_manuscript_figures.R: done")
