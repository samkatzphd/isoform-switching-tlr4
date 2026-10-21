# ---------------------------------------------------------------------------------------
# Shared helpers for the standalone manuscript figure scripts in explore/.
#
# Sourced by explore/manuscript_fig[1-4]_explore.R. This is the ONE dependency those scripts
# have beyond dplyr/ggplot2 -- it is a sibling file in the same directory, not the pipeline
# bootstrap. Nothing here reads config, the context objects, or the external drive.
#
# Run every figure script from the PROJECT ROOT.
# ---------------------------------------------------------------------------------------

MANUSCRIPT_TBL <- file.path("results", "tables", "manuscript")

#' Read one panel CSV by id, with a useful error if 21 has not been run.
panel_reader <- function(dir = MANUSCRIPT_TBL) {
  if (!dir.exists(dir)) {
    stop("Missing ", dir, " -- run: Rscript scripts/21_manuscript_data.R\n",
         "  (and run from the project root; current wd is ", getwd(), ")")
  }
  function(id) {
    f <- file.path(dir, paste0("panel_", id, ".csv"))
    if (!file.exists(f)) stop("Missing ", f, " -- run scripts/21_manuscript_data.R")
    utils::read.csv(f, stringsAsFactors = FALSE)
  }
}

# ---- shared scales ----------------------------------------------------------------------

TIER_LEVELS <- c("reference", "pb_redundant", "pb_contained", "pb_novel")
TIER_LABELS <- c("Reference model (ENST)", "PacBio, chain already in Ensembl",
                 "PacBio, contained fragment (c)", "PacBio, novel structure")
TIER_COLS <- c(reference = "#4F7CAC", pb_redundant = "#9DBEDC",
               pb_contained = "#F2B134", pb_novel = "#C1443C")
names(TIER_LABELS) <- TIER_LEVELS

LPS_COLS <- c(minus = "#B8B8B8", plus = "#1F6F54")

# Tables written by different pipeline scripts key the dataset two ways: some by the full
# label (T_HT_anchor) and some by the short config key (T_HT). Map both, so every figure
# shows the same facet strip regardless of which table a panel came from.
DS_LEVELS <- c("THP-1 (T) \u00b1 LPS", "Primary MDM (H) \u00b1 LPS")
DS_PRETTY <- c(T_HT_anchor        = DS_LEVELS[1], T_HT = DS_LEVELS[1],
               H_HT_supplementary = DS_LEVELS[2], H_HT = DS_LEVELS[2])
pretty_ds <- function(x) {
  out <- unname(DS_PRETTY[as.character(x)])
  unmapped <- is.na(out)
  if (any(unmapped)) {
    stop("pretty_ds(): unmapped dataset value(s) -- is UT-space data leaking in? ",
         paste(unique(as.character(x)[unmapped]), collapse = ", "))
  }
  factor(out, levels = DS_LEVELS)
}

#' House theme. Matches the pipeline's theme_bw base so these figures sit alongside the
#' existing results/figures output rather than looking like a different project.
ms_theme <- function(base_size = 11) {
  theme_bw(base_size = base_size) +
    theme(
      plot.subtitle = element_text(face = "bold", size = base_size, hjust = 0),
      plot.title = element_text(face = "bold", size = base_size + 1),
      strip.text = element_text(face = "bold"),
      strip.background = element_rect(fill = "grey93", colour = "grey70"),
      legend.key.size = unit(0.9, "lines"),
      panel.grid.minor = element_blank()
    )
}

#' A labelled empty frame standing in for a panel whose input is not available.
#' Renders as an explicit PLACEHOLDER rather than an empty space, and lists the inputs
#' required, so the gap is legible in the figure itself and not only in the guide.
placeholder_panel <- function(title, needs, have = NULL, body = NULL) {
  txt <- paste0(
    "PLACEHOLDER\n\n", title, "\n",
    if (!is.null(have)) paste0("\nHave: ", have, "\n") else "",
    if (!is.null(body)) paste0("\n", body, "\n") else "",
    "\nNeeds:\n", paste0("  \u2022 ", needs, collapse = "\n")
  )
  ggplot() +
    annotate("rect", xmin = 0, xmax = 1, ymin = 0, ymax = 1,
             fill = "grey97", colour = "#C1443C", linewidth = 0.7, linetype = "22") +
    annotate("text", x = 0.035, y = 0.975, label = txt, hjust = 0, vjust = 1,
             size = 2.75, lineheight = 1.18, colour = "grey20") +
    scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
    labs(subtitle = title) +
    ms_theme() +
    theme(axis.text = element_blank(), axis.ticks = element_blank(),
          axis.title = element_blank(), panel.grid = element_blank(),
          panel.border = element_blank())
}

#' Assemble panels with patchwork and either save or print.
#' Saves only when options(manuscript.save_dir=) is set, which keeps the explore scripts
#' side-effect free when run by hand.
#' Hard-wrap a caption. plot.caption does NOT wrap, so a long line silently runs off the
#' right edge of the saved PNG -- invisible until you look at the render. Wrap each
#' explicitly-newlined segment separately so intentional breaks survive.
wrap_caption <- function(x, width = 150) {
  if (is.null(x)) return(NULL)
  segs <- strsplit(x, "\n", fixed = TRUE)[[1]]
  paste(vapply(segs, function(s) paste(strwrap(s, width = width), collapse = "\n"),
               character(1)), collapse = "\n")
}

assemble_figure <- function(panels, design, heights = NULL, title = NULL, caption = NULL,
                            file, width = 11, height = 11, dpi = 200,
                            caption_width = 150) {
  caption <- wrap_caption(caption, width = caption_width)
  save_dir <- getOption("manuscript.save_dir", NULL)
  if (!requireNamespace("patchwork", quietly = TRUE)) {
    message("patchwork not installed; printing panels individually.")
    for (p in panels) print(p)
    return(invisible(NULL))
  }
  lab <- LETTERS[seq_along(panels)]
  panels <- Map(function(p, l) p + labs(tag = l), panels, lab)
  g <- patchwork::wrap_plots(panels, design = design, heights = heights) +
    patchwork::plot_annotation(
      title = title, caption = caption,
      theme = theme(
        plot.title = element_text(face = "bold", size = 14, margin = margin(b = 8)),
        plot.caption = element_text(hjust = 0, size = 8.5, colour = "grey30",
                                    margin = margin(t = 10))
      )
    ) &
    theme(plot.tag = element_text(face = "bold", size = 14))
  if (is.null(save_dir)) {
    print(g)
    message("Printed ", file, " (set options(manuscript.save_dir=) to write it).")
  } else {
    dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)
    out <- file.path(save_dir, file)
    ggsave(out, plot = g, width = width, height = height, dpi = dpi, bg = "white")
    message("Wrote ", out)
  }
  invisible(g)
}
