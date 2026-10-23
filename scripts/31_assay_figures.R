#!/usr/bin/env Rscript
# 31_assay_figures.R -- three figures for the assay-design stage.
#
#   assay_design_funnel.png        how 44 candidates become verified reagents, named
#   assay_assignment_by_candidate.png  which assay each candidate needs, and why
#   oligo_specificity.png          where every oligo actually binds, against 296,613
#
# Reads only CSVs written by 27-30. Writes only PNGs.

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2); library(patchwork)
})

ROOT <- normalizePath(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), ".."))
setwd(ROOT)
dir.create("results/figures/assays", recursive = TRUE, showWarnings = FALSE)

BASE <- 8
theme_pub <- function() {
  theme_minimal(base_size = BASE) +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major.y = element_blank(),
          plot.title = element_text(size = BASE, hjust = 0, face = "plain"),
          plot.subtitle = element_text(size = BASE - 1, colour = "grey35", hjust = 0),
          axis.text = element_text(size = BASE - 2),
          legend.text = element_text(size = BASE - 1),
          legend.title = element_text(size = BASE - 1),
          strip.text = element_text(size = BASE - 1))
}
INK <- "#1f2d3d"; FOCAL <- "#0b6fa4"; WARN <- "#c2531a"; MUTE <- "#9aa5b1"

assign  <- read.csv("results/tables/assay_assignment.csv", stringsAsFactors = FALSE)
order_s <- read.csv("results/tables/assay_order_sheet.csv", stringsAsFactors = FALSE)
prep    <- read.csv("results/tables/pair_selection_report.csv", stringsAsFactors = FALSE)
spec    <- read.csv("results/tables/assay_specificity.csv", stringsAsFactors = FALSE)
regions <- read.csv("results/tables/discriminating_regions.csv", stringsAsFactors = FALSE)

# ------------------------------------------------------------------ Figure A
# The funnel has to end in names, not a count: the point of the stage is that
# specific genes now have specific, checked reagents waiting to be ordered.
qpcr_pass <- prep %>% filter(grepl("^qPCR", pair), passes == "True") %>% pull(gene_name)
race_ok <- order_s %>%
  filter(grepl("^RACE5", oligo_role), suppressWarnings(as.integer(offtarget_loci)) == 0) %>%
  count(gene_name) %>% filter(n >= 2) %>% pull(gene_name)      # both GSPs clean

stages <- tibble::tibble(
  stage = factor(c("switching genes\n(donor-blocked refit)", "direction-supported\nshortlist",
                   "survive locus-fusion\nscreen", "oligos designed\n(5 alternatives each)",
                   "oligos selected after\ntranscriptome scan", "zero off-target loci"),
                 levels = rev(c("switching genes\n(donor-blocked refit)", "direction-supported\nshortlist",
                                "survive locus-fusion\nscreen", "oligos designed\n(5 alternatives each)",
                                "oligos selected after\ntranscriptome scan", "zero off-target loci"))),
  # 177 = 86 T-private + 83 H-private + 8 shared, from the donor-blocked refit
  n = c(177, 63, 44, nrow(spec), nrow(order_s),
        sum(suppressWarnings(as.integer(order_s$offtarget_loci)) == 0, na.rm = TRUE)),
  unit = c("genes", "genes", "genes", "oligos", "oligos", "oligos"))

funnel_panel <- function(d, fill, title, xlab) {
  ggplot(d, aes(n, stage)) +
    geom_col(fill = fill, width = 0.6) +
    geom_text(aes(label = format(n, big.mark = ",")), hjust = -0.2, size = 2.5, colour = INK) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.22))) +
    labs(title = title, x = xlab, y = NULL) +
    theme_pub()
}
# Genes and oligos are different units; they do not share an axis.
pA  <- funnel_panel(droplevels(filter(stages, unit == "genes")), FOCAL,
                    "Candidates", "genes")
pA2 <- funnel_panel(droplevels(filter(stages, unit == "oligos")), MUTE,
                    "Reagents", "oligos")

ready <- sort(union(qpcr_pass, race_ok))
lab <- paste(strwrap(paste(ready, collapse = "  "), width = 64), collapse = "\n")
pB <- ggplot() +
  annotate("text", x = 0, y = 1, hjust = 0, vjust = 1, size = 2.5, colour = INK,
           label = paste0(length(ready), " candidates have a primary assay whose oligos\n",
                          "passed the scan with no off-target locus:\n\n", lab)) +
  xlim(0, 1) + ylim(0, 1) + theme_void()

figA <- (pA | pA2) / pB +
  plot_layout(heights = c(1, 0.55)) +
  plot_annotation(
    title = "From switching genes to reagents that survived a specificity check",
    theme = theme(plot.title = element_text(size = BASE + 1, hjust = 0)))
ggsave("results/figures/assays/assay_design_funnel.png", figA,
       width = 7.6, height = 4.6, dpi = 300, bg = "white")

# ------------------------------------------------------------------ Figure B
ev <- assign %>%
  mutate(gene_name = reorder(gene_name, abs(own_dIF)),
         event = sub(" \\(TSS\\)", "", primary_event),
         assay = sub(" \\(\\+ 5' RACE to fix the TSS\\)", "", primary_assay),
         nested = ifelse(is_nested == "TRUE", "nested in partner", "resolvable alone"))

pC <- ggplot(ev, aes(abs(own_dIF), gene_name)) +
  geom_segment(aes(x = 0, xend = abs(own_dIF), yend = gene_name), colour = MUTE, linewidth = 0.3) +
  geom_point(aes(colour = assay, shape = nested), size = 1.9) +
  scale_colour_manual(values = c("5' RACE" = FOCAL, "unique-region qPCR / ddPCR" = "#1b7f5f",
                                 "3'-region qPCR / ddPCR" = "#6f5fa8",
                                 "junction-spanning qPCR / ddPCR" = "#b08900",
                                 "intron-spanning qPCR / ddPCR" = WARN, "3' RACE" = "#8a5a44"),
                      name = "primary assay") +
  scale_shape_manual(values = c("nested in partner" = 1, "resolvable alone" = 16), name = NULL) +
  labs(title = "Most candidates need an assay that reads the 5' end, not a short amplicon",
       subtitle = "29 of 44 switches are alternative-first-exon events; 10 isoforms are contained within their partner",
       x = "|dIF| in the discovering dataset", y = NULL) +
  theme_pub() + theme(panel.grid.major.y = element_line(colour = "grey92", linewidth = 0.3),
                      legend.position = "right", legend.key.height = unit(8, "pt"))
ggsave("results/figures/assays/assay_assignment_by_candidate.png", pC,
       width = 7.6, height = 6.4, dpi = 300, bg = "white")

# ------------------------------------------------------------------ Figure C
sp <- spec %>%
  mutate(n = pmax(as.integer(transcripts_full_match), 0L),
         role = dplyr::recode(oligo_role,
                              qPCR_F = "qPCR forward", qPCR_R = "qPCR reverse",
                              longread_F = "long-read forward", longread_R = "long-read reverse",
                              RACE5_GSP2_inner = "5' RACE inner GSP",
                              RACE5_GSP1_outer = "5' RACE outer GSP",
                              RACE3_GSP2_inner = "3' RACE inner GSP",
                              RACE3_GSP1_outer = "3' RACE outer GSP")) %>%
  filter(n > 0)
sel <- order_s %>% transmute(key = paste(isoform_id, oligo_role, alt), selected = TRUE)
sp <- sp %>% mutate(key = paste(isoform_id, oligo_role, alt)) %>%
  left_join(sel, by = "key") %>% mutate(selected = !is.na(selected))

extremes <- sp %>% filter(selected) %>% slice_max(n, n = 4)

sp <- sp %>% mutate(off = suppressWarnings(as.integer(offtarget_loci)),
                    status = ifelse(!selected, "alternative not chosen",
                             ifelse(off > 0, "selected: binds another locus",
                                    "selected: own locus only")))
extremes <- sp %>% filter(selected, off > 0) %>% slice_max(off, n = 4)

pD <- ggplot(sp, aes(n, role)) +
  geom_vline(xintercept = 1, colour = INK, linewidth = 0.4, linetype = "22") +
  geom_jitter(data = filter(sp, status == "alternative not chosen"),
              height = 0.18, size = 0.5, colour = MUTE, alpha = 0.45) +
  geom_jitter(data = filter(sp, status == "selected: own locus only"),
              height = 0.18, size = 1.2, colour = FOCAL) +
  geom_jitter(data = filter(sp, status == "selected: binds another locus"),
              height = 0.18, size = 1.6, colour = WARN) +
  geom_text(data = extremes, aes(label = gene_name), hjust = 1.2, size = 2.2, colour = WARN) +
  annotate("text", x = 1, y = 8.7, label = "binds one transcript\nin 296,613", hjust = -0.08,
           vjust = 1, size = 2.2, colour = INK) +
  scale_x_log10(breaks = c(1, 10, 100, 1000, 10000),
                labels = c("1", "10", "100", "1k", "10k")) +
  coord_cartesian(clip = "off") +
  labs(title = "Binding more than one transcript is expected; binding another gene is not",
       subtitle = paste0("Grey, alternatives not chosen. Blue, selected and confined to its own locus (",
                         sum(sp$status == "selected: own locus only"), " oligos). Orange, selected but binds a second locus (",
                         sum(sp$status == "selected: binds another locus"), ")."),
       x = "transcripts bound (full-length exact match, log scale)", y = NULL) +
  theme_pub() + theme(panel.grid.major.y = element_line(colour = "grey92", linewidth = 0.3),
                      plot.margin = margin(6, 26, 6, 6))
ggsave("results/figures/assays/oligo_specificity.png", pD,
       width = 7.4, height = 4.4, dpi = 300, bg = "white")

cat("figures written\n")
cat("  ready candidates:", length(ready), "\n")
cat("  oligos plotted  :", nrow(sp), "\n")
