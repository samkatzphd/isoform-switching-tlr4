#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------------------
# MANUSCRIPT FIGURE 4 -- candidate selection. STANDALONE.
#
#   Rscript explore/manuscript_fig4_explore.R        # prints the figure, writes nothing
#
# Reads only results/tables/manuscript/panel_f4*.csv. See explore/_manuscript_common.R.
#
# TWO OF THESE FOUR PANELS ARE PLACEHOLDERS, AND THE REASON IS THE SAME
#   Everything about the ANATOMY of a switch -- which AS event it is, how many bases
#   discriminate the two isoforms, which annotated protein features change -- was computed by
#   scripts 13/15/16 from U_T.gff3 and transcripts_U-T.fa, i.e. in the UT annotation space.
#   All 84 candidate isoform ids resolve in T_UT/U_UT but only 42 of 80 in T_HT, so that work
#   cannot be carried into an HT-only manuscript. Reproducing it here needs H_T.gff3 and
#   transcripts_H-T.fa. Panels C and D stay empty until those are in hand.
#
# WHAT PANEL B IS AND IS NOT
#   q1_candidates has three arms: anchor = T_HT, reannot = T_UT, cross = H_HT. Script 21
#   drops every *_reannot column, leaving the two HT arms. H_HT enters as PRESENCE/ABSENCE
#   corroboration in an independent cell system. Its dIF is deliberately NOT plotted against
#   the anchor's: H is batch-corrected primary cells and the two effect sizes are not on a
#   common scale (AGENT_CONTEXT fact 4).
# ---------------------------------------------------------------------------------------

suppressPackageStartupMessages({ library(dplyr); library(ggplot2) })

source(file.path("explore", "_manuscript_common.R"))
P <- panel_reader()

# --- 4A: the selection funnel ------------------------------------------------------------
f4a <- P("f4a") |>
  mutate(dataset = pretty_ds(dataset_label),
         stage = stats::reorder(stage, stage_order))
f4a_real <- f4a |> filter(status == "real")
f4a_ph   <- f4a |> filter(status == "placeholder")
p_a <- ggplot(f4a_real, aes(x = stage, y = n_genes)) +
  geom_col(width = 0.68, fill = "#1F4E6B") +
  geom_text(aes(label = scales::label_comma()(n_genes)), vjust = -0.35, size = 2.8) +
  # the un-run stage drawn as an explicitly empty, dashed bar rather than omitted
  geom_blank(data = f4a_ph, aes(x = stage, y = 1)) +
  # ymin = 1, not 0: the y axis is log10, where 0 maps to -Inf and ggplot warns about it.
  geom_rect(data = f4a_ph,
            aes(xmin = as.numeric(stage) - 0.34, xmax = as.numeric(stage) + 0.34,
                ymin = 1, ymax = Inf),
            inherit.aes = FALSE, fill = "grey96", colour = "#C1443C",
            linetype = "22", linewidth = 0.5) +
  geom_text(data = f4a_ph, aes(x = stage, y = 1,
                               label = "not run\n(needs HT-space\nassayability screen)"),
            inherit.aes = FALSE, vjust = -0.3, size = 2.5, colour = "#C1443C") +
  facet_wrap(~dataset, nrow = 1) +
  scale_y_log10(labels = scales::label_comma(), expand = expansion(mult = c(0, 0.14))) +
  labs(subtitle = "Candidate selection funnel (genes)", x = NULL, y = "genes (log scale)") +
  ms_theme() + theme(axis.text.x = element_text(angle = 24, hjust = 1))

# --- 4B: cross-system corroboration ------------------------------------------------------
f4b <- P("f4b")
CORR_COLS <- c(
  "both T_HT and H_HT" = "#1F6F54",
  "T_HT only (testable in H_HT, not called)" = "#F2B134",
  "T_HT only (not testable in H_HT)" = "#9E9E9E",
  "H_HT only" = "#9DBEDC",
  "neither" = "#E0E0E0"
)
# The full category names are too long for a legend; keep them in the CSV and the axis of
# panel B, and use these in the colour key so nothing gets clipped off the figure edge.
CORR_SHORT <- c(
  "both T_HT and H_HT" = "both systems",
  "T_HT only (testable in H_HT, not called)" = "T only (H testable)",
  "T_HT only (not testable in H_HT)" = "T only (H not testable)",
  "H_HT only" = "H only",
  "neither" = "neither"
)
corr_counts <- f4b |>
  count(corroboration, name = "n") |>
  mutate(corroboration = factor(corroboration, levels = names(CORR_COLS))) |>
  filter(!is.na(corroboration))
p_b <- ggplot(corr_counts, aes(x = stats::reorder(corroboration, n), y = n,
                               fill = corroboration)) +
  geom_col(width = 0.68) +
  geom_text(aes(label = n), hjust = -0.25, size = 2.9) +
  coord_flip() +
  scale_fill_manual(values = CORR_COLS, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.16))) +
  scale_x_discrete(labels = function(x) sub(" \\(", "\n(", x)) +
  labs(subtitle = "Cross-system corroboration of anchor candidates",
       x = NULL, y = "genes") +
  ms_theme()

# Top anchor genes, ordered by the anchor's own effect size and coloured by whether primary
# macrophages corroborate them. dif_cross is deliberately not drawn (see header).
top_n_show <- 18
f4b_top <- f4b |>
  filter(switches_anchor) |>
  arrange(desc(abs(dif_anchor))) |>
  slice_head(n = top_n_show) |>
  mutate(gene_name = stats::reorder(gene_name, abs(dif_anchor)),
         corroboration = factor(corroboration, levels = names(CORR_COLS)))
p_b2 <- ggplot(f4b_top, aes(x = abs(dif_anchor), y = gene_name, colour = corroboration)) +
  geom_segment(aes(x = 0, xend = abs(dif_anchor), y = gene_name, yend = gene_name),
               linewidth = 1.1) +
  geom_point(size = 2.4) +
  # drop = TRUE: only categories actually present in the top-N get a key, otherwise the
  # legend shows colourless entries for categories with no points.
  scale_colour_manual(values = CORR_COLS, labels = CORR_SHORT, name = NULL,
                      drop = TRUE) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.06))) +
  labs(subtitle = paste0("Top ", top_n_show, " anchor candidates"),
       x = "|dIF| in T_HT (anchor)", y = NULL) +
  ms_theme() +
  theme(legend.position = "bottom", legend.text = element_text(size = 7.2),
        legend.key.size = unit(0.75, "lines")) +
  guides(colour = guide_legend(ncol = 1))

# --- 4C / 4D: placeholders ----------------------------------------------------------------
f4c <- P("f4c")
p_c <- placeholder_panel(
  title = "Switch anatomy by AS event class",
  needs = c("H_T.gff3 (exon coordinates), on /Volumes/Expansion",
            "re-run of scripts/15_isoform_characterization.R in HT space"),
  have = "event classes enumerated, counts not derivable",
  body = paste(
    "Script 15 answered this in the UT space: 54 of 63",
    "candidate pairs differ by a distinct FIRST EXON, i.e.",
    "most are promoter switches, not splicing switches.",
    "That number cannot be carried over -- it was computed",
    "against U_T.gff3, and the candidate ids do not resolve",
    "in the HT annotation.",
    sep = "\n")
)
p_d <- placeholder_panel(
  title = "Candidate shortlist with platform class",
  needs = c("HT-space re-run of 13 (assayability), 15 (anatomy), 16 (functional sites)",
            "H_T.gff3 + transcripts_H-T.fa, on /Volumes/Expansion"),
  have = "selection through |dIF| and cross-system corroboration (panels A, B)",
  body = paste(
    "The existing shortlist -- RBKS, ATRN, IRAK3, NCOA7,",
    "IL23A, WDR33 and the platform split of 16 Illumina-",
    "findable vs 18 PacBio-required -- is UT-space output.",
    "All 84 candidate isoform ids resolve in T_UT/U_UT but",
    "only 42 of 80 in T_HT, so it is not admissible here.",
    sep = "\n")
)

# --- assemble ------------------------------------------------------------------------------
assemble_figure(
  panels = list(p_a, p_b, p_b2, p_c, p_d),
  # patchwork requires a RECTANGULAR design string -- every row the same number of
  # characters. A = funnel (wide), B = corroboration counts, C = top anchor genes (tall,
  # spans two rows), D/E = the two placeholders.
  design = "AAAC\nBBBC\nDDEE",
  heights = c(1.0, 0.85, 1.05),
  title = "Figure 4. Candidate selection on the HT reference",
  caption = paste0(
    "Panels D and E are placeholders. Everything describing the anatomy of a switch was ",
    "computed in the UT annotation space (scripts 13/15/16, from U_T.gff3 and ",
    "transcripts_U-T.fa) and does not transfer: all 84 candidate isoform ids resolve in ",
    "T_UT/U_UT but only 42 of 80 in T_HT.\n",
    "In B and C, H_HT contributes presence/absence corroboration only; its effect sizes are ",
    "not plotted against the anchor's because the two are not on a common scale. No UBL5 ",
    "knockout data enters this figure."),
  file = "fig4_candidates.png", width = 13.0, height = 13.6
)
