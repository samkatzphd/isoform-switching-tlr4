#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------------------
# MANUSCRIPT FIGURE 3 -- what the PacBio assembly actually added. STANDALONE.
#
#   Rscript explore/manuscript_fig3_explore.R        # prints the figure, writes nothing
#
# Reads only results/tables/manuscript/panel_f3*.csv. See explore/_manuscript_common.R.
#
# THE ARGUMENT THIS FIGURE MAKES
#   The pipeline's `is_novel_pacbio` flag is set from the transcript ID prefix, so it records
#   which assembler emitted the model, not whether its structure is new. gffcompare's
#   class_code is the structural statement: "=" means the intron chain exactly matches a named
#   Ensembl transcript, so short-read Illumina quantified against Ensembl would have had that
#   transcript in its index. Counting provenance as discovery overstates structural novelty
#   by 34-43 percentage points. Panel B is the correction; A and C show where it comes from.
#   Detail in docs/REVIEW_CHANGES.md §0j.
# ---------------------------------------------------------------------------------------

suppressPackageStartupMessages({ library(dplyr); library(ggplot2) })

source(file.path("explore", "_manuscript_common.R"))
P <- panel_reader()

# --- 3A: switching isoforms by detectability tier ---------------------------------------
f3a <- P("f3a") |>
  mutate(dataset = pretty_ds(dataset_label),
         detect_tier = factor(detect_tier, levels = TIER_LEVELS))
totals <- f3a |> group_by(dataset) |> summarize(tot = sum(n), .groups = "drop")
p_a <- ggplot(f3a, aes(x = dataset, y = n, fill = detect_tier)) +
  geom_col(width = 0.6) +
  geom_text(data = totals, aes(x = dataset, y = tot, label = paste0("n = ", tot)),
            inherit.aes = FALSE, vjust = -0.45, size = 2.9) +
  scale_fill_manual(values = TIER_COLS, labels = TIER_LABELS, name = NULL, drop = FALSE) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(subtitle = "Switching isoforms by what platform could have seen them",
       x = NULL, y = "switching isoforms") +
  ms_theme() + theme(legend.position = "right")

# --- 3B: the three novelty definitions ---------------------------------------------------
f3b <- P("f3b") |>
  mutate(dataset = pretty_ds(dataset_label),
         definition = factor(definition,
                             levels = c("published (ID prefix)", "structural",
                                        "strict (c excluded)", "Illumina-findable")))
over <- P("f3b") |>
  distinct(dataset_label, overcount_pp) |>
  mutate(dataset = pretty_ds(dataset_label),
         lab = paste0("overcount: ", round(overcount_pp, 1), " pp"))
p_b <- ggplot(f3b, aes(x = definition, y = pct, fill = definition)) +
  geom_col(width = 0.66) +
  geom_text(aes(label = paste0(round(pct, 1), "%")), vjust = -0.35, size = 2.8) +
  geom_text(data = over, aes(x = 2.5, y = Inf, label = lab), inherit.aes = FALSE,
            vjust = 1.4, size = 2.9, colour = "#C1443C", fontface = "bold") +
  facet_wrap(~dataset, nrow = 1) +
  scale_fill_manual(
    values = c("published (ID prefix)" = "#9E9E9E", "structural" = "#F2B134",
               "strict (c excluded)" = "#C1443C", "Illumina-findable" = "#4F7CAC"),
    guide = "none") +
  scale_y_continuous(limits = c(0, 92), expand = expansion(mult = c(0, 0))) +
  labs(subtitle = "\"Novel\" under three definitions",
       x = NULL, y = "% of switching isoforms") +
  ms_theme() + theme(axis.text.x = element_text(angle = 22, hjust = 1))

# --- 3C: class codes, all quantified vs switching ---------------------------------------
# If the assembly contributed something the LPS response actually uses, novel class codes
# should be over-represented among switching isoforms relative to the annotation at large.
f3c <- P("f3c") |>
  mutate(dataset = pretty_ds(dataset_label),
         set = factor(set, levels = c("all quantified", "switching")),
         class_code = stats::reorder(class_code, -pct, FUN = max))
p_c <- ggplot(f3c, aes(x = class_code, y = pct, fill = set)) +
  geom_col(width = 0.72, position = position_dodge2(preserve = "single", padding = 0.1)) +
  facet_wrap(~dataset, nrow = 1) +
  scale_fill_manual(values = c("all quantified" = "#C9D6E3", "switching" = "#1F4E6B"),
                    name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.06))) +
  labs(subtitle = "Class-code composition: the annotation at large vs the switching set",
       x = "class code", y = "% of set") +
  ms_theme() +
  theme(legend.position = "top", panel.grid.major.x = element_blank())

# --- 3D: effect size, structurally novel vs reference-matched ---------------------------
# Lead with the effect size, not a p-value: with hundreds of matched isoforms a p-value
# tracks n. Compare WITHIN a dataset only -- never across the facet boundary.
f3d <- P("f3d") |>
  mutate(dataset = pretty_ds(dataset_label),
         group = factor(ifelse(structurally_novel, "structurally novel", "reference-matched"),
                        levels = c("reference-matched", "structurally novel")))
p_d <- ggplot(f3d, aes(x = group, y = median_abs_dIF, colour = group)) +
  geom_errorbar(aes(ymin = q25, ymax = q75), width = 0.16, linewidth = 0.6) +
  geom_point(size = 3) +
  geom_text(aes(label = paste0("n = ", n)), vjust = -1.4, hjust = -0.18, size = 2.8,
            show.legend = FALSE) +
  facet_wrap(~dataset, nrow = 1) +
  scale_colour_manual(values = c("reference-matched" = "#4F7CAC",
                                 "structurally novel" = "#C1443C"), guide = "none") +
  labs(subtitle = "Effect size by structural novelty (median, IQR)",
       x = NULL, y = "|dIF| of switching isoforms") +
  ms_theme() + theme(axis.text.x = element_text(angle = 12, hjust = 1))

# --- assemble ----------------------------------------------------------------------------
assemble_figure(
  panels = list(p_a, p_b, p_c, p_d),
  design = "AB\nCD",
  heights = c(1, 1),
  title = "Figure 3. What the PacBio assembly added to the LPS isoform response",
  caption = paste0(
    "\"Published\" in B is the pipeline's own novelty rate, set from the transcript ID ",
    "prefix: it counts which assembler emitted the model, not whether the structure is in ",
    "the reference. \"Structural\" requires class_code != \"=\"; \"strict\" also excludes ",
    "class c (contained fragments), which are novel only if their 5' ends are real.\n",
    "Panel D compares groups within a dataset only. No UBL5 knockout data enters this ",
    "figure."),
  file = "fig3_novelty.png", width = 13.0, height = 9.8
)
