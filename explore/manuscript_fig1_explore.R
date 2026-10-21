#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------------------
# MANUSCRIPT FIGURE 1 -- experimental design and sequencing. STANDALONE.
#
#   Rscript explore/manuscript_fig1_explore.R        # prints the figure, writes nothing
#
# Reads only results/tables/manuscript/panel_f1*.csv. No bootstrap, no config, no .rds, no
# external drive. Needs dplyr + ggplot2 (+ patchwork to assemble; falls back to printing the
# panels one at a time without it).
#
# To save instead of print, set the option before sourcing -- this is what
# scripts/22_manuscript_figures.R does:
#   options(manuscript.save_dir = "results/figures/manuscript")
#
# SCOPE: HT reference, H and T libraries only. No UBL5 knockout data anywhere in this figure.
# ---------------------------------------------------------------------------------------

suppressPackageStartupMessages({ library(dplyr); library(ggplot2) })

source(file.path("explore", "_manuscript_common.R"))
P <- panel_reader()

# --- 1A: design schematic. PLACEHOLDER. --------------------------------------------------
f1a <- P("f1a")
p_a <- placeholder_panel(
  title = "Experimental design",
  needs = c("PMA differentiation protocol", "LPS dose", "LPS timepoint",
            "library chemistry (polyA vs ribo-depleted)", "5' cap verification"),
  have = paste0(sum(f1a$known), " of ", nrow(f1a), " design elements known"),
  body = paste(
    "THP-1 (T): one clonal line; T1-T3 matured,",
    "treated and extracted separately; paired -/+ LPS.",
    "Primary MDM (H): 3 donors, 3 days; paired -/+ LPS.",
    "H is the only batch-corrected dataset.",
    sep = "\n")
)

# --- 1B: per-library depth and detection -------------------------------------------------
f1b <- P("f1b")
p_b <- ggplot(f1b, aes(x = total_expected_count / 1e6, y = library, fill = lps)) +
  geom_col(width = 0.72) +
  facet_grid(rows = vars(cell_system), scales = "free_y", space = "free_y") +
  scale_fill_manual(values = LPS_COLS, name = "LPS") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(subtitle = "Library depth (12 HT libraries)",
       x = "RSEM expected counts (millions)", y = NULL) +
  ms_theme() + theme(legend.position = "top")

# --- 1C: merged HT annotation composition by class code ----------------------------------
# "=" dominates by two orders of magnitude, so this is log-scaled. The point of the panel is
# what the PacBio assembly CONTRIBUTED, so fill is the detectability tier, not the dataset.
f1c <- P("f1c") |>
  mutate(
    detect_tier = factor(detect_tier, levels = TIER_LEVELS),
    class_code = stats::reorder(class_code, -n, FUN = sum),
    dataset_label = pretty_ds(dataset_label)
  )
# NEVER stack on a log scale: geom_col's default position adds the TRANSFORMED values, so
# two stacked bars multiply instead of summing and the axis runs to ~1e9 for a quantity whose
# true maximum is ~1.5e5. Only class "=" carries two tiers (carried ENST models plus PacBio
# models whose chain matches), so dodging shows exactly the comparison that matters.
p_c <- ggplot(f1c, aes(x = class_code, y = n, fill = detect_tier)) +
  geom_col(width = 0.74, position = position_dodge2(preserve = "single", padding = 0.12)) +
  facet_wrap(~dataset_label, nrow = 1) +
  scale_fill_manual(values = TIER_COLS, labels = TIER_LABELS, name = NULL, drop = FALSE) +
  scale_y_log10(labels = scales::label_comma(),
                expand = expansion(mult = c(0, 0.06))) +
  labs(subtitle = "Merged HT annotation: all quantified isoforms by gffcompare class code",
       x = "class code", y = "isoforms (log scale)") +
  ms_theme() +
  theme(legend.position = "right",
        panel.grid.major.x = element_blank())

# --- 1D: the primary objects are pre-reduced ---------------------------------------------
f1d <- P("f1d") |>
  select(dataset_label, n_genes, n_genes_unfiltered,
         pct_primary = pct_genes_gene_q_significant,
         pct_unfiltered = pct_genes_significant_unfiltered) |>
  tidyr::pivot_longer(c(pct_primary, pct_unfiltered),
                      names_to = "layer", values_to = "pct") |>
  mutate(
    layer = factor(layer, levels = c("pct_primary", "pct_unfiltered"),
                   labels = c("primary ISA object", "unfiltered context")),
    dataset_label = pretty_ds(dataset_label)
  )
p_d <- ggplot(f1d, aes(x = layer, y = pct, fill = layer)) +
  geom_col(width = 0.6) +
  geom_text(aes(label = paste0(round(pct), "%")), vjust = -0.35, size = 3.1) +
  facet_wrap(~dataset_label, nrow = 1) +
  scale_fill_manual(values = c("primary ISA object" = "#9E9E9E",
                               "unfiltered context" = "#4F7CAC"), guide = "none") +
  scale_y_continuous(limits = c(0, 118), breaks = c(0, 25, 50, 75, 100),
                     expand = expansion(mult = c(0, 0))) +
  labs(subtitle = "Genes passing the gene-level q cutoff, by input layer",
       x = NULL, y = "% of genes significant") +
  ms_theme() +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))

# --- assemble ----------------------------------------------------------------------------
assemble_figure(
  panels = list(p_a, p_b, p_c, p_d),
  design = "AABB\nCCCC\nDDDD",
  heights = c(1.15, 1.0, 0.95),
  title = "Figure 1. Experimental design and sequencing of the HT reference",
  caption = paste0(
    "Panel A is a placeholder: the protocol parameters it needs are the standing wet-lab ",
    "questions. Panels B-D are derived from committed data.\n",
    "Scope is the HT reference transcriptome and the H and T libraries; no UBL5 knockout ",
    "data enters this figure."),
  file = "fig1_design.png", width = 11.5, height = 12.2
)
