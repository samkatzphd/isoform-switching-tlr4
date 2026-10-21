#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------------------
# MANUSCRIPT FIGURE 2 -- the isoform-usage response to LPS. STANDALONE.
#
#   Rscript explore/manuscript_fig2_explore.R        # prints the figure, writes nothing
#
# Reads only results/tables/manuscript/panel_f2*.csv. See explore/_manuscript_common.R.
#
# THE FRAMING RULE THIS FIGURE ENFORCES
#   H_HT is primary macrophages from three donors on three days, and it is the ONLY dataset
#   importRdata batch-corrected. That halves its replicate IF noise, which is most of why it
#   shows a 2.07% switching rate against T_HT's 0.87%. Those rates are NOT two measurements
#   of one quantity. Every panel here facets; nothing overlays the two datasets on a shared
#   axis, and the rate is reported per dataset with that caveat stated on the panel.
# ---------------------------------------------------------------------------------------

suppressPackageStartupMessages({ library(dplyr); library(ggplot2) })

source(file.path("explore", "_manuscript_common.R"))
P <- panel_reader()

# --- 2A: how much of the transcriptome switches -----------------------------------------
f2a <- P("f2a")
stopifnot(!any(f2a$comparable_across_datasets))
f2a_long <- f2a |>
  transmute(
    dataset = pretty_ds(dataset),
    `genes quantified` = n_genes,
    `switching genes` = n_switching_genes,
    `isoforms quantified` = n_isoforms_quantified,
    `switching isoforms` = n_switching_isoforms,
    pct_switching_genes
  ) |>
  tidyr::pivot_longer(c(`genes quantified`, `switching genes`,
                        `isoforms quantified`, `switching isoforms`),
                      names_to = "measure", values_to = "n") |>
  mutate(
    level = ifelse(grepl("gene", measure), "gene", "isoform"),
    kind = ifelse(grepl("switching", measure), "switching", "quantified"),
    level = factor(level, levels = c("gene", "isoform"))
  )

p_a <- ggplot(f2a_long, aes(x = level, y = n, fill = kind)) +
  geom_col(width = 0.66, position = position_dodge2(padding = 0.1)) +
  geom_text(aes(label = scales::label_comma()(n)),
            position = position_dodge2(width = 0.66, padding = 0.1),
            vjust = -0.3, size = 2.7) +
  facet_wrap(~dataset, nrow = 1) +
  scale_fill_manual(values = c(quantified = "#C9D6E3", switching = "#1F4E6B"), name = NULL) +
  scale_y_log10(labels = scales::label_comma(), expand = expansion(mult = c(0, 0.13))) +
  labs(subtitle = "Scale of the isoform-usage response",
       x = NULL, y = "count (log scale)") +
  ms_theme() + theme(legend.position = "top")

# --- 2B: the distribution of effect sizes -----------------------------------------------
f2b <- P("f2b") |> mutate(dataset = pretty_ds(dataset))
rates <- f2a |> transmute(dataset = pretty_ds(dataset),
                          lab = paste0(n_switching_isoforms, " switching isoforms\n",
                                       pct_switching_genes, "% of genes switch"))
p_b <- ggplot(f2b, aes(x = dIF)) +
  geom_histogram(binwidth = 0.025, fill = "#4F7CAC", colour = "white", linewidth = 0.18) +
  geom_vline(xintercept = c(-0.15, 0.15), linetype = "22", colour = "#C1443C") +
  geom_text(data = rates, aes(x = 0, y = Inf, label = lab), inherit.aes = FALSE,
            vjust = 1.15, size = 2.8, lineheight = 1.1, colour = "grey25") +
  facet_wrap(~dataset, nrow = 1) +
  scale_x_continuous(breaks = seq(-0.8, 0.8, 0.4)) +
  # headroom so the per-dataset count annotation clears the tallest bars
  scale_y_continuous(expand = expansion(mult = c(0, 0.22))) +
  labs(subtitle = "Effect size of switching isoforms (dashed line = |dIF| threshold)",
       x = "dIF (isoform fraction, +LPS \u2212 \u2212LPS)", y = "isoforms") +
  ms_theme()

# --- 2C: why there is an abundance floor -------------------------------------------------
# dIF is a ratio, so it is the GENE's abundance that sets how precisely the isoform fraction
# can be estimated. Below ~12 a threshold-sized switch sits within ~2 SD of replicate noise.
f2c <- P("f2c") |> mutate(dataset = pretty_ds(dataset))
f2c_floor <- P("f2c_floor") |> mutate(dataset = pretty_ds(dataset))
applied <- unique(f2c_floor$consensus_floor_applied)
p_c <- ggplot(f2c, aes(x = expr_median, y = if_sd_median)) +
  geom_line(colour = "#1F4E6B", linewidth = 0.7) +
  geom_point(aes(size = snr_at_threshold), colour = "#1F4E6B", alpha = 0.85) +
  geom_vline(xintercept = applied, linetype = "22", colour = "#C1443C") +
  annotate("text", x = applied, y = Inf, label = paste0(" floor = ", applied),
           hjust = 0, vjust = 1.6, size = 2.8, colour = "#C1443C") +
  facet_wrap(~dataset, nrow = 1) +
  scale_x_log10(labels = scales::label_comma()) +
  scale_size_continuous(range = c(0.9, 3.4), name = "SNR at\n|dIF| = 0.15") +
  labs(subtitle = "Replicate isoform-fraction noise falls with gene abundance",
       x = "gene expression, bin median (log scale)",
       y = "within-condition IF SD (median)") +
  ms_theme() + theme(legend.position = "right")

# --- 2D: sensitivity to the threshold ----------------------------------------------------
f2d <- P("f2d") |>
  mutate(dataset = pretty_ds(dataset)) |>
  tidyr::pivot_longer(c(n_sig_isoforms, n_sig_genes),
                      names_to = "level", values_to = "n") |>
  mutate(level = factor(level, levels = c("n_sig_genes", "n_sig_isoforms"),
                        labels = c("genes", "isoforms")))
chosen <- unique(P("f2d")$min_abs_dif[P("f2d")$is_chosen])
p_d <- ggplot(f2d, aes(x = min_abs_dif, y = n, colour = level)) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 1.5) +
  geom_vline(xintercept = chosen, linetype = "22", colour = "#C1443C") +
  annotate("text", x = chosen, y = Inf, label = paste0(" chosen = ", chosen),
           hjust = 0, vjust = 1.6, size = 2.8, colour = "#C1443C") +
  facet_wrap(~dataset, nrow = 1) +
  scale_colour_manual(values = c(genes = "#1F4E6B", isoforms = "#F2B134"), name = NULL) +
  labs(subtitle = "Calls retained across the |dIF| threshold",
       x = "min |dIF| threshold", y = "significant calls") +
  ms_theme() + theme(legend.position = "top")

# --- assemble ----------------------------------------------------------------------------
assemble_figure(
  panels = list(p_a, p_b, p_c, p_d),
  design = "AB\nCD",
  heights = c(1, 1),
  title = "Figure 2. The isoform-usage response to LPS on the HT reference",
  caption = paste0(
    "All panels facet by dataset and none overlays them: H_HT is batch-corrected primary ",
    "macrophages from three donors, so its switching rate and dIF distribution are not ",
    "comparable to T_HT's (AGENT_CONTEXT fact 4).\n",
    "The abundance floor in C is applied to GENE expression, not isoform expression, and ",
    "removes 14-36% of raw switching calls by design. No UBL5 knockout data enters this ",
    "figure."),
  file = "fig2_response.png", width = 13.0, height = 9.4
)
