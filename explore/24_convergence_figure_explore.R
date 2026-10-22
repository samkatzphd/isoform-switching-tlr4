#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------------------
# How we get from per-system hits to a named follow-up shortlist. STANDALONE.
#
#   Rscript explore/24_convergence_figure_explore.R                  # prints, writes nothing
#   options(convergence.save_dir = "results/figures/refit")          # ... to write the PNG
#
# Reads only results/tables/refit_*.csv (script 23). Needs dplyr + ggplot2 (+ patchwork).
# Run from the project root. Writes nothing unless the option above is set.
#
# THE ARGUMENT, PANEL BY PANEL
#   A  Each system finds ~90 switching genes, but only 8 are called in both. Read as a
#      binary that looks like near-total failure to replicate.
#   B  It is not. Graded by what the other system actually shows on the same isoform, the
#      largest single class is "all three replicates agree in direction but the effect sits
#      below the |dIF| threshold" -- real support that a significance cutoff at n=3 discards.
#   C  And the direction agreement is not wishful thinking: ~69% against a 50% null, and for
#      T discoveries all-three-donor agreement in H beats the genome-wide background.
#   D  So the follow-up set is the direction-supported union, with names.
#
# TWO CAVEATS CARRIED ONTO THE FIGURE ITSELF
#   - Genes are XLOC loci from the merged assembly. The locus-fusion screen (script 18) ran
#     in UT space and has NOT been re-run here, so symbols confirmed fused there are flagged:
#     a "fused locus" is two neighbouring genes merged, and its dIF is a ratio of two genes,
#     not an isoform switch.
#   - Rows whose "symbol" is a bare ENSG/XLOC identifier are not resolved gene symbols and
#     are excluded from the named panel rather than printed as if they were genes.
# ---------------------------------------------------------------------------------------

suppressPackageStartupMessages({ library(dplyr); library(ggplot2) })

TBL <- "results/tables"
rd <- function(stem) {
  f <- file.path(TBL, paste0(stem, ".csv"))
  if (!file.exists(f)) stop("Missing ", f, " -- run: Rscript scripts/23_donor_blocked_refit.R")
  utils::read.csv(f, stringsAsFactors = FALSE)
}
dual <- rd("refit_dual_anchor")
summ <- rd("refit_summary")
conc <- rd("refit_selconf_concordance")
# Fusion screen (script 25). Required: without it the named panel would advertise genes the
# screen has since disqualified, which is worse than showing a shorter list.
scr <- rd("followup_shortlist_fusion_screened")

PRIMARY <- "blocked (nested replicate)"
PRETTY  <- c(H_HT_supplementary = "Primary MDM (H)", T_HT_anchor = "THP-1 (T)")
COL_SYS <- c("Primary MDM (H)" = "#C1443C", "THP-1 (T)" = "#4F7CAC")
COL_SET <- c("shared" = "#1F6F54", "T-private" = "#4F7CAC", "H-private" = "#C1443C")

# Symbols confirmed in script 18 (UT space) to be two neighbouring genes merged into one
# XLOC. The equivalent screen has not been run on the HT assembly, so these are flagged as
# suspect rather than silently kept or silently dropped.
FUSION_SUSPECT <- c("PLD6", "DDX58", "IPCEF1", "THAP2", "MREG", "IFITM2", "CCL3")

# A real gene symbol, in this project's sense: not a bare accession.
is_real_symbol <- function(x) {
  !is.na(x) & nzchar(x) & !grepl("^(ENSG|ENST|XLOC|TCONS|LOC[0-9])", x)
}

th <- theme_bw(base_size = 11) +
  theme(plot.subtitle = element_text(face = "bold", size = 10.5, hjust = 0),
        strip.text = element_text(face = "bold"),
        panel.grid.minor = element_blank(),
        legend.key.size = unit(0.85, "lines"))

# "Direction-supported": all three replicates in the OTHER system move the same way as the
# discovery system. Tiers 3, 4 and 5 are exactly that, in increasing order of strength.
dual <- dual |>
  mutate(direction_supported = grepl("^[345] ", .data$tier),
         discovery_pretty = PRETTY[.data$discovery])

# --- A: hits per system, and the naive overlap -------------------------------------------
hits <- summ |>
  filter(.data$model == PRIMARY) |>
  transmute(system = PRETTY[.data$dataset_label], genes = .data$n_switching_genes)
shared_n <- dual |> filter(.data$set == "shared") |> distinct(.data$gene_id) |> nrow()
a_df <- bind_rows(
  hits |> mutate(cat = "switching genes"),
  tibble(system = c("Primary MDM (H)", "THP-1 (T)"), genes = shared_n,
         cat = "called in BOTH")
) |>
  mutate(cat = factor(.data$cat, levels = c("switching genes", "called in BOTH")))

p_a <- ggplot(a_df, aes(x = .data$system, y = .data$genes, fill = .data$system,
                        alpha = .data$cat)) +
  geom_col(width = 0.62, position = position_dodge2(padding = 0.1)) +
  geom_text(aes(label = .data$genes),
            position = position_dodge2(width = 0.62, padding = 0.1),
            vjust = -0.35, size = 3.1, show.legend = FALSE) +
  scale_fill_manual(values = COL_SYS, guide = "none") +
  scale_alpha_manual(values = c("switching genes" = 1, "called in BOTH" = 0.32), name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.2))) +
  labs(subtitle = "A naive overlap looks like failure",
       x = NULL, y = "genes") +
  th + theme(legend.position = "top")

# --- B: the evidence ladder ----------------------------------------------------------------
b_df <- dual |> count(.data$discovery_pretty, .data$tier, .data$direction_supported,
                      name = "genes")
p_b <- ggplot(b_df, aes(x = .data$genes, y = .data$tier, fill = .data$discovery_pretty)) +
  geom_col(width = 0.68, position = position_dodge2(padding = 0.12)) +
  # Bracket the three direction-supported tiers. NOTE the y positions: scale_y_discrete
  # (limits = rev) reverses the factor, so with six tiers "5" sits at position 1 and "0" at
  # position 6. Tiers 3/4/5 are therefore positions 1-3, NOT 4-6 -- shading 2.5-5.5 covers
  # the unsupported tiers instead, which is the opposite of the point.
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0.5, ymax = 3.5,
           fill = "#1F6F54", alpha = 0.07) +
  annotate("text", x = Inf, y = 3.5, hjust = 1.03, vjust = -0.5, size = 2.9,
           colour = "#1F6F54", fontface = "bold",
           label = "direction-supported: all 3 replicates agree") +
  scale_fill_manual(values = COL_SYS, name = "discovered in") +
  scale_y_discrete(limits = rev) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(subtitle = "What the other system actually shows for each hit",
       x = "genes", y = NULL) +
  th + theme(legend.position = "top")

# --- C: is the concordance more than chance? ----------------------------------------------
c_df <- bind_rows(
  conc |> transmute(discovery = PRETTY[.data$discovery],
                    metric = "same direction", observed = .data$pct_same_direction,
                    null = 50, p = .data$p_direction_vs_half),
  conc |> transmute(discovery = PRETTY[.data$discovery],
                    metric = "all 3 replicates agree", observed = .data$pct_all_three_agree,
                    null = .data$background_pct, p = .data$p_vs_background)
) |>
  mutate(metric = factor(.data$metric,
                         levels = c("same direction", "all 3 replicates agree")),
         star = ifelse(.data$p < 0.05, "*", "n.s."))

p_c <- ggplot(c_df, aes(x = .data$discovery, y = .data$observed, fill = .data$discovery)) +
  geom_col(width = 0.6) +
  geom_errorbar(aes(ymin = .data$null, ymax = .data$null), width = 0.6,
                linetype = "22", colour = "grey20") +
  geom_text(aes(label = paste0(round(.data$observed), "% ", .data$star)),
            vjust = -0.4, size = 3) +
  facet_wrap(~metric, nrow = 1) +
  scale_fill_manual(values = COL_SYS, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.22))) +
  labs(subtitle = "Concordance against its null (dashed = chance / genome-wide background)",
       x = NULL, y = "% of hits") +
  th + theme(axis.text.x = element_text(angle = 12, hjust = 1))

# --- D: the locus-fusion screen -----------------------------------------------------------
# Two isoforms of one gene always share exonic sequence. Where the candidate and the isoform
# moving hardest against it share NONE, the dIF is the relative expression of two neighbouring
# genes the assembler merged into one XLOC -- not an isoform switch. 19 of 63 fail.
scr <- scr |> mutate(screen = ifelse(.data$passes_screen, "genuine switch",
                                     "fused locus (0 shared exonic bp)"))
n_pass <- sum(scr$passes_screen); n_fail <- sum(!scr$passes_screen)

p_d <- ggplot(scr, aes(x = .data$set, fill = .data$screen)) +
  geom_bar(width = 0.66) +
  scale_fill_manual(values = c("genuine switch" = "#1F6F54",
                               "fused locus (0 shared exonic bp)" = "#C1443C"),
                    name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(subtitle = paste0("Locus screen: ", n_fail, " of ", nrow(scr), " are fused loci"),
       x = NULL, y = "genes") +
  th + theme(legend.position = "top", legend.text = element_text(size = 7.6),
             axis.text.x = element_text(angle = 15, hjust = 1)) +
  guides(fill = guide_legend(ncol = 1))

# --- E: the named shortlist, survivors only -------------------------------------------------
named <- scr |>
  filter(.data$passes_screen, is_real_symbol(.data$gene_name)) |>
  arrange(dplyr::desc(abs(.data$own_dIF))) |>
  slice_head(n = 22) |>
  mutate(label = stats::reorder(.data$gene_name, abs(.data$own_dIF)))

seg <- named |>
  transmute(.data$label, .data$set, x = abs(.data$own_dIF), xend = abs(.data$other_dIF))
pts <- bind_rows(
  named |> transmute(.data$label, .data$set, dIF = abs(.data$own_dIF),
                     which = "discovery system"),
  named |> transmute(.data$label, .data$set, dIF = abs(.data$other_dIF),
                     which = "other system")
)

p_e <- ggplot() +
  geom_segment(data = seg, aes(x = .data$x, xend = .data$xend, y = .data$label,
                               yend = .data$label, colour = .data$set),
               linewidth = 0.9, alpha = 0.55) +
  geom_point(data = pts, aes(x = .data$dIF, y = .data$label, colour = .data$set,
                             shape = .data$which), size = 2.4) +
  scale_colour_manual(values = COL_SET, name = "set") +
  scale_shape_manual(values = c("discovery system" = 16, "other system" = 1), name = NULL) +
  scale_x_continuous(expand = expansion(mult = c(0.02, 0.06))) +
  labs(subtitle = paste0("Follow-up shortlist: ", n_pass,
                         " genes surviving both screens, top 22 named"),
       x = "|dIF| (filled = discovery system, open = other system)", y = NULL) +
  th + theme(legend.position = "right")

# --- assemble --------------------------------------------------------------------------------
save_dir <- getOption("convergence.save_dir", NULL)
if (!requireNamespace("patchwork", quietly = TRUE)) {
  for (p in list(p_a, p_b, p_c, p_d, p_e)) print(p)
} else {
  g <- patchwork::wrap_plots(
    list(p_a, p_b, p_c, p_d, p_e),
    # Every patch's OWN area must be a rectangle, not merely the design string. Five panels:
    # A naive overlap, B evidence ladder, C concordance vs null, D locus screen, E named list.
    design = "AABBBB\nCCBBBB\nDDEEEE\nDDEEEE",
    heights = c(1, 1, 1.05, 1.05)
  ) +
    patchwork::plot_annotation(
      title = "From per-system hits to a named follow-up shortlist",
      subtitle = paste0(
        "Both systems refit from raw counts with the replicate blocked (script 23), so H and ",
        "T are measured the same way.\nNeither is treated as the anchor: a hit discovered in ",
        "either is graded by what the other shows on the same isoform, then screened for ",
        "fused loci."),
      tag_levels = "A",
      theme = theme(plot.title = element_text(face = "bold", size = 14),
                    plot.subtitle = element_text(size = 9.5, colour = "grey30"))
    ) &
    theme(plot.tag = element_text(face = "bold", size = 13))
  if (is.null(save_dir)) {
    print(g); message("Printed. Set options(convergence.save_dir=) to write the PNG.")
  } else {
    dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)
    ggsave(file.path(save_dir, "convergence_to_shortlist.png"), plot = g,
           width = 14.0, height = 15.5, dpi = 200, bg = "white")
    message("Wrote ", file.path(save_dir, "convergence_to_shortlist.png"))
  }
}

# --- the shortlist as a table ------------------------------------------------------------
cat("\n== Follow-up shortlist after BOTH screens ==\n")
cat("  direction-supported: ", nrow(scr), "\n", sep = "")
cat("  fused loci removed:  ", n_fail, "\n", sep = "")
cat("  surviving:           ", n_pass, "\n\n", sep = "")
print(as.data.frame(scr |> count(.data$set, .data$screen, name = "genes")), row.names = FALSE)

cat("\nTop 22 survivors with resolved symbols:\n")
print(as.data.frame(
  named |> transmute(gene = .data$gene_name,
                     discovered_in = .data$discovered_in, set = .data$set,
                     tier = substr(.data$tier, 1, 32),
                     own_dIF = round(.data$own_dIF, 3),
                     other_dIF = round(.data$other_dIF, 3),
                     shared_exonic_bp = .data$shared_exonic_bp)), row.names = FALSE)

cat("\nRemoved as fused loci: ",
    paste(sort(scr$gene_name[!scr$passes_screen]), collapse = ", "), "\n", sep = "")

if (!is.null(save_dir)) {
  utils::write.csv(
    scr |> arrange(.data$passes_screen, dplyr::desc(abs(.data$own_dIF))),
    file.path(TBL, "followup_shortlist_final.csv"), row.names = FALSE,
    fileEncoding = "UTF-8", na = "")
  message("Wrote ", file.path(TBL, "followup_shortlist_final.csv"))
}
