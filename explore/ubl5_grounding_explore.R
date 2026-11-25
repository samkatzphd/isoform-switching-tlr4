#!/usr/bin/env Rscript
# Standalone explorer for the UBL5 (U vs T) grounding analyses.
#
#   Rscript explore/ubl5_grounding_explore.R
#
# Reads ONLY CSVs from results/tables/. Writes nothing. Every object it builds is
# left in the global environment so you can poke at it interactively:
#   source("explore/ubl5_grounding_explore.R")
#
# Covers five grounded conclusions, each with the table that backs it:
#   A  UT reference is clean spliced cDNA      ut_reference_provenance.csv
#   B  UBL5 knockdown is ~30% at mRNA level    ubl5_knockdown_verification_libraries.csv
#   C  no shRNA escape isoform                 ubl5_shrna_target_sites.csv
#   D  KD carries NO elevated ISG tone         ut_isg_baseline_tone_camera.csv
#   E  genotype-level composition difference   ut_library_composition.csv

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(ggplot2)
})

ROOT <- Sys.getenv("ISOFORM_PROJECT_ROOT", unset = ".")
tbl  <- function(f) file.path(ROOT, "results", "tables", f)
theme_set(theme_bw(base_size = 11))

# ---------------------------------------------------------------- A. reference
prov <- read_csv(tbl("ut_reference_provenance.csv"), show_col_types = FALSE)

cat("\n=== A. UT reference provenance ===\n")
cat("All 315,349 FASTA records matched their GFF3 exon sums exactly.\n")
cat("Unlike the HT FASTA (12.2% unspliced genomic intervals), UT needs no rebuild.\n\n")
print(prov |> group_by(discovery_source) |>
        summarise(n = sum(n_transcripts), .groups = "drop") |>
        mutate(pct = round(100 * n / sum(n), 1)))
cat("\nPlatform class (structural novelty vs the reference annotation):\n")
print(prov |> group_by(platform_class) |>
        summarise(n = sum(n_transcripts), .groups = "drop") |>
        mutate(pct = round(100 * n / sum(n), 1)))

# The two axes are NOT the same question. 'discovery_source' says which assembly
# emitted the model; 'platform_class' says whether its structure is new. 8,366
# models are PB-emitted but class_code "=", i.e. PacBio-discovered yet
# Illumina-findable. Collapsing these two axes is the REVIEW_CHANGES 0j trap.
p_prov <- ggplot(prov, aes(reorder(class_code, -n_transcripts), n_transcripts,
                           fill = discovery_source)) +
  geom_col() + scale_y_log10() +
  labs(x = "gffcompare class_code", y = "transcripts (log10)",
       fill = "discovery source",
       title = "UT models by structural class and emitting assembly",
       subtitle = "class_code '=' means the intron chain already exists in the reference")

# ------------------------------------------------------------- B. knockdown
kd <- read_csv(tbl("ubl5_knockdown_verification_libraries.csv"), show_col_types = FALSE)

kd_summary <- kd |>
  group_by(unit, lps, genotype) |>
  summarise(mean = mean(value), sd = sd(value), n = n(), .groups = "drop") |>
  pivot_wider(names_from = genotype, values_from = c(mean, sd, n)) |>
  mutate(residual_fraction = round(mean_U_KD / mean_T_WT, 3),
         pct_knockdown     = round(100 * (1 - mean_U_KD / mean_T_WT), 1))

cat("\n=== B. UBL5 knockdown magnitude in the RNA-seq ===\n")
cat("CPM is primary; TPM is compositional and the arms differ in composition (see E).\n\n")
print(kd_summary |> select(unit, lps, mean_T_WT, mean_U_KD, residual_fraction, pct_knockdown))

kd_tests <- kd |> filter(unit == "CPM") |> group_by(lps) |>
  summarise(p_welch = t.test(value ~ genotype)$p.value, .groups = "drop")
cat("\nWelch t-test on CPM, n=3 per arm:\n"); print(kd_tests)
cat("\nINTERPRETATION: ~30% mRNA knockdown means ~70% residual UBL5 transcript.\n")
cat("A null result in this system is therefore weakly informative. Protein was not\n")
cat("measured; miR-30-based shRNAs can repress translation beyond mRNA loss, so the\n")
cat("protein-level knockdown is UNKNOWN and may be larger. Do not assume either way.\n")

p_kd <- ggplot(kd |> filter(unit == "CPM"),
               aes(interaction(genotype, lps), value, colour = genotype)) +
  geom_point(size = 2.5, position = position_jitter(width = .08, height = 0)) +
  stat_summary(fun = mean, geom = "crossbar", width = .4, colour = "grey30") +
  labs(x = NULL, y = "UBL5 CPM", title = "UBL5 knockdown, per library",
       subtitle = "each point is one sequencing library (n=3 per arm)")

# ------------------------------------------------------------- C. target site
site <- read_csv(tbl("ubl5_shrna_target_sites.csv"), show_col_types = FALSE)
iso  <- read_csv(tbl("ubl5_knockdown_verification_isoforms.csv"), show_col_types = FALSE)

cat("\n=== C. shRNA target-site coverage ===\n")
cat("Probe searched against all 315,349 UT transcripts. The full 22-mer\n")
cat("CGGGATGAACCTGGAGCTTTAT has ZERO matches; the leading CGG is a cloning\n")
cat("overhang, not target-complementary. The 19-mer GATGAACCTGGAGCTTTAT hits 14.\n\n")
print(site |> count(gene_name, name = "transcripts_with_site"))
cat("\nOnly off-target in the transcriptome is the pseudogene UBL5P2 (<1 TPM).\n")

cat("\nUBL5 output carried by site-bearing vs site-free models:\n")
print(iso |> group_by(lps, has_shrna_site) |>
        summarise(sum_T_WT = round(sum(T_WT), 1), sum_U_KD = round(sum(U_KD), 1),
                  n_models = n(), .groups = "drop"))
cat("\nThe single site-free model (TCONS_00145328) carries ~0.2% of UBL5 output,\n")
cat("so the knockdown ceiling is NOT explained by isoform escape.\n")

# -------------------------------------------------------------------- D. ISGs
cam <- read_csv(tbl("ut_isg_baseline_tone_camera.csv"), show_col_types = FALSE)
q2  <- read_csv(tbl("ut_q2_genelevel_KD_vs_WT.csv.gz"), show_col_types = FALSE)

cat("\n=== D. Is there chronic interferon tone in the transduced KD line? ===\n")
cat("camera() corrects for inter-gene correlation; ISGs are heavily co-regulated\n")
cat("so an uncorrected rank test here would be badly anti-conservative.\n\n")
print(cam |> select(set, contrast, NGenes, Direction, PValue, FDR) |> arrange(contrast, PValue))
cat("\nDirection is DOWN, not up. The feared artefact -- retroviral transduction\n")
cat("inducing ISGs and making the KD look inflamed -- is ABSENT. That removes the\n")
cat("specific worry, but note this is a COMPETITIVE test run against a background\n")
cat("in which ~56% of genes are called DE, so read it as a relative shift.\n")

p_isg <- ggplot(q2 |> mutate(cls = ifelse(is_isg, "ISG", "other")),
                aes(logFC, colour = cls)) +
  stat_ecdf(linewidth = .8) + facet_wrap(~contrast) + xlim(-3, 3) +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey50") +
  labs(x = "log2FC (KD - WT)", y = "cumulative fraction", colour = NULL,
       title = "ISG distribution shift, KD vs WT",
       subtitle = "ISG curve left of 'other' = ISGs lower in the knockdown")

# ------------------------------------------------------------- E. composition
comp <- read_csv(tbl("ut_library_composition.csv"), show_col_types = FALSE)

cat("\n=== E. Library composition differs systematically by genotype ===\n")
print(comp |> arrange(group, sample))
cat("\nTMM norm factors separate PERFECTLY by genotype: all six T libraries <= 1.013,\n")
cat("all six U libraries >= 1.059, no overlap. The WT libraries are more concentrated\n")
cat("(top-50 share ~18-20% vs ~15-16%), and this is present at BASELINE, so it is not\n")
cat("merely an LPS effect.\n\n")
cat("CONSEQUENCE: TMM's 'most genes unchanged' assumption is strained here, so every\n")
cat("ABSOLUTE abundance claim (Q2) rests on a normalisation the data contest. Measures\n")
cat("of composition WITHIN a gene (dIF, hence Q1) are invariant to library scaling and\n")
cat("are not exposed to this. Prefer Q1 phrasing wherever both are available.\n")

p_comp <- ggplot(comp, aes(nf_TMM, top50, colour = group, label = sample)) +
  geom_point(size = 3) +
  labs(x = "TMM normalisation factor", y = "% of counts in top 50 genes",
       title = "Genotype separates on library composition",
       subtitle = "T libraries cluster low-nf/high-concentration; U libraries the reverse")

cat("\n--- plots built: p_prov, p_kd, p_isg, p_comp (print any of them) ---\n")
