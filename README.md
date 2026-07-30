# Isoform switching — downstream analysis pipeline

Reproducible, modular R workflow for **isoform-level switching** results from [IsoformSwitchAnalyzeR](https://bioconductor.org/packages/IsoformSwitchAnalyzeR/) (or equivalent `SwitchList` / `SwitchAnalyzeRlist` objects). This repository starts **after** differential isoform analysis: you provide saved list objects, not raw RNA-seq.

**Remote:** [samkatzphd/isoform-switching-tlr4](https://github.com/samkatzphd/isoform-switching-tlr4)

**Local path:** `/Users/samkatz/projects/isoform-switch-pipeline`

## Biological framing

Two reference transcriptomes and four IsoformSwitchAnalyzeR (ISA) objects:

| Dataset | Biological role | Reference |
|---------|-----------------|-----------|
| `T_HT` | WT THP1 ± LPS (primary HT analysis) | HT |
| `H_HT` | H genotype ± LPS (supplementary HT) | HT |
| `T_UT` | WT THP1 ± LPS on UT annotation | UT |
| `U_UT` | UBL5 knockout ± LPS | UT |

Configured analysis questions (`config/config.yml` → `analysis_questions`):

1. **Q1 (HT):** WT pre/post LPS isoform differences anchored on T; H is supplementary.
2. **Q2 (UT):** Role of UBL5 in post-LPS isoform selection — compare T vs U on the UT reference.

## Design

- **Modular scripts** in `scripts/` (no monolithic notebook; Quarto reports are for exploration and sharing).
- **Config-driven** paths and thresholds: `config/config.yml` (no hard-coded project paths in code).
- **Explicit outputs**: intermediate tables as `.rds` (round-trip R objects) and `.csv` (spreadsheets).
- **Helper functions** in `utils/helper_functions.R` (isoform extraction, novel flags, gene-level collapse).

## Repository layout

| Path | Role |
|------|------|
| `config/config.yml` | Dataset paths, significance cutoffs, novel-id rule |
| `data/processed/` | Isoform tibbles from step `01` |
| `results/tables/` | Gene summaries, QC, novel, HT, and UT comparison tables |
| `results/figures/` | QC, novel, HT switchPlots, UT T-vs-U figures |
| `scripts/` | Numbered analysis steps (see below) |
| `utils/` | Bootstrap + shared helpers |
| `reports/*.qmd` | Quarto sources; HTML companions committed where rendered |
| `_quarto.yml` | Shared Quarto defaults (`embed-resources: true`) |
| `environment.yml` | Optional Conda stack |

## Significance definition

An isoform (and its gene) is treated as *switching* when:

- isoform q-value &lt; `significance.isoform_q` (default **0.05**),
- gene q-value (if present) &lt; `significance.gene_q` (default **0.05**),
- and **`|dIF| >= significance.min_abs_dif`** (currently **0.15**).

The `|dIF|` floor was chosen after a sensitivity scan (`00`): 0.15–0.20 keeps larger, more functionally plausible fraction changes than an unfiltered q-only call.

## Pipeline steps (what each script does)

Run from the **project root** unless noted. Order matters for dependents of `01`/`02`.

### `00_threshold_scan_abs_dif.R` — optional sensitivity helper

- **Purpose:** After `01`, explore how switching gene/isoform counts change as `min |dIF|` increases (q-filters held fixed).
- **Does not** rewrite `config.yml`; used to justify the current `min_abs_dif: 0.15`.
- **Outputs:** Console summaries + scan figures/tables under `results/` (when run).

### `01_load_data.R` — load ISA objects and annotate isoforms

- Load each configured ISA object (`.rds` / `.RData`).
- Extract isoform-level features (`extract_isoform_features()`).
- Join final reference **GFF3** (`annotation_path`) for `oId` / `cmp_ref` / `class_code` and gene symbols.
- Tag PacBio novel isoforms (`oId` prefix `PB` → `is_novel_pacbio`).
- **Writes:** `data/processed/isoformFeatures_<dataset>.{rds,csv}` for `T_HT`, `H_HT`, `T_UT`, `U_UT`.

### `02_gene_level_summary.R` — collapse isoforms to genes

- Read all processed isoform tables.
- Apply config q and `|dIF|` cutoffs via `summarize_genes_from_isoform_table()`.
- **Writes:** `results/tables/gene_level_summary_<dataset>.{rds,csv}`.
- Typical columns: `gene_id`, `gene_name`, `n_isoforms`, `n_switching_isoforms`, q / `|dIF|` summaries, `novel_involved`.

### `02b_qc_snapshot.R` — tabular QC

- One-row-per-dataset QC counts and top genes.
- **Writes:** QC snapshot tables including `results/tables/qc_top10_genes_per_dataset.csv`.

### `02c_qc_figures.R` — QC plots

- Gene- and isoform-level QC graphics.
- **Writes:** PNGs under `results/figures/`.

### `03_novel_isoform_analysis.R` — PacBio novelty + HT vs UT contrast

- Contrast novel vs known isoforms (effect sizes, SQANTI-like `class_code`, top novel switching genes/isoforms).
- HT vs UT novelty rates and shared novel-involved **gene-symbol** overlap.
- **Writes:** `results/tables/novel_*`, `results/figures/novel/`.
- **Report:** `reports/novel_isoform_analysis.{qmd,html}`.

### `04_comparison_T_vs_U.R` — UT WT (T) vs knockout (U)

- Shared / T-only / U-only switching genes (by clean gene symbol) and isoforms (by `isoform_id`).
- Concordance and ΔdIF-style stats on genes present in both analyses.
- Overlap figures, per-gene explorer panels, ISA `switchPlot`s where available.
- **Writes:** `results/tables/ut_*`, `results/figures/ut_t_vs_u/`.
- **Report:** `reports/ut_t_vs_u_comparison.{qmd,html}`.

### `05_pathway_enrichment.R` — scaffold only

- Placeholder for `clusterProfiler` (GO / KEGG / Reactome).
- Currently stops with “Not yet implemented.”

### `06_visualization.R` — HT top switches + T↔H overlap

- Focus on `T_HT` and `H_HT`.
- Rank genes by `max_abs_dif` among switching isoforms; flag genes with ≥2 isoforms in opposing directions.
- Top-30 tables; annotate ISA with GFF-derived gene symbols and `(PB)` on novel isoforms; draw IsoformSwitchAnalyzeR `switchPlot`s.
- T↔H overlap: isoforms significant in T shown in H even when H is non-significant.
- **Writes:** HT rank/overlap tables under `results/tables/`, plots under `results/figures/ht_top_switch_plots/`.
- **Report:** `reports/ht_switch_analysis.{qmd,html}`.

## Quick start

1. **Clone / open** the project at `/Users/samkatz/projects/isoform-switch-pipeline` and use that as the working directory, or set:

   ```bash
   export ISOFORM_PROJECT_ROOT="/Users/samkatz/projects/isoform-switch-pipeline"
   ```

2. **Edit** `config/config.yml` — `datasets:` paths (`isa_path`, `annotation_path`), and optionally `significance.min_abs_dif`.

3. **Install packages** (minimum for `01`–`04`; add ggplot2/scales for figures; IsoformSwitchAnalyzeR for `06` switchPlots):

   ```r
   install.packages(c("yaml", "tibble", "dplyr", "ggplot2", "scales"))
   # BiocManager::install("IsoformSwitchAnalyzeR")  # for switchPlot in 06 / 04
   ```

4. **Run the core pipeline:**

   ```bash
   Rscript scripts/01_load_data.R
   Rscript scripts/02_gene_level_summary.R
   Rscript scripts/02b_qc_snapshot.R
   Rscript scripts/02c_qc_figures.R
   Rscript scripts/03_novel_isoform_analysis.R
   Rscript scripts/04_comparison_T_vs_U.R
   Rscript scripts/06_visualization.R
   # Optional: Rscript scripts/00_threshold_scan_abs_dif.R
   ```

5. **Render HTML reports** ([Quarto](https://quarto.org/docs/get-started/)):

   ```bash
   quarto render reports/isoform_switching_overview.qmd
   quarto render reports/ht_switch_analysis.qmd
   quarto render reports/novel_isoform_analysis.qmd
   quarto render reports/ut_t_vs_u_comparison.qmd
   ```

   Committed HTML (self-contained via `embed-resources`):

   - `reports/isoform_switching_overview.html`
   - `reports/ht_switch_analysis.html`
   - `reports/novel_isoform_analysis.html`
   - `reports/ut_t_vs_u_comparison.html`

## Novel isoforms (PacBio)

- Default rule: `oId` (or `config$novel$id_column`) with prefix `PB` → `is_novel_pacbio`, rolled up to `novel_involved` at gene level.
- Step `03` ranks novel switching events and contrasts HT vs UT.
- Report: `reports/novel_isoform_analysis.qmd`.

## Reproducibility: `renv` (optional)

```r
if (!requireNamespace("renv", quietly = TRUE)) install.packages("renv")
renv::init()
# commit renv.lock when satisfied; restore elsewhere with renv::restore()
```

A starter `renv.lock` is not committed so platform/R versions are not over-constrained.

## License / contact

Replace `LICENSE` with your group’s license. Document ownership and data policies here or in `CONTRIBUTING.md` as needed.
