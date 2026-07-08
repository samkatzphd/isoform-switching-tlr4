# Isoform switching — downstream analysis pipeline

Reproducible, modular R workflow for **isoform-level switching** results from [IsoformSwitchAnalyzeR](https://bioconductor.org/packages/IsoformSwitchAnalyzeR/) (or equivalent `SwitchList` / `SwitchAnalyzeRlist` objects). This repository starts **after** differential isoform analysis: you provide saved list objects, not raw RNA-seq.

## Design

- **Modular scripts** in `scripts/` (no monolithic notebook; notebooks/reports are optional for exploration).
- **Config-driven** paths and thresholds: `config/config.yml` (no hard-coded project paths in code).
- **Explicit outputs**: intermediate tables as `.rds` (round-trip R objects) and `.csv` (interchange, spreadsheets).
- **Helper functions** in `utils/helper_functions.R` (isoform table extraction, novel flags, gene-level collapse).

## Repository layout

| Path | Role |
|------|------|
| `data/raw/` | Placeholder for copies of inputs if you do not use absolute paths in `config` |
| `data/processed/` | Isoform tibbles from `01_load_data.R` |
| `results/tables/` | Gene tables, QC snapshots, novel summaries, HT comparisons (`02`+) |
| `results/figures/` | QC, novel, and HT switch plots |
| `reports/isoform_switching_overview.qmd` | QC overview report |
| `reports/ht_switch_analysis.qmd` | HT top-switch + T↔H overlap report |
| `reports/novel_isoform_analysis.qmd` | Novel PacBio isoform report (per-dataset + HT vs UT) |
| `_quarto.yml` | Shared Quarto defaults for reports |
| `config/config.yml` | Input paths, labels (T/U/H), q-value cutoffs, novel id column/prefix |
| `utils/helper_functions.R` + `utils/bootstrap.R` | Shared I/O, extraction, gene summary |
| `environment.yml` | Optional Conda stack; or use R + `renv` (see below) |
| `reports/exploratory.qmd` | Optional Quarto scratchpad for QC and narrative |

## Quick start

1. **Clone** and open the project in RStudio (or any editor) with the **working directory** at the project root, or set:

   ```bash
   export ISOFORM_PROJECT_ROOT="/absolute/path/to/isoform-switch-pipeline"
   ```

2. **Set dataset inputs** in `config/config.yml` under `datasets:` (e.g. `T_HT`, `H_HT`, `U_UT`, `T_UT`).  
   Inputs can be `.rds` or `.RData/.Rdata`. For `.RData` with multiple objects, set `object_name` per dataset.  
   Use the final reference annotation `.gff3` per dataset (`annotation_path`) to map `isoform_id -> oId/cmp_ref/class_code`.

3. **Install R packages** (minimum for scripts 01–02; add **ggplot2** and **scales** for `02c` figures):

   ```r
   install.packages(c("yaml", "tibble", "dplyr", "ggplot2", "scales"))
   ```

   [Bioconductor](https://bioconductor.org/) packages (e.g. `IsoformSwitchAnalyzeR`, `clusterProfiler`) are only needed for upstream object creation or for planned scripts 05+.

4. **Run the pipeline** (from project root):

   ```bash
   Rscript scripts/01_load_data.R
   Rscript scripts/02_gene_level_summary.R
   Rscript scripts/02b_qc_snapshot.R
   Rscript scripts/02c_qc_figures.R
   Rscript scripts/03_novel_isoform_analysis.R
   ```

   If you run from `scripts/`, the bootstrap still finds the project as long as `config/config.yml` is discoverable (see `utils/bootstrap.R`).

4. **HTML reports (optional, requires [Quarto](https://quarto.org/docs/get-started/))**:

   ```bash
   quarto render reports/isoform_switching_overview.qmd
   quarto render reports/ht_switch_analysis.qmd
   quarto render reports/novel_isoform_analysis.qmd
   ```

   Outputs (with `embed-resources: true`):
   - `reports/isoform_switching_overview.html`
   - `reports/ht_switch_analysis.html`
   - `reports/novel_isoform_analysis.html`

## What each step does (implemented)

- **`01_load_data.R`**: `readRDS()` each configured list → `extract_isoform_features()` → `tag_novel_isoforms()` (default: `oId` starting with `PB`) → writes `data/processed/isoformFeatures_<label>.{rds,csv}`.
- **`02_gene_level_summary.R`**: Reads all `data/processed/isoformFeatures_*.rds`, runs `summarize_genes_from_isoform_table()` with thresholds from `config/config.yml` → `results/tables/gene_level_summary_<label>.{rds,csv}`.
- **`02b_qc_snapshot.R`**: One-row-per-dataset QC and `results/tables/qc_top10_genes_per_dataset.csv`.
- **`02c_qc_figures.R`**: Writes PNGs under `results/figures/` (gene + isoform QC plots).
- **`03_novel_isoform_analysis.R`**: Novel vs known contrasts (effect sizes, class codes, top novel switching genes/isoforms), plus **HT vs UT** novelty rates and shared novel-involved gene-symbol overlap → `results/tables/novel_*` and `results/figures/novel/`.
- **`06_visualization.R`**: HT-focused top-30 switch ranks, ISA `switchPlot`s, and T↔H overlap tables/report.

**Gene table columns** (typical): `gene_id`, `gene_name`, `n_isoforms`, `n_switching_isoforms`, `min_isoform_switch_q`, `min_gene_switch_q`, `max_abs_dif`, `novel_involved`.

*Switching* is defined in helpers as: isoform q and gene q (if present) below the configured FDR, and an optional `min_abs_dif` on |dIF| when `dIF` is available.

## Novel isoforms (PacBio)

- Default rule: `oId` (or `config$novel$id_column`) with prefix `PB` (see `config$novel$pb_prefix`) → `is_novel_pacbio` in isoform tables, then aggregated to `novel_involved` at gene level.
- Script `03` ranks novel switching isoforms/genes, summarizes annotation `class_code`, and contrasts HT vs UT novelty (rates + gene-symbol overlap).
- Report: `reports/novel_isoform_analysis.qmd` (render after running `03`).

## Reproducibility: `renv` (optional)

To pin package versions, from R in the project root:

```r
if (!requireNamespace("renv", quietly = TRUE)) install.packages("renv")
renv::init()  # captures DESCRIPTION / lockfile; commit renv.lock when satisfied
```

To restore on another machine: `renv::restore()`.

A starter `renv.lock` is **not** committed here so your platform/R version are not over-constrained; generate it in your environment.

## Python (optional)

`environment.yml` can install R and Python. Use Python only for extras (e.g. jupyter, custom plotting); the core pipeline is R.

## GitHub

Initialize remote and push (replace URL):

```bash
cd isoform-switch-pipeline
git remote add origin https://github.com/<org>/isoform-switch-pipeline.git
git push -u origin main
```

## What is scaffold only (for your feedback, then we extend)

- `04_comparison_T_vs_U.R` — shared/unique switching genes, Δ(ΔIF) style metrics.
- `05_pathway_enrichment.R` — `clusterProfiler` and pathway plots.

## License

Replace `LICENSE` with the license your group uses for shared work.

## Contact

This scaffold is designed for a shared computational biology project; document ownership and data policies in this README or a separate `CONTRIBUTING.md` as your team requires.
