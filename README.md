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
| `data/raw/` | Your `T_isa_list.rds`, `U_isa_list.rds`, optional `H_isa_list.rds` (as named in config) |
| `data/processed/` | Isoform tibbles from `01_load_data.R` |
| `results/tables/` | Gene-level and downstream tables (e.g. from `02` onward) |
| `results/figures/` | Plots (to be written by `06` and report chunks) |
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

3. **Install R packages** (minimum for scripts 01–02):

   ```r
   install.packages(c("yaml", "tibble", "dplyr"))
   ```

   [Bioconductor](https://bioconductor.org/) packages (e.g. `IsoformSwitchAnalyzeR`, `clusterProfiler`) are only needed for upstream object creation or for planned scripts 05+.

4. **Run the pipeline** (from project root):

   ```bash
   Rscript scripts/01_load_data.R
   Rscript scripts/02_gene_level_summary.R
   ```

   If you run from `scripts/`, the bootstrap still finds the project as long as `config/config.yml` is discoverable (see `utils/bootstrap.R`).

## What each step does (implemented)

- **`01_load_data.R`**: `readRDS()` each configured list → `extract_isoform_features()` → `tag_novel_isoforms()` (default: `oId` starting with `PB`) → writes `data/processed/isoformFeatures_<label>.{rds,csv}`.
- **`02_gene_level_summary.R`**: Reads all `data/processed/isoformFeatures_*.rds`, runs `summarize_genes_from_isoform_table()` with thresholds from `config/config.yml` → `results/tables/gene_level_summary_<label>.{rds,csv}`.

**Gene table columns** (typical): `gene_id`, `gene_name`, `n_isoforms`, `n_switching_isoforms`, `min_isoform_switch_q`, `min_gene_switch_q`, `max_abs_dif`, `novel_involved`.

*Switching* is defined in helpers as: isoform q and gene q (if present) below the configured FDR, and an optional `min_abs_dif` on |dIF| when `dIF` is available.

## Novel isoforms (PacBio)

- Default rule: `oId` (or `config$novel$id_column`) with prefix `PB` (see `config$novel$pb_prefix`) → `is_novel_pacbio` in isoform tables, then aggregated to `novel_involved` at gene level.

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

- `03_novel_isoform_analysis.R` — distribution / contrast plots for novel genes.
- `04_comparison_T_vs_U.R` — shared/unique switching genes, Δ(ΔIF) style metrics.
- `05_pathway_enrichment.R` — `clusterProfiler` and pathway plots.
- `06_visualization.R` — publication figures (volcano, per-gene isoform panels).

## License

Replace `LICENSE` with the license your group uses for shared work.

## Contact

This scaffold is designed for a shared computational biology project; document ownership and data policies in this README or a separate `CONTRIBUTING.md` as your team requires.
