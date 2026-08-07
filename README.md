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
| `config/config.yml` | Dataset paths, significance cutoffs, `analysis:` knobs, output policy |
| `data/processed/` | Isoform tibbles from step `01` (`.rds`; CSV twins suppressed) |
| `results/tables/` | Gene summaries, QC, novel, HT, and UT comparison tables |
| `results/figures/` | QC, novel, HT switchPlots, UT T-vs-U figures |
| `results/run_manifest.json` | Per-script provenance: git SHA, R/package versions, thresholds |
| `scripts/` | Numbered analysis steps (see below) |
| `utils/` | Bootstrap + shared helpers (single switching rule, symbol rule, writers) |
| `reports/_setup.R` | Shared report setup (root discovery, path helpers, caveat callout) |
| `reports/*.qmd` | Quarto sources; HTML companions committed where rendered |
| `docs/REVIEW_CHANGES.md` | Review of the original code and every correction made |
| `docs/AGENT_CONTEXT.md` | Source for the Agent Context system-prompt field |
| `docs/STATISTICAL_METHODS.md` | Every method used, why it was chosen, and how it can mislead |
| `docs/OPEN_QUESTIONS.md` | Tests still to run, ranked by what would change a conclusion |
| `docs/external_review/` | Second external review with its own figures, tables and gene sets |
| `_quarto.yml` | Shared Quarto defaults (`embed-resources: true`) |
| `environment.yml` | Optional Conda stack |

## ⚠️ Input caveat: the ISA objects are pre-reduced

All four saved `switchAnalyzeRlist` objects contain **only genes that already pass the
gene-level q cutoff** — they were written after
`isoformSwitchTestDEXSeq(reduceToSwitchingGenes = TRUE)`, the ISA default. Verified per
run and recorded in `results/tables/input_object_reduction_check.csv`.

This means gene counts are *genes retained*, not *genes tested*, and every
percentage-of-genes has a denominator already selected on the outcome.

**All four datasets now work around this via the context layer** (below), which recovers the
real tested denominator without re-running any analysis. See
[`docs/REVIEW_CHANGES.md`](docs/REVIEW_CHANGES.md) §0 and §0b.

## Unfiltered context layer

`isa_unfiltered_path` in `config/config.yml` points a dataset at an *unreduced* export of
the same analysis. Script `01` then writes a slim per-isoform context table plus
per-replicate isoform fractions to `data/processed/isoformContext*_<label>.rds`.

Configured for **all four datasets** — 10,907 / 8,160 / 11,582 / 11,969 genes tested for
`T_HT` / `H_HT` / `T_UT` / `U_UT`, versus 216 / 385 / 248 / 127 retained in the reduced
objects. `01` and `06` still handle a missing context explicitly, and `06` uses whichever
side is available rather than requiring both.

**Significance calls now come FROM the context layer**, via `load_scoring_table()`. They
previously came from the reduced objects while presence came from the context — two
different FDR universes, which misclassified genes (see `docs/REVIEW_CHANGES.md` §0c). What
the context enables:

- Gene explorer panels show **both genotypes for any gene**, whether or not it reached
  significance in each (101 of 102 T-only genes have a plottable U side; 42 of 44 the other way).
- Overlap classes distinguish *tested in the other dataset and not switching* (180
  isoforms) from *not detected there* (6).
- A co-occurrence test against a genuine background of 11,126 genes quantified in both
  (OR 36.1; the HT pair gives OR 5.4 over 7,671 genes).
- A real denominator for every rate, and the abundance floor below.

## Abundance floor for switching calls

`significance.min_gene_expression: 12`. dIF is a ratio, so at low gene abundance the isoform
fraction is estimated from few reads and swings between replicates.
`02d_expression_diagnostics.R` measures this: within-condition IF noise (condition fixed, so
pure noise) against gene expression. Below ~12 a threshold-sized shift sits within ~2 SD of
replicate noise, and switching calls are ~6x enriched there — a pattern that holds within
every isoform-count stratum, so it is abundance rather than transcript complexity.

The floor removes roughly a third to a half of prior calls, which is the intent. Excluded
isoforms are flagged (`low_expression`, `is_switching_before_expr_floor`) rather than
discarded from the tables, so the cost stays auditable.

## Getting input files onto the drive without corruption

Two ISA exports have arrived damaged, and neither was visible from a file listing: one was
truncated (130 MB of an expected 615 MB), the other was full size but corrupt mid-stream and
only failed 20 minutes into a pipeline run. Both are caught in seconds by the steps below.

> **The drive is exFAT.** `diskutil info /Volumes/Expansion` → *File System Personality:
> ExFAT*. exFAT has **no journaling**: if a write is interrupted by an unplug, a sleep, or a
> cable knock, the file is left partially written with no recovery and no error at the time.
> This is the most likely cause of the mid-stream corruption. Treat every write to this drive
> as unverified until checked, and always eject properly.

### The one command that gates a run

```bash
Rscript scripts/99_verify_inputs.R && Rscript scripts/01_load_data.R
```

`99_verify_inputs.R` reads `config/config.yml`, checks every `isa_path`,
`isa_unfiltered_path` and `annotation_path` for existence and compression integrity, writes
`results/tables/input_file_verification.csv`, and **exits non-zero** if anything fails — so
the `&&` stops the pipeline before it wastes time on a bad file. Add `--load` to also
deserialise each object (slow, but catches damage that survives the integrity check).

### Download and unzip checklist

**1. Verify the zip before trusting it.** The zip format stores a CRC32 for every member, so
this detects a truncated or corrupted download without extracting anything:

```bash
unzip -t archive.zip
```

Expect `No errors detected`. If it reports errors, the download is bad — re-download; do not
extract.

**2. Extract to the local disk, not directly onto the external drive.** Extracting straight
onto exFAT doubles the exposure to interrupted writes, and a failure midway leaves files that
look complete.

```bash
mkdir -p ~/isa_staging && cd ~/isa_staging
unzip /path/to/archive.zip
```

Prefer `unzip` over double-clicking in Finder — Archive Utility reports errors poorly and can
leave partial output. For very large archives, `ditto -x -k archive.zip .` is a good
alternative.

**3. Verify each extracted `.Rdata` independently.** R data files are gzip streams, so gzip
can check them end to end. This catches a bad extraction even if the zip was fine:

```bash
for f in **/*.Rdata; do
  printf '%s: ' "$f"; gzip -t "$f" 2>&1 && echo OK
done
```

**4. Copy to the drive with verification.** `rsync -c` compares checksums rather than size and
timestamp, so it catches a silently corrupted copy and re-sends only what is wrong:

```bash
rsync -avh --checksum --progress \
  ~/isa_staging/T_minus-T_plus_Unfiltered_2025-Oct-12 \
  /Volumes/Expansion/IsoformSwitchAnalyzer/H-T_Comparisons/
```

**5. Eject properly, then re-verify on the drive.** This is the step that matters most on
exFAT — the copy is not necessarily flushed to disk when `rsync` returns.

```bash
diskutil eject /Volumes/Expansion     # wait for it to disappear, then reconnect
Rscript scripts/99_verify_inputs.R    # confirms the files are intact where they now live
```

If you want belt and braces, compare hashes across the copy instead of trusting `rsync -c`:

```bash
shasum -a 256 ~/isa_staging/**/*.Rdata
shasum -a 256 /Volumes/Expansion/IsoformSwitchAnalyzer/H-T_Comparisons/**/*.Rdata
```

**Never copy an `.Rdata` into the repository.** Raw inputs belong on the drive; `.gitignore`
blocks `*.Rdata` and `*_Unfiltered_*/` because a truncated 130 MB partial once got swept in by
`git add -A` and was rejected by GitHub's 100 MB limit.

## Significance definition

An isoform (and its gene) is treated as *switching* when:

- isoform q-value &lt; `significance.isoform_q` (default **0.05**),
- gene q-value (when the column exists) &lt; `significance.gene_q` (default **0.05**),
- and **`|dIF| >= significance.min_abs_dif`** (currently **0.15**).

With `significance.require_finite_dif: true` (default) an isoform with a missing or
non-finite `dIF` or gene q **fails** the filter rather than passing on the q-value alone.

This rule has exactly one implementation — `score_isoforms()` in
`utils/helper_functions.R`. Every script calls it; none re-implements it.

The `|dIF|` floor was chosen after a sensitivity scan (`00`): 0.15–0.20 keeps larger, more functionally plausible fraction changes than an unfiltered q-only call. Note that the scan operates inside the already-significant set described above, so it shows how the floor *prunes*, not a power curve.

## Pipeline steps (what each script does)

Run from the **project root** unless noted. Order matters for dependents of `01`/`02`.

### `00_threshold_scan_abs_dif.R` — optional sensitivity helper

- **Purpose:** After `01`, explore how switching gene/isoform counts change as `min |dIF|` increases (q-filters held fixed).
- **Does not** rewrite `config.yml`; used to justify the current `min_abs_dif: 0.15`.
- **Outputs:** Console summaries + scan figures/tables under `results/` (when run).

### `01_load_data.R` — load ISA objects and annotate isoforms

- Load each configured ISA object (`.rds` / `.RData`), using the **pinned** `object_name`
  (each `.Rdata` holds 13 objects, so auto-picking was order-dependent).
- Extract isoform-level features (`extract_isoform_features()`).
- Join final reference **GFF3** (`annotation_path`) for `oId` / `cmp_ref` / `class_code` and gene symbols — each column joined only if actually missing.
- Tag PacBio novel isoforms (`oId` prefix `PB` → `is_novel_pacbio`).
- Check whether each object was pre-reduced to switching genes.
- **Writes:** `data/processed/isoformFeatures_<dataset>.rds`,
  `results/tables/input_object_reduction_check.{rds,csv}`, `results/run_manifest.json`.

### `02_gene_level_summary.R` — collapse isoforms to genes

- Read all processed isoform tables.
- Apply config q and `|dIF|` cutoffs via `summarize_genes_from_isoform_table()`.
- **One row per gene** — gene symbols are collapsed before grouping, and the script asserts
  `rows == distinct gene_id`.
- **Writes:** `results/tables/gene_level_summary_<dataset>.{rds,csv}`.
- Columns: `gene_id`, `gene_name`, `n_isoforms`, `n_switching_isoforms`, `n_novel_isoforms`,
  `min_isoform_switch_q`, `min_gene_switch_q`, **`max_abs_dif_all_isoforms`**,
  **`max_abs_dif_switching`**, `novel_involved`. The two `max_abs_dif_*` columns are
  different quantities and are named apart on purpose — ranking uses the switching one.

### `02b_qc_snapshot.R` — tabular QC

- One-row-per-dataset QC counts and top genes.
- **Writes:** QC snapshot tables including `results/tables/qc_top10_genes_per_dataset.csv`.

### `02c_qc_figures.R` — QC plots

- Gene- and isoform-level QC graphics.
- **Writes:** PNGs under `results/figures/`.

### `02d_expression_diagnostics.R` — abundance floor and between-group expression

- Relates within-condition IF noise to gene expression to derive the floor (see above).
- Compares baseline expression between the genotypes sharing a reference, to rule out a
  global power difference behind any "genotype X switches less" claim.
- **Writes:** `expression_if_stability`, `expression_floor_recommendation`,
  `expression_between_groups`, figures under `results/figures/expression/`.

### `02e_response_magnitude_control.R` — is a switching deficit just a smaller LPS response?

- Stratifies each gene by its own |log2FC| and asks whether a switching-rate difference
  survives. **Changes the Q2 conclusion:** the KO's LPS response runs at ~half amplitude
  (slope 0.52), and the crude switching OR of 0.51 becomes MH OR 0.77 (p = 0.11).
- Run on both pairs. The HT pair *strengthens* on adjustment (OR 2.29 → 2.49), which is the
  control showing stratification is not flattening everything.
- **Writes:** `ut_response_magnitude_{strata,control}.*`, figures under `results/figures/response/`.

### `02f_reference_concordance.R` — the same-sample technical ceiling

- `T_HT` and `T_UT` are **the same six libraries** under two references (verified from
  replicate sample names; gene-level r ≈ 0.99). All disagreement is annotation + thresholding.
- Switching-gene Jaccard **0.35** — so a single-dataset call has roughly a one-in-three chance
  of not replicating under re-annotation. Judge replication rates against this, not 100%.
- Note `TCONS_*` ids are **not** comparable across references; compare on gene symbols.
- **Writes:** `reference_{concordance,discordance_reasons,discordant_genes}_*.*`.

### `03_novel_isoform_analysis.R` — PacBio novelty + HT vs UT contrast

- Contrast novel vs known isoforms (effect sizes, SQANTI-like `class_code`, top novel switching genes/isoforms).
- HT vs UT novelty rates and shared novel-involved **gene-symbol** overlap.
- **Writes:** `results/tables/novel_*`, `results/figures/novel/`.
- **Report:** `reports/novel_isoform_analysis.{qmd,html}`.

### `04_comparison_T_vs_U.R` — UT WT (T) vs knockout (U)

- Shared / T-only / U-only switching genes (by clean gene symbol) and isoforms (by `isoform_id`).
- Retention accounting, ΔdIF effect sizes, and **one-way concordance** (condition on
  significance in one dataset, measure the effect in the other).
- **No enrichment test.** A Fisher test previously reported here had a background of 25
  already-significant genes; it was removed rather than caveated. Columns prefixed
  `selconf_` are conditioned on significance in both datasets and are descriptive only.
- Overlap figures, per-gene explorer panels, ISA `switchPlot`s where available.
- Skips plotting (rather than deleting existing figures) when the ISA objects are unreachable.
- **Writes:** `results/tables/ut_*`, `results/figures/ut_t_vs_u/`.
- **Report:** `reports/ut_t_vs_u_comparison.{qmd,html}`.

### `05_pathway_enrichment.R` — GO / Reactome / KEGG over switching gene sets

- **The universe is the point.** Enrichment runs against the genes *quantified and tested in
  the same experiment* (from the context layer), never the whole genome and never the
  pre-reduced object. A dataset without a context table is **skipped**, not silently given a
  whole-genome background.
- Gene sets: per-dataset switching genes, plus the UT overlap classes (shared / T-only /
  U-only) against genes tested in both UT arms.
- **Writes:** `results/tables/pathway_enrichment_{summary,results,top_terms,driver_genes,driver_summary}.*`,
  dotplots under `results/figures/pathway/`.
- The driver tables exist because a term count alone is misleading: GO and Reactome are
  nested hierarchies, so a few genes in one module generate dozens of overlapping terms.
  Result: 82 terms from the 15 shared UT genes, of which **95% of gene-term hits come from
  five chemokines/cytokines**; the two largest genotype-specific sets return nothing.

### `07_isg_analysis.R` — interferon-response gene profile

- Tests a **curated** interferon gene set once per dataset against the same tested
  background, rather than relying on open GO enrichment (which scatters interferon biology
  across dozens of nested terms).
- Default sets from GO via `org.Hs.eg.db` (offline): type I = GO:0034340 + GO:0060337,
  type II = GO:0034341. Swap in Schoggins/Interferome via `isg.custom_set_path`.
- **Result:** switching is ~6–7× over-represented among interferon-response genes in `T_HT`,
  `T_UT` and `U_UT` (all q < 0.01), driven by the type II set. `H_HT` shows no enrichment,
  though the difference from `T_HT` is not itself significant.
- **Writes:** `results/tables/isg_{enrichment,gene_profile,switching_genes,pair_profile,pair_summary,effect_size}.*`,
  figures under `results/figures/isg/`.

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
   install.packages(c("yaml", "tibble", "dplyr", "ggplot2", "scales", "jsonlite"))
   # BiocManager::install("IsoformSwitchAnalyzeR")  # for switchPlot in 06 / 04
   ```

   `jsonlite` is only needed for `results/run_manifest.json`; the pipeline degrades
   gracefully without it.

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
   quarto render reports/seminar_summary.qmd
   quarto render reports/isoform_switching_overview.qmd
   quarto render reports/ht_switch_analysis.qmd
   quarto render reports/novel_isoform_analysis.qmd
   quarto render reports/ut_t_vs_u_comparison.qmd
   ```

   Committed HTML (self-contained via `embed-resources`):

   - `reports/seminar_summary.html` — narrative summary of the analysis and findings
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
