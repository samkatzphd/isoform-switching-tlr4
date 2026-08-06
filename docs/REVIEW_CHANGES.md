# Review of the Cursor-era pipeline, and what changed

Branch: `fix/review-corrections-2026-07`
Review date: 2026-07-30
Baseline reviewed: `c15a108` (all work up to that commit was done in Cursor)

This document exists so the changes can be audited later. Each item states **what was
wrong**, **the evidence**, and **what was done**. Items are ordered by how much they
affect reported numbers.

Nothing here was found by reading alone — every claim was checked against the actual
data in `data/processed/`, `results/tables/`, and the source `.Rdata` objects.

---

## 0. The finding that changes interpretation, not code

**The input ISA objects were already reduced to significant switching genes.**

Verified directly against the source object:

```
T_minus-T_plus_IsoformSwitchAnalyzeR.Rdata  ->  isa_list
  isoformFeatures rows: 1769   genes: 248
  genes with gene_switch_q_value < 0.05: 248 of 248
  max gene_switch_q_value in the entire object: 0.0487
```

All four datasets behave the same way (216/216, 385/385, 127/127, 248/248). This is
`isoformSwitchTestDEXSeq(reduceToSwitchingGenes = TRUE)` — the ISA default — applied
before the object was saved.

Everything downstream that counts genes or assumes a background inherits this:

- Gene counts are **genes retained**, not **genes tested**.
- "Present in the other dataset" means **retained in that saved object**, not **tested
  and non-significant**.
- Only **25 gene symbols** and **172 isoforms** appear in both UT objects at all, yet the
  old summary reported 138 T-only and 60 U-only genes. Those T-only genes are
  overwhelmingly absent from the U object, not tested-and-negative.
- The threshold scan in `00` confirms the ceiling: shared genes never exceed 25 at *any*
  `min_abs_dif`, because 25 is the entire T∩U symbol overlap.

**What was done (this repo):** `01_load_data.R` now runs `detect_isa_reduction()` per
dataset and writes `results/tables/input_object_reduction_check.csv`. Reports render a
callout built from that file rather than a hard-coded claim. Overlap language throughout
was changed from *missing/present* to *not retained/retained*.

**Resolved for the UT datasets — see §0b.** The original conclusion here was that the fix
required re-exporting the objects. It turned out unfiltered exports already existed on the
drive for both UT datasets, which is a better outcome: no re-analysis was needed.

---

## 0b. The unfiltered context layer (added 2026-08-05)

`*_unfilteredR.Rdata` exports exist on the drive for **both UT datasets** — the same
analyses without the reduction step:

| | isoforms | genes | genes passing gene q<0.05 |
|---|---|---|---|
| T_UT unfiltered | 166,394 | 11,582 | 2,879 (24.9%) |
| U_UT unfiltered | 173,016 | 11,969 | 1,966 (16.4%) |
| T_UT reduced (primary) | 1,769 | 248 | 248 (100%) |

Verified consistent: all 1,769 isoforms of the reduced T_UT object appear in its
unfiltered counterpart, and gene q spans 1.4e-192 to 0.9999 rather than stopping at
0.0487. **11,126 gene_ids are shared between the two unfiltered objects, against 25
symbols in the reduced ones.**

**Design: context layer, not replacement.** Significance calls still come from the primary
(reduced) objects — they are identical values either way. The unfiltered data is used for
display and classification only. This matters because script `03` compares novelty rates
between HT and UT; had UT switched to a 166k-isoform denominator while HT stayed at 1.5k,
that figure would have quietly become meaningless.

`01_load_data.R` writes, per dataset with `isa_unfiltered_path` configured:

- `data/processed/isoformContext_<label>.rds` — per-isoform IF1/IF2/dIF, expression, q-values
- `data/processed/isoformContextRepIF_<label>.rds` — per-replicate isoform fractions

What it enables:

1. **Gene explorer panels show both genotypes for any gene.** Of 145 genes switching in T
   but not U, **127 (88%) now have U data to plot**; 60 of 66 in the other direction. The
   panels were rebuilt to show isoform fraction before and after LPS with an arrow for the
   shift and replicate points overlaid, rather than a single dIF bar — a ratio that moves
   in one genotype and holds in the other is now directly visible.
2. **Real overlap classification.** *Tested in the other dataset and not switching* is now
   distinguishable from *not detected*:

   | isoform class | n |
   |---|---|
   | T-significant / tested in U, not switching | 161 |
   | U-significant / tested in T, not switching | 77 |
   | Shared significant | 21 |
   | T-significant / not detected in U | 22 |
   | U-significant / not detected in T | 8 |

   238 isoforms are present in both with a ratio change in only one; only 30 are genuinely
   absent. Under the reduced objects all 268 were indistinguishable.
3. **A legitimate enrichment test**, reported as `ut_T_vs_U_background_enrichment`: over
   11,126 genes quantified in both, 18 switch in both against 1.02 expected under
   independence — OR 25.7 (95% CI 13.9–45.7), p = 4.4e-18. This is *not* the statistic
   removed in §1.1; that one had a 25-gene background selected on the outcome. Both are
   kept in the outputs under clearly different names so they cannot be confused.

**Note on §2.9:** the reduced objects contain no NA q-values, so the `require_finite_dif`
fix was described there as changing nothing. The unfiltered objects contain 186 and 142 NA
isoform q-values, so that fix is load-bearing for anything touching the context layer.

### HT context (added 2026-08-06)

Unfiltered HT exports were added under
`H-T_Comparisons/*_Unfiltered_2025-Oct-12/`. Status is split:

- **`H_HT` — in place and working.** 128,491 isoforms over 8,160 genes, 20.5% of genes
  significant, gene q spanning 7.8e-45 to 1, all 3,102 reduced isoforms present.
- **`T_HT` — corrupt, not usable.** `gzip -t` reports a data stream error; the file
  decompresses ~269 MB of an expected ~2.6 GB (the H file, of near-identical compressed
  size, yields 2.6 GB). The header is fine and the writing R version matches files that
  load correctly, so this is a bad transfer, not a compatibility problem. Needs re-copying
  from source. `config` keeps `isa_unfiltered_path: null` for `T_HT` with that noted.

`06` was relaxed so it uses whichever side is available — statements about H need only the
H context, and requiring both would have withheld the entire result over one bad file.

**This changed a headline number.** Judged by what survived `H_HT`'s reduction, 7 of 15
T-significant isoforms were also significant in H — an apparent replication rate of ~47%.
Against the 115 isoforms **actually tested** in H, the rate is **6.1%**: the biased figure
was inflated roughly eight-fold, for exactly the reason described in §0. The T-significant
switches are not absent from H, though — 74.8% move in the same direction (86 of 115,
binomial p = 9.8e-8) with a median |dIF| of 0.07, below the 0.15 calling threshold. Same
attenuation signature as the UBL5 arm.

`extract_t_h_overlap()` now also joins H's measured `dIF`/`IF1`/`IF2`/q from the context,
so the table can say what happened in H rather than only whether it was tested there.

**Still outstanding:** the corrupt `T_HT` export. Until it is replaced, the reciprocal
question — are H-significant switches merely sub-threshold in T? — cannot be asked. RSEM
count matrices for both references are on the drive under `gtf_files/*/[HU]_T_mapped/` if a
re-analysis is ever wanted, but that was explicitly out of scope.

---

## 1. Statistics that were not interpretable as reported

### 1.1 Fisher test removed, not caveated

`04_comparison_T_vs_U.R` ran `fisher.test()` on a 2×2 of *significant in T* × *significant
in U* over symbols present in both datasets, and reported the odds ratio and p in the
headline summary (`OR = 15.0, p = 0.031`).

Because "present" already means "significant", the background was **25 genes** and the
table was **3 / 1 / 3 / 18**. It measured association inside a set selected on the
outcome.

**Changed:** the test is gone. `ut_T_vs_U_retention_accounting.csv` and
`ut_T_vs_U_retention_contingency.csv` report the same counts descriptively, with an
explicit note. `fisher_odds_ratio_*` / `fisher_p_*` no longer appear in
`ut_T_vs_U_stats_summary`.

### 1.2 The `|dIF|`-by-class comparison was biased by construction

The Shared group was summarised with `pmax(T_abs_dIF, U_abs_dIF)` while T-only and U-only
used a single dataset's value. A maximum of two draws is upward-biased against one draw,
so "Shared isoforms have larger |dIF|" (Kruskal p = 9.3e-5) was partly an artefact of the
summary statistic.

**Changed:** every isoform is now measured in exactly one dataset (Shared and T-only in T,
U-only in U), recorded in a new `effect_measured_in` column.

### 1.3 Concordance statistics conditioned on the outcome twice

Spearman ρ = 0.91 and 21/21 same-direction (binomial p = 9.5e-7) were computed on
isoforms selected for passing the same threshold in *both* datasets. Selecting on
agreement and then measuring agreement is circular.

**Changed:** those columns are retained but prefixed **`selconf_`** and documented as
selection-conditioned. A new `ut_T_vs_U_one_way_concordance.csv` provides the defensible
version: condition on significance in one dataset, then measure median |dIF|, direction
agreement, and threshold-passing rate in the other.

### 1.4 Effect sizes now lead, p-values follow

With thousands of matched isoforms the Wilcoxon p-values are driven by n (median ΔdIF was
0.0029). `ut_T_vs_U_stats_summary` now lists medians and counts before p-values, and the
report says so in a "what can and cannot be concluded" section placed above the numbers.

---

## 2. Correctness bugs

### 2.1 Gene summaries had more rows than genes

`summarize_genes_from_isoform_table()` grouped by `(gene_id, gene_name)` without
collapsing symbols first, so a gene whose isoforms disagreed on the symbol was split
across rows. 11–18% of gene_ids are affected.

| dataset | old rows | actual genes |
|---|---|---|
| H_HT_supplementary | 447 | 385 |
| T_HT_anchor | 273 | 216 |
| T_UT | 293 | 248 |
| U_UT | 157 | 127 |

QC "n_genes" and every `pct_*` denominator were inflated, and joining anything on
`gene_id` fanned out (joining 30 top genes against the summary returned 41 rows).
`03_novel_isoform_analysis.R` had the same pattern; `04`/`06` did not, because they
collapsed symbols first — which is how the scripts had silently diverged.

**Changed:** `collapse_gene_symbols()` runs inside `score_isoforms()`, so all four
scripts share one behaviour. `02` and `02b` now `stopifnot(rows == unique gene_ids)`.
Corrected counts: 385 / 216 / 248 / 127.

### 2.2 Two different effect sizes under one name

The old `max_abs_dif` in `gene_level_summary_*` was the max over **all** isoforms, while
`03`/`04`/`06` ranked on max over **switching** isoforms under the near-identical name
`max_abs_dif_switching`. They disagreed for 11 of 41 joined top-30 rows (SH3BGRL2: 0.321
vs 0.039). `02b` ranked its QC top-10 on one and `06` on the other, so the two "top gene"
lists were not comparable.

**Changed:** the summary emits both `max_abs_dif_all_isoforms` and
`max_abs_dif_switching`. `02b` ranks on the switching version; `02c` plots the
all-isoform version with a title that says so.

### 2.3 `switchPlot` device leak could overwrite a good plot with a failed one

`on.exit()` inside `tryCatch({...})` registers on the **enclosing function** frame, not
per iteration. Confirmed by test: devices accumulated and closed LIFO at function return.
Because both candidates wrote to the same path, when the `gene_name` attempt failed and
the `gene_id` attempt succeeded, the failed device — opened first, therefore closed last —
overwrote the good plot. Separately, the `file.info()$size > 0` check ran before the PNG
had been flushed.

Latent rather than active: 58 of 60 committed HT plots are real 2400×1500 switchPlots.

**Changed (`04` and `06`):** each attempt renders to its own temp file, closes its device
immediately, and is promoted to the destination only after being verified non-empty.

### 2.4 The HT report contained no images at all

`ht_switch_analysis.html` had **0** embedded images, against 5 / 9 / 36 in the other three
reports. The cause: `include_graphics()` called inside a `for` loop returns a value that
is never printed. All 60 switchPlots generated by `06` were missing from the report.

**Changed:** `include_graphics()` is called once on a vector.

### 2.5 Running `04`/`06` without the external drive destroyed figures

Both scripts read the ISA objects from `/Volumes/Expansion`. Without it mounted, `06`
deleted the existing rank plots and wrote fallback bar charts over them, and `04`'s
`if (!ok && file.exists(out_png)) unlink(out_png)` deleted good plots on every failed
attempt.

**Changed:** both scripts now detect the missing objects, warn, and skip plotting
entirely, leaving existing figures in place. `04` also skips rewriting
`ut_switch_plot_index` so it cannot be emptied.

### 2.6 The annotation join was gated on the wrong condition

`01_load_data.R` skipped the whole GFF3 join when the ISA table already had `oId`, which
would have silently dropped `class_code`, `cmp_ref`, and the XLOC→symbol repair too.

**Changed:** each annotation column is joined only if actually missing.

### 2.7 Which object gets analysed was decided by storage order

`object_name: null` let `load_isa_input()` auto-pick the first list/data.frame in the
`.RData`. Each file contains **13 objects**, including `args` (list),
`genome_wide_splicing`, `top_gene_switches`, and `top_isoform_switches` (data frames) —
any of which the picker could have selected. It happened to pick correctly.

**Changed:** `object_name: "isa_list"` pinned for all four datasets in config.

### 2.8 Inconsistent gene-symbol rule

`03` filtered only `^XLOC_`, so ENSG/ENST identifiers leaked into its "gene symbol"
overlap as if they were symbols; `04`/`06` filtered `^(XLOC_|ENS[GTFP]\d)`.

**Changed:** one `is_real_gene_symbol()` in `utils/`, used everywhere including `01`'s
symbol repair and `00`'s overlap scan.

### 2.9 The switching rule failed open on missing values

`is.na(dIF) | abs(dIF) >= min_abs_dif` meant an isoform with **no measured effect size**
passed the effect-size filter on its q-value alone. Same for a missing gene q.

Currently harmless — there are zero NAs in all four datasets — so this changes no number
today.

**Changed:** `significance.require_finite_dif: true` (new, default on) makes such rows
fail. Set it to `false` to restore the old behaviour.

---

## 3. Structural changes

### 3.1 One definition of "switching"

The rule existed in four places (`utils` + local `score_isoforms` in `03`/`04`/`06`) plus
an ad-hoc fifth in `00`. That is how 2.1 and 2.8 happened.

**Changed:** `score_isoforms(iso, cfg, dataset_key, dataset_label)` in
`utils/helper_functions.R` is the only implementation. `is_real_gene_symbol()`,
`pick_gene_symbol()`, `collapse_gene_symbols()`, `sanitize()`, `min_finite()`,
`max_abs_finite()` are shared too. Roughly 200 duplicated lines removed.

### 3.2 Provenance is now recorded

The pipeline advertised reproducibility but recorded nothing about code version, package
versions, or the thresholds behind a given table.

**Added:** `write_run_manifest()` writes `results/run_manifest.json` per script — git SHA,
timestamp, R version, platform, key package versions (including
IsoformSwitchAnalyzeR), the `significance` and `analysis` blocks in force, and for `01`
the per-dataset reduction flags.

### 3.3 Hard-coded parameters moved to config

`top_n_genes`, `plot_top_n`, `gene_panels_per_class`, `switch_plots_per_class`,
`explorer_max_isoforms` and the threshold-scan grid were literals inside scripts despite
the config-driven design. They are now in `config/config.yml` under `analysis:`, read via
`analysis_param()`.

### 3.4 One report setup instead of three strategies

The four `.qmd` files each carried their own project-root discovery, in three different
styles (absolute `proj_path()`, `knitr` `root.dir`, and `normalizePath("..")`).

**Added:** `reports/_setup.R` — `PROJECT_ROOT`, `proj_path()`, `tab()`, `figs()`,
`make_fig()`, `require_outputs()`, `read_tab()`, `existing_figs()`, `reduction_caveat()`.
Absolute paths throughout; `root.dir` is deliberately **not** set, because it breaks
`include_graphics()` combined with `embed-resources` (the original reason for the
divergence). Verified after the change: 5 and 9 images still embed correctly.

### 3.5 CSV twins dropped for wide isoform-level tables

`isoformFeatures_*.csv` (1.5 MB for H_HT) and the `*_isoform_overlap_*_all.csv` tables
were re-committed on every regeneration; nothing reads them, since reports read the
summary tables.

**Changed:** `output.csv_twin_isoform_level: false` in config. Those tables are `.rds`
only; every summary table keeps its CSV twin.

---

## 4. Deliberately not done

- **No `.Rprofile`.** The ~10-line bootstrap block at the top of each script is a
  chicken-and-egg minimum: it has to locate `utils/bootstrap.R` before any helper exists.
  An `.Rprofile` would only fire when R starts in the project root — the case that
  already works — while adding a hidden file that also affects interactive sessions and
  Quarto renders. The duplication is left in place as the lesser cost.
- **`DESCRIPTION` was not turned into a real package.** It remains a dependency manifest,
  not something installable; `jsonlite` was added to `Suggests`. Converting `utils/` into
  a package is a larger change than this branch should carry.
- **`05_pathway_enrichment.R` is still a stub** that stops with "Not yet implemented".
  Untouched deliberately — implementing enrichment on a gene set drawn from pre-reduced
  objects would inherit the item-0 problem and needs the upstream fix first.

---

## 5. How to verify this branch

```bash
# needs /Volumes/Expansion mounted for 01, 04, 06
Rscript scripts/01_load_data.R
Rscript scripts/02_gene_level_summary.R
Rscript scripts/02b_qc_snapshot.R
Rscript scripts/02c_qc_figures.R
Rscript scripts/03_novel_isoform_analysis.R
Rscript scripts/04_comparison_T_vs_U.R
Rscript scripts/06_visualization.R
Rscript scripts/00_threshold_scan_abs_dif.R

quarto render reports/isoform_switching_overview.qmd
quarto render reports/ht_switch_analysis.qmd
quarto render reports/novel_isoform_analysis.qmd
quarto render reports/ut_t_vs_u_comparison.qmd
```

Checks worth repeating:

- `results/tables/input_object_reduction_check.csv` — all four rows should show
  `looks_reduced_to_switching_genes = TRUE` until the objects are re-exported upstream.
- `gene_level_summary_*.csv` row count must equal its distinct `gene_id` count
  (`02`/`02b` assert this).
- `ht_switch_analysis.html` should now contain embedded images (it had none).
- `ut_T_vs_U_stats_summary.csv` should contain no `fisher_*` columns.

### Verified after full regeneration (2026-07-30)

| check | result |
|---|---|
| Embedded images per report | overview 5, novel 9, **HT 24 (was 0)**, UT 36 |
| HT switchPlots | 58 real (2400×1500), 2 fallback — unchanged, so the device fix cost nothing |
| `fisher_*` columns in stats summary | 0 |
| Gene summary rows == genes | asserted in `02`/`02b`, passes for all four datasets |
| Stale `ut_T_vs_U_fisher_gene_contingency.*` | deleted from the working tree and from git |

The one-way concordance that replaced the removed test is also the more informative
result. Of 35 isoforms significant in T and retained in U, **94% change in the same
direction** in U and 66% would pass the |dIF| threshold there; of 28 significant in U,
89% agree in direction and 79% would pass. Directional agreement is high — but note the
small n, which is itself a consequence of item 0: only 25 symbols survive in both
objects, so this is all the overlap that exists to measure.
