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

> **Superseded in part by §0c.** The counts in this section were produced while significance
> still came from the reduced objects and before the abundance floor existed. They are kept
> as the record of what changed at the time; for current numbers see §0c and the generated
> tables, which every report reads live.

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

### HT context (added 2026-08-06, completed 2026-08-06)

Unfiltered HT exports were added under `H-T_Comparisons/*_Unfiltered_2025-Oct-12/`. The
first copy of `T_HT` was corrupt (`gzip -t` data stream error, ~269 MB of an expected
~2.6 GB); it was re-copied and now passes. **All four datasets now have a working context.**

| dataset | isoforms tested | genes tested | % of genes significant |
|---|---|---|---|
| T_HT | 148,951 | 10,907 | 23.3 |
| H_HT | 128,491 | 8,160 | 20.5 |
| T_UT | 166,394 | 11,582 | 24.9 |
| U_UT | 173,016 | 11,969 | 16.4 |

`06` uses whichever context side is available rather than requiring both, and
`extract_t_h_overlap()` joins each side's measured `dIF`/`IF1`/`IF2`/q so the tables can say
what happened in the other dataset, not merely whether it was tested there.

**This changed a headline number.** Judged by what survived `H_HT`'s reduction, 7 of 15
T-significant isoforms were also significant in H — an apparent replication of ~47%. Against
the 115 isoforms **actually tested** in H, the rate is **6.1%**: inflated roughly eight-fold,
for exactly the reason described in §0.

Both directions are now measurable, each conditioned on significance only once:

| direction | tested in the other | % also significant | % same direction | median \|dIF\| there | % clearing 0.15 |
|---|---|---|---|---|---|
| T-significant → measured in H | 115 | 6.1 | 74.8 (p = 9.8e-8) | 0.07 | 18.3 |
| H-significant → measured in T | 289 | 2.4 | 63.3 (p = 6.9e-6) | 0.05 | 8.0 |

Direction is preserved well above chance in both directions while magnitude collapses below
the calling threshold — the same signature as the UBL5 arm, more extreme.

`ht_T_vs_H_background_enrichment` mirrors the UT background test: over **7,671 genes
quantified in both** HT objects, 11 switch in both against 2.3 expected — OR 4.8, p = 6.4e-5.
Notably weaker than the UT pair (OR 25.7), which is the expected ordering: `T_UT` and `U_UT`
differ only by the UBL5 knockout, whereas `T_HT` and `H_HT` are different genotypes.

**Operational note.** Two exports arrived damaged and neither was visible from a file
listing. `scripts/99_verify_inputs.R` now checks every configured input for existence and
compression integrity (detected from magic bytes, not the extension) and exits non-zero on
failure, so it can gate a run. The drive is **exFAT**, which has no journaling — an
interrupted write leaves a partially written file with no error at the time. The README
carries a download/unzip/copy checklist built around that.

---

## 0c. Significance source and the abundance floor (2026-08-07)

Two corrections that change every switching count, both surfaced by following up a single
gene (CD86) that looked misclassified.

### The pipeline was mixing two multiple-testing universes

Significance was taken from the reduced objects while presence and effect were taken from
the unfiltered context. Across the 966 shared U_UT isoforms:

- **dIF identical in all 966** — same model fit
- **isoform q identical in only 132**; 13 isoforms significant in the reduced object only,
  16 in the context only

CD86 was the visible symptom: in U it passes every criterion (isoform q 3.6e-3, gene q
3.6e-3, |dIF| 0.26) yet was called **T-only**, because U's reduced object never retained it.

**Fixed:** `load_scoring_table()` prefers the unfiltered context, and `02`/`03`/`04`/`06`
all use it, so significance and presence share one FDR universe. The source is reported per
dataset at run time and falls back to the reduced object with a loud warning.

### Switching calls were concentrated where the estimator is unreliable

dIF is a ratio, so at low gene abundance the isoform fraction is estimated from few reads.
`02d_expression_diagnostics.R` measures this directly: the within-condition SD of IF across
replicates (condition fixed, so pure noise) against gene expression.

| gene expression | median within-condition IF SD | 0.15 / SD | % isoforms called switching |
|---|---|---|---|
| 2–6 | 0.086 | 1.7 | 1.51 |
| 6–12 | 0.062 | 2.4 | 1.14 |
| **12–20** | **0.050** | **3.0** | 0.74 |
| 31–44 | 0.038 | 3.9 | 0.26 |
| >274 | 0.017 | 8.6 | 0.24 |

Below gene expression ~12 a threshold-sized "switch" sits within ~2 SD of replicate noise,
and calls are ~6x enriched there. The pattern holds within every isoform-count stratum, so
it is abundance and not transcript complexity.

**Floor applied: `significance.min_gene_expression: 12`**, the point where |dIF| = 0.15
reaches 3x replicate noise in all four datasets. Excluded calls are also flagged
(`low_expression`, plus `is_switching_before_expr_floor`) so the cost stays visible.

Cost — roughly a third to a half of prior calls, which is the intent:

| dataset | switching isoforms / genes before | after |
|---|---|---|
| T_HT | 213 / 169 | 125 / 95 |
| H_HT | 288 / 258 | 192 / 169 |
| T_UT | 264 / 219 | 143 / 118 |
| U_UT | 153 / 128 | 79 / 62 |

CD86 illustrates both fixes: gene expression 15.2 in T but **3.0 in U**, so its U call is now
excluded as unreliable rather than silently absent. It is still T-only, but for a measured
reason instead of an artefact of object retention.

### Global power is not the confounder

Worth ruling out, since it is the obvious alternative explanation for "the knockout switches
less". On the 11,126 genes tested in both UT arms the median log2 expression ratio is
**0.003** and 50.3% of genes are lower in U — a coin flip. (Paired Wilcoxon p = 7e-9 is
sample size, not effect.) The HT pair is similar at 0.064. So the between-genotype
differences are not a depth or library-size artefact.

---

## 0d. External review follow-up (2026-08-07)

An external review (`REVIEW_HANDOFF.md`) raised three findings. All were independently
reproduced from the committed context tables before anything was changed.

### The Q2 conclusion did not survive its confounder

`02d` ruled out a *baseline abundance* confounder, which was the wrong one. The relevant
confounder for a post-LPS switching claim is the amplitude of the LPS response itself, and
it is large: over 8,595 genes tested in both UT arms the knockout's response has SD ratio
0.731, IQR ratio 0.577, regression slope **0.516** and Spearman 0.727 against wildtype. The
top 200 wildtype responders retain 57.9% of their magnitude in the knockout. Direction is
preserved; amplitude is roughly halved.

Stratifying each gene by its own response magnitude, the crude switching odds ratio of
**0.508** becomes a Mantel-Haenszel **0.773 (95% CI 0.565-1.058, p = 0.114)** — not
significant, and in the two highest-response strata the knockout switches slightly more.

**Changed:** `scripts/02e_response_magnitude_control.R`, run on both pairs. Reports reframed
to the supported claim: UBL5 loss blunts the LPS response globally by about half, and
reduced isoform switching follows from that.

Two things the review did not note, both recorded in the script:

- Response magnitude is plausibly a **mediator**, not a confounder. Adjusting for it removes
  part of the effect being measured, so a non-significant adjusted OR shows the design cannot
  separate a splicing-specific effect from a globally blunted response — not that none
  exists.
- Running the same control on the HT pair gives the opposite behaviour: crude OR 2.29 becomes
  MH OR **2.49 (p = 5.8e-13)**, i.e. the contrast *strengthens* on adjustment. That is a
  useful internal control — the stratification is not mechanically flattening everything.

### T_HT and T_UT are the same libraries

Replicate columns carry identical sample names including S-numbers, and gene-level expression
correlates at **r = 0.991-0.992** for matched samples (versus 0.943 same-library /
different-condition). Nothing in the project had noted this.

A caution for anyone repeating the check: correlating the raw matrices by `isoform_id` gives
**r ≈ 0**, which looks like different libraries but is an id-collision artefact — `TCONS_*`
ids are assigned per reference and denote different transcripts in HT and UT. All
cross-reference comparison must use gene symbols.

This is the project's best technical control: zero biological variability, so all
disagreement is annotation plus thresholding. Over 10,592 symbols testable in both, 95 and
116 switching genes share 55 — Jaccard **0.353** against 1.04 expected by chance, with
gene-level max |dIF| Spearman 0.771. Of 101 discordant genes, median |dIF| in the
non-calling reference is 0.132 against a 0.15 cutoff and 72.3% reach at least 0.10, so this
is threshold brittleness, not contradiction.

**Changed:** `scripts/02f_reference_concordance.R`, which asserts library identity from the
replicate column names before reporting. Reports now state that a single-dataset switching
call has roughly a one-in-three chance of not replicating under re-annotation alone, and
that cross-dataset replication rates should be judged against ~35% rather than 100%.

### ISG replication was partly one dataset counted twice

`07` corrected across all 12 dataset x set rows as if independent, and the seminar described
`T_HT` and `T_UT` agreement as replication. They are the same libraries; their switching ISGs
overlap 7 of 9 in the union.

**Changed:** `07` derives a `library_group` from the replicate sample names and reports both
`q_across_all_datasets` (optimistic, treats all rows as independent) and
`q_within_library_group`. Reports describe HT/UT agreement as technical reproducibility, and
show switching ISG counts (8, 8, 4, 3) beside every odds ratio.

### Smaller items

- **Switching is tail behaviour.** Above the abundance floor at mid-range fractions, the 99th
  percentile of |dIF| is 0.147 / 0.153 / 0.134 / 0.201 across the four datasets, so the 0.15
  cutoff sits near the 99th percentile of the whole distribution. Stated in the reports.
- **Possible paired design.** Sample names suggest three paired replicates but no blocking
  factor is used anywhere. Flagged in the reports as an open question with its power cost;
  **not** changed, since inferring a design from filenames is not sound.
- **Pathway term counts** were left in place rather than removed. The driver-gene summary is
  reported alongside every count, which addresses the misreading risk while keeping the
  count informative.

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
- **`05_pathway_enrichment.R` was a stub** while the only available gene lists came from
  pre-reduced objects, since enrichment against them would have inherited the §0 problem.
  Implemented 2026-08-06 once all four datasets had a context layer to supply a real tested
  universe. The result is largely negative: four of seven gene sets return no enriched terms,
  and the positive ones all reduce to the same five chemokines/cytokines (CCL3, CCL4, CCL22,
  IL23A, TNFSF4), which account for 95% of gene-term hits in the shared UT set. A
  `pathway_enrichment_driver_summary` table is written alongside the term counts so the
  ~150-terms-from-18-genes figure cannot be read as 150 independent findings.

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
