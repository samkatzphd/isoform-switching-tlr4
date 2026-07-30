# Resume note — review-corrections work in progress

**Branch:** `fix/review-corrections-2026-07` (off `main` @ c15a108)
**Last updated:** 2026-07-30
**Delete this file once the remaining work below is done** — it is superseded by
[`docs/REVIEW_CHANGES.md`](docs/REVIEW_CHANGES.md), which is the permanent record.

## State of this branch

All code changes from the review are **complete and parse clean**. What is *not* complete
is the regeneration of outputs, because scripts `04` and `06` read the ISA objects from
the external drive `/Volumes/Expansion`, which was unmounted partway through.

### Regenerated and consistent with the new code

- `data/processed/*.rds` (script `01`)
- `results/tables/gene_level_summary_*`, `qc_*`, `novel_*`, `min_abs_dif_threshold_scan`,
  `ut_overlap_vs_min_abs_dif`, `input_object_reduction_check`
- `results/figures/` QC, novel and threshold-scan figures
- `results/run_manifest.json`
- `reports/isoform_switching_overview.html`, `reports/novel_isoform_analysis.html`

### Stale — produced by the OLD code, do not cite

- `results/tables/ut_*` (all of them, from script `04`)
- `results/tables/top_switching_genes_*`, `isoform_overlap_T_HT_vs_H_HT_*`,
  `isoform_overlap_T_significant_in_H_context` (script `06`)
- `results/figures/ut_t_vs_u/**`, `results/figures/ht_top_switch_plots/**`
- `reports/ht_switch_analysis.html`, `reports/ut_t_vs_u_comparison.html`

`reports/ut_t_vs_u_comparison.qmd` will **fail to render** until `04` reruns, because it
now requires `ut_T_vs_U_retention_accounting.csv`, `ut_T_vs_U_retention_contingency.csv`
and `ut_T_vs_U_one_way_concordance.csv`, which the old `04` never wrote. That failure is
expected, not a bug.

## To finish

Mount `/Volumes/Expansion`, then from the project root:

```bash
Rscript scripts/04_comparison_T_vs_U.R
Rscript scripts/06_visualization.R
quarto render reports/ht_switch_analysis.qmd
quarto render reports/ut_t_vs_u_comparison.qmd
rm RESUME_NOTE.md
```

Expected differences from the committed (old) outputs:

- Gene counts drop where symbols were previously split across rows.
- `ut_T_vs_U_stats_summary.csv` loses every `fisher_*` column and gains `selconf_*`,
  one-way concordance and retention columns.
- Isoform overlap class labels change from *missing/present* to *not retained/retained*.
- `ht_switch_analysis.html` should contain embedded images — it previously had **zero**,
  because `include_graphics()` was being called inside a `for` loop.
- `04`/`06` now refuse to plot at all if the drive is missing, rather than overwriting
  good figures with fallbacks.

## Open item that cannot be fixed in this repo

The ISA objects were saved after
`isoformSwitchTestDEXSeq(reduceToSwitchingGenes = TRUE)`, so every overlap and background
statistic is conditioned on that selection. The real fix is upstream: re-export the
`switchAnalyzeRlist`s with `reduceToSwitchingGenes = FALSE`, which needs the original
count matrices (not in this repo). Section 0 of `docs/REVIEW_CHANGES.md` has the evidence.
