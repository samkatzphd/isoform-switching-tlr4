# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this repo is

A **downstream-only** R analysis pipeline for IsoformSwitchAnalyzeR (ISA) results. It starts *after* differential isoform analysis: inputs are saved `SwitchAnalyzeRlist` objects, not raw RNA-seq. There is no build system, no test suite, and no linter config — the deliverables are the numbered scripts in `scripts/`, the tables/figures they write under `results/`, and the Quarto reports in `reports/`.

Four datasets across two reference transcriptomes (config keys): `T_HT`, `H_HT` (HT reference, Q1: WT ± LPS), `T_UT`, `U_UT` (UT reference, Q2: UBL5 knockout vs WT). See `README.md` for the biological framing and a per-script description.

## Commands

Run everything from the project root (or set `ISOFORM_PROJECT_ROOT`).

```bash
Rscript scripts/99_verify_inputs.R        # gate: checks every configured input, exits non-zero if bad
Rscript scripts/01_load_data.R            # must run first; writes data/processed/
Rscript scripts/02_gene_level_summary.R   # must run after 01
Rscript scripts/02b_qc_snapshot.R
Rscript scripts/02c_qc_figures.R
Rscript scripts/02d_expression_diagnostics.R   # abundance floor + between-group expression
Rscript scripts/02e_response_magnitude_control.R  # LPS-response stratification; see the caveat below — NOT a confounder control
Rscript scripts/02f_reference_concordance.R       # same-sample concordance across references
Rscript scripts/03_novel_isoform_analysis.R
Rscript scripts/04_comparison_T_vs_U.R
Rscript scripts/06_visualization.R
Rscript scripts/07_isg_analysis.R         # curated interferon gene-set test
Rscript scripts/00_threshold_scan_abs_dif.R   # optional sensitivity scan
Rscript scripts/05_pathway_enrichment.R  # GO/Reactome/KEGG; needs context tables

# Every script takes an optional alternate config as argv[1]:
Rscript scripts/02_gene_level_summary.R config/other.yml

quarto render reports/isoform_switching_overview.qmd   # also: ht_switch_analysis,
                                                       # novel_isoform_analysis,
                                                       # ut_t_vs_u_comparison, seminar_summary
```

Rendered HTML is committed alongside each `.qmd` (self-contained via `embed-resources: true`).

### External data dependency

`config/config.yml` points `isa_path` / `annotation_path` at `/Volumes/Expansion/...` — an external drive. Scripts `01`, `04`, and `06` read those files directly. **Without the drive mounted, `01` fails and `04`/`06` warn and skip all plotting**, deliberately leaving existing figures in place rather than overwriting them with fallbacks. Scripts `00`, `02`–`03` and report rendering only need `data/processed/` and `results/`, which are committed, so most work can proceed offline.

## Architecture

### Script skeleton (repeated verbatim in every script — copy it for new ones)

1. Locate and `source()` `utils/bootstrap.R`, trying `utils/bootstrap.R` then `../utils/bootstrap.R`.
2. `bootstrap.R` resolves the project root (env `ISOFORM_PROJECT_ROOT`, else walk up ≤20 levels looking for `config/config.yml`), sources `utils/helper_functions.R` into the global env, and assigns `PROJECT_ROOT`.
3. Read `commandArgs(trailingOnly = TRUE)[[1]]` as an optional config path, defaulting to `config/config.yml`.
4. `resolve_path(...)` + `ensure_dir(...)` for every output directory; nothing is hard-coded to an absolute path.

`%||%` (null-coalesce) is defined in `helper_functions.R` and used pervasively for config defaults.

### Config keys vs. labels — the main filename gotcha

Config datasets are keyed `T_HT`, `H_HT`, `T_UT`, `U_UT`, but each carries a `label` (`T_HT_anchor`, `H_HT_supplementary`, `T_UT`, `U_UT`), and **output filenames use the sanitized label, not the key**: `data/processed/isoformFeatures_T_HT_anchor.rds`, `results/tables/gene_level_summary_H_HT_supplementary.csv`. When adding code that looks up a dataset's files, go through the config `label` (see `load_processed()` in `04_comparison_T_vs_U.R`) rather than assuming key == filename.

### Output convention

Every result table goes through `write_table_pair(x, dir, stem, cfg)` in `utils/`: `.rds` always, `.csv` twin unless suppressed. Wide isoform-level tables pass `csv = FALSE` (config `output.csv_twin_isoform_level`) because their CSV twins added megabytes to git on each regeneration. Reports read the `.csv` side of summary tables. Figures are PNGs under `results/figures/<topic>/`. Each script ends by calling `write_run_manifest()`, which appends provenance to `results/run_manifest.json`.

### The "switching" definition lives in exactly one place

`score_isoforms(iso, cfg, dataset_key, dataset_label)` in `utils/helper_functions.R` is the only implementation of the criterion — isoform q < `significance.isoform_q`, gene q < `significance.gene_q` (when the column exists), `|dIF| >= significance.min_abs_dif`, and finite values required unless `significance.require_finite_dif: false`. It returns the per-isoform columns (`q_i`, `q_g`, `dIF_n`, `abs_dIF`, `is_novel`, `is_switching`) that `02`/`03`/`04`/`06` all consume. Scripts used to carry their own copies and had drifted apart; do not reintroduce a local copy. `config/config.yml` is the only place thresholds are numerically defined.

### Gene summaries are one row per gene

`score_isoforms()` calls `collapse_gene_symbols()` first, so `group_by(gene_id, gene_name)` cannot split a gene across rows when its isoforms disagree on the symbol (11–18% of gene_ids here). `02` and `02b` assert `rows == distinct gene_id`. Two effect-size columns exist deliberately: `max_abs_dif_all_isoforms` (all isoforms) and `max_abs_dif_switching` (switching only, used for ranking). They are different quantities — don't merge them back into one name.

### Gene symbol handling

ISA objects carry `XLOC_*` placeholder gene names. `01_load_data.R` joins the reference GFF3 transcript map (`build_transcript_map_from_gtf()`) and overwrites placeholder `gene_name` values with GFF-derived symbols, preserving the original in `gene_name_original`. `is_real_gene_symbol()` / `pick_gene_symbol()` in `utils/` filter out `XLOC_*` and `ENSG/ENST/...` identifiers. Cross-dataset overlaps in `03`/`04` join on **clean gene symbols**, not `gene_id`, since gene ids differ between the HT and UT references.

### The inputs are pre-reduced — and the context layer is how that's handled

All four primary ISA objects were saved after `isoformSwitchTestDEXSeq(reduceToSwitchingGenes = TRUE)`, so every gene in them already passes the gene-level q cutoff (recorded per run in `results/tables/input_object_reduction_check.csv`).

`isa_unfiltered_path` in config points a dataset at an unreduced export; `01` turns it into `data/processed/isoformContext_<label>.rds` (+ `isoformContextRepIF_`). **Configured for all four datasets.** Load it with `load_context_table(processed_dir, label, "features"|"rep_if")`, which returns `NULL` when absent — always handle that branch, and prefer using whichever side is available over requiring both.

Rules that must hold:

- **Significance comes FROM the context layer** via `load_scoring_table()`, which prefers the context and falls back to the reduced object with a warning. Never take q-values from one object and presence from the other — they are different FDR universes (dIF matched across 966 U_UT isoforms but only 132 q-values did), and that mismatch misclassified genes.
- **`T_HT` and `T_UT` are the SAME six libraries** under two references (verified from replicate sample names; gene-level r ≈ 0.99). Their agreement is technical reproducibility, never biological replication, and statistics must not treat them as independent. `TCONS_*` ids are assigned per reference and are NOT comparable across them — compare on gene symbols. Same-sample switching Jaccard is 0.35 — but **do not call this a "reproducibility ceiling"**: HT merges T+H PacBio transcripts and UT merges T+U, so the two are not one transcript space viewed twice, and part of that disagreement is genuine annotation-space difference rather than measurement noise. The decomposition of the discordance still holds; the ceiling framing does not.
- **`02e` is NOT a confounder control, and its result must not be used to explain away the UT deficit.** The Mantel-Haenszel stratification on LPS-response magnitude (UT: crude OR 0.51 → MH 0.77, 95% CI 0.57–1.06, p = 0.11) is arithmetically correct but **causally over-adjusted**. UBL5/Hub1 is a spliceosome-associated modifier, so if splicing regulation is upstream of the LPS response then response magnitude is a **mediator**, not a confounder — and conditioning on a mediator removes the very effect being measured. Both causal models predict the same correlation, so this design cannot distinguish them. Keep the numbers as a record (`ut_response_magnitude_control.csv`, `review_switching_by_response_stratum.csv`); never present them as evidence that the deficit is an artefact. The HT arm (2.29 → 2.49) is separately invalid — see the `H_HT` rule below. Full withdrawal in `docs/external_review/HANDOFF.md` §"Superseded — do not implement".
- **`H` and `T`/`U` are different kinds of experiment.** All four have real replication — `T1`–`T3` were matured, LPS-treated and RNA-extracted *separately*, not split from one flask. But `T`/`U` are one clonal line processed on a single day, while `H` is human **primary** macrophages from **three different donors**, each prepped on a **different day**. `H` therefore carries two variance components (donor genotype, prep day) that `T`/`U` structurally lack — and in `H` those two are perfectly confounded, one donor per day, so genotype and batch can never be separated there.
- **`H_HT` alone is batch-corrected, so its IF/dIF are not on the same scale as the other three.** Following from the above, `importRdata` found 2 surrogate variables for it and none for `T_HT`/`T_UT`/`U_UT`, so it alone had `limma::removeBatchEffect` applied to its expression and all its IF values derived from the corrected matrix. This halves its replicate IF noise (0.0062 → 0.0030, inverting its noise ranking) and accounts for its 2.07%-vs-0.87% switching rate. **Never compare an `H_HT` rate or dIF distribution to another dataset's** — and note that the H-vs-T switching contrast is not rescued by reprocessing, since it confounds genotype with design; retire it rather than recompute it. The other three reproduce exactly from the raw RSEM counts (dIF r = 1.0000); `H_HT` reproduces only once the correction is reapplied (r = 0.777 → 0.99915). The fix is an **explicit** donor blocking factor (`~ donor + condition`), *not* `detectUnwantedEffects = FALSE`, which would leave real donor variation uncorrected. Full detail in `docs/REVIEW_CHANGES.md` §0g.
- **An abundance floor applies**: `significance.min_gene_expression` (12), derived in `02d` from replicate IF noise versus gene expression. `score_isoforms()` applies it, flags `low_expression`, and keeps `is_switching_before_expr_floor` so the cost is auditable. The floor is on GENE expression — IF is isoform/gene, so gene abundance sets the precision.
- **With context**, overlap classes say *tested in X, not switching* vs *not detected in X*. **Without it**, they must say *retained / not retained* — the reduced object cannot distinguish tested-and-negative from dropped.
- **Pathway enrichment (`05`) uses the context layer as its universe.** Never substitute a whole-genome background or the pre-reduced object — the script skips a dataset without a context rather than falling back. Report term counts alongside `pathway_enrichment_driver_summary.csv`: GO/Reactome nesting means a handful of genes can produce ~150 terms.
- Two kinds of enrichment statistic exist and mean different things: `ut_T_vs_U_background_enrichment` / `ht_T_vs_H_background_enrichment` (valid — 11,126 and 7,671 tested genes) versus the retention accounting (descriptive, 25 already-significant genes, no test). Don't merge them.

Full detail in `docs/REVIEW_CHANGES.md`. `docs/STATISTICAL_METHODS.md` explains every method and the traps this project has hit; `docs/OPEN_QUESTIONS.md` lists what is unresolved. Read both before proposing a new analysis. `docs/AGENT_CONTEXT.md` holds the
condensed version pasted into the project's Agent Context field — if you change an invariant
here, update that file too so the two do not drift.

### Novel isoform flag

`tag_novel_isoforms()` marks rows whose `config$novel$id_column` (`oId`) starts with `config$novel$pb_prefix` (`PB`) as `is_novel_pacbio`; this rolls up to `novel_involved` at gene level.

### Reports

All four `.qmd` files start with `source("_setup.R")`, which supplies `PROJECT_ROOT`, `proj_path()`, `tab()`, `figs()`, `make_fig()`, `require_outputs()`, `read_tab()`, `existing_figs()` and `reduction_caveat()`. Paths are **absolute**; do not set `knitr` `root.dir` — it breaks `include_graphics()` combined with `embed-resources`, which is why the reports originally diverged into three different setups.

Call `knitr::include_graphics()` **once on a vector**, never inside a `for` loop — the looped form returns a value that is never printed, which silently produced an HT report containing zero images.

Reports are read-only consumers of `results/tables/*.csv` and `results/figures/`; each guards with `require_outputs(..., "Run: Rscript scripts/NN_...")`. `reports/exploratory.qmd` is a scratch file with no real content.

## Style

Base R + `dplyr`/`tibble`/`ggplot2` with the native pipe `|>` and `.data$` pronouns inside dplyr verbs. `library()` calls are wrapped in `suppressPackageStartupMessages({...})`; optional packages (`ggplot2`, `scales`, `IsoformSwitchAnalyzeR`) are gated behind `requireNamespace(..., quietly = TRUE)` so scripts degrade rather than fail. Progress is reported with `message()`, and scripts end with `message("NN_name.R: done")`.
