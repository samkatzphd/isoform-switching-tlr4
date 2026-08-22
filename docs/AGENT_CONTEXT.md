# Isoform switching pipeline — agent context

> This file is the source for the "Agent Context" field that is injected into every agent's
> system prompt for this project. Keep it tight: it is a system prompt, not documentation.
> Detail belongs in `CLAUDE.md` (architecture) and `docs/REVIEW_CHANGES.md` (what changed and
> why). If you edit this file, paste the updated content back into the Agent Context field.

## What this is

A **downstream-only** R analysis of IsoformSwitchAnalyzeR (ISA) results. It starts *after*
differential isoform testing: inputs are saved `switchAnalyzeRlist` objects, not raw
RNA-seq. No build system, no test suite, no linter. Deliverables are the numbered scripts
in `scripts/`, the tables/figures under `results/`, and the Quarto reports in `reports/`.

Four datasets, two reference transcriptomes:
- **Q1 (HT):** `T_HT` (WT THP-1 ± LPS, anchor), `H_HT` (H genotype, supplementary)
- **Q2 (UT):** `T_UT` (WT) vs `U_UT` (UBL5 knockout), both ± LPS

Run from the project root. Every script takes an optional config path as `argv[1]`.
`config/config.yml` is the only place thresholds and paths are defined.

## Non-obvious facts about the data — read before trusting any count

1. **The primary ISA objects are pre-reduced.** They were saved after
   `isoformSwitchTestDEXSeq(reduceToSwitchingGenes = TRUE)`, so every gene in them already
   passes the gene-level q cutoff. "Genes present" ≠ "genes tested". Recorded per run in
   `results/tables/input_object_reduction_check.csv`.

2. **Significance comes from the unfiltered context layer**, via `load_scoring_table()`.
   `isa_unfiltered_path` in config points each dataset at an unreduced export; `01` turns it
   into `data/processed/isoformContext_<label>.rds`. **Never take q-values from one object
   and presence from the other** — they are different FDR universes. Across 966 shared
   U_UT isoforms, dIF matched in all 966 but isoform q matched in only 132, and that
   mismatch silently misclassified genes.

3. **An abundance floor applies:** `significance.min_gene_expression: 12`. dIF is a ratio, so
   at low *gene* abundance the isoform fraction rests on few reads. `02d` derives the floor
   from within-condition replicate IF noise; below ~12 a threshold-sized switch is within
   ~2 SD of noise. It removes 33–52% of raw calls by design. The floor is on GENE
   expression, not isoform expression.

4. **`H` and `T`/`U` are different kinds of experiment, and `H_HT` alone is batch-corrected.**
   All four have real replication (`T1`–`T3` were matured, treated and extracted separately,
   not split from one flask), but `T`/`U` are one clonal line processed on a single day while
   `H` is human **primary** macrophages from **three donors** on **three different days** —
   two extra variance components, perfectly confounded with each other. So `importRdata` found 2 surrogate
   variables for `H` and none for the others, and applied `limma::removeBatchEffect` to H's
   expression, deriving all its IF/dIF from the corrected matrix. That halves its replicate IF
   noise and drives its 2.07%-vs-0.87% switching rate. **Never compare an `H_HT` rate or dIF
   distribution to another dataset's** — the H-vs-T contrast in `02e` confounds genotype with
   design and should be retired, not recomputed. The fix for H is an **explicit** donor
   blocking factor, *not* `detectUnwantedEffects = FALSE`. `REVIEW_CHANGES.md` §0g.

5. **The ISA objects have been reconciled against raw RSEM counts** at
   `/Volumes/Expansion/IsoformSwitchAnalyzer/Counts/`. `T_HT`, `T_UT` and `U_UT` reproduce
   exactly (dIF r = 1.0000); `H_HT` only after its correction is reapplied (0.777 → 0.99915).
   The counts cover all twelve libraries per reference, so a model-based DTU refit can test
   `genotype x treatment` in one 2x2 model.

6. **Raw inputs live on an external drive at `/Volumes/Expansion`**, which is **exFAT and has
   no journaling** — interrupted writes leave silently corrupt files. Two exports have
   arrived damaged. Always run `Rscript scripts/99_verify_inputs.R` before the pipeline; it
   exits non-zero on failure. Scripts `01`, `04`, `06` need the drive; the rest work from
   committed data.

## Invariants — do not break these

- **One switching rule.** `score_isoforms()` in `utils/helper_functions.R` is the only
  implementation. Scripts used to carry their own copies and drifted apart. Never
  reintroduce a local copy.
- **Gene summaries are one row per gene.** Symbols are collapsed before grouping; `02`/`02b`
  assert `rows == distinct gene_id`.
- **`max_abs_dif_all_isoforms` and `max_abs_dif_switching` are different quantities.** Don't
  merge them back into one name.
- **Never commit `.Rdata`.** `.gitignore` blocks them; a truncated 130 MB partial once got
  swept in by `git add -A`. Check `git diff --cached --name-only` after any bulk add.
- **Call `knitr::include_graphics()` once on a vector**, never inside a `for` loop — the
  looped form returns a value that is never printed and silently renders zero images.
- **Reports use `reports/_setup.R`** with absolute paths. Do not set `knitr` `root.dir`; it
  breaks `include_graphics()` + `embed-resources`.
- Tables go through `write_table_pair()` (`.rds` always, `.csv` per config). Every script
  ends with `write_run_manifest()`.

## Statistical posture

This project has repeatedly found that plausible-looking results were artefacts. Assume the
same of new ones and check before reporting.

- **Never test against a whole-genome or pre-selected background.** The universe is the genes
  quantified and tested in that same experiment. A dataset without a context table is
  skipped, not given a fallback background.
- **Lead with effect sizes.** With thousands of matched isoforms, p-values track n. A
  paired Wilcoxon p of 7e-9 on a median log2 ratio of 0.003 is sample size, not biology.
- **Label statistics conditioned on the outcome.** Anything selecting on significance in both
  datasets is prefixed `selconf_` and is descriptive only.
- **A term count is not a finding.** GO/Reactome nesting means a few genes generate dozens of
  terms; always report `pathway_enrichment_driver_summary` alongside.
- **Check small-n patterns before claiming them.** Several compelling-looking contrasts have
  failed Fisher tests at p ≈ 0.24.
- **Do not condition on a mediator.** `02e`'s Mantel-Haenszel stratification on LPS-response
  magnitude is arithmetically right and causally over-adjusted: UBL5 is a spliceosome-
  associated modifier, so response magnitude is plausibly downstream of the very effect being
  measured. It cannot distinguish the two causal models.

## Retracted — do not cite these numbers

Each was reproducible from its table and still wrong. Reproducing one is a check on the code,
not a licence to quote it. Detail in `docs/external_review/HANDOFF.md` §"Superseded".

- **The layer asymmetry** (64.0% expression vs 37.7% splicing retained) — an artefact of
  comparing two estimators biased in opposite directions. `scripts/10_layer_estimator_calibration.R`
  calibrates it; the corrected gap is −0.04. The mediation split (36.3% / 63.7%) goes with it.
- **"The UT deficit does not survive its confounder"** — withdrawn, mediator over-adjustment
  (above). The deficit is not explained away.
- **ISG/LPS "sparing"** (63% / 48% retained vs 37%) — a ratio-scale artefact. On the absolute
  scale those genes lose *more* response, and ISG status does not survive conditioning on
  expression (p = 0.27). The transcriptome-wide headline stands; the sparing reading does not.
- **The "35% reproducibility ceiling"** — HT and UT are different annotation spaces, not one
  space viewed twice. The Jaccard is real; the ceiling framing is not.
- **The H-vs-T switching contrast** (OR 2.29/2.49) — retire, do not recompute (fact 4 above).
- **AFE percentages from the Ensembl subset** (47.6% coverage) — superseded by full GFF3
  coverage (143/143). Pure promoter switch is **21.7%**, not 10.3%. Also retracted: the TSS
  bimodality and its 2.7 kb antimode cutoff (Hartigan dip p = 0.99, unimodal).

## Conventions

Base R + `dplyr`/`tibble`/`ggplot2`, native pipe `|>`, `.data$` inside dplyr verbs.
`library()` wrapped in `suppressPackageStartupMessages()`; optional packages gated behind
`requireNamespace()`. Progress via `message()`; scripts end with `message("NN_name.R: done")`.

**Numbers in reports are computed inline from the results tables at render time, never typed
as literals.** Regeneration changes every number, and hardcoded prose goes stale silently.

`docs/REVIEW_CHANGES.md` is the permanent record of what was wrong, the evidence, and what
changed. Read it before re-litigating a decision; add to it rather than editing history when
something is superseded.
