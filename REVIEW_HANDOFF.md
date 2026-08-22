# Review handoff: isoform-switch-pipeline

> **SUPERSEDED — historical record, 2026-08-07.** This is the first review round. It is kept
> unedited (the append-don't-edit rule below applies to it too), but do not act on it
> directly. The current document is `docs/external_review/HANDOFF.md`.
>
> **Finding 1 below is withdrawn.** It concluded that the Q2 / UBL5 switching deficit "does
> not survive its confounder". Response magnitude is plausibly a **mediator**, not a
> confounder — UBL5 is a spliceosome-associated modifier — and conditioning on a mediator
> removes the effect being measured. The MH arithmetic was right, the causal logic was not.
> `02e` was built on this recommendation and its output is kept as a record only.
>
> **Finding 2's "technical ceiling" framing is retired.** HT and UT merge different PacBio
> transcript sets, so part of the 0.35 Jaccard is annotation-space difference, not noise.
>
> Full ledger of what has since been retracted: `docs/AGENT_CONTEXT.md`
> §"Retracted — do not cite", and `docs/external_review/HANDOFF.md` §"Superseded".

External review, 2026-08-07. Every number below was recomputed from the committed
`data/processed/isoformContext_*.rds` tables using the repo's own `score_isoforms()` and
`config/config.yml`. No repo files were modified (grant was read-only).

Read `docs/REVIEW_CHANGES.md` before acting on any of this. Nothing here re-litigates a
decision recorded there; these are new findings. Follow the project's append-don't-edit
rule and record the outcome as a new section `§0d`, not by editing history.

## Conventions any new code must follow

- Base R + `dplyr`/`tibble`/`ggplot2`, native pipe `|>`, `.data$` inside dplyr verbs.
- `library()` wrapped in `suppressPackageStartupMessages()`; optional packages behind
  `requireNamespace()`.
- Reuse `score_isoforms()` from `utils/helper_functions.R`. **Do not write a local copy** —
  that is the invariant that caused the drift documented in §3.1.
- Significance and presence must come from the same FDR universe: use
  `load_scoring_table()`, never mix the reduced object with the context layer.
- Thresholds live in `config/config.yml` only. Add new knobs there, read via
  `analysis_param()`.
- Tables via `write_table_pair()`; every script ends with `write_run_manifest()` and
  `message("NN_name.R: done")`.
- Scripts take an optional config path as `argv[1]`.
- Reports compute numbers inline from results tables at render time. Never type a literal.
- `knitr::include_graphics()` once on a vector, never in a loop.
- New scripts here work from committed `data/processed/` and do **not** need
  `/Volumes/Expansion`.

---

## Finding 1 (highest priority): the Q2 / UBL5 conclusion does not survive its confounder

`02d` rules out a *baseline expression* confounder and does so correctly. But the relevant
confounder for a post-LPS switching claim is the **magnitude of the LPS response itself**,
which was never tested. It is large.

Gene-level LPS response, 8,595 genes tested in both UT arms (from `gene_value_1`/
`gene_value_2`, verified constant within `gene_id`):

| metric | WT (T_UT) | KO (U_UT) | ratio KO/WT |
|---|---|---|---|
| SD of log2FC | 0.6827 | 0.4990 | 0.731 |
| IQR of log2FC | 0.6077 | 0.3508 | 0.577 |
| mean abs log2FC | 0.4507 | 0.2883 | 0.640 |
| regression slope KO~WT | — | — | **0.516** |
| Spearman WT vs KO | — | — | 0.727 |

Top-200 LPS responders in WT retain 57.9% of their response magnitude in the KO. Direction
is preserved (rho 0.73); amplitude is roughly halved.

Isoform switching requires a condition-driven shift in isoform fraction against a fixed
`min_abs_dif: 0.15`. If the whole perturbation runs at half amplitude, fewer shifts clear
the threshold irrespective of any splicing-specific role for UBL5.

Stratifying each gene by its own LPS response magnitude, measured in its own dataset:

| abs log2FC bin | n WT | n KO | switching WT | switching KO | % WT | % KO |
|---|---|---|---|---|---|---|
| <0.1 | 1566 | 2537 | 8 | 8 | 0.51 | 0.32 |
| 0.1-0.25 | 2039 | 2877 | 12 | 7 | 0.59 | 0.24 |
| 0.25-0.5 | 2327 | 2049 | 22 | 9 | 0.95 | 0.44 |
| 0.5-1 | 1859 | 775 | 29 | 14 | 1.56 | 1.81 |
| >1 | 804 | 357 | 28 | 15 | 3.48 | 4.20 |

- Crude odds ratio, KO vs WT switching: **0.532**
- Mantel-Haenszel OR adjusted for response stratum: **0.784, 95% CI 0.558-1.100, p = 0.173**

The deficit is not significant once response magnitude is matched, and in the two
highest-response strata the KO switches slightly more. Note also that the KO is depleted of
high responders (357 vs 804 genes with abs log2FC > 1), which is the mechanism.

Supporting negative control already checked: replicate IF noise is *not* worse in the KO.
Matched on gene-expression strata the KO is marginally quieter than WT (noise ratio
0.92-0.96 across bins below 250), so the deficit is not a KO-noise artefact either.

### What to build

Add `scripts/02e_response_magnitude_control.R`:

1. Read the context tables for `T_UT` and `U_UT` via `load_context_table()`.
2. Compute per-gene LPS log2 fold-change as `log2((gene_value_2 + 1) / (gene_value_1 + 1))`,
   taking one row per `gene_id` (values are constant within gene — assert this), filtered to
   `gene_value_1 + gene_value_2 >= significance.min_gene_expression`.
3. Join per-gene switching status from `score_isoforms()` on the context table.
4. Write `results/tables/ut_response_magnitude_strata.csv` (the stratified table above) and
   `results/tables/ut_response_magnitude_control.csv` (crude OR, MH OR, CI, p, plus the
   SD/IQR/slope summary).
5. Add a `response_strata` grid to `config/config.yml` under `analysis:` rather than
   hardcoding the bin edges.
6. Same treatment for the HT pair (`T_HT` vs `H_HT`) so the two arms are handled alike.

### What to change in the reports

`reports/seminar_summary.qmd` and `reports/ut_t_vs_u_comparison.qmd` currently frame reduced
switching in the KO as evidence for a UBL5 role in post-LPS isoform selection. Revise to the
supported claim: **UBL5 loss blunts the LPS response globally by roughly half, and reduced
isoform switching follows from that.** That is still a real and interesting phenotype — it is
a different claim. All numbers must be read inline from the new tables.

`ut_T_vs_U_background_enrichment` (OR 36.1) is unaffected — it measures co-switching among
genes tested in both, not the deficit. Leave it as is.

---

## Finding 2: an unused same-sample control that sets the technical ceiling

`T_HT` and `T_UT` are the **same six libraries** — replicate IF column names are
`T1_minus_S15, T1_plus_S16, T2_minus_S17, T2_plus_S18, T3_minus_S19, T3_plus_S20` in both —
quantified against two reference transcriptomes. Nothing in the code, reports, or docs notes
this. It is the best technical control the project has: zero biological variability, so all
disagreement is annotation plus thresholding.

Over 10,592 gene symbols testable in both:

- switching in `T_HT`: 95; in `T_UT`: 116; in **both: 55**; Jaccard **0.353**
- expected overlap under independence: 1.04 (so agreement is far above chance, OR 233)
- gene-level max abs dIF Spearman across all shared genes: **0.771**

Of the 101 discordant genes, status in the reference that did *not* call it:

| reason | n |
|---|---|
| near-miss, abs dIF 0.10-0.15, above expression floor | 53 |
| below the expression floor | 24 |
| clearly absent (abs dIF < 0.10, above floor) | 24 |

Median abs dIF in the non-calling reference is 0.132 against a 0.15 cutoff; 72.3% of
discordant genes reach at least 0.10. So the disagreement is threshold brittleness, not
contradiction — which is reassuring, but it means a single-dataset switching call carries
roughly a 1-in-3 chance of not replicating under re-annotation.

### What to build

Add `scripts/02f_reference_concordance.R` producing
`results/tables/reference_concordance_T_HT_vs_T_UT.csv` with the counts above and the
discordance decomposition. Drive the sample-identity check off the replicate IF column names
rather than hardcoding it, and have the script assert the libraries actually match before
reporting concordance.

### What to change in the reports

State the ~35% same-sample agreement wherever cross-dataset replication is reported, as the
technical ceiling. This reframes an existing headline: the `T_HT` to `H_HT` replication rate
of 6.1% should be judged against ~35%, not against 100%. Same for the UT one-way concordance
table.

---

## Finding 3: ISG replication is partly one dataset counted twice

`07_isg_analysis.R` line ~154 applies `p.adjust(p_value, "BH")` across all 12 dataset x set
rows as if independent. `T_HT` and `T_UT` are the same libraries (Finding 2), and their
switching ISGs overlap 7 of 9 union:

- `T_HT` (8): CCL3, CCL4, GBP1, IFIT1, IFITM2, RAB7B, SHFL, USP18
- `T_UT` (8): CCL22, CCL3, CCL4, GBP1, IFITM2, RAB7B, SHFL, USP18
- `U_UT` (4): CCL1, CCL22, CCL3, CCL4
- `H_HT` (3): IFITM2, JAK2, USP18

`seminar_summary.qmd` describes enrichment "in the wildtype THP-1 datasets on both
references" as replication. It is one biological result observed twice.

The enrichment itself looks real and is the project's strongest positive: OR ~6-8 against a
properly tested background, carried by type II, with `U_UT` showing independent type II
enrichment (OR 9.8). But n is small (8, 8, 4, 3 switching ISGs).

### What to change

1. Correct within reference-transcriptome families rather than across all 12 rows; keep the
   across-all column too if useful, but rename so the two cannot be confused (the project
   already does this well with `selconf_`).
2. Describe wildtype HT/UT agreement as technical reproducibility, not biological
   replication.
3. Report the switching ISG counts alongside every OR so the small n is visible.

---

## Smaller items

1. **A paired design appears to be analysed as unpaired.** Sample names (`T1_minus`/`T1_plus`)
   indicate three paired biological replicates. `02d`'s noise model pools within-condition SD
   without using pairing, and nothing passes a blocking factor to the test. If the design is
   genuinely paired, a paired analysis would materially increase power at n=3 — which matters
   given the floor removes 33-59% of calls. Confirm the design with the experimentalist
   first; if paired, either note the power cost explicitly in the report or re-test with a
   blocking factor. Do not change the model on the basis of filename inference alone.
2. **Switching calls are tail events; say so.** At gene expression >= 12 and mid-IF
   (0.05-0.95), the 99th percentile of abs dIF is 0.147 (T_HT), 0.153 (T_UT), 0.134 (U_UT),
   0.201 (H_HT). The 0.15 cutoff sits near the 99th percentile of the entire distribution, so
   calls come from the extreme tail where estimator error dominates. That is consistent with
   the 35% cross-reference agreement. The threshold choice is defensible; the reports should
   state that switching is tail behaviour, not typical.
3. **`H_HT` is an outlier and needs an annotation-bias check.** It has 2.07% switching genes
   vs 0.87% for `T_HT`, 3.11% of isoforms over abs dIF 0.15 vs 0.91%, the *lowest* replicate
   noise (median IF SD 0.0241), and loses only 14.2% of calls to the expression floor vs
   35.7% / 48.4% / 58.7% elsewhere. More calls plus less noise plus fewer low-expression
   calls may be real, but H is a different genotype on a reference built partly from it. Check
   for annotation bias before using it as supplementary support.
4. **Pathway enrichment: drop term counts from the headline.** `top5_share_of_gene_hits` is
   75-95%; 82 terms in the shared UT set come from 10 genes. The `driver_summary` table
   already prevents misreading — go further and report driver genes only, omitting term
   counts from prose entirely.

---

## Suggested order of work

1. `02e_response_magnitude_control.R` and the Q2 reframing (changes a conclusion).
2. `02f_reference_concordance.R` and the replication-ceiling caveat (changes how several
   existing numbers are read).
3. ISG correction scope in `07` and the seminar wording.
4. Confirm the pairing question with the experimentalist.
5. `H_HT` annotation-bias check.
6. Pathway prose cleanup.
7. Append all of it to `docs/REVIEW_CHANGES.md` as `§0d`, with evidence, in the existing
   what-was-wrong / evidence / what-changed format.

## Acceptance checks

- New scripts run from committed `data/processed/` without the external drive mounted.
- `Rscript scripts/99_verify_inputs.R` still exits 0 with the drive mounted.
- No new copy of the switching rule: `grep -rn "score_isoforms <- function" scripts/` returns
  nothing.
- Gene summary row counts still equal distinct `gene_id` counts (`02`/`02b` assertions pass).
- Reproduce these values exactly before changing prose: MH OR 0.784 (CI 0.558-1.100,
  p = 0.173); same-sample Jaccard 0.353 with 55 shared genes; LPS slope 0.516.
- Every regenerated report renders with images embedded (HT report had 0 once — regression
  risk).
- `git diff --cached --name-only` after any bulk add: no `.Rdata`.

## Caveat on this review's own method

Per-gene LPS fold-change was derived from ISA's condition means (`gene_value_1`/
`gene_value_2`), not a re-fit differential expression model. The blunting is large enough
that a proper DE model is unlikely to reverse it, but confirming Finding 1 against the count
matrices would make the argument airtight, and is worth doing before the claim goes into a
manuscript.
