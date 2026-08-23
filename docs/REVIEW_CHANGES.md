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

## 0e. Paired design, and what the layer question can actually answer (2026-08-07)

The design was confirmed paired: `T1_minus`/`T1_plus` are the same biological sample before
and after LPS, likewise T2/T3, U1-U3, H1-H3. `08_layer_decoupling.R` was rebuilt around a
paired permutation null (sign-flips within pairs).

Two implementation points that are easy to get wrong:

- The paired null has 2^3 = 8 sign patterns, collapsing to 4 up to a global flip, so only
  **3 null patterns** exist. That cannot support a per-gene SD, so the z-score standardisation
  used for the unpaired null is unavailable; excess over the null mean is used instead.
- Per-pair TVD is **invariant** to swapping a pair, being already an absolute value. The
  splicing statistic must be the TVD of the mean per-pair dIF vector, not the mean of
  per-pair TVDs, or the null is degenerate and identical to the observed value.

### What is robust

Retention of the wildtype LPS response in the knockout, across all four specifications
(paired/unpaired x z/excess): **expression 0.59-0.66, splicing 0.64-0.67**. Both layers are
reduced, by similar amounts, with splicing if anything slightly better retained. The claim
that splicing is disproportionately lost (16% vs 58%) does not hold under any specification
tried.

Baseline isoform composition shows no constitutive defect: unstimulated KO-vs-WT TVD is
0.095 against 0.142 within replicates (ratio 0.67, i.e. below replicate noise), while
H_HT vs T_HT gives 1.31 -- a positive control that the measure detects real composition
differences.

### What is not robust, and why that matters

The response-matched comparison -- the analysis that would show a splicing-specific effect --
depends on the normalisation, not on the pairing:

| paired | statistic | n matched | splicing ratio KO/WT |
|---|---|---|---|
| yes | z (divide by null SD) | 128 | 0.40 |
| no | z (divide by null SD) | 59 | 0.40 |
| yes | excess (subtract null mean) | 82 | 1.00 |
| no | excess (subtract null mean) | 81 | 0.78 |

Dividing by the null SD produces a deficit; subtracting the null mean does not. With 3 null
patterns the SD is essentially unestimable, so the z-based result is the less trustworthy of
the two -- but the honest conclusion is that **this design cannot resolve whether UBL5 has a
splicing-specific role**. Reporting either number alone would overstate the evidence. Written
out as `layer_specification_sensitivity.csv` so the ambiguity is visible rather than buried.

This supersedes the earlier matched-response result (ratio 0.49, p = 4.1e-9) reported from
the unpaired z-based specification alone.

### The hidden regulatory layer (09_hidden_layer.R)

Independent of the UBL5 question, and robust:

- **Cryptic switchers** -- switching genes whose total output barely moves (|log2FC| < 0.25):
  19.5% of switchers in T_UT, 29.0% U_UT, 26.3% T_HT, 32.0% H_HT. Includes IRAK3 (gene
  -0.077, dIF 0.238), SOCS4, SP1, RAB7B, SORT1, LIG1.
- **Compositional rescue** -- gene output falls while the baseline-dominant isoform holds and
  minor isoforms collapse: 0.7-2.9% of dropping genes. PENK: gene -0.39, dominant +0.21,
  minor -1.11, dominant IF 0.48 -> 0.76. Dominance is defined on unstimulated data only;
  defining it post-LPS selects on the outcome.
- Both survive re-annotation of the same libraries (cryptic 41% vs 0.31% base rate, 134x;
  rescue 18% vs 0.25%, 74x). That is robustness to annotation, not biological replication.

### Transcriptome design

HT = T+H transcripts, UT = T+U transcripts, one PacBio run per sample type. So HT and UT are
different transcript spaces and the earlier "35% technical ceiling" framing was wrong --
part of that discordance is genuine annotation difference. Detection is symmetric within the
UT reference (73.4% vs 70.9% of expression on PacBio-novel transcripts in T and U, no
undetected isoforms in either arm), so the Q2 comparison is not biased by the design.
Notably, 69-78% of expression sits on PacBio-novel transcripts across all four datasets from
only 33-39% of isoforms.

---

## 0f. Why the expression-vs-splicing asymmetry is not supported (2026-08-07)

The second external review (`docs/external_review/HANDOFF.md`) makes the layer asymmetry its
load-bearing claim: the knockout retains **64% of the expression response but only 38% of the
splicing response**, with non-overlapping CIs, and concludes the layers are partially
decoupled. It also survives an ENST-only annotation control, so it is not an artefact of the
PacBio reference.

That claim is **not rejected on the biology. It is rejected on the measurement**, and the
reason is specific and testable.

### The two layers are measured with estimators that are biased in opposite directions

Expression response is `|log2FC|` between condition means of three replicates. Averaging
before taking the absolute value gives high signal-to-noise, but the absolute value still
inflates small values, so the **ratio is pulled toward 1** — deficits are understated.

Splicing response is a total variation distance: a sum of absolute values over isoforms.
Absolute values do not average out noise. For signal small relative to noise,
`E|s + n| ~ E|n| + O(s^2)`, so any noise-subtracted TVD responds roughly **quadratically**
to signal near zero and **compresses ratios toward 0** — deficits are overstated.

Comparing an upward-biased number against a downward-biased one manufactures an asymmetry
even when both layers are reduced by exactly the same factor.

### Calibration against known ground truth

`scripts/10_layer_estimator_calibration.R` simulates a known signal ratio and measures what
each estimator reports. Estimator A is the review's (pairwise across-condition TVD minus
within-condition TVD); B is this pipeline's (TVD of the mean per-pair dIF vector, minus a
paired permutation null).

| true ratio | expression | splicing A | splicing B |
|---|---|---|---|
| 1.0 | 1.00 | 1.01 | 0.99 |
| 0.8 | 0.84 | 0.68 | 0.72 |
| 0.6 | 0.68 | 0.42 | 0.47 |
| 0.4 | 0.56 | 0.20 | 0.24 |
| 0.2 | 0.46 | 0.05 | 0.07 |

Both splicing estimators are compressed; the expression estimator is inflated. Inverting the
calibration on the observed values:

| source | observed expr / splicing | corrected expr / splicing | corrected gap |
|---|---|---|---|
| external review (estimator A) | 0.64 / 0.38 | **0.53 / 0.57** | -0.04 |
| this pipeline (estimator B) | 0.60 / 0.65 | **0.48 / 0.75** | -0.27 |

The review's observed gap of **+0.26 becomes -0.04** after correction — i.e. the two layers
are reduced by indistinguishable amounts, and if anything splicing is retained slightly
better. The two independent analyses, which disagree sharply on the raw numbers, **agree
after correction**.

### What this does and does not establish

- It does **not** show that UBL5 has no splicing-specific role. It shows this measurement
  cannot demonstrate one, in either direction.
- It does **not** overturn the finding that the KO's LPS response is globally blunted; that
  rests on expression alone and is unaffected.
- The mediation decomposition (36% indirect / 64% direct) inherits the same compressed
  splicing estimator on both sides and should not be quoted until recomputed on a calibrated
  scale.
- The simulation uses Gaussian homoscedastic noise, which real count data is not. The
  correction is approximate. The **direction** of the two biases is not approximate — it
  follows from the algebra of averaging before versus after an absolute value.

### What survives from the review and is adopted

- **Baseline composition is normal in the KO** (excess 0.0024, ratio 1.04) against a working
  positive control (H vs T, ratio 1.57). This does not involve a cross-layer comparison and
  reproduces here. It is genuine evidence against a constitutive splicing defect.
- **The withdrawal of the Mantel-Haenszel over-adjustment** is correct and was already acted
  on in section 0e: response magnitude is plausibly a mediator, and conditioning on it removes
  the effect being measured. `02e_response_magnitude_control.R` is retained but its output
  must be read as "cannot separate the models", not as evidence of no effect.
- **The 16% retention figure was a ratio of medians near zero** and the review corrects it to
  ~38%. Correct diagnosis, and the same instability is why script 08 reports means.
- **The "35% technical ceiling" retraction** is correct: HT and UT are genuinely different
  transcript spaces, so part of that discordance is annotation, not noise.
- **Class `j` is not enriched among switching isoforms** relative to the expressed background
  (OR 0.61). This corrects an impression left by the earlier class-code section, which
  compared switching classes to each other rather than to the background.
- **Class B reproduction: IRF3 and BCL7B do not reproduce** across annotations. Any writeup
  naming them must say so.

---

## 0g. The ISA objects reconciled against raw counts, and why `H_HT` is not comparable (2026-08-10)

RSEM matrices for all twelve libraries of each reference arrived at
`/Volumes/Expansion/IsoformSwitchAnalyzer/Counts/{H-T,U-T}/{Genes,Isoforms}/` (expected
counts, TPM and FPKM). Every ISA object was reconciled against them by recomputing isoform
fraction from scratch: IF per library = isoform TPM / sum over the gene's isoforms, then
condition means and `dIF`.

### Three datasets reproduce exactly; one does not

| dataset | dIF correlation | median abs difference | switching genes ISA -> recomputed | Jaccard |
|---|---|---|---|---|
| `T_HT_anchor` | **1.0000** | 1.4e-05 | 95 -> 106 | 0.90 |
| `T_UT` | **1.0000** | 1.4e-05 | 118 -> 112 | 0.95 |
| `U_UT` | **1.0000** | 1.5e-05 | 62 -> 64 | 0.94 |
| `H_HT_supplementary` | **0.777** | 1.9e-02 | 169 -> 35 | 0.17 |

The 1.4e-05 residual is exactly the four-decimal rounding of the stored IF values, so those
three are reproduced to the precision the objects carry. The Jaccard below 1 is threshold
ties at `|dIF| = 0.15`, not disagreement.

Gene-level totals agree for all four (r = 0.9999), including `H_HT`. Whatever differs about
`H_HT` changes how a gene's expression is **apportioned across its isoforms**, not how much
the gene has.

### The cause: a confounder correction applied to one dataset out of four

Three hypotheses were tested and rejected before the right one. It is not a mislabelled or
foreign count file — cross-correlating each ISA `H` column against all twelve RSEM columns,
every column matches its own library best (r ~ 0.987) against 0.75-0.82 for any `T` column.
It is not an annotation-version mismatch — the `gene_id` map is 100% identical and all
128,491 isoform ids are present. It is not a filtered-denominator artefact — restricting the
gene sum to the isoforms retained in the object changes nothing.

The `args` recorded in each `.Rdata` show all four runs read the same RSEM directory and the
same GTF. The design matrices differ:

| run | design columns | ISA version |
|---|---|---|
| `T_HT` | `sampleID, condition` | 2.6.0 |
| `T_UT` | `sampleID, condition` | 2.6.0 |
| `U_UT` | `sampleID, condition` | 2.6.0 |
| `H_HT` | `sampleID, condition, `**`sv1, sv2`** | 2.4.0 |

`IsoformSwitchAnalyzeR::importRdata` runs `sva::num.sv` on every import. When it finds
surrogate variables it applies
`limma::removeBatchEffect(log2(TPM+1), design = condition, covariates = SVs, method = "robust")`,
back-transforms with `2^x - 1`, clips negatives to zero, and **derives gene expression and
every IF/dIF value from the corrected matrix** (source lines ~1170-1205). It found two
surrogate variables for `H_HT` and zero for the other three.

Reapplying that exact transform to the RSEM TPM reproduces the object:

| | raw RSEM vs ISA | corrected RSEM vs ISA |
|---|---|---|
| per-replicate IF | r = 0.978-0.987 | **r = 0.9997-0.9999** |
| `dIF` | r = 0.777 | **r = 0.99915** |
| within-condition replicate IF SD (`H_minus`) | 0.0062 | 0.0030 (object holds 0.0031) |

So the counts are complete and correctly labelled, and all four objects are now verified
against them. The difference is processing, not data.

### Why only `H_HT` has surrogate variables — the sample provenance

Confirmed by the experimentalist 2026-08-10, and it makes the asymmetry expected rather than
accidental:

- **`T`** is the THP-1 **cell line** differentiated to macrophages, and **`U`** the UBL5
  knockout line. `T1`-`T3` are **independent replicates through the whole workflow** —
  matured separately, treated with LPS separately, RNA extracted separately. They are not
  splits of one flask, and they do capture maturation, treatment and extraction variability.
  What they do not span is **genotype** (one clonal line) or **day** (all processed
  together).
- **`H`** is **h**uman **primary** macrophages donated by volunteers. `H1`, `H2`, `H3` are
  **three different donors**, so genotype varies, and because scheduling followed donor
  availability **each was prepped on a different day**.

So all four datasets have real replication. `H` simply carries **two additional variance
components** — donor genotype and prep day — layered on top of the same process variability
`T`/`U` have. `sva` finding two surrogate variables in `H` and none in the others is
**correct behaviour, not a defect**: it is detecting structure that is genuinely present in
one design and genuinely absent from the other.

Note that in `H`, **donor and prep day are perfectly confounded** — one donor per day. A
single blocking term absorbs both, and no analysis of this data can attribute H's between-
replicate variation to genotype rather than batch.

The ISA version difference is **not** the cause: the H *unfiltered* export was built with
2.6.0 and still acquired `sv1`/`sv2`, while `T` under that same 2.6.0 acquired none. The
detection is data-driven.

### What this invalidates

This resolves Tier-2 open question #6. The answer is not "biological outlier" and not
"annotation bias" — it is that two datasets with fundamentally different replicate structures
were processed by one pipeline that adapted to each, and then their outputs were compared as
though they were on one scale. It explains `H_HT`'s whole profile:

- **Lowest replicate IF noise of the four datasets.** The correction halves it
  (0.0062 -> 0.0030). Uncorrected, H's donor variation would make it the *noisiest*, not the
  quietest — the ranking is inverted by the correction.
- **Its recommended expression floor is the lowest of the four.**
  `expression_floor_recommendation.csv` gives `T_UT` 20.16, `U_UT` 11.96, `T_HT` 10.07,
  **`H_HT` 6.02**. `02d` takes the *median* of the four, so H pulls the consensus down only
  modestly (~11 vs ~12.5 if H were on the same footing) — the median of four is set by the
  middle two. Real, but a second-order effect, not a distortion of the floor.
- **2.07% switching genes against 0.87% for `T_HT`.** Less replicate noise, plus `sv1`/`sv2`
  carried in the design into the downstream test, yields more calls at the same threshold.

Consequently the **H-vs-T switching contrast in `02e` (crude OR 2.29, MH 2.49) is not
interpretable as a biological difference** — and, given the provenance above, it would not be
even if both were processed identically. A switching rate is a count of genes whose effect is
consistent across that dataset's replicates, so the two rates answer **different questions**:
in `H`, "consistent across three donors on three days"; in `T`, "consistent across three
process replicates of one genotype on one day". H's is the strictly harder test. Their ratio
is therefore not a genotype effect. Re-scoring `H_HT` on raw-count `dIF` while holding its
ISA q-values fixed drops it from 169 to 35 switching genes and inverts the OR to 0.44; that
number is **not** the corrected answer (the q-values still come from the corrected fit), it
only bounds how much of the contrast rides on the correction.

One real concern for whoever reruns this: `H`'s `sv1` is **not orthogonal to condition**
(`H1_plus` +0.56, `H2_plus` +0.57, but `H3_plus` -0.46, against -0.21 to -0.23 for all three
minus samples). A latent variable partly collinear with the contrast of interest absorbs real
condition signal along with the unwanted variation. `sv2` looks like a cleaner H3-vs-H1/H2
split — i.e. donor. This is the specific hazard of letting `sva` *rediscover* a factor that
is already known.

### What is unaffected

Q1's anchor (`T_HT`) and both Q2 datasets (`T_UT`, `U_UT`) reproduce exactly from raw counts.
Nothing resting on those three is touched by this. `H_HT` is supplementary throughout, and
the ISG result it was already the odd one out on (no enrichment) now has a candidate
explanation.

### The fix

**Not** `detectUnwantedEffects = FALSE`. Given the provenance, that would leave genuine
donor and batch variation in `H` uncorrected — inflating its noise and discarding real power
— and it would not make `H` comparable to `T` in any meaningful sense, because the two
designs are not comparable to begin with.

The right change is that **donor is a known factor, so it belongs in the design explicitly**:
fit `~ donor + condition` for `H` (a paired donor design, which is what the experiment
actually is) rather than letting `sva` rediscover it latently. That removes the same
variation, is reproducible, and eliminates the collinearity hazard `sv1` exhibits.

> **Scope correction (2026-08-23, §0i).** This is not a change to how the upstream script is
> *called*. `NCBR-40-main/data/sample_sheet_isa_*.tsv` carries **only `sampleID` and
> `condition`** — there is no donor, batch, or replicate column anywhere, so the design matrix
> handed to `importRdata` is `~ condition` and sva had no declarable alternative. Implementing
> this needs a new column in the sample sheet **and** a covariate threaded through
> `isa_import()`. Do not plan it as a one-line edit.

Apply the analogous term uniformly — `~ replicate + condition` for `T` and `U`. Their
replicates are independent through maturation, treatment and extraction, and the design is
paired (§0e), so replicate is a legitimate blocking factor there too, not a formality. `sva`
finds nothing in those datasets, so they should barely move; running the same model
everywhere is easier to defend than blocking only the dataset that forced the issue.

Then re-import all four **under one ISA version** — a consistency fix in its own right, since
`H_HT` was built with 2.4.0 and the other three with 2.6.0 — but note that the version is not
what produced the asymmetry, so this will not by itself change H's behaviour.

Everything needed is on the drive: the RSEM matrices, the GTFs named in `args`
(`{H-T,U-T}_Comparisons/reference/GRCh38_Gencode_CHR_v40_plus_*_Isoseq.gtf`) and the
transcript fastas (`transcripts_{H-T,U-T}.fa`). The installed build is 2.4.0, so upgrade
before rerunning.

**What will not be fixed by any of this:** the H-vs-T switching-rate contrast. Correct
per-dataset modelling makes each dataset's own calls trustworthy; it does not make a
three-donor design and a one-day cell-line design yield comparable rates. That comparison
should be retired rather than recomputed.

---

## 0h. The model-based 2x2 refit, and what it says about UBL5 (2026-08-22)

`scripts/12_dtu_refit_2x2.R`. The blocking step from `OPEN_QUESTIONS.md` #1 is done. Because
the U-T RSEM matrix carries **all twelve libraries** (T1-T3 and U1-U3, each ± LPS), the layer
question was fitted as **one `genotype x treatment` interaction** in satuRn rather than as a
comparison of two separately-fit ISA runs — removing the cross-run normalisation problem and
the estimator problem (§0f) together.

### The headline: the interaction is almost empty

| contrast | what it asks | genes called (scaledTPM) | (expected_count) |
|---|---|---|---|
| `WT_LPS` | does LPS remodel isoform usage in wildtype? | **174** | 290 |
| `KO_LPS` | does it in the knockout? | **37** | 47 |
| **`interaction`** | **does the LPS effect DIFFER between genotypes?** | **2** | **4** |

Out of ~13,000 genes tested. The two scaledTPM interaction genes are **ATF5** and
`ENSG00000288684`, both at q = 0.020.

**This is the first direct test of the Q2 claim.** Every previous version compared two
separately-fit arms and read the gap between them; none tested the difference itself. The
174-vs-37 marginal gap looks dramatic and is the source of the "the KO switches less"
framing — but it is the textbook trap of treating a **difference of significance as the
significance of a difference**. Tested directly, almost nothing separates the arms.

**Do not over-read this as proof of no effect.** The interaction is intrinsically the least
powered contrast: it is a difference of differences, so variance adds, and the median
standard error is 0.655 against 0.470 for the WT main effect — roughly 1.4x. At n = 3 per
group, a real but modest interaction would not be detected. The honest reading is that the
data **cannot support** a transcriptome-wide claim that UBL5 loss reshapes LPS-driven isoform
selection, not that it has been ruled out.

### Retention, finally on one scale

Both arms now come from the same model and the same estimator, so a KO/WT ratio is
like-for-like — which §0f's retracted 64%/38% never was. Slope through the origin,
KO = beta x WT, on WT-significant isoforms:

| scale | model | retention (lower bound) | 95% CI |
|---|---|---|---|
| scaledTPM | primary | **0.557** | 0.487-0.627 |
| scaledTPM | + rep3 | 0.551 | 0.464-0.638 |
| expected_count | primary | 0.471 | 0.420-0.521 |
| expected_count | + rep3 | 0.414 | 0.352-0.476 |

It is a **lower bound**: selecting isoforms on WT significance gives the WT estimates a
winner's curse, inflating the denominator and biasing the slope down. The KO estimates for
those isoforms are not selected on and stay unbiased.

So retention is **at least ~56%** on the primary scale — consistent with the old *expression*
figure (~64%) and **not** with the retracted *splicing* figure (38%). On a single estimator
the two layers degrade at similar rates, which is what §0f predicted once the calibration was
inverted and the ordering disappeared.

**The unselected, transcriptome-wide slope is not usable and is kept only as a diagnostic.**
Regressing across all 127,315 isoforms returns 0.068, but that is regression dilution, not
biology: the reliability of the WT estimate — the share of its observed variance that is
signal rather than measurement error — is **0.0066**. About 99% of the transcriptome-wide
spread is noise, because most isoforms have no LPS effect at n = 3. The standard correction
(divide by reliability) is unstable at that magnitude and returns beta > 10. Any future
transcriptome-wide ratio must report its reliability alongside it.

### Two method findings worth carrying

**satuRn's empirical FDR fails at this sample size.** satuRn offers BH on the theoretical
null (`regular_FDR`) and an empirically re-estimated null (`empirical_FDR`, via locfdr).
The empirical version emits `f(z) misfit` warnings on every contrast and never drops below
q ≈ 0.30 — it returns **zero** genes even for the wildtype LPS effect, the least
controversial signal in the dataset, which ISA independently calls for 118 genes. A
correction that finds nothing where the positive control is strongest is broken, not
conservative: with n = 3 there are too few informative z-statistics to estimate a null, so
locfdr inflates it and absorbs the signal. **BH is reported**; both are written to
`dtu_refit_fdr_comparison.csv` so the failure stays on the record.

**The paired design cannot be blocked in this model.** A 6-level pair factor is rank
deficient in the 2x2 (ncol 9, rank 8) because pair nests entirely inside genotype. The
primary model is therefore unblocked `~ 0 + group`; a `~ 0 + group + rep3` sensitivity fit
assumes T1 and U1 were processed as matched batches, which has not been confirmed. The two
agree on every qualitative conclusion (interaction 2 vs 0, retention 0.557 vs 0.551).

### The count-scale decision, recorded

`OPEN_QUESTIONS.md` #1 asked for this to be chosen deliberately. **scaledTPM is primary**
(`config: dtu.count_scale`), because raw expected counts are biased by effective-length
differences between isoforms of the *same* gene — precisely the comparison DTU makes.
`expected_count` is fitted alongside every time as a sensitivity check. It is uniformly more
liberal (290 vs 174 WT genes) but changes no conclusion.

> **Incomplete as written — see §0i item 5.** Reading the upstream script confirms ISA's
> counts are abundance-derived, as claimed above, but they are also **TMM-normalised**
> (`importIsoformExpression(normalizationMethod = 'TMM', calculateCountsFromAbundance = TRUE)`).
> `12_dtu_refit_2x2.R` builds its scaledTPM from raw `colSums(expected_count)` library sizes,
> with no TMM, so the refit's primary scale is **not identical** to ISA's `isoformCountMatrix`.
> This is untested — the drive was unmounted before it could be checked. The conclusions are
> unlikely to move (the interaction is empty by a wide margin and `expected_count` agrees
> qualitatively), but the correspondence should be measured rather than assumed.

### What this does and does not settle

Settled: the retention comparison is no longer estimator-dependent, and the layers degrade
at similar rates. The blunted **expression** response (§0e, slope 0.52) is untouched and
still the best-supported Q2 phenotype.

Not settled: whether UBL5 has a splicing-specific role. The direct test is near-empty but
underpowered, so the verdict stays **unresolved** rather than moving to rejected. What would
resolve it is more replicates or a perturbation with a larger effect — not another estimator.

---

## 0i. The upstream ISA scripts, read against ours (2026-08-23)

`NCBR-40-main/` was added to the project folder: the code that *produced* the four `.Rdata`
objects. It is one CLI script, `scripts/IsoformSwitchAnalyzeR.R`, plus three sample sheets.
This repo starts where that script stops, so the two do not overlap in scope — but the
boundary between them carries assumptions that had never been checked against the source.
They are recorded here.

### The pipeline, end to end

```
NCBR-40-main                                   this repo
  RSEM per-sample isoform quant
  importIsoformExpression(TMM, countsFromAbundance)
  importRdata(gtf, fasta, designMatrix)         01_load_data.R
  preFilter(geneExpr>1, IF>=0.01, ...)          02..07  summaries, QC, comparisons
  isoformSwitchTestDEXSeq(reduce=TRUE)          08..10  layers, hidden layer, calibration
  analyzeORF                                    11,13   candidates, assayability
  [analyzeAlternativeSplicing]  <- NOT run      12      model-based 2x2 refit
  extractTopSwitches -> TSV
  save.image() -> .Rdata                        (reads the .Rdata)
```

One invocation per contrast, driven by `-c1`/`-c2`. Four runs for our datasets, each repeated
with `-p 1.0 -f 0.0 -u 0.0` to produce the unfiltered context exports.

### 1. The dIF threshold does not match, and ours is stricter

Upstream filters at **|dIF| >= 0.1** — both in `preFilter` and in `isoformSwitchTestDEXSeq`,
and stated in the README ("abs(dIF) >= 0.1"). `config/config.yml` uses
`significance.min_abs_dif: 0.15`.

Stricter downstream is safe and nothing needs to change numerically. What needs correcting is
how the pre-reduction caveat is phrased: the reduced objects were reduced against **q < 0.05
AND |dIF| >= 0.1**, not against the gene-level q alone. Genes whose best isoform sits between
0.1 and 0.15 are present in the objects and are discarded by our own scoring. So "every gene
in the object already passes the cutoff" is true of *their* cutoff, not ours, and the two
differ.

### 2. `geneExpressionCutoff = 1` is hardcoded, and it applies to the context layer too

`isa_import()` calls `preFilter` with `geneExpressionCutoff = 1` and
`isoformExpressionCutoff = 0` as literals. Neither is exposed as a CLI argument, so the
`-p 1.0 -f 0.0 -u 0.0` "unfiltered" runs **still have them applied**.

Consequence for the rule that pathway enrichment uses the context layer as its universe: that
universe is *genes quantified with expression above 1*, not every quantified gene. The cutoff
is low enough that no conclusion is likely to move, but the claim as written is slightly
stronger than the data supports and should be qualified.

### 3. Why sva had no alternative — the sample sheets have no batch column

`importRdata` is called with no `detectUnwantedEffects` argument, so it takes the package
default and sva runs. That confirms §0g's mechanism. The sample sheets add the reason:

```
sampleID        condition
H1_minus_S21    H_minus
H1_plus_S22     H_plus
...
```

**`sampleID` and `condition` are the only columns.** There is nowhere to declare donor, prep
day, or replicate. H's donor structure could only ever have been discovered latently, because
the design matrix passed to `importRdata` is `~ condition` and nothing else.

This enlarges §0g's fix. "Refit with an explicit donor blocking factor" is not a change to how
the script is *called* — it needs a new column in the sample sheet **and** a covariate threaded
through `isa_import()` into the design matrix. Recorded so nobody plans that work as a
one-line edit.

The same sheets independently confirm that `T_HT` and `T_UT` are the same six libraries:
`T1_minus_S15` through `T3_plus_S20` appear identically in `sample_sheet_isa_H-T.tsv` and
`sample_sheet_isa_U-T.tsv`. That was previously inferred from filenames; it is now confirmed
from the inputs the runs actually consumed.

### 4. Each contrast was fitted on six libraries, separately

The design matrix is subset to `condition_1`/`condition_2` before import, so every run sees
exactly 6 libraries. `T_UT` and `U_UT` were never in one model, which is what
`12_dtu_refit_2x2.R` changed — the refit is doing something new, not repeating existing work.

It also means **H's two surrogate variables were estimated from 6 samples**. That is a small
basis for a latent-variable correction, and an additional reason to prefer an explicit term.

### 5. The count scale: §0h is right in mechanism, incomplete in detail

`importIsoformExpression(calculateCountsFromAbundance = TRUE, normalizationMethod = 'TMM')`
confirms §0h's claim that ISA's `isoformCountMatrix` is abundance-derived rather than RSEM
expected counts.

But it is also **TMM-normalised**, and `12_dtu_refit_2x2.R` builds its scaledTPM as
`TPM * colSums(expected_count) / 1e6` — raw library sums, no TMM. So the refit's primary scale
is *not* identical to ISA's counts, and §0h implies a closer correspondence than holds.

**Unverified:** the drive was unmounted before this could be tested. The check to run when it
is back is a per-library correlation of `isa_list$isoformCountMatrix` against both the raw and
TMM-scaled versions. The refit's conclusions are unlikely to move — the interaction is empty
by a wide margin and `expected_count` already agrees qualitatively — but the discrepancy is
real and should be closed rather than assumed away.

### 6. `analyzeORF` ran; `analyzeAlternativeSplicing` did not

**ORF analysis was run on every dataset**, and its output reached our committed files.
`isoformFeatures_T_UT.rds` carries `PTC` (192 TRUE, 1,465 FALSE, 112 NA) and `IR` (293
isoforms with at least one retained intron). This is the annotation Tier-2 open question #4
asks for — "are the lost minor isoforms NMD targets or truncated, and the retained dominant
ones full-length coding" — and it is already on disk. **No drive needed.**

**Alternative-splicing analysis was NOT run.** `analyzeAlternativeSplicing` is gated behind
`--run_extra_analysis`, and none of the four commands in the README pass that flag. Two
consequences: the GFF3 event-structure work in the external review was necessary rather than
duplicative, and ISA's native ATSS classification — an independent check on the 21.7%
pure-promoter figure — is one flag and one rerun away. That one does need the drive.

### 7. Reproducibility posture differs sharply

Upstream: `module load R/4.4`, `packages.R` installs whatever BiocManager currently serves with
no version pin, and `save.image()` dumps the entire workspace. The absence of pinning explains
the ISA version inconsistency §0g found across the archived objects (2.4.0 for the `H_HT`
filtered run, 2.6.0 elsewhere) — the runs happened at different times and picked up whatever
was current. The `save.image()` habit is why `args` was recoverable from inside the `.Rdata`
at all, which is what made §0g diagnosable.

Downstream records git SHA, R version, key package versions and thresholds per script in
`results/run_manifest.json`.

### What agrees, and is worth stating

- **Direction convention.** `dIF = condition_2 - condition_1` with `c1 = X_minus`, so positive
  dIF means increased usage after LPS. Downstream treats it the same way. No sign inversion.
- **The unfiltered exports really are unreduced.** `reduceToSwitchingGenes = TRUE` still runs,
  but at `alpha = 1.0` and `dIFcutoff = 0.0` it retains everything (subject to item 2).
- **q-value columns.** `score_isoforms()` reads `isoform_switch_q_value` /
  `gene_switch_q_value`, the same DEXSeq outputs `extractTopSwitches` sorts on upstream.

### Not ours

`data/sample_sheet_isa_NCBR-390_H-T.tsv` belongs to a different project: 48 samples in groups
`grpA`-`grpP`, patient/healthy, sharing only the H-T transcriptome. It is not part of this
analysis and should not be mistaken for a fifth dataset.

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
