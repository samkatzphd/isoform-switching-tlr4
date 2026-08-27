# External review — findings, corrections, and work to do

External analysis of the isoform-switch-pipeline results, 2026-08-07, revised 2026-08-18.

> **Revision 2026-08-18 — read this before anything else.** Four items in this document are
> now superseded; they are listed in "Superseded — do not implement" as items 3-6, and the
> `FIGURES.md` walkthrough (`figures/F1`-`F5`) reflects the corrected state and is the fastest
> way in.
>
> **The most important one is item 6.** The expression-vs-splicing layer asymmetry — the
> load-bearing claim of the UBL5 argument in Finding 2 — is **retracted**. Claude Code's
> `scripts/10_layer_estimator_calibration.R` (commit 36690a2) showed it rests on comparing two
> estimators biased in opposite directions; I reproduced the calibration, tested whether
> matching the estimator structure rescues it (it does not — the layers differ 20x in
> signal-to-noise), and the rebuttal holds. The Q2 conclusion is **unresolved** pending a
> formal differential-splicing model, not supported.
>
> Items 3-5 follow from the GFF3 and transcript FASTA now being available locally (`gff3/`,
> `fa/` — both gitignored, 3.9 GB, do not commit), which allowed the event-structure analysis
> to run at 100% coverage rather than the 47.6% reachable via Ensembl: the ISG-sparing claim
> is withdrawn (a ratio artefact), the AFE percentages are replaced, and a reported bimodality
> in TSS distance turned out to be an artefact of mixing two annotations.

**This document supersedes the earlier `REVIEW_HANDOFF.md`.** If that file is in the repo,
read the "Superseded" section below before implementing anything from it — its headline
recommendation was based on an analysis I subsequently retracted.

Every number here was recomputed from the committed `data/processed/isoformContext_*.rds`
tables using the repo's own `score_isoforms()` and `config/config.yml`. No repo files were
modified. Each figure and claim traces to a CSV in `tables/`.

---

## Contents

```
HANDOFF.md                       this file — findings and work to do
FIGURES.md                       narrated walkthrough of F1-F5; start here
EXPERIMENTS.md                   proposed wet-lab follow-up, with primer designs
literature_review.md             cited literature grounding both questions

figures/
  F1_scope.png                   what the data is; reproducibility; cryptic switching
  F2_ubl5_mechanism.png          baseline composition + the retracted layer asymmetry (item 6)
  F3_specificity.png             what the deficit is NOT (three negative controls)
  F4_event_structure.png         what the switches structurally are (GFF3, 100% coverage)
  F5_implications.png            assay feasibility and a verdict per question

  fig_gff3_event_structure.png   full event classification + nesting detail
  fig_ubl5_mechanism.png         earlier version — panels b/c superseded, see item 3
  fig_ubl5_causality.png         annotation control, mediation, pathway concentration
  fig_hidden_layer.png           cryptic switchers and Class B
  fig_network.png                STRING distance to TLR4, degree-matched null
  fig_sparing_reanalysis.png     the ISG-sparing reanalysis (see item 3)
  fig_class_codes.png            class-code composition of switching isoforms
  fig_event_structure.png        SUPERSEDED — Ensembl-only, 47.6% coverage
  fig_tss_bimodality.png         SUPERSEDED — the retracted bimodality (item 5)

tables/                          one CSV per claim (see "Which table backs which claim")
gene_sets/                       gene symbol lists + the tested universe, for enrichment
```

---

## Conventions any new code must follow

- Base R + `dplyr`/`tibble`/`ggplot2`, native pipe `|>`, `.data$` inside dplyr verbs.
- `library()` wrapped in `suppressPackageStartupMessages()`; optional packages behind
  `requireNamespace()`.
- Reuse `score_isoforms()` from `utils/helper_functions.R`. **Do not write a local copy** —
  that is the invariant behind the drift documented in `REVIEW_CHANGES.md` §3.1.
- Significance and presence must come from one FDR universe: use `load_scoring_table()`.
- Thresholds live in `config/config.yml` only; read via `analysis_param()`.
- Tables via `write_table_pair()`; scripts end with `write_run_manifest()` and
  `message("NN_name.R: done")`.
- Scripts take an optional config path as `argv[1]`.
- Reports compute numbers inline from results tables at render time. Never type a literal.
- All analyses below run from committed `data/processed/` and do **not** need
  `/Volumes/Expansion`.
- Record the outcome in `docs/REVIEW_CHANGES.md` as a new section (`§0d`), appending rather
  than editing history.

---

## Superseded — do not implement

An earlier version of this review concluded that the UBL5 knockout's reduced isoform
switching was **explained away** by its blunted LPS expression response, based on a
Mantel-Haenszel analysis stratifying on response magnitude (crude OR 0.53 → MH OR 0.784,
95% CI 0.558-1.100, p = 0.173). That recommended building
`02e_response_magnitude_control.R` to implement the stratification and reframing the Q2
conclusion around it.

**That reasoning was wrong and the recommendation is withdrawn.** UBL5/Hub1 is a
spliceosome-associated modifier. If splicing regulation is upstream of the LPS response,
then LPS-response magnitude is a **mediator**, not a confounder — and conditioning on a
mediator removes the very effect being measured. The stratification was arithmetically
correct and causally over-adjusted; it cannot distinguish the two models, because both
predict the same correlation.

The MH numbers are retained in `tables/review_findings.csv` and
`tables/review_switching_by_response_stratum.csv` for the record, but must not be
presented as evidence that the splicing deficit is an artefact.

Two further corrections to that earlier document:

1. **A retention figure was misstated.** It reported the KO retaining "16%" of the wildtype
   splicing response. That was a ratio of medians, and the median splicing excess sits near
   zero, making the ratio numerically unstable. The correct mean-based figure is **~38%**.
   The qualitative conclusion (splicing degrades much more than expression) holds; the
   magnitude does not. `tables/ubl5_retention_by_layer.csv` carries the corrected values.
2. **The "35% same-sample reproducibility ceiling" framing was overstated.** It assumed
   `T_HT` and `T_UT` were two annotations of the same transcript space. Per the
   experimental design (below), they are not: HT merges T+H PacBio transcripts and UT
   merges T+U, so part of that disagreement is genuine annotation-space difference rather
   than measurement noise. The decomposition of the discordance still holds and is still
   useful; the "ceiling" language should not be used.

### Added 2026-08-18

3. **The ISG-sparing claim is withdrawn.** Finding 3 below reports that interferon-response
   and LPS-responsive genes retain *more* splicing response in the KO than the transcriptome
   average (63% and 48% vs 37%), and reads this as those pathways being spared. That is an
   artefact of expressing everything as a ratio. On the **absolute** scale the same genes
   lose *more* response, not less (ISG 0.0121, LPS-responsive 0.0102, transcriptome 0.0081 —
   `tables/sparing_relative_vs_absolute.csv`). Conditioning on expression, LPS-responsiveness
   survives (p = 0.002) but **ISG status does not** (p = 0.27) —
   `tables/sparing_conditioned_on_expression.csv`. Report both scales or neither; do not
   claim an interferon-specific effect. `figures/F3_specificity.png` panels a-b show this
   directly; `fig_sparing_reanalysis.png` is the full version.

4. **The AFE percentages are replaced.** An earlier pass classified switching events using
   Ensembl exon structures for the `=`-class isoforms only — 68 of 143 events, 47.6%. With
   `gff3/U_T.gff3` every event is resolvable (143/143, and 78/78 in the KO). The corrected
   numbers, from `tables/event_structure_summary_gff3.csv`:

   | | Ensembl subset | **full GFF3** |
   |---|---|---|
   | naive "first exon differs" | 97.1% | **98.6%** |
   | true alternative promoter (disjoint, >1 kb) | 51.5% | **69.2%** |
   | **pure promoter switch (no internal change)** | 10.3% | **21.7%** |
   | involves internal junction change | 69.1% | **59.4%** |

   Four-way classification of all 143 WT events: 47.6% promoter change **plus** internal
   splicing, 21.7% pure promoter switch, 18.9% other terminal, 11.9% internal splicing only.
   So Robinson's alternative-first-exon concern applies to ~22% of events, not ~50%, and the
   majority still involve genuine internal junction changes. Per-event flags are in
   `tables/event_structure_WT_gff3.csv` / `..._KO_gff3.csv`.

   The two Ensembl-derived tables (`event_structure_WT.csv`, `..._KO.csv`, and the `_modes`
   variants) are kept for the record only. **Use the `_gff3` tables.**

5. **A reported bimodality in TSS distance is retracted.** On the Ensembl subset the
   distribution of TSS distance between switch partners appeared bimodal (Hartigan dip
   p = 0.008, modes at ~85 bp and ~37 kb, antimode 2.7 kb), and a 2.7 kb cutoff was proposed
   as a data-derived threshold. On complete annotation the distribution is **unimodal**
   (dip p = 0.99, median 12,679 bp). The cause is not the added events — those have identical
   mode composition (Fisher p = 1.0) — but disagreement between the two annotations:
   Spearman 0.649 on TSS distance, 14 of 68 events differing more than 5-fold, and the same
   switch partner selected in only 57% of cases. The merged reference redefines 5' ends from
   PacBio reads and contains isoforms Ensembl lacks.

   **Do not implement the 2.7 kb threshold.** Use `>1 kb` with disjoint first exons as the
   alternative-promoter criterion, which is annotation-independent. Also retracted: a
   claimed association between UBL5-dependence and promoter proximity (Fisher p = 0.31 once
   classified by mode rather than tested on the continuous scale), and any claim that UBL5
   loss preferentially removes internal-splicing events (59% vs 51%, Fisher p = 0.26).

6. **The expression-vs-splicing layer asymmetry is retracted — Findings 2b and 2c.** This was
   the load-bearing claim of the whole UBL5 argument: 64% of the expression response retained
   in the KO against 38% of the splicing response, non-overlapping CIs, robust to an ENST-only
   annotation control. `scripts/10_layer_estimator_calibration.R` (commit 36690a2) rejects it
   on the measurement, and **the rebuttal is correct**. I reran the script and reproduced the
   calibration table exactly.

   The mechanism: expression response is |log2FC| between condition means, so averaging three
   replicates *before* the absolute value tracks truth almost linearly. Splicing response is a
   total variation distance — a sum of absolute values, which does not average out noise, so
   for signal small relative to noise E|s+n| ≈ E|n| + O(s²). The splicing estimator is
   therefore roughly **quadratic** near zero and compresses ratios toward 0, while the
   expression estimator does not. Comparing them manufactures an asymmetry even when both
   layers are reduced by the same factor:

   | true ratio | expression reports | splicing reports |
   |---|---|---|
   | 1.0 | 1.00 | 1.01 |
   | 0.8 | 0.84 | 0.68 |
   | 0.6 | 0.68 | 0.42 |
   | 0.4 | 0.56 | 0.20 |

   Inverting on the observed values turns **64% / 38% into 53% / 57%** — the observed gap of
   **+0.26 becomes −0.04**. The sign of the effect therefore depends entirely on the
   estimator, which is sufficient to retract the claim.

   **The two corrected analyses agree only on direction, not magnitude.** Do not repeat the
   "indistinguishable amounts" framing from commit 36690a2's message — the numbers beside it
   do not support it:

   | source | observed | corrected | corrected gap |
   |---|---|---|---|
   | external review (estimator A) | 0.640 / 0.377 | 0.535 / 0.572 | −0.038 |
   | this pipeline (estimator B) | 0.603 / 0.652 | 0.476 / 0.747 | −0.271 |

   Both corrected gaps are ≤ 0, so neither supports "splicing degrades more than expression".
   But they differ by **23 points**, and estimator B implies splicing is retained ~1.6x
   *better* than expression while estimator A implies near-parity. The defensible conclusion
   is that **after correction neither estimator yields a reliable cross-layer ordering** — the
   two do not converge on a shared answer, they merely both fail to support the original one.
   Backing: `tables/layer_estimator_calibration.csv`, `tables/layer_retention_corrected.csv`.

   **I checked whether matching the estimator structure rescues it. It does not.** Applying an
   identical across-minus-within pairwise form to *both* layers narrows the gap from 26 to 14
   points (expression 53.7%, splicing 39.9%) but cannot close it, because compression depends
   on **signal-to-noise**, not only estimator form, and the layers differ there by **20x**:

   | layer | arm | signal | noise | SNR |
   |---|---|---|---|---|
   | expression | WT | 0.3053 | 0.1686 | 1.81 |
   | expression | KO | 0.1638 | 0.1537 | 1.07 |
   | splicing | WT | 0.0127 | 0.1407 | 0.09 |
   | splicing | KO | 0.0051 | 0.1328 | 0.04 |

   At SNR 0.09 the splicing estimator is deep in the quadratic regime while expression is not.
   Matched structure is necessary but not sufficient — so this is not fixable by tweaking the
   distance measure.
   Backing: `tables/layer_matched_estimator.csv`, `tables/layer_snr.csv`.

   **Scope of the retraction.** It does *not* show UBL5 has no splicing-specific role; it
   shows this measurement cannot demonstrate one in either direction. Finding 2a survives
   untouched (within-layer, no cross-layer ratio). The globally blunted LPS expression
   response survives (expression alone). **Finding 2c must not be quoted** — the mediation
   decomposition inherits the same compressed splicing estimator.

   **What to build instead.** A formal differential-splicing model — DRIMSeq or satuRn, with
   the pairing included if the design is paired — fitted on the count matrices, returning a
   per-gene effect size on a scale comparable between layers. Until that exists, the Q2
   conclusion is **unresolved**, not supported and not refuted. Any future claim of the form
   "layer X is affected more than layer Y" must carry an estimator calibration alongside it.

---

## Experimental design facts that constrain interpretation

Supplied by the experimentalist; these are not inferable from the files and change how
several results read.

- **One PacBio full-length run per sample type** (one H, one T, one U) — budget-limited.
- **HT reference** = transcripts observed in T and H, merged. Used for the LPS-response
  question (Q1).
- **UT reference** = transcripts observed in T and U, merged. Used for the WT/KO UBL5
  question (Q2), deliberately including U so that KO-specific transcripts can be detected.

Consequences: HT and UT are **different transcript spaces**, not interchangeable
annotations. Cross-reference comparisons measure annotation sensitivity, not pure technical
noise. Within the UT pair, both arms are quantified against the same merged reference, so
the WT-vs-KO comparison is the clean one.

- **`T_HT` and `T_UT` are the same six libraries** (`T1_minus_S15` … `T3_plus_S20`),
  quantified against the two different references. Nothing in the code or reports notes
  this. It matters for any analysis treating them as independent (see Finding 4).
- **Sample names imply a paired design** (`T1_minus`/`T1_plus`). `02d`'s noise model pools
  within-condition SD without using pairing and no blocking factor is passed to any test.
  **Confirm the design with the experimentalist before changing any model** — this is
  filename inference, not established fact. If paired, a paired analysis would materially
  increase power at n=3, which matters given the floor removes 33-59% of calls.

---

## Finding 1 — The reference-sensitivity audit passes

Novel PacBio transcripts carry most of the expression in every arm (69-78%), so annotation
choice is consequential. But detection between the two UT arms is close to symmetric: of
154,360 shared isoforms, 0.29% are expressed only in U and 0.22% only in T (asymmetry
1.33x). No meaningful bias toward either genotype's own discoveries.

Decisive control — rerun the key result using **ENST-only transcripts**, excluding every
PacBio-novel transcript so both arms use identical reference annotation:

| transcript set | expression retained in KO | splicing retained in KO |
|---|---|---|
| all transcripts | 64.0% (62.5-65.6) | 37.7% (32.7-43.1) |
| ENST-only (annotation-neutral) | 63.9% (62.3-65.5) | 35.0% (25.3-46.2) |

The layer asymmetry is annotation-independent. Backing: `tables/ubl5_retention_by_layer.csv`,
`tables/reference_detection_symmetry.csv`. Figure: `fig_ubl5_causality.png` panel a.

---

## Finding 2 — The UBL5 result, restated correctly

> **2b and 2c are RETRACTED — see item 6 in "Superseded".** The estimator calibration in
> `scripts/10_layer_estimator_calibration.R` shows the expression-vs-splicing asymmetry is an
> artefact of comparing two estimators biased in opposite directions. I reproduced that
> calibration and confirmed it. **2a stands** — it is a within-layer comparison with no
> cross-layer ratio. Do not cite 2b or 2c.

Three constraints, each from an independent analysis:

**2a. Baseline isoform composition in the KO is normal.** Comparing unstimulated cells only
(per-gene total variation distance over shared isoforms, renormalised, calibrated against
within-genotype replicate pairs):

| contrast | within-genotype (noise) | between-genotype | excess | ratio |
|---|---|---|---|---|
| UBL5 KO vs WT | 0.1295 | 0.1343 | 0.0024 | 1.04 |
| H vs T (positive control) | 0.1087 | 0.1711 | 0.0496 | 1.57 |

The KO is indistinguishable from WT at rest; a genuine genotype difference (H vs T) shows
~21x more divergence by the same measure. **No constitutive splicing defect.**
Backing: `tables/baseline_composition_test.csv`. Figure: `fig_ubl5_mechanism.png` panel a.

**2b. Splicing response degrades about twice as much as expression response.** See the
Finding 1 table: 38% of splicing response retained vs 64% of expression response,
non-overlapping CIs, and the same in the annotation-neutral subset. Under a pure
"splicing-is-downstream" model the two layers should degrade together. They do not.

**2c. Most of the splicing deficit is not explained by the expression deficit.** Mediation
decomposition with 1,000 gene-level bootstrap resamples:

| component | share of total KO effect on splicing | 95% CI |
|---|---|---|
| explained by reduced expression response | 36.3% | 29.3-44.5% |
| not explained (direct) | 63.7% | 55.5-70.7% |

Backing: `tables/mediation_decomposition.csv`. Figure: `fig_ubl5_causality.png` panel b.

**Causal caveat, which must survive into any report.** "Direct" means *not explained by the
measured mediator*. Reverse causation — splicing failure degrading the expression response —
fits this decomposition identically. With two conditions and n=3 there is no design-based
way to break the symmetry. What the data establish is that the two layers are **partially
decoupled**, not which is upstream.

---

## Finding 3 — The deficit is transcriptome-wide, not LPS-pathway-specific

> **Partly superseded — see item 3 in "Superseded".** The headline conclusion (the deficit
> is transcriptome-wide, not LPS-pathway-specific) stands and is strengthened. The
> "spared / retain more" reading of the table below does **not**: it is an artefact of the
> ratio scale. On the absolute scale these genes lose more response, and the ISG-specific
> effect vanishes once expression is controlled for (p = 0.27).

If UBL5 were required specifically for an LPS-induced splicing program, the deficit should
concentrate in LPS-responsive genes. It does the opposite:

| gene class | mean splicing excess WT | KO | retained |
|---|---|---|---|
| not LPS-responsive | 0.00958 | 0.00258 | 26.9% |
| LPS-responsive (\|log2FC\|>0.5) | 0.01973 | 0.00949 | 48.1% |
| all other genes | 0.01271 | 0.00465 | 36.6% |
| interferon-response genes | 0.03265 | 0.02050 | 62.8% |

LPS-responsive genes retain **more** than the transcriptome average; ISGs retain the most.
The genotype x responsiveness interaction is significant but in the direction of protection,
not selective vulnerability. Backing: `tables/splicing_deficit_by_gene_class.csv`. Figure:
`fig_ubl5_causality.png` panel c.

**Supported model:** UBL5 loss produces a broad, condition-dependent reduction in splicing
responsiveness — not a constitutive defect (2a), and not one targeted at the inflammatory
program (Finding 3). A "splicing regulation gates the LPS response" framing is **not**
supported, because the genes executing that response are the least affected.

### What to build

`scripts/02e_splicing_response_control.R` — note the name change from the withdrawn
`02e_response_magnitude_control.R`:

1. Read context tables and replicate IF tables via `load_context_table()`.
2. Compute per-gene isoform-composition TVD: within-condition replicate pairs (noise) and
   across-condition pairs (LPS response); the response measure is `lps - noise`. Same
   function applied to unstimulated-only replicates gives the baseline contrast.
3. Write `results/tables/ubl5_retention_by_layer.csv`,
   `results/tables/baseline_composition_test.csv`,
   `results/tables/mediation_decomposition.csv`,
   `results/tables/splicing_deficit_by_gene_class.csv`.
4. Repeat the ENST-only control (filter on `oId` not matching `^PB`) as a column, not a
   separate script.
5. Put bin edges and the bootstrap count in `config/config.yml` under `analysis:`.
6. **Report means, not medians**, for the splicing-excess ratios — the median sits near zero
   and the ratio is unstable there. This is the mistake corrected above.

### What to change in the reports

`reports/seminar_summary.qmd` and `reports/ut_t_vs_u_comparison.qmd` currently frame reduced
switching in the KO as evidence for a UBL5 role in post-LPS isoform selection. Replace with
the three-constraint statement above, computed inline from the new tables. State the causal
caveat explicitly rather than in a footnote.

`ut_T_vs_U_background_enrichment` (OR 36.1) is unaffected — it measures co-switching among
genes tested in both, not the deficit. Leave it.

---

## Finding 4 — Cross-reference comparison and ISG correction scope

`T_HT` and `T_UT` are the same libraries on the two references. Over 10,592 gene symbols
testable in both: 95 switching in T_HT, 116 in T_UT, **55 in both** (Jaccard 0.353), against
1.04 expected under independence. Gene-level max \|dIF\| Spearman across all shared genes is
0.771.

Decomposition of the 101 discordant genes, by status in the reference that did not call them:

| reason | n |
|---|---|
| near-miss, \|dIF\| 0.10-0.15, above expression floor | 53 |
| below the expression floor | 24 |
| clearly absent (\|dIF\| < 0.10, above floor) | 24 |

Median \|dIF\| in the non-calling reference is 0.132 against a 0.15 cutoff; 72.3% reach at
least 0.10. Threshold brittleness dominates. Given the design, this is a mix of annotation
sensitivity and threshold effects — report it as such, **not** as a technical reproducibility
ceiling.

**ISG correction scope.** `07_isg_analysis.R` (~line 154) applies `p.adjust(p_value, "BH")`
across all 12 dataset x set rows as if independent. T_HT and T_UT are the same RNA, and their
switching ISGs overlap 7 of 9 union:

- `T_HT` (8): CCL3, CCL4, GBP1, IFIT1, IFITM2, RAB7B, SHFL, USP18
- `T_UT` (8): CCL22, CCL3, CCL4, GBP1, IFITM2, RAB7B, SHFL, USP18
- `U_UT` (4): CCL1, CCL22, CCL3, CCL4
- `H_HT` (3): IFITM2, JAK2, USP18

Correct within reference-transcriptome families; describe wildtype HT/UT agreement as
cross-annotation reproducibility on shared RNA, not biological replication. Report the
switching-ISG counts alongside every OR so the small n stays visible. The enrichment itself
(OR ~6-8, carried by type II, with independent type II enrichment in U_UT) looks real and is
the project's strongest positive.

---

## Finding 5 — A hidden regulatory layer that gene-level DE misses

Both classes below were checked for reproduction against the independent T_HT annotation of
the same libraries. Figure: `fig_hidden_layer.png`.

**Cryptic switchers.** Of 118 switching genes in T_UT (expression >= 12), **23 (19.5%)** have
\|gene log2FC\| < 0.25 — median 0.136 while shifting isoform fraction by a median of 0.171.
One in five switching genes is effectively invisible to differential expression. The list is
biologically pointed: IRAK3 (IRAK-M, negative TLR regulator), SOCS4, SP1, RAB7B (TLR4
trafficking), LIG1. **64.7% are called switching again** in the independent annotation,
against a 59.8% baseline for all switchers — so they are not preferentially artefacts.
Backing: `tables/cryptic_switchers_T_UT.csv` (carries a `reproduces_in_T_HT` flag per gene).

**Class B — "gene down, dominant isoform preserved".** Dominance is defined at **baseline**
to avoid selecting on the outcome; defining it post-LPS inflates the class. Of 2,611 genes
falling >0.3 log2, **56 (2.1%)** have their baseline-dominant isoform hold flat or rise while
minor isoforms collapse (median: gene -0.362, dominant -0.04, minor -0.532; dominant IF rises
0.25 → 0.32).

Reproduction, with denominators stated explicitly because they differ: of those 56 genes,
**47 are testable in T_HT** (9 have no T_HT counterpart) and **11 reproduce as Class B**
there. That is **23.4% of the 47 testable** (or 19.6% of all 56, if untestable genes are
counted as failures — the conservative reading). Against a 0.52% base rate of Class B among
all testable genes, 23.4% is a **45x enrichment, Fisher OR 78, p = 3e-16**; on the
conservative 19.6% denominator the enrichment is 37.5x. The rate quoted elsewhere in this
document is the 23.4% figure.

The 11 reproducing genes are **AMT, CD163L1, CKAP5, COQ6, DEDD2, FBXL19, HHLA3, NRF1, RBL2,
RTEL1, VPS26C**. Backing: `tables/class_b_purification_T_UT.csv`.

> **Correction on the record.** An earlier draft of this review named IRF3 and BCL7B as
> reproducing Class B genes and singled out IRF3 as the illustrative example. Both have
> `reproduces_in_T_HT = FALSE`. IRF3 remains an interesting single-reference observation and
> must be labelled as not reproducing across annotations. Given that HT and UT are genuinely
> different transcript spaces, non-reproduction is weaker evidence against it than it would
> be for a pure annotation swap — but it is a failed check, not a pass.

**Structural characterisation.** Class B genes are significantly more isoform-complex than
background (median 19.5 vs 13 isoforms, Wilcoxon p = 1.6e-5); cryptic switchers are not
(p = 0.53). That is mechanistically sensible: "gene falls while the dominant isoform holds"
requires enough minor isoforms to absorb the decline. Class B looks like a property of
complex loci, not a distinct functional program. Backing:
`tables/candidate_structural_profile.csv`.

Class B at 2.1% is a real but minority phenomenon. It matters for interpreting specific
genes, not as a claim that DE is broadly misleading.

### What to build

`scripts/03b_hidden_layer.R` writing `results/tables/cryptic_switchers_<label>.csv` and
`results/tables/class_b_purification_<label>.csv` for every dataset with a context table,
each carrying a cross-annotation reproduction flag where a second reference exists. Define
dominance at baseline; assert this in a comment, since the post-LPS definition is the
tempting bug.

### Open work — functional annotation

GO enrichment on these candidate sets could **not** be run here: Bioconductor mirrors are
unreachable from the analysis sandbox (`GenomeInfoDbData` failed against every release), so
`clusterProfiler` would not load. At these set sizes (11-56 genes) GO would likely be
underpowered anyway.

`gene_sets/` contains everything needed to run it where Bioconductor is available:
`universe_tested_genes.txt` (8,916 symbols, the T_UT tested universe) plus
`set_cryptic.txt`, `set_class_b.txt`, `set_class_b_repro.txt`, `set_all_switching.txt`.
Feed these through `05_pathway_enrichment.R`'s existing machinery so the universe convention
matches the rest of the project. Report `pathway_enrichment_driver_summary` alongside any
term counts, per the project's standing rule.

A more informative alternative given the set sizes: annotate the **isoforms** rather than the
genes — are the lost minor isoforms NMD-targets or truncated, and the retained dominant ones
full-length coding? That would convert "purification" from a pattern into a mechanism, and
needs CDS/NMD annotation the pipeline does not currently carry.

---

## Finding 6 — Class `j`: the premise is inverted

Class `j` (novel junction combination) is **not** enriched among switching isoforms — it is
slightly **depleted**:

| class | % of expressed isoforms | % of switching isoforms | Fisher OR | p |
|---|---|---|---|---|
| `=` matches reference | 63.08 | 69.23 | 1.32 | 0.14 |
| `c` contained in reference | 14.99 | 17.48 | 1.20 | 0.41 |
| `j` novel junction combination | 19.09 | 12.59 | **0.61** | 0.055 |

Where the impression comes from: `j` is 46% of all novel (PacBio) transcripts, and novel
transcripts are ~50% of switching isoforms, so `j` is the largest *novel* class among
switchers (18 of 71). But novel isoforms overall are not significantly enriched relative to
their abundance (OR 1.37, p = 0.062). The annotation is `j`-rich; the switching calls are not
`j`-rich beyond that. Backing: `tables/class_code_enrichment_T_UT.csv`. Figure:
`fig_class_codes.png`.

One marginal directional signal, **reported as marginal and not to be built on**: among WT
switching isoforms, reference (`=`) isoforms are predominantly lost on LPS while novel
(`c`/`j`) isoforms are gained (OR 2.13, p = 0.045). It fails to reproduce in T_HT (p = 0.13)
and reverses in the KO (OR 0.24, p = 0.013).

Worth adding to `03_novel_isoform_analysis.R`: the class-code composition of switching
isoforms *against the expressed background*, since the raw composition invites exactly this
misreading.

---

## Smaller items

1. **Switching calls are tail events; say so.** At gene expression >= 12 and mid-IF
   (0.05-0.95), the 99th percentile of \|dIF\| is 0.147 (T_HT), 0.153 (T_UT), 0.134 (U_UT),
   0.201 (H_HT). The 0.15 cutoff sits near the 99th percentile of the whole distribution, so
   calls come from the extreme tail where estimator error dominates. The threshold is
   defensible; the reports should state that switching is tail behaviour, not typical.
2. **`H_HT` is an outlier and needs an annotation-bias check.** 2.07% switching genes vs
   0.87% for T_HT, 3.11% of isoforms over \|dIF\| 0.15 vs 0.91%, the *lowest* replicate noise
   (median IF SD 0.0241), and it loses only 14.2% of calls to the expression floor vs
   35.7/48.4/58.7% elsewhere. More calls plus less noise plus fewer low-expression calls may
   be real, but H is a different genotype on a reference built partly from it.
3. **Pathway enrichment: drop term counts from prose.** `top5_share_of_gene_hits` is 75-95%;
   82 terms in the shared UT set come from 10 genes. Report driver genes only.
4. **The KO is not noisier.** Matched on gene-expression strata, replicate IF noise in U is
   marginally *lower* than T (ratio 0.92-0.96 below expression 250). The splicing deficit is
   not a KO-noise artefact. Worth stating, since it is the obvious first objection.

---

## Which table backs which claim

| table | claim |
|---|---|
| `ubl5_retention_by_layer.csv` | 64% expression vs 38% splicing retained; ENST-only control |
| `baseline_composition_test.csv` | KO baseline composition normal; H vs T positive control |
| `mediation_decomposition.csv` | 36% indirect / 64% direct, with bootstrap CIs |
| `splicing_deficit_by_gene_class.csv` | deficit is transcriptome-wide, spares ISGs |
| `splicing_response_by_expression_stratum.csv` | deficit persists at matched expression response |
| `reference_detection_symmetry.csv` | UT arms detect shared transcripts symmetrically |
| `cryptic_switchers_T_UT.csv` | 23 genes invisible to DE, with reproduction flags |
| `class_b_purification_T_UT.csv` | 56 Class B genes, dominant/minor fold-changes, reproduction flags |
| `candidate_structural_profile.csv` | Class B genes are isoform-complex; cryptic are not |
| `class_code_enrichment_T_UT.csv` | class `j` depleted among switchers |
| `review_findings.csv` | superseded MH analysis + concordance decomposition (record only) |
| `review_switching_by_response_stratum.csv` | superseded MH strata (record only) |

Note `candidate_structural_profile.csv` has `NA` in the switching-isoform columns for the
"class B (reproducing)" row: none of those 11 genes is itself a called switching gene, so
those columns are undefined rather than zero.

---

## Suggested order of work

1. `02e_splicing_response_control.R` + the Q2 reframing (Findings 1-3). Changes a conclusion.
2. `03b_hidden_layer.R` (Finding 5). New result, well supported.
3. ISG correction scope in `07` and the seminar wording (Finding 4).
4. Class-code background comparison in `03` (Finding 6).
5. Confirm the pairing question with the experimentalist.
6. `H_HT` annotation-bias check; pathway prose cleanup.
7. Append to `docs/REVIEW_CHANGES.md` as `§0d` in the what-was-wrong / evidence /
   what-changed format, including the two corrections recorded above.

## Acceptance checks

- New scripts run from committed `data/processed/` without the external drive mounted.
- `Rscript scripts/99_verify_inputs.R` still exits 0 with the drive mounted.
- No new copy of the switching rule: `grep -rn "score_isoforms <- function" scripts/` returns
  nothing.
- Gene summary row counts still equal distinct `gene_id` counts (`02`/`02b` assertions pass).
- Reproduce before changing prose: baseline excess 0.0024 (KO) vs 0.0496 (H vs T); Class B
  reproduction 11 of 47 testable (23.4%) vs 0.52% base rate, 11 of 56 total; class `j` OR 0.61.
  (The retention figures 64.0% / 37.7% and the mediation split 36.3% / 63.7% are still
  reproducible from their tables, but are **retracted** — item 6. Reproducing them is a check
  on the code, not a licence to quote them.)
- Item 6 checks: `Rscript scripts/10_layer_estimator_calibration.R` reproduces the calibration
  grid (true 0.6 → expression 0.68, splicing_A 0.42) and the corrected gap −0.04 for the
  external-review row. Matched-structure control: expression 53.7% / splicing 39.9% with WT
  SNR 1.81 vs 0.09.
- Added 2026-08-18. Event structure, from `gff3/U_T.gff3` at 143/143 WT events: pure promoter
  switch 21.7%, internal junction change 59.4%, true alternative promoter 69.2%, median TSS
  distance 12,679 bp, Hartigan dip p = 0.99 (unimodal). Absolute-scale sparing: ISG loss
  0.0121 vs transcriptome 0.0081. Nesting: 5 of 7 candidate isoforms nested, 6 of 7
  unassayable by qPCR, IRAK3 `TCONS_00065141` amplified by 1 of 315,349 transcripts.
- No stale threshold: `grep -rn "2700\|2\.7 kb" scripts/ reports/` returns nothing (the
  retracted antimode cutoff — item 5).
- `gff3/` and `fa/` stay untracked: `git check-ignore -q gff3/U_T.gff3 fa/transcripts_U-T.fa`
  exits 0. These are 3.9 GB and must never enter the repo.
- Every regenerated report renders with images embedded (the HT report had 0 once).
- `git diff --cached --name-only` after any bulk add: no `.Rdata`.

## Method caveat on this review

All analyses use ISA's condition means (`gene_value_1`/`gene_value_2`, verified constant
within `gene_id`) and replicate isoform fractions — not a re-fit differential expression or
differential splicing model. The layer-asymmetry result is the load-bearing claim and
deserves confirmation against the count matrices before it anchors a manuscript. The
composition measure (total variation distance over renormalised isoform fractions, calibrated
against within-condition replicate pairs) is a reasonable summary but is not a formal
differential-splicing test; a purpose-built method (for example DRIMSeq or satuRn with the
pairing included) would put the central claim on firmer statistical ground.
