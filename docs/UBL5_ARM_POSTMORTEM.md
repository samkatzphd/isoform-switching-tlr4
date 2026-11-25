# Why the UBL5 arm yields almost no candidates

*Attached: `ut_q1_attrition.png` (where the candidates went), `ut_q1_why_it_fails.png` (why
there was never enough signal). Tables: `ut_q1_power_vs_n.csv`, `ut_q1_power_audit.csv.gz`,
`ut_q1_locus_fusion_screen.csv`, `ubl5_sart1_knockdown_stats.csv`,
`ut_response_magnitude_control.csv`, `ut_library_composition.csv`.*

## The headline

The UBL5 arm asks whether LPS drives a different isoform choice in wildtype THP-1 than in the
UBL5 knockdown. Of **13,086 genes tested for that interaction, 2 reached significance, and 1
survived quality control.** That is the whole candidate list.

This is not a thresholding accident that a different cutoff would fix. It is the predictable
result of four limits that compound, plus three design features that would flag the arm even
if the power were there. Both figures are attached; this document is the argument.

**What is not broken.** The wildtype THP-1 data from this experiment are fine — 174 genes
switch on LPS in wildtype, comparable to the main analysis, and the aggregate layer-retention
work is sound and remains interpretable (see "Still worth testing"). What fails specifically
is **gene-level candidate discovery from the genotype × treatment interaction.** Please do not
read this as "the UBL5 experiment was bad"; read it as "this design cannot deliver a
gene-level candidate list, and here is exactly why."

---

## Four compounding limits

### 1. The perturbation is weak

Measured from the sequencing itself, UBL5 mRNA in the knockdown is **67–70% of wildtype** — a
30–33% knockdown with roughly 70% of the transcript still present (unstimulated ratio 0.671,
95% CI 0.483–0.931; LPS 0.696, CI 0.648–0.747).

The reagent is not at fault: the shRNA is on-target, it hits 13 of UBL5's 14 transcript models
plus one unexpressed pseudogene and nothing else in 315,349 transcripts, and the one model
lacking the target site carries 0.2% of UBL5 output. There is no escape isoform. This is
simply how much knockdown the line achieves.

**Protein was never measured.** miR-30-based shRNAs can repress translation beyond transcript
degradation, so protein loss may be deeper. That single missing number governs how everything
downstream should be read, which is why it tops the list of things still worth doing.

### 2. The LPS response the knockdown must act on is itself blunted

Across genes tested in both arms, the knockdown's LPS response runs at roughly two-thirds of
wildtype amplitude: mean |log2FC| 0.451 → 0.288 (×0.64), SD ×0.73, IQR ×0.58, and only 57.9%
of the top-200 response magnitude retained.

This matters mechanically, not just descriptively. Differential transcript usage needs a
condition-driven shift in isoform fraction to clear a threshold. If the whole perturbation
runs at reduced amplitude, fewer shifts clear it **for reasons that have nothing to do with
splicing.** Switching genes fall 118 → 62 (crude OR 0.508, p = 1.2 × 10⁻⁵).

### 3. The interaction test carries a variance penalty by construction

This is the step where the list actually collapses, and it is the least intuitive point, so it
is worth being precise.

A genotype × treatment interaction is a difference of differences — it combines four group
means rather than two, so its standard error is inflated by about √2 relative to a simple
contrast. Measured in this dataset, the median SE ratio of interaction to wildtype-LPS is
**1.35**, against a theoretical 1.41.

The consequence is that fewer than 1% of isoforms had an interaction effect large enough to be
detectable at all: median detection threshold 2.09 in log2 units versus 1.50 for
LPS-in-wildtype, giving 0.82% detectable against 1.99%.

The observed interaction effects do look larger — median |log2 effect| 0.424 versus 0.341 —
but **that is standard-error scaling, not signal.** When estimates are noisier their absolute
values are larger even when nothing is happening, and the two ratios track each other (SE
ratio 1.39, |effect| ratio 1.24). The right comparison is against what pure noise predicts.
Under a true effect of zero, |estimate| is half-normal with median 0.674 × SE:

| contrast | median \|effect\| observed | expected from noise alone | ratio | nominally p < 0.05 |
|---|---|---|---|---|
| LPS in wildtype | 0.341 | 0.317 | **1.075** | **6.37%** |
| LPS in knockdown | 0.277 | 0.277 | 1.000 | 4.66% |
| **genotype × treatment** | **0.424** | **0.442** | **0.959** | **4.05%** |

This is the sharpest statement in the document, and it is stronger than a power argument.
**The interaction contrast shows no excess over the null anywhere in its distribution** — not
merely at the FDR-corrected tail. Only 4.05% of isoforms reach nominal p < 0.05 against the
5% expected if nothing were real, and the median effect sits slightly *below* the noise
expectation. By contrast the wildtype LPS contrast shows a genuine excess (6.37% versus 5%,
median effect 7.5% above noise), so the data are perfectly capable of carrying signal — the
interaction simply does not contain any that this design can see.

The two genes that surfaced are therefore consistent with what 13,086 tests produce by chance.
That is also why the empirical null found nothing (q = 0.20) and why the pipeline's own
diagnostic marks it unusable.

A caution on reading the 4.05%: a slightly conservative figure can also mean the variance
estimates are slightly inflated rather than that the signal is exactly zero. The honest
statement is "no detectable excess", not "proven absence".

### 4. n = 3 is the binding constraint

Holding the observed effects fixed and scaling only the replicate count:

| libraries per arm | % of isoforms detectable |
|---|---|
| 3 | 0.82% |
| 4 | 2.16% |
| 6 | 6.09% |
| 8 | 10.50% |
| 12 | 18.62% |
| 20 | 30.68% |

Going from 3 to 6 flasks would multiply the detectable fraction roughly sevenfold. This
calculation is deliberately **optimistic** — it assumes the observed point estimates are the
true effects, which overstates things because the observed estimates are noisy and the large
ones are inflated by selection. Treat the shape of the curve as the message, not the exact
percentages.

---

## Why n = 3 was not a problem for the HT screen

The obvious objection is that the main (HT) screen ran on the same three-replicate design and
produced a usable candidate list, so n = 3 cannot be the whole story. That objection is
correct, and resolving it is the clearest way to see what went wrong here.

**The HT screen never asked an interaction question.** Its candidates come from the LPS effect
*within* one arm — a difference between two group means — with the second cell system used
afterwards as a direction-agreement check rather than as a formal interaction partner. Script
11's header states this explicitly: "the wildtype LPS switching response, on its own terms."
Only script 12 (this arm) and script 08 fit an interaction anywhere in the project.

Set the two side by side, same pipeline, same FDR, same n = 3:

| arm | question asked | contrast | switching genes | % of tested |
|---|---|---|---|---|
| H_HT, primary MDM | LPS effect within one arm | 2 group means | 61 | 0.49% |
| T_HT, THP-1 anchor | LPS effect within one arm | 2 group means | 75 | 0.61% |
| T_UT, UBL5 wildtype | LPS effect within one arm | 2 group means | **174** | **1.33%** |
| U_UT, UBL5 knockdown | LPS effect within one arm | 2 group means | 37 | 0.28% |
| **UBL5 genotype × treatment** | **difference of differences** | **4 group means** | **2** | **0.015%** |

The wildtype arm of the UBL5 experiment found **more** switching genes than either HT arm —
174 against 75 and 61. The data are not worse. The same libraries, analysed with a marginal
contrast, outperform the screen that produced the project's main candidate list.

Everything collapses at the one step where the question changes. Within a single model, on a
single dataset, moving from "does LPS change isoform usage in wildtype" (174 genes) to "does
LPS change isoform usage *differently* in the knockdown" (2 genes) costs a factor of 87.

So n = 3 was always the binding constraint; the HT screen simply never asked a question
expensive enough for it to bite. Three replicates buy enough precision for a marginal
contrast and not enough for an interaction — and no amount of care in the UBL5 analysis could
have changed that, because the cost is in the contrast, not in the execution.

---

## Three design features that flag the arm independently of power

### The genotype control is untransduced wildtype, not sh-Luciferase

The sequenced comparison is parental THP-1 versus a retrovirally transduced, G418-selected
knockdown line. Every wildtype-versus-knockdown difference therefore carries transduction and
selection alongside UBL5 loss, and those cannot be separated — in the design matrix, "UBL5
knocked down" and "transduced and selected" are the same column. No amount of sequencing depth
helps a confound that is perfectly collinear with the variable of interest.

The interaction is partly protected, because a constant offset cancels in a difference of
differences. Abundance comparisons are not protected at all. The sh-Luc line exists and was
used for the original qPCR validation, so this is fixable at the bench even though it cannot
be fixed in the data.

### Library composition separates perfectly by genotype

TMM normalisation factors show no overlap between arms — all six wildtype libraries ≤ 1.013,
all six knockdown libraries ≥ 1.059 — and the wildtype libraries are more concentrated
(top-50 gene share 18.0–19.8% versus 15.6–15.9%). The difference is present at baseline, so it
is not an LPS effect.

TMM assumes most genes are unchanged. Here that assumption is strained, and 56% of genes are
called differentially expressed, which is a symptom rather than a finding. Within-gene
composition measures such as dIF are invariant to library scaling and are not exposed — which
is a further reason to prefer the interaction framing over abundance.

### The reference is dominated by single-library transcript models

89% of the 315,349 UT transcript models carry `num_samples = 1`, and it is not currently
possible to say what that field counts — reference-derived models are 95.5% singletons against
70.3% for PacBio-derived ones, which suggests the reference annotation is itself one of the
merge inputs. Until that is resolved, the replicate-support covariate that the discovery-bias
mitigation calls for cannot be computed, and neither can a "detected in 2 of 3 knockdown
replicates" criterion.

---

## Three of the four hits were artefacts

The locus-fusion screen asks whether a "switch" is really one gene. It takes the switching
isoform and the isoform in the same locus moving hardest the other way — the pair whose ratio
the effect actually describes — and counts shared exonic bases. Two isoforms of one gene always
share some exonic sequence; **zero shared bases means the comparison is between two different
genes the assembler merged into one locus.**

| locus | shared exonic bases | verdict |
|---|---|---|
| ENSG00000288684 (DDX58 + TOPORS) | 0 | disqualified |
| SCAF8 (+ TIAM2) | 0 | disqualified |
| ENSG00000267174 (+ RAB3D, TMEM205) | 0 | disqualified |
| *ATF5* | 455 | passes |

A 75% artefact rate in a four-item list is itself evidence that the list sits at the noise
floor: at the margin of significance, structural artefacts are over-represented because they
can generate large apparent effects.

---

## The result is not robust to arbitrary analysis choices

The interaction gene count depends on decisions that have no principled answer:

| scale and model | genes significant |
|---|---|
| expected counts, primary (~0 + group) | 4 |
| scaledTPM, primary (~0 + group) | 2 |
| expected counts, + replicate term | 0 |
| scaledTPM, + replicate term | 0 |
| either scale, empirical null | 0 |

The empirical-null q-value for the surviving hits is **0.20**, and the pipeline's own
diagnostic marks the empirical null unusable on this dataset. A finding that moves between 0
and 4 depending on count scale and whether a replicate term is included is a finding to test,
not to report.

---

## The deeper problem: the design cannot answer the question

Beyond power, there is a structural issue worth stating plainly. Adjusting the switching
deficit for response magnitude moves the odds ratio from 0.508 (p = 1.2 × 10⁻⁵) to 0.773
(p = 0.114) — significance disappears.

But response magnitude is plausibly a **mediator**, not a confounder. If losing UBL5 blunts the
LPS response, and a blunted response yields fewer switches, then adjusting for response
magnitude removes part of the very effect being measured. So the non-significant adjusted
result does **not** show there is no splicing-specific role — it shows **this design cannot
separate a splicing-specific effect from a globally blunted response.** A different experiment
is required, not a different statistic.

---

## Still worth testing, and why

### 1. Protein-level UBL5 knockdown — do this first

A western blot on the existing line. It is cheap, it needs no new reagents, and it determines
how every result in this arm should be worded. If protein loss is deep, the arm is a genuine
hypomorph and negative results carry some weight. If protein tracks mRNA at ~30%, then every
null in this arm is uninformative rather than negative, and that is the honest framing for the
whole chapter. No other single measurement changes as much.

### 2. The one surviving candidate: a novel *ATF5* isoform

The reason to spend bench time on a one-gene list is not the statistics — it is what the gene
turned out to be. ATF5 is the mammalian functional ortholog of *C. elegans* ATFS-1, the master
UPR^mt transcription factor; UBL5 is the human homolog of UBL-5, which acts with DVE-1 as the
coregulator enabling ATFS-1-dependent transcription. ATFS-1, DVE-1 and UBL-5 are the three
nuclear components of one pathway. An unbiased transcriptome-wide scan in a UBL5 knockdown
returned UBL5's own direct pathway partner as its only survivor.

The switching species is a novel isoform that only long-read sequencing could have found, and
translating it gives an 88-aa product against 282 aa for full-length ATF5 — retaining the
N-terminal targeting region, losing the C-terminal bZIP, and not predicted to be degraded by
nonsense-mediated decay. Dominant-negative architecture.

The first assay is a 514-bp endpoint RT-PCR that either shows the band or does not. Designs
and falsification criteria are in `docs/ATF5_Q1_EXPERIMENTAL_PLAN.md`. If the band is absent,
the arm is closed cleanly and cheaply.

### 3. The aggregate layer result, which does not need a candidate list

This is the part of the arm that survives, and it is worth emphasising to collaborators
because it is a real finding rather than an absence. After calibrating both estimators against
simulated ground truth, UBL5 knockdown retains **47.6% of the wildtype expression response but
74.7% of the splicing response** — and the splicing-response-per-unit-expression-response slope
differs between genotypes (0.0202 versus 0.0140, interaction p = 4.2 × 10⁻⁴), which a purely
downstream model cannot produce.

Note that this runs **opposite** to the naive spliceosome prediction, and that the correction
is what reversed it. Two estimators were in play. The one used in external review gave raw
values that appeared to show the splicing layer collapsing — 0.377 retained versus 0.640 for
expression — which is exactly the headline a spliceosome story wants; after calibration those
become 0.572 and 0.535, i.e. the layers are degraded about equally. This pipeline's own
estimator never showed a splicing collapse even raw (0.652 splicing versus 0.603 expression)
and after calibration points the other way (0.747 versus 0.476).

So the two estimators disagree on magnitude but agree on sign after correction, and **neither
supports the claim that UBL5 loss preferentially destroys splicing responses.** The arm's
defensible aggregate claim is that UBL5 loss blunts the expression response at least as much
as the splicing response, with a small but detectable splicing-specific component on top. The
raw estimator-A numbers should not be quoted on their own.

### 4. The cryptic-switcher class

Genes whose total output barely moves while isoform composition shifts are a larger share of
surviving switching in the knockdown (29.0%) than in wildtype (19.5%). This is a class-level
statement that does not depend on any individual gene clearing FDR, and it is the cleanest
place to look for a splicing-specific effect because output change cannot explain it.

### 5. What a rescued experiment would need

If the genotype × treatment question is worth answering properly:

1. **n ≥ 6 per arm, ideally 12.** THP-1 flasks are cheap; this is the single largest lever.
2. **sh-Luciferase as the genotype control**, so UBL5 loss is separated from transduction.
3. **Deeper, faster depletion** — a degron or CRISPRi line rather than stable shRNA — which
   also avoids the cassette silencing that may explain the modest knockdown.
4. **Protein confirmation at the time of harvest**, not only after selection.
5. **A targeted readout** rather than a transcriptome-wide scan, if specific loci are the
   question — the genome-wide multiple-testing burden is what makes a 1% detectable fraction
   fatal.

---

## What should not be done

- **Do not loosen the FDR threshold to produce a longer list.** At an empirical-null q of 0.20,
  a longer list is a longer list of noise, and the 75% fusion-artefact rate shows what
  accumulates at the margin.
- **Do not define candidates as "switching in wildtype but not knockdown."** That comparison
  returns 174 versus 37 and looks compelling, but it is a power artefact: the knockdown has a
  smaller response and therefore fewer genes clearing a fixed threshold. Direction and
  significance must be graded together.
- **Do not report wildtype-versus-knockdown abundance differences as UBL5-dependent.** With an
  untransduced control and genotype-separated normalisation factors, those differences are
  between two cell lines, not attributable to UBL5.
