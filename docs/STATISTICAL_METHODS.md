# Statistical methods, and why each one

Written for a working scientist who reads and produces quantitative results but has no formal
statistics training. Each section says what the method does, why it was chosen over the
obvious alternative, and how it can mislead. Where this project has actually been misled, that
is stated — those are the most useful parts.

Companion documents: `OPEN_QUESTIONS.md` (what is still unresolved),
`REVIEW_CHANGES.md` (the audit trail).

---

## Part 1 — What we measure

### 1.1 Isoform fraction (IF) and dIF

A gene makes several transcript isoforms. **Isoform fraction** is the share of that gene's
output carried by one isoform: if a gene is expressed at 100 units and one isoform accounts
for 30, its IF is 0.30. **dIF** is the change in that share between conditions — LPS-minus
versus LPS-plus.

**Why fractions rather than isoform expression.** A gene can double its total output with
every isoform rising proportionally. Isoform-level expression would show large changes;
nothing about *which* transcript the cell chooses has changed. Fractions isolate the choice
from the volume, which is the question this project asks.

**The cost of that choice.** A fraction is a ratio, and ratios estimated from small
denominators are unstable. If a gene is barely expressed, its isoform fractions bounce around
between replicates for reasons that have nothing to do with the condition. This is the single
most important practical consequence in the whole project, and §4.1 explains what we do about
it.

**Compositional constraint.** Isoform fractions of a gene sum to 1, so the dIF values sum to
zero. One isoform cannot rise without another falling. This is why per-gene splicing change is
summarised as a **total variation distance** (below) rather than by adding up dIF values,
which would always give zero.

### 1.2 Total variation distance (TVD)

A single number for "how much did this gene's isoform usage change": half the sum of the
absolute dIF values across the gene's isoforms. It runs from 0 (no change) to 1 (complete
switch to different isoforms).

**Why halved.** Because every unit of fraction that leaves one isoform arrives at another, the
raw sum double-counts. Halving makes TVD interpretable as "the fraction of the gene's output
that moved between isoforms".

**The trap, which cost this project two rounds of analysis.** TVD is built from **absolute
values**, and absolute values do not average away noise — they accumulate it. Every
measurement error, positive or negative, adds to the total. So TVD is always biased upward,
and the bias is larger when measurements are noisier. §4.4 covers the consequences.

### 1.3 log2 fold-change

The standard expression measure: `log2(after / before)`. Log scale so that a doubling
(+1) and a halving (−1) are symmetric, which raw ratios are not (2× versus 0.5×).

A small constant is added before dividing (`log2((after+1)/(before+1))`) so that genes with
near-zero expression do not produce infinite or wildly unstable values.

### 1.4 Odds ratio

Used whenever we ask "is X over-represented among Y". If 5% of interferon genes switch and
0.8% of other genes switch, the odds ratio is roughly 6.5 — interferon genes have ~6.5× the
odds of switching.

**Why odds rather than a simple ratio of percentages.** Odds ratios have properties that make
them work correctly with the exact test below, and they behave sensibly when the outcome is
rare. For rare outcomes the odds ratio and the risk ratio are numerically close anyway.

**How to read one.** OR = 1 means no association. The confidence interval matters more than
the point estimate: an OR of 15 with an interval spanning 0.9 to 250 tells you almost nothing.
This project reports intervals for every odds ratio for that reason.

### 1.5 Jaccard index

Overlap between two sets: the size of the intersection divided by the size of the union. If
95 genes switch in one dataset, 116 in another, and 55 in both, the Jaccard is
55 / (95 + 116 − 55) = 0.35.

**Why not just "55 of 95 replicated"** — that denominator ignores the genes found only in the
second dataset, and the number changes depending on which dataset you call the reference.
Jaccard is symmetric.

---

## Part 2 — Deciding that something happened

### 2.1 p-values and q-values (FDR)

A **p-value** answers: if nothing were really going on, how often would I see a result this
extreme by chance? A p of 0.01 means once in a hundred.

The problem is scale. Testing 12,000 genes at p < 0.05 gives about 600 false positives even
if nothing is happening anywhere. The fix is the **false discovery rate**, and the
Benjamini–Hochberg procedure converts p-values into **q-values**. A q of 0.05 means: among all
the results you are calling significant, about 5% are expected to be wrong.

**Why FDR rather than the stricter Bonferroni correction.** Bonferroni controls the chance of
*even one* false positive, which is the right goal when a single false claim is catastrophic.
For genome-scale discovery it is far too conservative — it would discard most real findings.
FDR accepts a controlled proportion of errors in exchange for finding things.

**The critical detail.** A q-value depends on *the whole set of tests it was computed with*.
Take the same gene, the same p-value, and correct it against 12,000 tests instead of 250, and
the q-value changes. This is not a technicality — it caused a real error in this project (§4.2).

### 2.2 Fisher's exact test

Used for every "is X enriched among Y" question: build a 2×2 table (switching / not, in the
gene set / not) and compute the probability of seeing an association this strong by chance.

**Why "exact" and why not chi-squared.** The chi-squared test relies on an approximation that
degrades when any cell of the table has few observations. Several tables here have single-digit
cells — 3 switching interferon genes in one dataset. Fisher's test computes the probability
directly and stays valid at any count.

**What it cannot fix.** The test assumes your 2×2 table describes the right population. If the
background is wrong, the test is precisely and confidently wrong. See §3.1 — this is the error
that produced an odds ratio of 15 from a background of 25 genes.

### 2.3 Wilcoxon rank-sum and signed-rank tests

Used to compare distributions of effect sizes — for example, is |dIF| larger among interferon
genes than others.

**Why not a t-test.** A t-test assumes roughly normal distributions. Effect-size distributions
here are heavily skewed: most genes near zero, a long tail of large values. The Wilcoxon test
uses only the *ranks*, so it makes no distributional assumption.

**The paired version (signed-rank)** is used when the same gene is measured in two conditions
or genotypes, which is more sensitive than treating the two groups as unrelated.

**The trap at large n.** With thousands of genes, the Wilcoxon test detects differences far too
small to matter. This project has a paired Wilcoxon with p = 7 × 10⁻⁹ on a median difference of
0.003 — a 0.2% difference, statistically overwhelming and biologically meaningless. **Report
the effect size; the p-value at this scale is a statement about sample size.** Every summary
table here puts medians before p-values for that reason.

### 2.4 Binomial test

Used for direction agreement: of 105 isoforms changing in one genotype, 86 change the same way
in the other — is that more than the 50% expected by coin flip?

Simple, appropriate, and its main limitation is worth stating: it tells you the direction is
shared, not that the magnitude is comparable. Direction can agree perfectly while one genotype
shows one-fifth the effect.

### 2.5 Test for two proportions

Used to compare rates directly — 24.9% of tested genes switch in wildtype versus 16.4% in the
knockout. Reported with a confidence interval on the *difference*, which is the interpretable
quantity.

### 2.6 Mantel–Haenszel stratified odds ratio

Combines evidence across strata while holding a third variable fixed. Used to ask: if we
compare wildtype and knockout genes *matched on how strongly each responds to LPS*, is there
still a switching difference?

**Why not just adjust in a regression.** Same idea; the MH estimator is transparent — you can
see each stratum's table — which matters when strata are small.

**The serious trap, and this project fell into it.** Adjusting for a variable is only valid if
that variable is a **confounder** (something that influences both the cause and the effect),
not a **mediator** (something on the causal path between them). If losing UBL5 blunts the LPS
response, and a blunted response produces less switching, then response magnitude is a
mediator. Adjusting for it removes the very effect you are trying to measure, guaranteeing a
null result.

The arithmetic was correct; the causal logic was not. Both models predict the same
correlation, so no amount of stratification distinguishes them. See §4.3.

---

## Part 3 — What we compare against

This is where most of the real errors have been, and it gets the most space.

### 3.1 The background (universe) for enrichment

Every enrichment test needs a comparison population. Get it wrong and the test is
meaningless regardless of how carefully it is computed.

**Three candidate backgrounds and what each implies:**

*All human genes.* Wrong here. It treats genes never expressed in these cells as candidates
that failed to switch, which inflates every enrichment. Any immune gene set will look
"enriched" simply because immune genes are expressed in macrophages and olfactory receptors
are not.

*The genes in the saved analysis object.* Also wrong, and this was a live error. The
IsoformSwitchAnalyzeR objects were saved with `reduceToSwitchingGenes = TRUE`, meaning every
gene inside them had already been called significant. Testing "are switching genes enriched
for X" against a background of already-significant genes compares a significant set to a
significant set. One test in this project had a background of **25 genes**, all significant by
construction, and reported an odds ratio of 15.

*Genes actually quantified and tested in that experiment.* Correct, and what the pipeline now
uses throughout — 8,000–12,000 genes per dataset, recovered from unfiltered exports of the
same analyses.

**The rule adopted:** a dataset with no proper background is **skipped**, never given a
fallback. Silently substituting a whole-genome background is worse than reporting nothing.

### 3.2 Permutation nulls

Sometimes there is no standard test for the quantity you care about — "how much did this
gene's isoform composition shift" has no textbook distribution. A permutation null builds the
comparison from your own data: shuffle the condition labels, recompute, and repeat. That tells
you how large the statistic gets when nothing real is going on.

**Why exhaustive rather than random sampling.** With 6 samples in two groups of three there
are only 10 distinct ways to split them, so every one can be enumerated. No random sampling, no
seed, and the result is identical every run. Reproducibility for free whenever the design is
small enough.

**The averaging trap.** The null must be built the *same way* as the observed statistic. A
first attempt here compared condition means (each averaging 3 replicates) against
replicate-versus-replicate differences (no averaging). Averaging reduces noise by roughly √3,
so the null was structurally noisier than the signal — the "null" came out **larger** than the
observed effect (0.146 vs 0.104), and corrections built on it reversed the sign of the answer.
The rule: whatever averaging the real statistic uses, the null must use it too.

### 3.3 Paired designs

In this project each biological sample was measured before and after LPS — `T1_minus` and
`T1_plus` are the same cells. That pairing is information, and using it removes
sample-to-sample variation from the comparison.

**How the null changes.** For paired data the correct shuffle is to swap the two conditions
*within* each pair, not to reshuffle all samples freely. With 3 pairs there are 2³ = 8 such
patterns, collapsing to 4 once you account for a global flip having no effect on a magnitude,
leaving **3 null patterns**.

**What that costs.** Three null values cannot support a reliable standard deviation. So the
paired analysis reports **excess over the null mean** (observed minus average of the nulls)
rather than a z-score (which would divide by a standard deviation estimated from three
numbers). This matters: the choice between those two normalisations changed a headline result
in this project from "splicing halved" to "splicing unaffected".

**A subtle failure worth knowing.** Per-pair TVD is *unchanged* by swapping that pair, because
it is already an absolute value — so a naive paired permutation produces a null identical to
the observed value. The statistic has to be the TVD of the **averaged** difference vector, not
the average of per-pair TVDs.

### 3.4 Bootstrap confidence intervals

Used for quantities with no clean formula — the ratio of two means, for instance. Resample the
genes with replacement, recompute, repeat 2,000 times, and take the middle 95% of the results.

**Why here.** "What fraction of the wildtype response does the knockout retain" is a ratio of
means, whose sampling distribution has no simple closed form. The bootstrap sidesteps that.

**What it does not fix.** The bootstrap gives an honest interval *for the estimator you chose*.
If that estimator is biased, the bootstrap gives you a tight interval around the wrong number.
That is exactly what happened in §4.4.

---

## Part 4 — Traps this project has actually hit

Each of these produced a wrong result that was believed for a while.

### 4.1 Ratios estimated from small denominators

Isoform fraction is a ratio, so when a gene is barely expressed the fraction is estimated from
very few reads and swings between replicates. The symptom was visible without any modelling:
switching calls were **four times more common in the lowest expression quintile** than the
highest — backwards from what statistical power predicts, and the signature of an unstable
estimator rather than more biology.

**The fix.** Measure the noise directly. For each isoform, compute the standard deviation of
its fraction across replicates *within* a condition — condition is fixed, so that variation is
pure noise. Relate it to gene expression, and set a floor where the effect threshold sits at a
stated multiple of the noise.

At gene expression ≈ 12 the |dIF| threshold of 0.15 equals about 3 standard deviations of
replicate noise; below that, a "switch" is within about 2 SD of noise. The floor removes 33–52%
of raw calls, which is the intent.

**Generalisable lesson.** Whenever an effect size is a ratio, check whether the effect
correlates with the denominator. If small denominators show bigger effects, you are measuring
instability.

### 4.2 Mixing two multiple-testing universes

Significance was taken from one version of the data and presence from another. Because
q-values depend on the whole set of tests they were computed with, the two versions disagreed:
across 966 shared isoforms the effect sizes were **identical in all 966**, but the q-values
matched in only **132**, with 29 isoforms flipping across the significance threshold.

The visible symptom was a gene (CD86) that passed every criterion in the knockout yet was
classified as wildtype-only.

**Generalisable lesson.** Effect sizes are properties of the data; q-values are properties of
the *analysis*. Never take them from different sources.

### 4.3 Conditioning on the outcome

Several statistics in this project were computed on isoforms selected for being significant in
*both* datasets, then used to argue the datasets agree. Selecting on agreement and then
measuring agreement is circular — it will always look impressive.

These are retained but renamed with a `selconf_` prefix and labelled descriptive-only. The
defensible version conditions on significance in **one** dataset and measures the effect in the
other.

**The related trap: regression to the mean.** Select genes with a high value in both groups,
then compare a *second* measurement, and the group with generally lower values will look worse —
in either direction. The guard is a **reciprocal control**: run the matching both ways. If
matching on A shows a deficit in B *and* matching on B shows a deficit in A, most of what you
are seeing is regression to the mean. Only the *asymmetry* between the two directions is
evidence.

### 4.4 Comparing two estimators with different biases

The subtlest error, and the one that resolved a disagreement between two independent analyses.

The claim under test: the knockout retains 64% of its expression response but only 38% of its
splicing response, so the layers are decoupled.

Both numbers are ratios, both have tight bootstrap intervals, both are computed correctly. But
they use **different estimators with biases in opposite directions**:

- Expression (`|log2FC|` from 3-replicate means): averaging before the absolute value gives
  high signal-to-noise, but the absolute value still inflates small values, pulling the ratio
  **toward 1**. Deficits are understated.
- Splicing (TVD): a sum of absolute values, where noise accumulates rather than averaging out.
  For signal small relative to noise the response is roughly **quadratic**, compressing ratios
  **toward 0**. Deficits are overstated.

**How we established this rather than asserting it.** Simulate data with a *known* ratio and
see what each estimator reports (`10_layer_estimator_calibration.R`):

| true ratio | expression reports | splicing reports |
|---|---|---|
| 1.0 | 1.00 | 1.01 |
| 0.6 | 0.68 | 0.42 |
| 0.4 | 0.56 | 0.20 |

At a true ratio of 0.6 — both layers reduced *identically* — the two estimators report 0.68 and
0.42. **A gap of 0.26 appears from nothing.** That is almost exactly the gap in the original
claim. Inverting the calibration on the real data turns the observed +0.26 gap into −0.04.

**Generalisable lesson.** Before comparing two numbers, ask whether they are on the same scale.
Absolute values, ratios, bounded quantities and rank statistics all distort differently. If you
cannot answer analytically, simulate: generate data with a known answer and check that your
pipeline returns it. This is cheap and catches things nothing else will.

### 4.5 Counting nested terms as independent findings

GO and Reactome are hierarchies — one tight group of genes generates dozens of overlapping
terms. An analysis here returned **82 enriched terms from 15 genes**, of which **five genes
accounted for 95%** of all gene-term hits. The term count says nothing; the driver genes say
everything.

Every enrichment output therefore ships with a driver-gene summary alongside the term count.

### 4.6 Treating non-independent datasets as replicates

Two datasets in this project are the **same six sequencing libraries** quantified against two
different reference transcriptomes (gene-level correlation r ≈ 0.99 for matched samples). They
were being corrected for multiple testing as though independent, and their agreement described
as replication.

Agreement between them is **technical reproducibility on shared RNA**, not biological
replication. Corrections now report both a within-library-group and an across-all version, and
the library grouping is derived from the sample names rather than hardcoded.

**A related identifier trap.** Checking whether they were the same libraries by correlating
expression matrices on transcript ID gave r ≈ 0 — apparently conclusive evidence they were
different. That was wrong: `TCONS_*` identifiers are assigned *per reference*, so the same ID
denotes different transcripts in the two annotations. At gene-symbol level the correlation is
0.99. **Cross-reference comparisons must use stable identifiers**, and a near-zero correlation
where you expect a high one should prompt "is my key valid" before "is my hypothesis wrong".

---

## Part 5 — Reporting conventions adopted

1. **Effect size before p-value.** At genome scale p-values measure sample size. Medians,
   ratios and their intervals go first.
2. **Confidence intervals on every odds ratio.** A point estimate without one is unreadable.
3. **Denominators stated.** "23.4% reproduced" means nothing without knowing whether the
   denominator is all candidates or only testable ones. Both are given where they differ.
4. **Statistics conditioned on the outcome are labelled.** The `selconf_` prefix.
5. **Term counts travel with driver genes.**
6. **Nulls are exhaustive where the design allows**, so results do not depend on a random seed.
7. **Numbers in reports are computed from result tables at render time**, never typed in.
   Regeneration changes every number, and hardcoded prose goes stale silently.
8. **Specification sensitivity is reported when a result depends on an analysis choice.** Where
   a conclusion flips between two defensible normalisations, both are shown and the conclusion
   is that the data cannot resolve it — rather than picking one.

---

## Part 6 — The honest limits of this design

- **Three biological replicates per condition.** Enough to estimate an average, marginal for
  estimating a per-gene variance. Several analyses here are limited by that rather than by
  sample size in genes.
- **Two conditions, one time point.** No design-based way to establish direction of causation.
  "Splicing failure degrades the expression response" and "a blunted expression response
  produces less measured splicing change" fit the data identically.
- **Effect estimates inherited from upstream.** This pipeline consumes fitted results, so it
  inherits the upstream model's assumptions. Re-fitting from raw counts (`OPEN_QUESTIONS.md` #1)
  is the single change that would most improve what can be claimed.
- **Switching is tail behaviour.** The |dIF| ≥ 0.15 threshold sits near the 99th percentile of
  the whole distribution, so calls come from the extreme tail where estimator error dominates.
  The threshold is defensible; it should not be described as capturing typical behaviour.
