# Open questions and tests to revisit

Status as of 2026-08-07. Each entry states the question, what is blocking it, what would
answer it, and why it matters. Ordered by what would change a conclusion.

Companion documents: `STATISTICAL_METHODS.md` (why each method was chosen),
`REVIEW_CHANGES.md` (what was wrong before and what changed).

---

## Tier 1 — would change a stated conclusion

### 1. Re-fit both layers from raw counts

**Question.** Is the LPS response reduced more in one layer than the other in the UBL5
knockout?

**Blocked on.** Raw per-isoform counts for all genes (being supplied). Everything so far
uses IsoformSwitchAnalyzeR's condition means and replicate isoform fractions, not a fitted
model.

**Why the current answer is unsatisfactory.** The comparison depends on the estimator rather
than the biology. Expression response (`|log2FC|`) is biased *toward* no-difference; splicing
response (total variation distance) is biased *away* from it. Calibration against simulated
ground truth (`10_layer_estimator_calibration.R`) shows a true ratio of 0.6 reads as 0.68 for
expression and 0.42 for splicing — a fabricated gap of 0.26 from identical underlying change.

**What would answer it.** A beta-binomial or Dirichlet-multinomial model on isoform counts
gives a likelihood-based splicing effect with genuine standard errors, on the same footing as
a negative-binomial expression model. That removes the problem rather than calibrating around
it. `DRIMSeq`, `satuRn` or a custom `DirichletMultinomial` fit are all reasonable.

**Matters because.** This is the central Q2 claim. Two independent analyses currently
disagree on the raw numbers and agree only after an approximate correction.

---

### 2. Is the splicing deficit transcriptome-wide or targeted at the LPS programme?

**Question.** The external review reports that the deficit *spares* LPS-responsive genes
(62.8% retained for interferon-response genes vs 26.9% for non-responsive), which would argue
against "splicing gates the LPS response".

**Blocked on.** Nothing technically — but it uses the compressed splicing estimator, so the
magnitudes are unreliable. Needs redoing on a calibrated or model-based scale (see #1).

**What would answer it.** Recompute the per-gene-class retention with a linear estimator, then
test the genotype x responsiveness interaction.

**Matters because.** It is the single most direct test of the "splicing regulation is part of
how the LPS response is executed" model, and it currently points against it. If it survives
calibration, that is a real constraint on the hypothesis.

---

### 3. Mediation decomposition on a calibrated scale

**Question.** How much of the knockout's splicing deficit is explained by its reduced
expression response?

**Blocked on.** The published decomposition (36% indirect / 64% direct) uses the compressed
splicing estimator on both sides of the ratio, so both components inherit the bias.

**What would answer it.** Recompute once #1 lands. Note that even then, "direct" means *not
explained by the measured mediator* — reverse causation (splicing failure degrading the
expression response) fits the same decomposition. With two conditions and n=3 there is no
design-based way to break that symmetry.

---

## Tier 2 — would strengthen or qualify a result

### 4. Isoform-level functional consequence

**Question.** Do these switches change what protein is made, or only which transcript
carries it?

**Blocked on.** The pipeline carries no CDS, NMD or protein-domain annotation.

**What would answer it.** ISA's `analyzeORF` / `analyzeNMD` outputs, or re-deriving from the
GFF3. Then ask, for the compositional-rescue class specifically: are the lost minor isoforms
NMD targets or truncated, and the retained dominant one full-length coding? That converts
"purification" from a pattern into a mechanism.

**Matters because.** Gene-level enrichment is structurally blind here — both isoforms of a
switching gene sit in the same GO terms. This is the analysis that would say what switching
*does*.

---

### 5. Targeted validation of a cryptic switcher

**Question.** Do the switches invisible to differential expression hold up outside RNA-seq?

**Candidate.** **IRAK3** (IRAK-M, negative regulator of TLR signalling): gene-level log2FC
−0.077 — flat — while isoform fraction shifts 0.238. Also SOCS4, SP1, RAB7B, SORT1, LIG1.

**What would answer it.** Isoform-specific qPCR in WT vs KO, before and after LPS.

**Matters because.** It tests the hidden-layer claim directly, and for a gene where the
standard analysis would report nothing. It also tests the UBL5 splicing hypothesis without
depending on any of the normalisation arguments in #1.

---

### 6. `H_HT` annotation-bias check

**Question.** Is `H_HT` a genuine biological outlier or an artefact of its annotation?

**Evidence it is odd.** 2.07% switching genes vs 0.87% for `T_HT`; 3.11% of isoforms over
|dIF| 0.15 vs 0.91%; the *lowest* replicate noise (median IF SD 0.0241); loses only 14.2% of
calls to the expression floor against 35.7 / 48.4 / 58.7% elsewhere.

**The concern.** H is a different genotype quantified against a reference built partly from
its own PacBio run.

**What would answer it.** Compare the H-derived and T-derived transcript contributions to the
HT reference and ask whether H's expression concentrates on transcripts discovered in H.
Structurally the same test already run for the UT pair, which passed
(`02f_reference_concordance.R` machinery applies).

**Matters because.** `H_HT` is used as supplementary support for Q1, and it is also the one
dataset showing no ISG enrichment.

---

### 7. Curated ISG list instead of GO terms

**Question.** Does the interferon enrichment strengthen with a tighter gene-set definition?

**Current state.** Sets come from GO "response to type I/II interferon", which include
upstream signalling components (CDC37, FADD, HDAC4) alongside induced effectors. Enrichment
is OR ~6–8 with the loose definition.

**What would answer it.** A functional list — Schoggins et al. or Interferome. Drop a CSV at
`isg.custom_set_path` in config; nothing in `07_isg_analysis.R` changes.

**Matters because.** ISG enrichment is the project's strongest positive result, and a tighter
set would likely raise the odds ratios. Deliberately not reproduced from memory — that is a
place to introduce silent errors.

---

### 8. GO enrichment on the hidden-layer gene sets

**Question.** Are cryptic switchers or compositional-rescue genes functionally coherent?

**Blocked on.** Nothing — `docs/external_review/gene_sets/` has the symbol lists and the
tested universe; `05_pathway_enrichment.R` has the machinery.

**Caveat.** At 11–56 genes these sets are underpowered for GO. Expect a negative result, and
report driver genes alongside any term count per the project's standing rule.

---

## Tier 3 — housekeeping, but each is a real gap

### 9. Class-code composition against the expressed background

`03_novel_isoform_analysis.R` reports the class-code composition of switching isoforms, and
`04` compares composition *between* outcome classes. Neither compares against the **expressed
background**, which is what determines whether a class is enriched. The external review shows
class `j` (novel junction combination) is in fact slightly *depleted* among switchers
(OR 0.61, p = 0.055), not enriched — the impression of enrichment comes from `j` being 46% of
novel transcripts. Add the background comparison so the raw composition cannot be misread.

### 10. Re-examine the abundance floor under a paired noise model

The floor (gene expression ≥ 12) was derived from within-condition replicate IF noise pooled
without using the pairing. A paired noise estimate should be smaller, which would justify a
lower floor and recover some of the 33–52% of calls currently removed.

### 11. Confirm the `PENK` and Class B examples

`PENK` (gene −0.39, dominant isoform +0.21, minor −1.11, dominant IF 0.48 → 0.76) is the
clearest compositional-rescue example. The external review notes that `IRF3` and `BCL7B` —
named in an earlier draft — do **not** reproduce across annotations. Any example quoted in a
writeup needs its reproduction flag checked first.

---

## Resolved, recorded so they are not re-opened

| question | resolution | where |
|---|---|---|
| Are the ISA inputs a full tested background? | No — pre-reduced; context layer added | `REVIEW_CHANGES.md` §0, §0b |
| Should significance come from reduced or unfiltered objects? | Unfiltered, one FDR universe | §0c |
| Is there an abundance floor, and where? | Yes, gene expression ≥ 12 | §0c |
| Is the switching deficit explained by response magnitude? | Cannot be answered by stratification — mediator, not confounder | §0d, §0e |
| Are `T_HT` and `T_UT` the same libraries? | Yes, but different transcript spaces | §0d, §0f |
| Is the design paired? | Yes, confirmed by the experimentalist | §0e |
| Is the 64%-vs-38% layer asymmetry real? | Not supported — estimator artefact | §0f |
| Is baseline composition abnormal in the KO? | No, and the positive control works | §0e, §0f |
