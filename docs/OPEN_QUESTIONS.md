# Open questions and tests to revisit

Status as of 2026-08-10. Each entry states the question, what is blocking it, what would
answer it, and why it matters. Ordered by what would change a conclusion.

Companion documents: `STATISTICAL_METHODS.md` (why each method was chosen),
`REVIEW_CHANGES.md` (what was wrong before and what changed).

---

## Tier 1 — would change a stated conclusion

### 1. Re-fit both layers from raw counts

**Question.** Is the LPS response reduced more in one layer than the other in the UBL5
knockout?

**No longer blocked (2026-08-10).** The counts arrived at
`/Volumes/Expansion/IsoformSwitchAnalyzer/Counts/{H-T,U-T}/Isoforms/` — RSEM expected counts
for all twelve libraries of each reference, with `gene_id`/`GeneName`, isoform ids matching
the ISA objects 100%. Verified against the ISA objects in `REVIEW_CHANGES.md` §0g.

Because the U-T file carries **all twelve** libraries, the layer question becomes a
`genotype x treatment` interaction in **one** 2x2 model rather than a comparison of two
separately-fit ISA runs. That removes the cross-run normalisation problem as well as the
estimator problem. Counts are fractional (RSEM EM output) — fine for DRIMSeq/satuRn, round
for DEXSeq. `satuRn`, `DEXSeq`, `edgeR` and `limma` are installed; `DRIMSeq` is not.

**Choose the count scale deliberately.** ISA's own `isoformCountMatrix` is *not* RSEM
expected counts: its column totals match (40.5M vs 40.5M) but per-isoform it differs by a
factor of 0.18-2.3, because ISA derives counts from abundance (scaledTPM: TPM x libsize/1e6)
rather than carrying RSEM's effective-length-weighted expected counts. For DTU that is the
*recommended* scale — raw expected counts are biased by effective-length differences between
isoforms of the same gene. Both are derivable from the supplied files, so pick one on
purpose and record it.

Everything so far uses IsoformSwitchAnalyzeR's condition means and replicate isoform
fractions, not a fitted model.

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

### 6. `H_HT` — ANSWERED 2026-08-10, and it is neither option

**Question as posed.** Is `H_HT` a genuine biological outlier or an artefact of its
annotation?

**Answer: neither. `H_HT` has a different experimental design, and the pipeline adapted to
it.** `importRdata` detected two surrogate variables for `H_HT` and none for the other three,
so it applied `limma::removeBatchEffect` to H's expression and derived all its IF/dIF from
the corrected matrix. Its IFs are batch-corrected; the other three are raw. Full evidence and
reproduction from raw counts in `REVIEW_CHANGES.md` §0g.

**Why H and only H:** all four datasets have real replication — `T1`-`T3` were matured,
LPS-treated and RNA-extracted **separately**, not split from one flask. But `T`/`U` are one
clonal line processed on a single day, whereas `H` is human **primary** macrophages from
**three different donors**, each prepped on a **different day**. H therefore carries two extra
variance components (donor genotype, prep day) that T/U structurally lack, and sva detecting
them is correct behaviour. In H those two are **perfectly confounded** — one donor per day —
so no analysis of this data can separate genotype from batch.

Every listed oddity follows: the correction halves within-condition replicate IF noise
(0.0062 -> 0.0030) — uncorrected, H's donor variance would make it the *noisiest* dataset,
so the correction inverts the ranking — and the reduced noise plus `sv1`/`sv2` in the design
explains the elevated switching rate.

**What is now open instead.** Refit with donor as an **explicit** blocking factor
(`~ donor + condition`) rather than a latent one, and re-import all four under a single ISA
version. Do **not** simply set `detectUnwantedEffects = FALSE` — that would leave real donor
variation uncorrected. Details in §0g.

**Watch for.** `H`'s `sv1` is not orthogonal to condition (`H1_plus` +0.56, `H2_plus` +0.57,
`H3_plus` -0.46 against -0.21 to -0.23 for all minus samples), so the latent correction may
absorb real condition signal. An explicit donor term does not have this failure mode.

**Retire rather than recompute:** the H-vs-T switching-rate contrast in `02e` (OR 2.29/2.49).
A switching rate counts genes consistent across that dataset's replicates, so the two rates
answer different questions — "consistent across three donors on three days" versus
"consistent across three process replicates of one genotype on one day". H's is the harder
test. No modelling choice makes the ratio a genotype effect.

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

Now also entangled with #6: `H_HT`'s replicate IF noise was halved by a batch correction the
others did not get, so its recommended floor (6.02) is the lowest of the four — against
`T_HT` 10.07, `U_UT` 11.96, `T_UT` 20.16. `02d` takes the *median*, so the effect on the
applied value is modest (~11 vs ~12.5 if H were on the same footing), not a distortion. Still,
redo `02d` after the re-import rather than before.

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
| Do the ISA objects match the raw counts? | Yes for `T_HT`/`T_UT`/`U_UT` (dIF r = 1.0000); `H_HT` only after its batch correction is reapplied | §0g |
| Why is `H_HT` an outlier? | sva confounder correction applied to it alone | §0g, #6 |
