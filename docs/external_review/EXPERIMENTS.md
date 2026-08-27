# Proposed wet-lab follow-up

Organised around the two questions as you framed them, which are separable:

- **Q1. What is the role of alternative splicing in the LPS/TLR4 response?**
- **Q2. Is UBL5 an interesting splicing-regulation candidate, and what is it doing?**

Experiments are ordered so that the cheapest result that could kill a hypothesis comes
first. Each states the design, the readout, and what outcome would falsify the hypothesis.

---

## BLOCKER — information I do not have

Three things are needed before primers can be designed. Everything below marked
**[NEEDS COORDINATES]** is specified down to the level of "which junction to span", but the
actual oligo sequences require:

1. **The merged GFF3 annotations** — `/Volumes/Expansion/gtf_files/GTF_to_GFF311_3_2022/U_T/U_T.gff3`
   and the `H_T` equivalent. The drive was not mounted during this analysis. Without exon
   coordinates I cannot give you primer sequences, and — more importantly — I cannot
   distinguish alternative-first-exon events from internal splicing changes (see the AFE
   problem below, which is currently the single biggest threat to the UBL5 interpretation).
2. **The UBL5 knockout allele.** Is this a full knockout, a frameshift, a hypomorph, a
   degron, or a knockdown? Published work reports human UBL5 as essential for viability
   (Ammon et al. 2014), so a viable complete null needs explaining. If it is a hypomorph or
   partial knockdown, that changes the interpretation substantially and reconciles a
   literature contradiction (see `literature_review.md`, Q2).
3. **Whether the RNA-seq libraries are polyA-selected or ribo-depleted, and whether the
   three replicates are paired batches** (the `T1/T2/T3` naming suggests paired). Intron
   retention — the readout most diagnostic for UBL5 — is poorly captured by polyA selection.

---

## The AFE problem — read this before designing anything

The literature review surfaced a finding that, if it applies here, undercuts the central
UBL5 inference: **roughly half of splicing changes after LPS in human and murine macrophages
are alternative-first-exon (AFE) events** (Robinson et al. 2021). AFE is promoter/TSS choice,
not spliceosome catalysis. IsoformSwitchAnalyzeR scores AFE as isoform switching.

If a large share of the "splicing response" in this dataset is AFE, then attributing its loss
in the UBL5 knockout to a spliceosome modifier is mechanistically incoherent — a spliceosome
protein should not control promoter selection.

I could not resolve this from the committed data: `class_code` gives only a weak proxy
(17.5% of switching isoforms are `c`/truncation-like, 12.6% are `j`/internal-junction, and
69.2% are `=` reference transcripts whose AFE status is unresolvable without exon
coordinates).

**Recommended first analysis, before any bench work:** with the GFF3 in hand, classify every
switching event as AFE / alternative-last-exon / exon skipping / intron retention /
alternative 5' or 3' site (SUPPA2 will do this from the GTF), then recompute the 38%
retention figure within the catalytic classes only. If the UBL5 deficit lives in intron
retention and exon skipping, the hypothesis strengthens sharply. If it lives in AFE, the
mechanism must be rethought. This is a day of compute that determines whether the wet-lab
program below is worth running.

---

## Q1 — Splicing in the LPS response

### Q1.1 Validate the "cryptic switcher" class by isoform-specific RT-qPCR

**Rationale.** 19.5% of switching genes show almost no gene-level expression change. If real,
this is the project's most transferable claim: standard DE misses a fifth of the regulated
genes. It needs orthogonal validation on independent RNA.

**Design.** THP-1 (PMA-differentiated, matching the original protocol), ± LPS at the same
dose/timepoint as the RNA-seq, n = 4 independent biological replicates, plus a 0/2/4/8/24 h
time course on the top 4 targets. New RNA, not the sequenced material.

**Targets and the specific discrimination required** — all **[NEEDS COORDINATES]** for final
oligo design, but the junction to target is determined:

| Gene | Isoform that FALLS | Isoform that RISES | Assay strategy |
|---|---|---|---|
| IRAK3 | `ENST00000261233.9` (IF 0.85→0.61) | `PB.4682.2`, a `c`-class transcript contained within ENST00000261233 (IF 0.03→0.21) | The rising isoform is a truncation of the falling one. A junction-spanning primer pair cannot discriminate a contained isoform — you need a primer in the region **unique to the long form** (present in ENST00000261233, absent from PB.4682.2) for the falling species, and a boundary-spanning pair at the PB.4682.2-specific 5' or 3' terminus. If PB.4682.2 differs only by truncation with no novel junction, **use 3' RACE or long-read targeted sequencing instead of qPCR** — qPCR cannot resolve nested isoforms. |
| SPRING1 (C12orf49) | `PB.4684.4` (IF 0.245→0.057) | `PB.5028.2` (IF 0.106→0.391) | Both are `c`-class relative to ENST00000261318. Same nesting problem; check junction structure first. |
| RAB7B | `PB.1663.8` (IF 0.32→0.16) | `PB.1802.12` (IF 0.08→0.26) | Two distinct PB transcripts, most likely resolvable by junction-spanning primers. Best qPCR candidate of the four. |
| SOCS4 | `ENST00000339298.2` (IF 0.27→0.13) | `ENST00000395472.2` (IF 0.14→0.32) | Both are annotated Ensembl transcripts — design from Ensembl directly, no GFF3 needed. **Start here**; it is the only one designable today. |

**Controls.** (i) A gene-level primer pair per gene in a constitutive region, to confirm total
expression does not change — that is the whole point of "cryptic". (ii) Two housekeepers not
in any switching set. (iii) A positive-control switcher with a large gene-level change
(from the >1 log2FC bin) to show the assay detects switching when DE also does.

**Falsified if** isoform ratios do not shift while gene-level expression stays flat.

### Q1.2 Is the switching functional, or is it transcriptional noise?

**Rationale.** The strongest version of the Q1 claim is that splicing is a *regulatory* layer,
not a byproduct. The cleanest published architecture is dominant-negative isoform feedback
(Lee & Alper 2022), and **IRAK3/IRAK-M is already a negative regulator of TLR signalling** —
making its isoform switch the highest-value functional target in the dataset.

**Design.** Isoform-specific overexpression in THP-1: lentiviral constructs for the long
(ENST00000261233) and short (PB.4682.2) IRAK3 forms, plus empty vector, in both wildtype and
IRAK3-knockout backgrounds. Stimulate with LPS; read out TNF, IL6, IL1B, CXCL8 by qPCR and
ELISA at 2, 6, 24 h, plus NF-κB reporter and IκBα degradation by western.

**Prediction.** If the short form is dominant-negative (it lacks part of the region present in
the long form), overexpression should *increase* inflammatory output — the opposite of the
long form.

**Falsified if** the two isoforms behave identically.

### Q1.3 Splicing inhibition — does blocking splicing blunt the LPS response?

**Design.** Pladienolide B or herboxidiene (SF3B1 inhibitors), titrated to sub-lethal doses
that perturb splicing without killing cells; pre-treat 2 h, then LPS. Read out the same
cytokine panel plus RNA-seq at one dose.

**Interpretation.** This tests the general claim (splicing required for a full LPS response)
without touching UBL5. **Important caveat:** SF3B1 inhibition is catastrophic and pleiotropic,
so a blunted response is weak evidence. Include a dose-response and a viability/translation
control; treat a positive result as consistent-with rather than evidence-for.

---

## Q2 — UBL5

### Q2.1 The essential control: is the KO a hypomorph? (do this first, it is cheap)

**Rationale.** Two literature findings must be reconciled before anything else: human UBL5 is
reported essential (Ammon et al. 2014), and UBL5 depletion is reported to increase intron
retention constitutively (Oka et al. 2014) — whereas this dataset shows **normal baseline
isoform composition**. A hypomorphic allele explains all three observations at once.

**Design.** (i) Sequence the edited locus (amplicon sequencing across the guide site) and
report the exact allele. (ii) Western blot for UBL5 protein in WT and KO, quantified against
a dilution series of WT lysate so residual protein can be estimated as a percentage, not just
"absent". (iii) qPCR for UBL5 mRNA including a primer pair spanning the edit site.

**Why it matters.** If residual UBL5 protein is, say, 20% of WT, then "no baseline defect but
impaired stimulus response" becomes the expected phenotype of a partial loss-of-function under
increased splicing load — a coherent and publishable model. If protein is truly absent, the
contradiction with the essentiality literature must be addressed head-on.

### Q2.2 Direct test of the central claim: intron retention under stimulation

**Rationale.** This is the single most diagnostic experiment for the whole UBL5 hypothesis.
UBL5/Hub1 acts on the spliceosome via SART1/Snu66 — the STRING interactome confirms this
independently (51 of 82 UBL5 physical partners are core spliceosome/RBP, including a direct
SART1 edge). If UBL5 supports splicing under load, **intron retention should rise in the KO
specifically after LPS, not at baseline.**

**Design.** WT and UBL5-KO THP-1, ± LPS, n = 4, **ribo-depleted total RNA-seq** (not polyA —
retained-intron transcripts are frequently non-polyadenylated and would be lost). If budget
allows, add a 2 h and 6 h timepoint, since the model predicts the defect emerges with load.

**Readout.** Intron retention by IRFinder or SUPPA2; the key statistic is the
**genotype × stimulation interaction**, not the main effect of genotype.

**Prediction (and the falsification condition).** IR increases in KO+LPS above the additive
expectation. **If IR is elevated in the KO at baseline too, that contradicts this dataset's
baseline result and supports the published constitutive-defect model instead** — an
informative outcome either way, and the reason to run this before the more expensive work.

**Targeted confirmation:** RT-PCR across specific retained introns with primers in flanking
exons, resolved on a gel or Bioanalyzer so the retained and spliced species are visible as
separate products. **[NEEDS COORDINATES]** — intron choice depends on the RNA-seq above.

### Q2.3 Does UBL5 physically engage the spliceosome more under LPS?

**Design.** Co-IP of endogenous UBL5 (or a knock-in tagged allele, to avoid overexpression
artefacts) from WT THP-1 ± LPS; blot for SART1, SF3B1, PRPF8, SNRNP200. Quantify the ratio of
co-precipitated partner to precipitated UBL5.

**Prediction.** If UBL5's role is load-dependent, association should increase after
stimulation. A flat result does not kill the hypothesis (regulation could be by activity
rather than binding) but a positive is strong support.

**Extension.** UBL5 conjugation is non-canonical; if you have antibodies or an epitope-tagged
line, test whether SART1 modification changes with LPS.

### Q2.4 Rescue — the specificity control

**Design.** Re-express wildtype UBL5 in the KO, alongside (i) empty vector and (ii) a
conjugation-deficient C-terminal mutant (the Hub1 C-terminal double-tyrosine region;
**[NEEDS the exact residue numbering for human UBL5]** — I can specify this from UniProt
P0CG47 if you want it, but the construct should be designed against the published Hub1
structure-function work rather than my inference).

**Readout.** Does wildtype UBL5 restore the LPS splicing response (by the same TVD metric
computed on new RNA-seq), and does the conjugation-deficient mutant fail to?

**Why it matters.** This is what converts a correlation between genotype and splicing
phenotype into a statement about UBL5's molecular activity. Without a rescue, a reviewer will
reasonably attribute the phenotype to clonal drift in the KO line.

### Q2.5 Independent perturbation — clonal artefact control

**Design.** siRNA or shRNA knockdown of UBL5 in wildtype THP-1 (two independent
sequences), plus a second independently derived KO clone. Repeat the Q2.2 readout.

**Rationale.** A single CRISPR clone carries clonal history. The observed phenotype — reduced
LPS response amplitude — is exactly what a slow-growing or partially differentiated clone
would show. **This control is not optional for publication.**

---

## Cross-cutting: what the sequencing-side analysis should be, given n=3

The current analysis uses ISA's condition means and a composition distance I computed. Neither
is a formal differential-transcript-usage test. Before any of the above anchors a manuscript:

- Re-test differential transcript usage with **DRIMSeq or satuRn**, including the replicate
  pairing if the design is paired (see BLOCKER 3).
- Classify events by type (SUPPA2) and recompute the headline retention figure within
  catalytic classes (the AFE problem).
- The abundance floor (gene expression ≥ 12) removes 33–59% of calls; report results with and
  without it.

---

## Suggested sequence

1. **Q2.1** (allele characterisation) — cheap, and its outcome reframes everything else.
2. **AFE reclassification** (compute, needs GFF3) — determines whether Q2 is coherent.
3. **Q2.2** (IR under stimulation, ribo-depleted) — the central falsification test.
4. **Q1.1** starting with SOCS4 (designable today) — validates the cryptic-switcher class.
5. **Q2.5** (independent perturbation) in parallel with Q2.2.
6. **Q2.4** (rescue) and **Q1.2** (IRAK3 function) — only once 1–3 support the model.

## What I need from you to go further

- The two GFF3 files (or just the transcript/exon records for the ~15 candidate genes) —
  then I can produce actual primer sequences, amplicon sizes, and specificity checks.
- The KO allele description.
- Library prep chemistry and whether replicates are paired batches.
- Confirmation of the LPS dose, timepoint, and PMA differentiation protocol used, so the
  validation experiments match the original conditions.
