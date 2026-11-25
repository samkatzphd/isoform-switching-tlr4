# Q1 candidate and validation plan: a novel *ATF5* isoform

*Figure: `results/figures/ut_q1_atf5_candidate.png`. Oligos to order:
`results/tables/ut_q1_atf5_assay_order_sheet.csv`. Supporting tables, all under
`results/tables/`: `ut_q1_locus_fusion_screen.csv`, `ut_q1_atf5_isoform_fractions.csv`,
`ut_q1_atf5_effect_decomposition.csv`, `ut_q1_atf5_pair_discrimination.csv`,
`ut_q1_atf5_primer_specificity.csv`.*

> **Primer specificity of record is `results/tables/atf5_oligo_order_table.csv`**, not
> `ut_q1_atf5_primer_specificity.csv`. The latter is an earlier 8-oligo iteration covering only
> the N0/N3/T0/T3 qPCR pairs. The order table supersedes it: 16 oligos, each with a
> transcriptome-wide off-target count, including the D2 pair that gives the 514 bp endpoint
> RT-PCR product. Order from the order table.
>
> Provenance: `ut_q1_atf5_candidate.png`, `atf5_biology.png`,
> `ut_q1_atf5_pair_discrimination.csv` and `ut_q1_atf5_primer_specificity.csv` were built
> interactively on 2026-10-25 and committed on 2026-11-25. No script under `scripts/` produces
> the ATF5 tables, so rerunning the pipeline will not regenerate them.

## Summary

The Q1 question — does LPS drive a different isoform choice in wildtype THP-1 than in the
UBL5 knockdown — produces **one candidate**, not a list. Four genes reached significance for
the genotype × treatment interaction; three of them turned out to be two neighbouring genes
merged into one locus by the assembler and were disqualified. The survivor is *ATF5*, and the
switching species is a **novel isoform that only long-read sequencing could have found**.

I want to be blunt about the size of the list: this is one gene out of 13,086 tested, and the
honest framing is a single well-characterised candidate rather than a screen. The reason to
take it seriously is not the statistics alone — it is what the gene turned out to be.

## Why this candidate is worth bench time

ATF5 is the mammalian functional ortholog of *C. elegans* ATFS-1, the master transcription
factor of the mitochondrial unfolded protein response (UPR^mt). UBL5 is the human homolog of
*C. elegans* UBL-5, which acts with DVE-1 as the chromatin-level coregulator that enables
ATFS-1-dependent transcription of UPR^mt genes. ATFS-1, DVE-1 and UBL-5 are the three nuclear
components of the same pathway.

So an unbiased transcriptome-wide interaction scan in a UBL5 knockdown returned, as its only
surviving hit, UBL5's own direct pathway partner. That is not what noise looks like. It also
points at the UPR^mt arm of UBL5 biology rather than the spliceosome arm — though since the
event *is* an isoform switch, the interesting possibility is both at once: UBL5's splicing
function acting on the transcript of its own pathway partner.

## What was measured

Novel isoform **TCONS_00159376**, isoform fraction of the *ATF5* locus (n = 3 libraries per arm):

| arm | isoform fraction | replicates |
|---|---|---|
| WT unstimulated | 0.112 | 0.114, 0.114, 0.107 |
| WT + LPS | 0.135 | 0.143, 0.132, 0.129 |
| KD unstimulated | 0.185 | 0.187, 0.181, 0.188 |
| KD + LPS | 0.107 | 0.109, 0.095, 0.117 |

Decomposed:

| quantity | value |
|---|---|
| baseline difference, KD − WT | **+0.074** |
| LPS effect in WT | +0.023 |
| LPS effect in KD | −0.078 |
| interaction | **−0.101** (q = 0.020) |
| difference after LPS, KD − WT | −0.028 |

Read that as a pattern rather than four numbers: **the knockdown starts with more of the novel
isoform, and LPS erases the difference.** UBL5 appears to restrain this isoform under resting
conditions, and LPS overrides the restraint. Note that the single largest effect is the
*baseline* difference, not the LPS response.

## Structure and predicted consequence

The novel isoform differs from the dominant reference isoform (TCONS_00149857) in two ways:
it starts from an alternative first exon 378 bp upstream, and it terminates early, reading 93
bp into what is an intron in the reference and stopping there instead of splicing to the large
terminal exon.

Translating it from the canonical ATG gives **88 amino acids terminating in a genuine stop
codon**, against 282 aa for full-length ATF5 (which my translation reproduced exactly against
Ensembl, so the frame is right). The first 62 residues are shared with full-length ATF5.

That architecture matters: ATF5's mitochondrial targeting sequence is N-terminal and its
bZIP DNA-binding/dimerisation domain is C-terminal. An 88-aa product retains the targeting
region and loses the bZIP entirely — the classic shape of a dominant-negative or a decoy.
**And it should not be degraded**: the stop codon sits in the final exon, 383 nt downstream of
the last exon–exon junction, so it is not predicted to be an NMD substrate. If it is made, it
should persist.

I should flag that the domain assignment is from the literature, not measured here; the
88-aa length and the absence of the terminal exon are measured.

## The assays

### One structural caution that shaped the design

The novel isoform contains **no unique sequence at all** — every 20-mer in it also appears in
some other *ATF5* isoform. Its novelty is combinatorial: it is the only isoform that has both
the upstream first exon *and* the 93-bp read-through. So no single primer can be
isoform-specific, and the discrimination has to come from a pair that spans both features.
Two reference isoforms each carry one of the features: TCONS_00149856 has the upstream exon,
TCONS_00149858 has the read-through. Neither has both.

All oligos below were screened against all 315,349 transcripts in the long-read annotation:
**zero off-target loci**, and every primer position was verified back against the spliced cDNA.

### Assay 1 — does the isoform exist? (do this first)

Endpoint RT-PCR with **pair D2**, then Sanger the product.

| | sequence | position |
|---|---|---|
| forward | `CACAAGGAGAGGGCTGGATC` | chr19:49,928,968–49,928,987 (upstream exon) |
| reverse | `CACCTTCCACCTGCCCTTAC` | chr19:49,931,029–49,931,048 (read-through) |

Product **514 bp**, and it amplifies this isoform and no other — verified against all seven
*ATF5* models. A single band of the right size, with Sanger confirming the exon-1/exon-2
junction and the read-through terminus, establishes the structure. This is the assay that
decides whether anything else is worth doing.

### Assay 2 — quantify the switch

RT-qPCR, ratio of two pairs. Run all four arms (genotype × LPS).

- **Numerator, pair N0** — `AATCACAGGAGGAAGCAGGC` / `CAAAGGGGTCGAGGCTACAG`, 187 bp.
- **Denominator, pair T0** — `TCTCCTGTTGCCATCAGTGC` / `AGTCTACGAGCCATCCCAGC`, 119 bp, both
  primers inside one contiguous block shared by all six expressed isoforms, so it measures
  total *ATF5* output.

**Be aware of what N0 actually measures.** It is not novel-isoform-specific: it also picks up
TCONS_00149856, at the same product size, because that isoform shares the upstream exon. At
baseline the numerator is roughly 91% novel isoform and 9% TCONS_00149856, and TCONS_00149856
barely moves between arms (−0.001 in WT, −0.005 in KD), so it dilutes the signal slightly
without reversing it. Call this pair "upstream-TSS isoforms", not "the novel isoform". The
514-bp D2 assay is what attributes the signal to the novel species.

A fraction is the readout, so the denominator is not optional — LPS moves total *ATF5* output
substantially (WT 240 → 330 TPM; KD 128 → 223 TPM), and an un-normalised isoform Ct would be
uninterpretable.

### Assay 3 — confirm the 5′ end

5′ RACE with a gene-specific primer in the shared exon 2, reading back into the upstream exon.
The alternative start is called from long-read models quantified with short reads, and a TSS
deserves independent confirmation before it is described as an alternative promoter.

### Assay 4 — is the 88-aa protein made?

Western blot with an **N-terminal** ATF5 antibody (a C-terminal antibody cannot see this
product by construction). Expect full-length ATF5 near 32 kDa and, if the isoform is
translated, an additional band near 10 kDa. This is the cheapest test of the whole hypothesis
and the one that converts a transcript observation into a protein claim.

### Assay 5 — is it a dominant negative? (only if 1 and 4 succeed)

Express the 88-aa form in wildtype THP-1 and measure induction of UPR^mt targets — HSPD1,
HSPA9, LONP1, CLPP — under mitochondrial stress. If the truncated protein competes with
full-length ATF5, induction should be blunted.

## Controls

**The genotype control must be the sh-Luciferase line, not parental THP-1.** The sequenced
comparison used untransduced wildtype, which means every WT-versus-KD difference in the
sequencing carries retroviral transduction and G418 selection alongside UBL5 loss. The
baseline difference (+0.074) is the part of this result most exposed to that confound, and it
is the largest single effect — so the bench assay is where it gets resolved. Running sh-Luc
fixes at the bench what the sequencing cannot.

Also required: no-RT control on every RNA (pair T0 has both primers in one exon and will
amplify genomic DNA); no-template controls; standard curves with 90–110% efficiency for both
qPCR pairs, since a ratio between pairs of differing efficiency is not valid; and a melt curve
showing a single peak per pair.

## Falsification criteria

`docs/ASSAY_DESIGNS.md` uses a project-wide floor of |dIF| ≥ 0.16 at n = 3. **That floor is
the wrong yardstick here** — it is a genome-wide figure, and this locus is far less variable
than average: the pooled within-arm SD of the isoform fraction is 0.0073, giving a
locus-specific detection floor of ±0.027 for the interaction at n = 3. The observed
interaction of 0.101 clears that by about 3.8×.

Two honest cautions on that. First, the locus was *selected* because it was significant, which
selects for small standard error — so 0.0073 is biased downward and the true variance is
probably larger. Second, qPCR adds technical variance that sequencing does not have. The truth
sits between 0.027 and the project-wide 0.16, and the only way to settle it is to measure
assay variance directly.

So: **run a pilot of n = 3 per arm to estimate the qPCR variance, then power the real
experiment on that number.** At the sequencing-derived variance, n = 3 would suffice; at the
project-wide floor it would need roughly n = 16. Plan for n = 6 as a working compromise and
revise once the pilot variance is in hand — THP-1 flasks are cheap enough that this is not a
painful contingency.

**Supports the result if:** the novel isoform's share of the locus is higher in unstimulated
knockdown than in unstimulated sh-Luc, and falls with LPS in the knockdown while holding or
rising in sh-Luc, with the interaction exceeding the pilot-derived detection floor.

**Falsified if:** the direction reverses, or the 514-bp product is absent, or Sanger shows a
different structure.

**Uninformative rather than negative if:** the shift is below the pilot-derived floor. That
means the experiment was underpowered, and the conclusion is to add replicates.

## What this does not establish

The interaction is significant under Benjamini-Hochberg but the empirical-null q-value for
this contrast is 0.20, and the script's own diagnostic marks the empirical null unusable on
this dataset. The 2×2 sensitivity model with a shared replicate factor returns zero
interaction genes. Both are reasons to treat *ATF5* as a hypothesis to test at the bench
rather than a finding — which is exactly what the assays above are for.
