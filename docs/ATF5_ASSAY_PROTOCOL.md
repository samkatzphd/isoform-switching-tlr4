# Bench protocol: validating the novel *ATF5* isoform

*Companion to `docs/ATF5_DEEP_DIVE.md` (what the transcript is and why it matters) and
`docs/ATF5_Q1_EXPERIMENTAL_PLAN.md` (how it was found). Oligos: `atf5_oligo_order_table.csv`.*

## 0. Where each parameter comes from

The upstream half is **not new** — cells, differentiation, LPS and RNA are reproduced from the
original experiment so that a validation run cannot confuse a wrong candidate with a wrong
experiment.

| step | source | status |
|---|---|---|
| THP-1 culture and PMA differentiation | Thesis §2.5.1 (see `PROTOCOLS.md` §1.1) | reproduced verbatim |
| Primary human macrophage isolation | Thesis §2.5.6 (`PROTOCOLS.md` §1.2) | reproduced verbatim |
| LPS source, dose, duration | Thesis §2.6.1, §7.2, Fig. 7.1 (`PROTOCOLS.md` §1.3) | reproduced verbatim |
| RNA extraction | Thesis §2.6.8 (`PROTOCOLS.md` §1.4) | reproduced verbatim |
| Reverse transcription, qPCR cycling | Thesis §2.6.7 | reproduced, with the priming change in §2.1 |
| Endpoint RT-PCR on gel | Thesis §2.6.6 | reproduced verbatim |
| **Assays 1–6 below, all oligos, reference-gene panel** | — | **NEW** |

Key upstream parameters, repeated here so nobody has to cross-reference mid-bench: THP-1
differentiated with **5 ng/ml PMA for 72 h** (do not substitute a higher dose — PMA dose
changes the LPS response); primary MDM with **10 ng/ml GM-CSF for 7 days**, one donor = one
biological replicate; LPS **R595 rough chemotype at 100 ng/ml for 4 h**, paired untreated
subset, both arms of a replicate handled together.

## Genotype arms — read before planning

Use **sh-Luciferase as the knockdown's comparator, not parental THP-1.** The sequencing
compared the knockdown against untransduced wildtype, so every genotype difference in the data
carries retroviral transduction and G418 selection alongside UBL5 loss. The bench can fix this;
the sequencing cannot. Four arms per experiment: sh-Luc ± LPS, UBL5-KD ± LPS.

The LPS-induction readout *within* each line is internally controlled and is the primary
result (§2.6), so it survives even if the genotype comparison is contested.

---

## Assay 1 — Does the transcript exist? (do this first)

**Endpoint RT-PCR with pair D2, then Sanger the product.** Decides whether anything else is
worth doing. One afternoon.

1. RT 500 ng total RNA with **random hexamers** (see §2.1 for why not oligo-dT).
2. PCR: ATF5_D2_F + ATF5_D2_R, a proofreading polymerase, 35 cycles, Ta 60 °C, extension 35 s.
3. Run 15 µl on a 1.5% agarose gel alongside a 100 bp ladder.
4. Gel-extract the band and Sanger from both primers.

**Expected: a single 514 bp product.** It amplifies TCONS_00159376 and no other *ATF5* model —
verified against all seven models and all 315,349 transcripts in the long-read annotation (0
off-target loci).

**Sanger must confirm two things:** the exon 1 → exon 2 junction at chr19:49,929,006 |
49,929,243, and continuous read-through past chr19:49,931,028 into the 93-bp extension. Both
are required; the first alone is also present in TCONS_00149856.

**Genomic DNA gives 750 bp**, not 514, because the amplicon spans a 236-bp intron. A 750-bp
band means DNA contamination, not a new isoform.

**Decision point.** No 514 bp band in wildtype THP-1 ± LPS, with positive controls working →
stop. The transcript is an assembly artefact and the arm closes here, cheaply.

---

## Assay 2 — Quantify the induction (the primary quantitative readout)

The headline result to reproduce is **absolute LPS induction of the isoform in wildtype
(1.66-fold, p = 0.004) and its complete absence in the knockdown (1.00-fold, p = 0.96)** — not
a change in isoform fraction. Report both, but power the experiment on the induction.

### 2.1 Reverse transcription — one deliberate change from Thesis §2.6.7

**Use random hexamers, not oligo-dT.** This is the one real deviation and it matters: the
numerator amplicon sits ~190 nt from the 5′ end of a 728-nt transcript while the denominator
amplicon sits mid-transcript in isoforms up to 2,007 nt. Oligo-dT priming biases coverage
toward 3′ ends and would distort a ratio between amplicons at different distances from the
poly(A) tail. Random priming is approximately even.

500 ng total RNA per 20 µl reaction; dilute cDNA 1:5 before use; run every arm of a replicate
in the same RT batch.

### 2.2 Reactions per sample

Five pairs, triplicate technical replicates:

| pair | measures | product |
|---|---|---|
| ATF5_N0_F / ATF5_N0_R | upstream-TSS isoforms (numerator) | 187 bp |
| ATF5_T0_F / ATF5_T0_R | total *ATF5*, all 6 expressed isoforms (denominator) | 119 bp |
| HPRT1_F / HPRT1_R | reference | 145 bp |
| TBP_F / TBP_R | reference | 80 bp |
| YWHAZ_F / YWHAZ_R | reference | 145 bp |

**What N0 actually measures.** It is *not* novel-isoform-specific: it also amplifies
TCONS_00149856, at the identical product size, because that isoform shares the upstream first
exon. At baseline the numerator is ~91% novel isoform and ~9% TCONS_00149856, and
TCONS_00149856 barely moves between arms (26.8→44.5 TPM for the novel isoform versus 2.6→3.1
for TCONS_00149856), so it dilutes slightly without reversing anything. **Call this pair
"upstream-TSS isoforms" in any figure or text.** Assay 1 is what attributes the signal to the
novel species. There is no way around this: the novel isoform contains no unique sequence at
all — every 20-mer in it occurs in some other *ATF5* isoform — so the only discriminating
design is the 514-bp D2 pair spanning both distinguishing features.

### 2.3 Reference genes — chosen from these data, not from convention

HPRT1, TBP and YWHAZ, geometric mean of the three. Selected by stability across the twelve
THP-1 libraries and checked for LPS responsiveness:

| gene | CV across libraries | LPS fold (WT) | p | off-target loci |
|---|---|---|---|---|
| HPRT1 | 0.037 | 1.00 | 0.97 | 0 |
| TBP | 0.069 | 0.99 | 0.85 | 0 |
| YWHAZ | 0.071 | 1.08 | 0.23 | 2 (0.43% of locus output) |

**Do not use GAPDH.** It falls to 0.80× with LPS in wildtype (p = 0.003) and 0.86× in the
knockdown, which would inflate apparent induction by about 25%. ACTB is borderline (0.90,
p = 0.14). **PPIA was designed and rejected** — its pair amplifies four pseudogene loci and
PPIAP22 alone sits at 145 TPM, so it cannot be measured cleanly by PCR despite looking stable
in RNA-seq.

### 2.4 Validation before any biological comparison

- Standard curve over ≥5 logs for **all five pairs**, efficiency 90–110%. A ratio between two
  pairs of differing efficiency is not valid, and this assay is a ratio.
- Melt curve on every reaction: a single peak. These primers were screened against the
  transcriptome, but a second peak is a second product regardless.
- **No-RT control is mandatory, not optional.** ATF5_T0 has both primers inside one exon
  (chr19:49,930,811–49,930,929) and amplifies genomic DNA at the identical 119 bp. D2 and N0
  span an intron and are protected; T0 is not.
- No-template control on every plate.

### 2.5 Cycling — Thesis §2.6.7

As published. SYBR chemistry; 60 °C anneal suits all five pairs (Tm 58.9–61.4 °C).

### 2.6 Calculation and call

1. **Primary:** normalise N0 to the geometric mean of the three reference genes → relative
   abundance of the upstream-TSS class. Compare +LPS vs −LPS **within each line**. Expect
   ≈1.7-fold in sh-Luc and ≈1.0-fold in the knockdown.
2. **Secondary:** N0/T0 as an isoform fraction, to connect back to the sequencing dIF.
3. **Context:** T0 normalised to reference genes, to show the gene itself is induced in both
   lines (and harder in the knockdown) while the isoform is not.

**Replicates.** Run **n = 3 as a pilot to estimate assay variance, then power the real
experiment on that number.** The project-wide |dIF| ≥ 0.16 floor is the wrong yardstick for
this locus — its pooled within-arm SD of isoform fraction is 0.0073, giving a locus-specific
n = 3 interaction floor of ±0.027 — but that SD is biased low because the locus was selected
for significance, and qPCR adds variance that sequencing does not. Plan n = 6 as a working
compromise; the induction readout is far better powered than the fraction readout.

---

## Assay 3 — 3′ RACE: is the intronic poly(A) site real?

This is the assay for the most novel structural claim. The transcript carries a canonical
**AATAAA at chr19:49,931,106**, 15 nt from its 3′ end — a textbook poly(A) signal — so the
prediction is that transcription genuinely terminates inside intron 2 rather than splicing on.

1. RT with an **adapter-anchored oligo-dT primer** (kit-supplied; e.g. SMARTer/GeneRacer
   3′ primer).
2. First PCR: **ATF5_3RACE_GSP1** + adapter primer. Touchdown 70 → 64 °C; GSP Tm is 67.7 °C.
3. Nested PCR: 1:50 dilution of round 1 as template, **ATF5_3RACE_GSP2** + nested adapter.
4. Resolve on a 1.5% gel; gel-extract every band and Sanger.

**Expected ladder from GSP2** (all three products are real transcripts, and their sizes
distinguish them):

| 3′ end used | isoform | product from GSP2 |
|---|---|---|
| chr19:49,931,121 — intronic poly(A) | **novel TCONS_00159376** | **≈214 bp** + adapter |
| chr19:49,931,229 | TCONS_00149858 | ≈322 bp + adapter |
| chr19:49,933,935 — canonical, spliced | dominant TCONS_00149857 | ≈1,635 bp + adapter |

A 214-bp product whose Sanger read ends at 49,931,121 followed by poly(A) **establishes the
intronic polyadenylation event.** Its absence, with the 1,635-bp canonical band present, says
the 3′ end is an assembly artefact.

**Control:** a GSP against a gene with a single annotated poly(A) site, same reaction batch, to
show the adapter ligation and dT priming worked.

---

## Assay 4 — 5′ RACE: is the alternative TSS real, or an assembly truncation?

Worth being precise about what this tests. The novel isoform is annotated as starting at
chr19:49,928,827, while TCONS_00149856 — which shares the same first-exon 3′ end — starts
**125 nt further upstream at 49,928,702.** So there are two possibilities, and they matter:

- The two isoforms genuinely use different start sites, and the novel isoform has its own TSS.
- They share one TSS and the novel model is simply truncated at the 5′ end by the assembler —
  in which case **the only real difference between them is the 3′ end**, and the story
  simplifies to a pure intronic-polyadenylation event.

1. RT with **ATF5_5RACE_GSP1** (gene-specific, Tm 66.5 °C), then template-switch or
   terminal-transferase tailing per kit.
2. Nested PCR with **ATF5_5RACE_GSP2** (in exon 1, Tm 67.7 °C) + adapter primer.
3. Clone or Sanger directly; **read ≥8 clones** — TSS usage is a distribution, not a point.

**Expected ≈143 bp + adapter** if the annotated start is correct; ≈268 bp if the TSS is shared
with TCONS_00149856 at 49,928,702.

**Controls:** no-RT on the same RNA (a band means genomic DNA and every TSS call from that
sample is void); no-template in both rounds; a positive-control GSP against a single-TSS gene
in the same batch.

> Note for prioritisation: if several alternative-TSS candidates accumulate across the project,
> CAGE/nanoCAGE resolves them all in one experiment per condition and also settles whether
> these 5′ ends are genuinely capped — which per-gene RACE cannot. That remains the lead
> recommendation in `COLLABORATOR_UPDATE.md`.

---

## Assay 5 — Is the 88-aa protein made?

The single experiment that separates the leading hypotheses (deep dive §5). Cheap.

1. Lyse paired ± LPS macrophages in RIPA with protease inhibitors.
2. **A 15–16% gel or a Tris-tricine gel.** A standard 10% gel runs a 10 kDa product off the
   bottom — this is the most common way to miss the band.
3. Blot with an **N-terminal / full-length ATF5 antibody**. A C-terminal antibody cannot
   detect this product by construction, because the epitope is in the lost bZIP region. Check
   the datasheet's immunogen before ordering.
4. Expected: full-length ATF5 ≈32 kDa (282 aa) and, if translated, a second band ≈10 kDa
   (88 aa).

**Interpretation.** A ~10 kDa band supports a protein-level hypothesis (H2/H3) and justifies
Assay 6. No band, with full-length clearly detected and a high-percentage gel, supports the
transcript acting as a translational-control device (H1) — in which case the isoform's role is
the diversion of *ATF5* output, and no further protein work is warranted.

---

## Assay 6 — Functional test (only if Assay 5 shows the band)

The literature gives a directional prediction: in BV-2 microglia, siRNA knockdown of ATF5
**increased** IL-6, IL-1β and TNF-α (PMID 36834738, doi:10.3390/ijms24043322). ATF5 is a brake
on inflammatory output acting through the mitochondrial UPR.

Express the 88-aa form in wildtype THP-1 (and separately knock down the endogenous isoform with
an shRNA targeting the 93-bp read-through region, which is the only isoform-selective target
available), then measure:

- IL-6, IL-1β, TNF-α by ELISA and RT-qPCR at 4 h LPS;
- UPRmt targets HSPD1, HSPA9, LONP1, CLPP by RT-qPCR;
- full-length ATF5 protein, to check the construct is not simply reducing it.

**Prediction if the isoform limits functional ATF5:** more IL-6/IL-1β/TNF-α and blunted UPRmt
target induction. Do not test CEBPB/CEBPD sequestration as the primary mechanism — the
published dominant-negative ATF5 retains the leucine zipper, and this isoform does not
(deep dive §5, H3).

---

## Oligonucleotides

All 16 oligos are in `results/tables/atf5_oligo_order_table.csv` with sequence, length, Tm, GC,
genomic coordinates, product size and specificity result. Order at **25 nmol, desalted**;
nothing here needs HPLC purification.

Every *ATF5* oligo was screened against **all 315,349 transcripts** in the UT long-read
annotation and every position verified back against the spliced cDNA. **Zero off-target loci**
for all 14 *ATF5* and reference oligos except YWHAZ, which picks up two pseudogenes amounting
to 0.43% of locus output.

| oligo | 5′→3′ | assay |
|---|---|---|
| ATF5_D2_F | CACAAGGAGAGGGCTGGATC | 1 — structure, 514 bp |
| ATF5_D2_R | CACCTTCCACCTGCCCTTAC | 1 |
| ATF5_N0_F | AATCACAGGAGGAAGCAGGC | 2 — numerator, 187 bp |
| ATF5_N0_R | CAAAGGGGTCGAGGCTACAG | 2 |
| ATF5_T0_F | TCTCCTGTTGCCATCAGTGC | 2 — denominator, 119 bp |
| ATF5_T0_R | AGTCTACGAGCCATCCCAGC | 2 |
| ATF5_3RACE_GSP1 | TCGCCTTGGTGCCTGTCTTCGC | 3 — outer |
| ATF5_3RACE_GSP2 | GGGCTGGGATGGCTCGTAGACT | 3 — nested |
| ATF5_5RACE_GSP1 | GGGCACAAAGGGGTCGAGGCTA | 4 — outer |
| ATF5_5RACE_GSP2 | TGGCTCCGGGGTCTGGTTTCCT | 4 — nested |
| HPRT1_F | CCTGGCGTCGTGATTAGTGA | 2 — reference, 145 bp |
| HPRT1_R | TCACATCTCGAGCAAGACGT | 2 |
| TBP_F | GCACCACTCCACTGTATCCC | 2 — reference, 80 bp |
| TBP_R | CCCAGAACTCTCCGAAGCTG | 2 |
| YWHAZ_F | GCAGGCTGAGCGATATGATG | 2 — reference, 145 bp |
| YWHAZ_R | GACGACCCTCCAAGATGACC | 2 |

RACE GSPs are deliberately long and hot (Tm 65.8–67.7 °C) so they tolerate the touchdown
programme; the qPCR and RT-PCR primers all sit at 58.9–61.4 °C and share a 60 °C anneal.

---

## How to call a result

1. **Assay 1 negative** (no 514 bp band, controls working) → the transcript is an assembly
   artefact. Stop. Total cost: one RT and one PCR.
2. **Assay 1 positive, Assay 2 shows induction in sh-Luc but not in the knockdown** → the core
   claim reproduces. Proceed to 3 and 4 to pin the structure, and to 5 for the protein.
3. **Assay 1 positive, Assay 2 shows induction in both lines** → the transcript is real but
   UBL5-independent. Still a novel LPS-induced macrophage isoform worth reporting; the UBL5
   link is dropped.
4. **Assay 3 negative** (no intronic poly(A) product) → the 3′ end is an artefact and the
   coding prediction collapses, even if the 5′ structure holds.
5. **Any shift below the pilot-derived detection floor** → record as underpowered and
   uninformative, not as a negative. Add replicates.

A reminder on framing throughout: the interaction statistic that first surfaced this gene
carries **no genome-wide excess over the null** (4.05% of isoforms at nominal p < 0.05 against
5% expected). This candidate rests on its absolute LPS induction (p = 0.004), its presence in
all three primary donors, and the UPR^mt pathway convergence — not on the interaction q-value.

---

## Sources

- Upstream protocol: Thesis §2.5.1, §2.5.6, §2.6.1, §2.6.6, §2.6.7, §2.6.8, §7.2, Fig. 7.1.
- ATF5 as a brake on macrophage cytokine output: PMID 36834738, doi:10.3390/ijms24043322.
- ATF5 as the mammalian UPR^mt regulator: PMID 27426517, doi:10.1016/j.cub.2016.06.002.
- ATF5 uORF / eIF2α translational control: PMID 18195013, doi:10.1074/jbc.M708530200.
- ATF5 5′-UTR and NMD: PMID 23876217, doi:10.1111/febs.12440.
- Published dominant-negative ATF5 architecture: PMID 31676720,
  doi:10.1158/1541-7786.MCR-19-0631.
- UBL5/Hub1 and non-canonical splice-site usage: PMID 21614000, doi:10.1038/nature10143;
  PMID 28712727, doi:10.1016/j.molcel.2017.06.021.

Literature retrieved via PubMed.
