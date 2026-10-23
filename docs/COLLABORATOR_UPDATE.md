# Isoform switching in the LPS response: candidates, assays, and what to run

Prepared from the long-read (PacBio) isoform analysis of THP-1 macrophages and primary human monocyte-derived macrophages, ± LPS. This is the bench-facing summary; the complete per-candidate reference is `docs/ASSAY_DESIGNS.md`, and every number below comes from the tables in `results/tables/`.

## The short version

We have **44 candidate isoform switches** that survived a direction-of-effect filter and a locus-fusion screen. Of these, **28 are visible only because we used long reads** and 16 could in principle have been found by standard short-read RNA-seq.

The central practical finding is that **most of these are not qPCR problems.** 29 of the 44 switches are alternative-first-exon events — the two isoforms differ in where transcription *starts*, and share everything downstream. For those, the shorter isoform is contained inside the longer one, so no short amplicon can be specific to it. They need an assay that reads the 5' end: RACE, or CAGE.

I am proposing a **first batch of 10 candidates** (below) rather than all 44, chosen by an explicit rule: clean reagents, no quality flags, and either reproducible across both cell systems or supported by a complete qPCR assay. The other 34 are ready in the same format whenever you want them.

## 1. The candidates

**dIF** is the change in that isoform's share of its gene's total output on LPS. A dIF of +0.30 means the isoform went from, say, 20% to 50% of the gene. It is a *fraction*, not an expression change — which is why every qPCR below needs a total-gene denominator as well as an isoform-specific amplicon.

| gene | dIF | direction on LPS | seen in | needs long reads? | what differs | assay |
|---|---|---|---|---|---|---|
| **KDM6A** | +0.705 | up | both systems | yes | alternative first exon | 5' RACE |
| **NCOA7** | +0.651 | up | both systems | yes | alternative first exon | 5' RACE |
| **BCL2L13** | +0.391 | up | both systems | yes | alternative first exon | 5' RACE |
| **RSAD2** | +0.341 | up | THP-1 (T) | yes | alternative first exon | unique-region qPCR / ddPCR |
| **MERTK** | +0.307 | up | both systems | yes | both ends differ | unique-region qPCR / ddPCR |
| **OR52K3P** | -0.245 | down | Primary MDM (H) | yes | alternative first exon | unique-region qPCR / ddPCR |
| **ZC3HAV1** | -0.211 | down | both systems | yes | both ends differ | unique-region qPCR / ddPCR |
| **URM1** | -0.201 | down | Primary MDM (H) | yes | exon skipping | junction-spanning qPCR / ddPCR |
| **SPHK1** | +0.179 | up | THP-1 (T) | no | alternative first exon | unique-region qPCR / ddPCR |
| **CHMP5** | +0.166 | up | Primary MDM (H) | yes | alternative splice site | junction-spanning qPCR / ddPCR |

Direction is the prediction you are testing. *Up* means that isoform's share of the gene should be higher in the LPS-treated sample.

**KDM6A, NCOA7, BCL2L13, MERTK, ZC3HAV1** switch in both THP-1 and primary MDM with all three donors agreeing in direction. Those are the ones I would bet on.

## 2. Which assay, and why that one

| assay | candidates in batch 1 | what it reads out | why this one |
|---|---|---|---|
| 5' RACE (nested) | 3 | the exact 5' end of each transcript, and the ratio between two start sites in one reaction | the two isoforms differ only in where they start, so the short one is contained in the long one and no amplicon can be specific to it |
| Isoform-specific qPCR / ddPCR | 7 | the abundance of one isoform, expressed as a fraction of the gene's total | these have a stretch of sequence unique to one isoform that a primer can sit in — the quickest readout in the set |
| Targeted long-read cDNA | 10 | the whole transcript as one molecule, so every difference is read at once | confirmatory for everything here, and the only route at all for an isoform fully nested inside a sibling |

**The single highest-value experiment is CAGE (or nanoCAGE).** It maps transcription start sites genome-wide, so one experiment per condition would answer all 29 alternative-first-exon switches at once, instead of running that many separate RACE reactions. It would also settle a standing question about whether these alternative 5' ends are genuinely capped. If there is any route to CAGE through a core facility or a collaborator, that is where I would put the effort first; the RACE designs below are the per-candidate fallback and are worth running regardless for the handful we care most about.

## 3. What to order

### Oligos

**70 oligos** for the 10 batch-1 candidates, listed in `results/tables/order_batch1.csv` in a form you can paste straight into an order. Standard desalted, 25 nmol scale is fine for all of them — these are 20–25-mers and none needs HPLC purification.

Every oligo has been checked against all 296,613 transcripts in our reference by exact match on its 3'-terminal 15 nt, extended to full length. In the batch-1 set, every primary-assay oligo binds only within its own gene.

Two columns in the order sheet report specificity and they mean different things. `offtarget_loci` counts binding sites for that **single primer**. The `selection` column reports the verdict for the **pair**, which only amplifies where both primers bind in the right orientation and spacing — so a primer with one stray site can still belong to a pair that amplifies nothing else. Where they appear to disagree, the pair verdict is the one that governs whether you get a spurious band.

Each candidate gets up to four kinds of oligo, and they do different jobs:

| role | job |
|---|---|
| `qPCR_F` / `qPCR_R` | isoform-specific amplicon — the numerator |
| `qPCR_total_F` / `qPCR_total_R` | sits in sequence common to both isoforms — the denominator, so you can report a fraction rather than an expression change |
| `RACE5_GSP1_outer` / `RACE5_GSP2_inner` | nested antisense primers for 5' RACE; outer for the first PCR, inner for the nested PCR |
| `longread_F` / `longread_R` | outer pair spanning the whole transcript, for the confirmatory long-read amplicon |

### Kits and reagents

- **RNA prep**: your usual column prep plus an on-column or in-solution DNase step. The no-RT control below is the check that this worked.

- **Reverse transcription**: a standard RT kit for the qPCR work. For 5' RACE you need a *template-switching* RT kit specifically — the SMARTer RACE 5'/3' kit (Takara) is the common one, or any equivalent that adds a defined adapter at the 5' end. A conventional RT kit will not work for RACE.

- **qPCR**: SYBR-based master mix. If you would rather have absolute numbers and skip standard curves, ddPCR with EvaGreen works and is arguably better suited to a ratio measurement.

- **Gel / cleanup**: 1.5–2% agarose for sizing RACE products, plus a gel-extraction kit — you will want to sequence the bands, not just photograph them.

- **Sequencing**: Sanger for individual gel-extracted RACE bands. For the long-read confirmation, a PacBio or ONT amplicon run through a core facility; the amplicons are 0.5–9 kb.

### What to ask a core facility for

- CAGE or nanoCAGE on the same RNA, ± LPS (see above — highest value per experiment).

- PacBio or ONT amplicon sequencing of the long-read outer products, barcoded by condition.

## 4. What to run

### Stimulation and RNA — read this first

> **The validation has to use the same stimulation conditions as the sequencing.** The PMA differentiation protocol, the LPS dose and the timepoint are not recorded anywhere in the analysis repository, so they need to come from the original experiment records before anyone starts. If the validation runs at a different dose or timepoint and a switch fails to reproduce, there is no way to tell whether the switch is wrong or the conditions were.

Three independent biological replicates per condition, −LPS and +LPS. For the primary-MDM candidates that means three donors. Process both conditions of a replicate in the same batch: these are ratio measurements, and a ratio compared across batches is not interpretable.

### Protocol A — isoform-specific qPCR (the quickest answer)

1. DNase-treat RNA; reverse-transcribe with random hexamers or oligo-dT.

2. For each candidate run **two** assays on the same cDNA: the isoform-specific pair and the total-gene pair.

3. Validate both pairs first: standard curve over at least 5 logs, efficiency 90–110%, and a single melt peak. **A fraction built from two pairs of different efficiency is not a fraction** — this step is not optional.

4. Compute the isoform fraction as 2^−(Ct_isoform − Ct_total) per sample. Do **not** normalise to a housekeeping gene: the claim is about the ratio of two isoforms within a gene, and total-gene normalisation is what makes it a fraction.

5. Compare the fraction between −LPS and +LPS across the three replicates.

**Controls**: no-RT and no-template for every pair; melt curve every run.

### Protocol B — nested 5' RACE (for the start-site switches)

1. Template-switching RT on 1 µg total RNA, both conditions, same batch.

2. First PCR: universal adapter primer + the **outer** gene-specific primer (`GSP1`).

3. Nested PCR: nested universal primer + the **inner** primer (`GSP2`), on a 1:50 dilution of PCR 1.

4. Run on 1.5–2% agarose. The per-candidate entries give the expected size from each isoform — for most of these the two start sites differ by kilobases, so the two products are far apart and easy to resolve.

5. **Gel-extract both bands and sequence them.** A size shift is suggestive; only the sequence tells you which start site was used. If you want a real ratio rather than band intensity, sequence the un-separated nested product as an amplicon library and count reads per start site — that is a much better measurement than densitometry.

**Controls**: no-RT (a band here means genomic DNA and voids the TSS call for that sample); no-template in both PCRs; a positive-control GSP against a single-TSS gene to confirm adapter ligation worked.

### Protocol C — targeted long-read cDNA (the confirmation, and the only option for some)

1. Amplify with the outer pair, both conditions, barcoded.

2. Check the product on a Bioanalyzer/TapeStation before sequencing — a dominant short product will eat the run.

3. Sequence on PacBio or ONT and count full-length reads per isoform.

4. Reads shorter than the expected amplicon are PCR or library artefacts, not short isoforms — do not interpret them as evidence for a new transcript.

This is the only assay that can see an isoform that is fully nested inside a sibling, and it is also the cleanest way to get a per-isoform ratio for anything in the list.

### How to call a result

For each candidate the per-candidate entry gives a specific threshold, but the rule is the same throughout:

- **Supports the switch**: the isoform's share moves in the predicted direction by at least 0.16 (THP-1) or 0.17 (primary MDM), across three replicates.

- **Falsifies it**: the share moves in the *opposite* direction.

- **Uninformative, not negative**: a shift smaller than that threshold. Those numbers are the minimum detectable change at n = 3 with 80% power, so a smaller result means the experiment could not have resolved it either way. The honest response is more replicates, not a negative conclusion.

## 5. What is not ready, and what I would not trust yet

- **Six of the 15 isoform-specific qPCR designs do not work as discriminating assays** and are excluded from batch 1. RBKS, TMEM177 and KDM6B amplify the partner isoform as well, so they measure the gene rather than the switch (TMEM177 gives the same 133 bp product from both). ZP3, PHLDA3 and POLD4 have primers that bind a second gene, and that survived all five alternatives I tried.

- **ACP3's long-read reverse primer binds 4,355 loci** — a repeat element, not a bad pick. Its qPCR assay is fine; only the long-read confirmation needs a hand-designed primer.

- **Two candidates have structural quality flags.** KDM6B's candidate transcript is a 219 nt model that shares only 72 nt with its 6,728 nt partner, which looks like an assembly fragment rather than a real isoform. A few others (JDP2, SPRED1, GATD3) share so little sequence with their partner that whether they are two isoforms of one gene is a fair question. All are flagged in the per-candidate entries.

- **Antibody and mass-spec routes are specified but not validated.** Where a per-candidate entry calls for an isoform-resolving western, it states where the epitope has to fall for the antibody to tell the two isoforms apart. Whether such an antibody actually exists has not been checked for any of them.

- **These oligos are designed and computationally screened, not bench-tested.** The specificity numbers are exact-match counts against a transcriptome, which is a much better check than a thermodynamic score, but it is not a gel.

