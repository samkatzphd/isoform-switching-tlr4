# Which assay can resolve which switch

Built 2026-10-22 by `scripts/26_ht_switch_anatomy.R` and `scripts/27_assay_assignment.R`.
Per-candidate table: `results/tables/assay_assignment.csv`.

## The finding that drives everything below

Of the 44 candidates surviving the cross-system and locus-fusion screens:

| what structurally differs | candidates |
|---|---|
| alternative first exon (TSS) | **29** |
| both ends differ | 7 |
| alternative last exon (3′) | 4 |
| alternative splice site | 2 |
| exon skipping | 1 |
| intron retention | 1 |

**Only 13 of 44 are qPCR-addressable** — not nested inside a sibling isoform *and* carrying
≥60 bp of candidate-only sequence. Ten candidates are nested outright. The median TSS
distance between the two isoforms of a pair is **8,680 bp**.

This is not bad luck, and `scripts/13_assayability_screen.R` already said so: an isoform that
differs only by where it *starts* is contained inside its longer partner by construction, so
there is no unique stretch to put a primer in. The external review hit the same wall and
eliminated 6 of 7 Q2 candidates on it.

**The previous plan — design qPCR primers for everything — would have covered 30% of the list
and silently omitted the dominant biology, including KDM6A and NCOA7, the two largest
effects.** Nesting disqualifies one assay, not the candidate.

## Resulting primary assignment

| primary assay | candidates |
|---|---|
| 5′ RACE | 28 |
| unique-region qPCR/ddPCR (+ 5′ RACE to fix the TSS) | 8 |
| 3′-region qPCR/ddPCR | 3 |
| junction-spanning qPCR/ddPCR | 3 |
| 3′ RACE | 1 |
| intron-spanning qPCR/ddPCR | 1 |

Consequence-driven routes, layered on top: **2** PTC-positive rising isoforms (NMD
inhibition), **35** protein-changing by ORF length (isoform-resolving western or targeted MS),
**8** with unchanged ORF length (ribosome/polysome profiling), **10** nested (targeted
long-read is the only way to see them as a distinct molecule).

## What each assay measures, and what it cannot

**qPCR / ddPCR over a unique region.** Relative (qPCR) or absolute (ddPCR) abundance of one
isoform. Cheap, quantitative, the only option here that scales to many samples and timepoints.
It reads *abundance, not structure*: it confirms the isoform is present and moves, and says
nothing about where it starts. **Prefer ddPCR for these effect sizes** — most survivors sit at
|dIF| 0.15–0.35, and a 2^−ΔΔCt ratio resting on a shifting denominator is fragile in exactly
that range.

**Junction- or intron-spanning qPCR.** Same readout, but the discriminator is a junction
rather than a stretch of sequence — appropriate for exon skipping and alternative splice
sites. An intron-spanning pair is the standard readout for intron retention, and gives the
retained and spliced forms as two product sizes.

**5′ RACE.** Sequences the 5′ end directly, returning the actual TSS rather than inferring it.
*The* assay for an alternative first exon, which is two-thirds of this list. Low throughput and
only semi-quantitative; template-switching chemistries can generate artefactual ends, so a
known-TSS positive control and a no-RT control are not optional.

**CAGE / nanoCAGE.** Genome-wide, cap-dependent TSS quantification. Not per-candidate and not
cheap, but it would answer the entire shortlist at once *and* settle the standing
cap-selection question — whether these 5′ ends are real caps or assembly truncation. For a
list that is 29/44 promoter switches, **this is arguably the single highest-value experiment
available**, and it would retire a caveat that currently qualifies every 5′ claim in the
project.

**Targeted long-read cDNA (PacBio or ONT amplicon).** Amplify the full transcript with
flanking primers, sequence the whole molecule. The universal fallback: it resolves nested
isoforms, multi-difference pairs, and phasing that no short amplicon can reach. Caveats that
matter here: quantification is crude, and long-template PCR biases toward shorter products —
which is awkward when the shorter isoform is often the one being measured. Use it to establish
*that* the isoform exists with the claimed structure, then quantify by another route.

**3′ RACE / 3′-end sequencing (QuantSeq, PAS-seq).** The 3′ equivalent, for the 4 alternative
last-exon candidates.

**NMD inhibition (UPF1 knockdown, or cycloheximide).** For the 2 PTC-positive rising isoforms:
if the isoform is an NMD target, blocking NMD should raise it specifically relative to its
partner. Tests the predicted *fate*, not the structure. CHX is blunt and pleiotropic; UPF1
knockdown is cleaner but slower and has its own transcriptome-wide effects.

**Isoform-resolving western or targeted MS.** For the 35 protein-changing switches. **Requires
an antibody whose epitope lies in the differing region** — N-terminal for a 5′ switch,
C-terminal for a 3′ switch. I cannot validate antibodies or confirm availability; these are
specified with the epitope requirement stated, for a human to check against catalogues.
Targeted MS (PRM/SRM on a discriminating peptide) avoids the antibody problem but needs the
peptide to be unique and detectable.

**Ribosome / polysome profiling.** For the 8 switches with unchanged ORF length, where no
protein length changes and the plausible consequence is translational efficiency rather than
composition. Note the caveat below.

**Minigene splicing reporter.** Asks whether the event is *cis*-encoded and LPS-responsive,
rather than merely correlated with LPS. The right orthogonal test for the splicing-class
events.

**Alternative-promoter luciferase.** For promoter switches: is the alternative promoter itself
LPS-inducible?

**CRISPRi or promoter deletion.** The functional test — silence the alternative promoter and
ask whether the locus's LPS response changes.

## Caveats carried in the table

- **ORF comparison is by LENGTH only.** `orfTransciptStart` is transcript-relative, so two
  isoforms with different 5′ ends have different ORF starts even when the CDS is identical —
  comparing starts made "ORF identical" zero for all 44, which is impossible with 29
  alternative-first-exon events. Equal length implies equal protein *length*, not equal
  *sequence*; confirming identity needs the translated ORF, which the context layer does not
  carry. Treat the 8 "same ORF length" candidates as *candidates* for a UTR-only
  interpretation, to be confirmed before committing to ribosome profiling.
- **PTC is absent for 3 of 44** candidates, so the NMD route may be under-called.
- **Event classification comes from exon coordinates in the HT assembly**, so it inherits
  whatever the assembly got wrong — which is exactly why the fusion screen ran first.
- **The cap-selection question is unresolved and bears directly on 29 of these 44.** If the
  alternative first exons are assembly truncations rather than real caps, those switches are
  artefacts. 5′ RACE answers it per candidate; CAGE answers it for all of them at once.
