# Next session — start here

Updated 2026-08-22 (second pass, drive mounted). Delete or rewrite once the task below is done.

## The task: primer design + empirical specificity for the 25 survivors

`scripts/13_assayability_screen.R` has run. It says **where** on each surviving isoform a
discriminating primer can sit; it has **not** designed one, and it has **not** run the
empirical specificity match against all 315,349 transcripts. That match is what confirmed the
single surviving Q2 assay at 1-of-315,349, and it is the remaining gap.

### Start here

`results/tables/q1_assayability.csv` — the 25 rows with `assayable == TRUE`. Each carries
`unique_region` as `chr:start-end`, so design can start from a coordinate.

**Do the five both-sides genes first**, since only these let the switch be measured *as a
switch* rather than inferred from one side:

| gene | list | longest unique (bp) | note |
|---|---|---|---|
| **RBKS** | response | 273 | best candidate — dIF −0.33 and +0.44, both sides designable |
| KDM6B | response | 147 | both sides |
| ATRN | cryptic | 4,396 | both sides, lots of room |
| PLD6 | cryptic + response | 1,654 | both sides |
| SAMD14 | cryptic | 104 | both sides, tightest |

### Method to follow

The review's step 3, in `docs/external_review/CANDIDATES.md` and `tables/primer_designs.csv`:
design into the unique interval, then match **both** primers against every transcript in
`fa/transcripts_U-T.fa` and count amplified products. Do not trust predicted specificity —
the review's 5′-anchored alternatives amplified 4–10 transcripts each despite looking fine.

A new `scripts/14_primer_design.R` following the standard skeleton. Guard the FASTA read the
way `13` guards the GFF3: warn and skip, do not fail.

## Done since the last note

- **`12_dtu_refit_2x2.R`** — the blocking step. Interaction test: **2 genes of 13,086**. The
  174-vs-37 marginal gap between arms does not survive a direct test. Retention **≥56%** on a
  single estimator, matching the expression figure, not the retracted splicing one. Q2 stays
  *unresolved* (underpowered), not rejected. `REVIEW_CHANGES.md` §0h.
- **`13_assayability_screen.R`** — 25 of 63 Q1 isoforms designable, 5 genes both-sides.
  Validated against the review's Q2 verdicts (all five checkable isoforms reproduce exactly).
- Report and artifact both updated with these two sections.

## State

- Branch `fix/review-corrections-2026-07`, pushed, **~39 commits ahead of `main`, still
  unmerged**. This is the longest-standing loose end in the project.
- Artifact (private): <https://claude.ai/code/artifact/2a915114-23fa-4b47-8d5d-15ed4ce6e674>
- External drive was mounted for this pass; `99_verify_inputs.R` exits 0 on all 12 inputs.
  Primer design needs only local `fa/` and `gff3/`, so it runs without the drive.

## Queued: PTC / IR annotation of the hidden-layer sets — NO DRIVE NEEDED

Deferred by request 2026-08-23. Reading the upstream scripts (`NCBR-40-main`, §0i) showed
`analyzeORF` **was** run, and its output is already in committed files:
`isoformFeatures_T_UT.rds` carries `PTC` (192 TRUE / 1,465 FALSE / 112 NA) and `IR` (293
isoforms with at least one retained intron). Nothing needs the external drive.

This is the annotation Tier-2 open question #4 asks for. The informative version, per
`HANDOFF.md` §"Open work": annotate the **isoforms**, not the genes — are the lost minor
isoforms NMD targets (`PTC == TRUE`) or intron-retaining, and are the retained dominant ones
full-length coding? That converts Class B "purification" from a pattern into a mechanism.

Apply to: `hidden_layer_cryptic_switchers`, `hidden_layer_compositional_rescue`, and the
Class B set. Watch the 112 NAs — `PTC` is undefined where no ORF was called, which is not the
same as "not a target".

**Separately, and this one DOES need the drive:** `analyzeAlternativeSplicing` was never run
(gated behind `--run_extra_analysis`, which none of the four README commands pass). Re-running
with that flag gives ISA's native ATSS classification — an independent cross-check on the
21.7% pure-promoter figure from the GFF3 work.

## New question for the experimentalist, added 2026-08-23

**Was the PacBio/IsoSeq library cap-selected?** (TeloPrime, or IsoSeq with 5'-cap
verification.) 23,044 of 166,394 quantified transcripts — **13.8%** — are gffcompare class
`c`, contained inside another isoform. Those are either genuine alternative-internal-promoter
isoforms or 5'-incomplete assemblies, and without cap selection the second class is produced
in quantity.

It bears on the biology, not just assay design: 6 cryptic candidates (PCSK7, CCDC117,
SPRING1, DCUN1D4, SORT1, ITSN1) rest on `c`-class transcripts **alone**. All candidate tables
now carry `class_codes` / `n_switching_contained` so the exposure is visible per gene.

Mild reassurance already checked: the nested forms do not move as a class (8 rise, 5 fall,
p = 0.58), so there is no global 5'-coverage signature. Not conclusive.

## Still waiting on the wet-lab side

Unchanged, and now the main brake on Q2. All three are in the lab meeting update §6: the
UBL5 allele (full KO vs hypomorph — cheapest, reframes everything), polyA vs ribo-depleted
libraries (intron retention is the diagnostic readout and polyA captures it poorly), and the
LPS dose / timepoint / PMA protocol.

**One thing the refit changed about priorities.** Q2 is now limited by *power*, not by
analysis method — the interaction test is real but underpowered at n = 3. No further
reanalysis will settle it. More replicates or a larger perturbation will.
