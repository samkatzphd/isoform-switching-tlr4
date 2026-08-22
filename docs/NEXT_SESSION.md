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

## Still waiting on the wet-lab side

Unchanged, and now the main brake on Q2. All three are in the lab meeting update §6: the
UBL5 allele (full KO vs hypomorph — cheapest, reframes everything), polyA vs ribo-depleted
libraries (intron retention is the diagnostic readout and polyA captures it poorly), and the
LPS dose / timepoint / PMA protocol.

**One thing the refit changed about priorities.** Q2 is now limited by *power*, not by
analysis method — the interaction test is real but underpowered at n = 3. No further
reanalysis will settle it. More replicates or a larger perturbation will.
