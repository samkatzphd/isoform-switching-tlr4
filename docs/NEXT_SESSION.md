# Next session — start here

Written 2026-08-22. Delete or rewrite this file once the task below is done.

## The task: assayability screen on the two Q1 candidate lists

`scripts/11_q1_candidates.R` now produces two Q1 candidate lists. **Neither has been through
the nesting / primer-specificity screen** that eliminated 6 of the 7 Q2 candidates. That is
the next job.

Everything it needs is on local disk. **The external drive is NOT required.**

### Inputs, all present

| What | Where | Note |
|---|---|---|
| Candidate lists | `results/tables/q1_response_candidates.csv` (20)<br>`results/tables/q1_cryptic_candidates.csv` (37) | from script 11 |
| Full ranked union | `results/tables/q1_candidates.csv` (157) | if you want a wider net |
| Exon structures | `gff3/U_T.gff3` | 1.1 GB, gitignored, already local |
| Transcript sequences | `fa/transcripts_U-T.fa` | 2.8 GB, gitignored, already local |
| Per-isoform detail pattern | `data/processed/isoformContext_T_UT.rds` | `isoform_id`, `dIF`, `IF1`/`IF2`, `class_code` |

Note the lists are **gene**-level. The screen is **isoform**-level, so step one is expanding
each gene to its switching isoforms plus their rising/falling partners — the same rule
script 11's provenance describes and the review used: `(is_switching OR |dIF| > 0.15) AND
expressed`. You cannot validate a switch by measuring only one side of it.

### What the Q2 screen did, to replicate

From `docs/external_review/CANDIDATES.md` (§"assay feasibility") and
`docs/external_review/tables/primer_designs.csv`:

1. **Nesting test.** For each candidate isoform, ask whether every exon and every junction it
   contains is also present in a longer isoform of the same gene. If so it is **fully
   nested** and *no primer pair can distinguish it* — 5 of the 7 Q2 candidates died here.
2. **Unique-sequence test.** Even when not nested, an isoform can be tiled over by the
   union of the others. Require a contiguous unique stretch of **≥ 60 bp**. SOCS4
   `TCONS_00086549` failed this despite a distinct structure.
3. **Empirical primer specificity.** Design against the unique region, then match both
   primers against **every** transcript in `fa/transcripts_U-T.fa` and count how many are
   amplified. Do not trust predicted specificity — the review's one surviving assay (IRAK3
   `TCONS_00065141`) was verified this way at 1 of 315,349 transcripts, while the
   5′-anchored alternatives amplified 4–10 each.

### What to expect

Expect the screen to remove several candidates, and expect the 5′-end problem to recur: if
~22% of events are pure promoter switches and much of the rest is 5′ truncation
(`event_structure_summary_gff3.csv`), nested isoforms are the **norm** here, not the
exception. If most of List A also fails, the honest conclusion is the same one the review
reached for Q2 — **5′ RACE, CAGE, or targeted long-read cDNA**, not qPCR — and long-read
cDNA resolves everything at once.

### Where to put it

A new `scripts/12_assayability_screen.R`, following the standard skeleton (copy the header
of `11_q1_candidates.R`). Write `results/tables/q1_assayability{,_primer_designs}.*` via
`write_table_pair()` and end with `write_run_manifest()`. Guard the GFF3/FASTA reads the way
`04`/`06` guard the drive: **warn and skip, do not fail**, since those two files are
gitignored and will be absent on a fresh clone.

## State as of this note

- Branch `fix/review-corrections-2026-07`, pushed, **35 commits ahead of `main`, still
  unmerged**. Working tree clean.
- Reports: `reports/lab_meeting_update.{qmd,html}` is the current whole-project document and
  covers both Q1 lists. `seminar_summary` is kept for its `00`–`07` detail and points at it.
- Published artifact (private):
  <https://claude.ai/code/artifact/2a915114-23fa-4b47-8d5d-15ed4ce6e674> — republish to the
  same URL after the screen so the candidate sections stay in sync.
- R is **4.4.1** with ISA 2.4.0; `satuRn`, `DEXSeq`, `edgeR`, `limma` present, `DRIMSeq`
  absent. A transient R 4.5.3 that wrote a null-SHA manifest entry in August is gone;
  `write_run_manifest()` now warns if that recurs.

## Still blocked on the external drive

The model-based 2×2 refit (`OPEN_QUESTIONS.md` #1) — the step that would settle the
splicing-specific UBL5 claim — needs the RSEM counts at
`/Volumes/Expansion/IsoformSwitchAnalyzer/Counts/`. Unchanged by anything above.

## Still waiting on the wet-lab side

Three answers gate the validation design, all listed in the lab meeting update §6: the UBL5
allele (full KO vs hypomorph — cheapest, and it reframes everything), polyA vs ribo-depleted
libraries (intron retention is the diagnostic readout and polyA captures it poorly), and the
LPS dose / timepoint / PMA protocol.
