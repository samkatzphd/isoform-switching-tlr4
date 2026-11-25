# Next session: start here

Updated 2026-11-24. Rewritten from the commit history and the Oct 23-25 docs; the previous
version (2026-08-27) predated scripts 21-35 and the whole assay-design phase, so anything it
said about "primer design not done" or the RBKS/ATRN/IRAK3 shortlist is stale.

## State of the repo

- `main` now carries everything. `fix/review-corrections-2026-07` was fast-forwarded into it
  on 2026-11-24 (70 commits, no conflicts), together with the Oct 23-25 ATF5 and UBL5 work
  that had been left uncommitted. Work on `main` from here.
- `core.fileMode` is set to `false` for this repo: something (Claude Science or the folder
  mount) had flipped 33 files to executable with no content change.
- Four large UT tables are on disk but ignored (`ut_q1_power_audit.csv.gz`,
  `ut_tx2gene.tsv`, `ut_reference_audit.csv.gz`, `ut_q2_genelevel_KD_vs_WT.csv.gz`), as is
  `counts/` (295 MB). None of these has a second physical copy yet.

## Where the work stopped (2026-10-25)

Two arms, at different stages.

### HT arm (WT +/- LPS): candidates designed, ready for the bench

- `25` screened the 63-gene shortlist for fused loci in HT space: 19 disqualified, 44 survive.
- `27` routes each candidate to the assay that can resolve it (28 go to 5' RACE).
- `28`-`30` designed 942 oligos, searched every one against all 296,613 HT transcripts, and
  selected on the measured binding. `28` also designs a total-gene pair per gene, so the
  isoform *fraction* is measurable and not just the isoform's Ct.
- `32`/`33` write `docs/ASSAY_DESIGNS.md` with a platform call per candidate (28 of 44 need
  long reads).
- `docs/COLLABORATOR_UPDATE.md` is the bench-facing summary: first batch of 10 (KDM6A, NCOA7,
  BCL2L13, RSAD2, MERTK, OR52K3P, ZC3HAV1, URM1, SPHK1, CHMP5), 70 oligos.
- `docs/PROTOCOLS.md` reproduces the thesis conditions with section citations. The lab qPCR
  is probe-based, the designs here are primer-only; 3.2 gives both routes.

### UBL5 arm (U vs T): one candidate, and a postmortem

Read in this order:

1. `docs/UBL5_HANDOVER.md`: decisions already made. Work in **UT** space (corrected in the
   last commit of Oct 23; the HT-space recommendation was wrong). **U is an shRNA knockdown,
   not a knockout.**
2. `docs/UBL5_KD_EFFICACY_NOTE.md`: the knockdown is ~30% at the mRNA level.
3. `docs/UBL5_ARM_POSTMORTEM.md`: why the genotype x treatment interaction yields almost
   nothing (2 of 13,086 genes; 1 survives QC) and what is still worth testing.
4. `docs/ATF5_DEEP_DIVE.md` then `docs/ATF5_Q1_EXPERIMENTAL_PLAN.md` and
   `docs/ATF5_ASSAY_PROTOCOL.md`: the one survivor, a novel ATF5 isoform (TCONS_00159376).
   The deep dive §1 reframes the plan: lead with the within-WT LPS induction, which is
   internally controlled, not the genotype contrast.

## Loose ends, roughly in order

1. **Knockout vs knockdown.** `CLAUDE.md`, `docs/AGENT_CONTEXT.md` and `README.md` (two
   places) still say "knockout". Fix, and paste the updated AGENT_CONTEXT into the project
   instructions.
2. **"Q1" means two different things.** In `CLAUDE.md`, Q1 is WT +/- LPS in HT. In the Oct 25
   docs and `ut_q1_*` tables, "Q1" is the genotype x treatment interaction in UT. Pick one
   naming and apply it.
3. **ATF5 provenance gap.** The ATF5 tables have no generating script under `scripts/`; only
   `explore/ubl5_grounding_explore.R` exists. The plan also cites
   `ut_q1_atf5_pair_discrimination.csv`, `ut_q1_atf5_primer_specificity.csv` and
   `ut_q1_atf5_candidate.png`, which are not in the repo. Regenerate them from a numbered
   script before anyone orders the ATF5 oligos.
4. **Scripts 25-35 hardcode HT paths.** The method is verified; the reference needs
   parameterising before the UT arm can reuse it.
5. **Questions only the lab can answer:** UBL5 protein level in the knockdown (governs how
   every negative is read); whether the sequenced WT is untransduced or sh-Luciferase; the
   Iso-Seq chemistry version. Library type is settled: TruSeq Stranded mRNA, polyA.
6. **Housekeeping:** `Rplots.pdf` is still tracked despite its ignore entry
   (`git rm --cached Rplots.pdf`). The two status artifacts linked in the August version of
   this file have not been updated since August.

## Still standing from before

- Q2 is limited by power, not analysis method. No further reanalysis will settle it.
- 29 of the HT candidates are TSS switches; CAGE would answer all of them in one experiment,
  and whether the libraries were 5'-cap verified still bears on most of the list.
