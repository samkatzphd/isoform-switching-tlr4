# Next session — start here

Updated 2026-08-27 (drive mounted). Delete or rewrite once the task below is done.

## The task: primer design + empirical specificity, on the CORRECTED shortlist

`scripts/13_assayability_screen.R` says **where** on each surviving isoform a discriminating
primer can sit. It has not designed one, and it has not run the empirical specificity match
against all 315,349 transcripts. That match is what confirmed the single surviving Q2 assay at
1-of-315,349, and it is still the remaining gap.

**Read `results/tables/candidate_profiles.csv` first, not `q1_assayability.csv`.** The
shortlist changed this session — see "What changed" below.

### The order to work in

| gene | why it is here | unique bp |
|---|---|---|
| **RBKS** | biggest both-sides effect (dIF 0.44); loses the catalytic proton acceptor and 13 of 14 binding sites; rising isoform is PTC+; **zero** isoform-level prior art | 273 |
| **ATRN** | membrane→secreted switch, loses only the TM helix (1280–1300); the secreted/membrane pair was described in PMID 10811918, so there is a predicted answer to check | 4,396 |
| **IRAK3** | loses the Death domain, Ser110, the PIN1 site and 78% of the pseudokinase domain | 3,465 |
| **KDM6B** | both sides, dIF 0.30 — but the rising isoform has **no ORF called**; inspect it manually before designing | 147 |
| **NCOA7** | best single-sided candidate: internal-promoter isoform gaining the TLDc domain, and it is the known interferon-inducible NCOA7-B (PMID 30478388) | 81 |

**PLD6 is off the list** — see below. It was previously ranked as the second-best both-sides
gene, so anything written before 2026-08-27 that names it is stale.

### Method to follow

Unchanged: design into the `unique_region` interval, then match **both** primers against every
transcript in `fa/transcripts_U-T.fa` and count amplified products. Do not trust predicted
specificity — the review's 5′-anchored alternatives amplified 4–10 transcripts each despite
looking fine. A new `scripts/19_primer_design.R` following the standard skeleton (14–18 are
taken). Guard the FASTA read the way `13` guards the GFF3: warn and skip, do not fail.

## What changed 2026-08-27

### Seven candidates are not switches at all

Two isoforms of one gene always share some exonic sequence. Sixteen candidate pairs share
**none** — zero exonic bases, sides 5.9–136 kb apart — because the assembler merged two
neighbouring genes into one `XLOC`. Confirmed against Ensembl in `18`:

`PLD6`(=PLD6+FLCN) · `DDX58`(=TOPORS+RIGI) · `IPCEF1`(=IPCEF1+CNKSR3) ·
`THAP2`(=THAP2+TMEM19) · `MREG`(=PECR+MREG) · `IFITM2`(=IFITM1+IFITM2) · `CCL3`(=CCL3L3+CCL3)

For these, dIF is the relative expression of two genes. **PLD6 was a both-sides priority
gene; MREG and THAP2 were in the assayable 25.** `NCOA7` is the one suspect that clears —
both sides are `ENSG00000111912`, a real internal-promoter pair.

### Class B is closed (OPEN_QUESTIONS #4)

Carrying `PTC` + `orfAnalysis` from the unfiltered exports into the context layer took Class B
coverage from 5–21% to 86–93%. The purification reading fails three ways: the minor-vs-dominant
PTC asymmetry is universal (background OR 75–319), the Class B minor-isoform PTC *rate* is not
enriched (depleted in T_HT/T_UT, null in H_HT/U_UT), and PTC does not predict which minors are
shed (null in all four). The first version of that last test looked like a clean four-of-four
result in the opposite direction and was a pseudocount artifact — read the comment in `14` §4c
before touching it.

### Four new scripts

- `15_isoform_characterization.R` — pairwise structure, reading frame, coding consequence,
  splice sites. Offline from `gff3/` + `fa/` + the context layer.
- `16_functional_sites.R` — maps each switch onto UniProt features. Cached, offline after
  first run.
- `17_candidate_literature.R` — isoform vs immunity vs both, per gene. Cached.
- `18_candidate_profiles.R` — one row per gene joining all of the above.

### What the candidate set actually looks like

- **54 of 63 pairs differ by a distinct first exon.** These are promoter switches, not
  splicing switches: 14 exon skips, 5 alt donors, 5 alt acceptors, 1 intron retention. That
  sharpens the cap-selection question rather than settling it.
- **14 of 63 change no protein at all** (UTR only): AGPAT3, R3HCC1L, RGS19, RSAD2, SAMD14,
  SOCS4, SPRING1, TMEM63B, ZDHHC21. SOCS4 and SPRING1 were named review candidates.
- **All 458 differential junctions are canonical** (452 GT-AG, 6 GC-AG).
- **12 candidates are `c`-class AND differ only at the 5′ end** — indistinguishable from a
  5′-incomplete assembly until the cap-selection question is answered.

## State

- Branch `fix/review-corrections-2026-07`, **~47 commits ahead of `main`, still unmerged**.
  This remains the longest-standing loose end in the project.
- Profile artifact (private): <https://claude.ai/code/artifact/a08054fc-2b30-43f5-a4aa-f9470ab79cfa>
- Project-status artifact (private): <https://claude.ai/code/artifact/2a915114-23fa-4b47-8d5d-15ed4ce6e674>
  — **not yet updated with any of the above.**
- The `reports/` Quarto tree does not yet carry sections 15–18 either.

## Still waiting on the wet-lab side

Unchanged, and still the main brake on Q2. All in the lab meeting update §6: the UBL5 allele
(full KO vs hypomorph — cheapest, reframes everything), polyA vs ribo-depleted libraries, and
the LPS dose / timepoint / PMA protocol. Plus the cap-selection question from 2026-08-23,
which this session made more load-bearing, not less: with 54 of 63 candidate switches turning
on where the transcript starts, whether the library was 5′-cap-verified now bears on most of
the shortlist.

**Q2 is limited by power, not by analysis method.** No further reanalysis will settle it.
