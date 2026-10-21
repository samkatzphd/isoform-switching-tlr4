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

| gene | platform | why it is here | unique bp |
|---|---|---|---|
| **RBKS** | Illumina-findable | biggest both-sides effect (dIF 0.44); loses the catalytic proton acceptor and 13 of 14 binding sites; rising isoform is PTC+; **zero** isoform-level prior art | 273 |
| **ATRN** | Illumina-findable | membrane→secreted switch, loses only the TM helix (1280–1300); the secreted/membrane pair was described in PMID 10811918, so there is a predicted answer to check | 4,396 |
| **IRAK3** | **PacBio required** | loses the Death domain, Ser110, the PIN1 site and 78% of the pseudokinase domain. Both sides map to `ENST00000261233` — Ensembl holds one transcript where PacBio resolves two | 3,465 |
| **KDM6B** | Illumina-findable | both sides, dIF 0.30 — but the rising isoform has **no ORF called**; inspect it manually before designing | 147 |
| **NCOA7** | **PacBio required** | best single-sided candidate: internal-promoter isoform gaining the TLDc domain, and it is the known interferon-inducible NCOA7-B (PMID 30478388) | 81 |
| **IL23A** | **PacBio required** | added 2026-10-21. Priority 7, assayable, frameshift consequence, rising side is a class-`j` novel junction structure. Was not on the previous list and is the most immunologically direct target on it | 75 |
| **WDR33** | **PacBio required** | added 2026-10-21. Priority 6, assayable, C-terminal change, 5,526 unique bp — the largest discriminable region on the list. A cleavage/polyadenylation factor, which bears on the promoter/cap story | 5,526 |

Platform column from `results/tables/candidate_detectability_by_switch.csv` (script 19).
"PacBio required" means no annotation-based short-read quantifier could have seen the switch
— either a side is a structure absent from Ensembl, or both sides collapse onto one Ensembl
transcript. Explore it with `explore/19_novelty_partition_explore.R`.

**PLD6 is off the list** — see below. It was previously ranked as the second-best both-sides
gene, so anything written before 2026-08-27 that names it is stale.

### Method to follow

Unchanged: design into the `unique_region` interval, then match **both** primers against every
transcript in `fa/transcripts_U-T.fa` and count amplified products. Do not trust predicted
specificity — the review's 5′-anchored alternatives amplified 4–10 transcripts each despite
looking fine. A new **`scripts/20_primer_design.R`** following the standard skeleton (14–19
are taken; 19 is the novelty partition added 2026-10-21). Guard the FASTA read the way `13`
guards the GFF3: warn and skip, do not fail.

**Design for both platform classes, and mark which is which on the order sheet.** For a
PacBio-required switch the isoform-specific assay is the *only* evidence available — no
short-read reanalysis of a public dataset can corroborate it — so those assays carry the
weight of the claim and need the stricter specificity check. For an Illumina-findable switch,
a public macrophage RNA-seq dataset can serve as an orthogonal check before any bench work.

## What changed 2026-10-21

### "Novel" was counting provenance, not discovery

`is_novel_pacbio` comes from the transcript ID prefix (`oId` starting `PB`), so it says which
assembler emitted the model — not whether the structure is new. 89 `PB` transcripts in the
T_HT object are class `=`, i.e. their intron chain exactly matches a named Ensembl transcript.
`scripts/19_novelty_partition.R` recounts it structurally. The published
`pct_switching_that_are_novel` overstates structural novelty by **34–43 percentage points**;
genuinely novel junction structures are **11–16%** of switching isoforms, not 46–58%. The
published figure is not retracted — it is mislabelled. Call it "isoforms contributed by the
PacBio assembly". Full detail in `REVIEW_CHANGES.md` §0j.

### The shortlist splits roughly in half by platform

Of 34 distinct candidate switches (locus fusions excluded): **16 Illumina-findable, 18
PacBio-required**. Nine of the 18 are PacBio-required because *both sides map to the same
Ensembl transcript* — a failure mode the isoform-level view cannot see at all. Six are
PacBio-required *and* assayable: **IL23A, IRAK3, NCOA7, RPUSD2, WDR33, ZDHHC21.** IL23A and
WDR33 were not on the previous priority list and should be.

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
