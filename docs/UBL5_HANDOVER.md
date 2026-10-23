# Handover: the UBL5 arm (U vs T)

Written at the end of the session that completed the HT assay-design chain, to start the
UBL5 analysis in a fresh session with clean context. Memory and artifacts are
project-scoped and carry over automatically; what follows is the part that does not.

---

## 1. What is being asked for

Two candidate lists, each taken all the way to collaborator-ready reagents and protocols —
the same end state the HT arm reached in `COLLABORATOR_UPDATE.md` + `PROTOCOLS.md`:

**Q1 — genotype-dependent isoform switching.** Genes where an LPS-driven change in isoform
usage is present in wildtype but absent in UBL5 knockdown, or present in knockdown and
absent in wildtype. This is a **genotype × treatment interaction**, not two separately-fit
analyses compared by eye.

**Q2 — genotype-dependent abundance.** Transcripts significantly different in abundance
between WT and KD **at baseline**, and separately **under LPS**. Two contrasts, reported
separately, not pooled.

Final deliverables, mirroring the HT arm:
- a screened candidate list with novel-PacBio vs Illumina-findable separated
- an assay per candidate, with designed and specificity-screened oligos
- a collaborator update and a bench protocol document
- a standalone R/dplyr/ggplot explorer per conclusion, reading only CSVs

---

## 2. Decisions already made — do not relitigate these

### Work in HT annotation space, not UT

This is the single most important thing in this document. **Both RSEM quantifications
contain all 18 libraries.** `counts/H_T_pipe/RSEM.isoforms.*.all_samples.txt` has
`U1_minus_S9_001` … `U3_plus_S14_001` alongside the H and T libraries — verified by reading
the header. The U arm is therefore quantified against the HT reference, and the whole UBL5
question can be answered there.

Why that matters: `gff3/H_T.gff3` and `fa/transcripts_H-T.fa` are already local and
verified, and **the entire downstream chain is built and checked against HT** — the
spliced-cDNA set (`results/sequences/ht_candidates_spliced.json`), the platform caller
(`33`), the 296,613-transcript specificity scan (`29`). Running in UT space would mean
pulling two more large files off `/Volumes/Expansion`, re-deriving all of it against a
second annotation, and reintroducing the annotation-space mixing this project already
decided against for the manuscript.

Use the **HT-space** matrices. Keep `counts/U_T_pipe/` only as a sensitivity check if a
result looks annotation-dependent.

### It is a knockdown, not a knockout

The thesis (§2.5.4) is explicit: shRNA, sequence `CGGGATGAACCTGGAGCTTTAT`, cloned into
`pEN_miRc2` (Shin et al., 2006), subcloned to `pDS_FBneo`, retrovirally delivered, stable
line selected, knockdown confirmed by qPCR (Fig. 2.4, p < 0.01). Several scripts in this
repo — including `04_comparison_T_vs_U.R` and the header of `12_dtu_refit_2x2.R` — call the
U arm a "knockout". **That is wrong and should be corrected as you touch each file.** It
matters for interpretation: residual UBL5 protein means a null result is not evidence of
no requirement.

### Count scale for DTU

`scaledTPM`, not RSEM expected counts. Script 12's header documents why: raw expected counts
are biased by effective-length differences between isoforms *of the same gene*, which is
exactly the comparison DTU makes. Script 12 already writes the expected-count result
alongside as a sensitivity check. Keep that.

---

## 3. What already exists

| script | what it does | state |
|---|---|---|
| `12_dtu_refit_2x2.R` | **Q1 as a single 2×2.** Builds contrasts `WT_LPS`, `KO_LPS`, `interaction` by hand from `~ 0 + group`. The interaction is exactly the Q1 question. | written, needs its input path repointed |
| `04_comparison_T_vs_U.R` | T_UT vs U_UT switch overlap: shared / T-only / U-only, concordance, delta-dIF, per-gene panels | written, UT-space |
| `25`–`35` | fusion screen → structural classification → assay routing → reagent design → specificity scan → selection → figures → bench docs → protocols | **built and verified for HT; reusable with a different candidate list** |

**Script 12 needs one fix before it runs.** It reads
`<counts_dir>/U-T/Isoforms/RSEM.isoforms.expected_count.all_samples.txt`, with `counts_dir`
defaulting to `/Volumes/Expansion/IsoformSwitchAnalyzer/Counts`, and exits 0 with a warning
when that is absent. The local copies are laid out differently —
`counts/U_T_pipe/RSEM.isoforms.*.all_samples.txt`, no `U-T/Isoforms` nesting. Repoint it
via `dtu.counts_dir` in `config/config.yml` plus a path adjustment, or add a local-layout
fallback. Do not copy files around to satisfy a hardcoded path.

**Script 12 does not answer Q2.** It is a DTU model; Q2 is differential abundance. That is
new work — see step 3 below.

**One caution inside script 12.** Its sensitivity design is `~ 0 + group + rep3`, a shared
3-level replicate factor. For the H/T comparison that construction was rejected in this
project because it asserts that donor H1 pairs with flask T1, which is meaningless. For U
vs T it is less obviously wrong — both are THP-1 flasks — but it still asserts that U1 pairs
with T1, which is only true if they were cultured and harvested together. **Check the lab
records before trusting the blocked fit**; treat `~ 0 + group` as primary unless pairing is
confirmed.

---

## 4. The caveat that needs a human answer first

**Which THP-1 line went into the T libraries?**

The knockdown was validated against an **sh-Luciferase control** (Fig. 2.4). But Fig. 7.1
labels the four sequenced arms as human PBMC, mouse BMDM, "THP1 cells", and "THP1 cells with
UBL5 knockdown" — which reads as **plain wildtype THP-1 versus a transduced, G418-selected
knockdown line**.

If T is plain WT rather than sh-Luc, then every WT-vs-KD difference carries retroviral
transduction and antibiotic selection along with UBL5 loss. Consequences:

- **Q2 is affected directly.** A transcript differing in abundance between WT and KD may be
  responding to transduction, not to UBL5. Interferon-stimulated genes are the obvious
  worry — retroviral transduction induces them, and this is an LPS/innate-immunity system
  where ISGs are exactly what you are looking at.
- **Q1 is more robust.** The interaction asks whether the *LPS response* differs between
  genotypes, and the LPS effect is internally controlled within each line. A constant
  transduction offset cancels.

So: **resolve this before reporting Q2**, and if T is plain WT, Q2 ships with the
confound stated plainly and ISGs flagged. Do not bury it in a methods note.
(Note also that Fig. 7.1's legend says "knockdown by siRNA" while §2.5.4 describes shRNA —
an internal inconsistency in the thesis worth resolving at the same time.)

---

## 5. Suggested step plan

1. **Repoint and run script 12** on the HT-space matrix. Report how many genes are
   significant for `interaction`, and for `WT_LPS` and `KO_LPS` separately. This is the
   diagnostic that sizes the whole job — if the interaction returns very few genes, say so
   plainly rather than loosening thresholds to manufacture a list.
2. **Classify the Q1 candidates into the two directions the question actually asks**:
   switching in WT and not in KD (UBL5-dependent), and switching in KD and not in WT
   (UBL5-restrained). A significant interaction alone does not say which; it needs the two
   marginal effects read alongside it. Grade on direction and significance together, the
   way the HT arm's evidence ladder did — never on "called in one and not the other" alone,
   which is a power artefact.
3. **Build the Q2 contrasts** as a new script: transcript-level differential abundance,
   `U_minus vs T_minus` and `U_plus vs T_plus`, limma-voom or edgeR-QLF on the same HT-space
   matrix, reported as two separate tables. Keep it transcript-level, as asked, but also
   report the gene-level result beside it — a transcript that moves only because its whole
   gene moved is a different finding from one that moves against its gene.
4. **Fusion-screen the survivors** with `25_ht_locus_fusion_screen.R`. Non-negotiable: the
   HT arm lost 19 of 63 candidates to fused loci, including the single largest effect.
5. **Structural classification and assay routing** with `26` and `27`.
6. **Sequence extraction, reagent design, specificity scan, selection** with `28`–`30`.
7. **Platform split** with `33` — novel-PacBio vs Illumina-findable, computed from the GFF3
   for *this* candidate set.
8. **Figures, collaborator update, protocols** with `31`, `34`, `35`, plus a standalone
   explorer in `explore/`.

Steps 4–8 are existing, verified code pointed at a new candidate list. The genuinely new
work is steps 1–3.

---

## 6. Gotchas that must not be rediscovered

These cost real time in the HT arm. All are in project memory too, but they are the ones
most likely to bite again here.

1. **`fa/transcripts_H-T.fa` is not uniformly spliced cDNA.** 36,115 of its 296,613 records
   (12.2%) store the unspliced genomic interval. Always check extracted sequence length
   against the summed exon length from the GFF3, and rebuild offenders from hg38 (UCSC
   `getData/sequence`, 0-based half-open, UCSC chrom naming) using GFF3 exon coordinates.
   `results/tables/reference_fasta_audit.csv.gz` already records which records are affected —
   check against it rather than re-auditing.
2. **Specificity is off-target *loci*, not transcripts bound.** A RACE GSP is deliberately
   placed in sequence shared with the partner and a long-read outer pair is meant to amplify
   both isoforms; counting same-gene isoforms as failures rejects designs that are working.
   `results/tables/tx2gene.tsv` maps transcript → locus.
3. **Every qPCR candidate needs two pairs**: isoform-specific (numerator) and total-gene
   with both primers inside one contiguous shared block (denominator). dIF is a fraction,
   and LPS moves total output hard at these loci. This is doubly true here — a genotype
   comparison on top of an LPS comparison makes an un-normalised isoform Ct meaningless.
4. **Check `oligo_role` → strand mapping** in `29_assay_specificity.py`. A role missing from
   `SENSE_ROLES` is silently searched as its own reverse complement and reports zero binding
   sites, which reads as a perfectly clean primer. The script now asserts every role is
   classified; keep that assert if you add a role.
5. **Do not join on `candidate_detectability_by_switch.csv`.** It was built on an earlier
   candidate list, coincidentally also has 44 rows, and shares no `isoform_id` with the
   current shortlist — joining on it yields silently empty platform columns. Use `33` to
   recompute platform for whatever candidate set you end up with, and assert coverage.
6. **Nesting is an assay-routing signal, not a disqualification.** An isoform differing only
   at its 5′ end is contained in its partner by construction; that rules out a unique-region
   qPCR, not the candidate.
7. **ORF comparison by length only.** `orfTransciptStart` is transcript-relative, so two
   isoforms with different 5′ ends have different ORF starts even when the CDS is identical.
8. **`git push` does not work from the analysis sandbox** (no keychain). Commit locally and
   let Sam push.

---

## 7. What "done" looks like

The same bar the HT arm met:

- Every conclusion backed by a CSV plus standalone R/dplyr/ggplot code that reads only CSVs
  and writes nothing.
- Novel-PacBio versus Illumina-findable separated at every stage, computed for this
  candidate set rather than inherited.
- Oligos designed, screened against all 296,613 transcripts, and verified back against the
  spliced cDNA at their recorded coordinates before anything is reported.
- A protocol document whose upstream half is reproduced from the thesis with sections cited
  — for this arm that includes §2.5.4 for the knockdown line and the sh-Luciferase control
  question above.
- Design-vs-measurement contradictions printed where they occur, not averaged away.

## 8. Protocol source for this arm

`/Users/samkatz/Documents/SIGNAL2.0/Thesis_final_sk.pdf` (read-only grant). Methods in
chapter 2 (pp. 95–103), design in chapter 7 (pp. 212–215). `docs/PROTOCOLS.md` already
reproduces THP-1 culture, PMA differentiation, LPS, RNA extraction and RT-qPCR with sections
cited; the UBL5 arm needs §2.5.4 (shRNA line construction) added, and the sh-Luc control
question resolved, but otherwise the upstream protocol is identical — same cells, same
stimulation, same extraction.
