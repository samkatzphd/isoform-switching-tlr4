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

### Work in UT annotation space. HT is the sensitivity check, not the primary

**An earlier draft of this handover said the opposite. It was wrong, and the reasoning is
worth recording so nobody repeats it.** The argument for HT space was pipeline reuse: the
whole verified assay-design chain is built against HT, and `counts/H_T_pipe/` does carry all
18 libraries including `U1_minus_S9_001` … `U3_plus_S14_001`, so the U arm *is* quantified
there. All of that is true and none of it matters, because it optimises for convenience
instead of for the question.

The references are built from different long-read pools. **HT = H + T long reads. UT = U + T
long reads.** So the UT assembly contains transcript models discovered in the UBL5-knockdown
cells, and the HT assembly cannot contain them — no KD read ever contributed to it. The
sizes show it directly:

| reference | transcripts | genes |
|---|---|---|
| HT | 296,613 | 57,502 |
| UT | **315,349** | 57,841 |

That is **18,736 transcript models present in UT and absent from HT**. Quantifying the U
libraries against HT does not leave those transcripts unmeasured — it is worse than that.
They are not in the index, so reads arising from them are assigned to whichever HT model
they best match, diluted into other isoforms of the locus. A KD-specific isoform does not
show up as missing; it shows up as noise on its neighbours.

This matters most for exactly the question being asked. "Are there novel transcripts
differentially abundant between WT and KD" is **unanswerable in HT space by construction**,
and UBL5 is the reason to expect such transcripts in the first place — it is the human Hub1
homolog, and Hub1 acts on the spliceosome, so a knockdown plausibly produces splice forms
that wildtype cells do not make.

**Use the UT-space matrices as primary.** Because the same six WT THP-1 libraries were
quantified against both references, HT space remains available as a genuine sensitivity
check for any candidate whose model exists in both — and that check is worth running, but it
is a check, not the analysis.

### The bias this introduces, and how to handle it

Switching to UT buys the right question and brings one problem with it. It is **not** that
KD reads preferentially map to KD-discovered models — once the index is fixed, RSEM
quantifies every transcript in all twelve libraries symmetrically, and a WT read matching a
KD-discovered structure is assigned to it normally. The problem is upstream of
quantification: **a structure only enters the index if some sample's long reads supported
it.** A model admitted on the strength of a handful of KD reads enters with weak evidence,
and its apparent KD-specificity is then partly an artefact of which library was in the
discovery pool.

This is the same trap this project already fell into once — `REVIEW_CHANGES.md` §0j, where
`is_novel_pacbio` turned out to encode *which assembler emitted the model*, not whether the
structure was new. Treat discovery provenance as a covariate, not a result:

1. Check what provenance the UT GFF3 carries per transcript. The HT GFF3 has `oId`,
   `class_code`, `cmp_ref`, `num_samples` and `tss_id` on every transcript line; if
   `num_samples` or an equivalent records how many libraries supported the model, that is
   the covariate you need.
2. Require replicate support, not just significance: a KD-specific novel transcript should
   be detected in at least 2 of 3 KD replicates before it is a candidate. The direction-
   agreement machinery from the HT arm (`23`, `24`) already does this kind of grading.
3. Report novel-transcript candidates with their discovery provenance attached, and treat a
   model supported by one library as a hypothesis the assay is there to test — which is
   fine, because producing validation assays is the whole point.
4. Run `33`-equivalent platform classification in UT space. Do not reuse the HT platform
   calls; `class_code` and `cmp_ref` are assembly-specific.

### What this costs, practically

Two files need to come off `/Volumes/Expansion`: the UT GFF3 and the UT transcript FASTA
(by analogy with `gtf_files/GTF_to_GFF311_3_2022/H_T/H_T.gff3` and
`IsoformSwitchAnalyzer/fa/transcripts_H-T.fa` — confirm the actual names). Expect roughly
600 MB and ~2.5 GB respectively, a little larger than the HT pair given the extra models.

**Check free space first.** The boot volume was at ~97% during the HT work. The HT copies
can now be deleted to make room — `gff3/H_T.gff3` and `fa/transcripts_H-T.fa` are gitignored,
every HT result derived from them is committed, and the two things that would be expensive to
regenerate (`results/sequences/ht_candidates_spliced.json` and
`results/tables/reference_fasta_audit.csv.gz`) are in the repo. Both HT files are
recoverable from the drive if an HT sensitivity check is needed later.

Also note: the UT space has its own fused-locus problem, already partly characterised by
`scripts/18` — PLD6+FLCN, DDX58 (=TOPORS+RIGI), IPCEF1+CNKSR3, THAP2+TMEM19, MREG(+PECR),
IFITM2+IFITM1, CCL3+CCL3L3. The fusion screen in step 4 is not optional here either.

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
| `25`–`35` | fusion screen → structural classification → assay routing → reagent design → specificity scan → selection → figures → bench docs → protocols | **built and verified for HT; the logic transfers, the paths do not — see below** |

> **These scripts hardcode the HT inputs.** `gff3/H_T.gff3`, `fa/transcripts_H-T.fa`,
> `results/sequences/ht_candidates_spliced.json`, `results/tables/tx2gene.tsv` and the HT
> assay tables are referenced by literal path throughout `25`–`35`. Running the UBL5 arm in
> UT space means parameterising the reference (a `--ref UT` argument or a config key) rather
> than reusing them as-is. Budget for that: the method is verified, the plumbing is not.
> Specifically, these HT artefacts **cannot** be reused and must be regenerated in UT space:
> `tx2gene.tsv` (transcript→locus map), `reference_fasta_audit.csv.gz` (spliced/unspliced
> audit), the spliced-cDNA sequence set, and all `class_code`/`cmp_ref` platform calls.

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

0. **Get the UT reference onto disk** — GFF3 and transcript FASTA, after checking free space
   (see §2). Verify the FASTA against the GFF3 exon sums before anything is designed against
   it; 12.2% of the HT FASTA records were unspliced genomic intervals and there is no reason
   to assume UT is cleaner. Build the UT `tx2gene` map from the FASTA headers in the same
   pass.
1. **Repoint and run script 12** on the UT-space matrix. Report how many genes are
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
   `U_minus vs T_minus` and `U_plus vs T_plus`, limma-voom or edgeR-QLF on the UT-space
   matrix, reported as two separate tables. Keep it transcript-level, as asked, but also
   report the gene-level result beside it — a transcript that moves only because its whole
   gene moved is a different finding from one that moves against its gene. **This is the
   step that most needs UT space**, and the one where discovery provenance (§2) must be
   carried as a column: partition the results into transcripts annotated in both references
   and transcripts that exist only in UT, and report those two groups separately rather than
   pooled. The second group is the novel-isoform claim the long-read design was for, and it
   is also the group most exposed to discovery bias — both things are true at once.
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
   **The existing `reference_fasta_audit.csv.gz` is HT-specific and does not transfer** — run
   the same audit against the UT FASTA before designing anything (step 0).
2. **Specificity is off-target *loci*, not transcripts bound.** A RACE GSP is deliberately
   placed in sequence shared with the partner and a long-read outer pair is meant to amplify
   both isoforms; counting same-gene isoforms as failures rejects designs that are working.
   `results/tables/tx2gene.tsv` maps transcript → locus for HT; build the UT equivalent from
   the UT FASTA headers.
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
