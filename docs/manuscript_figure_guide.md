# Manuscript figures — HT-only build

Built 2026-10-21. Covers experimental design through candidate selection, i.e. everything
**before** assay design. Four figures, 17 panels, 13 real and 4 placeholder.

## Scope, and why it is narrow

This build must stand **without the UBL5 knockout**. Scope is therefore the **HT reference
transcriptome** and the **H and T libraries** only.

`T_UT` is excluded even though it is the same six wildtype libraries as `T_HT` — it is
quantified against the UT annotation, and mixing annotation spaces is the error
`REVIEW_CHANGES.md` §0b documents. `U_UT` is excluded as the knockout arm. The exclusion is
**asserted, not merely intended**: `scripts/21_manuscript_data.R` stops if any UT-space
label reaches a panel table, and `pretty_ds()` in the figure scripts errors on any dataset
value it cannot map.

The check scans **every** character column, not only those named `dataset*`. The first
version looked only at dataset columns and missed `called_in` from `q1_candidates`, whose
values (`"T_UT only"`, `"both annotations"`) describe agreement between the T_HT and T_UT
*annotations* and are therefore UT-dependent despite the column name. It is now dropped
alongside the `*_reannot` columns. If you add a panel, re-run the scan:

```
grep -lE 'T_UT|U_UT' results/tables/manuscript/panel_*.csv   # must return nothing
```

## The framing constraint that governs Figures 2–4

`H_HT` is human **primary** macrophages, three donors, three days, and it is the only dataset
`importRdata` batch-corrected (2 surrogate variables, `limma::removeBatchEffect`). That halves
its replicate IF noise, which is most of why it shows a 2.07% switching rate against `T_HT`'s
0.87%. **Those rates are not two measurements of one quantity** (AGENT_CONTEXT fact 4).

Consequences enforced in the figures:

- Every panel carrying both datasets **facets**; nothing overlays them on a shared axis.
- Panel tables carry a `comparable_across_datasets` column, set `FALSE` wherever it applies.
- In Figure 4, `H_HT` contributes **presence/absence corroboration only**. `dif_cross` is
  deliberately never plotted against `dif_anchor`.

H's value to the argument is that it is an independent cell system, not a second replicate.

## How to run

```
Rscript scripts/21_manuscript_data.R      # writes results/tables/manuscript/panel_*.csv
Rscript scripts/22_manuscript_figures.R   # writes results/figures/manuscript/fig*.png
```

To explore a figure by hand — reads only the CSVs, writes nothing, prints to your device:

```
Rscript explore/manuscript_fig1_explore.R     # ... fig2, fig3, fig4
```

`scripts/22` is a thin wrapper: it sets `options(manuscript.save_dir=)` and sources the same
four explore scripts. There is one copy of the plotting code, so the version a colleague runs
by hand is the version that made the figure.

> Running an explore script from the repo root with no graphics device open writes
> `Rplots.pdf` into the project root, which is a **tracked file**. Either open a device first
> (`pdf(tempfile())`) or restore it afterwards (`git checkout -- Rplots.pdf`). Consider
> gitignoring it.

## Panel inventory

| panel | figure | status | rows | title | source |
|---|---|---|---|---|---|
| f1a | 1A | **placeholder** | 10 | Experimental design schematic | AGENT_CONTEXT fact 4 + lab meeting update sec 6 |
| f1b | 1B | real | 12 | Per-library depth and detection | `counts/H_T_pipe/RSEM.isoforms.expected_count.all_samples.txt` |
| f1c | 1C | real | 26 | Merged HT annotation by class code | `isoformContext_{T_HT_anchor,H_HT_supplementary}.rds` |
| f1d | 1D | real | 2 | Primary objects are pre-reduced | `input_object_reduction_check.rds` |
| f2a | 2A | real | 2 | Switching genes and isoforms | `qc_snapshot_gene_level.rds` + context |
| f2b | 2B | real | 317 | dIF distribution of switching isoforms | context layers, `score_isoforms()` |
| f2c | 2C | real | 20 | Replicate IF noise vs gene abundance | `expression_if_stability.rds` (02d) |
| f2c_floor | 2C | real | 2 | Abundance floor recommendation | `expression_floor_recommendation.rds` (02d) |
| f2d | 2D | real | 14 | Sensitivity to the \|dIF\| threshold | `min_abs_dif_threshold_scan.rds` (00) |
| f3a | 3A | real | 8 | Switching isoforms by detectability tier | 19's tier logic, recomputed HT-only |
| f3b | 3B | real | 8 | Novelty under three definitions | `novelty_detectability_summary.rds` (19) |
| f3c | 3C | real | 34 | Class-code composition, all vs switching | context layers |
| f3d | 3D | real | 4 | Effect size by structural novelty | context layers |
| f4a | 4A | **placeholder** | 10 | Candidate selection funnel | context + config thresholds |
| f4b | 4B, 4C | real | 157 | Anchor vs cross-system corroboration | `q1_candidates.rds` (11), `*_reannot` + `called_in` dropped |
| f4c | 4D | **placeholder** | 6 | Switch anatomy by AS event class | not derivable in HT space |
| f4d | 4E | **placeholder** | 0 | Shortlist with platform class | not derivable in HT space |

## The placeholders, and exactly what each needs

**1A — experimental design schematic.** Five of ten design elements are known. Needs: PMA
differentiation protocol, LPS dose, LPS timepoint, library chemistry (polyA vs
ribo-depleted), and 5′ cap verification. These are the standing wet-lab questions; the last
two are load-bearing, because intron retention is poorly captured by polyA selection and
because most candidate switches turn on where the transcript starts.

**4A, stage 5 — "assayable at the bench".** Stages 1–4 are real (genes quantified → gene
q < 0.05 → above the abundance floor of 12 → |dIF| ≥ 0.15). Stage 5 needs an HT-space
assayability screen: `H_T.gff3` **and** `transcripts_H-T.fa`.

**4D — switch anatomy by AS event class.** Needs `H_T.gff3` (exon coordinates) and a re-run of
`scripts/15_isoform_characterization.R` in HT space.

**4E — candidate shortlist with platform class.** Needs an HT-space re-run of 13
(assayability), 15 (anatomy) and 16 (functional sites), i.e. `H_T.gff3` +
`transcripts_H-T.fa`.

### Why these four cannot be filled from the existing results

Everything describing the **anatomy** of a switch — which AS event it is, how many bases
discriminate the two isoforms, which annotated protein features change, whether the locus is a
two-gene fusion — was computed by scripts 13/15/16 from `U_T.gff3` and `transcripts_U-T.fa`,
i.e. in **UT annotation space**. The candidate transcript ids do not transfer: all 84 resolve
in `T_UT`/`U_UT`, but only **42 of 80** in `T_HT`.

So the existing shortlist — RBKS, ATRN, IRAK3, NCOA7, IL23A, WDR33, and the 16
Illumina-findable vs 18 PacBio-required split from `candidate_detectability_by_switch` — is UT
output and is **not admissible in an HT-only manuscript**. It has to be re-derived, not
relabelled.

All four placeholders resolve with the same two files, both on `/Volumes/Expansion`:

```
/Volumes/Expansion/gtf_files/GTF_to_GFF311_3_2022/H_T/H_T.gff3
transcripts_H-T.fa   (HT equivalent of fa/transcripts_U-T.fa)
```

Note that `/Volumes` is outside the host grant used for this build, so the drive must be
mounted *and* the path granted before the HT-space screens can run.

## What each figure argues

1. **Design and sequencing.** Twelve HT libraries, balanced ±LPS, with T and H structured
   differently on purpose. Panel D states the pre-reduction caveat up front: the primary ISA
   objects are 100% significant by construction, which is why significance is taken from the
   unfiltered context layer (23% / 21% of genes).
2. **The isoform-usage response.** A small, well-powered response — 95 switching genes in T,
   169 in H — with effect sizes thresholded at |dIF| ≥ 0.15 and resting on a gene-level
   abundance floor derived from replicate noise, not chosen by eye.
3. **What PacBio added.** The correction: counting assembly provenance as discovery overstates
   structural novelty by 34.4 pp (T_HT) and 42.7 pp (H_HT). Genuinely novel junction
   structures are 11.2% and 15.6% of switching isoforms. 74.4% / 49.0% of switching isoforms
   would have been findable by short reads against Ensembl.
4. **Candidate selection.** The funnel from ~10,900 quantified genes to 95 switching, and
   cross-system corroboration: of 157 ranked candidates, 8 switch in both T_HT and H_HT, 69
   are T-only despite being testable in H.

The weakest link in the current narrative is Figure 4B/4C: only **8** genes are corroborated
in both systems. That number, not the placeholder panels, is what most constrains how strong a
claim the paper can make before assay validation.
