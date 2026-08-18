# Literature review: alternative splicing in the LPS/TLR4 response, and UBL5 as a candidate splicing regulator

**Scope.** Prepared to support external review of a THP-1 macrophage isoform-switching study
(±LPS, PacBio long-read-augmented reference, IsoformSwitchAnalyzeR, WT vs UBL5 knockout).
Every factual claim carries a citation verified against Crossref or PubMed. Where a claim could
not be sourced, it is marked as a gap rather than asserted. Evidence level is stated explicitly:
**[human macrophage]**, **[human cells]**, **[mouse]**, **[worm/plant]**, **[yeast]**, **[inferred]**.

A companion machine-readable table is in `literature_citations.csv` (64 rows).

---

## 1. UBL5 / Hub1 molecular function

**Bottom line.** UBL5/Hub1 is a well-characterised, non-covalent spliceosome modifier whose loss
impairs splicing of *weak, non-canonical* introns disproportionately to canonical ones — which is
mechanistically the right shape for a substrate-selective, transcriptome-wide splicing deficit. But
its reported condition-dependence rests almost entirely on yeast, and the one mammalian
loss-of-function study reports a *global baseline* intron-retention phenotype that this dataset does
not see.

### 1.1 Mechanism

UBL5 (yeast Hub1) is atypical among ubiquitin-like proteins: it lacks the C-terminal di-glycine
motif required for conjugation and acts by binding partners non-covalently
(Chanarat 2021, *Int J Mol Sci*, 10.3390/ijms22179384). Two distinct surfaces carry two distinct
functions:

- **Asp22 surface → SART1/Snu66 (HIND motif).** Hub1 binds non-covalently to a conserved element
  termed HIND in Snu66 (SART1 in mammals) and Prp38 in plants. Spliceosomes lacking Hub1, or
  defective in the Hub1–HIND interaction, cannot use certain non-canonical 5′ splice sites and fail
  to perform alternative splicing of *SRC1*, while general splicing in *S. cerevisiae* is barely
  affected **[yeast]** (Mishra 2011, *Nature*, 10.1038/nature10143).
- **His63 surface → DDX46/Prp5.** Hub1 binds the DEAD-box helicase Prp5, a regulator of early
  spliceosome assembly, and stimulates its ATPase activity. High Hub1 enhances splicing efficiency
  but *relaxes fidelity*, tolerating suboptimal splice sites and branchpoint sequences; a
  Hub1-dependent cryptic intron in *PRP5* itself provides negative feedback that curbs excessive
  mis-splicing **[yeast]** (Karaduman 2017, *Mol Cell*, 10.1016/j.molcel.2017.06.021).

The practical consequence of the Karaduman model is that **UBL5 level sets how permissive the
spliceosome is toward suboptimal splice sites**. This is a dose-dependent knob, not a binary switch
— which is the correct shape for a KO that retains a substantial fraction of a splicing response
rather than abolishing it.

### 1.2 Conservation and essentiality — a caveat for this dataset

| System | Essential? | Phenotype | Source |
|---|---|---|---|
| *S. cerevisiae* | No | Mating defects; *SRC1* AS fails; general splicing largely intact | Mishra 2011 (10.1038/nature10143) |
| *S. pombe* | **Yes** | Cell-cycle defects, inefficient pre-mRNA splicing | Wilkinson 2004 (10.1016/j.cub.2004.11.058); Yashiroda 2004 (10.1111/j.1365-2443.2004.00807.x) |
| *C. elegans* | Yes (dev.) | Arrest and death at larval stage 3, splicing defects | Kolathur 2022 (10.1002/1873-3468.14555) |
| Human cells | **Reported essential** | Speckle abnormalities, partial nuclear mRNA retention, mitotic catastrophe, apoptosis on prolonged depletion | Ammon 2014 (10.1093/jmcb/mju026) |

> **This is a direct challenge to the experimental system.** Ammon et al. report human Hub1/UBL5 is
> *essential for viability* **[human cells]**. A viable UBL5 knockout THP-1 line with an
> indistinguishable baseline isoform composition is therefore not the expected result. The review
> should ask the authors to document (a) the exact lesion and whether it is a true null or a
> hypomorph, (b) residual UBL5 protein by western blot, and (c) whether any compensating paralogue
> or partial-length product exists. A hypomorphic allele would in fact fit the data better than a
> null — and would fit the dose-dependent Prp5 model above.

### 1.3 Knockdown phenotypes in human cells — a partial contradiction

UBL5 depletion in human cells **decreases pre-mRNA splicing efficiency and causes globally enhanced
intron retention**; the specific downregulation of the cohesion factor Sororin through intron
retention explains a sister-chromatid cohesion defect **[human cells]** (Oka 2014, *EMBO Rep*,
10.15252/embr.201439478). Notably, a UBL5 D22A mutant — the mutation that abrogates Hub1–Snu66
binding in yeast — still bound SART1 and still rescued the cohesion defect, indicating the SART1
interface is not the only functional surface in human cells.

> **This CONTRADICTS the dataset's finding of an indistinguishable baseline.** The published
> mammalian loss-of-function phenotype is a *constitutive, global* intron-retention defect. Two
> reconciliations are available and both are testable: (i) the KO is hypomorphic and retains enough
> UBL5 for baseline but not for surge demand; (ii) siRNA depletion in Oka et al. was acute and deep
> whereas a stable KO line has had time to adapt. This tension should be stated in the manuscript,
> not elided.

### 1.4 Is the role constitutive or condition-dependent?

This is the crux of Q2, and the honest answer is: **in mammals, essentially unstudied. The
condition-dependent evidence is yeast, and the stress is not TLR4.**

- **[yeast]** Hub1 is transcriptionally upregulated via the yeast AP-1 (Yap1) regulon upon oxidative
  and heavy-metal stress and promotes efficient splicing of introns with non-canonical splice sites;
  *hub1* mutants become cadmium-sensitive when metallothionein defence is impaired, and other
  splicing-factor mutants showed similar cadmium sensitivity (Chanarat 2020, *BBA Mol Cell Res*,
  10.1016/j.bbamcr.2019.118565). This is the single closest published precedent for
  "UBL5 becomes important specifically under stimulus."
- **[human cells / mouse]** UBL5 is an ER-stress-responsive protein, undergoing rapid depletion in
  mammalian cells and mouse liver via proteasome-dependent but ubiquitin-independent proteolysis;
  PERK activation is necessary and sufficient for this degradation (Wang 2023, *JBC*,
  10.1016/j.jbc.2023.104915). This is the strongest *mammalian* evidence that UBL5 abundance is
  itself under stimulus control — though the direction (stress *depletes* UBL5) is opposite to the
  yeast oxidative-stress result.
- A review explicitly frames UBL5/Hub1 as a **stress-responsive** regulator across splicing, the
  Fanconi anaemia pathway and UPRmt (Chanarat 2021, 10.3390/ijms22179384) — but it is a review, not
  primary evidence, and the splicing-under-stress primary data it rests on are the yeast work above.

> **Gap, stated plainly:** there is **no published study of UBL5 in an LPS/TLR4 response, in
> macrophages, or in any mammalian innate-immune stimulation.** A condition-dependent mammalian
> splicing requirement for UBL5 has not been demonstrated by anyone. If this dataset shows one, it
> is a genuinely novel claim — which raises rather than lowers the evidentiary bar.

### 1.5 The mitochondrial-UPR role — real, but does not transfer cleanly

The UPRmt link is genuine and primary, not folklore:

- *ubl-5* RNAi compromises UPRmt signalling in *C. elegans*, perturbs mitochondrial morphology and
  multi-subunit complex assembly, and mitochondrial stress drives nuclear accumulation of GFP-tagged
  UBL-5 **[worm]** (Benedetti 2006, *Genetics*, 10.1534/genetics.106.061580).
- Mitochondrial unfolded-protein stress drives complex formation between the homeodomain
  transcription factor DVE-1 and UBL-5, which binds mitochondrial chaperone gene promoters
  **[worm]** (Haynes 2007, *Dev Cell*, 10.1016/j.devcel.2007.07.016).

**Caveat the user should apply:** this axis is defined in *C. elegans*, DVE-1 has no clean human
orthologue, and the mammalian UPRmt is routed through ATF4/ATF5/CHOP instead. Do **not** present
"UBL5 has a mitochondrial-UPR role" as a mammalian fact. It is a worm result with an unresolved
mammalian counterpart. It is, however, a legitimate reason to check whether the KO has a
mitochondrial-stress signature that could confound the splicing interpretation.

### 1.6 Splicing-independent UBL5 functions (confounders)

UBL5 promotes functional integrity of the Fanconi anaemia DNA-repair pathway via FANCI
**[human cells]** (Oka 2015, *EMBO J*, 10.15252/embj.201490376). Any UBL5 KO phenotype attributed to
splicing must contend with the fact that UBL5 also has DNA-repair and (in worms) transcriptional
co-regulator activities.

### 1.7 Immunity and inflammation

The only retrieved UBL5-immunity literature is a review of UBL5 in viral infection, whose one
concrete antiviral result is UBL5 binding and mediating degradation of rice stripe virus NS3
(Xia 2024, *Viruses*, 10.3390/v16121922) — a plant-virus system. **No mammalian innate-immune,
inflammatory, or macrophage function for UBL5 has been established.** Searches for UBL5 with
macrophage, LPS, TLR, interferon and inflammation terms returned nothing on point.

---

## 2. Alternative splicing in TLR4 / LPS / innate immune signalling

**Bottom line.** Splicing is an established regulatory layer in TLR signalling — specifically a
negative-feedback layer producing dominant-negative isoforms — and multiple splicing factors
(including core spliceosome components) demonstrably tune innate immune output. Crucially, at least
two studies show splicing regulators whose function is **dispensable at baseline and required only
upon stimulation**, which is the same architecture this dataset proposes for UBL5.

### 2.1 Splicing as regulation, not noise

Alternative splicing of at least 13 genes in the TLR signalling pathway can produce dominant-negative
inhibitors of signalling; many are induced by immune challenge, making their production a negative
feedback loop (Lee & Alper 2022, *Front Immunol*, 10.3389/fimmu.2022.1023567). This review is the
best single entry point for section 2.

**Documented nodes:**

| Node | Isoform | Function | Evidence level | Source |
|---|---|---|---|---|
| MyD88 | MyD88-S (exon 2 skip) | Dominant-negative; fails to recruit IRAK4, so IRAK1 is not phosphorylated | [human/mouse cells] | Janssens 2002 (10.1016/S0960-9822(02)00712-1); Burns 2003 (10.1084/jem.20021790) |
| IRAK1 | IRAK1c | Negative regulator of Toll/IL-1R signalling | [cells] | Rao 2005 (10.1128/MCB.25.15.6521-6532.2005) |
| TLR4 | smTLR4 (soluble) | Inhibits LPS-induced NF-κB and TNF; LPS-induced in RAW264.7 | **[mouse only]** | Iwami 2000 (10.4049/jimmunol.165.12.6682) |
| MD-2 (LY96) | MD-2s | Negatively regulates LPS-induced TLR4 signalling | **[human]** | Gray 2010 (10.4049/jimmunol.0903543) |
| NLRP3 | Δexon 5 | Lacks NEK7 interaction surface → inactive; stochastically regulated | **[human MDM and THP-1]** | Hoss 2019 (10.1038/s41467-019-11076-1) |

> **Species caution the review should flag.** The soluble-TLR4 story is **mouse**. There is no direct
> homologue of the mouse smTLR4 alternative exon in the human *TLR4* gene, and although human *TLR4*
> transcript variants predicted to produce a similar truncated protein exist, this has not been
> tested experimentally (Lee & Alper 2022, 10.3389/fimmu.2022.1023567). If the dataset reports a
> human TLR4 isoform switch, it cannot be interpreted through the mouse smTLR4 literature.
>
> The NLRP3 result is the most directly transferable: it is human, includes **THP-1 cells**, and
> demonstrates a functionally consequential isoform difference in an innate immune sensor.

### 2.2 Splicing factors implicated in macrophage LPS responses

- **Core spliceosome (SF3A/SF3B).** When SF3A levels or activity are diminished, MyD88-S levels rise
  through skipping of MyD88 exon 2, limiting the innate immune response **[macrophage]**
  (De Arras & Alper 2013, *PLoS Genet*, 10.1371/journal.pgen.1003855). This is the closest published
  analogue to perturbing UBL5: a *core spliceosome* component, not an SR/hnRNP regulator, tuning
  immune output.
- **SRSF1 / hnRNP U.** SRSF1 and HNRNPU binding to MyD88 pre-mRNA is regulated by LPS, LPS exposure
  decreases SRSF1 expression, and both regulate the extent of the innate immune response
  **[macrophage]** (Lee 2024, *J Mol Biol*, 10.1016/j.jmb.2024.168497).
- **hnRNP M — the key precedent.** hnRNP M represses innate immune gene expression, and its function
  is regulated by pathogen-sensing cascades. Mutating specific serines on hnRNP M had **little effect
  on pre-mRNA splicing or housekeeping transcript levels in resting macrophages, but greatly
  impacted its ability to dampen induction of specific innate immune transcripts following pathogen
  sensing** **[macrophage]** (West 2019, *Cell Rep*, 10.1016/j.celrep.2019.09.078).
- **SR/hnRNP panel.** Knockdown of ten SR/hnRNP proteins in RAW264.7 macrophages showed these factors
  influence *different genes* in uninfected versus *Salmonella*-infected cells and drive differential
  isoform usage for thousands of transcripts, differently in each state **[mouse macrophage]**
  (Wagner 2021, *Front Immunol*, 10.3389/fimmu.2021.656885).
- **SRSF3.** Upregulated in LPS-stimulated RAW264.7; knockdown suppresses inflammatory cytokines and
  raises the short MD2B variant of the TLR4 co-receptor **[mouse macrophage]** (Fu 2024,
  10.3390/cimb46060372).

> **West 2019 is the single most important comparator for this manuscript.** It establishes, in
> macrophages, exactly the architecture claimed here — a splicing regulator dispensable at baseline
> and required on stimulation. The *direction* differs (hnRNP M loss causes hyper-induction; UBL5
> loss here causes blunted response), so it is a precedent for the *architecture*, not for the
> *phenotype*. That distinction should be made explicitly rather than citing it as support for the
> result itself.

### 2.3 Timing, kinetics and the transcription-versus-splicing balance

- Kinetic analysis of chromatin, nucleoplasmic and cytoplasmic RNA fractions in LPS-stimulated
  macrophages showed the inflammatory response is governed **primarily at transcription initiation**,
  with intron-containing full-length nascent transcripts accumulating on chromatin before release
  **[mouse macrophage]** (Bhatt 2012, *Cell*, 10.1016/j.cell.2012.05.043).
- However, for a subset of Lipid A-response genes, ligation of certain exon pairs is *delayed*
  relative to synthesis of the complete transcript, so **splicing kinetics and chromatin release
  limit the rate of induced gene expression** **[mouse macrophage]** (Pandya-Jones 2013, *RNA*,
  10.1261/rna.039081.113).
- Macrophage development and activation involve **coordinated intron retention in key inflammatory
  regulators** (Green 2020, *NAR*, 10.1093/nar/gkaa435). Detained introns are a widespread class of
  post-transcriptionally spliced introns (Boutz 2015, *Genes Dev*, 10.1101/gad.247361.114) — a
  plausible mechanism by which a splicing-efficiency deficit blunts the *speed and amplitude* of
  induction without changing steady state.

> **A methodological warning that matters a great deal here.** About **50% of splicing changes after
> LPS activation in both human and murine macrophages are alternative-first-exon (AFE) events**
> (Robinson 2021, *eLife*, 10.7554/eLife.69431). AFE events are promoter-choice events, not
> spliceosome-catalysed exon-selection events. IsoformSwitchAnalyzeR will score them as isoform
> switches. **If half of this dataset's "splicing response" is AFE, then attributing its loss to a
> spliceosome modifier is mechanistically incoherent** — UBL5 acts at 5′ splice-site selection and
> Prp5-dependent A-complex formation, not at promoter choice. This is the single highest-value
> analysis the authors could add: stratify the LPS isoform-switching response by
> `AlternativeSplicingAnalysis` event class (ATSS/AFE vs ES/IR/A5/A3), and re-compute the WT-vs-KO
> retention figure *within* the spliceosome-plausible classes. If the 38% retention figure holds
> specifically in IR and A5/A3 events — the classes Hub1 mechanistically governs — the argument
> becomes strong. If the deficit is concentrated in ATSS/AFE, it points somewhere other than UBL5's
> known mechanism.

### 2.4 Transcriptome-wide LPS splicing surveys

A deep RNA-seq study of LPS-stimulated human PBMCs from three donors reported 490 differentially
expressed genes and differential alternative splicing affecting TLR signalling, PI3K/AKT signalling
and pro-inflammatory macrophage polarisation (Chavez-Iglesias 2025, 10.1101/2025.11.09.687437).
**This is a preprint and has not been peer reviewed** — use only as a directional comparison, and
note it is PBMCs (mixed populations) rather than a macrophage line.

---

## 3. Specific candidate genes

**Bottom line.** Of the four priority genes, only **RAB7B** and **IRAK3** have strong LPS/macrophage
literature, and *neither* has a functionally characterised alternative isoform. **SOCS4** and
**SPRING1** have essentially no isoform-level or macrophage-specific literature. For every gene on
this list, an isoform switch would be a novel observation with no prior functional anchor — which
means the manuscript cannot lean on "known biology" to validate these calls.

### Priority genes

**IRAK3 (IRAK-M)** — *strong gene-level prior, no isoform-level support.*
IRAK-M is an LPS-inducible pseudokinase and negative regulator of TLR signalling that prevents
dissociation of IRAK1/IRAK4 from MyD88 and formation of the IRAK–TRAF6 complex **[mouse/human
macrophage]** (Kobayashi 2002, *Cell*, 10.1016/S0092-8674(02)00827-9); in humans its expression is
restricted to the myeloid lineage (Pereira & Gazzinelli 2023, *Front Immunol*,
10.3389/fimmu.2023.1133354). RefSeq annotation records that alternative splicing produces multiple
transcript variants, including one with a deletion spanning the death-domain region
(NCBI Gene ID 11213). **However: no functional characterisation of a distinct or dominant-negative
*human* IRAK3 splice isoform was retrieved.** The only systematic functional test of IRAK3 splice
variants found is in **rainbow trout**, and it was largely negative — overexpressed full-length and
truncated variants had only modest, non-dose-dependent effects on TLR-stimulated NF-κB, though one
C-terminal-domain variant did quench IL-1β and IL-8 production (Rebl 2019, *Front Immunol*,
10.3389/fimmu.2019.02246). **Do not assume a dominant-negative IRAK3 isoform; the one direct test
argues against it.**

**RAB7B** — *strong macrophage prior, zero isoform literature.*
Rab7b negatively regulates TLR4 signalling in macrophages by promoting lysosomal degradation of
TLR4, suppressing LPS-induced TNF-α, IL-6, NO and IFN-β **[mouse macrophage]** (Wang 2007, *Blood*,
10.1182/blood-2007-01-066027). It also suppresses TLR9-initiated cytokine and type I IFN production
by promoting TLR9 degradation, and TLR9 ligation inhibits Rab7b expression via ERK and p38
(Yao 2009, *J Immunol*, 10.4049/jimmunol.0900249). **No documented functionally distinct RAB7B splice
isoform.** A RAB7B switch here would be new.

**SOCS4** — *thin literature, phenotype is not macrophage-intrinsic.*
SOCS4-deficient mice succumb rapidly to pathogenic H1N1 influenza with dysregulated pulmonary
cytokine and chemokine production, delayed viral clearance and impaired trafficking of
influenza-specific CD8 T cells linked to defects in T-cell receptor activation — described as the
first functional phenotype reported for SOCS4 **[mouse]** (Kedzierski 2014, *PLoS Pathog*,
10.1371/journal.ppat.1004134). Note the phenotype is largely **T-cell**, not macrophage-intrinsic.
No isoform-level literature. Interpret a SOCS4 switch cautiously.

**SPRING1 / C12orf49** — *well-defined function, entirely outside immunology.*
Identified by haploid genetic screens as a Golgi-resident determinant of SREBP signalling; its loss
reduces SCAP levels and mislocalises SCAP, attenuating SREBP-dependent programmes, and *Spring*
deletion is embryonic lethal in mice (Loregger 2020, *Nat Commun*, 10.1038/s41467-020-14811-1).
Independently identified by metabolic coessentiality mapping (Bayraktar 2020, *Nat Metab*,
10.1038/s42255-020-0206-9), and shown to act by promoting maturation of site-1 protease
(Xiao 2020, *Protein Cell*, 10.1007/s13238-020-00753-3). Relevant *only* by the indirect route that
LPS reprograms macrophage lipid metabolism. **No immune function and no isoform-level data.**

### Remaining genes

**PLD6** — MitoPLD is anchored in the mitochondrial outer membrane, hydrolyses cardiolipin to
phosphatidic acid to drive mitochondrial fusion, and separately has nuclease activity in piRNA
biogenesis (Frohman 2015, *Trends Pharmacol Sci*, 10.1016/j.tips.2015.01.001). No innate-immune or
isoform-level role. Weak candidate — though note the coincidence that both UBL5 (via UPRmt) and PLD6
touch mitochondrial biology, which is worth a sanity check rather than a claim.

**SP1, LIG1, IPCEF1, ATRN, SAMD14, CTTNBP2NL** — searches returned **no characterisation of
functionally distinct alternative isoforms** in an immune or macrophage context for any of these.
SP1 and LIG1 have annotated variants; no macrophage isoform-function literature was found. These
should be reported as uncharacterised hits, not interpreted.

---

## 4. Methodological literature

**Bottom line.** The project's abundance floor on *gene* expression is well justified by the
statistics of ratio estimators; the field's own guidance is that n=3 is at the low end for DTU FDR
control, and that DRIMSeq specifically is the method most prone to FDR inflation. Long-read
augmentation is supported and does change which switches are detectable — which is a reason the
results cannot be compared naively to short-read-reference studies.

### 4.1 dIF/PSI estimation at low expression

- PSI is a ratio of reads supporting inclusion to all reads spanning the event, so at low read
  support PSI estimates become **highly variable**; unlike expression, splicing ratios cannot be
  treated as zero when missing (Jiang 2026, *iScience*, 10.1016/j.isci.2026.116090). Framed for
  single cells, but the estimator mathematics is identical for low-abundance bulk genes.
- Robust PSI estimation requires an explicit measure of estimate uncertainty, because
  position-specific bias can dramatically influence estimates, particularly for transcripts with
  minimal coverage (Kakaradov 2012, *BMC Bioinformatics*, 10.1186/1471-2105-13-S6-S11).
- Beta distributions modelled from inclusion/exclusion counts yield PSI estimates whose width narrows
  as coverage rises, making precision explicit rather than relying on **context-independent
  predefined cutoffs** (Ascensão-Ferreira 2024, *RNA*, 10.1261/rna.079764.123).

> **Supports this project's design decision.** The `min_gene_expression: 12` floor, derived
> empirically from within-condition replicate IF noise, is the same idea betAS implements
> analytically: a fixed dIF threshold is not equally meaningful at all coverages. The published
> literature does *not* supply a canonical numeric floor, so the empirical derivation in script `02d`
> is the right approach and should be presented as a contribution, not a nuisance parameter. One
> refinement worth noting: the field's framing (Ascensão-Ferreira 2024) favours *coverage-aware
> uncertainty per event* over a *hard global threshold* — a per-gene weighting is the more principled
> version of the same instinct.

### 4.2 Differential transcript usage with n=3

| Finding | Source |
|---|---|
| Standard workflow filters lowly expressed genes/transcripts before DTU testing (too little power); stage-wise testing via stageR recommended for gene-then-transcript inference | Love 2018, *F1000Res*, 10.12688/f1000research.15398.1 |
| DRIMSeq: Dirichlet-multinomial model with likelihood-ratio test on relative transcript abundance | Nowicka 2016, *F1000Res*, 10.12688/f1000research.8900.1 |
| **DEXSeq and DoubleExpSeq had the highest performance; satuRn was on par with them; DRIMSeq's performance varied strongly between datasets and most methods — DRIMSeq in particular — failed to control FDR at nominal levels** | Gilis 2022, *F1000Res*, 10.12688/f1000research.51749.2 |
| SUPPA2 and RATs always controlled FDR but with consistently low sensitivity (~50%); DRIMSeq and DEXSeq had higher sensitivity while sometimes exceeding target FDR; FDR control improved with larger sample size | Erdogdu 2024, *Cell Rep Methods*, 10.1016/j.crmeth.2024.100736 |
| **Performance of all tools improves with eight replicates**; edgeR and DRIMSeq achieved low FDR across DTU event types, while satuRn and DEXSeq maintained low FDR with higher recall in some scenarios | Lio 2025, *NAR Genom Bioinform*, 10.1093/nargab/lqaf117 |
| Propagating quantification uncertainty from bootstrap replicates (CompDTU/CompDTUme) improves sensitivity and specificity vs RATs, DRIMSeq, DEXSeq and SUPPA2, and reduces FPR vs DRIMSeq | Young 2023, *Biostatistics*, 10.1093/biostatistics/kxad008 |
| SUPPA2 provides fast, uncertainty-aware differential splicing across multiple conditions | Trincado 2018, *Genome Biol*, 10.1186/s13059-018-1417-1 |

> **Actionable for this project.** (1) If the ISA runs used `isoformSwitchTestDEXSeq`, that is the
> better-supported choice; if any analysis used DRIMSeq, the FDR is likely optimistic (Gilis 2022).
> (2) At n=3 the design sits *below* every benchmark's comfortable regime — Lio 2025 shows all tools
> improve at eight replicates, and Erdogdu 2024 shows FDR control specifically improves with sample
> size. This argues for leading with effect sizes over p-values (which the project already does) and
> for treating the KO-vs-WT retention percentages as estimates with wide intervals rather than point
> facts. (3) Salmon/kallisto bootstrap replicates are already available in most such pipelines;
> propagating them (Young 2023) would materially strengthen the low-abundance calls that the
> expression floor is currently handling bluntly.
>
> **On paired designs:** no benchmark retrieved gives specific guidance for paired/blocked DTU
> designs. DEXSeq and satuRn both accept general design formulae so blocking is implementable, but
> the FDR behaviour of paired DTU designs at n=3 appears **unbenchmarked** in the published
> literature. That is a genuine gap, not an oversight in the search.

### 4.3 Long-read-augmented references

- LRGASP consortium: libraries with longer, more accurate reads produce more accurate transcript
  models, while greater read depth improves quantification accuracy; reference-based tools performed
  best in well-annotated genomes; **detecting novel transcripts was harder than recovering annotated
  ones**, and quantifying complex and lowly expressed transcripts remained challenging
  (Pardo-Palacios 2024, *Nat Methods*, 10.1038/s41592-024-02298-3).
- Short and long reads gave comparable *gene-level* estimates but **differed substantially for
  individual isoforms**, and several novel long-read-discovered transcripts participated in
  cell-line-specific isoform switching events (Chen 2025, *Nat Methods*, 10.1038/s41592-025-02623-4).

> **Direct answer to the question "do long-read-augmented references change isoform-switch calls?"
> — yes, demonstrably** (Chen 2025). The PacBio-for-annotation / short-read-for-quantification design
> is the configuration LRGASP recommends. Two caveats follow: novel PacBio-derived models are the
> *least* reliable class, particularly at low abundance, so novel-isoform-driven switches deserve
> separate reporting from annotated-isoform switches; and results are not directly comparable to
> GENCODE-only studies.

### 4.4 IsoformSwitchAnalyzeR-specific

The tool identifies isoform switches and predicts functional consequences genome-wide
(Vitting-Seerup & Sandelin 2019, *Bioinformatics*, 10.1093/bioinformatics/btz247). In the TCGA
landscape paper, isoform switches with predicted functional consequences affected approximately
**19% (N=2,352) of multi-isoform genes** across >5,500 patients (Vitting-Seerup & Sandelin 2017,
*Mol Cancer Res*, 10.1158/1541-7786.MCR-16-0459).

> **Warning — a numerical coincidence that will mislead readers.** That 19% is the fraction of
> *multi-isoform genes showing consequential switches*. This dataset's ~19.5% is the fraction of
> *switching genes with little gene-level expression change*. **These are different denominators and
> different quantities.** They should never be juxtaposed in the manuscript without explicit
> distinction; a reader will otherwise read a replication where none exists.

---

## 5. Spliceosome perturbation and blunted stimulus responses

**Bottom line.** The general phenomenon — perturb a core splicing component, get a *reduced-amplitude*
response to a stimulus while cells remain viable — is published, and the strongest example is
quantitatively close to this dataset's effect size. But the same literature shows spliceosome
perturbation can *activate* interferon signalling, which provides a sharp falsifiable control.

- **The closest precedent.** siRNA knockdown of SF3B1 in HeLa cells reduced heat-shock induction of
  *HSPA6* to **18%** and *DNAJB1* to **31%** of control — nearly as severe as knocking down HSF1
  itself — and SF3B1 regulates both HSF1 concentration and activity **[human cells]**
  (Kim Guisbert 2017, *PLOS ONE*, 10.1371/journal.pone.0176382). A core spliceosome perturbation
  blunting an inducible stress response to roughly a third of normal is directly comparable in
  magnitude to a KO retaining ~38% of the splicing response.
- **Graded, not binary.** Titrating SF3B1 activity with different pladienolide B concentrations
  produces distinct, dose-dependent effects on the transcriptome and cell physiology rather than a
  single all-or-nothing splicing failure (Kim Guisbert 2020, *Int J Mol Sci*, 10.3390/ijms21249641).
  Supports a partial-loss model.
- **Splicing is rate-limiting for induced genes.** Delayed exon ligation and chromatin release limit
  the rate of Lipid A-induced gene expression (Pandya-Jones 2013, 10.1261/rna.039081.113), and
  detained introns provide a regulated post-transcriptional reservoir (Boutz 2015,
  10.1101/gad.247361.114). Together these give a mechanism for amplitude loss without baseline change.
- **The counter-evidence, and a control the authors should run.** Pharmacological modulation of SF3B1
  with pladienolide B induces aberrant RNA species and a **robust type I interferon response** via
  RIG-I and IRF3 (Chang 2021, *JBC*, 10.1016/j.jbc.2021.101277).

> **This last point is the most useful experimental suggestion in this review.** If gross spliceosome
> dysfunction were occurring in the UBL5 KO, the expected signature is **elevated baseline type I
> interferon signalling** from aberrant RNA species. The dataset reports the opposite pattern —
> interferon-response genes *retain more* splicing response than average, and baseline composition is
> normal. **Explicitly checking that baseline ISG expression is not elevated in the KO would convert
> this from an absence of evidence into a positive control** supporting a specific, non-catastrophic
> UBL5 function rather than general spliceosome breakage.

---

## 6. Direct relevance to this dataset's findings

### Finding 1 — KO retains ~64% of the expression response but only ~38% of the splicing response

**Supported in principle.** The magnitude is comparable to the SF3B1/heat-shock precedent, where core
spliceosome knockdown reduced induction of stress genes to 18–31% of control while cells stayed
viable (Kim Guisbert 2017, 10.1371/journal.pone.0176382). The mechanistic route — splicing kinetics
limiting the rate of induced gene expression (Pandya-Jones 2013, 10.1261/rna.039081.113) — is
established in macrophages specifically.

**Not supported directly.** No published work measures expression-response and splicing-response
retention as separable quantities in a splicing-factor mutant, so the 64%/38% *ratio* has no
literature comparator. It is a novel measurement.

**Caution.** The ratio is sensitive to how each response is normalised and to the abundance floor
applied to only one of the two. A reviewer will ask whether the same filtering universe underlies
both percentages.

### Finding 2 — Baseline isoform composition indistinguishable from WT

**Partially CONTRADICTED by the mammalian literature.** UBL5 depletion in human cells causes
**globally enhanced intron retention** at baseline (Oka 2014, 10.15252/embr.201439478), and human
UBL5 is reported essential (Ammon 2014, 10.1093/jmcb/mju026). A viable KO with a normal baseline is
not the published expectation.

**Supported by yeast.** Hub1 binding "barely affects general splicing" in *S. cerevisiae* while
selectively enabling non-canonical 5′ splice-site usage (Mishra 2011, 10.1038/nature10143) — exactly
a normal-baseline, selective-deficit phenotype. But that is yeast, and yeast Hub1 is non-essential.

**Architecturally supported in macrophages.** hnRNP M serine mutations had little effect on splicing
in *resting* macrophages but greatly impaired dampening of induction after pathogen sensing
(West 2019, 10.1016/j.celrep.2019.09.078). This is the strongest support that "no baseline defect,
stimulus-specific deficit" is a real architecture for a splicing regulator — while noting the
direction of the hnRNP M effect is opposite.

**Recommended:** report residual UBL5 protein and allele structure; a hypomorph reconciles the
contradiction cleanly and fits the dose-dependent Prp5 model (Karaduman 2017,
10.1016/j.molcel.2017.06.021).

### Finding 3 — Deficit is transcriptome-wide and *spares* LPS-responsive and interferon-response genes

**Consistent with, and partly explained by, the transcription-dominance of the LPS response.** The
inflammatory response is governed primarily at transcription initiation (Bhatt 2012,
10.1016/j.cell.2012.05.043), so the most strongly induced genes are the ones least dependent on
splicing regulation for their output — they would be expected to be *less* sensitive to a splicing
deficit, not more.

**A serious alternative explanation the authors must exclude.** LPS-responsive and interferon genes
are also the most highly *expressed* genes after stimulation. Since the analysis applies a gene-level
abundance floor and dIF is a ratio estimator whose variance collapses with coverage
(Kakaradov 2012, 10.1186/1471-2105-13-S6-S11; Ascensão-Ferreira 2024, 10.1261/rna.079764.123), high-expression
genes will retain *measurable* splicing response more reliably than low-expression genes **regardless
of biology**. The "sparing" of LPS/IFN genes may be an expression-abundance artefact.
**Test: regress splicing-response retention on gene expression level, and ask whether the LPS/IFN
sparing survives conditioning on abundance.** Given this project's stated posture that
plausible-looking results have repeatedly been artefacts, this is the finding most in need of that check.

**Additionally:** confirm baseline ISG levels are *not* elevated in the KO (Chang 2021,
10.1016/j.jbc.2021.101277), which would otherwise indicate spliceosome dysfunction and confound the
interferon-sparing interpretation.

### Finding 4 — ~19.5% "cryptic switchers" (isoform switch without gene-level expression change)

**Conceptually well supported.** This is the founding rationale for isoform-level analysis: shifts in
isoform usage cannot be detected at gene level (Vitting-Seerup & Sandelin 2019,
10.1093/bioinformatics/btz247). The macrophage literature supplies concrete instances — the human
NLRP3 Δexon 5 isoform differs functionally without requiring a gene-level change
(Hoss 2019, 10.1038/s41467-019-11076-1), and intron retention is coordinated during macrophage
activation independently of expression (Green 2020, 10.1093/nar/gkaa435).

**Explicit warning against a false comparison.** The ~19% figure in Vitting-Seerup & Sandelin 2017
(10.1158/1541-7786.MCR-16-0459) is a *different quantity* (multi-isoform genes with consequential
switches, in cancer). The numerical agreement is coincidental.

**Caution.** A gene with no net expression change is, on average, at lower total abundance than a
strongly induced one — so the cryptic-switcher class is enriched for exactly the regime where dIF is
least stable. Report the abundance distribution of cryptic switchers against all switchers.

### Finding 5 — ~2% where gene expression falls, dominant isoform is preserved, minor isoforms collapse

**Weakest literature support of the five.** No published class matching this description was
retrieved.

**A mundane explanation must be excluded first.** When total gene expression falls, minor isoforms —
being low-count by definition — lose read support fastest, so their apparent "collapse" and the
dominant isoform's apparent "preservation" is the expected behaviour of a ratio estimator under
declining coverage (Kakaradov 2012, 10.1186/1471-2105-13-S6-S11; Jiang 2026,
10.1016/j.isci.2026.116090). This pattern can arise from sampling alone.

**If it survives that check,** the nearest mechanistic frames are NMD-coupled degradation of minor
isoforms and the detained-intron reservoir (Boutz 2015, 10.1101/gad.247361.114) — both testable
against the ISA NMD-sensitivity predictions already available in the pipeline. Given the project's
own record of small-n patterns failing Fisher tests around p ≈ 0.24, a 2% class deserves an explicit
power statement before it is described as a class at all.

---

## 7. Summary of gaps and contradictions

**Genuine gaps (nothing published):**

1. UBL5 in macrophages, LPS, TLR signalling, or any mammalian innate-immune stimulation — **nothing**.
2. A condition-dependent mammalian splicing requirement for UBL5 — **not demonstrated by anyone**.
3. Functionally characterised alternative isoforms for RAB7B, SOCS4, SPRING1, IPCEF1, ATRN, PLD6,
   SAMD14, CTTNBP2NL — **none found**; for human IRAK3, annotated but functionally uncharacterised.
4. FDR behaviour of **paired** DTU designs at n=3 — **unbenchmarked**.
5. Separable expression-response vs splicing-response retention metrics in a splicing mutant — **no
   comparator exists**.

**Contradictions the manuscript must address:**

1. Human UBL5 is reported **essential** (Ammon 2014, 10.1093/jmcb/mju026) — a viable KO needs
   documenting.
2. Human UBL5 knockdown causes **global baseline intron retention** (Oka 2014,
   10.15252/embr.201439478) — directly against "no constitutive splicing defect."
3. The one direct functional test of IRAK3 splice variants was **largely negative** (Rebl 2019,
   10.3389/fimmu.2019.02246) — argues against a dominant-negative IRAK3 isoform.
4. Spliceosome perturbation can **induce** type I interferon (Chang 2021, 10.1016/j.jbc.2021.101277)
   — makes baseline ISG levels a decisive control.
5. About **half** of macrophage LPS splicing changes are alternative-first-exon events
   (Robinson 2021, 10.7554/eLife.69431) — promoter choice, not spliceosome catalysis; undermines a
   UBL5 attribution unless the response is stratified by event class.

**Highest-value additional analyses, in order:**

1. Stratify the LPS isoform-switch response by AS event class; re-compute WT-vs-KO retention within
   spliceosome-plausible classes (IR, A5/A3, ES) versus ATSS/AFE.
2. Regress splicing-response retention on gene abundance; test whether LPS/IFN "sparing" survives.
3. Report baseline ISG expression in the KO.
4. Document the UBL5 allele and residual protein.
5. Report the abundance distribution of cryptic switchers versus all switchers.
