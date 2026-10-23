#!/usr/bin/env python3
"""
35_protocols_doc.py -- generate docs/PROTOCOLS.md: step-by-step bench protocols for the
designed assays, anchored to the protocol actually used for the original RNA-seq.

Source of the upstream conditions is the thesis, Katz S., "Thesis_final_sk.pdf", sections
cited inline per step. Nothing upstream is invented here: cell culture, differentiation,
LPS source/dose/duration and RNA extraction are reproduced from the original so that a
validation run is comparable to the sequencing it is validating. Where a protocol step has
NO precedent in the thesis (RACE, long-read amplicons, SYBR detection) that is stated
explicitly and the source is named.

Oligo tables are injected from results/tables/order_batch1.csv so the document cannot drift
from the order sheet.
"""
import csv, os
from collections import defaultdict

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
os.chdir(ROOT)

oligos = list(csv.DictReader(open("results/tables/order_batch1.csv")))
by_gene = defaultdict(dict)
for o in oligos:
    by_gene[o["gene_name"]][o["oligo_role"]] = o
genes = sorted(by_gene)

out = []
w = out.append

# ---------------------------------------------------------------- header
w("# Bench protocols for isoform-switch validation\n")
w("Step-by-step protocols for the assays designed against the long-read candidate list. "
  "Written to be run alongside, and to be directly comparable with, the RNA-seq that "
  "generated the candidates.\n")
w("**Companion documents.** `COLLABORATOR_UPDATE.md` says which candidates and why; "
  "`ASSAY_DESIGNS.md` gives the per-candidate rationale, expected product sizes and "
  "falsification criteria; this document is how to run it. Oligo sequences are in "
  "`results/tables/order_batch1.csv`.\n")

w("## 0. Where each parameter comes from\n")
w("The single most important property of this protocol is that **the upstream half is not "
  "new**. Cell source, differentiation, LPS and RNA extraction are reproduced from the "
  "original experiment, because a validation run under different conditions cannot "
  "distinguish a wrong candidate from a wrong experiment.\n")
w("| step | source | status |")
w("|---|---|---|")
w("| THP-1 culture and PMA differentiation | Thesis §2.5.1 | reproduced verbatim |")
w("| Primary human macrophage isolation and differentiation | Thesis §2.5.6 | reproduced verbatim |")
w("| LPS source, dose, duration | Thesis §2.6.1; §7.2; Fig. 7.1 | reproduced verbatim |")
w("| RNA extraction | Thesis §2.6.8 | reproduced verbatim |")
w("| Reverse transcription and qPCR cycling | Thesis §2.6.7 | reproduced, with one documented change (§3.2) |")
w("| Semi-quantitative RT-PCR on gel | Thesis §2.6.6; precedent in §6.8 | reproduced verbatim |")
w("| TNF-α / HPRT stimulation control | Thesis §2.6.7, Table 2.1; Fig. 7.2 | reproduced verbatim, primers included below |")
w("| Isoform-specific and total-gene primer pairs | this project, `scripts/28`–`30` | NEW — designed and screened *in silico*, not bench-tested |")
w("| Nested 5′ RACE | no precedent in the thesis | NEW — kit manufacturer's protocol, cited in §4 |")
w("| Targeted long-read amplicon | no precedent in the thesis | NEW — see §5 |")
w("")
w("> **Precedent worth knowing about.** The thesis already contains an experiment of exactly "
  "this type. In §6.8, primer sets spanning versus flanking the alternatively spliced region "
  "of *IRAK1* (exon 11; Rao et al., 2005) and of *TMED7-TICAM2* (Doyle et al., 2012) were run "
  "by RT-PCR on THP-1 RNA at 0, 2 and 4 h of LPS and resolved on a gel (Figs. 6.13, 6.14). "
  "The spanning primers lost signal with LPS while the flanking primers did not. Protocol D "
  "below is that same design, and the primer-naming convention used here (S = sense, "
  "A = antisense) follows Table 2.1.\n")

# ---------------------------------------------------------------- 1 upstream
w("## 1. Cells, stimulation and RNA\n")
w("### 1.1 THP-1–derived macrophages — Thesis §2.5.1\n")
w("1. Maintain THP-1 in RPMI-1640 + 10% FBS + 2 mM L-glutamine.\n")
w("2. Differentiate to a macrophage-like state with **5 ng/ml phorbol-12-myristate-13-acetate "
  "(PMA; Sigma P1585) for 72 h**.\n")
w("   *Note: this is a deliberately low PMA dose and a long differentiation. Do not "
  "substitute a higher-dose/shorter protocol — PMA dose changes the LPS response.*\n")
w("3. Plate for RNA at **4.0 × 10⁷ cells per 10 cm dish** (Thesis §2.6.8). For qPCR-only "
  "validation this can be scaled down; see §1.4.\n")

w("### 1.2 Primary human monocyte-derived macrophages — Thesis §2.5.6\n")
w("1. Obtain peripheral blood monocytes from screened healthy donors. In the original work "
  "these came from the NIH Department of Transfusion Medicine under NIH Clinical Center "
  "**IRB-approved protocol 99-CC-0168**. *Any new collection needs its own approved protocol "
  "— cite yours here before running.*\n")
w("2. Propagate in RPMI + 10% FBS + 10 mM HEPES + β-mercaptoethanol, 37 °C, 5% CO₂.\n")
w("3. Differentiate with **10 ng/ml GM-CSF (R&D Systems) for 7 days**.\n")
w("4. **One donor is one biological replicate.** The original design used three donors, "
  "each collected on a different day. This is the dominant source of variance in the primary "
  "arm, and it is why the detection threshold for primary cells (§7) is higher than for "
  "THP-1.\n")

w("### 1.3 LPS stimulation — Thesis §2.6.1, §7.2, Fig. 7.1\n")
w("1. LPS: *Salmonella minnesota* R595, TLRgrade, **Enzo Life Sciences ALX-581-008-L002** "
  "(formerly Alexis Biochemicals).\n")
w("2. Treat at **100 ng/ml for 4 h**. Leave a paired subset untreated.\n")
w("3. **Both arms of a replicate must be handled together.** Every readout here is a ratio "
  "within a gene; a ratio compared across handling batches is not interpretable.\n")
w("> R595 is a rough-chemotype LPS and is not interchangeable with smooth *E. coli* "
  "O111:B4 or with a different supplier's preparation — potency and TLR4/MD-2 engagement "
  "differ. If the original lot is gone, note the substitution in the record; it is a "
  "candidate explanation for any failure to reproduce.\n")

w("### 1.4 RNA extraction — Thesis §2.6.8\n")
w("1. At the timepoint, remove media and wash twice with PBS.\n")
w("2. Lyse in **1 ml Buffer RLT : β-mercaptoethanol at 100:1** per 10 cm dish.\n")
w("3. Plates may be held at **−80 °C** until extraction.\n")
w("4. Extract total RNA with the **RNeasy Mini Kit (Qiagen 74106)** per the manufacturer's "
  "instructions.\n")
w("5. **Add a DNase step.** The original protocol as written does not include one. It is "
  "required here because the no-RT control is load-bearing for both RACE and the "
  "junction-spanning assays — on-column RNase-free DNase (Qiagen 79254) or an equivalent "
  "in-solution treatment.\n")
w("6. Check integrity. The sequencing used **RIN ≥ 9.6** (Thesis §2.6.10); for RACE, which "
  "depends on intact 5′ ends, aim for the same. For qPCR alone RIN ≥ 8 is adequate.\n")
w("> **Scale.** 4 × 10⁷ cells per plate yields far more RNA than these assays need. qPCR "
  "needs ~1 µg per RT reaction; RACE needs ~1 µg. A 6-well format is sufficient if you are "
  "not also banking RNA — but keep the cell density per cm² the same, because macrophage "
  "density alters the LPS response.\n")

# ---------------------------------------------------------------- 2 QC
w("## 2. Stimulation control — run this first, every time\n")
w("Before any isoform assay, confirm the LPS worked, using the assay the original work used "
  "for exactly this purpose (Thesis §2.6.7, Table 2.1; Fig. 7.2 shows the expected result "
  "across all three replicates of each cell type).\n")
w("Measure **TNF-α induction normalised to HPRT**, treated versus untreated, per replicate.\n")
w("| target | oligonucleotide (5′→3′) |")
w("|---|---|")
w("| Human TNF-α, forward (sense) | `CCAGGGACCTCTCTCTAATCA` |")
w("| Human TNF-α, reverse (antisense) | `TCAGCTTGAGGGTTTGCTAC` |")
w("| Human TNF-α, FAM probe | `AGTGACAAGCCTGTAGCCCATGTT` |")
w("| Human HPRT, forward (sense) | `CTGGAAAGAATGTCTTGATTGTGG` |")
w("| Human HPRT, reverse (antisense) | `CTTGCGACCTTGACCATCTT` |")
w("| Human HPRT, FAM probe | `AGACTTTGCTTTCCTTGGTCAGGCA` |")
w("| Human β-actin, forward (sense) | `TGAAGTCTGACGTGGACATC` |")
w("| Human β-actin, reverse (antisense) | `ACTCGTCATACTCCTGCTTG` |")
w("")
w("*Sequences transcribed from Thesis Table 2.1 (§2.6.7).*\n")
w("**Gate:** proceed only for replicates showing clear TNF-α induction over the paired "
  "untreated sample. A replicate that did not respond to LPS cannot inform an isoform-ratio "
  "question, and including it will only add variance.\n")

# ---------------------------------------------------------------- 3 qPCR
w("## 3. Protocol A — isoform-specific qPCR / ddPCR\n")
w("**Applies to:** the candidates whose primary assay is a qPCR in `COLLABORATOR_UPDATE.md` "
  "§1. **Reads out:** the candidate isoform's share of its gene's total output.\n")
w("### 3.1 Reverse transcription — Thesis §2.6.7\n")
w("1. Reverse-transcribe **1 µg total RNA** with the **Bio-Rad reverse transcription kit**. "
  "The thesis specifies the supplier but not which Bio-Rad kit — confirm the exact product "
  "from the lab records so the RT chemistry matches the original.\n")
w("2. Use **200 ng RNA-equivalent of cDNA per reaction**.\n")
w("3. Prime with **random hexamers**, not oligo-dT alone. Several of these amplicons sit far "
  "from the 3′ end — ARHGAP26's transcript is 9.2 kb — and oligo-dT priming will "
  "under-represent them. This is a change from a 3′-biased default and is deliberate.\n")

w("### 3.2 Detection chemistry — the one real deviation\n")
w("> The original qPCR was **probe-based**: gene-specific primers with FAM-conjugated probes "
  "(IDT) and **Solaris qPCR mix** (Dharmacon/Life Technologies), on a Bio-Rad CFX "
  "(Thesis §2.6.7). **The primer pairs designed in this project are primer-only and assume "
  "an intercalating dye.** Two ways forward:\n")
w("> **(a) Run them as SYBR/EvaGreen assays** (recommended to start). Requires an "
  "intercalating-dye master mix rather than Solaris, and makes the melt curve in §3.4 "
  "mandatory rather than optional. Same instrument, same cycling.\n")
w("> **(b) Keep the probe workflow** and have a FAM probe designed per amplicon. For the "
  "isoform-specific assays the probe must sit over the discriminating sequence or junction, "
  "or it will not add specificity over the primers. This is the better long-term format and "
  "matches the lab's existing reagents — ask for the probes to be designed if you want to go "
  "this route.\n")

w("### 3.3 Reactions to set up — two per candidate, not one\n")
w("For each candidate run **both** assays on the same cDNA:\n")
w("1. **Isoform-specific pair** (`qPCR_F` / `qPCR_R`) — the numerator. At least one primer "
  "sits in sequence unique to the candidate isoform.\n")
w("2. **Total-gene pair** (`qPCR_total_F` / `qPCR_total_R`) — the denominator. Both primers "
  "sit inside one block of sequence common to both isoforms, so it reports the locus's total "
  "output.\n")
w("> **Why both.** dIF is a *fraction*. LPS moves total output at several of these loci hard "
  "(*RSAD2* and *USP18* are interferon-induced). An isoform-specific Ct on its own cannot "
  "distinguish \"this isoform's share rose\" from \"the whole gene was induced\" — and only "
  "the first is the claim being tested. Do **not** substitute HPRT or β-actin "
  "normalisation here: a housekeeping gene answers a different question. HPRT belongs in §2, "
  "as the stimulation control.\n")

w("### 3.4 Assay validation — before any biological comparison\n")
w("1. **Standard curve** for each pair over ≥5 logs of input cDNA. Require efficiency "
  "90–110% and R² ≥ 0.99.\n")
w("2. **Melt curve** on every reaction: a single peak. These primers were screened against "
  "all 296,613 transcripts in the long-read reference, but a second melt peak means a second "
  "product regardless of what the screen said.\n")
w("3. **Efficiency matching is not optional.** The fraction is a ratio of two assays; if the "
  "numerator and denominator pairs differ in efficiency, the ratio is a systematic artefact "
  "that will look like a switch. If a pair fails, use one of the alternatives in "
  "`results/tables/assay_reagents.csv` (up to five were designed and screened per assay) "
  "rather than proceeding with a poor curve.\n")

w("### 3.5 Cycling — Thesis §2.6.7\n")
w("```\n95 °C  15 min\n40 ×  ( 95 °C  15 s  →  60 °C  60 s )\n[+ melt curve 65→95 °C if using an "
  "intercalating dye]\n```\n")
w("All designed amplicons are 70–200 bp, compatible with this two-step programme.\n")

w("### 3.6 Calculation and call\n")
w("1. Per sample: `isoform fraction ∝ 2^−(Ct_isoform − Ct_total)`.\n")
w("2. Compare the fraction between paired −LPS and +LPS samples within each replicate.\n")
w("3. Apply the thresholds in §7.\n")
w("**ddPCR alternative.** If a ddPCR platform is available, it is the better fit: it gives "
  "absolute counts for numerator and denominator without standard curves, and a ratio of two "
  "absolute counts avoids the efficiency-matching problem in §3.4 entirely. Same primers, "
  "same amplicon sizes.\n")

# ---------------------------------------------------------------- 4 RACE
w("## 4. Protocol B — nested 5′ RACE\n")
w("**Applies to:** the alternative-first-exon candidates — the largest group, and the ones "
  "no qPCR can address, because the shorter isoform is wholly contained within the longer "
  "one. **Reads out:** the actual transcription start site used, and the ratio between two "
  "start sites in a single reaction.\n")
w("> **No precedent in the thesis.** The original work did not run RACE. The steps below "
  "follow the **template-switching RACE** method as implemented in the Takara SMARTer RACE "
  "5′/3′ kit; confirm the current part number with the supplier. Any equivalent kit that adds "
  "a defined adapter at the 5′ end during first-strand synthesis will work. A conventional "
  "RT kit — including the Bio-Rad kit used in §3 — **cannot** be used for RACE: there is no "
  "adapter to prime from.\n")

w("### 4.1 Why the primers are where they are\n")
w("Each candidate has two antisense gene-specific primers, placed in sequence **shared** with "
  "the partner isoform rather than in the candidate-specific region. That is deliberate: a "
  "shared-sequence GSP amplifies from *both* start sites in one reaction, so the two products "
  "differ in length by the distance between the TSSs and the readout is a ratio rather than a "
  "presence call. Expected sizes for both isoforms are in `ASSAY_DESIGNS.md` per candidate.\n")
w("`RACE5_GSP1_outer` is used in the first PCR; `RACE5_GSP2_inner` sits 5′ of it and is used "
  "in the nested PCR. Both were designed at Tm 64–71 °C and 22–30 nt, which satisfies the "
  "kit's GSP requirements (≥65 °C, 23–28 nt) without redesign.\n")

w("### 4.2 Steps\n")
w("1. **First-strand synthesis** on ~1 µg total RNA per condition with the kit's 5′-CDS "
  "primer and template-switching oligo. Process paired −LPS/+LPS samples side by side.\n")
w("2. **First PCR**: kit universal primer mix + **`RACE5_GSP1_outer`**. Use the kit's "
  "touchdown programme, setting extension time for the *longest* expected product (up to "
  "~3 kb for several of these):\n")
w("```\n 5 ×  ( 94 °C 30 s → 72 °C 3 min )\n 5 ×  ( 94 °C 30 s → 70 °C 30 s → 72 °C 3 min )\n"
  "25 ×  ( 94 °C 30 s → 68 °C 30 s → 72 °C 3 min )\n```\n")
w("3. **Nested PCR**: kit nested universal primer + **`RACE5_GSP2_inner`**, on a 1:50 "
  "dilution of the first PCR. 20–25 cycles at the final annealing temperature.\n")
w("4. **Resolve** on a 1.5–2% agarose gel. For most of these candidates the two start sites "
  "are kilobases apart — the median TSS separation across the candidate list is 8,680 bp — so "
  "the two products are far apart and easy to call.\n")
w("5. **Sequence the products.** Gel-extract each band and Sanger sequence it. A size shift "
  "is suggestive; only the sequence establishes which TSS was used. **This step is not "
  "optional** — a band of roughly the right size is the most common way a RACE result is "
  "over-interpreted.\n")
w("6. **For a quantitative ratio**, do not use densitometry. Take the un-separated nested "
  "product into an amplicon library and count reads per start site. Band intensity on a gel "
  "is a poor estimator of a ratio, and the ratio is the measurement.\n")

w("### 4.3 Controls — Thesis §2.6.8 for the RNA, kit manual for the rest\n")
w("1. **No-RT control** on the same RNA. A product here means genomic DNA carry-through and "
  "**voids every TSS call from that sample** — RACE products are primed from an adapter, so "
  "a gDNA template produces plausible-looking artefacts. This is why the DNase step in §1.4 "
  "was added.\n")
w("2. **No-template control** in both the first and the nested PCR.\n")
w("3. **Positive-control GSP** against a gene with a single well-annotated TSS, run in the "
  "same batch, to show adapter addition and template switching worked. Without it, a failed "
  "RACE and a genuinely absent 5′ end look identical.\n")
w("4. Both conditions of a replicate in the same batch, as in §1.3.\n")

# ---------------------------------------------------------------- 5 long read
w("## 5. Protocol C — targeted long-read cDNA amplicon\n")
w("**Applies to:** every candidate, as confirmation — and it is the **only** assay that can "
  "see an isoform fully nested inside a sibling as a distinct molecule. **Reads out:** the "
  "entire transcript in one read, so every difference between the isoforms is observed at "
  "once rather than inferred from one amplicon.\n")
w("1. Amplify from the same cDNA with the outer pair (`longread_F` / `longread_R`), which "
  "sits in the two termini of the transcript. Use a high-fidelity, long-range polymerase; "
  "products range from ~0.5 to ~9 kb depending on the candidate.\n")
w("2. Barcode by condition and replicate so one run covers the whole set.\n")
w("3. **Check the product distribution on a Bioanalyzer or TapeStation before sequencing.** "
  "A dominant short product will consume the run.\n")
w("4. Sequence on PacBio or ONT through a core facility. Count full-length reads per isoform.\n")
w("5. **Reads shorter than the expected amplicon are PCR or library artefacts, not short "
  "isoforms.** Do not interpret them as evidence for a new transcript — this is the failure "
  "mode that makes amplicon long-read data look more interesting than it is.\n")
w("> The original long-read data is described only as \"PacBio Whole-Genome De Novo "
  "Sequencing\" analysed with SQANTI (Thesis §2.6.11, §2.4.6; Tardaguila et al., 2018). The "
  "specific Iso-Seq chemistry and version are not recorded there. For a targeted amplicon run "
  "this does not need to match the original, but it should be recorded this time.\n")

# ---------------------------------------------------------------- 6 RT-PCR
w("## 6. Protocol D — semi-quantitative RT-PCR on gel\n")
w("**Applies to:** a fast, low-cost first look, and the only protocol here with a direct "
  "precedent in the thesis (§6.8, Figs. 6.13–6.14).\n")
w("The design logic is the one already used for *IRAK1* and *TMED7-TICAM2*: run **one primer "
  "set that spans the variable region** and **one that sits outside it**, on the same cDNA. "
  "The flanking set reports the locus's total output and acts as the internal reference; the "
  "spanning set reports the isoform. A change in the spanning product with no change in the "
  "flanking product is the signature of a splicing change rather than an expression change — "
  "which is exactly what was seen for *IRAK1* exon 11 at 2 and 4 h of LPS.\n")
w("Here, the `qPCR_F`/`qPCR_R` pair is the spanning set and `qPCR_total_F`/`qPCR_total_R` is "
  "the flanking set.\n")
w("Cycling — Thesis §2.6.6, Applied Biosystems Veriti 96-well:\n")
w("```\n95 °C  3 min\n35 ×  ( 95 °C 15 s → 55 °C 15 s → 72 °C 30 s )\n```\n")
w("Run products on a **2% agarose gel** (Thesis §2.6.6).\n")
w("> Treat this as qualitative. 35 cycles is deep into plateau for an induced transcript, so "
  "band intensity is not a quantitative ratio. Use it to see whether the event is there at "
  "all and whether the direction looks right, then quantify with Protocol A or C. The thesis "
  "used it the same way — as a demonstration that the event occurs within the LPS timeframe, "
  "not as a measurement.\n")

# ---------------------------------------------------------------- 7 calling
w("## 7. How to call a result\n")
w("Per-candidate thresholds are in `ASSAY_DESIGNS.md`; the rule is uniform:\n")
w("| outcome | criterion | interpretation |")
w("|---|---|---|")
w("| **Supports the switch** | isoform share moves in the predicted direction by ≥ **0.156** "
  "(THP-1) or ≥ **0.172** (primary MDM), across 3 replicates | reproduces the sequencing |")
w("| **Falsifies it** | share moves in the *opposite* direction | the call does not hold |")
w("| **Uninformative** | shift smaller than the threshold | underpowered, *not* negative |")
w("")
w("Those thresholds are the minimum detectable |dIF| at n = 3 with 80% power, computed from "
  "the replicate-level isoform-fraction variance in each arm of the original data. A result "
  "below them means the experiment could not have resolved the effect either way; the honest "
  "response is more replicates, not a negative conclusion. The primary-cell threshold is "
  "higher because donor-to-donor variance exceeds flask-to-flask variance — consistent with "
  "the concordance analysis in Thesis §7.3, where replicate agreement was highest in the "
  "most tightly controlled cell source.\n")

# ---------------------------------------------------------------- 8 oligos
w("## 8. Oligonucleotides\n")
w("Full table with expected product sizes, measured off-target counts and selection notes: "
  "`results/tables/order_batch1.csv`. Standard desalted, 25 nmol; no HPLC purification "
  "needed for any of them.\n")
for g in genes:
    d = by_gene[g]
    w(f"### {g}\n")
    w("| role | sequence (5′→3′) | nt | Tm °C | GC% | product, candidate | product, partner |")
    w("|---|---|---|---|---|---|---|")
    for role in ("qPCR_F", "qPCR_R", "qPCR_total_F", "qPCR_total_R",
                 "RACE5_GSP1_outer", "RACE5_GSP2_inner",
                 "RACE3_GSP1_outer", "RACE3_GSP2_inner",
                 "longread_F", "longread_R"):
        o = d.get(role)
        if not o:
            continue
        w(f"| `{role}` | `{o['sequence']}` | {o['length']} | {o['tm_C']} | {o['gc_pct']} | "
          f"{o['expected_product_candidate_bp']} | {o['expected_product_partner_bp']} |")
    w("")

# ---------------------------------------------------------------- 9 refs
w("## 9. Sources\n")
w("**Primary protocol source.** Katz, S. *Thesis_final_sk.pdf*. Sections cited inline: "
  "§2.4.6 (SQANTI), §2.5.1 (THP-1 culture and PMA), §2.5.6 (primary human macrophages), "
  "§2.6.1 (LPS), §2.6.6 (PCR), §2.6.7 and Table 2.1 (RT-qPCR, TNF-α/HPRT primers and probes), "
  "§2.6.8 (RNA extraction), §2.6.10 (short-read RNA-seq), §2.6.11 (long-read RNA-seq), "
  "§6.8 and Figs. 6.13–6.14 (RT-PCR validation of *IRAK1* and *TMED7-TICAM2* splicing), "
  "§7.2 and Figs. 7.1–7.2 (experimental design, LPS confirmation), §7.3 (replicate "
  "concordance).\n")
w("**Reagents named in the original protocol.** PMA, Sigma P1585. LPS, *S. minnesota* R595 "
  "TLRgrade, Enzo Life Sciences ALX-581-008-L002. RNeasy Mini Kit, Qiagen 74106. GM-CSF, "
  "R&D Systems. Solaris qPCR mix, Dharmacon/Life Technologies. FAM probes, IDT. Bio-Rad "
  "reverse transcription kit and CFX Real-Time System. Applied Biosystems Veriti 96-well "
  "thermocycler.\n")
w("**Reagents introduced here.** RNase-free DNase (Qiagen 79254 or equivalent) — §1.4. "
  "Template-switching 5′/3′ RACE kit (Takara SMARTer RACE or equivalent; confirm part number) "
  "— §4. Intercalating-dye qPCR master mix, if following §3.2(a). Long-range high-fidelity "
  "polymerase — §5.\n")
w("**Literature cited.** Rao, N. et al. (2005) — *IRAK1* exon 11 variants. Doyle, S. et al. "
  "(2012) — *TMED7-TICAM2* splice variants. Tardaguila, M. et al. (2018) — SQANTI. "
  "Shen, S. et al. (2014) — rMATS. Li, N. et al. (2015) — THP-1 B5 and RAW264.7 G9 reporter "
  "lines. All as cited in the thesis.\n")
w("**Designed in this project.** `scripts/28_assay_reagent_design.py` (primer3-py 2.4.0), "
  "`scripts/29_assay_specificity.py` (binding sites across all 296,613 transcripts of the "
  "H_T reference), `scripts/30_select_oligos.py` (selection on the measurement), "
  "`scripts/33_platform_for_shortlist.py` (novel-vs-known from the GFF3).\n")
w("> **Not validated at the bench.** Every oligo below was designed *in silico* and screened "
  "for specificity against the project's own long-read transcriptome. None has been run. The "
  "standard curves and melt curves in §3.4 and the controls in §4.3 are where that gets "
  "established.\n")

os.makedirs("docs", exist_ok=True)
open("docs/PROTOCOLS.md", "w").write("\n".join(out) + "\n")
print(f"docs/PROTOCOLS.md written: {len(out)} lines")
print(f"oligo tables for {len(genes)} genes, {len(oligos)} oligos")
