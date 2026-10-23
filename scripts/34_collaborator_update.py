#!/usr/bin/env python3
"""
34_collaborator_update.py -- the document that goes to the bench.

docs/ASSAY_DESIGNS.md is the complete reference, 1,269 lines, one entry per candidate.
This is the thing a collaborator actually reads: what the candidates are, which assay each
needs, what to order, and what to run. It names a first batch rather than handing over 44
candidates at once, and the batch is chosen by an explicit rule, not by taste.

Outputs:
  docs/COLLABORATOR_UPDATE.md
  results/tables/order_batch1.csv      oligos for the first batch, ready to paste into an order
"""
import csv, os
from collections import defaultdict, Counter

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
os.chdir(ROOT)
rd = lambda p: list(csv.DictReader(open(p)))
B = lambda v: str(v).strip().lower() == "true"
F = lambda v: float(v) if v not in ("", None) else 0.0

assign = rd("results/tables/assay_assignment.csv")
plat = {r["isoform_id"]: r for r in rd("results/tables/shortlist_platform.csv")}
regions = {r["isoform_id"]: r for r in rd("results/tables/discriminating_regions.csv")}
order = defaultdict(list)
for o in rd("results/tables/assay_order_sheet.csv"):
    order[o["isoform_id"]].append(o)
prep = defaultdict(dict)
for p in rd("results/tables/pair_selection_report.csv"):
    prep[p["isoform_id"]][p["pair"]] = p

# ------------------------------------------------------------------ batch rule
# A candidate enters batch 1 when all three hold:
#   1. the oligos for its PRIMARY assay carry no off-target locus
#   2. it has no QC flag (not a short model, not a barely-overlapping pair)
#   3. it is either switching in BOTH THP-1 and primary MDM, or it has a complete
#      discriminating qPCR assay (isoform-specific pair AND a valid total-gene denominator)
for a in assign:
    iso = a["isoform_id"]
    ol = order.get(iso, [])
    primary = [o for o in ol if not o["oligo_role"].startswith(("longread", "qPCR_total"))]
    a["_clean"] = bool(primary) and all(str(o["offtarget_loci"]) == "0" for o in primary)
    r = regions.get(iso, {})
    a["_flag"] = B(r.get("low_overlap_flag")) or B(r.get("short_model_flag"))
    q = prep[iso].get("qPCR_F/qPCR_R")
    t = prep[iso].get("qPCR_total_F/qPCR_total_R")
    a["_qpcr_complete"] = bool(q and q["passes"] == "True"
                               and t and t["passes"] == "True" and B(t["amplifies_partner"]))
    a["_batch1"] = a["_clean"] and not a["_flag"] and (a["set"] == "shared" or a["_qpcr_complete"])

batch1 = sorted([a for a in assign if a["_batch1"]], key=lambda x: -abs(F(x["own_dIF"])))
rest = sorted([a for a in assign if not a["_batch1"]], key=lambda x: -abs(F(x["own_dIF"])))

out = []
w = out.append

w("# Isoform switching in the LPS response: candidates, assays, and what to run\n")
w("Prepared from the long-read (PacBio) isoform analysis of THP-1 macrophages and primary "
  "human monocyte-derived macrophages, ± LPS. This is the bench-facing summary; the "
  "complete per-candidate reference is `docs/ASSAY_DESIGNS.md`, and every number below "
  "comes from the tables in `results/tables/`.\n")

w("## The short version\n")
n_pb = sum(1 for a in assign if plat[a["isoform_id"]]["platform_required"].startswith("PacBio"))
w(f"We have **{len(assign)} candidate isoform switches** that survived a direction-of-effect "
  f"filter and a locus-fusion screen. Of these, **{n_pb} are visible only because we used "
  f"long reads** and {len(assign)-n_pb} could in principle have been found by standard "
  f"short-read RNA-seq.\n")
w(f"The central practical finding is that **most of these are not qPCR problems.** "
  f"{Counter(a['primary_event'] for a in assign)['alternative first exon (TSS)']} of the "
  f"{len(assign)} switches are alternative-first-exon events — the two isoforms differ in "
  f"where transcription *starts*, and share everything downstream. For those, the shorter "
  f"isoform is contained inside the longer one, so no short amplicon can be specific to it. "
  f"They need an assay that reads the 5' end: RACE, or CAGE.\n")
w(f"I am proposing a **first batch of {len(batch1)} candidates** (below) rather than all "
  f"{len(assign)}, chosen by an explicit rule: clean reagents, no quality flags, and either "
  f"reproducible across both cell systems or supported by a complete qPCR assay. The other "
  f"{len(rest)} are ready in the same format whenever you want them.\n")

# ------------------------------------------------------------------ 1. candidates
w("## 1. The candidates\n")
w("**dIF** is the change in that isoform's share of its gene's total output on LPS. "
  "A dIF of +0.30 means the isoform went from, say, 20% to 50% of the gene. It is a "
  "*fraction*, not an expression change — which is why every qPCR below needs a total-gene "
  "denominator as well as an isoform-specific amplicon.\n")
w("| gene | dIF | direction on LPS | seen in | needs long reads? | what differs | assay |")
w("|---|---|---|---|---|---|---|")
for a in batch1:
    p = plat[a["isoform_id"]]
    d = F(a["own_dIF"])
    w(f"| **{a['gene_name']}** | {d:+.3f} | {'up' if d>0 else 'down'} | "
      f"{'both systems' if a['set']=='shared' else a['discovered_in']} | "
      f"{'yes' if not B(p['illumina_findable']) else 'no'} | "
      f"{a['primary_event'].replace(' (TSS)','').replace(' (3 end)','')} | "
      f"{a['primary_assay'].split(' (')[0]} |")
w("")
w("Direction is the prediction you are testing. *Up* means that isoform's share of the gene "
  "should be higher in the LPS-treated sample.\n")

shared = [a for a in batch1 if a["set"] == "shared"]
if shared:
    w(f"**{', '.join(a['gene_name'] for a in shared)}** switch in both THP-1 and primary MDM "
      f"with all three donors agreeing in direction. Those are the ones I would bet on.\n")

# ------------------------------------------------------------------ 2. assays
w("## 2. Which assay, and why that one\n")
w("| assay | candidates in batch 1 | what it reads out | why this one |")
w("|---|---|---|---|")
rows = [
    ("5' RACE (nested)", "the exact 5' end of each transcript, and the ratio between two "
     "start sites in one reaction",
     "the two isoforms differ only in where they start, so the short one is contained in "
     "the long one and no amplicon can be specific to it"),
    ("Isoform-specific qPCR / ddPCR", "the abundance of one isoform, expressed as a "
     "fraction of the gene's total",
     "these have a stretch of sequence unique to one isoform that a primer can sit in — "
     "the quickest readout in the set"),
    ("3' RACE", "the exact 3' end / polyA site used",
     "the isoforms differ at the 3' end and the unique region is too short for a qPCR"),
    ("Targeted long-read cDNA", "the whole transcript as one molecule, so every difference "
     "is read at once", "confirmatory for everything here, and the only route at all for an "
     "isoform fully nested inside a sibling"),
]
assay_counts = Counter(a["primary_assay"].split(" (")[0] for a in batch1)
for name, reads, why in rows:
    key = next((k for k in assay_counts if k.lower().startswith(name.split()[0].lower())), None)
    n = assay_counts.get(key, 0) if key else 0
    if name.startswith("Isoform-specific"):
        n = sum(1 for a in batch1 if "qPCR" in a["primary_assay"])
    if name.startswith("Targeted"):
        n = len(batch1)
    if n == 0:
        continue          # don't list an assay no batch-1 candidate needs
    w(f"| {name} | {n} | {reads} | {why} |")
w("")
w("**The single highest-value experiment is CAGE (or nanoCAGE).** It maps transcription "
  "start sites genome-wide, so one experiment per condition would answer all "
  f"{Counter(a['primary_event'] for a in assign)['alternative first exon (TSS)']} "
  "alternative-first-exon switches at once, instead of running that many separate RACE "
  "reactions. It would also settle a standing question about whether these alternative 5' "
  "ends are genuinely capped. If there is any route to CAGE through a core facility or a "
  "collaborator, that is where I would put the effort first; the RACE designs below are the "
  "per-candidate fallback and are worth running regardless for the handful we care most about.\n")

# ------------------------------------------------------------------ 3. order
w("## 3. What to order\n")
w("### Oligos\n")
b1iso = {a["isoform_id"] for a in batch1}
b1_oligos = []
for a in batch1:
    for o in sorted(order[a["isoform_id"]], key=lambda x: x["oligo_role"]):
        b1_oligos.append(o)
w(f"**{len(b1_oligos)} oligos** for the {len(batch1)} batch-1 candidates, listed in "
  "`results/tables/order_batch1.csv` in a form you can paste straight into an order. "
  "Standard desalted, 25 nmol scale is fine for all of them — these are 20–25-mers and none "
  "needs HPLC purification.\n")
w("Every oligo has been checked against all 296,613 transcripts in our reference by exact "
  "match on its 3'-terminal 15 nt, extended to full length. In the batch-1 set, every "
  "primary-assay oligo binds only within its own gene.\n")
w("Two columns in the order sheet report specificity and they mean different things. "
  "`offtarget_loci` counts binding sites for that **single primer**. The `selection` column "
  "reports the verdict for the **pair**, which only amplifies where both primers bind in the "
  "right orientation and spacing — so a primer with one stray site can still belong to a pair "
  "that amplifies nothing else. Where they appear to disagree, the pair verdict is the one "
  "that governs whether you get a spurious band.\n")
w("Each candidate gets up to four kinds of oligo, and they do different jobs:\n")
w("| role | job |")
w("|---|---|")
w("| `qPCR_F` / `qPCR_R` | isoform-specific amplicon — the numerator |")
w("| `qPCR_total_F` / `qPCR_total_R` | sits in sequence common to both isoforms — the "
  "denominator, so you can report a fraction rather than an expression change |")
w("| `RACE5_GSP1_outer` / `RACE5_GSP2_inner` | nested antisense primers for 5' RACE; outer "
  "for the first PCR, inner for the nested PCR |")
w("| `longread_F` / `longread_R` | outer pair spanning the whole transcript, for the "
  "confirmatory long-read amplicon |")
w("")

w("### Kits and reagents\n")
w("- **RNA prep**: your usual column prep plus an on-column or in-solution DNase step. The "
  "no-RT control below is the check that this worked.\n")
w("- **Reverse transcription**: a standard RT kit for the qPCR work. For 5' RACE you need a "
  "*template-switching* RT kit specifically — the SMARTer RACE 5'/3' kit (Takara) is the "
  "common one, or any equivalent that adds a defined adapter at the 5' end. A conventional "
  "RT kit will not work for RACE.\n")
w("- **qPCR**: SYBR-based master mix. If you would rather have absolute numbers and skip "
  "standard curves, ddPCR with EvaGreen works and is arguably better suited to a ratio "
  "measurement.\n")
w("- **Gel / cleanup**: 1.5–2% agarose for sizing RACE products, plus a gel-extraction kit — "
  "you will want to sequence the bands, not just photograph them.\n")
w("- **Sequencing**: Sanger for individual gel-extracted RACE bands. For the long-read "
  "confirmation, a PacBio or ONT amplicon run through a core facility; the amplicons are "
  "0.5–9 kb.\n")

w("### What to ask a core facility for\n")
w("- CAGE or nanoCAGE on the same RNA, ± LPS (see above — highest value per experiment).\n")
w("- PacBio or ONT amplicon sequencing of the long-read outer products, barcoded by condition.\n")

# ------------------------------------------------------------------ 4. run
w("## 4. What to run\n")
w("### Stimulation and RNA — read this first\n")
w("> **The validation has to use the same stimulation conditions as the sequencing.** The PMA "
  "differentiation protocol, the LPS dose and the timepoint are not recorded anywhere in the "
  "analysis repository, so they need to come from the original experiment records before "
  "anyone starts. If the validation runs at a different dose or timepoint and a switch fails "
  "to reproduce, there is no way to tell whether the switch is wrong or the conditions "
  "were.\n")
w("Three independent biological replicates per condition, −LPS and +LPS. For the primary-MDM "
  "candidates that means three donors. Process both conditions of a replicate in the same "
  "batch: these are ratio measurements, and a ratio compared across batches is not "
  "interpretable.\n")

w("### Protocol A — isoform-specific qPCR (the quickest answer)\n")
w("1. DNase-treat RNA; reverse-transcribe with random hexamers or oligo-dT.\n")
w("2. For each candidate run **two** assays on the same cDNA: the isoform-specific pair and "
  "the total-gene pair.\n")
w("3. Validate both pairs first: standard curve over at least 5 logs, efficiency 90–110%, and "
  "a single melt peak. **A fraction built from two pairs of different efficiency is not a "
  "fraction** — this step is not optional.\n")
w("4. Compute the isoform fraction as 2^−(Ct_isoform − Ct_total) per sample. Do **not** "
  "normalise to a housekeeping gene: the claim is about the ratio of two isoforms within a "
  "gene, and total-gene normalisation is what makes it a fraction.\n")
w("5. Compare the fraction between −LPS and +LPS across the three replicates.\n")
w("**Controls**: no-RT and no-template for every pair; melt curve every run.\n")

w("### Protocol B — nested 5' RACE (for the start-site switches)\n")
w("1. Template-switching RT on 1 µg total RNA, both conditions, same batch.\n")
w("2. First PCR: universal adapter primer + the **outer** gene-specific primer (`GSP1`).\n")
w("3. Nested PCR: nested universal primer + the **inner** primer (`GSP2`), on a 1:50 dilution "
  "of PCR 1.\n")
w("4. Run on 1.5–2% agarose. The per-candidate entries give the expected size from each "
  "isoform — for most of these the two start sites differ by kilobases, so the two products "
  "are far apart and easy to resolve.\n")
w("5. **Gel-extract both bands and sequence them.** A size shift is suggestive; only the "
  "sequence tells you which start site was used. If you want a real ratio rather than band "
  "intensity, sequence the un-separated nested product as an amplicon library and count reads "
  "per start site — that is a much better measurement than densitometry.\n")
w("**Controls**: no-RT (a band here means genomic DNA and voids the TSS call for that "
  "sample); no-template in both PCRs; a positive-control GSP against a single-TSS gene to "
  "confirm adapter ligation worked.\n")

w("### Protocol C — targeted long-read cDNA (the confirmation, and the only option for some)\n")
w("1. Amplify with the outer pair, both conditions, barcoded.\n")
w("2. Check the product on a Bioanalyzer/TapeStation before sequencing — a dominant short "
  "product will eat the run.\n")
w("3. Sequence on PacBio or ONT and count full-length reads per isoform.\n")
w("4. Reads shorter than the expected amplicon are PCR or library artefacts, not short "
  "isoforms — do not interpret them as evidence for a new transcript.\n")
w("This is the only assay that can see an isoform that is fully nested inside a sibling, and "
  "it is also the cleanest way to get a per-isoform ratio for anything in the list.\n")

w("### How to call a result\n")
w("For each candidate the per-candidate entry gives a specific threshold, but the rule is the "
  "same throughout:\n")
w("- **Supports the switch**: the isoform's share moves in the predicted direction by at "
  "least 0.16 (THP-1) or 0.17 (primary MDM), across three replicates.\n")
w("- **Falsifies it**: the share moves in the *opposite* direction.\n")
w("- **Uninformative, not negative**: a shift smaller than that threshold. Those numbers are "
  "the minimum detectable change at n = 3 with 80% power, so a smaller result means the "
  "experiment could not have resolved it either way. The honest response is more replicates, "
  "not a negative conclusion.\n")

# ------------------------------------------------------------------ 5. caveats
w("## 5. What is not ready, and what I would not trust yet\n")
w("- **Six of the 15 isoform-specific qPCR designs do not work as discriminating assays** and "
  "are excluded from batch 1. RBKS, TMEM177 and KDM6B amplify the partner isoform as well, so "
  "they measure the gene rather than the switch (TMEM177 gives the same 133 bp product from "
  "both). ZP3, PHLDA3 and POLD4 have primers that bind a second gene, and that survived all "
  "five alternatives I tried.\n")
w("- **ACP3's long-read reverse primer binds 4,355 loci** — a repeat element, not a bad pick. "
  "Its qPCR assay is fine; only the long-read confirmation needs a hand-designed primer.\n")
w("- **Two candidates have structural quality flags.** KDM6B's candidate transcript is a "
  "219 nt model that shares only 72 nt with its 6,728 nt partner, which looks like an assembly "
  "fragment rather than a real isoform. A few others (JDP2, SPRED1, GATD3) share so little "
  "sequence with their partner that whether they are two isoforms of one gene is a fair "
  "question. All are flagged in the per-candidate entries.\n")
w("- **Antibody and mass-spec routes are specified but not validated.** Where a per-candidate "
  "entry calls for an isoform-resolving western, it states where the epitope has to fall for "
  "the antibody to tell the two isoforms apart. Whether such an antibody actually exists has "
  "not been checked for any of them.\n")
w("- **These oligos are designed and computationally screened, not bench-tested.** The "
  "specificity numbers are exact-match counts against a transcriptome, which is a much better "
  "check than a thermodynamic score, but it is not a gel.\n")

os.makedirs("docs", exist_ok=True)
open("docs/COLLABORATOR_UPDATE.md", "w").write("\n".join(out) + "\n")

cols = ["gene_name", "oligo_role", "sequence", "length", "tm_C", "gc_pct",
        "expected_product_candidate_bp", "expected_product_partner_bp",
        "offtarget_loci", "assay", "isoform_id", "partner", "selection"]
with open("results/tables/order_batch1.csv", "w", newline="") as fh:
    wtr = csv.DictWriter(fh, fieldnames=cols, extrasaction="ignore")
    wtr.writeheader()
    wtr.writerows(b1_oligos)

print(f"docs/COLLABORATOR_UPDATE.md written: {len(out)} lines")
print(f"batch 1: {len(batch1)} candidates, {len(b1_oligos)} oligos")
print("  " + ", ".join(a["gene_name"] for a in batch1))
print(f"deferred: {len(rest)} candidates")
