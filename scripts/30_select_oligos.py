#!/usr/bin/env python3
"""
30_select_oligos.py -- choose, per candidate, the alternative that the transcriptome
scan actually supports, and write the order sheet.

28 designs several alternatives per assay; 29 measures where each one binds across all
296,613 transcripts. This script applies the selection rule, which differs by assay
because the assays want different things:

  discriminating qPCR/ddPCR : must amplify the candidate, must NOT amplify the partner,
                              must not amplify anything outside the locus
  targeted long-read cDNA   : must not amplify anything outside the locus, and SHOULD
                              amplify the partner too -- the point is to capture both
                              isoforms in one amplicon pool and let the reads separate them
  5'/3' RACE GSP            : must not bind outside the locus; binding the partner is
                              preferred, because then one reaction sizes both ends

An on-target miss is ignored when the transcript's own reference record is one of the
36,115 unspliced ones: a junction-spanning primer cannot be found in a record that still
contains its introns, so that is a defect of the reference rather than of the primer.

Output: results/tables/assay_order_sheet.csv
"""
import csv, os
from collections import defaultdict

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
os.chdir(ROOT)

spec = list(csv.DictReader(open("results/tables/assay_specificity.csv")))
pair = list(csv.DictReader(open("results/tables/assay_pair_specificity.csv")))
reag = list(csv.DictReader(open("results/tables/assay_reagents.csv")))
assign = {r["isoform_id"]: r for r in csv.DictReader(open("results/tables/assay_assignment.csv"))}

B = lambda v: str(v).strip().lower() == "true"
I = lambda v: int(v) if str(v).strip().isdigit() else 0

# per-oligo specificity, keyed by the oligo's identity
sp = {}
for s in spec:
    sp[(s["isoform_id"], s["oligo_role"], s.get("alt", "0"))] = s
rg = {}
for r in reag:
    rg[(r["isoform_id"], r["oligo_role"], r.get("alt", "0"))] = r

chosen = []


def take(iso, role, alt, assay, rank_note):
    r = rg.get((iso, role, alt)); s = sp.get((iso, role, alt))
    if r is None:
        return
    a = assign.get(iso, {})
    chosen.append(dict(
        gene_name=r["gene_name"], isoform_id=iso, partner=r["partner"],
        assay=assay, oligo_role=role, alt=alt,
        sequence=r["sequence"], length=r["length"], tm_C=r["tm_C"], gc_pct=r["gc_pct"],
        chrom=r["chrom"], genomic_start=r["genomic_start"], genomic_end=r["genomic_end"],
        in_candidate_specific_seq=r["in_candidate_specific_seq"],
        expected_product_candidate_bp=r["expected_product_candidate_bp"],
        expected_product_partner_bp=r["expected_product_partner_bp"],
        transcripts_bound=(s or {}).get("transcripts_full_match", ""),
        offtarget_loci=(s or {}).get("offtarget_loci", ""),
        offtarget_locus_names=(s or {}).get("offtarget_locus_names", ""),
        primary_event=a.get("primary_event", ""), is_nested=a.get("is_nested", ""),
        pacbio_dependent=a.get("qpcr_addressable", ""),
        selection=rank_note, note=r["note"]))


# ----------------------------------------------------------- pairs
bykey = defaultdict(list)
for p in pair:
    bykey[(p["isoform_id"], p["pair"])].append(p)

pair_report = []
for (iso, pr), cands in bykey.items():
    is_qpcr = pr.startswith("qPCR")
    scored = []
    for c in cands:
        ref_bad = B(c.get("ref_record_unspliced"))
        on_target = B(c["locus_clean"]) or ref_bad
        off = I(c["offtarget_loci"])
        amp_partner = B(c["amplifies_partner"])
        if is_qpcr:
            ok = on_target and off == 0 and not amp_partner
            # prefer: valid, then fewest templates, then primer3 rank
            score = (0 if ok else 1, off, I(c["amplifiable_templates"]), int(c["alt"]))
        else:
            ok = on_target and off == 0
            score = (0 if ok else 1, off, 0 if amp_partner else 1, int(c["alt"]))
        scored.append((score, ok, c))
    scored.sort(key=lambda x: x[0])
    best_score, best_ok, best = scored[0]
    role_f, role_r = pr.split("/")
    if best_ok:
        # say what was actually measured for THIS pair, not what the class usually does
        note = "verified: no off-target locus"
        note += "; excludes the partner" if is_qpcr else (
            "; amplifies both isoforms" if B(best["amplifies_partner"])
            else "; amplifies the candidate only")
    else:
        note = f"BEST AVAILABLE - {I(best['offtarget_loci'])} off-target loci"
        if B(best["amplifies_partner"]) and is_qpcr:
            note += "; also amplifies the partner, so it does not discriminate"
    take(iso, role_f, best["alt"], best["assay"], note)
    take(iso, role_r, best["alt"], best["assay"], note)
    pair_report.append(dict(gene_name=best["gene_name"], isoform_id=iso, pair=pr,
                            alternatives_tested=len(cands), chosen_alt=best["alt"],
                            passes=best_ok, offtarget_loci=best["offtarget_loci"],
                            amplifies_partner=best["amplifies_partner"],
                            ref_record_unspliced=best.get("ref_record_unspliced", "")))

# ----------------------------------------------------------- single oligos (RACE)
single = defaultdict(list)
for s in spec:
    if s["oligo_role"].startswith("RACE"):
        single[(s["isoform_id"], s["oligo_role"])].append(s)
for (iso, role), cands in single.items():
    scored = sorted(cands, key=lambda c: (I(c["offtarget_loci"]),
                                          0 if B(c["hits_own_transcript"]) else 1,
                                          I(c["transcripts_full_match"]),
                                          int(c.get("alt", 0))))
    best = scored[0]
    ok = I(best["offtarget_loci"]) == 0
    take(iso, role, best.get("alt", "0"), best["assay"],
         "verified: no off-target locus" if ok
         else f"BEST AVAILABLE - {best['offtarget_loci']} off-target loci")

with open("results/tables/assay_order_sheet.csv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(chosen[0])); w.writeheader()
    w.writerows(sorted(chosen, key=lambda c: (c["gene_name"], c["oligo_role"])))
with open("results/tables/pair_selection_report.csv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(pair_report[0])); w.writeheader(); w.writerows(pair_report)

print(f"oligos on the order sheet : {len(chosen)}")
print(f"candidates covered        : {len(set(c['isoform_id'] for c in chosen))}")
clean = sum(1 for c in chosen if str(c["offtarget_loci"]) == "0")
print(f"with zero off-target loci : {clean} ({100*clean/len(chosen):.0f}%)")
print(f"pairs passing in full     : {sum(p['passes'] for p in pair_report)} of {len(pair_report)}")
from collections import Counter
print("by role:", dict(Counter(c["oligo_role"] for c in chosen)))
