#!/usr/bin/env python3
"""
33_platform_for_shortlist.py -- recompute novel-vs-known for the ACTUAL 44 candidates.

results/tables/candidate_detectability_by_switch.csv was built on the earlier 63-gene
candidate list. It also has 44 rows, which is a coincidence, and ZERO isoform_id overlap
with the fusion-screened shortlist -- so joining on it silently produced empty platform
columns. This script reads class_code / cmp_ref / oId straight out of gff3/H_T.gff3 for the
88 transcripts that make up the 44 switches and applies the project's established rule.

The unit is the SWITCH, not the isoform. A switch needs long reads when either
  (a) either side has a class code other than '=', i.e. its intron chain is not an exact
      match to a reference transcript, or
  (b) both sides map to the SAME cmp_ref -- the non-obvious case, because then no
      isoform-level short-read analysis can separate them at any depth.

Provenance (oId prefix: PB = emitted by the PacBio assembly, ENST = carried over from
Ensembl) is reported separately, because which assembler emitted a model is NOT the same
question as whether its structure is new -- 89 PB models in this annotation carry class
code '=' and match a named Ensembl transcript exactly.

Output: results/tables/shortlist_platform.csv
"""
import csv, os, urllib.parse
from collections import defaultdict

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
os.chdir(ROOT)

assign = list(csv.DictReader(open("results/tables/assay_assignment.csv")))
need = set()
for a in assign:
    need.add(a["isoform_id"]); need.add(a["partner"])
print(f"candidates: {len(assign)}   transcripts needed: {len(need)}")

attr = {}
with open("gff3/H_T.gff3") as fh:
    for line in fh:
        if line.startswith("#"):
            continue
        f = line.split("\t")
        if len(f) < 9 or f[2] != "transcript":
            continue
        d = {}
        for kv in f[8].rstrip().split(";"):
            if "=" in kv:
                k, v = kv.split("=", 1)
                d[k] = urllib.parse.unquote(v)
        tid = d.get("transcript_id") or d.get("ID")
        if tid in need:
            attr[tid] = d
print(f"found in GFF3: {len(attr)} of {len(need)}")

rows = []
for a in assign:
    c = attr.get(a["isoform_id"], {})
    p = attr.get(a["partner"], {})
    cc, pc = c.get("class_code", "?"), p.get("class_code", "?")
    cref, pref = c.get("cmp_ref", ""), p.get("cmp_ref", "")
    same_ref = bool(cref) and cref == pref
    novel_struct = [x for x in (cc, pc) if x != "="]
    pacbio = bool(novel_struct) or same_ref
    if same_ref and not novel_struct:
        reason = (f"both isoforms match the same reference transcript ({cref}), so an "
                  f"isoform-level short-read analysis cannot separate them at any depth")
    elif novel_struct and same_ref:
        reason = (f"class codes {cc}/{pc} are not exact reference matches, and both map to "
                  f"{cref}")
    elif novel_struct:
        reason = (f"class code {cc} (candidate) / {pc} (partner): "
                  f"{'the candidate' if cc != '=' else 'the partner'} is not an exact match "
                  f"to any reference transcript")
    else:
        reason = (f"both isoforms match distinct reference transcripts ({cref} and {pref}), "
                  f"so short-read Illumina against the standard annotation could have found "
                  f"this switch")
    rows.append(dict(
        gene_name=a["gene_name"], isoform_id=a["isoform_id"], partner=a["partner"],
        class_code_candidate=cc, class_code_partner=pc,
        cmp_ref_candidate=cref, cmp_ref_partner=pref,
        same_reference_transcript=same_ref,
        oId_candidate=c.get("oId", ""), oId_partner=p.get("oId", ""),
        assembly_contributed_candidate=c.get("oId", "").startswith("PB"),
        assembly_contributed_partner=p.get("oId", "").startswith("PB"),
        tss_id_candidate=c.get("tss_id", ""), tss_id_partner=p.get("tss_id", ""),
        shared_tss_id=(c.get("tss_id", "x") == p.get("tss_id", "y")),
        platform_required="PacBio (long reads)" if pacbio else "Illumina short reads suffice",
        illumina_findable=not pacbio, reason=reason,
        own_dIF=a["own_dIF"], set=a["set"], primary_assay=a["primary_assay"]))

with open("results/tables/shortlist_platform.csv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)

from collections import Counter
print("\nplatform requirement:", dict(Counter(r["platform_required"] for r in rows)))
print("  of the PacBio-required, via same reference transcript:",
      sum(1 for r in rows if r["same_reference_transcript"]))
print("class codes, candidate side:", dict(Counter(r["class_code_candidate"] for r in rows)))
print("class codes, partner side  :", dict(Counter(r["class_code_partner"] for r in rows)))
print("candidate model emitted by the PacBio assembly:",
      sum(r["assembly_contributed_candidate"] for r in rows), "of", len(rows))
print("\nIllumina-findable candidates:",
      ", ".join(sorted(r["gene_name"] for r in rows if r["illumina_findable"])))
