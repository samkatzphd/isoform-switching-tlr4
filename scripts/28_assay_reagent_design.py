#!/usr/bin/env python3
"""
28_assay_reagent_design.py -- design the oligos each assigned assay actually needs.

Inputs (all tracked, all produced by earlier steps):
  results/tables/assay_assignment.csv        one row per candidate, with primary_assay
  results/tables/candidate_exons.tsv         transcript_id, chrom, start, end, strand (GFF3)
  results/sequences/ht_candidates_spliced.json   verified spliced cDNA, keyed by transcript id

Output:
  results/tables/assay_reagents.csv          one row per OLIGO
  results/tables/discriminating_regions.csv  one row per candidate: what differs, where

The design is driven by one idea: an oligo is only useful if it can tell the two isoforms
of the switch apart. So every candidate is first reduced to a per-base mask of its own
cDNA saying "this base is also in the partner" / "this base is not". Primers are then
constrained against that mask, and -- the part that matters for interpreting a gel or a
Ct -- every primer is mapped back to genomic coordinates and re-located in the PARTNER
transcript, so the expected product size in BOTH isoforms is reported. A pair that gives
the same size in both is not a discriminating assay no matter what primer3 says about it.
"""
import json, csv, sys, os
from collections import defaultdict
import primer3

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
os.chdir(ROOT)

# ---------------------------------------------------------------- load inputs
exons = defaultdict(list)
for line in open("results/tables/candidate_exons.tsv"):
    t, c, s, e, st = line.rstrip("\n").split("\t")
    exons[t].append((c, int(s), int(e), st))
for t in exons:
    exons[t].sort(key=lambda x: x[1])

seqs = json.load(open("results/sequences/ht_candidates_spliced.json"))
cand = list(csv.DictReader(open("results/tables/assay_assignment.csv")))
print(f"candidates: {len(cand)}   sequences: {len(seqs)}")

COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")
rc = lambda s: s.translate(COMP)[::-1]


def tx_positions(t):
    """Genomic position of every base, in transcript (5'->3') order."""
    blocks = exons[t]
    strand = blocks[0][3]
    pos = []
    for _, a, b, _ in blocks:
        pos.extend(range(a, b + 1))
    if strand == "-":
        pos.reverse()
    return pos, strand


def runs(mask, want):
    """Maximal [start, end) runs where mask == want."""
    out, i, n = [], 0, len(mask)
    while i < n:
        if mask[i] == want:
            j = i
            while j < n and mask[j] == want:
                j += 1
            out.append((i, j))
            i = j
        else:
            i += 1
    return out


# ---------------------------------------------------------------- primer3 settings
QPCR = {
    "PRIMER_OPT_SIZE": 21, "PRIMER_MIN_SIZE": 18, "PRIMER_MAX_SIZE": 25,
    "PRIMER_OPT_TM": 60.0, "PRIMER_MIN_TM": 59.0, "PRIMER_MAX_TM": 62.0,
    "PRIMER_MIN_GC": 40.0, "PRIMER_MAX_GC": 60.0,
    "PRIMER_PRODUCT_SIZE_RANGE": [[70, 200]],
    "PRIMER_MAX_POLY_X": 4, "PRIMER_MAX_SELF_ANY_TH": 45.0,
    "PRIMER_MAX_SELF_END_TH": 35.0, "PRIMER_PAIR_MAX_COMPL_ANY_TH": 45.0,
    "PRIMER_SALT_MONOVALENT": 50.0, "PRIMER_DNA_CONC": 250.0,
    "PRIMER_NUM_RETURN": 5,
}
# gene-specific primers for RACE run hot and long: they must survive a high-stringency
# touchdown against a universal adapter primer.
RACE = {
    "PRIMER_OPT_SIZE": 25, "PRIMER_MIN_SIZE": 22, "PRIMER_MAX_SIZE": 30,
    "PRIMER_OPT_TM": 67.0, "PRIMER_MIN_TM": 64.0, "PRIMER_MAX_TM": 71.0,
    "PRIMER_MIN_GC": 40.0, "PRIMER_MAX_GC": 65.0,
    "PRIMER_MAX_POLY_X": 4, "PRIMER_MAX_SELF_ANY_TH": 45.0,
    "PRIMER_MAX_SELF_END_TH": 35.0,
    "PRIMER_SALT_MONOVALENT": 50.0, "PRIMER_DNA_CONC": 250.0,
    "PRIMER_NUM_RETURN": 5,
}
# outer pair for targeted long-read cDNA: long amplicon, only needs to be specific
LONGREAD = dict(QPCR)
LONGREAD.update({
    "PRIMER_PRODUCT_SIZE_RANGE": [[500, 10000]],
    "PRIMER_OPT_TM": 62.0, "PRIMER_MIN_TM": 60.0, "PRIMER_MAX_TM": 66.0,
    "PRIMER_OPT_SIZE": 23, "PRIMER_MAX_SIZE": 28,
})


def pick_pair(seq, ok_region, settings, name):
    sa = {"SEQUENCE_ID": name, "SEQUENCE_TEMPLATE": seq}
    if ok_region:
        sa["SEQUENCE_PRIMER_PAIR_OK_REGION_LIST"] = ok_region
    try:
        r = primer3.bindings.design_primers(sa, settings)
    except (OSError, ValueError) as e:
        return None, f"primer3 error: {e}"
    if r.get("PRIMER_PAIR_NUM_RETURNED", 0) < 1:
        return None, r.get("PRIMER_LEFT_EXPLAIN", "") + " | " + r.get("PRIMER_RIGHT_EXPLAIN", "")
    return r, None


def pick_right_only(seq, lo, hi, settings, name):
    """Antisense gene-specific primer inside transcript window [lo, hi)."""
    lo, hi = max(0, lo), min(len(seq), hi)
    if hi - lo < 60:
        return None, f"window too small ({hi-lo} nt)"
    s = dict(settings)
    s.update({"PRIMER_PICK_LEFT_PRIMER": 0, "PRIMER_PICK_RIGHT_PRIMER": 1,
              "PRIMER_PICK_INTERNAL_OLIGO": 0, "PRIMER_TASK": "generic"})
    sa = {"SEQUENCE_ID": name, "SEQUENCE_TEMPLATE": seq,
          "SEQUENCE_INCLUDED_REGION": [lo, hi - lo]}
    try:
        r = primer3.bindings.design_primers(sa, s)
    except (OSError, ValueError) as e:
        return None, f"primer3 error: {e}"
    if r.get("PRIMER_RIGHT_NUM_RETURNED", 0) < 1:
        return None, r.get("PRIMER_RIGHT_EXPLAIN", "no right primer")
    return r, None


def pick_left_only(seq, lo, hi, settings, name):
    """Sense gene-specific primer inside transcript window [lo, hi) -- for 3' RACE."""
    lo, hi = max(0, lo), min(len(seq), hi)
    if hi - lo < 60:
        return None, f"window too small ({hi-lo} nt)"
    s = dict(settings)
    s.update({"PRIMER_PICK_LEFT_PRIMER": 1, "PRIMER_PICK_RIGHT_PRIMER": 0,
              "PRIMER_PICK_INTERNAL_OLIGO": 0, "PRIMER_TASK": "generic"})
    sa = {"SEQUENCE_ID": name, "SEQUENCE_TEMPLATE": seq,
          "SEQUENCE_INCLUDED_REGION": [lo, hi - lo]}
    try:
        r = primer3.bindings.design_primers(sa, s)
    except (OSError, ValueError) as e:
        return None, f"primer3 error: {e}"
    if r.get("PRIMER_LEFT_NUM_RETURNED", 0) < 1:
        return None, r.get("PRIMER_LEFT_EXPLAIN", "no left primer")
    return r, None


# ---------------------------------------------------------------- main loop
oligos, regions = [], []

for row in cand:
    t, p = row["isoform_id"], row["partner"]
    gene, assay = row["gene_name"], row["primary_assay"]
    if t not in seqs or t not in exons:
        print(f"  SKIP {gene}: no verified sequence"); continue
    cseq = seqs[t]
    cpos, cstrand = tx_positions(t)
    ppos, _ = tx_positions(p) if p in exons else ([], cstrand)
    pindex = {g: i for i, g in enumerate(ppos)}
    mask = [g in pindex for g in cpos]           # True = base also present in partner

    uruns = runs(mask, False)                     # candidate-specific stretches
    sruns = runs(mask, True)                      # stretches shared with partner
    longest_u = max(uruns, key=lambda r: r[1] - r[0]) if uruns else None
    longest_s = max(sruns, key=lambda r: r[1] - r[0]) if sruns else None

    regions.append(dict(
        gene_name=gene, isoform_id=t, partner=p, primary_assay=assay,
        chrom=exons[t][0][0], strand=cstrand,
        candidate_len=len(cseq), partner_len=len(ppos),
        unique_bp=sum(b - a for a, b in uruns), shared_bp=sum(b - a for a, b in sruns),
        n_unique_blocks=len(uruns),
        longest_unique_tx_start=longest_u[0] if longest_u else "",
        longest_unique_tx_end=longest_u[1] if longest_u else "",
        longest_unique_bp=(longest_u[1] - longest_u[0]) if longest_u else 0,
        five_prime_shared_from=(sruns[0][0] if sruns else ""),
        shared_frac_of_candidate=round(sum(b - a for a, b in sruns) / max(1, len(cseq)), 3),
        # A pair that shares almost no exonic sequence is easy to design against but hard
        # to believe as one gene's two isoforms. The fusion screen only removed ZERO
        # overlap; this is the graded version of the same worry.
        low_overlap_flag=(sum(b - a for a, b in sruns) < 200),
        short_model_flag=(len(cseq) < 300),
        is_nested=row["is_nested"], primary_event=row["primary_event"]))

    def binds_partner(i0, i1):
        """Does a primer occupying candidate tx [i0, i1) have a contiguous partner site?"""
        idx = [pindex.get(g) for g in cpos[i0:i1]]
        if any(v is None for v in idx):
            return None
        if idx[-1] - idx[0] != len(idx) - 1:      # junction absent in partner
            return None
        return idx[0]

    def emit(role, seq_str, tm, gc, i0, i1, assay_label, prod_c="", prod_p="", note="", alt=0):
        oligos.append(dict(
            gene_name=gene, isoform_id=t, partner=p, assay=assay_label, oligo_role=role, alt=alt,
            sequence=seq_str, length=len(seq_str), tm_C=round(float(tm), 1),
            gc_pct=round(float(gc), 1),
            tx_start=i0, tx_end=i1, chrom=exons[t][0][0],
            genomic_start=min(cpos[i0:i1]), genomic_end=max(cpos[i0:i1]),
            in_candidate_specific_seq=not all(mask[i0:i1]),
            binds_partner=binds_partner(i0, i1) is not None,
            expected_product_candidate_bp=prod_c,
            expected_product_partner_bp=prod_p, note=note))

    # ---- pairs for anything PCR-like -------------------------------------
    if "qPCR" in assay:
        ok, used = None, ""
        if longest_u and (longest_u[1] - longest_u[0]) >= 90:
            a, b = longest_u
            ok = [[a, b - a, a, b - a]]; used = "both primers in candidate-specific sequence"
        elif longest_u:
            a, b = longest_u
            ok = [[a, b - a, b, len(cseq) - b], [0, a, a, b - a]]
            used = "one primer anchored in candidate-specific sequence"
        res, err = pick_pair(cseq, ok, QPCR, f"{gene}_qpcr")
        if res is None and ok:                    # relax: let primer3 place freely
            res, err = pick_pair(cseq, None, QPCR, f"{gene}_qpcr")
            used = "UNCONSTRAINED - verify discrimination manually" if res else used
        if res is None:
            print(f"  {gene}: no qPCR pair ({err[:70]})")
        else:
            # Emit every pair primer3 returns, not just its favourite. primer3 ranks on
            # thermodynamics against ONE template; which pair is actually usable is
            # decided later, empirically, by 29's transcriptome scan.
            for k in range(res.get("PRIMER_PAIR_NUM_RETURNED", 0)):
                ls, ll = res[f"PRIMER_LEFT_{k}"]; rs, rl = res[f"PRIMER_RIGHT_{k}"]
                prod = res[f"PRIMER_PAIR_{k}_PRODUCT_SIZE"]
                lp, rp = binds_partner(ls, ls + ll), binds_partner(rs - rl + 1, rs + 1)
                pprod = (rp + rl - lp) if (lp is not None and rp is not None and rp > lp) else "no product"
                emit("qPCR_F", res[f"PRIMER_LEFT_{k}_SEQUENCE"], res[f"PRIMER_LEFT_{k}_TM"],
                     res[f"PRIMER_LEFT_{k}_GC_PERCENT"], ls, ls + ll, assay, prod, pprod, used, k)
                emit("qPCR_R", res[f"PRIMER_RIGHT_{k}_SEQUENCE"], res[f"PRIMER_RIGHT_{k}_TM"],
                     res[f"PRIMER_RIGHT_{k}_GC_PERCENT"], rs - rl + 1, rs + 1, assay, prod, pprod, used, k)

    # ---- "total gene" pair, both primers in SHARED sequence ---------------
    # Without this the qPCR cannot answer the question. dIF is a FRACTION, and LPS changes
    # total output at these loci substantially -- RSAD2 and USP18 are interferon-induced.
    # An isoform-specific Ct alone cannot separate "the isoform's share rose" from "the
    # whole gene was induced". The isoform-specific amplicon is the numerator; this is the
    # denominator, and it has to sit in sequence common to both isoforms.
    if longest_s and (longest_s[1] - longest_s[0]) >= 90:
        a, b = longest_s
        res, err = pick_pair(cseq, [[a, b - a, a, b - a]], QPCR, f"{gene}_total")
        if res is None:
            print(f"  {gene}: no total-gene pair ({str(err)[:60]})")
        else:
            for k in range(res.get("PRIMER_PAIR_NUM_RETURNED", 0)):
                ls, ll = res[f"PRIMER_LEFT_{k}"]; rs, rl = res[f"PRIMER_RIGHT_{k}"]
                prod = res[f"PRIMER_PAIR_{k}_PRODUCT_SIZE"]
                lp, rp = binds_partner(ls, ls + ll), binds_partner(rs - rl + 1, rs + 1)
                pprod = (rp + rl - lp) if (lp is not None and rp is not None and rp > lp) else "no product"
                for role, s_, tm_, gc_, i0, i1 in (
                        ("qPCR_total_F", res[f"PRIMER_LEFT_{k}_SEQUENCE"], res[f"PRIMER_LEFT_{k}_TM"],
                         res[f"PRIMER_LEFT_{k}_GC_PERCENT"], ls, ls + ll),
                        ("qPCR_total_R", res[f"PRIMER_RIGHT_{k}_SEQUENCE"], res[f"PRIMER_RIGHT_{k}_TM"],
                         res[f"PRIMER_RIGHT_{k}_GC_PERCENT"], rs - rl + 1, rs + 1)):
                    emit(role, s_, tm_, gc_, i0, i1, "total-gene reference amplicon", prod, pprod,
                         "denominator for the isoform fraction; must amplify BOTH isoforms equally", k)

    # ---- nested gene-specific primers for RACE ---------------------------
    if "5' RACE" in assay:
        # Preference: put the GSP in sequence SHARED with the partner, so one reaction
        # reports both TSSs as two product sizes. If the pair shares too little sequence
        # to hold a primer, fall back to an isoform-SPECIFIC GSP -- still a valid
        # experiment, it just reports presence of one TSS rather than the ratio.
        win, how = None, ""
        for a, b in sruns:
            if b - a >= 140 and b > 70:
                win, how = (max(a, 70), b), "GSP in shared sequence: both TSSs amplify, sized apart"
                break
        if win is None and longest_s and longest_s[1] - longest_s[0] >= 140:
            win, how = (max(longest_s[0], 70), longest_s[1]), "GSP in shared sequence (internal block)"
        if win is None:
            win, how = (70, len(cseq)), "ISOFORM-SPECIFIC GSP: pair shares too little sequence to place a common primer"
        lo, hi = win
        res2, err2 = pick_right_only(cseq, lo, min(hi, lo + 330), RACE, f"{gene}_GSP2")
        if res2 is None:
            res2, err2 = pick_right_only(cseq, lo, hi, RACE, f"{gene}_GSP2")
        if res2 is None:
            print(f"  {gene} RACE5_GSP2: {str(err2)[:60]}")
        else:
            for k in range(res2.get("PRIMER_RIGHT_NUM_RETURNED", 0)):
                ks, kl = res2[f"PRIMER_RIGHT_{k}"]
                j0, j1 = ks - kl + 1, ks + 1
                pk = binds_partner(j0, j1)
                emit("RACE5_GSP2_inner", res2[f"PRIMER_RIGHT_{k}_SEQUENCE"],
                     res2[f"PRIMER_RIGHT_{k}_TM"], res2[f"PRIMER_RIGHT_{k}_GC_PERCENT"],
                     j0, j1, "5' RACE", ks + 1,
                     (pk + kl) if pk is not None else "no site in partner", how, k)
            rs, rl = res2["PRIMER_RIGHT_0"]
            # outer primer must sit 3' of the inner one for the nesting to mean anything
            res1, err1 = pick_right_only(cseq, rs + 15, min(rs + 470, hi), RACE, f"{gene}_GSP1")
            if res1 is None:
                res1, err1 = pick_right_only(cseq, rs + 15, min(rs + 700, len(cseq)), RACE, f"{gene}_GSP1")
            if res1 is None:
                print(f"  {gene} RACE5_GSP1: {str(err1)[:60]}")
            else:
                for k in range(res1.get("PRIMER_RIGHT_NUM_RETURNED", 0)):
                    r1s, r1l = res1[f"PRIMER_RIGHT_{k}"]
                    j0, j1 = r1s - r1l + 1, r1s + 1
                    pp1 = binds_partner(j0, j1)
                    emit("RACE5_GSP1_outer", res1[f"PRIMER_RIGHT_{k}_SEQUENCE"],
                         res1[f"PRIMER_RIGHT_{k}_TM"], res1[f"PRIMER_RIGHT_{k}_GC_PERCENT"],
                         j0, j1, "5' RACE", r1s + 1,
                         (pp1 + r1l) if pp1 is not None else "no site in partner",
                         f"outer primer for first-round PCR; {r1s - rs} nt 3' of the inner GSP", k)

    if "3' RACE" in assay:
        e0 = sruns[-1][0] if sruns else 0
        eend = sruns[-1][1] if sruns else len(cseq)
        for role, lo, hi in (("RACE3_GSP2_inner", max(e0, eend - 330), eend - 70),
                             ("RACE3_GSP1_outer", max(e0, eend - 700), eend - 260)):
            res, err = pick_left_only(cseq, lo, hi, RACE, f"{gene}_{role}")
            if res is None:
                print(f"  {gene} {role}: {err[:60]}"); continue
            ls, ll = res["PRIMER_LEFT_0"]
            pp = binds_partner(ls, ls + ll)
            emit(role, res["PRIMER_LEFT_0_SEQUENCE"], res["PRIMER_LEFT_0_TM"],
                 res["PRIMER_LEFT_0_GC_PERCENT"], ls, ls + ll, "3' RACE",
                 len(cseq) - ls, (len(ppos) - pp) if pp is not None else "no site in partner",
                 "sense GSP in sequence shared with partner; product length reports which 3' end was used")

    # ---- outer flanking pair for targeted long-read cDNA -----------------
    # The only way to see a nested isoform as a distinct molecule: amplify the whole
    # thing and sequence it, rather than asking a short amplicon to be isoform-specific.
    LR = dict(LONGREAD)
    LR["PRIMER_PRODUCT_SIZE_RANGE"] = [[min(200, max(100, len(cseq) - 80)), len(cseq)]]
    flank = min(180, max(40, len(cseq) // 4))
    res, err = pick_pair(cseq, [[0, flank, max(0, len(cseq) - flank), flank]], LR, f"{gene}_lr")
    if res is None:
        res, err = pick_pair(cseq, None, LR, f"{gene}_lr")
    if res is None:
        print(f"  {gene}: no long-read outer pair ({str(err)[:60]})")
    else:
        for k in range(res.get("PRIMER_PAIR_NUM_RETURNED", 0)):
            ls, ll = res[f"PRIMER_LEFT_{k}"]; rs, rl = res[f"PRIMER_RIGHT_{k}"]
            prod = res[f"PRIMER_PAIR_{k}_PRODUCT_SIZE"]
            lp, rp = binds_partner(ls, ls + ll), binds_partner(rs - rl + 1, rs + 1)
            pprod = (rp + rl - lp) if (lp is not None and rp is not None and rp > lp) else "no product"
            for role, s_, tm_, gc_, i0, i1 in (
                    ("longread_F", res[f"PRIMER_LEFT_{k}_SEQUENCE"], res[f"PRIMER_LEFT_{k}_TM"],
                     res[f"PRIMER_LEFT_{k}_GC_PERCENT"], ls, ls + ll),
                    ("longread_R", res[f"PRIMER_RIGHT_{k}_SEQUENCE"], res[f"PRIMER_RIGHT_{k}_TM"],
                     res[f"PRIMER_RIGHT_{k}_GC_PERCENT"], rs - rl + 1, rs + 1)):
                emit(role, s_, tm_, gc_, i0, i1, "targeted long-read cDNA", prod, pprod,
                     "outer pair; sequence the amplicon rather than sizing it", k)

# ---------------------------------------------------------------- write
os.makedirs("results/tables", exist_ok=True)
with open("results/tables/assay_reagents.csv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(oligos[0])); w.writeheader(); w.writerows(oligos)
with open("results/tables/discriminating_regions.csv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(regions[0])); w.writeheader(); w.writerows(regions)

print(f"\noligos designed : {len(oligos)}")
print(f"candidates covered: {len(set(o['isoform_id'] for o in oligos))} of {len(cand)}")
from collections import Counter
print("by role:", dict(Counter(o["oligo_role"] for o in oligos)))
print("discriminating (no partner product):",
      sum(1 for o in oligos if o["expected_product_partner_bp"] in ("no product", "no site in partner")))
