#!/usr/bin/env python3
"""
29_assay_specificity.py -- find out where every designed oligo ACTUALLY binds.

The external review's lesson was that primer3 scores a primer against the one template
you hand it, so a primer can look perfect and still sit in sequence shared by ten
transcripts. The fix is not a better score, it is a search: take every oligo and look
for its binding site in all 296,613 transcripts of the reference, then ask of each PAIR
how many templates could actually produce an amplicon.

Mispriming is governed by the 3' end -- a primer mismatched at its 5' end still extends,
one mismatched at its 3' terminus usually does not. So the probe is the 3'-TERMINAL 15 nt
of each oligo, written on the strand it anneals to, searched exactly with Aho-Corasick in
a single streaming pass. Every 15-mer hit is then extended to the full oligo length to
separate perfect sites from 3'-anchored partial ones.

The same pass audits the reference itself: each record's length is compared with that
transcript's summed exon length in the GFF3, because 14 of the 97 candidate sequences
turned out to be unspliced genomic intervals rather than cDNA.

Outputs:
  results/tables/assay_specificity.csv        per oligo
  results/tables/assay_pair_specificity.csv   per primer pair: amplifiable templates
  results/tables/reference_fasta_audit.csv    per transcript: spliced or not
"""
import csv, os, sys
from collections import defaultdict
import ahocorasick

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
os.chdir(ROOT)
FASTA = "fa/transcripts_H-T.fa"
PROBE = 15

COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")
rc = lambda s: s.translate(COMP)[::-1]

SENSE_ROLES = {"qPCR_F", "longread_F", "RACE3_GSP1_outer", "RACE3_GSP2_inner"}

oligos = list(csv.DictReader(open("results/tables/assay_reagents.csv")))
print(f"oligos: {len(oligos)}")

# template_seq = the sequence as it appears on the transcript (sense) strand
for i, o in enumerate(oligos):
    o["_id"] = i
    o["_sense"] = o["oligo_role"] in SENSE_ROLES
    o["_template"] = o["sequence"] if o["_sense"] else rc(o["sequence"])
    # 3'-terminal PROBE nt of the primer, expressed on the template strand
    o["_probe"] = o["_template"][-PROBE:] if o["_sense"] else o["_template"][:PROBE]

A = ahocorasick.Automaton()
probe_map = defaultdict(list)
for o in oligos:
    probe_map[o["_probe"]].append(o["_id"])
for k in probe_map:
    A.add_word(k, k)
A.make_automaton()
print(f"distinct probes: {len(probe_map)}")

exon_sum = {}
if os.path.exists("/tmp/all_exon_sums.tsv"):
    for line in open("/tmp/all_exon_sums.tsv"):
        t, s = line.split("\t")
        exon_sum[t] = int(s)
print(f"GFF3 exon sums loaded: {len(exon_sum)}")

by_oligo = defaultdict(list)      # oligo_id -> [(tx, pos, perfect)]
audit = []
n_rec = 0


def scan(tid, seq):
    global n_rec
    n_rec += 1
    if exon_sum:
        es = exon_sum.get(tid)
        audit.append((tid, len(seq), es if es is not None else "",
                      "" if es is None else ("spliced" if len(seq) == es else "NOT_spliced")))
    for end, key in A.iter(seq):
        start = end - len(key) + 1
        for oid in probe_map[key]:
            o = oligos[oid]
            tmpl = o["_template"]
            # extend the 15-mer hit to the full oligo footprint
            if o["_sense"]:
                s0 = end + 1 - len(tmpl)
            else:
                s0 = start
            perfect = 0 <= s0 and s0 + len(tmpl) <= len(seq) and seq[s0:s0 + len(tmpl)] == tmpl
            by_oligo[oid].append((tid, s0, perfect))


cur, buf = None, []
with open(FASTA) as fh:
    for line in fh:
        if line.startswith(">"):
            if cur is not None:
                scan(cur, "".join(buf).upper())
            cur = line[1:].split()[0]
            buf = []
            if n_rec and n_rec % 50000 == 0:
                print(f"  {n_rec} records...", flush=True)
        else:
            buf.append(line.strip())
    if cur is not None:
        scan(cur, "".join(buf).upper())
print(f"records scanned: {n_rec}")

# ------------------------------------------------------------------ per oligo
# Hitting other isoforms OF THE SAME LOCUS is expected and often intended -- a RACE GSP
# is deliberately placed in shared sequence, and a long-read outer pair is meant to
# amplify both isoforms. The number that decides whether an oligo is usable is how many
# OTHER LOCI it binds.
unspliced_ids = {a[0] for a in audit if a[3] == "NOT_spliced"}
print(f"unspliced reference records: {len(unspliced_ids)}")

tx2gene = {}
for line in open("results/tables/tx2gene.tsv"):
    a, b = line.rstrip("\n").split("\t")
    tx2gene[a] = b

rows = []
for o in oligos:
    h = by_oligo.get(o["_id"], [])
    perfect = [x for x in h if x[2]]
    tx_perfect = {x[0] for x in perfect}
    own_locus = tx2gene.get(o["isoform_id"], "?")
    loci = {tx2gene.get(t, "?") for t in tx_perfect}
    rows.append(dict(
        gene_name=o["gene_name"], isoform_id=o["isoform_id"], assay=o["assay"],
        oligo_role=o["oligo_role"], alt=o.get("alt", "0"),
        sequence=o["sequence"], length=o["length"],
        tm_C=o["tm_C"], gc_pct=o["gc_pct"],
        sites_3p_anchored=len(h), sites_full_match=len(perfect),
        transcripts_full_match=len(tx_perfect),
        loci_full_match=len(loci), offtarget_loci=len(loci - {own_locus}),
        offtarget_locus_names=";".join(sorted(loci - {own_locus})[:6]),
        hits_own_transcript=o["isoform_id"] in tx_perfect,
        locus_clean=(loci <= {own_locus}),
        unique_in_transcriptome=(len(tx_perfect) == 1 and o["isoform_id"] in tx_perfect)))
with open("results/tables/assay_specificity.csv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)

# ------------------------------------------------------------------ per pair
PAIRS = [("qPCR_F", "qPCR_R", 70, 200),
         ("longread_F", "longread_R", 200, 12000)]
pair_rows = []
byiso = defaultdict(dict)
for o in oligos:
    byiso[o["isoform_id"]][(o["oligo_role"], o.get("alt", "0"))] = o
for iso, d in byiso.items():
    alts = sorted({k[1] for k in d})
    for fr, rr, lo, hi in PAIRS:
      for a in alts:
        if (fr, a) not in d or (rr, a) not in d:
            continue
        F, R = d[(fr, a)], d[(rr, a)]
        fh_ = defaultdict(list); rh_ = defaultdict(list)
        for tx, pos, perf in by_oligo.get(F["_id"], []):
            fh_[tx].append((pos, perf))
        for tx, pos, perf in by_oligo.get(R["_id"], []):
            rh_[tx].append((pos, perf))
        amp_all, amp_perfect, amp_tx = 0, 0, set()
        for tx in set(fh_) & set(rh_):
            for fp, fperf in fh_[tx]:
                for rp, rperf in rh_[tx]:
                    size = rp + len(R["_template"]) - fp
                    if lo <= size <= hi and rp > fp:
                        amp_all += 1
                        amp_tx.add(tx)
                        if fperf and rperf:
                            amp_perfect += 1
        own_locus = tx2gene.get(iso, "?")
        amp_loci = {tx2gene.get(t, "?") for t in amp_tx}
        partner = F["partner"]
        pair_rows.append(dict(
            gene_name=F["gene_name"], isoform_id=iso, partner=partner,
            assay=F["assay"], pair=fr + "/" + rr, alt=a,
            # If the reference record for this transcript is one of the 36,115 unspliced
            # ones, a junction-spanning pair CANNOT be found in it. An on-target miss here
            # is a defect of the reference, not of the primer -- flagged, not scored.
            ref_record_unspliced=(iso in unspliced_ids),
            partner_record_unspliced=(partner in unspliced_ids),
            forward=F["sequence"], reverse=R["sequence"],
            product_on_candidate_bp=F["expected_product_candidate_bp"],
            product_on_partner_bp=F["expected_product_partner_bp"],
            amplifiable_templates=len(amp_tx),
            amplifiable_loci=len(amp_loci),
            offtarget_loci=len(amp_loci - {own_locus}),
            offtarget_locus_names=";".join(sorted(amp_loci - {own_locus})[:6]),
            amplifies_partner=(partner in amp_tx),
            amplicons_3p_anchored=amp_all, amplicons_full_match=amp_perfect,
            # usable = amplifies the candidate, nothing outside the locus. Whether it
            # should also amplify the partner depends on the assay: a discriminating
            # qPCR must NOT; a long-read outer pair is supposed to.
            locus_clean=(amp_loci <= {own_locus}) and iso in amp_tx,
            specific=(len(amp_tx) == 1 and iso in amp_tx)))
with open("results/tables/assay_pair_specificity.csv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(pair_rows[0])); w.writeheader(); w.writerows(pair_rows)

if audit:
    with open("results/tables/reference_fasta_audit.csv", "w", newline="") as fh:
        w = csv.writer(fh); w.writerow(["transcript_id", "fasta_len", "gff3_exon_sum", "status"])
        w.writerows(audit)
    ns = sum(1 for a in audit if a[3] == "NOT_spliced")
    sp = sum(1 for a in audit if a[3] == "spliced")
    print(f"\nreference audit: {sp} spliced, {ns} NOT spliced, {len(audit)-sp-ns} not in GFF3")

print(f"\noligos: {len(rows)}")
print(f"  bind only within their own locus : {sum(r['locus_clean'] for r in rows)}")
print(f"  unique transcript site           : {sum(r['unique_in_transcriptome'] for r in rows)}")
print(f"  median off-target loci           : "
      f"{sorted(r['offtarget_loci'] for r in rows)[len(rows)//2]}")
print(f"\npairs: {len(pair_rows)}")
print(f"  amplify nothing outside the locus: {sum(p['locus_clean'] for p in pair_rows)}")
print(f"  amplify exactly one template     : {sum(p['specific'] for p in pair_rows)}")
qp = [p for p in pair_rows if p["pair"].startswith("qPCR")]
print(f"  qPCR pairs that are locus-clean AND do not amplify the partner: "
      f"{sum(p['locus_clean'] and not p['amplifies_partner'] for p in qp)} of {len(qp)}")
