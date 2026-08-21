#!/usr/bin/env python3
import sys, os, glob, statistics
from collections import defaultdict

# Usage: python 03_call_gene_blocks.py work/hmmsearch work/maps work/consensus [min_taxa_frac=0.4]
if len(sys.argv) < 4:
    sys.exit(f"Usage: {sys.argv[0]} HMMSEARCH_DIR MAPS_DIR OUT_PREFIX [min_taxa_frac=0.4]")

hmm_dir = sys.argv[1]
maps_dir = sys.argv[2]
out_prefix = sys.argv[3]
min_taxa_frac = float(sys.argv[4]) if len(sys.argv) > 4 else 0.4

# Load maps: taxon -> 1-based array mapping ungapped_pos to supermatrix_col
maps = {}
for tsv in glob.glob(os.path.join(maps_dir, "*.tsv")):
    tax = os.path.basename(tsv).rsplit(".", 1)[0]
    arr = [None]  # 1-based
    with open(tsv) as f:
        next(f)  # header
        for line in f:
            up, col = line.strip().split("\t")
            arr.append(int(col))
    maps[tax] = arr

def best_domain_per_taxon(domtbl_path):
    best = {}
    with open(domtbl_path) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            cols = line.strip().split()
            target = cols[0]           # sequence id (taxon)
            evalue = float(cols[21])   # i-Evalue (domain)
            ali_start = int(cols[17])  # target start (ungapped)
            ali_end   = int(cols[18])  # target end   (ungapped)
            cur = best.get(target)
            if (cur is None) or (evalue < cur[2]):
                best[target] = (ali_start, ali_end, evalue)
    return best

gene_to_cols = {}
for domtbl in glob.glob(os.path.join(hmm_dir, "*.domtbl")):
    gene = os.path.basename(domtbl).replace(".domtbl", "")
    per_tax = best_domain_per_taxon(domtbl)
    spans = []
    for tax, (us, ue, ev) in per_tax.items():
        if tax not in maps: 
            continue
        mp = maps[tax]
        if us >= len(mp) or ue >= len(mp):
            continue
        cs, ce = mp[us], mp[ue]
        if cs is None or ce is None:
            continue
        if cs > ce: cs, ce = ce, cs
        spans.append((cs, ce))
    if not spans:
        continue
    starts = [s for s,_ in spans]
    ends   = [e for _,e in spans]
    med_s, med_e = int(statistics.median(starts)), int(statistics.median(ends))
    if med_s > med_e: med_s, med_e = med_e, med_s

    # keep columns supported by >= min_taxa_frac of taxa
    support = defaultdict(int)
    for s,e in spans:
        for c in range(s, e+1):
            support[c] += 1
    n = len(spans)
    kept = sorted([c for c,v in support.items() if v / n >= min_taxa_frac])
    if not kept:
        continue
    gene_to_cols[gene] = (kept[0], kept[-1])

# Resolve overlaps by sorted start, keep monotonic blocks
ordered = sorted(gene_to_cols.items(), key=lambda x: x[1][0])
nonoverlap = []
for g,(s,e) in ordered:
    if not nonoverlap:
        nonoverlap.append([g,s,e]); continue
    pg, ps, pe = nonoverlap[-1]
    if s <= pe:
        # overlap; pick the larger interval and/or trim current
        if e <= pe:
            prev_span = pe-ps
            cur_span  = e-s
            if cur_span > prev_span:
                nonoverlap[-1] = [g,s,e]
        else:
            s2 = pe+1
            if s2 < e:
                nonoverlap.append([g,s2,e])
    else:
        nonoverlap.append([g,s,e])

part_path = f"{out_prefix}.partitions.txt"
with open(part_path, "w") as out:
    for g,s,e in nonoverlap:
        out.write(f"{g} = {s}-{e}\n")

print(f"Wrote partition: {part_path}")

