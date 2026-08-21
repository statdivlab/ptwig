#!/usr/bin/env python3
import sys, re, os
from Bio import SeqIO

# Usage: python 04_extract_gene_alignments.py s97.msa.fa work/consensus.partitions.txt work/gene_alignments_all
if len(sys.argv) != 4:
    sys.exit(f"Usage: {sys.argv[0]} SUPERMATRIX.fasta PARTITIONS.txt OUT_DIR")

super_fa, part_file, out_dir = sys.argv[1], sys.argv[2], sys.argv[3]
os.makedirs(out_dir, exist_ok=True)

parts = []
with open(part_file) as f:
    for line in f:
        line = line.strip()
        if not line or line.startswith("#"): 
            continue
        m = re.match(r'^(\S+)\s*=\s*(\d+)\s*-\s*(\d+)\s*$', line)
        if not m:
            sys.stderr.write(f"Skip malformed partition line: {line}\n")
            continue
        gene, s, e = m.group(1), int(m.group(2)), int(m.group(3))
        if s > e: s, e = e, s
        parts.append((gene, s, e))

records = list(SeqIO.parse(super_fa, "fasta"))
for gene, s, e in parts:
    out_path = os.path.join(out_dir, f"{gene}.fasta")
    with open(out_path, "w") as out:
        for rec in records:
            subseq = str(rec.seq)[s-1:e]  # 1-based inclusive slice
            out.write(f">{rec.id}\n{subseq}\n")
    print(f"Wrote {out_path}")

