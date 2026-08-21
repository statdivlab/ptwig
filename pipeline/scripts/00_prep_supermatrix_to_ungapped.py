#!/usr/bin/env python3
import sys, os
from Bio import SeqIO

# Usage: python 00_prep_supermatrix_to_ungapped.py s97.msa.fa work
if len(sys.argv) != 3:
    sys.exit(f"Usage: {sys.argv[0]} SUPERMATRIX.fasta OUT_WORK_DIR")

super_fa = sys.argv[1]
work = sys.argv[2]
os.makedirs(work, exist_ok=True)
ungapped_fa = os.path.join(work, "supermatrix.ungapped.fasta")
maps_dir = os.path.join(work, "maps")
os.makedirs(maps_dir, exist_ok=True)

records = list(SeqIO.parse(super_fa, "fasta"))
with open(ungapped_fa, "w") as outfa:
    for rec in records:
        seq = str(rec.seq)
        ungapped = []
        ung2col = []  # 1-based mapping: ungapped_pos -> supermatrix_col
        ungapped_count = 0
        for col_idx, aa in enumerate(seq, start=1):
            if aa not in ('-', '.'):
                ungapped_count += 1
                ungapped.append(aa)
                ung2col.append((ungapped_count, col_idx))
        outfa.write(f">{rec.id}\n{''.join(ungapped)}\n")
        with open(os.path.join(maps_dir, f"{rec.id}.tsv"), "w") as m:
            m.write("ungapped_pos\tsupermatrix_col\n")
            for up, col in ung2col:
                m.write(f"{up}\t{col}\n")

print(f"Wrote: {ungapped_fa}")
print(f"Wrote maps for {len(records)} taxa to: {maps_dir}")

