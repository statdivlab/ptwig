#!/usr/bin/env python
import argparse
from pathlib import Path
from Bio import SeqIO

ALN_EXTS = (".fa", ".fasta", ".faa", ".fas", ".aln", ".fna")

def read_taxa_list(taxa_path, lowercase=False):
    taxa = set()
    with open(taxa_path) as fh:
        for line in fh:
            s = line.strip()
            if not s:
                continue
            taxa.add(s.lower() if lowercase else s)
    return taxa

def header_key(header, strip_after_space=False, lowercase=False):
    key = header.split(None, 1)[0] if strip_after_space else header
    return key.lower() if lowercase else key

def subset_alignment(aln_path, taxa_set, out_dir, strip_after_space, lowercase, min_seqs):
    recs = list(SeqIO.parse(aln_path, "fasta"))
    if not recs:
        return False, 0, 0

    kept = []
    for r in recs:
        key = header_key(r.id, strip_after_space=strip_after_space, lowercase=lowercase)
        if key in taxa_set:
            kept.append(r)

    if len(kept) >= min_seqs:
        out_path = out_dir / aln_path.name
        SeqIO.write(kept, out_path, "fasta")
        return True, len(recs), len(kept)
    else:
        return False, len(recs), len(kept)

def main():
    p = argparse.ArgumentParser(description="Subset each alignment to a specified organism list.")
    p.add_argument("in_dir", help="Directory with input FASTA alignments")
    p.add_argument("taxa_file", help="Text file with taxa IDs, one per line")
    p.add_argument("out_dir", help="Output directory for subset FASTAs")
    p.add_argument("--strip-after-space", action="store_true",
                   help="Match taxa using only the first token of the FASTA header")
    p.add_argument("--lowercase", action="store_true",
                   help="Lowercase both taxa list and headers before matching")
    p.add_argument("--min-seqs", type=int, default=2,
                   help="Write file only if >= this many sequences remain (default: 2)")
    args = p.parse_args()

    in_dir = Path(args.in_dir)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    taxa_set = read_taxa_list(args.taxa_file, lowercase=args.lowercase)

    total, written, skipped = 0, 0, 0
    for path in sorted(in_dir.iterdir()):
        if not path.is_file() or path.suffix.lower() not in ALN_EXTS:
            continue
        total += 1
        ok, n_in, n_out = subset_alignment(
            path, taxa_set, out_dir,
            strip_after_space=args.strip_after_space,
            lowercase=args.lowercase,
            min_seqs=args.min_seqs
        )
        if ok:
            written += 1
            print(f"[write] {path.name}: {n_out}/{n_in} kept")
        else:
            skipped += 1
            print(f"[skip ] {path.name}: {n_out}/{n_in} kept (<{args.min_seqs})")

    print(f"\nDone. Files scanned: {total}, written: {written}, skipped: {skipped}")

if __name__ == "__main__":
    main()
