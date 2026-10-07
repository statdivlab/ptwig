#!/usr/bin/env python
import argparse
import shutil
import subprocess
from pathlib import Path
from Bio import SeqIO

ALN_EXTS = (".fa", ".fasta", ".faa", ".fas", ".aln", ".fna")
# Treat these as "missing" (adjust if you like)
MISSING = set("-?NnXx")

def which_iqtree():
    for cand in ("iqtree3", "iqtree2", "iqtree"):
        p = shutil.which(cand)
        if p:
            return p, cand
    return None, None

def non_gap_prop(seq: str) -> float:
    if not seq:
        return 0.0
    non_gaps = sum(1 for c in seq if c not in MISSING)
    return non_gaps / len(seq)

def clean_alignment(aln_path: Path, cleaned_dir: Path, min_prop: float):
    """Return path to cleaned alignment, (orig_n, kept_n).
       Drops sequences with non-gap proportion < min_prop."""
    recs = list(SeqIO.parse(aln_path, "fasta"))
    if not recs:
        return None, 0, 0
    kept = [r for r in recs if non_gap_prop(str(r.seq)) >= min_prop]
    if not kept:
        return None, len(recs), 0
    cleaned_dir.mkdir(parents=True, exist_ok=True)
    out_path = cleaned_dir / aln_path.name
    SeqIO.write(kept, out_path, "fasta")
    return out_path, len(recs), len(kept)

def build_cmd(iqbin, aln, outdir, model, boot, alrt, threads, seed, redo, st):
    base = Path(aln).stem
    prefix = str((Path(outdir) / base))
    cmd = [
        iqbin,
        "-s", str(aln),
        "-pre", prefix,
        "-m", model,
        "-B", str(boot),
        "-alrt", str(alrt),
        "-nt", str(threads),
        "-seed", str(seed),
    ]
    if redo:
        cmd.append("-redo")
    if st:
        cmd.extend(["-st", st])  # DNA / AA / BIN
    return cmd

def main():
    ap = argparse.ArgumentParser(
        description="Clean alignments (drop all-gap/mostly-gap taxa) and run IQ-TREE (v3 preferred)."
    )
    ap.add_argument("in_dir", help="Directory with input FASTA alignments")
    ap.add_argument("out_dir", help="Directory to place IQ-TREE outputs (prefix directory)")
    ap.add_argument("--model", default="MFP", help="Substitution model (default: MFP)")
    ap.add_argument("--boot", type=int, default=1000, help="Ultrafast bootstrap replicates")
    ap.add_argument("--alrt", type=int, default=1000, help="SH-aLRT replicates")
    ap.add_argument("--threads", default="AUTO", help="Thread count or AUTO")
    ap.add_argument("--seed", type=int, default=12345, help="Random seed")
    ap.add_argument("--redo", action="store_true", help="Overwrite existing results")
    ap.add_argument("--st", default=None, help="Force data type: DNA/AA/BIN (optional)")
    ap.add_argument("--min-seqs", type=int, default=2, help="Skip if fewer sequences remain after cleaning")
    ap.add_argument("--min-non-gap-prop", type=float, default=0.05,
                    help="Minimum non-gap proportion to keep a sequence (default 0.05)")
    ap.add_argument("--cleaned-dir", default=None,
                    help="Directory for cleaned FASTAs (default: <out_dir>/cleaned)")
    args = ap.parse_args()

    iqbin, binname = which_iqtree()
    if not iqbin:
        raise SystemExit("ERROR: iqtree3/iqtree2/iqtree not found in PATH.")

    print(f"Using {binname}: {iqbin}")

    in_dir = Path(args.in_dir)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    cleaned_dir = Path(args.cleaned_dir) if args.cleaned_dir else (out_dir / "cleaned")
    cleaned_dir.mkdir(parents=True, exist_ok=True)

    processed = skipped = 0
    for aln in sorted(in_dir.iterdir()):
        if not (aln.is_file() and aln.suffix.lower() in ALN_EXTS):
            continue

        cleaned, orig_n, kept_n = clean_alignment(aln, cleaned_dir, args.min_non_gap_prop)
        if cleaned is None or kept_n < args.min_seqs:
            print(f"[skip ] {aln.name}: {kept_n}/{orig_n} usable seqs (<{args.min_seqs})")
            skipped += 1
            continue

        cmd = build_cmd(
            iqbin=iqbin,
            aln=cleaned,                 # <- run on the CLEANED file
            outdir=out_dir,
            model=args.model,
            boot=args.boot,
            alrt=args.alrt,
            threads=args.threads,
            seed=args.seed,
            redo=args.redo,
            st=args.st
        )

        print(f"[run  ] {aln.name} → {' '.join(cmd)}")
        subprocess.run(cmd, check=True)
        processed += 1

    print(f"\nDone. Processed: {processed}, skipped: {skipped}")

if __name__ == "__main__":
    main()
