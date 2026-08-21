#!/usr/bin/env bash
set -euo pipefail

# Usage: bash 01_build_hmms.sh gene_alignments work/hmms
GENE_DIR=${1:-gene_alignments}
OUT_DIR=${2:-work/hmms}
mkdir -p "$OUT_DIR"

shopt -s nullglob
found=0
for aln in "$GENE_DIR"/*.trimmed.faa "$GENE_DIR"/*.trimmed.fa; do
  [ -e "$aln" ] || continue
  found=1
  base=$(basename "$aln")
  gene=${base%%.trimmed.*}      # e.g., arCOG00081
  hmmbuild "$OUT_DIR/$gene.hmm" "$aln" > "$OUT_DIR/$gene.hmmbuild.log"
done

if [[ $found -eq 0 ]]; then
  echo "No *.trimmed.faa (or *.trimmed.fa) files found in $GENE_DIR" >&2
  exit 2
fi

echo "HMMs written to $OUT_DIR"

