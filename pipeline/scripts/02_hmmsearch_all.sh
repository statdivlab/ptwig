#!/usr/bin/env bash
set -euo pipefail

# Usage: bash 02_hmmsearch_all.sh work/hmms work/supermatrix.ungapped.fasta work/hmmsearch
HMM_DIR=${1:-work/hmms}
DB=${2:-work/supermatrix.ungapped.fasta}
OUT_DIR=${3:-work/hmmsearch}
mkdir -p "$OUT_DIR"

shopt -s nullglob
for hmm in "$HMM_DIR"/*.hmm; do
  gene=$(basename "$hmm" .hmm)
  hmmsearch --tblout "$OUT_DIR/$gene.tbl" --domtblout "$OUT_DIR/$gene.domtbl" \
            --noali "$hmm" "$DB" > "$OUT_DIR/$gene.search.log"
done
echo "tblout/domtbl written to $OUT_DIR"

