As detailed in the preprint, Zhang et al made their multiple sequence alignment available, but not individual alignments nor individual gene trees that include the Eukarya. There are, however, alignments and gene trees for the Archaea. Therefore, we used the archaeal alignments to train HMMs to detect these genes in the eukaryotic genomes, then re-aligned and re-estimated gene trees to obtain trees that include both archaeal and eukaryotic organisms. 

This directory contains scripts to estimate those trees. There are a number of steps. We use pixi to organize the `scripts` using only the `input` data. You will need to install pixi to reproduce this. (And you should! pixi is awesome and a billion times better than conda!!)

We have copied the output of this pipeline (the gene trees) to `../data/`, in keeping with typical practices for R packages. 

We have done our best to make this analysis reproducible, but please let us know if there is anything missing. Many thanks to Zhang et al for making their data available. 

# Overview

This pipeline takes a concatenated protein supermatrix and a set of seed
single-gene alignments, recovers which columns of the supermatrix belong to
which gene, and estimates a maximum-likelihood tree for each gene on a chosen
subset of taxa. The resulting gene trees are the input to the `ptwig` analyses
in the parent repository.

## Layout

```
pipeline/
├── input/                        given data (tracked)
│   ├── s97.msa.fa                supermatrix: 425 taxa x 20067 columns
│   ├── gene_alignments/          97 seed alignments, *.trimmed.faa
│   └── chosen_orgs.txt           15 taxa to retain for tree estimation
├── scripts/                      pipeline steps 00-06
├── work/                         intermediates (generated, not tracked)
├── output/final_trees/           the deliverable gene trees (generated)
├── pixi.toml                     environment + task definitions
└── pixi.lock                     pinned dependency versions
```

## Requirements

All software is pinned in `pixi.lock`; nothing else needs to be installed.

```sh
cd pipeline
pixi install
```

This provides Python 3.10-3.12, Biopython, HMMER >= 3.3, and IQ-TREE >= 3.0.
Every command below must be run from the `pipeline/` directory.

## Running

The whole pipeline:

```sh
pixi run all
```

Or one step at a time:

| Task | Script | Reads | Writes |
|---|---|---|---|
| `prep` | `00_prep_supermatrix_to_ungapped.py` | `input/s97.msa.fa` | `work/supermatrix.ungapped.fasta`, `work/maps/<taxon>.tsv` |
| `build-hmms` | `01_build_hmms.sh` | `input/gene_alignments/*.trimmed.faa` | `work/hmms/<gene>.hmm` |
| `hmmsearch-all` | `02_hmmsearch_all.sh` | `work/hmms`, ungapped supermatrix | `work/hmmsearch/<gene>.{tbl,domtbl}` |
| `call-blocks` | `03_call_gene_blocks.py` | `work/hmmsearch`, `work/maps` | `work/consensus.partitions.txt` |
| `extract` | `04_extract_gene_alignments.py` | `input/s97.msa.fa`, partitions | `work/gene_alignments_all/<gene>.fasta` |
| `subset` | `05_subset_alignments.py` | `work/gene_alignments_all`, `input/chosen_orgs.txt` | `work/gene_alignments_subset/` |
| `iqtree-batch` | `06_run_iqtree.py` | `work/gene_alignments_subset` | `work/trees/` |
| `collect-trees` | — | `work/trees/*.treefile` | `output/final_trees/` |

Two grouped tasks are also defined: `alignments` runs steps 0-4 (supermatrix to
per-gene alignments), and `all` runs steps 0-7. `pixi run clean` deletes `work/`
and `output/` without touching `input/` or `scripts/`.

## How it works

The supermatrix has already had its gene boundaries erased by concatenation, so
they are recovered by search rather than by bookkeeping.

1. **Prep.** Gaps are stripped from each taxon's row, and a map is written from
   each ungapped residue position back to its supermatrix column. HMMER searches
   ungapped sequence, so this map is what lets a hit be translated back into
   supermatrix coordinates.
2. **Build HMMs.** One profile HMM is built per seed alignment. The seeds cover
   a wider taxon set (538 sequences) than the supermatrix itself; they are used
   only as models of each gene family, not as data.
3. **Search.** Every HMM is searched against every taxon's ungapped sequence.
4. **Call blocks.** For each gene, the best-scoring domain per taxon gives an
   interval in ungapped coordinates, which the maps convert to supermatrix
   columns. Columns supported by at least 40% of the taxa hitting that gene are
   kept, and overlapping gene blocks are resolved into a non-overlapping,
   monotonic partition. Of the 97 seed genes, 85 survive this step.
5. **Extract.** Each partition is cut out of the supermatrix, giving one
   alignment per gene across all 425 taxa.
6. **Subset.** Each gene alignment is reduced to the 15 organisms in
   `chosen_orgs.txt`.
7. **Trees.** Sequences that are more than 95% gaps after subsetting are
   dropped, then IQ-TREE estimates a tree per gene with ModelFinder (`MFP`),
   1000 ultrafast bootstrap replicates, and 1000 SH-aLRT replicates, seeded at
   12345.

## Expected output

85 Newick files in `output/final_trees/`, one per recovered gene, each on the
15 chosen organisms.

## Notes

- Gene names carry the prefix of the family database they came from: 92 arCOG,
  2 Pfam, 2 TIGRFAM, 1 COG.
- Step 0 writes one map file per taxon and step 2 writes full HMMER output, so
  `work/` grows to roughly 140 MB. It is safe to delete after `output/` exists.
- The files in `input/gene_alignments/` came from `2. phylogenome files/2.S97/1.individual markers files/2.remove paralogs/<gene>.files/<gene>.trimmed.faa` in Zhang et al's supplementary data. 
- Similarly, `input/s97.msa.fa` came from `2. phylogenome files/2.S97/2.concatenated files`.
- `input/chosen_orgs.txt` was chosen by me (Amy) to balance having an interesting and well-motivated data analysis with the computational intensity of our consensus tree estimation method.  
- The seeds and the supermatrix do not sample the same taxa. The seeds cover 579
  taxa and the supermatrix 425, with 410 in common. The 169 seed-only taxa are
  non-Asgard archaeal outgroups (Thermoproteia, Thermoplasmata, Nitrososphaeria,
  Korarchaeia, and DPANN lineages). The 15 supermatrix-only taxa are all 14
  eukaryotes plus `DZ1B.157`, so the HMMs are trained on archaeal sequence alone
  and then used to delimit columns in a matrix that also contains eukaryotes.
  This is sound, since a partition is a property of the alignment columns and
  applies to every row regardless of which taxa trained the model, but note that
  5 of the 15 organisms in `chosen_orgs.txt` are eukaryotes unseen by any HMM.
