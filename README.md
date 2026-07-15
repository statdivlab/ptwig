# ptwig


**P**hylogenetic **T**rees **W**ith **I**nterpretable **G**uarantees

`ptwig` is an `R` package for **stable and FDR-controlled inference on collections of phylogenetic trees**.  
The package implements the algorithms described in

> **Consensus Tree Estimation with False Discovery Rate Control via Partially Ordered Sets**
> *arXiv:2511.23433*  
> https://doi.org/10.48550/arXiv.2511.23433

The package provides efficient implementations of:
- **Algorithm 1**: Stable search for phylogenetic trees  
- **Algorithm 3**: Subposet construction induced by stable trees  
- **Algorithm 2**: False discovery rate–controlled inference using sample splitting  

All core algorithms are implemented in **C++ via Rcpp**, with high-level R wrappers.

---

## Installation

To install the development version of `ptwig` from GitHub:

``` r
# install.packages("devtools")
devtools::install_github("statdivlab/ptwig")
library(ptwig)
```

## Use

### Stable Search (---)

The function `run_stable_search()` identifies stable trees from a collection 
of phylogenetic trees. Trees are provided via newick-formatted strings. 
The input may be given either as an array of strings or via a file (with a tree per line).


**Example: using a file of Newick trees**

``` r
res <- run_stable_search(
  file = "trees.nwk",
  alpha = 0.85
)
```

### FDR-Controlled Search (Algorithm 1)

`run_FDRcontrol_search()` runs the full pipeline — subposet construction followed by
FDR-controlled inference on a held-out portion of the sample — and returns a single
consensus tree whose splits carry a false-discovery-rate guarantee at level `q`.

The result is a length-1 character vector holding one Newick string. An empty result,
`"();"`, means no split could be reported at the requested FDR level. Parse it back with
`ape::read.tree(text = res)`.

**Choosing the subposet builder (`SPbuilder`)**

How the candidate subposet is built is controlled by `SPbuilder`. The default,
`"basic_score"`, is the recommended general-purpose option:

| `SPbuilder`               | What it does                                                            | Key parameters                                    |
|---------------------------|------------------------------------------------------------------------|---------------------------------------------------|
| `"basic_score"` (default) | Score-based construction, grown `"upwards"` by default                 | `q`, `top_width`, `bottom_width`, `orientation`   |
| `"basic_bifurcation"`     | Fixed bifurcation with a known number of maximal trees and anchor rank | `q`, `Mt`, `rb`                                   |
| `"basic_auto"`            | Bifurcation with automatic `Mt` / `rb` selection                       | `q`, `tau`                                        |
| `"stability"`             | Builds a stable tree first, then branches out (Algorithm 3)            | `q`, `alpha`, `tau`                               |

For the default `"basic_score"` builder the tuning knobs are `q` (the FDR level), the
subposet width (`top_width`, `bottom_width`), and the growth `orientation`
(`"upwards"` or `"downwards"`). `alpha` and `tau` only affect the stability/auto builders
and are ignored here.

**Minimal example (inline Newick strings, all defaults)**

``` r
library(ptwig)

trees <- c(
  "((A,B),(C,D),(E,F));",
  "((A,B),(C,D),(E,F));",
  "((A,B),(C,D),(E,F));",
  "((A,B),(C,E),(D,F));",
  "((A,C),(B,D),(E,F));"
)

res <- run_FDRcontrol_search(newicks = trees, q = 0.1)  # SPbuilder = "basic_score"
res
#> the FDR-controlled consensus tree, as a Newick string

ape::read.tree(text = res)   # parse it back into a phylo object
```

**From a file, with automatic sample splitting**

The pipeline splits the sample in two: one part builds the subposet, the other drives the
FDR-controlled inference. If `n1` is not given, the sample is split in half.

``` r
res <- run_FDRcontrol_search(
  file = "trees.nwk",   # one Newick tree per line
  q    = 0.1
)
```

Set the split size with `n1`, and randomize which trees go into each half with
`random_subsampling = TRUE`:

``` r
res <- run_FDRcontrol_search(
  file = "trees.nwk",
  q    = 0.1,
  n1   = 60,
  random_subsampling = TRUE
)
```

**Two independent tree samples**

If you already have separate samples for construction and testing, pass them directly
(as files or as Newick vectors):

``` r
res <- run_FDRcontrol_search(
  file1 = "trees_construct.nwk",   # builds the subposet
  file2 = "trees_test.nwk",        # FDR-controlled inference
  q     = 0.1
)
```

**Using a different builder**

``` r
# Stability-based construction (Algorithm 3): uses alpha and tau
res <- run_FDRcontrol_search(
  file      = "trees.nwk",
  SPbuilder = "stability",
  q         = 0.1,
  alpha     = 0.85,
  tau       = 0.80
)
```

## Minor / utility functions

A few functions to explore tree measures that the main algorithms are built on.
Two quantities recur throughout `ptwig` (both defined in the paper):

- **rank** — an integer expressing the complexity of a single tree: `0` for a fully unresolved
  star tree, growing as the tree becomes more resolved.
- **ρ (rho)** — the similarity between two trees, on the same scale as rank
  (`0 ≤ ρ(t1, t2) ≤ rank(t1)`): how much structure the two trees share in common.

All of these take Newick input (via `newicks*` character vectors or `file*` paths, one
tree per line) and are vectorized. For the two-tree helpers you supply exactly one of
`newicks1 + newicks2`, `file1 + file2`, `file1 + newicks2`, or `file2 + newicks1`; the two
lists are compared element-by-element, except that a length-1 second list is treated as a
single shared target compared against every tree in the first list.

### `run_compute_rank()` — tree resolution

Returns the rank of each input tree.

``` r
run_compute_rank(newicks = c("((A,B),(C,D),(E,F));", "(A,B,C,D,E,F);"))
#> 5 0        # a resolved 6-tip tree vs. an unresolved star
```

### `run_compute_similarity()` — ρ similarity

Returns `ρ(tree1, tree2)` for each pair.

``` r
run_compute_similarity(
  newicks1 = "((A,B),(C,D),(E,F));",
  newicks2 = "((A,B),(C,E),(D,F));"
)
#> 3
```

### `run_compute_FD()` — false discoveries

Returns the number of False Discoveries incurred by the first tree in comparison to the second (the
"target tree"): `rank(tree1) - ρ(tree1, tree2)`.

``` r
run_compute_FD(
  newicks1 = "((A,B),(C,D),(E,F));",
  newicks2 = "((A,B),(C,E),(D,F));"
)
#> 2          # rank 5 minus similarity 3
```

### `run_compute_FDP()` — false discovery proportion

Returns the false discovery proportion `FD / rank(tree1)`, a value in `[0, 1]`
(`0` when the first tree has rank 0).

``` r
run_compute_FDP(
  newicks1 = "((A,B),(C,D),(E,F));",
  newicks2 = "((A,B),(C,E),(D,F));"
)
#> 0.4        # 2 / 5
```

## Citation

If you use ptwig in your research, please cite:
> **Consensus Tree Estimation with False Discovery Rate Control via Partially Ordered Sets**
> *arXiv:2511.23433*  
> https://doi.org/10.48550/arXiv.2511.23433
 
## Bug Reports / Feature Requests
If you encounter a bug, unexpected behavior, or would like to request a new feature, please file an issue at:
https://github.com/statdivlab/ptwig/issues

