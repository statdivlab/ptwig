# ptwig 

**P** hylogenetic **T** rees **W**ith **I**nterpretable **G**uarantees
---

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

### Stable Search (Algorithm 1)

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

**Example: summarized tree sample**

If many trees are repeated in the sample, setting summarized = TRUE can substantially
improve performance:

``` r
res <- run_stable_search(
  file = "trees.nwk",
  alpha = 0.85,
  summarized = TRUE
)
```

### FDR-Controlled Search (Algorithms 1–3)

The function `run_FDRcontrol_search()` internally performs the stable search, subposet
construction, and FDR-controlled inference, returning the trees with FDR control guarantees.
It supports flexible input configurations, including automatic sample splitting.

**Example: single file with specified sample splitting**

The function allows for the user to specify how many trees are used for the stable search and subposet 
construction through the parameter `n1`.

``` r
res <- run_FDRcontrol_search(
  file = "trees.nwk",
  alpha = 0.85,
  q = 0.1,
  tau = 0.95,
  n1 = 60
)
```

**Example: single file with automatic sample splitting**

If `n1` is not provided, it automatically takes half the sample for the stable search and the other half for 
the FDR-controlled inference. 

``` r
res <- run_FDRcontrol_search(
  file = "trees.nwk",
  alpha = 0.85,
  q = 0.1,
  tau = 0.95
)
```

**Random splitting can be enabled via:**

``` r
random_subsampling = TRUE
```

**Example: two independent tree samples**

If separate tree samples are available for stable search and FDR-controlled
testing, they can be provided explicitly:

``` r
res <- run_FDRcontrol_search(
  file1 = "trees_stable.nwk",
  file2 = "trees_test.nwk",
  alpha = 0.85,
  q = 0.1,
  tau = 0.95
)
```
 
**Example: summarized samples for both stages**

``` r
res <- run_FDRcontrol_search(
  file1 = "trees_stable.nwk",
  file2 = "trees_test.nwk",
  alpha = 0.85,
  q = 0.1,
  tau = 0.95,
  summarized = TRUE
)
```
 
## Citation

If you use ptwig in your research, please cite:
> **Consensus Tree Estimation with False Discovery Rate Control via Partially Ordered Sets**
> *arXiv:2511.23433*  
> https://doi.org/10.48550/arXiv.2511.23433
 
## Bug Reports / Feature Requests
If you encounter a bug, unexpected behavior, or would like to request a new feature, please file an issue at:
https://github.com/statdivlab/ptwig/issues

