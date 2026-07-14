#' Compute stability of every leaf and split of a tree against a tree sample
#'
#' Given a tree `V` and a sample of trees, this computes, for every leaf and
#' every internal split of `V`, its stability: the (weighted) fraction of trees
#' `Z` in the sample for which removing that feature from `V` strictly decreases
#' the similarity rho(V, Z).
#'
#' @param tree_newick Single Newick string for the tree `V` (optional).
#' @param tree_file Path to a file containing a single Newick tree `V` (optional).
#' @param sample_newicks Character vector of Newick strings for the tree sample (optional).
#' @param sample_file Path to a file containing the tree sample (optional).
#' @param summarized Logical; if TRUE, identical trees in the sample are collapsed
#'   and counted (faster for samples with repeated topologies).
#'
#' @return A data.frame with columns `feature` (leaf name as a string, or the
#'   split printed with `{side1|side2}`), `type` (`"leaf"` or `"split"`) and
#'   `stability` (numeric in [0, 1]).
#' @export
#'
#' @importFrom ape read.tree write.tree
run_compute_stability <- function(tree_newick = NULL, tree_file = NULL,
                                  sample_newicks = NULL, sample_file = NULL,
                                  summarized = FALSE) {

  ## --- Argument validation -------------------------------------------------
  if (is.null(tree_newick) && is.null(tree_file)) {
    stop("You must provide either `tree_newick` or `tree_file`.")
  }
  if (!is.null(tree_newick) && !is.null(tree_file)) {
    stop("Provide ONLY one of `tree_newick` or `tree_file`, not both.")
  }
  if (is.null(sample_newicks) && is.null(sample_file)) {
    stop("You must provide either `sample_newicks` or `sample_file`.")
  }
  if (!is.null(sample_newicks) && !is.null(sample_file)) {
    stop("Provide ONLY one of `sample_newicks` or `sample_file`, not both.")
  }

  ## --- Helper: unroot and strip branch lengths / node labels ---------------
  normalize <- function(tr) {
    if (!is.null(tr$root.edge) || !ape::is.rooted(tr)) {
      tr <- ape::unroot(tr)
    }
    tr$edge.length <- NULL
    tr$node.label <- NULL
    tr
  }

  ## --- Read the tree V -----------------------------------------------------
  if (!is.null(tree_file)) {
    if (!file.exists(tree_file)) stop("File does not exist: ", tree_file)
    V <- ape::read.tree(tree_file)
    if (inherits(V, "multiPhylo")) {
      if (length(V) != 1) stop("`tree_file` must contain exactly one tree.")
      V <- V[[1]]
    }
  } else {
    if (length(tree_newick) != 1) stop("`tree_newick` must be a single Newick string.")
    V <- ape::read.tree(text = tree_newick)
  }
  V <- normalize(V)
  tree_clean <- ape::write.tree(V)

  ## --- Read the tree sample ------------------------------------------------
  if (!is.null(sample_file)) {
    if (!file.exists(sample_file)) stop("File does not exist: ", sample_file)
    trees <- ape::read.tree(sample_file)
  } else {
    trees <- lapply(sample_newicks, function(x) ape::read.tree(text = x))
  }
  trees <- lapply(trees, normalize)

  ## --- Compute stability ---------------------------------------------------
  if (summarized) {
    ## --- Collapse identical trees, counting multiplicities -----------------
    Unique_trees <- list()
    Count_trees <- c()

    for (tree in trees) {
      Found <- FALSE
      i <- 0
      for (Top in Unique_trees) {
        i <- i + 1
        if (isTRUE(all.equal(tree, Top))) {
          Found <- TRUE
          Count_trees[i] <- Count_trees[i] + 1
          break
        }
      }
      if (!Found) {
        Unique_trees[[length(Unique_trees) + 1]] <- tree
        Count_trees <- c(Count_trees, 1)
      }
    }

    sample_clean <- vapply(Unique_trees, ape::write.tree, FUN.VALUE = character(1))
    res <- computeStabilityRcppS(treeR = tree_clean,
                                 treeSampleR = sample_clean,
                                 nSampleR = Count_trees)
  } else {
    sample_clean <- vapply(trees, ape::write.tree, FUN.VALUE = character(1))
    res <- computeStabilityRcpp(treeR = tree_clean,
                                treeSampleR = sample_clean)
  }

  return(res)
}
