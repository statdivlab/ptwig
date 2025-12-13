#' Run stable search
#'
#' @param newicks Character vector of Newick strings (optional)
#' @param file Path to a file containing Newick trees (optional)
#' @param alpha Numeric value passed to stableSearchRcpp
#' @param summarized Boolean factor indicating if the function is to be runned with a summarized version of the sample
#'
#' @return Output of stableSearchRcpp
#' @export
#'
#' @importFrom ape read.tree
run_stable_search <- function(newicks = NULL, file = NULL, alpha, summarized = FALSE) {
  
  ## --- Argument validation -------------------------------------------------
  if (is.null(newicks) && is.null(file)) {
    stop("You must provide either `newicks` or `file`.")
  }
  if (!is.null(newicks) && !is.null(file)) {
    stop("Provide ONLY one of `newicks` or `file`, not both.")
  }
  
  ## --- Read trees from file (if provided) ---------------------------------
  if (!is.null(file)) {
    if (!file.exists(file)) stop("File does not exist: ", file)
    trees <- ape::read.tree(file)
  } else {
    ## Parse Newick strings directly
    trees <- lapply(newicks, ape::read.tree, text = TRUE)
  }
  
  ## --- Normalize: unroot and remove edge lengths ---------------------------
  trees <- lapply(trees, function(tr) {
    # Unroot the tree
    if (!is.null(tr$root.edge) || !ape::is.rooted(tr)) {
      tr <- ape::unroot(tr)
    }
    
    # Remove branch lengths (set to NULL so ape::write.tree does not print them)
    tr$edge.length <- NULL
    tr
  })
  
  ## --- Compute union of all tip labels -------------------------------------
  completeLeaveSet <- unique(unlist(lapply(trees, function(x) x$tip.label)))
  
 ## --- Preparing to run the final function -----------------------------------
  if(summarized){
    ## --- Joining trees that are identical if summarized = TRUE --------------
    Unique_trees = list()
    Count_trees = c()
    
    for (tree in trees){
      Found = FALSE;
      i = 0;
      for (Top in Unique_trees) {
        i = i+1
        if (all.equal(tree, Top)){
          Found = TRUE;
          Count_trees[i] = Count_trees[i] + 1;
        }
      }
      if (!Found){
        Unique_trees[[length(Unique_trees)+1]] = tree;
        Count_trees = c(Count_trees, 1);
      }
    }
    
    ## --- Write cleaned trees back to Newick strings ---------------------------
    cleaned_newicks <- vapply(Unique_trees, ape::write.tree, FUN.VALUE = character(1))
    
    ## --- Call your Rcpp backend ----------------------------------------------
    res <- stableSearchRcppS(treeSampleR = cleaned_newicks, nSampleR = Count_trees, compLeafSetR = completeLeaveSet, alphaR = alpha)
    
  } else {
    ## --- Write cleaned trees back to Newick strings ---------------------------
    cleaned_newicks <- vapply(trees, ape::write.tree, FUN.VALUE = character(1))
    
    ## --- Call your Rcpp backend ----------------------------------------------
    res <- stableSearchRcpp(treeSampleR = cleaned_newicks,compLeafSetR = completeLeaveSet, alphaR = alpha)
  }
  
  return(res)
}