#' Running a complete analysis of the properties of the subposet constructed from a set of samples
#'
#' @param newicks Character vector of Newick strings for subPoset Building(optional)
#' @param file Path to a file containing Newick trees for subPoset Building (optional)
#' @param Mt Numeric value representing number of maximal trees in subPoset
#' @param rb rank at which bifurcation in subposet starts
#' @param tau Extra value for subposet building.
#'
#' @return Output of computeNewSubposet
#' @export
#'
#' @importFrom ape read.tree
run_NewSubposet <- function(newicks = NULL, file = NULL, Mt, rb) {
  
  ## --- Argument validation -------------------------------------------------
  
  # Define the six allowed conditions as logical variables
  
  cond <- xor((!is.null(newicks) && is.null(file)), (is.null(newicks) && !is.null(file)))
  
  
  # Check if exactly one source for each tree/sample tree

  if (!cond){
    stop(
      "You must provide the Sample for SubPoset construction in EXACTLY ONE of the following arguments:\n",
      "1. as a collection of newick strings in newicks or 2. as a file containing the newick strings in file."
    )
  }
  
 
  ## --- Read trees from files (if provided) ---------------------------------
  
  if (!is.null(newicks)){
    trees <- lapply(newicks, function(x) ape::read.tree(text = x))
  } else {
    if (!file.exists(file)) stop("File does not exist: ", file)
    trees <- ape::read.tree(file)
  }
  
  
  ## --- Normalize: unroot and remove edge lengths ---------------------------
  
  
  trees <- lapply(trees, function(tr) {
    # Unroot the tree
    if (!is.null(tr$root.edge) || !ape::is.rooted(tr)) {
      tr <- ape::unroot(tr)
    }
    
    # Remove branch lengths (set to NULL so ape::write.tree does not print them)
    tr$edge.length <- NULL
    tr$node.label <- NULL
    tr
  })
  
  
  
  ## --- Compute union of all tip labels -------------------------------------
  completeLeaveSet <- unique(unlist(lapply(c(trees), function(x) x$tip.label)))
  
  ## --- Preparing to run the final function -----------------------------------
    ## --- Joining trees that are identical -----------------------------------
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
          break
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
    res <-  computeNewSubposet(treeSampleR = cleaned_newicks, nSampleR = Count_trees,
                               compLeafSetR = completeLeaveSet,
                                MtR = Mt, rbR = rb)
  
  
  return(res)
}