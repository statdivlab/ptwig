#' Compute ranks
#'
#' @param newicks Character vector of Newick strings (optional)
#' @param file Path to a file containing Newick trees (optional)
#'
#' @return Output of computeRank
#' @export
#'
#' @importFrom ape read.tree
run_compute_rank <- function(newicks = NULL,file = NULL) {

  ## --- Argument validation -------------------------------------------------
  
  # Define the six allowed conditions as logical variables
  cond1 <- !is.null(newicks) # (newicks != NULL)
  
  cond2 <- !is.null(file) # (file != NULL)
  
  # Create a vector of the conditions
  conditions <- c(cond1, cond2)
  
  # Count how many conditions are TRUE
  num_true_conditions <- sum(conditions)
  
  # Check if exactly one condition is TRUE
  if (num_true_conditions != 1) {
    # If zero conditions are true
    if (num_true_conditions == 0) {
      stop(
        "You must provide arguments satisfying EXACTLY ONE of the following conditions:\n",
        "1. `newicks` is provided\n",
        "2. `file` is provided\n"
      )
    } 
    # If more than one condition is true
    else {
      stop(
        "You have provided both newicks and file. You must provide EXACTLY ONE."
      )
    }
  }
  
  ## --- Read trees from files (if provided) ---------------------------------
  
  if (!is.null(file)) {
    if (!file.exists(file)) stop("File does not exist: ", file)
    trees <- ape::read.tree(file)
  } else {
    ## Parse Newick strings directly
    trees <- lapply(newicks, function(x) ape::read.tree(text = x))
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
  
  cleaned_newicks <- vapply(trees, ape::write.tree, FUN.VALUE = character(1))
  
  ## --- Computing the rank for every tree in the list -----------------------
  res <- lapply(cleaned_newicks, computeRank);
  
  return(res)
}
  
  