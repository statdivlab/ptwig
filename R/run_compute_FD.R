#' Run computation of FD 
#'
#' @param newicks1 Character vector of Newick strings for first tree(s) (optional).
#' @param newicks2 Character vector of Newick strings for target trees(s)  (optional).
#' @param file1 Path to a file containing Newick trees for first tree(s) (optional).
#' @param file2 Path to a file containing Newick trees for target trees(s) (optional).
#'
#' @return Output of computeFD
#' @export
#'
#' @importFrom ape read.tree
run_compute_FD <- function(newicks1 = NULL, newicks2 = NULL,
                                  file1 = NULL, file2 = NULL) {
  
  ## --- Argument validation -------------------------------------------------
  
  # Define the six allowed conditions as logical variables
  cond1 <- !is.null(newicks1) && !is.null(newicks2) # (newicks1 !=NULL && newicks2 != NULL)
  
  cond2 <- !is.null(file1) && !is.null(file2) # (file1 != NULL && file2 != NULL)
  
  cond3 <- !is.null(file1) && !is.null(newicks2) # (file1 != NULL && newicks2 != NULL)
  
  cond4 <- !is.null(file2) && !is.null(newicks1) # (file2 != NULL && newicks1 != NULL)
  
  # Create a vector of the conditions
  conditions <- c(cond1, cond2, cond3, cond4)
  
  # Count how many conditions are TRUE
  num_true_conditions <- sum(conditions)
  
  # Check if exactly one condition is TRUE
  if (num_true_conditions != 1) {
    # If zero conditions are true
    if (num_true_conditions == 0) {
      stop(
        "You must provide arguments satisfying EXACTLY ONE of the following conditions:\n",
        "1. `newicks1` AND `newicks2` are provided\n",
        "2. `file1` AND `file2` are provided\n",
        "3. `file1` AND `newicks2` are provided\n",
        "4. `file2` AND `newicks1` are provided"
      )
    } 
    # If more than one condition is true
    else {
      stop(
        "You have provided arguments satisfying MORE THAN ONE of the allowed combinations. ",
        "Please provide arguments satisfying EXACTLY ONE of the following conditions:\n",
        "1. `newicks1` AND `newicks2` are provided\n",
        "2. `file1` AND `file2` are provided\n",
        "3. `file1` AND `newicks2` are provided\n",
        "4. `file2` AND `newicks1` are provided"
      )
    }
  }
  
  ## --- Read trees from files (if provided) ---------------------------------
  if (cond1){
    trees1 <- lapply(newicks1, function(x) ape::read.tree(text = x))
    trees2 <- lapply(newicks2, function(x) ape::read.tree(text = x))
  } else if (cond2){
    if (!file.exists(file1)) stop("File does not exist: ", file1)
    trees1 <- ape::read.tree(file1)
    if (!file.exists(file2)) stop("File does not exist: ", file2)
    trees2 <- ape::read.tree(file2)
  } else if (cond3){
    if (!file.exists(file1)) stop("File does not exist: ", file1)
    trees1 <- ape::read.tree(file1)
    trees2 <- lapply(newicks2, function(x) ape::read.tree(text = x))
  } else if (cond4){
    if (!file.exists(file2)) stop("File does not exist: ", file2)
    trees2 <- ape::read.tree(file2)
    trees1 <- lapply(newicks1, function(x) ape::read.tree(text = x))
  }
  
  n1 <- length(trees1)
  n2 <- length(trees2)
  
  if (n1 != n2){
    if (n2 != 1){
      stop("Number of trees in both lists should be the same or there should be only one tree in the second list")
    }
  }
  
  
  ## --- Normalize: unroot and remove edge lengths ---------------------------
  trees1 <- lapply(trees1, function(tr) {
    # Unroot the tree
    if (!is.null(tr$root.edge) || !ape::is.rooted(tr)) {
      tr <- ape::unroot(tr)
    }
    
    # Remove branch lengths (set to NULL so ape::write.tree does not print them)
    tr$edge.length <- NULL
    tr
  })
  
  trees2 <- lapply(trees2, function(tr) {
    # Unroot the tree
    if (!is.null(tr$root.edge) || !ape::is.rooted(tr)) {
      tr <- ape::unroot(tr)
    }
    
    # Remove branch lengths (set to NULL so ape::write.tree does not print them)
    tr$edge.length <- NULL
    tr
  })
  
  cleaned_newicks1 <- vapply(trees1, ape::write.tree, FUN.VALUE = character(1))
  cleaned_newicks2 <- vapply(trees2, ape::write.tree, FUN.VALUE = character(1))
  
  ## --- Computing the similarity for every tree in the first list against tree(s) in second list ---
  
  if (n2 == 1){
    res <- lapply(cleaned_newicks1, function(x) computeFD(x, cleaned_newicks2[1]));
  } else {
    res <- lapply((1:n1), function(x) computeFD(cleaned_newicks1[x],cleaned_newicks2[x]))
  }
  
  return(res)
}