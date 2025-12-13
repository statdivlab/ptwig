#' Run stable search
#'
#' @param newicks Character vector of Newick strings (optional)
#' @param newicks1 Character vector of Newick strings for Stable Search (optional)
#' @param newicks2 Character vector of Newick strings for FDR control Search (optional)
#' @param file Path to a file containing Newick trees (optional)
#' @param file1 Path to a file containing Newick trees for Stable Search (optional)
#' @param file2 Path to a file containing Newick trees for FDR control Search (optional)
#' @param n1 Number of trees to be used for Stable Search (optional)
#' @param random_subsampling Boolean indicating if the subsampling is to be made random 
#' @param alpha Numeric value passed to Stable Search for stable threshold.
#' @param q Numeric value passed to FDR control purposes.
#' @param tau Extra value for subposet building.
#' @param summarized Boolean factor indicating if the function is to be runned with a summarized version of the sample
#'
#' @return Output of completeSearchRcpp
#' @export
#'
#' @importFrom ape read.tree
run_FDRcontrol_search <- function(newicks = NULL, newicks1 = NULL, newicks2 = NULL,
                                  file = NULL, file1 = NULL, file2 = NULL, 
                                  n1 = NULL,
                                  random_subsampling = FALSE,
                                  alpha = 0.85, q = 0.1, tau = 0.95, 
                                  summarized = FALSE) {
  
  ## --- Argument validation -------------------------------------------------
  
  # Define the six allowed conditions as logical variables
  cond1 <- !is.null(newicks) # (newicks != NULL)
  
  cond2 <- !is.null(newicks1) && !is.null(newicks2) # (newicks1 !=NULL && newicks2 != NULL)
  
  cond3 <- !is.null(file) # (file != NULL)
  
  cond4 <- !is.null(file1) && !is.null(file2) # (file1 != NULL && file2 != NULL)
  
  cond5 <- !is.null(file1) && !is.null(newicks2) # (file1 != NULL && newicks2 != NULL)
  
  cond6 <- !is.null(file2) && !is.null(newicks1) # (file2 != NULL && newicks1 != NULL)
  
  # Create a vector of the conditions
  conditions <- c(cond1, cond2, cond3, cond4, cond5, cond6)
  
  # Count how many conditions are TRUE
  num_true_conditions <- sum(conditions)
  
  # Check if exactly one condition is TRUE
  if (num_true_conditions != 1) {
    # If zero conditions are true
    if (num_true_conditions == 0) {
      stop(
        "You must provide arguments satisfying EXACTLY ONE of the following conditions:\n",
        "1. `newicks` is provided\n",
        "2. `newicks1` AND `newicks2` are provided\n",
        "3. `file` is provided\n",
        "4. `file1` AND `file2` are provided\n",
        "5. `file1` AND `newicks2` are provided\n",
        "6. `file2` AND `newicks1` are provided"
      )
    } 
    # If more than one condition is true
    else {
      stop(
        "You have provided arguments satisfying MORE THAN ONE of the allowed combinations. ",
        "Please provide arguments satisfying EXACTLY ONE of the following conditions:\n",
        "1. `newicks` is provided\n",
        "2. `newicks1` AND `newicks2` are provided\n",
        "3. `file` is provided\n",
        "4. `file1` AND `file2` are provided\n",
        "5. `file1` AND `newicks2` are provided\n",
        "6. `file2` AND `newicks1` are provided"
      )
    }
  }
  
  ## --- Read trees from files (if provided) ---------------------------------
  if (cond1 || cond3){
    if (!is.null(file)) {
      if (!file.exists(file)) stop("File does not exist: ", file)
      trees <- ape::read.tree(file)
    } else {
      ## Parse Newick strings directly
      trees <- lapply(newicks, ape::read.tree, text = TRUE)
    }
    
    if (is.null(n1)){
      n1 = floor(length(trees)/2)
    }
    
    if (random_subsampling) {
      # randomly sample indices
      idx1 <- sample(length(trees), n1)
    } else {
      # take the first n1 indices
      idx1 <- seq_len(n1)
    }
    
    # create the two subsets
    trees1 <- trees[idx1]
    trees2 <- trees[-idx1]
  } else if (cond2){
    trees1 <- lapply(newicks1, ape::read.tree, text = TRUE)
    trees2 <- lapply(newicks2, ape::read.tree, text = TRUE)
  } else if (cond4){
    if (!file.exists(file1)) stop("File does not exist: ", file1)
    trees1 <- ape::read.tree(file1)
    if (!file.exists(file2)) stop("File does not exist: ", file2)
    trees2 <- ape::read.tree(file2)
  } else if (cond5){
    if (!file.exists(file1)) stop("File does not exist: ", file1)
    trees1 <- ape::read.tree(file1)
    trees2 <- lapply(newicks2, ape::read.tree, text = TRUE)
  } else if (cond6){
    if (!file.exists(file2)) stop("File does not exist: ", file2)
    trees2 <- ape::read.tree(file2)
    trees1 <- lapply(newicks1, ape::read.tree, text = TRUE)
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
  
  ## --- Compute union of all tip labels -------------------------------------
  completeLeaveSet <- unique(unlist(lapply(c(trees1,trees2), function(x) x$tip.label)))
  
  ## --- Preparing to run the final function -----------------------------------
  if(summarized){
    ## --- Joining trees that are identical if summarized = TRUE --------------
    Unique_trees1 = list()
    Count_trees1 = c()
    
    Unique_trees2 = list()
    Count_trees2 = c()
    
    for (tree in trees1){
      Found = FALSE;
      i = 0;
      for (Top in Unique_trees1) {
        i = i+1
        if (all.equal(tree, Top)){
          Found = TRUE;
          Count_trees1[i] = Count_trees1[i] + 1;
        }
      }
      if (!Found){
        Unique_trees1[[length(Unique_trees1)+1]] = tree;
        Count_trees1 = c(Count_trees1, 1);
      }
    }
    
    for (tree in trees2){
      Found = FALSE;
      i = 0;
      for (Top in Unique_trees2) {
        i = i+1
        if (all.equal(tree, Top)){
          Found = TRUE;
          Count_trees2[i] = Count_trees2[i] + 1;
        }
      }
      if (!Found){
        Unique_trees2[[length(Unique_trees2)+1]] = tree;
        Count_trees2 = c(Count_trees2, 1);
      }
    }
    
    ## --- Write cleaned trees back to Newick strings ---------------------------
    cleaned_newicks1 <- vapply(Unique_trees1, ape::write.tree, FUN.VALUE = character(1))
    cleaned_newicks2 <- vapply(Unique_trees2, ape::write.tree, FUN.VALUE = character(1))
    
    ## --- Call your Rcpp backend ----------------------------------------------
    res <- completeSearchRcppS(treeSample1R = cleaned_newicks1, nSample1R = Count_trees1, 
                               treeSample2R = cleaned_newicks2, nSample2R = Count_trees2, 
                               compLeafSetR = completeLeaveSet, 
                               alphaR = alpha, qR = q, tauR = tau)
  } else {
    ## --- Write cleaned trees back to Newick strings ---------------------------
    cleaned_newicks1 <- vapply(trees1, ape::write.tree, FUN.VALUE = character(1))
    cleaned_newicks2 <- vapply(trees2, ape::write.tree, FUN.VALUE = character(1))
    
    ## --- Call your Rcpp backend ----------------------------------------------
    res <- completeSearchRcpp(treeSample1R = cleaned_newicks1, 
                               treeSample2R = cleaned_newicks2,
                               compLeafSetR = completeLeaveSet, 
                               alphaR = alpha, qR = q, tauR = tau)
  }
  
  return(res)
}