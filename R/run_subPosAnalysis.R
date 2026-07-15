#' Running a complete analysis of the properties of the subposet constructed from a set of samples
#'
#' @param tstar_newick Character vector of Newick string for the target tree (optional).
#' @param newicks1 Character vector of Newick strings for Stable Search (optional)
#' @param newicks2 Character vector of Newick strings for subposet Variability Analysis (optional)
#' @param bnewicks Character vector of Newick strings for probability computation (optional)
#' @param tstar_file Path to a file containing Newick string for the target tree (optional)
#' @param file1 Path to a file containing Newick trees for Stable Search (optional)
#' @param file2 Path to a file containing Newick trees for subposet Variability Analysis (optional)
#' @param bfile Path to a file containing Newick trees for probability computation (optional)
#' @param SPbuilder String indication what method to use to build subposet.
#' @param alpha Numeric value passed to Stable Search for stable threshold. (in stability-based SubPoset)
#' @param q Numeric value used in FDR control purposes, used in the construction of the subPoset in stability-based subposet
#' @param tau Extra value for subposet building in stability-based subposet
#' @param Mt Number of maximal trees for basic-bifurcation subposet building
#' @param rb Anchor rank for basic-bifurcation subposet building
#'
#' @return Output of SPAnalysisR
#' @export
#'
#' @importFrom ape read.tree
run_subPostAnalysis <- function(tstar_newick = NULL, newicks1 = NULL, newicks2 = NULL, bnewicks = NULL,
                                  tstar_file = NULL, file1 = NULL, file2 = NULL, bfile = NULL,
                                  SPbuilder = "stability",
                                  alpha = 0.85, q = 0.1, tau = 0.80,
                                  Mt = 1, rb = NULL, delta = NULL) {
  
  ## --- Argument validation -------------------------------------------------
  
  # Define the six allowed conditions as logical variables
  cond0 <- xor((!is.null(tstar_newick) && is.null(tstar_file)), (is.null(tstar_newick) && !is.null(tstar_file)))
  
  cond1 <- xor((!is.null(newicks1) && is.null(file1)), (is.null(newicks1) && !is.null(file1)))
  
  cond2 <- xor((!is.null(newicks2) && is.null(file2)), (is.null(newicks2) && !is.null(file2)))  
  
  cond3 <- xor((!is.null(bnewicks) && is.null(bfile)), (is.null(bnewicks) && !is.null(bfile)))  
  
  
  # Check if exactly one source for each tree/sample tree
  
  if(!cond0){
    stop(
      "You must provide the target tree in EXACTLY ONE of the following arguments:\n",
      "1. as a newick string in tstar_newick or 2. as a file containing the newick string in tstar_file."
    )
  }
  
  if (!cond1){
    stop(
      "You must provide the First Sample for SubPoset construction in EXACTLY ONE of the following arguments:\n",
      "1. as a collection of newick strings in newicks1 or 2. as a file containing the newick strings in file1."
    )
  }
  
  if (!cond2){
    stop(
      "You must provide the Second Sample for SubPoset variance analysis in EXACTLY ONE of the following arguments:\n",
      "1. as a collection of newick strings in newicks2 or 2. as a file containing the newick strings in file2."
    )
  }
  
  if (!cond2){
    stop(
      "You must provide a Large Sample for probability estimation in EXACTLY ONE of the following arguments:\n",
      "1. as a collection of newick strings in bnewicks or 2. as a file containing the newick strings in bfile."
    )
  }
  
  if ((SPbuilder != "stability") && (SPbuilder != "basic")){
    stop(
      "The subPoset builder method must be a valid option. It is either:\n",
      "SPbuilder = 'stability': Building a stable tree first and branching out from there\n",
      "SPbuilder = 'basic': For basic bifurcation with known number of maximal trees and bifurcation rank level."
    )
  }
  
  if (is.null(delta)){
    delta = q/2
  }
  
  if (delta >= q){
    stop(
      "The value of delta must be strictly lower than q \n"
    )
  }
  
  ## --- Read trees from files (if provided) ---------------------------------
  
  if (!is.null(tstar_newick)){
    tree_star <- lapply(tstar_newick, function(x) ape::read.tree(text = x))
  } else {
    if (!file.exists(tstar_file)) stop("File does not exist: ", tstar_file)
    tree_star <- ape::read.tree(tstar_file)
  }
  
  nt <- length(tree_star)
  if (nt>1){ 
    warning("You provided more than one target tree. Will use only the first in the list")
    tree_star1 <- tree_star[1]
  } else {
    tree_star1 <- tree_star[[1]]
  }
  
  if (!is.null(newicks1)){
    trees1 <- lapply(newicks1, function(x) ape::read.tree(text = x))
  } else {
    if (!file.exists(file1)) stop("File does not exist: ", file1)
      trees1 <- ape::read.tree(file1)
  }
  
  if (!is.null(newicks2)){
    trees2 <- lapply(newicks2, function(x) ape::read.tree(text = x))
  } else {
    if (!file.exists(file2)) stop("File does not exist: ", file2)
    trees2 <- ape::read.tree(file2)
  }
  
  if (!is.null(bnewicks)){
    btrees <- lapply(bnewicks, function(x) ape::read.tree(text = x))
  } else {
    if (!file.exists(bfile)) stop("File does not exist: ", bfile)
    btrees <- ape::read.tree(bfile)
  }
  
  
  ## --- Normalize: unroot and remove edge lengths ---------------------------
  tree_star1 <- ape::unroot(tree_star1)
  tree_star1$edge.length <- NULL
  tree_star1$node.label <- NULL
  cleaned_treeStar <- ape::write.tree(phy = tree_star1)
  
  trees1 <- lapply(trees1, function(tr) {
    # Unroot the tree
    if (!is.null(tr$root.edge) || !ape::is.rooted(tr)) {
      tr <- ape::unroot(tr)
    }
    
    # Remove branch lengths (set to NULL so ape::write.tree does not print them)
    tr$edge.length <- NULL
    tr$node.label <- NULL
    tr
  })
  
  trees2 <- lapply(trees2, function(tr) {
    # Unroot the tree
    if (!is.null(tr$root.edge) || !ape::is.rooted(tr)) {
      tr <- ape::unroot(tr)
    }
    
    # Remove branch lengths (set to NULL so ape::write.tree does not print them)
    tr$edge.length <- NULL
    tr$node.label <- NULL
    tr
  })
  
  btrees <- lapply(btrees, function(tr) {
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
  completeLeaveSet <- unique(unlist(lapply(c(trees1,trees2,btrees), function(x) x$tip.label)))
  
  ## --- Preparing to run the final function -----------------------------------
  ## --- Joining trees that are identical -----------------------------------
  Unique_trees1 = list()
  Count_trees1 = c()

  Unique_trees2 = list()
  Count_trees2 = c()

  Unique_btrees = list()
  Count_btrees = c()

  for (tree in trees1){
    Found = FALSE;
    i = 0;
    for (Top in Unique_trees1) {
      i = i+1
      if (all.equal(tree, Top)){
        Found = TRUE;
        Count_trees1[i] = Count_trees1[i] + 1;
        break
      }
    }
    if (!Found){
      Unique_trees1[[length(Unique_trees1)+1]] = tree;
      Count_trees1 = c(Count_trees1, 1);
    }
  }

  for (tree in btrees){
    Found = FALSE;
    i = 0;
    for (Top in Unique_btrees) {
      i = i+1
      if (all.equal(tree, Top)){
        Found = TRUE;
        Count_btrees[i] = Count_btrees[i] + 1;
        break
      }
    }
    if (!Found){
      Unique_btrees[[length(Unique_btrees)+1]] = tree;
      Count_btrees = c(Count_btrees, 1);
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
        break
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
  cleaned_newicksb <- vapply(Unique_btrees, ape::write.tree, FUN.VALUE = character(1))

  ## --- Call your Rcpp backend ----------------------------------------------
  if (SPbuilder == "stability"){
    res <-  SPAnalysisRS(treeStar = cleaned_treeStar,
                         treeSample1R = cleaned_newicks1, nSample1R = Count_trees1,
                         treeSample2R = cleaned_newicks2, nSample2R = Count_trees2,
                         bigTreeSampleR = cleaned_newicksb, nBSampleR = Count_btrees,
                         compLeafSetR = completeLeaveSet,
                         alphaR = alpha, qR = q, tauR = tau, deltaR = delta)
  } else if (SPbuilder == "basic") {
    if (is.null(rb)){
      rb = length(completeLeaveSet) - 4;
    }

    res <-  SPAnalysisR2S(treeStar = cleaned_treeStar,
                  treeSample1R = cleaned_newicks1, nSample1R = Count_trees1,
                  treeSample2R = cleaned_newicks2, nSample2R = Count_trees2,
                  bigTreeSampleR = cleaned_newicksb, nBSampleR = Count_btrees,
                  compLeafSetR = completeLeaveSet,
                  MtR = Mt, rbR = rb, qR = q, deltaR = delta)
  }
  
  return(res)
}