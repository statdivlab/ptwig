#include "pTree.h"
#include "mPhylo.h"
#include "rho.h"
#include "coverTrees.h"
#include <Rcpp.h>
#include <iostream>
#include <set>
#include <queue>
#include <vector>
#include <stack>
#include <fstream>
#include <iostream>
#include <string>
#include <chrono>
#include <random>
#include <variant>
#include <utility>
#include <filesystem> // C++17
#include <future>
#include <mutex>
#include <cmath>

using namespace Rcpp;
using namespace std;

// [[Rcpp::export]]
int computeFD(CharacterVector treeR1, CharacterVector treeR2) {
    
    pTree tree1 = pTree(as<std::string>(treeR1));
    pTree tree2 = pTree(as<std::string>(treeR2));
    
    return (tree1.rank - rho(tree1, tree2));
}