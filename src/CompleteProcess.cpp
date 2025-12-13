#include "pTree.h"
#include "mPhylo.h"
#include "rho.h"
#include "coverTrees.h"
#include "stableSearch.h"
#include "subPoset.h"
#include "FDRSearch.h"
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
#include <Rcpp.h>

using namespace Rcpp;
using namespace std;

// [[Rcpp::export]]
CharacterVector completeSearchRcpp(CharacterVector treeSample1R,
                                 CharacterVector treeSample2R,
                                 CharacterVector compLeafSetR,
                                 double alphaR, double qR, double tauR) {
    
    std::vector<pTree> treeSample1;
    treeSample1.reserve(treeSample1R.size());

    for (int i = 0; i < treeSample1R.size(); i++) {
        if (treeSample1R[i] == NA_STRING)
            stop("treeSample cannot contain NA.");
        treeSample1.emplace_back(
            pTree(as<std::string>(treeSample1R[i]))
        );
    }
    
    std::vector<pTree> treeSample2;
    treeSample2.reserve(treeSample2R.size());

    for (int i = 0; i < treeSample2R.size(); i++) {
        if (treeSample2R[i] == NA_STRING)
            stop("treeSample cannot contain NA.");
        treeSample2.emplace_back(
            pTree(as<std::string>(treeSample2R[i]))
        );
    }
    

    //
    // 2. Convert compLeafSetR → set<string>
    //
    std::set<std::string> compLeafSet;
    for (int i = 0; i < compLeafSetR.size(); i++) {
        if (compLeafSetR[i] == NA_STRING)
            stop("compLeafSet cannot contain NA.");
        compLeafSet.insert(as<std::string>(compLeafSetR[i]));
    }
    
    //
    // 2.1 Converting the remaining floats
    //
    
    float alpha = static_cast<float>(alphaR);
    float q = static_cast<float>(qR);
    float tau = static_cast<float>(tau);

    //
    // 3. Call C++ function
    //
    
    std::vector<pTree> stTrees = stableSearch(treeSample1, compLeafSet, alpha);
    
    
    int stM = static_cast<int>(stTrees.size());
    int B2 = static_cast<int>(treeSample2.size());
    int rank_anchor = floor((log(q) - log(stM) + B2*(tau-0.5)*(tau-0.5))/log(2));
    
    subPoset subPost = subPoset(stTrees, treeSample1, compLeafSet, rank_anchor);
    
    vector<pTree> FDRTrees = FDRSearch(treeSample2, subPost, q);

    //
    // 4. Convert vector<pTree> → CharacterVector
    //
    CharacterVector out(FDRTrees.size());
    for (size_t i = 0; i < FDRTrees.size(); i++) {
        mPhylo rP = mPhylo(FDRTrees[i]);
        out[i] = rP.toNewick();
    }

    return out;
    
    
}

// [[Rcpp::export]]
CharacterVector completeSearchRcppS(CharacterVector treeSample1R,
                                 IntegerVector nSample1R,
                                 CharacterVector treeSample2R,
                                 IntegerVector nSample2R,
                                 CharacterVector compLeafSetR,
                                 double alphaR, double qR, double tauR) {
    
    std::vector<pTree> treeSample1;
    treeSample1.reserve(treeSample1R.size());

    for (int i = 0; i < treeSample1R.size(); i++) {
        if (treeSample1R[i] == NA_STRING)
            stop("treeSample cannot contain NA.");
        treeSample1.emplace_back(
            pTree(as<std::string>(treeSample1R[i]))
        );
    }
    
    std::vector<pTree> treeSample2;
    treeSample2.reserve(treeSample2R.size());

    for (int i = 0; i < treeSample2R.size(); i++) {
        if (treeSample2R[i] == NA_STRING)
            stop("treeSample cannot contain NA.");
        treeSample2.emplace_back(
            pTree(as<std::string>(treeSample2R[i]))
        );
    }
    
    std::vector<int> nSample1;
    
    for (int i = 0; i < nSample1R.size(); i++){
        nSample1.push_back(static_cast<int>(nSample1R[i]));
    }
    
    std::vector<int> nSample2;
    
    for (int i = 0; i < nSample2R.size(); i++){
        nSample2.push_back(static_cast<int>(nSample2R[i]));
    }

    //
    // 2. Convert compLeafSetR → set<string>
    //
    std::set<std::string> compLeafSet;
    for (int i = 0; i < compLeafSetR.size(); i++) {
        if (compLeafSetR[i] == NA_STRING)
            stop("compLeafSet cannot contain NA.");
        compLeafSet.insert(as<std::string>(compLeafSetR[i]));
    }
    
    //
    // 2.1 Converting the remaining floats
    //
    
    float alpha = static_cast<float>(alphaR);
    float q = static_cast<float>(qR);
    float tau = static_cast<float>(tau);

    //
    // 3. Call C++ function
    //
    
    std::vector<pTree> stTrees = stableSearch(treeSample1, nSample1, compLeafSet, alpha);
    
    
    int stM = static_cast<int>(stTrees.size());
    int B2 = std::accumulate(nSample2.begin(), nSample2.end(), 0);
    int rank_anchor = floor((log(q) - log(stM) + B2*(tau-0.5)*(tau-0.5))/log(2));
    
    subPoset subPost = subPoset(stTrees, treeSample1, nSample1, compLeafSet, rank_anchor);
    
    vector<pTree> FDRTrees = FDRSearch(treeSample2, nSample2, subPost, q);

    //
    // 4. Convert vector<pTree> → CharacterVector
    //
    CharacterVector out(FDRTrees.size());
    for (size_t i = 0; i < FDRTrees.size(); i++) {
        mPhylo rP = mPhylo(FDRTrees[i]);
        out[i] = rP.toNewick();
    }

    return out;
    
}
