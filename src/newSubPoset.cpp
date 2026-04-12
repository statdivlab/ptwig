#include "pTree.h"
#include "mPhylo.h"
#include "rho.h"
#include "coverTrees.h"
#include "stableSearch.h"
#include "subPoset.h"
#include "FDRSearch.h"
#include "idNullCoveringPairsComputation.h"
#include "subPosetAnalysis.h"
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
Rcpp::List computeNewSubposet(CharacterVector treeSampleR,
                                 IntegerVector nSampleR,
                                 CharacterVector compLeafSetR,
                                 int MtR, int rbR) {
    
    std::vector<pTree> treeSample;
    treeSample.reserve(treeSampleR.size());

    for (int i = 0; i < treeSampleR.size(); i++) {
        if (treeSampleR[i] == NA_STRING)
            stop("treeSample cannot contain NA.");
        treeSample.emplace_back(
            pTree(as<std::string>(treeSampleR[i]))
        );
    }
    
    std::vector<int> nSample;
    
    for (int i = 0; i < nSampleR.size(); i++){
        nSample.push_back(static_cast<int>(nSampleR[i]));
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
    cout << "We are over here! \n";
    
    //
    // 2.1 Converting the integers
    //
    
    int Mt = static_cast<int>(MtR);
    
    
    int rb = static_cast<int>(rbR);

    //
    // 3. Call C++ function
    //
    
    cout<< "Right about to enter the subPoset builder \n";
    
    subPoset subPost = subPoset(treeSample, nSample, compLeafSet, Mt, rb);
    
    cout<< "Got out of the builder \n";
    
    std::vector<std::pair<int,int>> edges;
    
    
    //Finding the pairs with upper tree of rank 1 and forming the edges 
    
    int curIndx = subPost.firstRank.at(1);
    while (curIndx > -1){
        edges.push_back({-1, (curIndx+1)});
        curIndx = subPost.Poset.at(curIndx).next;
    }
    
    //Finding the pairs above
    
    for (int i = 0; i < subPost.Poset.size(); i++){
        for (int j : subPost.Poset.at(i).over){
            edges.push_back({(i+1),(j+1)});
        }
    }
    
    vector<string> subPosetTrees;
    vector<int> subPosetRanks;
    
    for (int k = 0; k < subPost.Poset.size(); k++){
        mPhylo rP = mPhylo(subPost.Poset.at(k).Tree);
        subPosetTrees.push_back(rP.toNewick());
        subPosetRanks.push_back(subPost.Poset.at(k).Tree.rank);
    }
    
    // Convert edges: split pairs into two parallel integer vectors
    int nEdges = edges.size();
    Rcpp::IntegerVector edgeFrom(nEdges), edgeTo(nEdges);
    for (int i = 0; i < nEdges; i++) {
        edgeFrom[i] = edges[i].first;
        edgeTo[i]   = edges[i].second;
    }
    
    
    return Rcpp::List::create(
        Rcpp::Named("subPosetTrees")       = Rcpp::wrap(subPosetTrees),
        Rcpp::Named("subPosetRanks")       = Rcpp::wrap(subPosetRanks),
        Rcpp::Named("CoveringPairsLower")  = edgeFrom,
        Rcpp::Named("CoveringPairsUpper")  = edgeTo
    );
    
}