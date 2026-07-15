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
CharacterVector completeSearch_stability(CharacterVector treeSample1R,
                                 IntegerVector nSample1R,
                                 CharacterVector treeSample2R,
                                 IntegerVector nSample2R,
                                 CharacterVector compLeafSetR,
                                 double alphaR, double qR, double tauR) {
    cout << "It entered the first Cpp function \n"<< std::flush;
    std::vector<pTree> treeSample1;
    treeSample1.reserve(treeSample1R.size());

    for (int i = 0; i < treeSample1R.size(); i++) {
        if (treeSample1R[i] == NA_STRING)
            stop("treeSample cannot contain NA.");
        treeSample1.emplace_back(
            pTree(as<std::string>(treeSample1R[i]))
        );
    }
    
    cout << "Created the first Sample \n"<< std::flush;
    std::vector<pTree> treeSample2;
    treeSample2.reserve(treeSample2R.size());

    for (int i = 0; i < treeSample2R.size(); i++) {
        if (treeSample2R[i] == NA_STRING)
            stop("treeSample cannot contain NA.");
        treeSample2.emplace_back(
            pTree(as<std::string>(treeSample2R[i]))
        );
    }
    
    cout << "Created the second Sample \n"<< std::flush;
    
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
    
    int B2 = std::accumulate(nSample2.begin(), nSample2.end(), 0);
    
    //
    // 2.1 Converting the remaining floats
    //
    
    float alpha = static_cast<float>(alphaR);
    float q = static_cast<float>(qR);
    float tau = static_cast<float>(tauR);

    //
    // 3. Call C++ function
    //
    
    cout<< "Computing Stable Trees for SubPoset build-up \n"<< std::flush;
    std::vector<pTree> stTrees = stableSearch(treeSample1, nSample1, compLeafSet, alpha);
    
    
    int stM = static_cast<int>(stTrees.size());
    int rank_anchor = floor((log(q) - log(stM) + B2*(tau-0.5)*(tau-0.5))/log(2));
    
    cout<< "About to build SubPoset, with Mt = " << stM << " and rb = "<< rank_anchor <<" \n"<< std::flush;
    subPoset subPost = subPoset(stTrees, treeSample1, nSample1, compLeafSet, rank_anchor);
    
    //----- Prepare bounds and the treeSample2 oRho cache ------------------//
    // Each node's zeta was computed inside the subposet constructor (over its
    // covers, from treeSample1). We only need the max-level antichain bounds
    // and the per-node oRho cache against treeSample2 that FDRSearchGreedy uses.
    computeAllMaxLevelBounds(subPost);

    cout << "Building oRho cache for treeSample2 \n"<< std::flush;

    int nNodes = static_cast<int>(subPost.Poset.size());
    int nS2 = static_cast<int>(treeSample2.size());

    pTree emptyTree = pTree("();");
    vector<set<string>> emptyLeaves;
    oRho emptyORho = oRho(0, emptyLeaves);

    vector<vector<oRho>> storedORho2(nNodes, vector<oRho>(nS2, emptyORho));

    // Rank-1 nodes: incremental rho from the empty tree.
    int spCurIndx = subPost.firstRank.at(1);

    while (spCurIndx > -1) {
        spNode& spN = subPost.Poset[spCurIndx];
        pTree Ta = spN.Tree;

        for (int k = 0; k < nS2; k++){
            pTree T = treeSample2.at(k);
            oRho oRhoTa = rho(Ta, emptyTree, emptyORho, T);
            storedORho2.at(spCurIndx).at(k) = oRhoTa;
        }
        spCurIndx = subPost.Poset.at(spCurIndx).next;
    }

    // Higher ranks: incremental rho from the node directly below.
    for (int crank = 2; crank < subPost.firstRank.size(); crank++){
        spCurIndx = subPost.firstRank.at(crank);

        while (spCurIndx > -1) {
            spNode& spN = subPost.Poset[spCurIndx];
            pTree Ta = spN.Tree;

            int underIdx = spN.under.at(0);
            pTree Tunder = subPost.Poset.at(underIdx).Tree;

            for (int k = 0; k < nS2; k++){
                pTree T = treeSample2.at(k);
                oRho baseORho = storedORho2.at(underIdx).at(k);
                oRho oRhoTa = rho(Ta, Tunder, baseORho, T);
                storedORho2.at(spCurIndx).at(k) = oRhoTa;
            }
            spCurIndx = subPost.Poset.at(spCurIndx).next;
        }
    }

    cout << "About to enter FDR-controlled tree search \n"<< std::flush;
    pTree FDRTree = FDRSearchGreedy(treeSample2, nSample2, storedORho2, subPost, q);

    CharacterVector out(1);
    if (FDRTree.rank == 0){
        out[0] = "();";
    } else {
        out[0] = mPhylo(FDRTree).toNewick();
    }

    return out;
    
}

// [[Rcpp::export]]
CharacterVector completeSearch_basic_bifurcation(CharacterVector treeSample1R,
                                 IntegerVector nSample1R,
                                 CharacterVector treeSample2R,
                                 IntegerVector nSample2R,
                                 CharacterVector compLeafSetR,
                                 int MtR, int rbR, double qR) {
    
    cout << "It entered the first Cpp function \n"<< std::flush;
    std::vector<pTree> treeSample1;
    treeSample1.reserve(treeSample1R.size());

    for (int i = 0; i < treeSample1R.size(); i++) {
        if (treeSample1R[i] == NA_STRING)
            stop("treeSample cannot contain NA.");
        treeSample1.emplace_back(
            pTree(as<std::string>(treeSample1R[i]))
        );
    }
    
    cout << "Created the first Sample \n"<< std::flush;
    
    std::vector<pTree> treeSample2;
    treeSample2.reserve(treeSample2R.size());

    for (int i = 0; i < treeSample2R.size(); i++) {
        if (treeSample2R[i] == NA_STRING)
            stop("treeSample cannot contain NA.");
        treeSample2.emplace_back(
            pTree(as<std::string>(treeSample2R[i]))
        );
    }
    
    cout << "Created the second Sample \n"<< std::flush;
    
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
    
    int Mt = static_cast<int>(MtR);
    
    int rb = static_cast<int>(rbR);
    
    float q = static_cast<float>(qR);

    //
    // 3. Call C++ function
    //
    
    cout<< "About to build SubPoset \n"<< std::flush;
    subPoset subPost = subPoset(treeSample1, nSample1, compLeafSet, Mt, rb);
    
    computeAllMaxLevelBounds(subPost);
    
    
    //----- Build the treeSample2 oRho cache -------------------------------//
    // Each node's zeta was computed inside the subposet constructor (over its
    // covers, from treeSample1). Here we only build storedORho2 - the per-node
    // oRho cache against treeSample2 - that FDRSearchGreedy needs.
    cout << "Building oRho cache for treeSample2 \n"<< std::flush;

    int nNodes = static_cast<int>(subPost.Poset.size());
    int nS2 = static_cast<int>(treeSample2.size());

    pTree emptyTree = pTree("();");
    vector<set<string>> emptyLeaves;
    oRho emptyORho = oRho(0, emptyLeaves);

    vector<vector<oRho>> storedORho2(nNodes, vector<oRho>(nS2, emptyORho));

    int spCurIndx = subPost.firstRank.at(1);

    while (spCurIndx > -1) {
        spNode& spN = subPost.Poset[spCurIndx];
        pTree Ta = spN.Tree;

        for (int k = 0; k < nS2; k++){
            pTree T = treeSample2.at(k);
            oRho oRhoTa = rho(Ta, emptyTree, emptyORho, T);
            storedORho2.at(spCurIndx).at(k) = oRhoTa;
        }
        spCurIndx = subPost.Poset.at(spCurIndx).next;
    }

    for (int crank = 2; crank < subPost.firstRank.size(); crank++){
        spCurIndx = subPost.firstRank.at(crank);

        while (spCurIndx > -1) {
            spNode& spN = subPost.Poset[spCurIndx];
            pTree Ta = spN.Tree;

            int underIdx = spN.under.at(0);
            pTree Tunder = subPost.Poset.at(underIdx).Tree;

            for (int k = 0; k < nS2; k++){
                pTree T = treeSample2.at(k);
                oRho baseORho = storedORho2.at(underIdx).at(k);
                oRho oRhoTa = rho(Ta, Tunder, baseORho, T);
                storedORho2.at(spCurIndx).at(k) = oRhoTa;
            }
            spCurIndx = subPost.Poset.at(spCurIndx).next;
        }
    }

    cout << "About to enter FDR-controlled tree search \n"<< std::flush;
    pTree FDRTree = FDRSearchGreedy(treeSample2, nSample2, storedORho2, subPost, q);

    CharacterVector out(1);
    if (FDRTree.rank == 0){
        out[0] = "();";
    } else {
        out[0] = mPhylo(FDRTree).toNewick();
    }

    return out;
    
}


// [[Rcpp::export]]
CharacterVector completeSearch_basic_score(CharacterVector treeSample1R,
                                 IntegerVector nSample1R,
                                 CharacterVector treeSample2R,
                                 IntegerVector nSample2R,
                                 CharacterVector compLeafSetR,
                                 double qR, int top_widthR, int bottom_widthR,
                                 std::string orientationR) {


    cout << "It entered the first Cpp function \n"<< std::flush;
    std::vector<pTree> treeSample1;
    treeSample1.reserve(treeSample1R.size());

    for (int i = 0; i < treeSample1R.size(); i++) {
        if (treeSample1R[i] == NA_STRING)
            stop("treeSample cannot contain NA.");
        treeSample1.emplace_back(
            pTree(as<std::string>(treeSample1R[i]))
        );
    }

    cout << "Created the first Sample \n"<< std::flush;

    std::vector<pTree> treeSample2;
    treeSample2.reserve(treeSample2R.size());

    for (int i = 0; i < treeSample2R.size(); i++) {
        if (treeSample2R[i] == NA_STRING)
            stop("treeSample cannot contain NA.");
        treeSample2.emplace_back(
            pTree(as<std::string>(treeSample2R[i]))
        );
    }

    cout << "Created the second Sample \n"<< std::flush;

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
    // 2.1 Converting the remaining parameters
    //

    float q = static_cast<float>(qR);

    int top_width    = static_cast<int>(top_widthR);

    int bottom_width = static_cast<int>(bottom_widthR);

    //
    // 3. Call C++ function
    //

    cout<< "About to build SubPoset \n"<< std::flush;

    subPoset subPost = subPoset(treeSample1, nSample1, compLeafSet,
                                top_width, bottom_width, orientationR);

    computeAllMaxLevelBounds(subPost);
    
    //subPost.print();

    //----- Build oRho cache for treeSample2 -------------//
    // The fixed-width constructor already computes each node's zeta exactly
    // (over all of its covers, from treeSample1), so the zeta re-estimation
    // that V3 did here is redundant and dropped. We still build storedORho2 —
    // the per-node oRho cache against treeSample2 — which FDRSearchGreedy needs.
    cout << "Building oRho cache for treeSample2 \n"<< std::flush;

    int nNodes = static_cast<int>(subPost.Poset.size());
    int nS2 = static_cast<int>(treeSample2.size());

    // Base oRho for the empty tree
    pTree emptyTree = pTree("();");
    vector<set<string>> emptyLeaves;
    oRho emptyORho = oRho(0, emptyLeaves);

    vector<vector<oRho>> storedORho2(nNodes, vector<oRho>(nS2, emptyORho));

    //Computing Rhos for the trees in Subposet at rank 1;
    int spCurIndx = subPost.firstRank.at(1);

    while (spCurIndx > -1) {
        spNode& spN = subPost.Poset[spCurIndx];
        pTree Ta = spN.Tree;

        // --- treeSample2 ---
        for (int k = 0; k < nS2; k++){
            pTree T = treeSample2.at(k);

            oRho oRhoTa = rho(Ta, emptyTree, emptyORho, T);

            storedORho2.at(spCurIndx).at(k) = oRhoTa;
        }
        spCurIndx = subPost.Poset.at(spCurIndx).next;
    }

    for (int crank = 2; crank < subPost.firstRank.size(); crank++){
        spCurIndx = subPost.firstRank.at(crank);

        while (spCurIndx > -1) {
            spNode& spN = subPost.Poset[spCurIndx];
            pTree Ta = spN.Tree;

            int underIdx = spN.under.at(0);
            pTree Tunder = subPost.Poset.at(underIdx).Tree;

            // --- treeSample2 ---
            for (int k = 0; k < nS2; k++){
                pTree T = treeSample2.at(k);
                oRho baseORho = storedORho2.at(underIdx).at(k);
                oRho oRhoTa = rho(Ta, Tunder, baseORho, T);

                storedORho2.at(spCurIndx).at(k) = oRhoTa;
            }
            spCurIndx = subPost.Poset.at(spCurIndx).next;
        }
    }

    cout << "About to enter FDR-controlled tree search \n"<< std::flush;
    pTree FDRTree = FDRSearchGreedy(treeSample2, nSample2, storedORho2, subPost, q);

    //
    // 4. Convert vector<pTree> → CharacterVector
    //
    CharacterVector out(1);
    
    
    if (FDRTree.rank == 0){
        out[0] = "();";
    } else {
        out[0] = mPhylo(FDRTree).toNewick();
    }

    return out;

}
