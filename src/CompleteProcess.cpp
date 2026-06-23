#include "pTree.h"
#include "mPhylo.h"
#include "rho.h"
#include "coverTrees.h"
#include "stableSearch.h"
#include "subPoset.h"
#include "FDRSearch.h"
#include <rcpptimer.h>
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
    
    cout << "Created the first Sample \n" << std::flush; 
    
    int B1 = treeSample1R.size();
    
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
    float tau = static_cast<float>(tauR);

    //
    // 3. Call C++ function
    //
    
    cout<< "Computing Stable Trees for SubPoset build-up \n"<< std::flush;
    std::vector<pTree> stTrees = stableSearch(treeSample1, compLeafSet, alpha);
    
    
    int stM = static_cast<int>(stTrees.size());
    int B2 = static_cast<int>(treeSample2.size());
    int rank_anchor = floor((log(q) - log(stM) + B2*(tau-0.5)*(tau-0.5))/log(2));
    
    cout<< "About to build SubPoset, with Mt = " << stM << " and rb = "<< rank_anchor <<" \n"<< std::flush;
    subPoset subPost = subPoset(stTrees, treeSample1, compLeafSet, rank_anchor);
    
    //----- Computing eta's lower bounds -------------//
    
     cout<< "Computing Etas' lower bounds \n"<< std::flush;
    
    vector<float> etasLowerBounds;
    
    int counterEtas = 0;
    for (spNode spN : subPost.Poset){
        pTree Ta = spN.Tree;
        
        vector<pTree> Tbs = coverTrees(Ta, compLeafSet);
        
        int nTbs = static_cast<int>(Tbs.size());
        
        int sumOfIndicators = 0;
        
        for (pTree T : treeSample1){
            int rhoTa = rho(Ta, T);
            for (pTree Tb : Tbs){
                sumOfIndicators += (int)(rho(Tb,T) > rhoTa);
            }
        }
        for (pTree T : treeSample2){
            int rhoTa = rho(Ta,T);
            for (pTree Tb : Tbs){
                sumOfIndicators += (int)(rho(Tb,T) > rhoTa);
            }
        }
        
        etasLowerBounds.push_back(std::max(1 - static_cast<float>(sumOfIndicators)/(static_cast<float>(nTbs*(B1+B2))), 0.5f));
        cout<< "New Lower Bound Added "<< counterEtas <<" \n"<< std::flush;
        counterEtas++;
    }
    
    {
        pTree Ta = pTree("();");
        
        vector<pTree> Tbs = coverTrees(Ta, compLeafSet);
        
        int nTbs = static_cast<int>(Tbs.size());
        
        int sumOfIndicators = 0;
        
        for (pTree T : treeSample1){
            int rhoTa = rho(Ta, T);
            for (pTree Tb : Tbs){
                sumOfIndicators += (int)(rho(Tb,T) > rhoTa);
            }
        }
        for (pTree T : treeSample2){
            int rhoTa = rho(Ta,T);
            for (pTree Tb : Tbs){
                sumOfIndicators += (int)(rho(Tb,T) > rhoTa);
            }
        }
        
        etasLowerBounds.push_back(std::max(1 - static_cast<float>(sumOfIndicators)/(static_cast<float>(nTbs*(B1+B2))),0.5f));
        cout<< "New Lower Bound Added "<< counterEtas <<" \n"<< std::flush;
        counterEtas++;
    }
    
    // Complete the FDR search//
    
    cout<< "Starting FDR-controled search"<< counterEtas <<" \n"<< std::flush;
    vector<pTree> FDRTrees = FDRSearch(treeSample2, subPost, etasLowerBounds, q);

    //
    // 4. Convert vector<pTree> → CharacterVector
    //
    CharacterVector out(FDRTrees.size());
    for (size_t i = 0; i < FDRTrees.size(); i++) {
        if (FDRTrees[i].rank == 0){
            out[i] = "();";
        } else {
            mPhylo rP = mPhylo(FDRTrees[i]);
            out[i] = rP.toNewick();
        }
    }

    return out;
    
    
}

// [[Rcpp::export]]
CharacterVector completeSearchRcpp_V2(CharacterVector treeSample1R,
                                 CharacterVector treeSample2R,
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
    subPoset subPost = subPoset(treeSample1, compLeafSet, Mt, rb);
    
    subPost.print();
    
    //----- Computing eta's lower bounds -------------//
    
    cout << "Computing etasLowerBounds \n"<< std::flush;
    
    // For each node index, store the oRho objects for each sample tree.
    // storedORho1[i][k] = oRho for (Poset[i].Tree, treeSample1[k])
    // storedORho2[i][k] = oRho for (Poset[i].Tree, treeSample2[k])
    int nNodes = static_cast<int>(subPost.Poset.size());
    int nS1 = static_cast<int>(treeSample1.size());
    int nS2 = static_cast<int>(treeSample2.size());
    
    // Base oRho for the empty tree
    pTree emptyTree = pTree("();");
    vector<set<string>> emptyLeaves;
    oRho emptyORho = oRho(0, emptyLeaves);

    vector<vector<oRho>> storedORho1(nNodes, vector<oRho>(nS1, emptyORho));
    vector<vector<oRho>> storedORho2(nNodes, vector<oRho>(nS2, emptyORho));
    
    //Computing Rhos for the trees in Subposet at rank 1;
    int spCurIndx = subPost.firstRank.at(1);
    
    while (spCurIndx > -1) {
        spNode& spN = subPost.Poset[spCurIndx];
        pTree Ta = spN.Tree;
        
        // --- treeSample1 ---
        for (int k = 0; k < nS1; k++){
            pTree T = treeSample1.at(k);
            oRho oRhoTa = rho(Ta, emptyTree, emptyORho, T);

            // Store for use by nodes above Ta
            storedORho1.at(spCurIndx).at(k) = oRhoTa;
        }

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

            // --- treeSample1 ---
            for (int k = 0; k < nS1; k++){
                pTree T = treeSample1.at(k);
                oRho baseORho = storedORho1.at(underIdx).at(k);
                oRho oRhoTa = rho(Ta, Tunder, baseORho, T);
                
                // Store for use by nodes above Ta
                storedORho1.at(spCurIndx).at(k) = oRhoTa;
            }

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
    
    int nodesCount = 0;
    vector<float> etasLowerBounds;
    for (int idx = 0; idx < nNodes; idx++){
        spNode& spN = subPost.Poset[idx];
        pTree Ta = spN.Tree;

        vector<pTree> Tbs = coverTrees(Ta, compLeafSet);
        int nTbs = static_cast<int>(Tbs.size());
        int sumOfIndicators = 0;

        // --- treeSample1 ---
        for (int k = 0; k < nS1; k++){
            pTree T = treeSample1.at(k);

            // Compute rho(Tb, T) for each cover tree, reusing oRhoTa
            for (pTree Tb : Tbs){
                oRho oRhoTb = rho(Tb, Ta, storedORho1.at(idx).at(k), T);
                sumOfIndicators += ((int)(oRhoTb.rho > storedORho1.at(idx).at(k).rho));
            }
        }

        // --- treeSample2 ---
        for (int k = 0; k < nS2; k++){
            pTree T = treeSample2.at(k);

            for (pTree Tb : Tbs){
                oRho oRhoTb = rho(Tb, Ta, storedORho2.at(idx).at(k), T);
                sumOfIndicators += ((int)(oRhoTb.rho > storedORho2.at(idx).at(k).rho));
            }
        }

        etasLowerBounds.push_back(std::max(1 - static_cast<float>(sumOfIndicators) /
                                  (static_cast<float>(nTbs * (nS1 + nS2))), 0.5f));

        cout << "Eta Lower Bound number " << nodesCount << " computed \n" << std::flush;
        nodesCount++;
    }

    // --- Empty tree case ---
    {
        vector<pTree> Tbs = coverTrees(emptyTree, compLeafSet);
        int nTbs = static_cast<int>(Tbs.size());
        int sumOfIndicators = 0;

        for (int k = 0; k < nS1; k++){
            pTree T = treeSample1.at(k);
            for (pTree Tb : Tbs){
                oRho oRhoTb = rho(Tb, emptyTree, emptyORho, T);
                sumOfIndicators += ((int)(oRhoTb.rho > 0));
            }
        }

        for (int k = 0; k < nS2; k++){
            pTree T = treeSample2.at(k);
            for (pTree Tb : Tbs){
                oRho oRhoTb = rho(Tb, emptyTree, emptyORho, T);
                sumOfIndicators += ((int)(oRhoTb.rho > 0));
            }
        }

        etasLowerBounds.push_back(std::max(1 - static_cast<float>(sumOfIndicators) /
                                  (static_cast<float>(nTbs * (nS1 + nS2))), 0.5f));
        cout << "Eta Lower Bound number " << nodesCount << " computed \n" << std::flush;
        nodesCount++;
    }
    
    // Complete the FDR search//
    cout<< "Starting FDR-controled search \n"<< std::flush;
    pTree FDRTree = FDRSearchGreedy(treeSample2, storedORho2, subPost, etasLowerBounds, q);

    //
    // 4. Convert vector<pTree> → CharacterVector
    //
    CharacterVector out(1);
    //for (size_t i = 0; i < FDRTrees.size(); i++) {
    //    mPhylo rP = mPhylo(FDRTrees[i]);
    //    out[i] = rP.toNewick();
    //}
    
    
    if (FDRTree.rank == 0){
        out[0] = "();";
    } else {
        out[0] = mPhylo(FDRTree).toNewick();
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
    
    int B1 = std::accumulate(nSample1.begin(), nSample1.end(), 0);
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
    
    //----- Computing eta's lower bounds -------------//
    
    cout<< "Computing Etas' lower bounds \n"<< std::flush;
    vector<float> etasLowerBounds;
    
    int counterEtas = 0;
    for (spNode spN : subPost.Poset){
        pTree Ta = spN.Tree;
        
        vector<pTree> Tbs = coverTrees(Ta, compLeafSet);
        
        int nTbs = static_cast<int>(Tbs.size());
        
        int sumOfIndicators = 0;
        
        for (int k = 0; k < treeSample1.size(); k++){
            pTree T = treeSample1.at(k);
            int rhoTa = rho(Ta, T);
            for (pTree Tb : Tbs){
                sumOfIndicators += nSample1.at(k)*((int)(rho(Tb,T) > rhoTa));
            }
        }
        for (int k = 0; k < treeSample2.size(); k++){
            pTree T = treeSample2.at(k);
            int rhoTa = rho(Ta,T);
            for (pTree Tb : Tbs){
                sumOfIndicators += nSample2.at(k)*((int)(rho(Tb,T) > rhoTa));
            }
        }
        
        etasLowerBounds.push_back(std::max(1 - static_cast<float>(sumOfIndicators)/(static_cast<float>(nTbs*(B1+B2))),0.5f));
        cout<< "New Lower Bound Added "<< counterEtas <<" \n"<< std::flush;
        counterEtas++;
    }
    
    {
        pTree Ta = pTree("();");
        
        vector<pTree> Tbs = coverTrees(Ta, compLeafSet);
        
        int nTbs = static_cast<int>(Tbs.size());
        
        int sumOfIndicators = 0;
        
        for (int k = 0; k < treeSample1.size(); k++){
            pTree T = treeSample1.at(k);
            int rhoTa = rho(Ta, T);
            for (pTree Tb : Tbs){
                sumOfIndicators += nSample1.at(k)*((int)(rho(Tb,T) > rhoTa));
            }
        }
        
        for (int k = 0; k < treeSample2.size(); k++){
            pTree T = treeSample2.at(k);
            int rhoTa = rho(Ta,T);
            for (pTree Tb : Tbs){
                sumOfIndicators += nSample2.at(k)*((int)(rho(Tb,T) > rhoTa));
            }
        }
        
        etasLowerBounds.push_back(std::max(1 - static_cast<float>(sumOfIndicators)/(static_cast<float>(nTbs*(B1+B2))),0.5f));
        cout<< "New Lower Bound Added "<< counterEtas <<" \n"<< std::flush;
        counterEtas++;
    }
    
    cout << "About to find FDR controlled trees \n"<< std::flush;
    //vector<pTree> FDRTrees = FDRSearch(treeSample2, nSample2, subPost, etasLowerBounds, q);
    pTree FDRTree = FDRSearchGreedy(treeSample2, nSample2, subPost, etasLowerBounds, q);
    

    //
    // 4. Convert vector<pTree> → CharacterVector
    //
    CharacterVector out(1);
    //for (size_t i = 0; i < FDRTrees.size(); i++) {
    //    mPhylo rP = mPhylo(FDRTrees[i]);
    //    out[i] = rP.toNewick();
    //}
    
    if (FDRTree.rank == 0){
        out[0] = "();";
    } else {
        out[0] = mPhylo(FDRTree).toNewick();
    }

    return out;
    
}

// [[Rcpp::export]]
CharacterVector completeSearchRcppS_V2(CharacterVector treeSample1R,
                                 IntegerVector nSample1R,
                                 CharacterVector treeSample2R,
                                 IntegerVector nSample2R,
                                 CharacterVector compLeafSetR,
                                 int MtR, int rbR, double qR) {
    
    Rcpp::Timer timer;
    cout << "It entered the first Cpp function \n"<< std::flush;
    timer.tic("TreeReading");
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
    
    int B1 = std::accumulate(nSample1.begin(), nSample1.end(), 0);
    int B2 = std::accumulate(nSample2.begin(), nSample2.end(), 0);
    timer.toc("TreeReading");
    //
    // 2. Convert compLeafSetR → set<string>
    //
    timer.tic("LeavesReading");
    std::set<std::string> compLeafSet;
    for (int i = 0; i < compLeafSetR.size(); i++) {
        if (compLeafSetR[i] == NA_STRING)
            stop("compLeafSet cannot contain NA.");
        compLeafSet.insert(as<std::string>(compLeafSetR[i]));
    }
    timer.toc("LeavesReading");
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
    timer.tic("Subposet1");
    subPoset subPost = subPoset(treeSample1, nSample1, compLeafSet, Mt, rb);
    timer.toc("Subposet1");
    
    timer.tic("Subposet2");
    computeAllMaxLevelBounds(subPost);
    timer.toc("Subposet2");
    
    
    //----- Computing eta's lower bounds -------------//
    
    cout << "Computing etasLowerBounds \n"<< std::flush;
    
    // For each node index, store the oRho objects for each sample tree.
    // storedORho1[i][k] = oRho for (Poset[i].Tree, treeSample1[k])
    // storedORho2[i][k] = oRho for (Poset[i].Tree, treeSample2[k])
    timer.tic("Etas");
    int nNodes = static_cast<int>(subPost.Poset.size());
    int nS1 = static_cast<int>(treeSample1.size());
    int nS2 = static_cast<int>(treeSample2.size());
    
    // Base oRho for the empty tree
    pTree emptyTree = pTree("();");
    vector<set<string>> emptyLeaves;
    oRho emptyORho = oRho(0, emptyLeaves);

    vector<vector<oRho>> storedORho1(nNodes, vector<oRho>(nS1, emptyORho));
    vector<vector<oRho>> storedORho2(nNodes, vector<oRho>(nS2, emptyORho));
    
    //Computing Rhos for the trees in Subposet at rank 1;
    int spCurIndx = subPost.firstRank.at(1);
    
    while (spCurIndx > -1) {
        spNode& spN = subPost.Poset[spCurIndx];
        pTree Ta = spN.Tree;
        
        // --- treeSample1 ---
        for (int k = 0; k < nS1; k++){
            pTree T = treeSample1.at(k);
            oRho oRhoTa = rho(Ta, emptyTree, emptyORho, T);

            // Store for use by nodes above Ta
            storedORho1.at(spCurIndx).at(k) = oRhoTa;
        }

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

            // --- treeSample1 ---
            for (int k = 0; k < nS1; k++){
                pTree T = treeSample1.at(k);
                oRho baseORho = storedORho1.at(underIdx).at(k);
                oRho oRhoTa = rho(Ta, Tunder, baseORho, T);
                
                // Store for use by nodes above Ta
                storedORho1.at(spCurIndx).at(k) = oRhoTa;
            }

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
    
    //int nodesCount = 0;
    vector<float> etasLowerBounds;
    for (int idx = 0; idx < nNodes; idx++){
        spNode& spN = subPost.Poset[idx];
        pTree Ta = spN.Tree;

        vector<pTree> Tbs = coverTrees(Ta, compLeafSet);
        int nTbs = max(static_cast<int>(Tbs.size()),1);
        int sumOfIndicators = 0;

        // --- treeSample1 ---
        for (int k = 0; k < nS1; k++){
            pTree T = treeSample1.at(k);

            // Compute rho(Tb, T) for each cover tree, reusing oRhoTa
            for (pTree Tb : Tbs){
                oRho oRhoTb = rho(Tb, Ta, storedORho1.at(idx).at(k), T);
                sumOfIndicators += nSample1.at(k) * ((int)(oRhoTb.rho > storedORho1.at(idx).at(k).rho));
            }
        }

        // --- treeSample2 ---
        for (int k = 0; k < nS2; k++){
            pTree T = treeSample2.at(k);

            for (pTree Tb : Tbs){
                oRho oRhoTb = rho(Tb, Ta, storedORho2.at(idx).at(k), T);
                sumOfIndicators += nSample2.at(k) * ((int)(oRhoTb.rho > storedORho2.at(idx).at(k).rho));
            }
        }

        etasLowerBounds.push_back(std::max(1 - static_cast<float>(sumOfIndicators) /
                                  (static_cast<float>(nTbs * (B1 + B2))), 0.5f));
        
        float tempZeta = (spN.zeta*B1 + static_cast<float>(sumOfIndicators)/static_cast<float>(nTbs))/(static_cast<float>(B1 + B2));
            
        subPost.Poset[idx].setZeta(tempZeta);
        //nodesCount++;
    }
    timer.toc("Etas");
    // --- Empty tree case ---
    /*{
        vector<pTree> Tbs = coverTrees(emptyTree, compLeafSet);
        int nTbs = static_cast<int>(Tbs.size());
        int sumOfIndicators = 0;

        for (int k = 0; k < nS1; k++){
            pTree T = treeSample1.at(k);
            for (pTree Tb : Tbs){
                oRho oRhoTb = rho(Tb, emptyTree, emptyORho, T);
                sumOfIndicators += nSample1.at(k) * ((int)(oRhoTb.rho > 0));
            }
        }

        for (int k = 0; k < nS2; k++){
            pTree T = treeSample2.at(k);
            for (pTree Tb : Tbs){
                oRho oRhoTb = rho(Tb, emptyTree, emptyORho, T);
                sumOfIndicators += nSample2.at(k) * ((int)(oRhoTb.rho > 0));
            }
        }

        etasLowerBounds.push_back(std::max(1 - static_cast<float>(sumOfIndicators) /
                                  (static_cast<float>(nTbs * (B1 + B2))), 0.5f));
        cout << "Eta Lower Bound number " << nodesCount << " computed \n" << std::flush;
        nodesCount++;
    }*/
    
    cout << "About to find FDR controlled trees \n"<< std::flush;
    //vector<pTree> FDRTrees = FDRSearch(treeSample2, nSample2, subPost, etasLowerBounds, q);
    //pTree FDRTree = FDRSearchGreedy(treeSample2, nSample2, storedORho2, subPost, etasLowerBounds, q);
    timer.tic("FDRsearch");
    pTree FDRTree = FDRSearchGreedy(treeSample2, nSample2, storedORho2, subPost, q);
    timer.toc("FDRsearch");
    //
    // 4. Convert vector<pTree> → CharacterVector
    //
    CharacterVector out(1);
    //for (size_t i = 0; i < FDRTrees.size(); i++) {
    //    mPhylo rP = mPhylo(FDRTrees[i]);
    //    out[i] = rP.toNewick();
    //}
    
    if (FDRTree.rank == 0){
        out[0] = "();";
    } else {
        out[0] = mPhylo(FDRTree).toNewick();
    }

    return out;
    
}

// [[Rcpp::export]]
CharacterVector completeSearchRcppS_V3(CharacterVector treeSample1R,
                                 IntegerVector nSample1R,
                                 CharacterVector treeSample2R,
                                 IntegerVector nSample2R,
                                 CharacterVector compLeafSetR,
                                 double qR, double qoR) {
    
    Rcpp::Timer timer;
    
    cout << "It entered the first Cpp function \n"<< std::flush;
    timer.tic("TreeReading");
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
    
    int B1 = std::accumulate(nSample1.begin(), nSample1.end(), 0);
    int B2 = std::accumulate(nSample2.begin(), nSample2.end(), 0);
    timer.toc("TreeReading");

    //
    // 2. Convert compLeafSetR → set<string>
    //
    timer.tic("LeavesReading");
    std::set<std::string> compLeafSet;
    for (int i = 0; i < compLeafSetR.size(); i++) {
        if (compLeafSetR[i] == NA_STRING)
            stop("compLeafSet cannot contain NA.");
        compLeafSet.insert(as<std::string>(compLeafSetR[i]));
    }
    timer.toc("LeavesReading");
    
    //
    // 2.1 Converting the remaining floats
    //
    
    float q = static_cast<float>(qR);
    
    float qo = static_cast<float>(qoR);

    //
    // 3. Call C++ function
    //
    
    cout<< "About to build SubPoset \n"<< std::flush;
    
    timer.tic("Subposet1");
    subPoset subPost = subPoset(treeSample1, nSample1, compLeafSet, qo);
    timer.toc("Subposet1");
    
    timer.tic("Subposet2");
    computeAllMaxLevelBounds(subPost);
    timer.toc("Subposet2");
    //subPost.print();
    
    
    //----- Computing eta's lower bounds -------------//
    
    cout << "Computing etasLowerBounds \n"<< std::flush;
    
    // For each node index, store the oRho objects for each sample tree.
    // storedORho1[i][k] = oRho for (Poset[i].Tree, treeSample1[k])
    // storedORho2[i][k] = oRho for (Poset[i].Tree, treeSample2[k])
    timer.tic("Etas");
    int nNodes = static_cast<int>(subPost.Poset.size());
    int nS2 = static_cast<int>(treeSample2.size());
    // RNG shared across all shuffles
    auto rd  = random_device{};
    auto rng = default_random_engine{rd()};
    int ZETA_SAMPLE = 30;
    
    // Base oRho for the empty tree
    pTree emptyTree = pTree("();");
    vector<set<string>> emptyLeaves;
    oRho emptyORho = oRho(0, emptyLeaves);

    vector<vector<oRho>> storedORho2(nNodes, vector<oRho>(nS2, emptyORho));
    
    //Computing Rhos for the trees in Subposet at rank 1;
    int spCurIndx = subPost.firstRank.at(1);
    
    //int countNodes = 0;
    while (spCurIndx > -1) {
        spNode& spN = subPost.Poset[spCurIndx];
        pTree Ta = spN.Tree;

        // --- treeSample2 ---
        for (int k = 0; k < nS2; k++){
            pTree T = treeSample2.at(k);

            oRho oRhoTa = rho(Ta, emptyTree, emptyORho, T);
            
            storedORho2.at(spCurIndx).at(k) = oRhoTa;
        }
        //countNodes++;
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
            //countNodes++;
            spCurIndx = subPost.Poset.at(spCurIndx).next;
        }
    }
    
    //cout << " The total nodes in the poset is " << nNodes << "\n" << flush;
    //cout << " The values ORho computed were " << countNodes << "\n" << flush;
    
    /*for (int idx = 0; idx < nNodes; idx++){
        spNode& spN = subPost.Poset[idx];
        pTree Ta = spN.Tree;
        
        //cout << "The tree in this node " << idx << "is " << mPhylo(Ta).toNewick() << "\n" << flush;
        
        for (int k = 0; k < nS2; k++){
            cout << "   The rho is " << storedORho2.at(idx).at(k).rho << "\n" << flush;
        }
    }*/
    
    
    //int nodesCount = 0;
    //vector<float> etasLowerBounds;
    for (int idx = 0; idx < nNodes; idx++){
        spNode& spN = subPost.Poset[idx];
        pTree Ta = spN.Tree;

        vector<pTree> tbCovers = coverTrees(Ta, compLeafSet);
        int nTbs = static_cast<int>(tbCovers.size());
        int nsamp = min(ZETA_SAMPLE, nTbs);
        
        //cout << "   The tbCovers is of size " << nTbs << "\n" << flush;
        //cout << "   The nsamp is " << nsamp << "\n" << flush;
        
        if (nsamp > 0){
            shuffle(tbCovers.begin(), tbCovers.end(), rng);
            int sumOfIndicators = 0;

            // --- treeSample2 ---
            for (int k = 0; k < nS2; k++){
                pTree T = treeSample2.at(k);

                for (int i = 0; i < nsamp; i++){
                    pTree Tb = tbCovers[i];
                    //cout << "      calling rho for idx=" << idx << " k=" << k << " i=" << i << "with tree " << mPhylo(Tb).toNewick() << "\n" << flush;
                    oRho oRhoTb = rho(Tb, Ta, storedORho2.at(idx).at(k), T);
                    //cout << "      rho returned " << oRhoTb.rho << "\n" << flush;
                    sumOfIndicators += nSample2.at(k) * ((int)(oRhoTb.rho > storedORho2.at(idx).at(k).rho));
                }
            }
            
            //cout << "      Value of zeta before adjusting is" << spN.zeta << "\n" << flush;

            float tempZeta = (spN.zeta*B1 + static_cast<float>(sumOfIndicators)/static_cast<float>(nsamp))/(static_cast<float>(B1 + B2));
            //cout << "      Value of zeta after adjusting is" << tempZeta << "\n" << flush;
            
            subPost.Poset[idx].setZeta(tempZeta);

        }
        
        //cout << "Eta Lower Bound number " << nodesCount << " computed \n" << std::flush;
        //nodesCount++;
    }
    timer.toc("Etas");
    // --- Empty tree case ---
    // For now, we will estimate that zeta for empty tree is 1/3 (check the map on quartets)
    // We could compute this directly by counting how many leaves (quartets?) each tree in the sample has. 
    
    
    cout << "About to find FDR controlled trees \n"<< std::flush;
    //vector<pTree> FDRTrees = FDRSearch(treeSample2, nSample2, subPost, etasLowerBounds, q);
    timer.tic("FDRsearch");
    pTree FDRTree = FDRSearchGreedy(treeSample2, nSample2, storedORho2, subPost, q);
    timer.toc("FDRsearch");
    
    subPost.print();

    //
    // 4. Convert vector<pTree> → CharacterVector
    //
    CharacterVector out(1);
    //for (size_t i = 0; i < FDRTrees.size(); i++) {
    //    mPhylo rP = mPhylo(FDRTrees[i]);
    //    out[i] = rP.toNewick();
    //}

    if (FDRTree.rank == 0){
        out[0] = "();";
    } else {
        out[0] = mPhylo(FDRTree).toNewick();
    }

    return out;

}

// [[Rcpp::export]]
CharacterVector completeSearchRcppS_V4(CharacterVector treeSample1R,
                                 IntegerVector nSample1R,
                                 CharacterVector treeSample2R,
                                 IntegerVector nSample2R,
                                 CharacterVector compLeafSetR,
                                 double qR, int top_widthR, int bottom_widthR,
                                 std::string orientationR) {

    Rcpp::Timer timer;

    cout << "It entered the first Cpp function \n"<< std::flush;
    timer.tic("TreeReading");
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
    timer.toc("TreeReading");

    //
    // 2. Convert compLeafSetR → set<string>
    //
    timer.tic("LeavesReading");
    std::set<std::string> compLeafSet;
    for (int i = 0; i < compLeafSetR.size(); i++) {
        if (compLeafSetR[i] == NA_STRING)
            stop("compLeafSet cannot contain NA.");
        compLeafSet.insert(as<std::string>(compLeafSetR[i]));
    }
    timer.toc("LeavesReading");

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

    timer.tic("Subposet1");
    subPoset subPost = subPoset(treeSample1, nSample1, compLeafSet,
                                top_width, bottom_width, orientationR);
    timer.toc("Subposet1");

    timer.tic("Subposet2");
    computeAllMaxLevelBounds(subPost);
    timer.toc("Subposet2");
    
    //subPost.print();

    //----- Build oRho cache for treeSample2 -------------//
    // The fixed-width constructor already computes each node's zeta exactly
    // (over all of its covers, from treeSample1), so the zeta re-estimation
    // that V3 did here is redundant and dropped. We still build storedORho2 —
    // the per-node oRho cache against treeSample2 — which FDRSearchGreedy needs.
    cout << "Building oRho cache for treeSample2 \n"<< std::flush;

    timer.tic("Rho2Cache");
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
    timer.toc("Rho2Cache");

    cout << "About to enter FDR-controlled tree search \n"<< std::flush;
    //cout << "With the subposet of size: " << subPost.Poset.size() << "\n" << std::flush;
    timer.tic("FDRsearch");
    pTree FDRTree = FDRSearchGreedy(treeSample2, nSample2, storedORho2, subPost, q);
    timer.toc("FDRsearch");

    //
    // 4. Convert vector<pTree> → CharacterVector
    //
    CharacterVector out(1);
    
    //cout << "Tranforming? \n" << std::flush;
    
    if (FDRTree.rank == 0){
        out[0] = "();";
    } else {
        out[0] = mPhylo(FDRTree).toNewick();
    }

    return out;

}
