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
        mPhylo rP = mPhylo(FDRTrees[i]);
        out[i] = rP.toNewick();
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
    
    int B2 = treeSample2R.size();

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
    
    out[0] = mPhylo(FDRTree).toNewick();

    return out;
    
}

// [[Rcpp::export]]
CharacterVector completeSearchRcppS_V2(CharacterVector treeSample1R,
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
    
    int B1 = std::accumulate(nSample1.begin(), nSample1.end(), 0);
    int B2 = std::accumulate(nSample2.begin(), nSample2.end(), 0);

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
    
    
    //----- Computing eta's lower bounds -------------//
    
    cout << "Computing etasLowerBounds \n"<< std::flush;
    
    vector<float> etasLowerBounds;
    
    int nodesCount = 0;
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
        
        cout << "Eta Lower Bound number " << nodesCount << " computed \n"<< std::flush;
        nodesCount++; 
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
        cout << "Eta Lower Bound number " << nodesCount << " computed \n"<< std::flush;
        nodesCount++;
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
    
    out[0] = mPhylo(FDRTree).toNewick();

    return out;
    
}
