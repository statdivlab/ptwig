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
Rcpp::List SPAnalysisR(CharacterVector treeStar,
                                  CharacterVector treeSample1R,
                                 CharacterVector treeSample2R,
                                 CharacterVector bigTreeSampleR,
                                 CharacterVector compLeafSetR,
                                 double alphaR, double qR, double tauR, double deltaR) {
    
    pTree tStar = pTree(as<std::string>(treeStar));
    
    int B1 = treeSample1R.size();
    
    std::vector<pTree> treeSample1;
    treeSample1.reserve(treeSample1R.size());

    for (int i = 0; i < treeSample1R.size(); i++) {
        if (treeSample1R[i] == NA_STRING)
            stop("treeSample cannot contain NA.");
        treeSample1.emplace_back(
            pTree(as<std::string>(treeSample1R[i]))
        );
    }
    
    int B2 = treeSample2R.size();
    
    std::vector<pTree> treeSample2;
    treeSample2.reserve(B2);

    for (int i = 0; i < B2; i++) {
        if (treeSample2R[i] == NA_STRING)
            stop("treeSample cannot contain NA.");
        treeSample2.emplace_back(
            pTree(as<std::string>(treeSample2R[i]))
        );
    }
    
    
    std::vector<pTree> bigTreeSample;
    bigTreeSample.reserve(bigTreeSampleR.size());

    for (int i = 0; i < bigTreeSampleR.size(); i++) {
        if (bigTreeSampleR[i] == NA_STRING)
            stop("treeSample cannot contain NA.");
        bigTreeSample.emplace_back(
            pTree(as<std::string>(bigTreeSampleR[i]))
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
    float tau = static_cast<float>(tauR);
    float delta = static_cast<float>(deltaR);

    //
    // 3. Call C++ function
    //
    
    std::vector<pTree> stTrees = stableSearch(treeSample1, compLeafSet, alpha);
    
    
    int stM = static_cast<int>(stTrees.size());
    int rank_anchor = floor((log(q) - log(stM) + B2*(tau-0.5)*(tau-0.5))/log(2));
    
    cout<< "Stable trees number is " << stM << "\n";
    cout<< "And B2 is " << B2 << "\n";
    cout<< "Just to be sure: q = "<< q << " and tau =" << tau << "\n"; 
    cout<< "The rank anchor is : " << rank_anchor << "\n";
    
    subPoset subPost = subPoset(stTrees, treeSample1, compLeafSet, rank_anchor);
    
    //----- Computing eta's lower bounds -------------//
    
    vector<float> etasLowerBounds;
    
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
        
        etasLowerBounds.push_back(1 - static_cast<float>(sumOfIndicators)/(static_cast<float>(nTbs*(B1+B2))));
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
        
        etasLowerBounds.push_back(1 - static_cast<float>(sumOfIndicators)/(static_cast<float>(nTbs*(B1+B2))));
    }
    
    //------------------------------------------------//
    
    vector<string> subPosetTrees;
    vector<int> subPosetRanks;
    vector<int> subPosetKappas;
    
    for (int k = 0; k < subPost.Poset.size(); k++){
        mPhylo rP = mPhylo(subPost.Poset.at(k).Tree);
        subPosetTrees.push_back(rP.toNewick());
        subPosetRanks.push_back(subPost.Poset.at(k).Tree.rank);
        subPosetKappas.push_back(subPost.Poset.at(k).kappa);
    }
    
    subPosetOutput allResults = SPanalisys(tStar, treeSample2, bigTreeSample, subPost, q, delta);
    
    // Convert edges: split pairs into two parallel integer vectors
    int nEdges = allResults.edges.size();
    Rcpp::IntegerVector edgeFrom(nEdges), edgeTo(nEdges);
    for (int i = 0; i < nEdges; i++) {
        edgeFrom[i] = allResults.edges[i].first;
        edgeTo[i]   = allResults.edges[i].second;
    }
    
    return Rcpp::List::create(
        Rcpp::Named("subPosetTrees")       = Rcpp::wrap(subPosetTrees),
        Rcpp::Named("subPosetRanks")       = Rcpp::wrap(subPosetRanks),
        Rcpp::Named("subPosetKappas")       = Rcpp::wrap(subPosetKappas),
        Rcpp::Named("EtasLowerBounds")  = Rcpp::wrap(etasLowerBounds),
        Rcpp::Named("CoveringPairsLower")  = edgeFrom,
        Rcpp::Named("CoveringPairsUpper")  = edgeTo,
        Rcpp::Named("nullCovering")     = Rcpp::wrap(allResults.nullCovering),
        Rcpp::Named("etaValues")        = Rcpp::wrap(allResults.etaValues),
        Rcpp::Named("coveringMeans")    = Rcpp::wrap(allResults.coveringMean),
        Rcpp::Named("coveringVariance") = Rcpp::wrap(allResults.coveringVariance),
        Rcpp::Named("nullCoveringProb") = allResults.nullCoveringProb,
        Rcpp::Named("minLower")         = allResults.minLower,
        Rcpp::Named("minUpper")         = allResults.minUpper,
        Rcpp::Named("RademacherComplexity") = allResults.RademacherComplex,
        Rcpp::Named("RademacherComplexity2") = allResults.RademacherComplex2,
        Rcpp::Named("kappaThresholds05") = Rcpp::wrap(allResults.kappa_Ts_05),
        Rcpp::Named("kappaThresholdsP") = Rcpp::wrap(allResults.kappa_Ts_p),
        Rcpp::Named("radThresholds05") = Rcpp::wrap(allResults.rad_Ts_05),
        Rcpp::Named("radThresholdsP") = Rcpp::wrap(allResults.rad_Ts_p)
    );
    
}

// [[Rcpp::export]]
Rcpp::List SPAnalysisR2(CharacterVector treeStar,
                                  CharacterVector treeSample1R,
                                 CharacterVector treeSample2R,
                                 CharacterVector bigTreeSampleR,
                                 CharacterVector compLeafSetR,
                                 int MtR, int rbR, double qR, double deltaR) {
    
    pTree tStar = pTree(as<std::string>(treeStar));
    
    std::vector<pTree> treeSample1;
    treeSample1.reserve(treeSample1R.size());
    
    int B1 = treeSample1R.size();

    for (int i = 0; i < treeSample1R.size(); i++) {
        if (treeSample1R[i] == NA_STRING)
            stop("treeSample cannot contain NA.");
        treeSample1.emplace_back(
            pTree(as<std::string>(treeSample1R[i]))
        );
    }
    
    int B2 = treeSample2R.size();
    
    std::vector<pTree> treeSample2;
    treeSample2.reserve(B2);

    for (int i = 0; i < B2; i++) {
        if (treeSample2R[i] == NA_STRING)
            stop("treeSample cannot contain NA.");
        treeSample2.emplace_back(
            pTree(as<std::string>(treeSample2R[i]))
        );
    }
    
    
    std::vector<pTree> bigTreeSample;
    bigTreeSample.reserve(bigTreeSampleR.size());

    for (int i = 0; i < bigTreeSampleR.size(); i++) {
        if (bigTreeSampleR[i] == NA_STRING)
            stop("treeSample cannot contain NA.");
        bigTreeSample.emplace_back(
            pTree(as<std::string>(bigTreeSampleR[i]))
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
    // 2.1 Converting the integers
    //
    
    int Mt = static_cast<int>(MtR);
    
    int rb = static_cast<int>(rbR);
    
    float q = static_cast<float>(qR);
    
    float delta = static_cast<float>(deltaR);

    //
    // 3. Call C++ function
    //
    
    subPoset subPost = subPoset(treeSample1, compLeafSet, Mt, rb);
    
    //----- Computing eta's lower bounds -------------//
    
    vector<float> etasLowerBounds;
    
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
        
        etasLowerBounds.push_back(1 - static_cast<float>(sumOfIndicators)/(static_cast<float>(nTbs*(B1+B2))));
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
        
        etasLowerBounds.push_back(1 - static_cast<float>(sumOfIndicators)/(static_cast<float>(nTbs*(B1+B2))));
    }
    
    //------------------------------------------------//
    
    vector<string> subPosetTrees;
    vector<int> subPosetRanks;
    vector<int> subPosetKappas;
    
    for (int k = 0; k < subPost.Poset.size(); k++){
        mPhylo rP = mPhylo(subPost.Poset.at(k).Tree);
        subPosetTrees.push_back(rP.toNewick());
        subPosetRanks.push_back(subPost.Poset.at(k).Tree.rank);
        subPosetKappas.push_back(subPost.Poset.at(k).kappa);
    }
    
    subPosetOutput allResults = SPanalisys(tStar, treeSample2, bigTreeSample, subPost, q, delta);
    
    // Convert edges: split pairs into two parallel integer vectors
    int nEdges = allResults.edges.size();
    Rcpp::IntegerVector edgeFrom(nEdges), edgeTo(nEdges);
    for (int i = 0; i < nEdges; i++) {
        edgeFrom[i] = allResults.edges[i].first;
        edgeTo[i]   = allResults.edges[i].second;
    }
    
    return Rcpp::List::create(
        Rcpp::Named("subPosetTrees")       = Rcpp::wrap(subPosetTrees),
        Rcpp::Named("subPosetRanks")       = Rcpp::wrap(subPosetRanks),
        Rcpp::Named("subPosetKappas")       = Rcpp::wrap(subPosetKappas),
        Rcpp::Named("EtasLowerBounds")  = Rcpp::wrap(etasLowerBounds),
        Rcpp::Named("CoveringPairsLower")  = edgeFrom,
        Rcpp::Named("CoveringPairsUpper")  = edgeTo,
        Rcpp::Named("nullCovering")     = Rcpp::wrap(allResults.nullCovering),
        Rcpp::Named("etaValues")        = Rcpp::wrap(allResults.etaValues),
        Rcpp::Named("coveringMeans")    = Rcpp::wrap(allResults.coveringMean),
        Rcpp::Named("coveringVariance") = Rcpp::wrap(allResults.coveringVariance),
        Rcpp::Named("nullCoveringProb") = allResults.nullCoveringProb,
        Rcpp::Named("minLower")         = allResults.minLower,
        Rcpp::Named("minUpper")         = allResults.minUpper,
        Rcpp::Named("RademacherComplexity") = allResults.RademacherComplex,
        Rcpp::Named("RademacherComplexity2") = allResults.RademacherComplex2,
        Rcpp::Named("kappaThresholds05") = Rcpp::wrap(allResults.kappa_Ts_05),
        Rcpp::Named("kappaThresholdsP") = Rcpp::wrap(allResults.kappa_Ts_p),
        Rcpp::Named("radThresholds05") = Rcpp::wrap(allResults.rad_Ts_05),
        Rcpp::Named("radThresholdsP") = Rcpp::wrap(allResults.rad_Ts_p)
    );
}


// [[Rcpp::export]]
Rcpp::List SPAnalysisRS(CharacterVector treeStar,
                                 CharacterVector treeSample1R,
                                 IntegerVector nSample1R,
                                 CharacterVector treeSample2R,
                                 IntegerVector nSample2R,
                                 CharacterVector bigTreeSampleR,
                                 IntegerVector nBSampleR,
                                 CharacterVector compLeafSetR,
                                 double alphaR, double qR, double tauR, double deltaR) {
    
    pTree tStar = pTree(as<std::string>(treeStar));
    
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
    
    std::vector<pTree> bigTreeSample;
    bigTreeSample.reserve(bigTreeSampleR.size());

    for (int i = 0; i < bigTreeSampleR.size(); i++) {
        if (bigTreeSampleR[i] == NA_STRING)
            stop("treeSample cannot contain NA.");
        bigTreeSample.emplace_back(
            pTree(as<std::string>(bigTreeSampleR[i]))
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
    
    std::vector<int> nBSample;
    
    for (int i = 0; i < nBSampleR.size(); i++){
        nBSample.push_back(static_cast<int>(nBSampleR[i]));
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
    
    float alpha = static_cast<float>(alphaR);
    float q = static_cast<float>(qR);
    float tau = static_cast<float>(tauR);
    float delta = static_cast<float>(deltaR);

    //
    // 3. Call C++ function
    //
    
    std::vector<pTree> stTrees = stableSearch(treeSample1, nSample1, compLeafSet, alpha);
    
    
    int stM = static_cast<int>(stTrees.size());
    int rank_anchor = floor((log(q) - log(stM) + B2*(tau-0.5)*(tau-0.5))/log(2));
    
    cout<< "Stable trees number is " << stM << "\n";
    cout<< "And B2 is " << B2 << "\n";
    cout<< "Just to be sure: q = "<< q << " and tau =" << tau << "\n"; 
    cout<< "The rank anchor is : " << rank_anchor << "\n";
    
    subPoset subPost = subPoset(stTrees, treeSample1, nSample1, compLeafSet, rank_anchor);
    
    //----- Computing eta's lower bounds -------------//
    
    vector<float> etasLowerBounds;
    
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
        
        etasLowerBounds.push_back(1 - static_cast<float>(sumOfIndicators)/(static_cast<float>(nTbs*(B1+B2))));
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
        
        etasLowerBounds.push_back(1 - static_cast<float>(sumOfIndicators)/(static_cast<float>(nTbs*(B1+B2))));
    }
    
    //------------------------------------------------//
    
    vector<string> subPosetTrees;
    vector<int> subPosetRanks;
    vector<int> subPosetKappas;
    
    for (int k = 0; k < subPost.Poset.size(); k++){
        mPhylo rP = mPhylo(subPost.Poset.at(k).Tree);
        subPosetTrees.push_back(rP.toNewick());
        subPosetRanks.push_back(subPost.Poset.at(k).Tree.rank);
        subPosetKappas.push_back(subPost.Poset.at(k).kappa);
    }
    
    subPosetOutput allResults = SPanalisys(tStar, treeSample2, nSample2, bigTreeSample, nBSample, subPost,  q, delta);
    
    // Convert edges: split pairs into two parallel integer vectors
    int nEdges = allResults.edges.size();
    Rcpp::IntegerVector edgeFrom(nEdges), edgeTo(nEdges);
    for (int i = 0; i < nEdges; i++) {
        edgeFrom[i] = allResults.edges[i].first;
        edgeTo[i]   = allResults.edges[i].second;
    }
    
    return Rcpp::List::create(
        Rcpp::Named("subPosetTrees")       = Rcpp::wrap(subPosetTrees),
        Rcpp::Named("subPosetRanks")       = Rcpp::wrap(subPosetRanks),
        Rcpp::Named("subPosetKappas")       = Rcpp::wrap(subPosetKappas),
        Rcpp::Named("EtasLowerBounds")  = Rcpp::wrap(etasLowerBounds),
        Rcpp::Named("CoveringPairsLower")  = edgeFrom,
        Rcpp::Named("CoveringPairsUpper")  = edgeTo,
        Rcpp::Named("nullCovering")     = Rcpp::wrap(allResults.nullCovering),
        Rcpp::Named("etaValues")        = Rcpp::wrap(allResults.etaValues),
        Rcpp::Named("coveringMeans")    = Rcpp::wrap(allResults.coveringMean),
        Rcpp::Named("coveringVariance") = Rcpp::wrap(allResults.coveringVariance),
        Rcpp::Named("nullCoveringProb") = allResults.nullCoveringProb,
        Rcpp::Named("minLower")         = allResults.minLower,
        Rcpp::Named("minUpper")         = allResults.minUpper,
        Rcpp::Named("RademacherComplexity") = allResults.RademacherComplex,
        Rcpp::Named("RademacherComplexity2") = allResults.RademacherComplex2,
        Rcpp::Named("kappaThresholds05") = Rcpp::wrap(allResults.kappa_Ts_05),
        Rcpp::Named("kappaThresholdsP") = Rcpp::wrap(allResults.kappa_Ts_p),
        Rcpp::Named("radThresholds05") = Rcpp::wrap(allResults.rad_Ts_05),
        Rcpp::Named("radThresholdsP") = Rcpp::wrap(allResults.rad_Ts_p)
    );
    
}

// [[Rcpp::export]]
Rcpp::List SPAnalysisR2S(CharacterVector treeStar,
                                 CharacterVector treeSample1R,
                                 IntegerVector nSample1R,
                                 CharacterVector treeSample2R,
                                 IntegerVector nSample2R,
                                 CharacterVector bigTreeSampleR,
                                 IntegerVector nBSampleR,
                                 CharacterVector compLeafSetR,
                                 int MtR, int rbR, double qR, double deltaR) {
    
    pTree tStar = pTree(as<std::string>(treeStar));
    
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
    
    std::vector<pTree> bigTreeSample;
    bigTreeSample.reserve(bigTreeSampleR.size());

    for (int i = 0; i < bigTreeSampleR.size(); i++) {
        if (bigTreeSampleR[i] == NA_STRING)
            stop("treeSample cannot contain NA.");
        bigTreeSample.emplace_back(
            pTree(as<std::string>(bigTreeSampleR[i]))
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
    
    std::vector<int> nBSample;
    
    for (int i = 0; i < nBSampleR.size(); i++){
        nBSample.push_back(static_cast<int>(nBSampleR[i]));
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
    // 2.1 Converting the integers
    //
    
    int Mt = static_cast<int>(MtR);
    
    int rb = static_cast<int>(rbR);
    
    float q = static_cast<float>(qR);
    
    float delta = static_cast<float>(deltaR);

    //
    // 3. Call C++ function
    //
    
    subPoset subPost = subPoset(treeSample1, nSample1, compLeafSet, Mt, rb);
    
    
    //----- Computing eta's lower bounds -------------//
    
    vector<float> etasLowerBounds;
    
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
        
        etasLowerBounds.push_back(1 - static_cast<float>(sumOfIndicators)/(static_cast<float>(nTbs*(B1+B2))));
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
        
        etasLowerBounds.push_back(1 - static_cast<float>(sumOfIndicators)/(static_cast<float>(nTbs*(B1+B2))));
    }
    
    //------------------------------------------------//
    
    vector<string> subPosetTrees;
    vector<int> subPosetRanks;
    vector<int> subPosetKappas;
    
    for (int k = 0; k < subPost.Poset.size(); k++){
        mPhylo rP = mPhylo(subPost.Poset.at(k).Tree);
        subPosetTrees.push_back(rP.toNewick());
        subPosetRanks.push_back(subPost.Poset.at(k).Tree.rank);
        subPosetKappas.push_back(subPost.Poset.at(k).kappa);
    }
    
    subPosetOutput allResults = SPanalisys(tStar, treeSample2, nSample2, bigTreeSample, nBSample, subPost,  q, delta);
    
    // Convert edges: split pairs into two parallel integer vectors
    int nEdges = allResults.edges.size();
    Rcpp::IntegerVector edgeFrom(nEdges), edgeTo(nEdges);
    for (int i = 0; i < nEdges; i++) {
        edgeFrom[i] = allResults.edges[i].first;
        edgeTo[i]   = allResults.edges[i].second;
    }
    
    return Rcpp::List::create(
        Rcpp::Named("subPosetTrees")       = Rcpp::wrap(subPosetTrees),
        Rcpp::Named("subPosetRanks")       = Rcpp::wrap(subPosetRanks),
        Rcpp::Named("subPosetKappas")       = Rcpp::wrap(subPosetKappas),
        Rcpp::Named("EtasLowerBounds")  = Rcpp::wrap(etasLowerBounds),
        Rcpp::Named("CoveringPairsLower")  = edgeFrom,
        Rcpp::Named("CoveringPairsUpper")  = edgeTo,
        Rcpp::Named("nullCovering")     = Rcpp::wrap(allResults.nullCovering),
        Rcpp::Named("etaValues")        = Rcpp::wrap(allResults.etaValues),
        Rcpp::Named("coveringMeans")    = Rcpp::wrap(allResults.coveringMean),
        Rcpp::Named("coveringVariance") = Rcpp::wrap(allResults.coveringVariance),
        Rcpp::Named("nullCoveringProb") = allResults.nullCoveringProb,
        Rcpp::Named("minLower")         = allResults.minLower,
        Rcpp::Named("minUpper")         = allResults.minUpper,
        Rcpp::Named("RademacherComplexity") = allResults.RademacherComplex,
        Rcpp::Named("RademacherComplexity2") = allResults.RademacherComplex2,
        Rcpp::Named("kappaThresholds05") = Rcpp::wrap(allResults.kappa_Ts_05),
        Rcpp::Named("kappaThresholdsP") = Rcpp::wrap(allResults.kappa_Ts_p),
        Rcpp::Named("radThresholds05") = Rcpp::wrap(allResults.rad_Ts_05),
        Rcpp::Named("radThresholdsP") = Rcpp::wrap(allResults.rad_Ts_p)
    );
    
}

// [[Rcpp::export]]
Rcpp::List simpleSPAnalysis(CharacterVector treeSample1R,
                                 IntegerVector nSample1R,
                                 CharacterVector compLeafSetR,
                                 int MtR, int rbR, double qR, double deltaR) {
    
    std::vector<pTree> treeSample1;
    treeSample1.reserve(treeSample1R.size());

    for (int i = 0; i < treeSample1R.size(); i++) {
        if (treeSample1R[i] == NA_STRING)
            stop("treeSample cannot contain NA.");
        treeSample1.emplace_back(
            pTree(as<std::string>(treeSample1R[i]))
        );
    }
    
    std::vector<int> nSample1;
    
    for (int i = 0; i < nSample1R.size(); i++){
        nSample1.push_back(static_cast<int>(nSample1R[i]));
    }
    
    int B1 = std::accumulate(nSample1.begin(), nSample1.end(), 0);
    
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
    // 2.1 Converting the integers
    //
    
    int Mt = static_cast<int>(MtR);
    
    int rb = static_cast<int>(rbR);
    
    float q = static_cast<float>(qR);
    
    float delta = static_cast<float>(deltaR);

    //
    // 3. Call C++ function
    //
    
    subPoset subPost = subPoset(treeSample1, nSample1, compLeafSet, 0.5f);
    
    cout << "Computing 1 \n" << std::flush;
    computeAllChainCounts(subPost);
    cout << "Computing 2 \n" << std::flush;
    computeAllMaxLevelBounds(subPost);
    
    cout << "Computing 3 \n" << std::flush;
    
    //------------------------------------------------//
    
    vector<string> subPosetTrees;
    vector<int> subPosetRanks;
    vector<int> subPosetKappas;
    
    for (int k = 0; k < subPost.Poset.size(); k++){
        mPhylo rP = mPhylo(subPost.Poset.at(k).Tree);
        subPosetTrees.push_back(rP.toNewick());
        subPosetRanks.push_back(subPost.Poset.at(k).Tree.rank);
        subPosetKappas.push_back(subPost.Poset.at(k).kappa);
    }
    
    
    // Convert edges: split pairs into two parallel integer vectors
    std::vector<std::pair<int,int>> edges;
    std::vector<int> antiChainEst;
    std::vector<int64_t> antiChainUBound;
    
    { int r = 1;
      int v = subPost.firstRank[r];
      while (v != -1){
          edges.push_back({0,(v+1)});
          antiChainEst.push_back(subPost.Poset[v].boundAntichain[0]);
          antiChainUBound.push_back(subPost.Poset[v].chainCountIe[0]);
          v = subPost.Poset[v].next;
      }
      
    }
    
    for (int r = 2; r < (int)subPost.firstRank.size(); ++r) {
        int v = subPost.firstRank[r];
        while (v != -1){
            for (int i = 0; i < subPost.Poset[v].under.size(); ++i) {
                int u = subPost.Poset[v].under[i];
                edges.push_back({(u+1),(v+1)});
                antiChainEst.push_back(subPost.Poset[v].boundAntichain[i]);
                antiChainUBound.push_back(subPost.Poset[v].chainCountIe[i]);
            }
            v = subPost.Poset[v].next;
        }
    }
    
    int nEdges = edges.size();
    
    Rcpp::IntegerVector edgeFrom(nEdges), edgeTo(nEdges);
    
    for (int i = 0; i < nEdges; i++) {
        edgeFrom[i] = edges[i].first;
        edgeTo[i]   = edges[i].second;
    }
    
    return Rcpp::List::create(
        Rcpp::Named("subPosetTrees")       = Rcpp::wrap(subPosetTrees),
        Rcpp::Named("subPosetRanks")       = Rcpp::wrap(subPosetRanks),
        Rcpp::Named("subPosetKappas")       = Rcpp::wrap(subPosetKappas),
        Rcpp::Named("subPosetNus_HighBound")     = Rcpp::wrap(antiChainUBound),
        Rcpp::Named("subPosetNus_Est")       = Rcpp::wrap(antiChainEst),
        Rcpp::Named("CoveringPairsLower")  = edgeFrom,
        Rcpp::Named("CoveringPairsUpper")  = edgeTo
    );
    
}

