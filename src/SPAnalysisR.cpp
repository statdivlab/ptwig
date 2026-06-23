#include "pTree.h"
#include "mPhylo.h"
#include "rho.h"
#include "coverTrees.h"
#include "stableSearch.h"
#include "subPoset.h"
#include "FDRSearch.h"
#include "idNullCoveringPairsComputation.h"
#include "subPosetAnalysis.h"
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
                                 int MtR, int rbR, double qR, double q0R,
                                 int top_width1, int bottom_width1,
                                 int top_width2, int bottom_width2) {
    
    Rcpp::Timer timer("times_main");
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
    
    //float q = static_cast<float>(qR);
    
    float q0 = static_cast<float>(q0R);

    //
    // 3. Call C++ function
    //
    
    timer.tic("UpwardsSubposet");
    subPoset subPost_upwards = subPoset(treeSample1, nSample1, compLeafSet, q0);
    timer.toc("UpwardsSubposet");
    
    
    computeAllMaxLevelBounds(subPost_upwards);
    
    //------------------------------------------------//
    
    vector<string> subPosetTrees_upwards;
    vector<int> subPosetRanks_upwards;
    
    for (int k = 0; k < subPost_upwards.Poset.size(); k++){
        mPhylo rP = mPhylo(subPost_upwards.Poset.at(k).Tree);
        subPosetTrees_upwards.push_back(rP.toNewick());
        subPosetRanks_upwards.push_back(subPost_upwards.Poset.at(k).Tree.rank);
    }
    
    
    // Convert edges: split pairs into two parallel integer vectors
    std::vector<std::pair<int,int>> edges_upwards;
    std::vector<int> antiChainEst_upwards;
    
    { int r = 1;
      int v = subPost_upwards.firstRank[r];
      while (v != -1){
          edges_upwards.push_back({0,(v+1)});
          antiChainEst_upwards.push_back(subPost_upwards.Poset[v].boundAntichain[0]);
          v = subPost_upwards.Poset[v].next;
      }
      
    }
    
    for (int r = 2; r < (int)subPost_upwards.firstRank.size(); ++r) {
        int v = subPost_upwards.firstRank[r];
        while (v != -1){
            for (int i = 0; i < subPost_upwards.Poset[v].under.size(); ++i) {
                int u = subPost_upwards.Poset[v].under[i];
                edges_upwards.push_back({(u+1),(v+1)});
                antiChainEst_upwards.push_back(subPost_upwards.Poset[v].boundAntichain[i]);
            }
            v = subPost_upwards.Poset[v].next;
        }
    }
    
    int nEdges_upwards = edges_upwards.size();
    
    Rcpp::IntegerVector edgeFrom_upwards(nEdges_upwards), edgeTo_upwards(nEdges_upwards);
    
    for (int i = 0; i < nEdges_upwards; i++) {
        edgeFrom_upwards[i] = edges_upwards[i].first;
        edgeTo_upwards[i]   = edges_upwards[i].second;
    }
    
    //
    // 3. Call C++ function
    //
    
    timer.tic("BasicSubposet");
    subPoset subPost_basic = subPoset(treeSample1, nSample1, compLeafSet, Mt, rb);
    timer.toc("BasicSubposet");
    
    
    computeAllMaxLevelBounds(subPost_basic);
    
    //------------------------------------------------//
    
    vector<string> subPosetTrees_basic;
    vector<int> subPosetRanks_basic;
    
    for (int k = 0; k < subPost_basic.Poset.size(); k++){
        mPhylo rP = mPhylo(subPost_basic.Poset.at(k).Tree);
        subPosetTrees_basic.push_back(rP.toNewick());
        subPosetRanks_basic.push_back(subPost_basic.Poset.at(k).Tree.rank);
    }
    
    
    // Convert edges: split pairs into two parallel integer vectors
    std::vector<std::pair<int,int>> edges_basic;
    std::vector<int> antiChainEst_basic;
    
    { int r = 1;
      int v = subPost_basic.firstRank[r];
      while (v != -1){
          edges_basic.push_back({0,(v+1)});
          antiChainEst_basic.push_back(subPost_basic.Poset[v].boundAntichain[0]);
          v = subPost_basic.Poset[v].next;
      }
      
    }
    
    for (int r = 2; r < (int)subPost_basic.firstRank.size(); ++r) {
        int v = subPost_basic.firstRank[r];
        while (v != -1){
            for (int i = 0; i < subPost_basic.Poset[v].under.size(); ++i) {
                int u = subPost_basic.Poset[v].under[i];
                edges_basic.push_back({(u+1),(v+1)});
                antiChainEst_basic.push_back(subPost_basic.Poset[v].boundAntichain[i]);
            }
            v = subPost_basic.Poset[v].next;
        }
    }
    
    int nEdges_basic = edges_basic.size();

    Rcpp::IntegerVector edgeFrom_basic(nEdges_basic), edgeTo_basic(nEdges_basic);

    for (int i = 0; i < nEdges_basic; i++) {
        edgeFrom_basic[i] = edges_basic[i].first;
        edgeTo_basic[i]   = edges_basic[i].second;
    }

    // ------------------------------------------------------------------
    // New fixed-width constructor, upwards orientation (top_width1, bottom_width1)
    // ------------------------------------------------------------------
    // This builder uses the same implicit-empty-bottom convention as the
    // upwards/basic builders (the empty tree is node 0, not stored in Poset),
    // so trees/edges are extracted exactly like the _U / _B blocks: rank-1
    // nodes connect to bottom 0, real nodes are labelled 1..n.
    timer.tic("FixedWidthUpSubposet");
    subPoset subPost_fwUp = subPoset(treeSample1, nSample1, compLeafSet,
                                     top_width1, bottom_width1,
                                     std::string("upwards"));
    timer.toc("FixedWidthUpSubposet");

    computeAllMaxLevelBounds(subPost_fwUp);

    vector<string> subPosetTrees_fwUp;
    vector<int>    subPosetRanks_fwUp;

    for (int k = 0; k < subPost_fwUp.Poset.size(); k++){
        mPhylo rP = mPhylo(subPost_fwUp.Poset.at(k).Tree);
        subPosetTrees_fwUp.push_back(rP.toNewick());
        subPosetRanks_fwUp.push_back(subPost_fwUp.Poset.at(k).Tree.rank);
    }

    std::vector<std::pair<int,int>> edges_fwUp;
    std::vector<int> antiChainEst_fwUp;

    { int r = 1;
      int v = subPost_fwUp.firstRank[r];
      while (v != -1){
          edges_fwUp.push_back({0,(v+1)});
          antiChainEst_fwUp.push_back(subPost_fwUp.Poset[v].boundAntichain[0]);
          v = subPost_fwUp.Poset[v].next;
      }
    }

    for (int r = 2; r < (int)subPost_fwUp.firstRank.size(); ++r) {
        int v = subPost_fwUp.firstRank[r];
        while (v != -1){
            for (int i = 0; i < subPost_fwUp.Poset[v].under.size(); ++i) {
                int u = subPost_fwUp.Poset[v].under[i];
                edges_fwUp.push_back({(u+1),(v+1)});
                antiChainEst_fwUp.push_back(subPost_fwUp.Poset[v].boundAntichain[i]);
            }
            v = subPost_fwUp.Poset[v].next;
        }
    }

    int nEdges_fwUp = edges_fwUp.size();
    Rcpp::IntegerVector edgeFrom_fwUp(nEdges_fwUp), edgeTo_fwUp(nEdges_fwUp);
    for (int i = 0; i < nEdges_fwUp; i++) {
        edgeFrom_fwUp[i] = edges_fwUp[i].first;
        edgeTo_fwUp[i]   = edges_fwUp[i].second;
    }

    // ------------------------------------------------------------------
    // New fixed-width constructor, downwards orientation (top_width2, bottom_width2)
    // ------------------------------------------------------------------
    timer.tic("FixedWidthDownSubposet");
    subPoset subPost_fwDn = subPoset(treeSample1, nSample1, compLeafSet,
                                     top_width2, bottom_width2,
                                     std::string("downwards"));
    timer.toc("FixedWidthDownSubposet");

    computeAllMaxLevelBounds(subPost_fwDn);

    vector<string> subPosetTrees_fwDn;
    vector<int>    subPosetRanks_fwDn;

    for (int k = 0; k < subPost_fwDn.Poset.size(); k++){
        mPhylo rP = mPhylo(subPost_fwDn.Poset.at(k).Tree);
        subPosetTrees_fwDn.push_back(rP.toNewick());
        subPosetRanks_fwDn.push_back(subPost_fwDn.Poset.at(k).Tree.rank);
    }

    std::vector<std::pair<int,int>> edges_fwDn;
    std::vector<int> antiChainEst_fwDn;

    { int r = 1;
      int v = subPost_fwDn.firstRank[r];
      while (v != -1){
          edges_fwDn.push_back({0,(v+1)});
          antiChainEst_fwDn.push_back(subPost_fwDn.Poset[v].boundAntichain[0]);
          v = subPost_fwDn.Poset[v].next;
      }
    }

    for (int r = 2; r < (int)subPost_fwDn.firstRank.size(); ++r) {
        int v = subPost_fwDn.firstRank[r];
        while (v != -1){
            for (int i = 0; i < subPost_fwDn.Poset[v].under.size(); ++i) {
                int u = subPost_fwDn.Poset[v].under[i];
                edges_fwDn.push_back({(u+1),(v+1)});
                antiChainEst_fwDn.push_back(subPost_fwDn.Poset[v].boundAntichain[i]);
            }
            v = subPost_fwDn.Poset[v].next;
        }
    }

    int nEdges_fwDn = edges_fwDn.size();
    Rcpp::IntegerVector edgeFrom_fwDn(nEdges_fwDn), edgeTo_fwDn(nEdges_fwDn);
    for (int i = 0; i < nEdges_fwDn; i++) {
        edgeFrom_fwDn[i] = edges_fwDn[i].first;
        edgeTo_fwDn[i]   = edges_fwDn[i].second;
    }

    return Rcpp::List::create(
        Rcpp::Named("subPosetTrees_U")       = Rcpp::wrap(subPosetTrees_upwards),
        Rcpp::Named("subPosetRanks_U")       = Rcpp::wrap(subPosetRanks_upwards),
        Rcpp::Named("subPosetNus_U")       = Rcpp::wrap(antiChainEst_upwards),
        Rcpp::Named("CoveringPairsLower_U")  = edgeFrom_upwards,
        Rcpp::Named("CoveringPairsUpper_U")  = edgeTo_upwards,
        Rcpp::Named("subPosetTrees_B")       = Rcpp::wrap(subPosetTrees_basic),
        Rcpp::Named("subPosetRanks_B")       = Rcpp::wrap(subPosetRanks_basic),
        Rcpp::Named("subPosetNus_B")       = Rcpp::wrap(antiChainEst_basic),
        Rcpp::Named("CoveringPairsLower_B")  = edgeFrom_basic,
        Rcpp::Named("CoveringPairsUpper_B")  = edgeTo_basic,
        Rcpp::Named("subPosetTrees_FU")      = Rcpp::wrap(subPosetTrees_fwUp),
        Rcpp::Named("subPosetRanks_FU")      = Rcpp::wrap(subPosetRanks_fwUp),
        Rcpp::Named("subPosetNus_FU")        = Rcpp::wrap(antiChainEst_fwUp),
        Rcpp::Named("CoveringPairsLower_FU") = edgeFrom_fwUp,
        Rcpp::Named("CoveringPairsUpper_FU") = edgeTo_fwUp,
        Rcpp::Named("subPosetTrees_FD")      = Rcpp::wrap(subPosetTrees_fwDn),
        Rcpp::Named("subPosetRanks_FD")      = Rcpp::wrap(subPosetRanks_fwDn),
        Rcpp::Named("subPosetNus_FD")        = Rcpp::wrap(antiChainEst_fwDn),
        Rcpp::Named("CoveringPairsLower_FD") = edgeFrom_fwDn,
        Rcpp::Named("CoveringPairsUpper_FD") = edgeTo_fwDn
    );

}

