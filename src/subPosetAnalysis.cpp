#include "pTree.h"
#include "mPhylo.h"
#include "subPoset.h"
#include "rho.h"
#include "coverTrees.h"
#include "subPosetAnalysis.h"
#include <iostream>
#include <set>
#include <queue>
#include <vector>
#include <stack>
#include <fstream>
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
#include <algorithm>
#include <functional>

using namespace Rcpp;
using namespace std;


float rad_threshold(float omegaT, float q, int N, float Rad, float delta){
    float Udelta = 2*Rad + 2*sqrt((log(2) - log(delta))/(2*N));
    float prelim = log(omegaT) - log((q - delta));
    float thrs1 = std::max(2*prelim/N, 0.0f);
    
    return(Udelta + sqrt(thrs1));
}

subPosetOutput SPanalisys(pTree tStar, vector<pTree> treeSample, vector<pTree> bigTreeSample, subPoset SP, float q, float delta){
    
    subPosetOutput result;
    
    //vector<array<int, 2>> NullPairs;
    
    //Finding the pairs in C_null with upper tree of rank 1 and forming the edges above
    
    int curIndx = SP.firstRank.at(1);
    //int numberPairsFirstLevel = 0;
    
    vector<vector<float>> covPairsIndicators;
    //int numPairs = 0;
    while (curIndx > -1){
        result.edges.push_back({-1, (curIndx+1)});
        vector<float> IndividualIndicators;
        if (rho(SP.Poset.at(curIndx).Tree,tStar) == 0){
            //NullPairs.push_back({curIndx,-1});
            result.nullCovering.push_back(true);
            //numberPairsFirstLevel++;
            //numPairs++;
            //cout << "Tree pair " << numPairs << " that belong to null pairs with 0 : \n";
            //SP.Poset.at(curIndx).Tree.print();
            //cout << "\n";
        } else {
            result.nullCovering.push_back(false);
        }
        
        float curmean = (rho(SP.Poset.at(curIndx).Tree, treeSample.at(0)) > 0);
        float curvar = 0;
        int curn = 1;
        
        IndividualIndicators.push_back(curmean);
        
        for (int k = 1; k < treeSample.size(); k++){
            float Xvalue = (rho(SP.Poset.at(curIndx).Tree, treeSample.at(k)) > 0);
            IndividualIndicators.push_back(Xvalue);
            curmean = (curn*curmean + Xvalue)/(1+curn);
            curvar = (curn*curvar)/(1+curn) + (curmean-Xvalue)*(curmean-Xvalue)/curn;
            curn++;
        }
        
        covPairsIndicators.push_back(IndividualIndicators);
        
        result.coveringMean.push_back(curmean);
        
        result.coveringVariance.push_back(curvar);
        
        curIndx = SP.Poset.at(curIndx).next;
    }
    
    //Finding the pairs above
    
    for (int i = 0; i < SP.Poset.size(); i++){
        for (int j : SP.Poset.at(i).over){
            result.edges.push_back({(i+1),(j+1)});
            vector<float> IndividualIndicators;
            if (rho(SP.Poset.at(i).Tree, tStar) == rho(SP.Poset.at(j).Tree, tStar)){
                //NullPairs.push_back({j,i});
                result.nullCovering.push_back(true);
                //numPairs++;
                //cout << "Tree pair " << numPairs << " that belong to null pairs: \n";
                //SP.Poset.at(j).Tree.print();
                //SP.Poset.at(i).Tree.print();
                //cout << "\n";
            } else {
                result.nullCovering.push_back(false);
            }
            
            float curmean = (rho(SP.Poset.at(j).Tree, treeSample.at(0)) > rho(SP.Poset.at(i).Tree, treeSample.at(0)));
            float curvar = 0;
            int curn = 1;
            
            IndividualIndicators.push_back(curmean);
        
            for (int k = 1; k < treeSample.size(); k++){
                float Xvalue = (rho(SP.Poset.at(j).Tree, treeSample.at(k)) > rho(SP.Poset.at(i).Tree, treeSample.at(k)));
                
                IndividualIndicators.push_back(Xvalue);
            
                curmean = (curn*curmean + Xvalue)/(1+curn);
                curvar = (curn*curvar)/(1+curn) + (curmean-Xvalue)*(curmean-Xvalue)/curn;
                curn++;
            }
        
            result.coveringMean.push_back(curmean);
        
            result.coveringVariance.push_back(curvar);
            
            covPairsIndicators.push_back(IndividualIndicators);
        }
    }
    
    // Estimating Rademacher Complexity bound
    
    int B = (int)treeSample.size();
    
    int T = 10000;
    int k = (int)result.edges.size();

    std::mt19937 rng(42);
    std::bernoulli_distribution coin(0.5);

    double total = 0.0;
    std::vector<int> sigma(B);

    for (int t = 0; t < T; ++t) {
        // Draw sigma_i in {-1, +1}
        for (int i = 0; i < B; ++i)
            sigma[i] = coin(rng) ? 1 : -1;

        // Compute (1/m) sum_i sigma_i f_i(a) for each a, take max
        double best = -std::numeric_limits<double>::infinity();
        
        for (int a = 0; a < k; ++a) {
            double val = 0.0;
            int sigmaCounter = 0;
            for (int i = 0; i < treeSample.size(); i++){
                val += sigma[sigmaCounter] * covPairsIndicators[a][i];
                sigmaCounter++;
            }
                
            best = std::max(best, val / B);
        }
        total += best;
    }
    result.RademacherComplex =  static_cast<float>(total / T);
    
    // Precompute le[i][j] = true means edge[i] <= edge[j]
std::vector<std::vector<bool>> le(k, std::vector<bool>(k, false));
for (int i = 0; i < k; ++i)
    for (int j = 0; j < k; ++j)
        if (i != j){
            if (result.edges[j].first > 0){
                le[i][j] = SP.Poset.at(result.edges[j].first -1 )
                             .Tree.over(SP.Poset.at(result.edges[i].second - 1).Tree);
            }
        }
            

// perTrialVal[t][a] = (1/B) sum_i sigma_i * covPairsIndicators[a][i]  for trial t
std::vector<std::vector<double>> perTrialVal(T, std::vector<double>(k, 0.0));

for (int t = 0; t < T; ++t) {
    for (int i = 0; i < B; ++i)
        sigma[i] = coin(rng) ? 1 : -1;
    for (int a = 0; a < k; ++a) {
        double val = 0.0;
        for (int i = 0; i < B; ++i)
            val += sigma[i] * covPairsIndicators[a][i];
        perTrialVal[t][a] = val / B;
    }
}

// For a given antichain (subset of edge indices), its score is:
// (1/T) * sum_t  max_{a in A} perTrialVal[t][a]
/*auto scoreAntichain = [&](const std::vector<int>& A) -> double {
    if (A.empty()) return 0.0;
    double total = 0.0;
    for (int t = 0; t < T; ++t) {
        double best = -std::numeric_limits<double>::infinity();
        for (int a : A)
            best = std::max(best, perTrialVal[t][a]);
        total += best;
    }
    return total / T;
};*/

// Greedy antichain construction seeded from every singleton.
// At each step, add the edge with highest marginal score gain,
// provided it is incomparable to all current members.
double bestScore = 0.0;  // empty antichain scores 0
std::vector<int> bestAntichain;

for (int seed = 0; seed < k; ++seed) {
    std::vector<int> antichain = {seed};
    std::vector<bool> compatible(k, true);
    compatible[seed] = false;
    for (int j = 0; j < k; ++j)
        if (j != seed && (le[seed][j] || le[j][seed]))
            compatible[j] = false;

    // Cache current per-trial maxima for the antichain
    std::vector<double> trialMax(T);
    for (int t = 0; t < T; ++t)
        trialMax[t] = perTrialVal[t][seed];

    double currentScore = 0.0;
    for (int t = 0; t < T; ++t) currentScore += trialMax[t];
    currentScore /= T;

    while (true) {
        int bestJ = -1;
        double bestGain = 0.0;  // only add if gain > 0

        for (int j = 0; j < k; ++j) {
            if (!compatible[j]) continue;
            // Marginal gain: replacing trialMax[t] with max(trialMax[t], perTrialVal[t][j])
            double gain = 0.0;
            for (int t = 0; t < T; ++t)
                gain += std::max(0.0, perTrialVal[t][j] - trialMax[t]);
            gain /= T;
            if (gain > bestGain) { bestGain = gain; bestJ = j; }
        }

        if (bestJ == -1) break;

        // Update trialMax and compatibility
        for (int t = 0; t < T; ++t)
            trialMax[t] = std::max(trialMax[t], perTrialVal[t][bestJ]);
        currentScore += bestGain;

        antichain.push_back(bestJ);
        compatible[bestJ] = false;
        for (int j = 0; j < k; ++j)
            if (compatible[j] && (le[bestJ][j] || le[j][bestJ]))
                compatible[j] = false;
    }

    if (currentScore > bestScore) {
        bestScore = currentScore;
        bestAntichain = antichain;
    }
}

result.RademacherComplex2 = static_cast<float>(bestScore);
    
    // Now, computing the probabilities estimates for each pair and keeping track of the minimum computed.
    
    int minTimesCorrectlyClassified = bigTreeSample.size(); // We will keep track of the number of times the pair was correctly classified to make sure we deal with the minimum probability found.
    int sampleSize = bigTreeSample.size();
    
    int minLower = -2;
    int minUpper = -2;
    
    for (int cp = 0 ; cp < result.edges.size() ; cp++){
        int lIndx = result.edges[cp].first - 1;
        int uIndx = result.edges[cp].second - 1;
        
        int curSum = 0;
        
        if (lIndx < 0){
            for (pTree T : bigTreeSample){
                if (rho(T, SP.Poset.at(uIndx).Tree) == 0){
                    curSum++;
                }
            }
            if (result.nullCovering[cp]){
                if (curSum < minTimesCorrectlyClassified){
                    minTimesCorrectlyClassified = curSum;
                    minUpper = uIndx;
                    minLower = lIndx;
                } 
            }
            
        } else {
             for (pTree T : bigTreeSample){
                 if (rho(T, SP.Poset.at(lIndx).Tree) == rho(T, SP.Poset.at(uIndx).Tree)){
                     curSum++;
                 }
             }
            if (result.nullCovering[cp]){
                if (curSum < minTimesCorrectlyClassified){
                    minTimesCorrectlyClassified = curSum;
                    minUpper = uIndx;
                    minLower = lIndx;
                } 
            }
        }
        
        result.etaValues.push_back(static_cast<float>(curSum)/sampleSize);
    }
    
    result.nullCoveringProb = static_cast<float>(minTimesCorrectlyClassified)/sampleSize;
    result.minLower = minLower+1;
    result.minUpper = minUpper+1;
    
    curIndx = SP.firstRank.at(1);
    int rmax = static_cast<int>(SP.firstRank.size()) - 1;
    
    float addedEta = std::min(0.95f, 1 - result.nullCoveringProb);
    
    while (curIndx > -1){
        float rad_t = rad_threshold(1, q, B, result.RademacherComplex2, delta);
        
        result.rad_Ts_05.push_back(0.5f + rad_t);
        result.rad_Ts_p.push_back(addedEta + rad_t);
        
         curIndx = SP.Poset.at(curIndx).next;
    }
    
    
    for (int i = 0; i < SP.Poset.size(); i++){
        for (int j : SP.Poset.at(i).over){
            float omegaTemp =  static_cast<float>(rmax - SP.Poset.at(j).Tree.rank + 1)/(static_cast<float>(rmax));
            
            float rad_t = rad_threshold(omegaTemp, q, B, result.RademacherComplex2, delta);
        
            result.rad_Ts_05.push_back(0.5f + rad_t);
            result.rad_Ts_p.push_back(addedEta + rad_t);
        }
    }
    
    return result;
}

subPosetOutput SPanalisys(pTree tStar, vector<pTree> treeSample, vector<int> nSample, vector<pTree> bigTreeSample, vector<int> nBSample, subPoset SP, float q, float delta){
    
    subPosetOutput result;
    
    //vector<array<int, 2>> NullPairs;
    
    //Finding the pairs in C_null with upper tree of rank 1 and forming the edges above.
    
    int curIndx = SP.firstRank.at(1);
    //int numberPairsFirstLevel = 0;
    //int numPairs = 0;
    
    vector<vector<float>> covPairsIndicators;
    
    while (curIndx > -1){
        result.edges.push_back({-1, (curIndx+1)});
        vector<float> IndividualIndicators;
        if (rho(SP.Poset.at(curIndx).Tree,tStar) == 0){
            //NullPairs.push_back({curIndx,-1});
            result.nullCovering.push_back(true);
            //numberPairsFirstLevel++;
            //numPairs++;
            //cout << "Tree pair " << numPairs << " that belong to null pairs with 0 : \n";
            //SP.Poset.at(curIndx).Tree.print();
            //cout << "\n";
        } else {
            result.nullCovering.push_back(false);
        }
        
        float curmean = (rho(SP.Poset.at(curIndx).Tree, treeSample.at(0)) > 0);
        float curvar = 0;
        int curn = nSample.at(0);
        
        IndividualIndicators.push_back(curmean);
        
        for (int k = 1; k < treeSample.size(); k++){
            float Xvalue = (rho(SP.Poset.at(curIndx).Tree, treeSample.at(k)) > 0);
            
            IndividualIndicators.push_back(Xvalue);
            
            curmean = (curn*curmean + nSample.at(k)*Xvalue)/(nSample.at(k) + curn);
            curvar = (curn*curvar)/(nSample.at(k) + curn) + (curmean-Xvalue)*(curmean-Xvalue)*nSample.at(k)/curn;
            curn += nSample.at(k);
        }
        
        covPairsIndicators.push_back(IndividualIndicators);
        result.coveringMean.push_back(curmean);
        
        result.coveringVariance.push_back(curvar);
        
        curIndx = SP.Poset.at(curIndx).next;
    }
    
    //Finding the pairs above
    
    for (int i = 0; i < SP.Poset.size(); i++){
        for (int j : SP.Poset.at(i).over){
            result.edges.push_back({(i+1),(j+1)});
            vector<float> IndividualIndicators;
            if (rho(SP.Poset.at(i).Tree, tStar) == rho(SP.Poset.at(j).Tree, tStar)){
                //NullPairs.push_back({j,i});
                result.nullCovering.push_back(true);
                //numPairs++;
                //cout << "Tree pair " << numPairs << " that belong to null pairs: \n";
                //SP.Poset.at(j).Tree.print();
                //SP.Poset.at(i).Tree.print();
                //cout << "\n";
            } else {
                result.nullCovering.push_back(false);
            }
            
            float curmean = (rho(SP.Poset.at(j).Tree, treeSample.at(0)) > rho(SP.Poset.at(i).Tree, treeSample.at(0)));
            float curvar = 0;
            int curn = nSample.at(0);
            
            IndividualIndicators.push_back(curmean);
        
            for (int k = 1; k < treeSample.size(); k++){
                float Xvalue = (rho(SP.Poset.at(j).Tree, treeSample.at(k)) > rho(SP.Poset.at(i).Tree, treeSample.at(k)));
                
                IndividualIndicators.push_back(Xvalue);
            
                curmean = (curn*curmean + nSample.at(k)*Xvalue)/(nSample.at(k) + curn);
                curvar = (curn*curvar)/(nSample.at(k) + curn) + (curmean-Xvalue)*(curmean-Xvalue)*nSample.at(k)/curn;
                curn += nSample.at(k);
            }
            
            covPairsIndicators.push_back(IndividualIndicators);
        
            result.coveringMean.push_back(curmean);
        
            result.coveringVariance.push_back(curvar);
        }
    }
    
    // Estimating Rademacher Complexity bound
    
    int B = std::accumulate(nSample.begin(), nSample.end(), 0);
    
    int T = 10000;
    int k = (int)result.edges.size();

    std::mt19937 rng(42);
    std::bernoulli_distribution coin(0.5);

    double total = 0.0;
    std::vector<int> sigma(B);

    for (int t = 0; t < T; ++t) {
        // Draw sigma_i in {-1, +1}
        for (int i = 0; i < B; ++i)
            sigma[i] = coin(rng) ? 1 : -1;

        // Compute (1/m) sum_i sigma_i f_i(a) for each a, take max
        double best = -std::numeric_limits<double>::infinity();
        
        for (int a = 0; a < k; ++a) {
            double val = 0.0;
            int sigmaCounter = 0;
            for (int i = 0; i < treeSample.size(); i++){
                for (int l = 0; l < nSample[i]; l++){
                    val += sigma[sigmaCounter] * covPairsIndicators[a][i];
                    sigmaCounter++;
                }
            }
                
            best = std::max(best, val / B);
        }
        total += best;
    }
    result.RademacherComplex =  static_cast<float>(total / T);
    
    // Precompute le[i][j] = true means edge[i] <= edge[j]
std::vector<std::vector<bool>> le(k, std::vector<bool>(k, false));
for (int i = 0; i < k; ++i)
    for (int j = 0; j < k; ++j)
        if (i != j){
            if (result.edges[j].first > 0){
                le[i][j] = SP.Poset.at(result.edges[j].first - 1)
                             .Tree.over(SP.Poset.at(result.edges[i].second - 1).Tree);
            }
        }
            

// perTrialVal[t][a] = (1/B) sum_i sigma_i * covPairsIndicators[a][i]  for trial t
    
std::vector<std::vector<double>> perTrialVal(T, std::vector<double>(k, 0.0));

for (int t = 0; t < T; ++t) {
    for (int i = 0; i < B; ++i)
        sigma[i] = coin(rng) ? 1 : -1;
    for (int a = 0; a < k; ++a) {
        double val = 0.0;
        int sigmaCounter = 0;
        for (int i = 0; i < treeSample.size(); i++){
            for (int l = 0; l < nSample[i]; l++){
                val += sigma[sigmaCounter] * covPairsIndicators[a][i];
                sigmaCounter++;
            }
        }
        perTrialVal[t][a] = val / B;
    }
}

// For a given antichain (subset of edge indices), its score is:
// (1/T) * sum_t  max_{a in A} perTrialVal[t][a]
 
/*auto scoreAntichain = [&](const std::vector<int>& A) -> double {
    if (A.empty()) return 0.0;
    double total = 0.0;
    for (int t = 0; t < T; ++t) {
        double best = -std::numeric_limits<double>::infinity();
        for (int a : A)
            best = std::max(best, perTrialVal[t][a]);
        total += best;
    }
    return total / T;
};*/

// Greedy antichain construction seeded from every singleton.
// At each step, add the edge with highest marginal score gain,
// provided it is incomparable to all current members.
double bestScore = 0.0;  // empty antichain scores 0
std::vector<int> bestAntichain;

for (int seed = 0; seed < k; ++seed) {
    std::vector<int> antichain = {seed};
    std::vector<bool> compatible(k, true);
    compatible[seed] = false;
    for (int j = 0; j < k; ++j)
        if (j != seed && (le[seed][j] || le[j][seed]))
            compatible[j] = false;

    // Cache current per-trial maxima for the antichain
    std::vector<double> trialMax(T);
    for (int t = 0; t < T; ++t)
        trialMax[t] = perTrialVal[t][seed];

    double currentScore = 0.0;
    for (int t = 0; t < T; ++t) currentScore += trialMax[t];
    currentScore /= T;

    while (true) {
        int bestJ = -1;
        double bestGain = 0.0;  // only add if gain > 0

        for (int j = 0; j < k; ++j) {
            if (!compatible[j]) continue;
            // Marginal gain: replacing trialMax[t] with max(trialMax[t], perTrialVal[t][j])
            double gain = 0.0;
            for (int t = 0; t < T; ++t)
                gain += std::max(0.0, perTrialVal[t][j] - trialMax[t]);
            gain /= T;
            if (gain > bestGain) { bestGain = gain; bestJ = j; }
        }

        if (bestJ == -1) break;

        // Update trialMax and compatibility
        for (int t = 0; t < T; ++t)
            trialMax[t] = std::max(trialMax[t], perTrialVal[t][bestJ]);
        currentScore += bestGain;

        antichain.push_back(bestJ);
        compatible[bestJ] = false;
        for (int j = 0; j < k; ++j)
            if (compatible[j] && (le[bestJ][j] || le[j][bestJ]))
                compatible[j] = false;
    }

    if (currentScore > bestScore) {
        bestScore = currentScore;
        bestAntichain = antichain;
    }
}
    

result.RademacherComplex2 = static_cast<float>(bestScore);
    
    // Now, computing the probabilities estimates for each pair and keeping track of the minimum computed.
    int sampleSize = std::accumulate(nBSample.begin(), nBSample.end(), 0);
    int minTimesCorrectlyClassified = sampleSize;
    
    int minLower = -2;
    int minUpper = -2;
    
    for (int cp = 0 ; cp < result.edges.size() ; cp++){
        int lIndx = result.edges[cp].first - 1;
        int uIndx = result.edges[cp].second - 1;
        
        int curSum = 0;
        
        if (lIndx < 0){
            for (int i = 0; i < treeSample.size(); i++){
                pTree T = treeSample.at(i);
                if (rho(T, SP.Poset.at(uIndx).Tree) == 0){
                    curSum += nBSample.at(i);
                }
            }
            if (result.nullCovering[cp]){
                if (curSum < minTimesCorrectlyClassified){
                    minTimesCorrectlyClassified = curSum;
                    minUpper = uIndx;
                    minLower = lIndx;
                } 
            }
            
        } else {
             for (int i = 0; i < treeSample.size(); i++){
                 pTree T = treeSample.at(i);
                 if (rho(T, SP.Poset.at(lIndx).Tree) == rho(T, SP.Poset.at(uIndx).Tree)){
                     curSum += nBSample.at(i);
                 }
             }
            if (result.nullCovering[cp]){
                if (curSum < minTimesCorrectlyClassified){
                    minTimesCorrectlyClassified = curSum;
                    minUpper = uIndx;
                    minLower = lIndx;
                } 
            }
        }
        
        result.etaValues.push_back(static_cast<float>(curSum)/sampleSize);
    }
    
    result.nullCoveringProb = static_cast<float>(minTimesCorrectlyClassified)/sampleSize;
    result.minLower = (minLower+1);
    result.minUpper = (minUpper+1);
    
    curIndx = SP.firstRank.at(1);
    int rmax = static_cast<int>(SP.firstRank.size()) - 1;
    
    float addedEta = std::min(0.95f, 1 - result.nullCoveringProb);
    
    while (curIndx > -1){
        float rad_t = rad_threshold(1, q, B, result.RademacherComplex2, delta);
        
        result.rad_Ts_05.push_back(0.5f + rad_t);
        result.rad_Ts_p.push_back(addedEta + rad_t);
        
         curIndx = SP.Poset.at(curIndx).next;
    }
    
    
    for (int i = 0; i < SP.Poset.size(); i++){
        for (int j : SP.Poset.at(i).over){
            float omegaTemp =  static_cast<float>(rmax - SP.Poset.at(j).Tree.rank + 1)/(static_cast<float>(rmax));
            
            float rad_t = rad_threshold(omegaTemp, q, B, result.RademacherComplex2, delta);
        
            result.rad_Ts_05.push_back(0.5f + rad_t);
            result.rad_Ts_p.push_back(addedEta + rad_t);
        }
    }
    
    
    return result;
}