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

float kappa_threshold(float omegaT, float q, int N, int kappa){
    //cout << "kappa" << kappa << "\n";
    //cout << "log kappa" << log(kappa) << "\n";
    //cout << "log omegaT" << log(omegaT) << "\n";
    //cout << "log q" << log(q) << "\n";
    float prelim = log(kappa) + log(omegaT) - log(q);
    //cout << "Prelim is " << prelim << "\n"; 
    float thrs1 = std::max(prelim/(2*N), 0.0f);
    
    //cout << "Thrs1 " << thrs1 << "\n"; 
    
    return(sqrt(thrs1));
    
}

float rad_threshold(float omegaT, float q, int N, float Rad, float delta){
    float Udelta = 2*Rad + 2*sqrt(-log(delta)/(2*N));
    float prelim = log(omegaT) - log((q - delta));
    float thrs1 = std::max(2*prelim/N, 0.0f);
    
    return(Udelta + sqrt(thrs1));
}

subPosetOutput SPanalisys(pTree tStar, vector<pTree> treeSample, vector<pTree> bigTreeSample, subPoset SP, float q, float delta){
    
    subPosetOutput result;
    
    vector<array<int, 2>> NullPairs;
    
    //Finding the pairs in C_null with upper tree of rank 1 and forming the edges above
    
    int curIndx = SP.firstRank.at(1);
    int numberPairsFirstLevel = 0;
    
    vector<vector<float>> covPairsIndicators;
    //int numPairs = 0;
    while (curIndx > -1){
        result.edges.push_back({-1, (curIndx+1)});
        vector<float> IndividualIndicators;
        if (rho(SP.Poset.at(curIndx).Tree,tStar) == 0){
            NullPairs.push_back({curIndx,-1});
            result.nullCovering.push_back(true);
            numberPairsFirstLevel++;
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
                NullPairs.push_back({j,i});
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
    
    int T = 1000;
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
    
    // Now, computing the probabilities estimates for each pair and keeping track of the minimum computed.
    
    int minTimesCorrectlyClassified = bigTreeSample.size(); // We will keep track of the number of times the pair was correctly classified to make sure we deal with the minimum probability found.
    int sampleSize = bigTreeSample.size();
    
    int minLower = -2;
    int minUpper = -2;
    
    vector<vector<float>> covPairsIndicatorsB;
    
    for (int k = 0; k < numberPairsFirstLevel; k++){
        int curSum = 0;
        bool itBroke = false;
        for (pTree T : bigTreeSample){
            if (rho(T, SP.Poset.at(NullPairs[k][0]).Tree) == 0){
                curSum++;
            }
            if (curSum >= minTimesCorrectlyClassified) {
                itBroke = true;
                break;
            }
        }
        if (!itBroke){
            minTimesCorrectlyClassified = curSum;
            minUpper = NullPairs[k][0];
            //cout << "Minimum changed with tree pair " << (k+1) << " with 0: \n";
            //SP.Poset.at(NullPairs[k][0]).Tree.print();
            //cout << "\n";
        }
        
    }
    
    for (int k = numberPairsFirstLevel; k < NullPairs.size(); k++){
        int curSum = 0;
        bool itBroke = false;
        for (pTree T : bigTreeSample){
            if (rho(T, SP.Poset.at(NullPairs[k][0]).Tree) == rho(T, SP.Poset.at(NullPairs[k][1]).Tree)){
                curSum++;
            }
            if (curSum >= minTimesCorrectlyClassified) {
                itBroke = true;
                break;
            }
        }
        if (!itBroke){
            minTimesCorrectlyClassified = curSum;
            minUpper = NullPairs[k][0];
            minLower = NullPairs[k][1];
            //cout << "Minimum changed with tree pair " << (k+1) << ": \n";
            //SP.Poset.at(NullPairs[k][0]).Tree.print();
            //SP.Poset.at(NullPairs[k][1]).Tree.print();
            //cout << "\n";
        }
    }
    
    result.nullCoveringProb = static_cast<float>(minTimesCorrectlyClassified)/sampleSize;
    result.minLower = minLower+1;
    result.minUpper = minUpper+1;
    
    curIndx = SP.firstRank.at(1);
    int rmax = static_cast<int>(SP.firstRank.size()) - 1;
    
    float addedEta = std::min(0.95f, 1 - result.nullCoveringProb);
    
    while (curIndx > -1){
        float kappa_t = kappa_threshold(1, q, B, SP.Poset.at(curIndx).kappa);
        float rad_t = rad_threshold(1, q, B, result.RademacherComplex, delta);
        
        result.kappa_Ts_05.push_back(0.5f + kappa_t);
        result.kappa_Ts_p.push_back(addedEta + kappa_t);
        result.rad_Ts_05.push_back(0.5f + rad_t);
        result.rad_Ts_p.push_back(addedEta + rad_t);
        
         curIndx = SP.Poset.at(curIndx).next;
    }
    
    
    for (int i = 0; i < SP.Poset.size(); i++){
        for (int j : SP.Poset.at(i).over){
            float omegaTemp =  static_cast<float>(rmax - SP.Poset.at(j).Tree.rank + 1)/(static_cast<float>(rmax));
            
            float kappa_t = kappa_threshold(omegaTemp, q, B, SP.Poset.at(j).kappa);
            float rad_t = rad_threshold(omegaTemp, q, B, result.RademacherComplex, delta);
        
            result.kappa_Ts_05.push_back(0.5f + kappa_t);
            result.kappa_Ts_p.push_back(addedEta + kappa_t);
            result.rad_Ts_05.push_back(0.5f + rad_t);
            result.rad_Ts_p.push_back(addedEta + rad_t);
        }
    }
    
    return result;
}

subPosetOutput SPanalisys(pTree tStar, vector<pTree> treeSample, vector<int> nSample, vector<pTree> bigTreeSample, vector<int> nBSample, subPoset SP, float q, float delta){
    
    subPosetOutput result;
    
    vector<array<int, 2>> NullPairs;
    
    //Finding the pairs in C_null with upper tree of rank 1 and forming the edges above.
    
    int curIndx = SP.firstRank.at(1);
    int numberPairsFirstLevel = 0;
    //int numPairs = 0;
    
    vector<vector<float>> covPairsIndicators;
    
    while (curIndx > -1){
        result.edges.push_back({-1, (curIndx+1)});
        vector<float> IndividualIndicators;
        if (rho(SP.Poset.at(curIndx).Tree,tStar) == 0){
            NullPairs.push_back({curIndx,-1});
            result.nullCovering.push_back(true);
            numberPairsFirstLevel++;
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
                NullPairs.push_back({j,i});
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
    
    
    
    // Now, computing the probabilities estimates for each pair and keeping track of the minimum computed.
    int sampleSize = std::accumulate(nBSample.begin(), nBSample.end(), 0);
    int minTimesCorrectlyClassified = sampleSize;
    
    int minLower = -2;
    int minUpper = -2;
    
    for (int k = 0; k < numberPairsFirstLevel; k++){
        int curSum = 0;
        bool itBroke = false;
        for (int i = 0; i < treeSample.size(); i++){
            pTree T = treeSample.at(i);
            if (rho(T, SP.Poset.at(NullPairs[k][0]).Tree) == 0){
                curSum += nBSample.at(i);
            }
            if (curSum >= minTimesCorrectlyClassified) {
                itBroke = true;
                break;
            }
        }
        if(!itBroke){
            minTimesCorrectlyClassified = curSum;
            minUpper = NullPairs[k][0];
        }
    }
    
    for (int k = numberPairsFirstLevel; k < NullPairs.size(); k++){
        int curSum = 0;
        bool itBroke = false;
        for (int i = 0; i < treeSample.size(); i++){
            pTree T = treeSample.at(i);
            if (rho(T, SP.Poset.at(NullPairs[k][0]).Tree) == rho(T, SP.Poset.at(NullPairs[k][1]).Tree)){
                curSum += nBSample.at(i);
            }
            if (curSum >= minTimesCorrectlyClassified) {
                itBroke = true;
                break;
            }
        }
        if (!itBroke){
            minTimesCorrectlyClassified = curSum;
            minUpper = NullPairs[k][0];
            minLower = NullPairs[k][1];
        }
    }
    
    result.nullCoveringProb = static_cast<float>(minTimesCorrectlyClassified)/sampleSize;
    result.minLower = (minLower+1);
    result.minUpper = (minUpper+1);
    
    curIndx = SP.firstRank.at(1);
    int rmax = static_cast<int>(SP.firstRank.size()) - 1;
    
    float addedEta = std::min(0.95f, 1 - result.nullCoveringProb);
    
    while (curIndx > -1){
        float kappa_t = kappa_threshold(1, q, B, SP.Poset.at(curIndx).kappa);
        float rad_t = rad_threshold(1, q, B, result.RademacherComplex, delta);
        
        result.kappa_Ts_05.push_back(0.5f + kappa_t);
        result.kappa_Ts_p.push_back(addedEta + kappa_t);
        result.rad_Ts_05.push_back(0.5f + rad_t);
        result.rad_Ts_p.push_back(addedEta + rad_t);
        
         curIndx = SP.Poset.at(curIndx).next;
    }
    
    
    for (int i = 0; i < SP.Poset.size(); i++){
        for (int j : SP.Poset.at(i).over){
            float omegaTemp =  static_cast<float>(rmax - SP.Poset.at(j).Tree.rank + 1)/(static_cast<float>(rmax));
            
            float kappa_t = kappa_threshold(omegaTemp, q, B, SP.Poset.at(j).kappa);
            float rad_t = rad_threshold(omegaTemp, q, B, result.RademacherComplex, delta);
        
            result.kappa_Ts_05.push_back(0.5f + kappa_t);
            result.kappa_Ts_p.push_back(addedEta + kappa_t);
            result.rad_Ts_05.push_back(0.5f + rad_t);
            result.rad_Ts_p.push_back(addedEta + rad_t);
        }
    }
    
    
    return result;
}