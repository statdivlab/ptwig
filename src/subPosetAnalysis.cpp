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

subPosetOutput SPanalisys(pTree tStar, vector<pTree> treeSample, vector<pTree> bigTreeSample, subPoset SP){
    
    subPosetOutput result;
    
    vector<array<int, 2>> NullPairs;
    
    //Finding the pairs in C_null with upper tree of rank 1 and forming the edges above
    
    int curIndx = SP.firstRank.at(1);
    int numberPairsFirstLevel = 0;
    //int numPairs = 0;
    while (curIndx > -1){
        result.edges.push_back({-1, (curIndx+1)});
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
        
        for (int k = 1; k < treeSample.size(); k++){
            float Xvalue = (rho(SP.Poset.at(curIndx).Tree, treeSample.at(k)) > 0);
            
            curmean = (curn*curmean + Xvalue)/(1+curn);
            curvar = (curn*curvar)/(1+curn) + (curmean-Xvalue)*(curmean-Xvalue)/curn;
            curn++;
        }
        
        result.coveringMean.push_back(curmean);
        
        result.coveringVariance.push_back(curvar);
        
        curIndx = SP.Poset.at(curIndx).next;
    }
    
    //Finding the pairs above
    
    for (int i = 0; i < SP.Poset.size(); i++){
        for (int j : SP.Poset.at(i).over){
            result.edges.push_back({(i+1),(j+1)});
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
        
            for (int k = 1; k < treeSample.size(); k++){
                float Xvalue = (rho(SP.Poset.at(j).Tree, treeSample.at(k)) > rho(SP.Poset.at(i).Tree, treeSample.at(k)));
            
                curmean = (curn*curmean + Xvalue)/(1+curn);
                curvar = (curn*curvar)/(1+curn) + (curmean-Xvalue)*(curmean-Xvalue)/curn;
                curn++;
            }
        
            result.coveringMean.push_back(curmean);
        
            result.coveringVariance.push_back(curvar);
        }
    }
    
    
    // Now, computing the probabilities estimates for each pair and keeping track of the minimum computed.
    
    int minTimesCorrectlyClassified = bigTreeSample.size(); // We will keep track of the number of times the pair was correctly classified to make sure we deal with the minimum probability found.
    int sampleSize = bigTreeSample.size();
    
    int minLower = -2;
    int minUpper = -2;
    
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
    
    return result;
}

subPosetOutput SPanalisys(pTree tStar, vector<pTree> treeSample, vector<int> nSample, vector<pTree> bigTreeSample, vector<int> nBSample, subPoset SP){
    
    subPosetOutput result;
    
    vector<array<int, 2>> NullPairs;
    
    //Finding the pairs in C_null with upper tree of rank 1 and forming the edges above
    
    int curIndx = SP.firstRank.at(1);
    int numberPairsFirstLevel = 0;
    //int numPairs = 0;
    while (curIndx > -1){
        result.edges.push_back({-1, (curIndx+1)});
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
        
        for (int k = 1; k < treeSample.size(); k++){
            float Xvalue = (rho(SP.Poset.at(curIndx).Tree, treeSample.at(k)) > 0);
            
            curmean = (curn*curmean + nSample.at(k)*Xvalue)/(nSample.at(k) + curn);
            curvar = (curn*curvar)/(nSample.at(k) + curn) + (curmean-Xvalue)*(curmean-Xvalue)*nSample.at(k)/curn;
            curn += nSample.at(k);
        }
        
        result.coveringMean.push_back(curmean);
        
        result.coveringVariance.push_back(curvar);
        
        curIndx = SP.Poset.at(curIndx).next;
    }
    
    //Finding the pairs above
    
    for (int i = 0; i < SP.Poset.size(); i++){
        for (int j : SP.Poset.at(i).over){
            result.edges.push_back({(i+1),(j+1)});
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
        
            for (int k = 1; k < treeSample.size(); k++){
                float Xvalue = (rho(SP.Poset.at(j).Tree, treeSample.at(k)) > rho(SP.Poset.at(i).Tree, treeSample.at(k)));
            
                curmean = (curn*curmean + nSample.at(k)*Xvalue)/(nSample.at(k) + curn);
                curvar = (curn*curvar)/(nSample.at(k) + curn) + (curmean-Xvalue)*(curmean-Xvalue)*nSample.at(k)/curn;
                curn += nSample.at(k);
            }
        
            result.coveringMean.push_back(curmean);
        
            result.coveringVariance.push_back(curvar);
        }
    }
    
    
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
    
    return result;
}