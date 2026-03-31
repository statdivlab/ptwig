#include "pTree.h"
#include "mPhylo.h"
#include "subPoset.h"
#include "rho.h"
#include "coverTrees.h"
#include <iostream>
#include <set>
#include <queue>
#include <vector>
#include <array>
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

float idNullCoveringPairsComputation(pTree tStar, vector<pTree> treeSample, subPoset SP){
    vector<array<int, 2>> NullPairs; //Vector containing C_null pairs in the subPoset for tStar
    
    //Finding the pairs in C_null with upper tree of rank 1
    
    int curIndx = SP.firstRank.at(1);
    int numberPairsFirstLevel = 0;
    int numPairs = 0;
    while (curIndx > -1){
        if (rho(SP.Poset.at(curIndx).Tree,tStar) == 0){
            NullPairs.push_back({curIndx,-1});
            numberPairsFirstLevel++;
            numPairs++;
            cout << "Tree pair " << numPairs << " that belong to null pairs with 0 : \n";
            SP.Poset.at(curIndx).Tree.print();
            cout << "\n";
        }
        curIndx = SP.Poset.at(curIndx).next;
    }
    
    //Finding the pairs above
    
    for (int i = 0; i < SP.Poset.size(); i++){
        for (int j : SP.Poset.at(i).over){
            if (rho(SP.Poset.at(i).Tree, tStar) == rho(SP.Poset.at(j).Tree, tStar)){
                NullPairs.push_back({j,i});
                numPairs++;
                cout << "Tree pair " << numPairs << " that belong to null pairs: \n";
                SP.Poset.at(j).Tree.print();
                SP.Poset.at(i).Tree.print();
                cout << "\n";
            }
        }
    }
    
    // Now, computing the probabilities estimates for each pair and keeping track of the minimum computed.
    
    int minTimesCorrectlyClassified = treeSample.size(); // We will keep track of the number of times the pair was correctly classified to make sure we deal with the minimum probability found.
    int sampleSize = treeSample.size();
    
    for (int k = 0; k < numberPairsFirstLevel; k++){
        int curSum = 0;
        bool itBroke = false;
        for (pTree T : treeSample){
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
            cout << "Minimum changed with tree pair " << (k+1) << " with 0: \n";
            SP.Poset.at(NullPairs[k][0]).Tree.print();
            cout << "\n";
        }
        
    }
    
    for (int k = numberPairsFirstLevel; k < NullPairs.size(); k++){
        int curSum = 0;
        bool itBroke = false;
        for (pTree T : treeSample){
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
            cout << "Minimum changed with tree pair " << (k+1) << ": \n";
            SP.Poset.at(NullPairs[k][0]).Tree.print();
            SP.Poset.at(NullPairs[k][1]).Tree.print();
            cout << "\n";
        }
    }
    
    return(static_cast<float>(minTimesCorrectlyClassified)/sampleSize);
    
}

float idNullCoveringPairsComputation(pTree tStar, vector<pTree> treeSample, vector<int> nSample, subPoset SP){
    vector<array<int, 2>> NullPairs; //Vector containing C_null pairs in the subPoset for tStar
    
    //Finding the pairs in C_null with upper tree of rank 1
    
    int curIndx = SP.firstRank.at(1);
    int numberPairsFirstLevel = 0;
    int numPairs = 0;
    while (curIndx > -1){
        if (rho(SP.Poset.at(curIndx).Tree,tStar) == 0){
            NullPairs.push_back({curIndx,-1});
            numberPairsFirstLevel++;
            numPairs++;
            cout << "Tree pair " << numPairs << " that belong to null pairs with 0 : \n";
            SP.Poset.at(curIndx).Tree.print();
            cout << "\n";
        }
        curIndx = SP.Poset.at(curIndx).next;
    }
    
    //Finding the pairs above
    
    for (int i = 0; i < SP.Poset.size(); i++){
        for (int j : SP.Poset.at(i).over){
            if (rho(SP.Poset.at(i).Tree, tStar) == rho(SP.Poset.at(j).Tree, tStar)){
                NullPairs.push_back({j,i});
                numPairs++;
                cout << "Tree pair " << numPairs << " that belong to null pairs: \n";
                SP.Poset.at(j).Tree.print();
                SP.Poset.at(i).Tree.print();
                cout << "\n";
            }
        }
    }
    
    // Now, computing the probabilities estimates for each pair and keeping track of the minimum computed.
    
    int sampleSize = std::accumulate(nSample.begin(), nSample.end(), 0);
    int minTimesCorrectlyClassified = sampleSize; // We will keep track of the number of times the pair was correctly classified to make sure we deal with the minimum probability found.
    
    for (int k = 0; k < numberPairsFirstLevel; k++){
        int curSum = 0;
        bool itBroke = false;
        for (int i = 0; i < treeSample.size(); i++){
            pTree T = treeSample.at(i);
            if (rho(T, SP.Poset.at(NullPairs[k][0]).Tree) == 0){
                curSum += nSample.at(i);
            }
            if (curSum >= minTimesCorrectlyClassified) {
                itBroke = true;
                break;
            }
        }
        if(!itBroke){
            minTimesCorrectlyClassified = curSum;
            cout << "Minimum changed with tree pair " << (k+1) << " with 0: \n";
            SP.Poset.at(NullPairs[k][0]).Tree.print();
            cout << "\n";
        }
    }
    
    for (int k = numberPairsFirstLevel; k < NullPairs.size(); k++){
        int curSum = 0;
        bool itBroke = false;
        for (int i = 0; i < treeSample.size(); i++){
            pTree T = treeSample.at(i);
            if (rho(T, SP.Poset.at(NullPairs[k][0]).Tree) == rho(T, SP.Poset.at(NullPairs[k][1]).Tree)){
                curSum += nSample.at(i);
            }
            if (curSum >= minTimesCorrectlyClassified) {
                itBroke = true;
                break;
            }
        }
        if (!itBroke){
            minTimesCorrectlyClassified = curSum;
            cout << "Minimum changed with tree pair " << (k+1) << ": \n";
            SP.Poset.at(NullPairs[k][0]).Tree.print();
            SP.Poset.at(NullPairs[k][1]).Tree.print();
            cout << "\n";
        }
    }
    
    return(static_cast<float>(minTimesCorrectlyClassified)/sampleSize);
    
}