#include "pTree.h"
#include "mPhylo.h"
#include "subPoset.h"
#include "rho.h"
#include "coverTrees.h"
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


float boundThreshold(float omegaT, float q, int N, int nBound){
    float prelim = log(nBound) + log(omegaT) - log(q);
    float thrs1 = std::max(prelim/(2*N), 0.0f);
    
    
    return(sqrt(thrs1));
    
}



pTree FDRSearchGreedy(vector<pTree> treeSample, vector<int> nSample, vector<vector<oRho>> storedORho, subPoset SP, float q){

    int B = std::accumulate(nSample.begin(), nSample.end(), 0);
    int rmax = static_cast<int>(SP.firstRank.size() - 1);
    int numTrees = static_cast<int>(treeSample.size());
    

    // Score for the bottom-level transition (rank-1 node standing alone):
    // fraction of trees for which rho(node, T) > 0, shifted by lbEta
    auto scoreBase = [&](int nodeIdx) -> float {
        float sum = 0;
        for (int t = 0; t < numTrees; ++t){
            if (storedORho.at(nodeIdx).at(t).rho > 0)
                sum += nSample.at(t);
        }
        return (sum / B);
    };

    // Score for the transition from parentIdx -> childIdx (moving up):
    // fraction of trees where rho increases, shifted by lbEta of the parent
    auto scoreTransition = [&](int TaIdx, int TbIdx) -> float {
        float sum = 0;
        for (int t = 0; t < numTrees; ++t){
            if ((storedORho.at(TbIdx).at(t).rho - storedORho.at(TaIdx).at(t).rho) > 0)
                sum += nSample.at(t);
        }
            
        return (sum / B);
    };

    // Threshold for a given node
    auto threshold = [&](int nodeIdx, int childIdx) -> float {
        float omega = static_cast<float>(rmax - SP.Poset.at(nodeIdx).Tree.rank + 1)
                    / static_cast<float>(rmax);
        float taZeta = 1.0f/(3.0f);
        if (SP.Poset.at(nodeIdx).Tree.rank > 1){
            taZeta = min(SP.Poset.at(SP.Poset.at(nodeIdx).under[childIdx]).zeta, 0.5f);
        } 
        
        return (boundThreshold(omega, q, B, SP.Poset.at(nodeIdx).boundAntichain[childIdx]) + taZeta);
    };

    // ----------------------------------------------------------------
    // Step 1: scan rank-1 nodes in random order, pick first that passes
    // ----------------------------------------------------------------
    int curIndx = SP.firstRank.at(1);
    

    // Collect all rank-1 node indices
    vector<int> rank1Nodes;
    while (curIndx > -1) {
        rank1Nodes.push_back(curIndx);
        curIndx = SP.Poset.at(curIndx).next;
    }
    

    // Shuffle for random order
    auto rd  = std::random_device{};
    auto rng = std::default_random_engine{ rd() };
    shuffle(rank1Nodes.begin(), rank1Nodes.end(), rng);

    int current = -1;
    for (int idx : rank1Nodes) {
        float s = scoreBase(idx);
        if (s >= threshold(idx, 0)) {
            current = idx;
            break;
        }
    }
    if (current == -1) {
        pTree emptyTree  = pTree("();");
        return {emptyTree};
    }

    // ----------------------------------------------------------------
    // Step 2: greedily climb upward
    // ----------------------------------------------------------------
    while (!SP.Poset.at(current).over.empty()) {
        const vector<int>& candidates = SP.Poset.at(current).over;
        // Shuffle candidates for random order
        vector<int> shuffled(candidates.begin(), candidates.end());
        auto rd2  = std::random_device{};
        auto rng2 = std::default_random_engine{ rd2() };
        shuffle(shuffled.begin(), shuffled.end(),rng2);

        int next = -1;
        for (int upIdx : shuffled) {
            float s = scoreTransition(current, upIdx);
            int chldIdx = static_cast<int>(find(SP.Poset.at(upIdx).under.begin(), SP.Poset.at(upIdx).under.end(), current) - SP.Poset.at(upIdx).under.begin());
            if (s >= threshold(upIdx, chldIdx)) {
                next = upIdx;
                break;
            }
        }

        if (next == -1) {
            break;  // stuck — report current as result
        }

        current = next;
    }

    return { SP.Poset.at(current).Tree };
}