#include "pTree.h"
#include "mPhylo.h"
#include "rho.h"
#include "coverTrees.h"
#include <Rcpp.h>
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


using namespace Rcpp;
using namespace std;

float Stability(vector<pTree> treeSample, float w, pTree V, Split s){
    
    int K = static_cast<int>(treeSample.size());
    pTree U = V.Remove(s);
    
    float Sum = 0;
    
    for (pTree T : treeSample){
        if ((rho(V, T) - rho(U, T)) > 0){
            Sum++;
        }
    }
    
    return Sum/K;
}

float Stability(vector<pTree> treeSample, pTree V, Split s, const vector<float>& rhoV, int K){
    
    pTree U = V.Remove(s);
    
    float Sum = 0;
    
    for (int i = 0; i < K; i++){
        if ((rhoV[i] - rho(U, treeSample[i])) > 0){
            Sum++;
        }
    }
    
    return Sum/K;
}

float Stability(vector<pTree> treeSample, vector<int> nSample, float w, pTree V, Split s){
    
    int K = std::accumulate(nSample.begin(), nSample.end(), 0);
    pTree U = V.Remove(s);
    
    float Sum = 0;
    
    int curIndx = 0;
    for (pTree T : treeSample){
        if ((rho(V, T) - rho(U, T)) > 0){
            Sum += nSample.at(curIndx);
        }
        curIndx++;
    }
    
    return Sum/K;
}

float Stability(vector<pTree> treeSample, vector<int> nSample, pTree V, Split s, const vector<float>& rhoV, int K){
    
    pTree U = V.Remove(s);
    
    float Sum = 0;
    
    for (int i = 0; i < treeSample.size(); i++){
        if ((rhoV[i] - rho(U, treeSample[i])) > 0){
            Sum += nSample.at(i);
        }
    }
    
    return Sum/K;
}

float Stability(vector<pTree> treeSample, float w, pTree V, string a){
    int K = static_cast<int>(treeSample.size());
    pTree U = V.Remove(a);
    
    float Sum = 0;
    
    for (pTree T : treeSample){
        if ((rho(V, T) - rho(U, T)) > 0){
            Sum++;
        }
    }
    
    return Sum/K;
}

float Stability(vector<pTree> treeSample, pTree V, string a, const vector<float>& rhoV, int K){
    pTree U = V.Remove(a);
    
    float Sum = 0;
    
    for (int i = 0; i < K; i++){
        if ((rhoV[i] - rho(U, treeSample[i])) > 0){
            Sum++;
        }
    }
    
    return Sum/K;
}

float Stability(vector<pTree> treeSample, vector<int> nSample, float w, pTree V, string a){
    int K = std::accumulate(nSample.begin(), nSample.end(), 0);
    pTree U = V.Remove(a);
    
    float Sum = 0;
    
    int curIndx = 0;
    for (pTree T : treeSample){
        if ((rho(V, T) - rho(U, T)) > 0){
            Sum += nSample.at(curIndx);
        }
        curIndx++;
    }
    
    return Sum/K;
}

float Stability(vector<pTree> treeSample, vector<int> nSample, pTree V, string a, const vector<float>& rhoV, int K){
    
    pTree U = V.Remove(a);
    
    float Sum = 0;
    
    for (int i = 0; i < treeSample.size(); i++){
        if ((rhoV[i] - rho(U, treeSample[i])) > 0){
            Sum += nSample.at(i);
        }
    }
    
    return Sum/K;
}

float MinimumStability(vector<pTree> treeSample, float w, pTree V){
    float Result = 1;
    
    int K = static_cast<int>(treeSample.size());
    
    for (string a : V.leafSet){
        pTree U = V.Remove(a);
        
        float Sum = 0;
        
        for (pTree T : treeSample){
            if ((rho(V, T) > rho(U, T)) > 0){
                Sum++;
            }
        }
        
        if (Result > Sum/K){
            Result = Sum/K;
        }
    }
    
    for (Split s : V.intSplits){
        pTree U = V.Remove(s);
        
        float Sum = 0;
        
        for (pTree T : treeSample){
            if ((rho(V, T) > rho(U, T)) > 0){
                Sum++;
            }
        }
        
        if (Result > Sum/K){
            Result = Sum/K;
        }
    }
    
    return Result;
}

float MinimumStability(vector<pTree> treeSample,
                       const vector<float>& rhoV, // precomputed rho(V, Z)
                       pTree V, int K){
    float Result = 1.0f;
    
    auto checkFeature = [&](pTree U_minus) {
        float Sum = 0;
        for (int i = 0; i < K; i++) {
            if (rhoV[i] > rho(U_minus, treeSample[i]))
                Sum++;
        }
        Result = min(Result, Sum / K);
    };
    
    for (const string& a : V.leafSet)    checkFeature(V.Remove(a));
    for (const Split&  s : V.intSplits)  checkFeature(V.Remove(s));
    
    return Result;
}

float MinimumStability(vector<pTree> treeSample, vector<int> nSample, float w, pTree V){
    float Result = 1;
    
    int K = std::accumulate(nSample.begin(), nSample.end(), 0);
    
    for (string a : V.leafSet){
        pTree U = V.Remove(a);
        
        float Sum = 0;
        
        int curIndx = 0;
        for (pTree T : treeSample){
            if ((rho(V, T) > rho(U, T)) > 0){
                Sum += nSample.at(curIndx);
            }
            curIndx++;
        }
        
        if (Result > Sum/K){
            Result = Sum/K;
        }
    }
    
    for (Split s : V.intSplits){
        pTree U = V.Remove(s);
        
        float Sum = 0;
        
        int curIndx = 0;
        for (pTree T : treeSample){
            if ((rho(V, T) > rho(U, T)) > 0){
                Sum += nSample.at(curIndx);
            }
            curIndx++;
        }
        
        if (Result > Sum/K){
            Result = Sum/K;
        }
    }
    
    return Result;
}

float MinimumStability(vector<pTree>& treeSample, vector<int>& nSample,
                       const vector<float>& rhoV,   // precomputed rho(V, Z)
                       pTree V, int K) {
    float Result = 1.0f;

    auto checkFeature = [&](pTree U_minus) {
        float Sum = 0;
        for (int i = 0; i < treeSample.size(); i++) {
            if (rhoV[i] > rho(U_minus, treeSample[i]))
                Sum += nSample[i];
        }
        Result = min(Result, Sum / K);
    };

    for (const string& a : V.leafSet)    checkFeature(V.Remove(a));
    for (const Split&  s : V.intSplits)  checkFeature(V.Remove(s));

    return Result;
}

vector<pTree> stableSearch(vector<pTree> treeSample, set<string> compLeafSet, float alpha){
    
    vector<pTree> CollectionTrees;
    pTree U = pTree("();");
    int B = static_cast<int>(treeSample.size());
    
    bool RecentlyRemoved = false;
    int Count = 0;
    int highRank = 0;
    
    // --- Fix 3: visited set to prevent cycles ---
    set<string> visited;

    // --- Fix 1: cache rho(U, Z) — recomputed only when U changes ---
    vector<float> rhoU(B,0.0f);
    
    while (true) {
        cout << "Entering cycle " << Count << "\n"<< std::flush;
        vector<pTree> AllV = coverTrees(U, compLeafSet);
        
        auto rd  = std::random_device{};
        auto rng = std::default_random_engine{ rd() };
        shuffle(begin(AllV), end(AllV), rng);

        // --- Fix 2: find best candidate, saving its rhoV as we go ---
        int   IndexMax  = -1;
        float MaxValue  = 0;
        vector<float> bestRhoV(B);

        // Separate trackers for large/small leaf set tiebreak (preserved from original)
        int   IndexMaxS = -1, IndexMaxL = -1;
        float MaxValueS = 0,  MaxValueL = 0;
        vector<float> bestRhoVS(B), bestRhoVL(B);
        
        for (int indx = 0; indx < (int)AllV.size(); indx++){
            const pTree& V = AllV[indx];

            // Compute rho(V, Z) once, reuse for score and later for MinStab
            vector<float> rhoV(B);
            double tempSum = 0;
            for (int i = 0; i < B; i++) {
                rhoV[i] = rho(V, treeSample[i]);
                if (rhoV[i] - rhoU[i] > 0)
                    tempSum++;
            }

            float score = (float)(tempSum / B);

            if ((V.leafSet.size() > U.leafSet.size()) && (V.leafSet.size() > 4)) {
                if (score > MaxValueL) {
                    MaxValueL = score;
                    IndexMaxL = indx;
                    bestRhoVL = rhoV;
                }
            } else {
                if (score > MaxValueS) {
                    MaxValueS = score;
                    IndexMaxS = indx;
                    bestRhoVS = rhoV;
                }
            }
        }
        if (MaxValueL > MaxValueS) {
            MaxValue  = MaxValueL;
            IndexMax  = IndexMaxL;
            bestRhoV  = bestRhoVL;
        } else {
            MaxValue  = MaxValueS;
            IndexMax  = IndexMaxS;
            bestRhoV  = bestRhoVS;
        }

        cout << "Best candidate found \n"<< std::flush;
        pTree V = AllV.at(IndexMax);

        if (RecentlyRemoved && MaxValue < alpha)
            break;

        RecentlyRemoved = false;
        
        // --- Sequential removal block ---
        // Each time V is modified, bestRhoV is recomputed for the new V
        if (MaxValue < alpha) {

            // Helper lambda: recompute rhoV for current V
            auto recomputeRhoV = [&]() {
                for (int i = 0; i < B; i++)
                    bestRhoV[i] = rho(V, treeSample[i]);
            };

            auto tryRemoveLeaves = [&](const set<string>& leaves) {
                for (const string& a : leaves) {
                    // Only attempt if leaf still present after prior removals
                    if (V.leafSet.count(a)) {
                        if (Stability(treeSample, V, a, bestRhoV, B) < alpha) {
                            V = V.Remove(a);
                            recomputeRhoV();   // V changed — cache is stale
                            RecentlyRemoved = true;
                        }
                    }
                }
            };

            auto tryRemoveSplits = [&](const set<Split>& splits) {
                for (Split s : splits) {
                    Split s2 = s.TDR(V.leafSet);
                    if (s2.isInternal()) {
                        if (Stability(treeSample, V, s2, bestRhoV, B) < alpha) {
                            V = V.Remove(s2);
                            recomputeRhoV();   // V changed — cache is stale
                            RecentlyRemoved = true;
                        }
                    }
                }
            };

            if ((V.leafSet.size() > U.leafSet.size()) && (V.returnRank() > 1)) {
                tryRemoveLeaves(U.leafSet);
                tryRemoveSplits(V.intSplits);   // use V.intSplits as in original
            } else {
                tryRemoveLeaves(U.leafSet);
                tryRemoveSplits(U.intSplits);   // use U.intSplits as in original
            }
        }

        // MinimumStability now reuses bestRhoV — no extra rho(V,Z) calls
        if (MinimumStability(treeSample, bestRhoV, V, B) < alpha)
            break;
        
        cout << "The candidate was finally selected \n"<< std::flush;
        U = V;
        
        // --- Fix 3: cycle detection ---
        mPhylo mpU = mPhylo(U);
        string newickU = mpU.toNewick();
        if (visited.count(newickU)) {
            cout << "WARNING: Cycle detected, stopping.\n"<< std::flush;
            break;
        }
        visited.insert(newickU);
        
        // Update rhoU since U just changed
        for (int i = 0; i < B; i++)
            rhoU[i] = rho(U, treeSample[i]);

        if (U.returnRank() > highRank) {
            CollectionTrees.clear();
            CollectionTrees.push_back(U);
            highRank = U.returnRank();
        } else if (U.returnRank() == highRank) {
            CollectionTrees.push_back(U);
        }

        Count++;
        
        if (U.returnRank() == 2*compLeafSet.size()-7){
            break;
        }
        
        if (Count > 4*compLeafSet.size()-7){
            cout<< "WARNING!: Forced to stop. \n"<< std::flush;
            break;
        }
        
    }
    
    return CollectionTrees;
    
}

vector<pTree> stableSearch(vector<pTree> treeSample, vector<int> nSample, set<string> compLeafSet, float alpha){
    
    vector<pTree> CollectionTrees;
    pTree U = pTree("();");
    int B = std::accumulate(nSample.begin(), nSample.end(), 0);
    int N = static_cast<int>(treeSample.size());
    
    bool RecentlyRemoved = false;
    int Count = 0;
    int highRank = 0;
    
    // --- Fix 3: visited set to prevent cycles ---
    set<string> visited;

    // --- Fix 1: cache rho(U, Z) — recomputed only when U changes ---
    vector<float> rhoU(N,0.0f);
    
    while (true) {
        cout << "Entering cycle " << Count << "\n"<< std::flush;
        vector<pTree> AllV = coverTrees(U, compLeafSet);
        
        auto rd  = std::random_device{};
        auto rng = std::default_random_engine{ rd() };
        shuffle(begin(AllV), end(AllV), rng);

        // --- Fix 2: find best candidate, saving its rhoV as we go ---
        int   IndexMax  = -1;
        float MaxValue  = 0;
        vector<float> bestRhoV(N);

        // Separate trackers for large/small leaf set tiebreak (preserved from original)
        int   IndexMaxS = -1, IndexMaxL = -1;
        float MaxValueS = 0,  MaxValueL = 0;
        vector<float> bestRhoVS(N), bestRhoVL(N);
        
        for (int indx = 0; indx < (int)AllV.size(); indx++){
            const pTree& V = AllV[indx];

            // Compute rho(V, Z) once, reuse for score and later for MinStab
            vector<float> rhoV(N);
            double tempSum = 0;
            for (int i = 0; i < N; i++) {
                rhoV[i] = rho(V, treeSample[i]);
                if (rhoV[i] - rhoU[i] > 0)
                    tempSum += nSample[i];
            }

            float score = (float)(tempSum / B);

            if ((V.leafSet.size() > U.leafSet.size()) && (V.leafSet.size() > 4)) {
                if (score > MaxValueL) {
                    MaxValueL = score;
                    IndexMaxL = indx;
                    bestRhoVL = rhoV;
                }
            } else {
                if (score > MaxValueS) {
                    MaxValueS = score;
                    IndexMaxS = indx;
                    bestRhoVS = rhoV;
                }
            }
        }
        if (MaxValueL > MaxValueS) {
            MaxValue  = MaxValueL;
            IndexMax  = IndexMaxL;
            bestRhoV  = bestRhoVL;
        } else {
            MaxValue  = MaxValueS;
            IndexMax  = IndexMaxS;
            bestRhoV  = bestRhoVS;
        }

        cout << "Best candidate was found \n"<< std::flush;
        pTree V = AllV.at(IndexMax);

        if (RecentlyRemoved && MaxValue < alpha)
            break;

        RecentlyRemoved = false;
        
        // --- Sequential removal block ---
        // Each time V is modified, bestRhoV is recomputed for the new V
        if (MaxValue < alpha) {

            // Helper lambda: recompute rhoV for current V
            auto recomputeRhoV = [&]() {
                for (int i = 0; i < N; i++)
                    bestRhoV[i] = rho(V, treeSample[i]);
            };

            auto tryRemoveLeaves = [&](const set<string>& leaves) {
                for (const string& a : leaves) {
                    // Only attempt if leaf still present after prior removals
                    if (V.leafSet.count(a)) {
                        if (Stability(treeSample, nSample, V, a, bestRhoV, B) < alpha) {
                            V = V.Remove(a);
                            recomputeRhoV();   // V changed — cache is stale
                            RecentlyRemoved = true;
                        }
                    }
                }
            };

            auto tryRemoveSplits = [&](const set<Split>& splits) {
                for (Split s : splits) {
                    Split s2 = s.TDR(V.leafSet);
                    if (s2.isInternal()) {
                        if (Stability(treeSample, nSample, V, s2, bestRhoV, B) < alpha) {
                            V = V.Remove(s2);
                            recomputeRhoV();   // V changed — cache is stale
                            RecentlyRemoved = true;
                        }
                    }
                }
            };

            if ((V.leafSet.size() > U.leafSet.size()) && (V.returnRank() > 1)) {
                tryRemoveLeaves(U.leafSet);
                tryRemoveSplits(V.intSplits);   // use V.intSplits as in original
            } else {
                tryRemoveLeaves(U.leafSet);
                tryRemoveSplits(U.intSplits);   // use U.intSplits as in original
            }
        }

        // MinimumStability now reuses bestRhoV — no extra rho(V,Z) calls
        if (MinimumStability(treeSample, nSample, bestRhoV, V, B) < alpha)
            break;
        
        cout << "The candidate was finally selected \n"<< std::flush;
        U = V;
        
        // --- Fix 3: cycle detection ---
        mPhylo mpU = mPhylo(U);
        string newickU = mpU.toNewick();
        if (visited.count(newickU)) {
            cout << "WARNING: Cycle detected, stopping.\n"<< std::flush;
            break;
        }
        visited.insert(newickU);
        
        // Update rhoU since U just changed
        for (int i = 0; i < N; i++)
            rhoU[i] = rho(U, treeSample[i]);

        if (U.returnRank() > highRank) {
            CollectionTrees.clear();
            CollectionTrees.push_back(U);
            highRank = U.returnRank();
        } else if (U.returnRank() == highRank) {
            CollectionTrees.push_back(U);
        }

        Count++;
        
        if (U.returnRank() == 2*compLeafSet.size()-7){
            break;
        }
        
        if (Count > 4*compLeafSet.size()-7){
            cout<< "WARNING!: Forced to stop. \n"<< std::flush;
            break;
        }
        
    }
    
    return CollectionTrees;
    
}

// [[Rcpp::export]]
CharacterVector stableSearchRcpp(CharacterVector treeSampleR,
                                 CharacterVector compLeafSetR,
                                 double alphaR) {
    //
    // 1. Convert treeSampleR → vector<pTree>
    //
    std::vector<pTree> treeSample;
    treeSample.reserve(treeSampleR.size());

    for (int i = 0; i < treeSampleR.size(); i++) {
        if (treeSampleR[i] == NA_STRING)
            stop("treeSample cannot contain NA.");
        treeSample.emplace_back(
            pTree(as<std::string>(treeSampleR[i]))
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
    // 3. Call C++ function
    //
    float alpha = static_cast<float>(alphaR);
    std::vector<pTree> result = stableSearch(treeSample, compLeafSet, alpha);

    //
    // 4. Convert vector<pTree> → CharacterVector
    //
    CharacterVector out(result.size());
    for (size_t i = 0; i < result.size(); i++) {
        mPhylo rP = mPhylo(result[i]);
        out[i] = rP.toNewick();
    }

    return out;
}

// Build the vector<pTree> sample from a CharacterVector of Newick strings.
static std::vector<pTree> buildTreeSample(CharacterVector treeSampleR) {
    std::vector<pTree> treeSample;
    treeSample.reserve(treeSampleR.size());
    for (int i = 0; i < treeSampleR.size(); i++) {
        if (treeSampleR[i] == NA_STRING)
            stop("treeSample cannot contain NA.");
        treeSample.emplace_back(pTree(as<std::string>(treeSampleR[i])));
    }
    return treeSample;
}

// [[Rcpp::export]]
DataFrame computeStabilityRcpp(CharacterVector treeR,
                               CharacterVector treeSampleR) {
    //
    // 1. Build the tree V
    //
    if (treeR.size() != 1 || treeR[0] == NA_STRING)
        stop("treeR must be a single non-NA Newick string.");
    pTree V = pTree(as<std::string>(treeR[0]));

    //
    // 2. Build the tree sample
    //
    std::vector<pTree> treeSample = buildTreeSample(treeSampleR);
    int K = static_cast<int>(treeSample.size());

    //
    // 3. Precompute rhoV = rho(V, Z) for every tree Z in the sample
    //
    std::vector<float> rhoV(K);
    for (int i = 0; i < K; i++)
        rhoV[i] = static_cast<float>(rho(V, treeSample[i]));

    //
    // 4. Stability of every leaf and every internal split of V
    //
    std::vector<std::string> features;
    std::vector<std::string> types;
    std::vector<double>      stabilities;

    for (const std::string& a : V.leafSet) {
        features.push_back(a);
        types.push_back("leaf");
        stabilities.push_back(Stability(treeSample, V, a, rhoV, K));
    }
    for (const Split& s : V.intSplits) {
        features.push_back(s.printSt());
        types.push_back("split");
        stabilities.push_back(Stability(treeSample, V, s, rhoV, K));
    }

    return DataFrame::create(
        _["feature"]          = features,
        _["type"]             = types,
        _["stability"]        = stabilities,
        _["stringsAsFactors"] = false
    );
}

// [[Rcpp::export]]
DataFrame computeStabilityRcppS(CharacterVector treeR,
                                CharacterVector treeSampleR,
                                IntegerVector nSampleR) {
    //
    // 1. Build the tree V
    //
    if (treeR.size() != 1 || treeR[0] == NA_STRING)
        stop("treeR must be a single non-NA Newick string.");
    pTree V = pTree(as<std::string>(treeR[0]));

    //
    // 2. Build the (summarized) tree sample and its multiplicities
    //
    std::vector<pTree> treeSample = buildTreeSample(treeSampleR);

    if (nSampleR.size() != (int)treeSample.size())
        stop("treeSample and nSample must have the same length.");

    std::vector<int> nSample;
    nSample.reserve(nSampleR.size());
    for (int i = 0; i < nSampleR.size(); i++)
        nSample.push_back(static_cast<int>(nSampleR[i]));

    int N = static_cast<int>(treeSample.size());          // # of unique trees
    int K = std::accumulate(nSample.begin(), nSample.end(), 0); // total count

    //
    // 3. Precompute rhoV = rho(V, Z) for every unique tree Z in the sample
    //
    std::vector<float> rhoV(N);
    for (int i = 0; i < N; i++)
        rhoV[i] = static_cast<float>(rho(V, treeSample[i]));

    //
    // 4. Stability of every leaf and every internal split of V
    //
    std::vector<std::string> features;
    std::vector<std::string> types;
    std::vector<double>      stabilities;

    for (const std::string& a : V.leafSet) {
        features.push_back(a);
        types.push_back("leaf");
        stabilities.push_back(Stability(treeSample, nSample, V, a, rhoV, K));
    }
    for (const Split& s : V.intSplits) {
        features.push_back(s.printSt());
        types.push_back("split");
        stabilities.push_back(Stability(treeSample, nSample, V, s, rhoV, K));
    }

    return DataFrame::create(
        _["feature"]          = features,
        _["type"]             = types,
        _["stability"]        = stabilities,
        _["stringsAsFactors"] = false
    );
}

// [[Rcpp::export]]
CharacterVector stableSearchRcppS(CharacterVector treeSampleR,
                                 IntegerVector nSampleR,
                                 CharacterVector compLeafSetR,
                                 double alphaR) {
    //
    // 1. Convert treeSampleR → vector<pTree>
    //
    std::vector<pTree> treeSample;
    treeSample.reserve(treeSampleR.size());

    for (int i = 0; i < treeSampleR.size(); i++) {
        if (treeSampleR[i] == NA_STRING)
            stop("treeSample cannot contain NA.");
        treeSample.emplace_back(
            pTree(as<std::string>(treeSampleR[i]))
        );
    }
    
    std::vector<int> nSample;
    
    for (int i = 0; i < nSampleR.size(); i++){
        nSample.push_back(static_cast<int>(nSampleR[i]));
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
    // 3. Call C++ function
    //
    float alpha = static_cast<float>(alphaR);
    std::vector<pTree> result = stableSearch(treeSample, nSample, compLeafSet, alpha);

    //
    // 4. Convert vector<pTree> → CharacterVector
    //
    CharacterVector out(result.size());
    for (size_t i = 0; i < result.size(); i++) {
        mPhylo rP = mPhylo(result[i]);
        out[i] = rP.toNewick();
    }

    return out;
}