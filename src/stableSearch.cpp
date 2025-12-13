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
        if ((rho(V, T) > rho(U, T)) > 0){
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
        if ((rho(V, T) > rho(U, T)) > 0){
            Sum += nSample.at(curIndx);
        }
        curIndx++;
    }
    
    return Sum/K;
}

float Stability(vector<pTree> treeSample, float w, pTree V, string a){
    int K = static_cast<int>(treeSample.size());
    pTree U = V.Remove(a);
    
    float Sum = 0;
    
    for (pTree T : treeSample){
        if ((rho(V, T) > rho(U, T)) > 0){
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
        if ((rho(V, T) > rho(U, T)) > 0){
            Sum += nSample.at(curIndx);
        }
        curIndx++;
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

vector<pTree> stableSearch(vector<pTree> treeSample, set<string> compLeafSet, float alpha){
    
    vector<pTree> CollectionTrees;
    pTree U = pTree("();");
    int B = static_cast<int>(treeSample.size());
    
    bool RecentlyRemoved = false;
    int Count = 0;
    int highRank = 0;
    
    while (true) {
        
        vector<pTree> AllV = coverTrees(U, compLeafSet);
        
        auto rd = std::random_device {};
        auto rng = std::default_random_engine { rd() };
        shuffle(begin(AllV), std::end(AllV), rng);
        
        //Version where maximal value is found.
        
        int IndexMaxS = -1;
        int IndexMaxL = -1;
        int IndexMax  = -1;
        float MaxValueS = 0;
        float MaxValueL = 0;
        float MaxValue  = 0;
        int indx = 0;
        for (pTree V : AllV){
            double tempSum = 0;
            for (pTree Z : treeSample){
                if ((rho(V, Z) - rho(U, Z)) > 0){
                    tempSum++;
                }
            }
            if ((V.leafSet.size() > U.leafSet.size()) && (V.leafSet.size()>4)){
                if ((tempSum/B) > MaxValueL){
                    MaxValueL = tempSum/B;
                    IndexMaxL = indx;
                }
            } else {
                if ((tempSum/B) > MaxValueS){
                    MaxValueS = tempSum/B;
                    IndexMaxS = indx;
                }
            }
            indx++;
        }
        if (MaxValueL > MaxValueS){
            MaxValue = MaxValueL;
            IndexMax = IndexMaxL;
        } else {
            MaxValue = MaxValueS;
            IndexMax = IndexMaxS;
        }
        
        pTree V = AllV.at(IndexMax);
        
        
        if ((RecentlyRemoved) && (MaxValue < alpha)){
            break;
        }
        
        RecentlyRemoved = false;
        
        if (MaxValue < alpha){
            if((V.leafSet.size() > U.leafSet.size()) && (V.returnRank()>1)){
                for (string a : U.leafSet){
                    if (Stability(treeSample, 1.0, V, a) < alpha){
                        V = V.Remove(a);
                        RecentlyRemoved = true;
                    }
                }
                set<Split> FirstSplits = V.intSplits;
                for (Split s : FirstSplits){
                    Split s2 = s.TDR(V.leafSet);
                    if (s2.isInternal()){
                        if (Stability(treeSample, 1.0, V, s2) < alpha){
                            V = V.Remove(s2);
                            RecentlyRemoved = true;
                        }
                    }
                }
            } else {
                for (string a : U.leafSet){
                    if (Stability(treeSample, 1.0, V, a) < alpha){
                        V = V.Remove(a);
                        RecentlyRemoved = true;
                    }
                }
                set<Split> FirstSplits = U.intSplits;
                for (Split s : FirstSplits){
                    Split s2 = s.TDR(V.leafSet);
                    if (s2.isInternal()){
                        if (Stability(treeSample, 1.0, V, s2) < alpha){
                            V = V.Remove(s2);
                            RecentlyRemoved = true;
                        }
                    }
                }
            }
        }
        
        if (MinimumStability(treeSample,1.0,V) < alpha){
            break;
        }
        
        U = V;
        
        if (U.returnRank() > highRank){
            CollectionTrees.clear();
            CollectionTrees.push_back(U);
            highRank = U.returnRank();
        } else if (U.returnRank() == highRank){
            CollectionTrees.push_back(U);
        }
        
        Count++;
        
        
        if (U.returnRank() == 2*compLeafSet.size()-7){
            break;
        }
        
        if (Count > 4*compLeafSet.size()-7){
            cout<< "WARNING!: Forced to stop. \n";
            break;
        }
        
    }
    
    return CollectionTrees;
    
}

vector<pTree> stableSearch(vector<pTree> treeSample, vector<int> nSample, set<string> compLeafSet, float alpha){
    
    vector<pTree> CollectionTrees;
    pTree U = pTree("();");
    int B = std::accumulate(nSample.begin(), nSample.end(), 0);
    
    bool RecentlyRemoved = false;
    int Count = 0;
    int highRank = 0;
    
    while (true) {
        
        vector<pTree> AllV = coverTrees(U, compLeafSet);
        
        auto rd = std::random_device {};
        auto rng = std::default_random_engine { rd() };
        shuffle(begin(AllV), std::end(AllV), rng);
        
        //Version where maximal value is found.
        
        int IndexMaxS = -1;
        int IndexMaxL = -1;
        int IndexMax  = -1;
        float MaxValueS = 0;
        float MaxValueL = 0;
        float MaxValue  = 0;
        int indx = 0;
        for (pTree V : AllV){
            double tempSum = 0;
            int curIndx = 0;
            for (pTree Z : treeSample){
                if ((rho(V, Z) - rho(U, Z)) > 0){
                    tempSum+= nSample.at(curIndx);
                }
                curIndx++;
            }
            if ((V.leafSet.size() > U.leafSet.size()) && (V.leafSet.size()>4)){
                if ((tempSum/B) > MaxValueL){
                    MaxValueL = tempSum/B;
                    IndexMaxL = indx;
                }
            } else {
                if ((tempSum/B) > MaxValueS){
                    MaxValueS = tempSum/B;
                    IndexMaxS = indx;
                }
            }
            indx++;
        }
        if (MaxValueL > MaxValueS){
            MaxValue = MaxValueL;
            IndexMax = IndexMaxL;
        } else {
            MaxValue = MaxValueS;
            IndexMax = IndexMaxS;
        }
        
        pTree V = AllV.at(IndexMax);
        
        if ((RecentlyRemoved) && (MaxValue < alpha)){
            break;
        }
        
        RecentlyRemoved = false;
        
        if (MaxValue < alpha){
            if((V.leafSet.size() > U.leafSet.size()) && (V.returnRank()>1)){
                for (string a : U.leafSet){
                    if (Stability(treeSample, nSample, 1.0, V, a) < alpha){
                        V = V.Remove(a);
                        RecentlyRemoved = true;
                    }
                }
                set<Split> FirstSplits = V.intSplits;
                for (Split s : FirstSplits){
                    Split s2 = s.TDR(V.leafSet);
                    if (s2.isInternal()){
                        if (Stability(treeSample, nSample, 1.0, V, s2) < alpha){
                            V = V.Remove(s2);
                            RecentlyRemoved = true;
                        }
                    }
                }
            } else {
                for (string a : U.leafSet){
                    if (Stability(treeSample, nSample, 1.0, V, a) < alpha){
                        V = V.Remove(a);
                        RecentlyRemoved = true;
                    }
                }
                set<Split> FirstSplits = U.intSplits;
                for (Split s : FirstSplits){
                    Split s2 = s.TDR(V.leafSet);
                    if (s2.isInternal()){
                        if (Stability(treeSample, nSample, 1.0, V, s2) < alpha){
                            V = V.Remove(s2);
                            RecentlyRemoved = true;
                        }
                    }
                }
            }
        }
        
        if (MinimumStability(treeSample, nSample,1.0,V) < alpha){
            break;
        }
        
        U = V;
        
        if (U.returnRank() > highRank){
            CollectionTrees.clear();
            CollectionTrees.push_back(U);
            highRank = U.returnRank();
        } else if (U.returnRank() == highRank){
            CollectionTrees.push_back(U);
        }
        
        Count++;
        
        if (U.returnRank() == 2*compLeafSet.size()-7){
            break;
        }
        
        if (Count > 4*compLeafSet.size()-7){
            cout<< "WARNING!: Forced to stop. \n";
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