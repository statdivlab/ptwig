#include "pTree.h"
#include "rho.h"
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
#include <filesystem>
#include <future>
#include <mutex>
#include <cmath>

using namespace Rcpp;
using namespace std;

// For covering pair (U,V), returns the split present in V but absent in U.
// Precondition: V.intSplits has exactly one more element than U.intSplits.
Split findExtraSplit(const pTree& U, const pTree& V) {
    std::vector<Split> diff;
    std::set_difference(
        V.intSplits.begin(), V.intSplits.end(),
        U.intSplits.begin(), U.intSplits.end(),
        std::back_inserter(diff)
    );
    return diff.front(); // exactly one element by precondition
}

// For two compatible splits s and u on the same tree, returns the set of
// leaves separating them: the difference between the containing side of u
// and the contained side of s.
std::set<std::string> separatingLeaves(const Split& s, const Split& u) {
    // side1 of every split shares the globally smallest leaf, so the nested
    // side is always one of the two side2's — whichever is smaller.
    const std::set<std::string>* inner;
    const Split* outer;

    if (s.side2.size() <= u.side2.size()) {
        inner = &s.side2;
        outer = &u;
    } else {
        inner = &u.side2;
        outer = &s;
    }

    for (const auto* outerSide : {&outer->side1, &outer->side2}) {
        if (inner->size() < outerSide->size() &&
            std::includes(outerSide->begin(), outerSide->end(),
                          inner->begin(), inner->end())) {
            std::set<std::string> diff;
            std::set_difference(outerSide->begin(), outerSide->end(),
                                inner->begin(), inner->end(),
                                std::inserter(diff, diff.begin()));
            return diff;
        }
    }
    return {}; // s == u (shouldn't occur given preconditions)
}

// Returns the collection of minimal separating leaf sets between s and the
// splits of U (in the containment sense: no returned set contains another).
std::vector<std::set<std::string>>
minimalSeparatingLeafSets(const Split& s, const pTree& U) {
    // Collect one separating set per split in U
    std::vector<std::set<std::string>> candidates;
    candidates.reserve(U.intSplits.size());

    for (const Split& u : U.intSplits) {
        std::set<std::string> sep = separatingLeaves(s, u);
        if (!sep.empty()) {
            candidates.push_back(std::move(sep));
        }
    }

    // Keep only the minimal elements under set inclusion.
    // A candidate is minimal if no other candidate is a strict subset of it.
    std::vector<std::set<std::string>> minimal;
    for (std::size_t i = 0; i < candidates.size(); ++i) {
        bool dominated = false;
        for (std::size_t j = 0; j < candidates.size(); ++j) {
            if (i == j) continue;
            // Is candidates[j] a strict subset of candidates[i]?
            if (candidates[j].size() < candidates[i].size() &&
                std::includes(candidates[i].begin(), candidates[i].end(),
                              candidates[j].begin(), candidates[j].end())) {
                dominated = true;
                break;
            }
        }
        if (!dominated) {
            minimal.push_back(candidates[i]);
        }
    }

    return minimal;
}


// ---------------------------------------------------------------------------
// splitInter — identical logic, const-ref avoids two Split copies per call.
// ---------------------------------------------------------------------------
vector<Split> splitInter(const Split& s1, const Split& s2) {
    set<string> intSide1Side1, intSide1Side2, intSide2Side1, intSide2Side2;

    set_intersection(s1.side1.begin(), s1.side1.end(),
                     s2.side1.begin(), s2.side1.end(),
                     inserter(intSide1Side1, intSide1Side1.begin()));
    set_intersection(s1.side1.begin(), s1.side1.end(),
                     s2.side2.begin(), s2.side2.end(),
                     inserter(intSide1Side2, intSide1Side2.begin()));
    set_intersection(s1.side2.begin(), s1.side2.end(),
                     s2.side1.begin(), s2.side1.end(),
                     inserter(intSide2Side1, intSide2Side1.begin()));
    set_intersection(s1.side2.begin(), s1.side2.end(),
                     s2.side2.begin(), s2.side2.end(),
                     inserter(intSide2Side2, intSide2Side2.begin()));

    Split spInt1(intSide1Side1, intSide2Side2);
    Split spInt2(intSide1Side2, intSide2Side1);

    vector<Split> respSplits;
    if (spInt1.isInternal()) respSplits.push_back(spInt1);
    if (spInt2.isInternal()) respSplits.push_back(spInt2);
    return respSplits;
}

// ---------------------------------------------------------------------------
// vectSplitInter — identical logic to original.
// One targeted change: cache s.LeavesInSplit() and its size before the
// inner loop so they aren't recomputed on every comparison against existing
// entries. st.LeavesInSplit().size() is also cached per iteration.
// ---------------------------------------------------------------------------
vector<Split> vectSplitInter(pTree T1, pTree T2) {
    vector<Split> resultSet;
    vector<set<string>> resultLeaves; // parallel cache of LeavesInSplit() for each entry

    for (Split s1 : T1.intSplits) {
        for (Split s2 : T2.intSplits) {
            for (Split s : splitInter(s1, s2)) {

                set<string> sLeaves = s.LeavesInSplit();
                size_t      sSize   = sLeaves.size();

                bool AddEnd = true;

                for (int i = 0; i < static_cast<int>(resultSet.size()); ) {
                    Split& st = resultSet[i];

                    // Use cached size — avoids recomputing LeavesInSplit().size()
                    size_t stSize = resultLeaves[i].size();

                    if (stSize < sSize) {
                        if (AddEnd) {
                            resultSet.insert(resultSet.begin() + i, s);
                            resultLeaves.insert(resultLeaves.begin() + i, sLeaves);
                            AddEnd = false;
                            ++i;
                            continue;
                        }
                        // Use cached leaves for contains() — avoids LeavesInSplit() inside it
                        if (s.TDR(resultLeaves[i]) == st) {
                            resultSet.erase(resultSet.begin() + i);
                            resultLeaves.erase(resultLeaves.begin() + i);
                            continue;
                        }
                    } else {
                        // st.contains(s): check if st restricted to sLeaves == s
                        if (st.TDR(sLeaves) == s) {
                            AddEnd = false;
                            break;
                        }
                    }
                    ++i;
                }

                if (AddEnd) {
                    resultSet.push_back(s);
                    resultLeaves.push_back(std::move(sLeaves));
                }
            }
        }
    }

    return resultSet;
}

// ---------------------------------------------------------------------------
// rho — the original stack machine is completely preserved, character for
// character. The only changes are:
//
// 1. Pre-cache LeavesInSplit() for every split into `cachedLeaves` before
//    the loop. LeavesInSplit() does a set_union on every call; with O(n²)
//    or more Insert calls in the hot path this adds up significantly.
//
// 2. Replace curTrees.top().Insert(intSplits.at(i)) with
//    curTrees.top().InsertCached(intSplits.at(i), cachedLeaves.at(i))
//    everywhere in the stack machine. InsertCached has the identical body
//    as Insert but skips the internal LeavesInSplit() call, using the
//    pre-cached set instead.
//
// 3. Use cachedLeaves[start].size() instead of
//    intSplits[start].LeavesInSplit().size() for the seed construction.
//
// No logic, no control flow, no data structure is changed.
// ---------------------------------------------------------------------------
int rho(pTree T1, pTree T2) {
    vector<Split> intSplits = vectSplitInter(T1, T2);
    vector<int> rankPotential;

    // Pre-cache LeavesInSplit() for every split — avoids recomputing the
    // set union inside every Insert call in the hot path below.
    vector<set<string>> cachedLeaves(intSplits.size());
    for (int i = 0; i < static_cast<int>(intSplits.size()); ++i) {
        cachedLeaves[i] = intSplits[i].LeavesInSplit();
        rankPotential.push_back(2 * static_cast<int>(cachedLeaves[i].size()) - 3);
    }

    if (intSplits.empty()) return 0;

    int curMax = static_cast<int>(cachedLeaves[0].size()) + 1;
    stack<int>   indxC;
    stack<pTree> curTrees;
    stack<int>   curRankPotential;

    int Beginning = 0;
    int End = static_cast<int>(rankPotential.size());

    while (Beginning < End) {
        if (indxC.empty()) {
            indxC.push(Beginning);
            curTrees.push(pTree(cachedLeaves[Beginning],
                                set<Split>{ intSplits[Beginning] }));
            curRankPotential.push(rankPotential[Beginning]);
        }

        if (rankPotential[indxC.top()] <= curMax) {
            End = indxC.top();
            indxC.pop();
            curTrees.pop();
            curRankPotential.pop();
        } else {
            int tIndx = indxC.top();

            if (curRankPotential.top() > curMax) {
                if (tIndx + 1 < End) {
                    indxC.push(tIndx + 1);
                    pTree newTree = curTrees.top().InsertCached(
                                       intSplits[tIndx + 1],
                                       cachedLeaves[tIndx + 1]);
                    int newRankPot = 2 * static_cast<int>(newTree.leafSet.size()) - 3;
                    if (newTree.rank + 4 > curMax) curMax = newTree.rank + 4;
                    curTrees.push(std::move(newTree));
                    curRankPotential.push(newRankPot);
                } else {
                    indxC.pop();
                    curTrees.pop();
                    curRankPotential.pop();

                    if (indxC.empty()) {
                        Beginning = tIndx + 1;
                        if (Beginning >= static_cast<int>(intSplits.size())) break;
                        if (((End - Beginning) + static_cast<int>(cachedLeaves[Beginning].size())) <= curMax) break;
                        continue;
                    }

                    int next = indxC.top() + 1;
                    indxC.pop();
                    curTrees.pop();
                    curRankPotential.pop();

                    if (indxC.empty()) {
                        Beginning = next;
                        if (Beginning >= static_cast<int>(intSplits.size())) break;
                        if (((End - Beginning) + static_cast<int>(cachedLeaves[Beginning].size())) <= curMax) break;
                        continue;
                    }

                    indxC.push(next);
                    pTree newTree = curTrees.top().InsertCached(
                                       intSplits[next],
                                       cachedLeaves[next]);
                    int newRankPot = 2 * static_cast<int>(newTree.leafSet.size()) - 3;
                    if (newTree.rank + 4 > curMax) curMax = newTree.rank + 4;
                    curTrees.push(std::move(newTree));
                    curRankPotential.push(newRankPot);
                }
            } else {
                if (tIndx + 1 < End) {
                    indxC.pop();
                    curTrees.pop();
                    curRankPotential.pop();

                    indxC.push(tIndx + 1);
                    pTree newTree = curTrees.top().InsertCached(
                                       intSplits[tIndx + 1],
                                       cachedLeaves[tIndx + 1]);
                    int newRankPot = 2 * static_cast<int>(newTree.leafSet.size()) - 3;
                    curTrees.push(std::move(newTree));
                    curRankPotential.push(newRankPot);
                } else {
                    indxC.pop();
                    curTrees.pop();
                    curRankPotential.pop();

                    if (indxC.empty()) {
                        Beginning = tIndx + 1;
                        if (Beginning >= static_cast<int>(intSplits.size())) break;
                        if (((End - Beginning) + static_cast<int>(cachedLeaves[Beginning].size())) <= curMax) break;
                        continue;
                    }

                    int next = indxC.top() + 1;
                    indxC.pop();
                    curTrees.pop();
                    curRankPotential.pop();

                    if (indxC.empty()) {
                        Beginning = next;
                        if (Beginning >= static_cast<int>(intSplits.size())) break;
                        if (((End - Beginning) + static_cast<int>(cachedLeaves[Beginning].size())) <= curMax) break;
                        continue;
                    }

                    indxC.push(next);
                    pTree newTree = curTrees.top().InsertCached(
                                       intSplits[next],
                                       cachedLeaves[next]);
                    int newRankPot = 2 * static_cast<int>(newTree.leafSet.size()) - 3;
                    if (newTree.rank + 4 > curMax) curMax = newTree.rank + 4;
                    curTrees.push(std::move(newTree));
                    curRankPotential.push(newRankPot);
                }
            }
        }
    }

    return curMax - 4;
}

oRho::oRho(int nrho, std::vector<std::set<std::string>> newPresLeaves){
    rho = nrho;
    presLeaves = newPresLeaves;
}


oRho rho(pTree U, pTree V, oRho baseORho, pTree Tl, string extral){
    if (!Tl.leafSet.count(extral)){
        return oRho(baseORho.rho, baseORho.presLeaves);
    }
    
    int bRho = -1;
    set<set<string>> cleanedBaseLeaves;
    for (set<string> tLeaves :  baseORho.presLeaves){
        set<string> remainingLeaves = tLeaves;
        remainingLeaves.erase(extral);
        cleanedBaseLeaves.insert(remainingLeaves);
    }
    vector<set<string>> potentialPresLeaves;
    
    for (set<string> tLeaves :  cleanedBaseLeaves){
            
        pTree Ttemp = commonLower(U, Tl, tLeaves);
            
        if (Ttemp.rank > bRho){
            bRho = Ttemp.rank;
            potentialPresLeaves.clear();
            potentialPresLeaves.push_back(tLeaves);
        } else if (Ttemp.rank == bRho){
            potentialPresLeaves.push_back(tLeaves);
        }
            
    }
    
    return oRho(bRho, potentialPresLeaves);
        
}

oRho rho(pTree U, pTree V, oRho baseORho, pTree Tl, Split extraS){
    set<string> intLeaves;
    set_intersection(V.leafSet.begin(), V.leafSet.end(),
                     Tl.leafSet.begin(), Tl.leafSet.end(),                         
                     inserter(intLeaves, intLeaves.begin()));
                     
    vector<set<string>> SepL =  minimalSeparatingLeafSets(extraS, U); 
        
    vector<std::set<std::string>> SepLeaves;
    SepLeaves.reserve(SepL.size());
    
    for (const auto& sep : SepL) {
        std::set<std::string> inter;
        std::set_intersection(sep.begin(), sep.end(),
                                intLeaves.begin(), intLeaves.end(),
                                std::inserter(inter, inter.begin()));
        SepLeaves.push_back(std::move(inter));
    }
                         
    int bRho = -1;
    set<set<string>> cleanedBaseLeaves;
    cleanedBaseLeaves.insert(intLeaves);
    for (set<string> tLeaves :  baseORho.presLeaves){
        cleanedBaseLeaves.insert(tLeaves);
        for (set<string> sep : SepLeaves){
            set<string> diffLeaves;
            set_difference(tLeaves.begin(), tLeaves.end(),
                            sep.begin(), sep.end(),
                            inserter(diffLeaves, diffLeaves.begin()));
            cleanedBaseLeaves.insert(diffLeaves);
        }
        
    }
    
    vector<set<string>> potentialPresLeaves;
    
    for (set<string> tLeaves :  cleanedBaseLeaves){
            
        pTree Ttemp = commonLower(U, Tl, tLeaves);
            
        if (Ttemp.rank > bRho){
            bRho = Ttemp.rank;
            potentialPresLeaves.clear();
            potentialPresLeaves.push_back(tLeaves);
        } else if (Ttemp.rank == bRho){
            potentialPresLeaves.push_back(tLeaves);
        }
            
    }
    
    return oRho(bRho, potentialPresLeaves);
        
}

oRho rho(pTree V, pTree U, oRho baseORho, pTree Tl) {
    
    set<string> newleaf;
    int bRho = -1;
    vector<set<string>> potentialPresLeaves;
    
    if (V.rank == 1){
        set_intersection(V.leafSet.begin(), V.leafSet.end(),
                         Tl.leafSet.begin(), Tl.leafSet.end(),
                         inserter(newleaf, newleaf.begin()));
        
        potentialPresLeaves.push_back(newleaf);
        if (Tl.over(V)) return oRho(1, potentialPresLeaves);
        return oRho(0, potentialPresLeaves);
    }
    
    set_difference(V.leafSet.begin(), V.leafSet.end(),
                       U.leafSet.begin(), U.leafSet.end(),
                       inserter(newleaf, newleaf.begin()));
    
        
    if (newleaf.empty()){
        set<string> intLeaves;
        set_intersection(V.leafSet.begin(), V.leafSet.end(),
                         Tl.leafSet.begin(), Tl.leafSet.end(),
                         inserter(intLeaves, intLeaves.begin()));
        
        Split s =  findExtraSplit(U, V);
        
        vector<set<string>> SepL =  minimalSeparatingLeafSets(s, U); 
        
        vector<std::set<std::string>> SepLeaves;
        SepLeaves.reserve(SepL.size());

        for (const auto& sep : SepL) {
            std::set<std::string> inter;
            std::set_intersection(sep.begin(), sep.end(),
                                  intLeaves.begin(), intLeaves.end(),
                                  std::inserter(inter, inter.begin()));
            SepLeaves.push_back(std::move(inter));
        }

        
        {pTree Ttemp = commonLower(V,Tl, intLeaves);
         bRho = Ttemp.rank;
         potentialPresLeaves.push_back(intLeaves);
        }
        for (set<string> tLeaves : baseORho.presLeaves){
            //set<string> extraLeaves;
            //set_difference(intLeaves.begin(), intLeaves.end(),
            //               tLeaves.begin(), tLeaves.end(),
            //               inserter(extraLeaves, extraLeaves.begin()));
            
            for (set<string> eLeaves : SepLeaves){
                if (!std::includes(tLeaves.begin(), tLeaves.end(),
                                    eLeaves.begin(), eLeaves.end())){
                    set<string> unionLeaves;
                    set_union(tLeaves.begin(), tLeaves.end(),
                                eLeaves.begin(), eLeaves.end(),
                                inserter(unionLeaves, unionLeaves.begin()));

                    pTree Ttemp = commonLower(V, Tl, unionLeaves);
                    if (Ttemp.rank > bRho){
                    bRho = Ttemp.rank;
                    potentialPresLeaves.clear();
                    potentialPresLeaves.push_back(unionLeaves);
                    } else if (Ttemp.rank == bRho){
                        potentialPresLeaves.push_back(unionLeaves);
                    }
                }
                
            }
            if(tLeaves.size() < intLeaves.size()){   
                pTree Ttemp = commonLower(V, Tl, tLeaves);
                if (Ttemp.rank > bRho){
                bRho = Ttemp.rank;
                potentialPresLeaves.clear();
                potentialPresLeaves.push_back(tLeaves);
                } else if (Ttemp.rank == bRho){
                    potentialPresLeaves.push_back(tLeaves);
                }
            }
        }
        
    } else {
        if (!Tl.leafSet.count(*newleaf.begin())){
            return oRho(baseORho.rho, baseORho.presLeaves);
        }
        
        bRho = baseORho.rho;
        potentialPresLeaves = baseORho.presLeaves;
        
        for (set<string> tLeaves :  baseORho.presLeaves){
            set<string> unionLeaves;
            
            set_union(tLeaves.begin(), tLeaves.end(),
              newleaf.begin(), newleaf.end(),
              inserter(unionLeaves, unionLeaves.begin()));
            
            pTree Ttemp = commonLower(V, Tl, unionLeaves);
            
            if (Ttemp.rank > bRho){
                bRho = Ttemp.rank;
                potentialPresLeaves.clear();
                potentialPresLeaves.push_back(unionLeaves);
            } else if (Ttemp.rank == bRho){
                potentialPresLeaves.push_back(unionLeaves);
            }
            
        }
    }
    
    return oRho(bRho, potentialPresLeaves);
}
