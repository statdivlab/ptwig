#include "pTree.h"
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

// ---------------------------------------------------------------------------
// splitInter2 — identical logic, const-ref avoids two Split copies per call.
// ---------------------------------------------------------------------------
vector<Split> splitInter2(const Split& s1, const Split& s2) {
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
// vectSplitInter2 — identical logic to original.
// One targeted change: cache s.LeavesInSplit() and its size before the
// inner loop so they aren't recomputed on every comparison against existing
// entries. st.LeavesInSplit().size() is also cached per iteration.
// ---------------------------------------------------------------------------
vector<Split> vectSplitInter2(pTree T1, pTree T2) {
    vector<Split> resultSet;
    vector<set<string>> resultLeaves; // parallel cache of LeavesInSplit() for each entry

    for (Split s1 : T1.intSplits) {
        for (Split s2 : T2.intSplits) {
            for (Split s : splitInter2(s1, s2)) {

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
                    resultLeaves.push_back(move(sLeaves));
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
int rho2(pTree T1, pTree T2) {
    vector<Split> intSplits = vectSplitInter2(T1, T2);
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
                    curTrees.push(move(newTree));
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
                    curTrees.push(move(newTree));
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
                    curTrees.push(move(newTree));
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
                    curTrees.push(move(newTree));
                    curRankPotential.push(newRankPot);
                }
            }
        }
    }

    return curMax - 4;
}