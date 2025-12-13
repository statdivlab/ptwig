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
#include <filesystem> // C++17
#include <future>
#include <mutex>
#include <cmath>

using namespace Rcpp;
using namespace std;

vector<Split> splitInter(Split s1, Split s2){
    set<string> intSide1Side1;
    set<string> intSide1Side2;
    set<string> intSide2Side1;
    set<string> intSide2Side2;
    
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
    
    Split spInt1 = Split(intSide1Side1, intSide2Side2);
    Split spInt2 = Split(intSide1Side2, intSide2Side1);
    
    vector<Split> respSplits;
    if (spInt1.isInternal()){
        respSplits.push_back(spInt1);
    }
    
    if (spInt2.isInternal()){
        respSplits.push_back(spInt2);
    }
    
    return respSplits;
};

vector<Split> vectSplitInter(pTree T1, pTree T2){
    vector<Split> resultSet;
    
    for (Split s1 : T1.intSplits){
        for (Split s2 : T2.intSplits){
            for (Split s : splitInter(s1,s2)){
                bool AddEnd = true;
               
                for (auto Indx = resultSet.begin(); Indx != resultSet.end(); ){
                    Split &st = *Indx;
                    if (st.LeavesInSplit().size() < s.LeavesInSplit().size()){
                        if (AddEnd){
                            Indx = resultSet.insert(Indx, s);
                            AddEnd = false;
                            ++Indx; // move past the inserted element
                            continue;
                        }
                        if (s.contains(st)){
                            Indx = resultSet.erase(Indx);
                            continue;
                        }
                    } else if (st.contains(s)){
                        AddEnd = false;
                        break;
                    }
                    ++Indx;
                }
                if (AddEnd){
                    resultSet.push_back(s);
                }
            }
        }
    }
    return resultSet;
}

int rho(pTree T1, pTree T2){
    vector<Split> intSplits = vectSplitInter(T1, T2);
    vector<int> rankPotential;
    
    for (Split s : intSplits){
        rankPotential.push_back(2*(static_cast<int>(s.LeavesInSplit().size())) - 3);
    }
    
    if (intSplits.empty()){
        return (0);
    }
    
    int curMax = (static_cast<int>(intSplits.at(0).LeavesInSplit().size())) + 1;
    stack<int> indxC;
    stack<pTree> curTrees;
    stack<int> curRankPotential;
    
    int Beginning = 0;
    int End = static_cast<int>(rankPotential.size());
    while(Beginning < End){
        if (indxC.empty()){
            indxC.push(Beginning);
            curTrees.push(pTree(intSplits.at(Beginning).LeavesInSplit(),set<Split> {intSplits.at(Beginning)}));
            curRankPotential.push(rankPotential.at(Beginning));
        }
        if (rankPotential.at(indxC.top()) <= curMax){
            End = indxC.top();
            indxC.pop();
            curTrees.pop();
            curRankPotential.pop();
        } else {
            int tIndx = indxC.top();
            if (curRankPotential.top() > curMax){
                if (tIndx + 1 < End){
                    indxC.push(tIndx+1);
                    pTree newTree = curTrees.top().Insert(intSplits.at(tIndx+1));
                    curTrees.push(newTree);
                    curRankPotential.push(2*(static_cast<int>(newTree.leafSet.size())) -3);
                    
                    if ((newTree.rank+4) > curMax){
                        curMax = (newTree.rank + 4);
                    }
                } else {
                    indxC.pop();
                    curTrees.pop();
                    curRankPotential.pop();
                    
                    if (indxC.empty()){
                        Beginning = tIndx + 1;
                        if (Beginning >= intSplits.size()){
                            break;
                        }
                        if (((End - Beginning) + static_cast<int>(intSplits.at(Beginning).LeavesInSplit().size())) <= curMax){
                            break;
                        }
                        continue;
                    }
                    
                    int next = indxC.top() + 1;
                    indxC.pop();
                    curTrees.pop();
                    curRankPotential.pop();
                    
                    if (indxC.empty()){
                        Beginning = next;
                        if (Beginning >= intSplits.size()){
                            break;
                        }
                        if (((End - Beginning) + static_cast<int>(intSplits.at(Beginning).LeavesInSplit().size())) <= curMax){
                            break;
                        }
                        continue;
                    }
                    
                    indxC.push(next);
                    pTree newTree = curTrees.top().Insert(intSplits.at(next));
                    curTrees.push(newTree);
                    curRankPotential.push(2*(static_cast<int>(newTree.leafSet.size())) -3);
                    
                    if ((newTree.rank + 4) > curMax){
                        curMax = newTree.rank + 4;
                    }
                }
            } else {
                if (tIndx + 1 < End){
                    indxC.pop();
                    curTrees.pop();
                    curRankPotential.pop();
                    
                    indxC.push(tIndx+1);
                    pTree newTree = curTrees.top().Insert(intSplits.at(tIndx+1));
                    curTrees.push(newTree);
                    curRankPotential.push(2*(static_cast<int>(newTree.leafSet.size())) -3);
                } else {
                    indxC.pop();
                    curTrees.pop();
                    curRankPotential.pop();
                    
                    if (indxC.empty()){
                        Beginning = tIndx + 1;
                        if (Beginning >= intSplits.size()){
                            break;
                        }
                        if (((End - Beginning) + static_cast<int>(intSplits.at(Beginning).LeavesInSplit().size())) <= curMax){
                            break;
                        }
                        continue;
                    }
                    
                    int next = indxC.top() + 1;
                    indxC.pop();
                    curTrees.pop();
                    curRankPotential.pop();
                    
                    if (indxC.empty()){
                        Beginning = next;
                        if (Beginning >= intSplits.size()){
                            break;
                        }
                        if (((End - Beginning) + static_cast<int>(intSplits.at(Beginning).LeavesInSplit().size())) <= curMax){
                            break;
                        }
                        continue;
                    }
                    
                    
                    indxC.push(next);
                    pTree newTree = curTrees.top().Insert(intSplits.at(next));
                    curTrees.push(newTree);
                    curRankPotential.push(2*(static_cast<int>(newTree.leafSet.size())) -3);
                    
                    if ((newTree.rank+4) > curMax){
                        curMax = (newTree.rank + 4);
                    }
                }
            }
        }
    }
    return curMax - 4;
}