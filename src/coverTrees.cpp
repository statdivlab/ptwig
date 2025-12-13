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

bool compareSplitsBySide2(const Split& a, const Split& b){
    return a.side2.size() > b.side2.size();
}

bool precepSplitsBySide2(const Split& a, const Split& b){
    set<string> Intercept;
    
    set_intersection(a.side2.begin(), a.side2.end(),
                     b.side2.begin(), b.side2.end(),
                     inserter(Intercept, Intercept.begin()));
    
    return Intercept == b.side2;
}

class SplitNode{
    public:
        Split spt;
        int parent;
        vector<int> children;
    
        SplitNode (Split nSpt){
            spt = nSpt;
            parent = -1;
        }
    
        void setParent(int nPar){
            parent = nPar;
        }
    
        void insertChild(int nCh){
            children.push_back(nCh);
        }
    
        void removeChild(int oCh){
            children.erase(find(children.begin(), children.end(), oCh));
        }
};

vector<set<string>> getAllSubsets(const set<string>& inputSet, int sizemin, int sizemax) {
    vector<string> elements(inputSet.begin(), inputSet.end());
    int n = static_cast<int>(elements.size());
    vector<set<string>> powerSet;

    // Loop through all possible subsets (2^n)
    for (int mask = 0; mask < (1 << n); ++mask) {
        set<string> subset;
        for (int i = 0; i < n; ++i) {
            if (mask & (1 << i)) {  // Check if the i-th element should be included
                subset.insert(elements[i]);
            }
        }
        if ((subset.size() >= sizemin) && (subset.size() <= sizemax)){
            powerSet.push_back(subset);
        }
    }

    return powerSet;
}

vector<vector<int>> getIndicesSubsets(int n, int sizemin, int sizemax){
    vector<vector<int>> powerSet;
    
    for (int mask = 0; mask < (1 << n); ++mask){
        vector<int> subset;
        for (int i = 0; i < n; ++i){
            if (mask & (1 << i)){
                subset.push_back(i);
            }
        }
        if ((subset.size() >= sizemin) && (subset.size() <= sizemax)){
            powerSet.push_back(subset);
        }
    }
    
    return powerSet;
}

vector<pTree> coverTrees(pTree T, set<string> compLeaves){
    
    vector<pTree> resultTrees;
    
    if (T.intSplits.size() == 0){
        vector<string> AllLeaves(compLeaves.begin(), compLeaves.end());
        
        for (int i1 = 0; i1 < AllLeaves.size(); i1++){
            for (int i2 = i1 + 1; i2 < AllLeaves.size();  i2++){
                set<string> newSide1 = {AllLeaves.at(i1),AllLeaves.at(i2)};
                for (int j1 = i1 + 1; j1 < AllLeaves.size(); j1++){
                    if (j1 != i2){
                        for (int j2 = j1 + 1; j2 < AllLeaves.size(); j2++){
                            if (j2 != i2){
                                set<string> newSide2 = {AllLeaves.at(j1),AllLeaves.at(j2)};
                                set<Split> NewSplits;
                                NewSplits.insert(Split(newSide1, newSide2));
                                resultTrees.push_back(pTree(set<string>{AllLeaves.at(i1),AllLeaves.at(i2),AllLeaves.at(j1),AllLeaves.at(j2)},
                                                            NewSplits));
                            }
                        }
                    }
                }
            }
        }
    } else {
        vector<Split> orderedVector(T.intSplits.begin(), T.intSplits.end());
        std::sort(orderedVector.begin(), orderedVector.end(), compareSplitsBySide2);
        
        vector<int> UpperSplitsIndices;
        
        set<string> newLeaves;
        
        set_difference(compLeaves.begin(), compLeaves.end(),
                       T.leafSet.begin(), T.leafSet.end(),
                       inserter(newLeaves, newLeaves.begin()));
        
        vector<SplitNode> SplitsOrdered;
        SplitsOrdered.push_back(SplitNode(orderedVector.at(0)));
        UpperSplitsIndices.push_back(0);
        if (orderedVector.size() > 1){
            for (int j = 1; j < orderedVector.size(); j++){
                SplitsOrdered.push_back(SplitNode(orderedVector.at(j)));
                bool FindingPosition = true;
                int k = 0;
                while (FindingPosition){
                    if (!precepSplitsBySide2(SplitsOrdered.at(k).spt,SplitsOrdered.at(SplitsOrdered.size()-1).spt)){
                        k++;
                        if (k == SplitsOrdered.size()-1){
                            FindingPosition = false;
                            UpperSplitsIndices.push_back(k);
                        }
                        continue;
                    }
                    if (SplitsOrdered.at(k).children.size()>0){
                        bool SomeChild = false;
                        for (int l = 0; l < SplitsOrdered.at(k).children.size(); l++){
                            int k2 = SplitsOrdered.at(k).children.at(l);
                            if (precepSplitsBySide2(SplitsOrdered.at(k2).spt, SplitsOrdered.at(SplitsOrdered.size()-1).spt)){
                                k = k2;
                                SomeChild = true;
                                break;
                            } else if (precepSplitsBySide2(SplitsOrdered.at(SplitsOrdered.size()-1).spt, SplitsOrdered.at(k2).spt)){
                                SplitsOrdered.at(k).removeChild(k2);
                                SplitsOrdered.at(k).insertChild(static_cast<int>(SplitsOrdered.size())-1);
                                SplitsOrdered.at(static_cast<int>(SplitsOrdered.size())-1).setParent(k);
                                SplitsOrdered.at(static_cast<int>(SplitsOrdered.size())-1).insertChild(k2);
                                SplitsOrdered.at(k2).setParent(static_cast<int>(SplitsOrdered.size())-1);
                                FindingPosition = false;
                                SomeChild = true;
                                break;
                            }
                        }
                        if (!SomeChild){
                            SplitsOrdered.at(k).insertChild(static_cast<int>(SplitsOrdered.size())-1);
                            SplitsOrdered.at(static_cast<int>(SplitsOrdered.size())-1).setParent(k);
                            FindingPosition = false;
                        }
                    } else {
                        SplitsOrdered.at(k).insertChild(static_cast<int>(SplitsOrdered.size())-1);
                        SplitsOrdered.at(static_cast<int>(SplitsOrdered.size())-1).setParent(k);
                        FindingPosition = false;
                    }
                }
            }
        }
        
        bool AddLeafScheme[static_cast<int>(T.intSplits.size()) +1][static_cast<int>(T.intSplits.size())] ;
        
        for (int rn = 0; rn < static_cast<int>(T.intSplits.size()) +1; rn++){
            for (int cn = 0; cn < static_cast<int>(T.intSplits.size()); cn++) {
                AddLeafScheme[rn][cn] = false;
            }
        }
        
        //Dealing with "root" node in the tree (r)
        
        set<string> UpperNodeLeaves = SplitsOrdered.at(UpperSplitsIndices.at(0)).spt.side1;
        for (int indx : UpperSplitsIndices){
            set<string> tempUpperNodeLeaves;
            set_intersection(SplitsOrdered.at(indx).spt.side1.begin(), SplitsOrdered.at(indx).spt.side1.end(),
                             UpperNodeLeaves.begin(), UpperNodeLeaves.end(),
                             inserter(tempUpperNodeLeaves, tempUpperNodeLeaves.begin()));
            UpperNodeLeaves = tempUpperNodeLeaves;
        }
        
        UpperNodeLeaves.erase(UpperNodeLeaves.begin());
        
        vector<set<string>> UpperSubsets = getAllSubsets(UpperNodeLeaves, 0, static_cast<int>(UpperNodeLeaves.size()));
        vector<vector<int>> UpperSplitsSubsets =  getIndicesSubsets(static_cast<int>(UpperSplitsIndices.size()),0, static_cast<int>(UpperSplitsIndices.size()));
        
        for (set<string> LooseLeaves : UpperSubsets){
    
            for(vector<int> indxSubsets : UpperSplitsSubsets){
    
                if (((LooseLeaves.size() + indxSubsets.size()) > 1) && ((LooseLeaves.size() + indxSubsets.size()) < UpperNodeLeaves.size() + UpperSplitsIndices.size())){
                    set<string> side2 = LooseLeaves;
                    for (int indx : indxSubsets){
                        set<string> tempSide2;
                        set_union(SplitsOrdered.at(UpperSplitsIndices.at(indx)).spt.side2.begin(),
                                  SplitsOrdered.at(UpperSplitsIndices.at(indx)).spt.side2.end(),
                                  side2.begin(), side2.end(), inserter(tempSide2, tempSide2.begin()));
                        side2 = tempSide2;
                    }
                    set<string> side1;
                    set_difference(T.leafSet.begin(), T.leafSet.end(),
                                   side2.begin(), side2.end(), inserter(side1, side1.begin()));
                    
                    set<Split> UnionSplitsUpper = T.intSplits;
    
                    UnionSplitsUpper.insert(Split(side1,side2));
                    resultTrees.push_back(pTree(T.leafSet,UnionSplitsUpper));
                }
            }
        }
        
        int count = 0;
        while (count < SplitsOrdered.size()){
            int parIndex = SplitsOrdered.at(count).parent;
            AddLeafScheme[count][count] = true;
            
            while(parIndex != -1){
                AddLeafScheme[count][parIndex] = true;
                parIndex = SplitsOrdered.at(parIndex).parent;
            }
            
            set<string> LeavesAtNode = SplitsOrdered.at(count).spt.side2;
            
            for (int indx : SplitsOrdered.at(count).children){
                set<string> tempLeavesAtNode;
                set_difference(LeavesAtNode.begin(), LeavesAtNode.end(),
                               SplitsOrdered.at(indx).spt.side2.begin(),
                               SplitsOrdered.at(indx).spt.side2.end(),
                               inserter(tempLeavesAtNode, tempLeavesAtNode.begin()));
                LeavesAtNode = tempLeavesAtNode;
            }
            
            vector<set<string>> LeavesSubsets = getAllSubsets(LeavesAtNode, 0, static_cast<int>(LeavesAtNode.size()));
            vector<vector<int>> SplitsSubsets =  getIndicesSubsets(static_cast<int>(SplitsOrdered.at(count).children.size()),0, static_cast<int>(SplitsOrdered.at(count).children.size()));
            
            for (set<string> LooseLeaves : LeavesSubsets){
                for(vector<int> indxSubsets : SplitsSubsets){
                    if (((LooseLeaves.size() + indxSubsets.size()) > 1) && ((LooseLeaves.size() + indxSubsets.size()) < LeavesAtNode.size() + SplitsOrdered.at(count).children.size())){
                        
                        set<string> side2 = LooseLeaves;
                        for (int indx : indxSubsets){
                            set<string> tempSide2;
                            set_union(SplitsOrdered.at(SplitsOrdered.at(count).children.at(indx)).spt.side2.begin(),
                                      SplitsOrdered.at(SplitsOrdered.at(count).children.at(indx)).spt.side2.end(),
                                      side2.begin(), side2.end(), inserter(tempSide2, tempSide2.begin()));
                            side2 = tempSide2;
                        }
                        set<string> side1;
                        set_difference(T.leafSet.begin(), T.leafSet.end(),
                                       side2.begin(), side2.end(), inserter(side1, side1.begin()));
                        
                        set<Split> UnionSplits = T.intSplits;
    
                        UnionSplits.insert(Split(side1,side2));
                        resultTrees.push_back(pTree(T.leafSet,UnionSplits));
                    }
                }
            }
            count++;
        }
        
        for (string nleaf : newLeaves){
            for (int rn = 0 ; rn < static_cast<int>(T.intSplits.size()) + 1; rn++){
                set<Split> NewSplits;
                int count = 0;
                for (SplitNode nSpt : SplitsOrdered){
                    Split tempSpt = nSpt.spt;
                    if (AddLeafScheme[rn][count]){
                        tempSpt.addLeaf(2, nleaf);
                    } else {
                        tempSpt.addLeaf(1, nleaf);
                    }
                    NewSplits.insert(Split(tempSpt.side1,tempSpt.side2));
                    count++;
                }
                set<string> LeafSetNew = T.leafSet;
                LeafSetNew.insert(nleaf);
                resultTrees.push_back(pTree(LeafSetNew,NewSplits));
            }
        }
    }
    
    return resultTrees;
}