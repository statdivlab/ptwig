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


// =====================
// Split implementations
// =====================

Split::Split(){
    side1 = {};
    side2 = {};
}
    
Split::Split(set<string> inputSet1, set<string> inputSet2){
    set<string> intersection_result;
    set_intersection(inputSet1.begin(), inputSet1.end(),
                     inputSet2.begin(), inputSet2.end(),
                     inserter(intersection_result, intersection_result.begin()));

    if (intersection_result.empty()){
        if (inputSet1 < inputSet2){
            side1 = inputSet1;
            side2 = inputSet2;
        } else {
            side1 = inputSet2;
            side2 = inputSet1;
        }
    } else {
        side1 = {};
        side2 = {};
        cout << "Not an actual split. It keeps it empty \n";
    }
}
    
bool Split::operator<(const Split& other) const {//Operation that defines which edge goes first in order
    if (side1 < other.side1){
        return true;
    } else if (side1 == other.side1){
        return side2 < other.side2;
    }
    return false;
}
    
bool Split::operator==(const Split& other) const{
    return (((side1 == other.side1) && (side2 == other.side2)) || ((side1 == other.side2) && (side2 == other.side1)));
}
    
void Split::addLeaf(int side, string nLeaf){
    if (side == 1){
        side1.insert(nLeaf);
    } else if (side == 2){
        side2.insert(nLeaf);
    }
}
    
const void Split::print(){
    cout << "{";
    bool first = true;
    for (string leaf : side1){
        if(!first) cout << ", ";
        cout << leaf;
        first = false;
    }
    cout << "|";
    first = true;
    for (string leaf : side2){
        if(!first) cout << ", ";
        cout << leaf;
        first = false;
    }
    cout << "} \n";
}

string Split::printSt() const{
    string Resp = "{";
    bool first = true;
    for (string leaf : side1){
        if (!first) Resp = Resp + ", ";
        Resp = Resp + leaf;
        first = false;
    }
    Resp = Resp + "|";
    first = true;
    for (string leaf : side2){
        if (!first) Resp = Resp + ", ";
        Resp = Resp + leaf;
        first = false;
    }
    Resp = Resp + "}";

    return Resp;
}

Split Split::TDR(set<string> L){
    set<string> intersection1;
    set<string> intersection2;

    set_intersection(side1.begin(), side1.end(),
                     L.begin(), L.end(),
                     inserter(intersection1, intersection1.begin()));

    set_intersection(side2.begin(), side2.end(),
                     L.begin(), L.end(),
                     inserter(intersection2, intersection2.begin()));

    return Split(intersection1, intersection2);
}

bool Split::isInternal(){
    return ((side1.size() > 1) && (side2.size() > 1));
}

set<string> Split::LeavesInSplit(){
    set<string> UnionLeaves;
    set_union(side1.begin(),side1.end(),
              side2.begin(),side2.end(),
              inserter(UnionLeaves, UnionLeaves.begin()));

    return UnionLeaves;
}

bool Split::contains(Split otherS){
    return ( TDR(otherS.LeavesInSplit()) == otherS); //Checks if the splits coincide in the presummably smaller set of leaves in otherS.
}


// =====================
// pTree implementations
// =====================
    
pTree::pTree(){
    leafSet = set<string>();
    intSplits = set<Split>();
    complexity = 0;
    rank = 0;
}

pTree::pTree(string newick){
    string workingNewick = newick.substr(1,newick.find_last_of(')')-1);

    vector<Split> runningSplits1;
    vector<Split> runningSplits2;
    Split tempSplit;

    int i = 0;
    while((i < workingNewick.size()) && (i > -1)){
        switch (workingNewick[i]) {
            case '(':
                
                tempSplit = Split();
                for (string cLeaf : leafSet){
                    tempSplit.addLeaf(1, cLeaf);
                }
                runningSplits2.push_back(tempSplit);
                i++;
                break;
            case ')':
                tempSplit = runningSplits2.back();
                runningSplits1.push_back(tempSplit);
                runningSplits2.pop_back();
                i++;
                break;
            case ',':
                i++;
                break;
            default:
                int pos1 = static_cast<int>(workingNewick.find(',',i));
                int pos2 = static_cast<int>(workingNewick.find(')',i));
                int nxtPos = static_cast<int>(workingNewick.length());
                if ((pos1 > -1) && (pos2 > -1)){
                    if (pos1 < pos2){
                        nxtPos = pos1;
                    } else {
                        nxtPos = pos2;
                    }
                } else if (pos1 > -1){
                    nxtPos = pos1;
                } else if (pos2 > -1){
                    nxtPos = pos2;
                }
                string nwLeaf = workingNewick.substr(i,nxtPos - i);
                leafSet.insert(nwLeaf);
                for (int k = 0; k < runningSplits1.size(); k++){
                    tempSplit = runningSplits1[k];
                    tempSplit.addLeaf(1, nwLeaf);
                    runningSplits1.at(k) = tempSplit;
                }
                for (int k = 0; k < runningSplits2.size(); k++){
                    tempSplit = runningSplits2[k];
                    tempSplit.addLeaf(2, nwLeaf);
                    runningSplits2.at(k) = tempSplit;
                }
                i = nxtPos;
                break;
        }
    }

    for (Split cSplit : runningSplits1){
        Split newSplit = Split(cSplit.side1, cSplit.side2);
        intSplits.insert(newSplit);
    }

    complexity = -1;
    rank = max((static_cast<int>(leafSet.size()) + (static_cast<int>(intSplits.size()) - 4)),0);

    if (static_cast<int>(intSplits.size()) == 0) {
        rank = 0;
    }
}

pTree::pTree(set<string> nwLeafSet, set<Split> nwIntSplits){
    if (nwIntSplits.size()>0){
        leafSet = nwLeafSet;
        intSplits = nwIntSplits;

        complexity = -1;
        rank = max((static_cast<int>(leafSet.size()) + (static_cast<int>(intSplits.size()) - 4)),0);
        if (static_cast<int>(intSplits.size()) == 0) {
            rank = 0;
        }
    } else {
        leafSet = set<string>();
        intSplits = set<Split>();
        complexity = 0;
        rank = 0;
    }
}

pTree::pTree(set<string> nwLeafSet, set<Split> nwIntSplits, int nwComplexity){
    leafSet = nwLeafSet;
    intSplits = nwIntSplits;

    complexity = nwComplexity;
    rank = max((static_cast<int>(leafSet.size()) + (static_cast<int>(intSplits.size()) - 4)),0);
}

void pTree::setComplexity(float w){
    complexity = (static_cast<int>(leafSet.size()) + w*(static_cast<int>(intSplits.size()) - 4));
}

bool pTree::operator==(const pTree& other) const{
    return (intSplits == other.intSplits);
}

void pTree::print(){
    cout << "Leaf set: \n";
    bool first = true;
    for (string sLeaf : leafSet){
        if(!first) cout << ", ";
        cout << sLeaf;
        first = false;
    }
    cout << "\n \n";

    cout << "The set of splits is: \n";

    for (Split sp1 : intSplits){
        sp1.print();
    }

    cout << "\n";
}

string pTree::printSt(){
    string Resp = "Leaf set: \n";

    bool first = true;
    for (string sLeaf : leafSet){
        if(!first) Resp = Resp + ", ";
        Resp = Resp + sLeaf;
        first = false;
    }
    Resp = Resp + "\n \n";

    Resp = Resp + "The set of splits is: \n";

    for (Split sp1 : intSplits){
        Resp = Resp + sp1.printSt() + "\n";
    }

    Resp = Resp + "\n";

    return Resp;
}

pTree pTree::TDR(set<string> L){
    set<string> leafIntersection;
    set_intersection(leafSet.begin(), leafSet.end(),
                     L.begin(), L.end(),
                     inserter(leafIntersection, leafIntersection.begin()));

    set<Split> newSplits;
    for (Split sp1 : intSplits){
        Split nwSp1 = sp1.TDR(L);
        if (nwSp1.isInternal()){
            newSplits.insert(nwSp1);
        }
    }

    return pTree(leafIntersection, newSplits);
}

pTree pTree::Remove(Split s){
    set<Split> Copy = intSplits;
    Copy.erase(s);
    return (pTree(leafSet,Copy));
}

pTree pTree::Remove(string a){
    set<string> Copy = leafSet;
    Copy.erase(a);
    return (TDR(Copy));
}

pTree pTree::Insert(Split s){
    set<string> commonLeaves;

    const auto& leavesInS = s.LeavesInSplit();
    set_intersection(leafSet.begin(), leafSet.end(),
                     leavesInS.begin(), leavesInS.end(),
                     inserter(commonLeaves, commonLeaves.begin()));

    set<Split> nwSplits;
    for (Split st : intSplits){
        Split nSt = st.TDR(commonLeaves);
        if (nSt.isInternal()){
            nwSplits.insert(nSt);
        }
    }
    Split nS = s.TDR(commonLeaves);
    if (nS.isInternal()){
        nwSplits.insert(nS.TDR(commonLeaves));
    }

    return pTree(commonLeaves, nwSplits);
}

int pTree::returnRank(){
    return (rank);
}

int pTree::returnComplexity(){
    return (complexity);
}

bool pTree::over(pTree tOther){
    if (includes(leafSet.begin(),leafSet.end(),
                 tOther.leafSet.begin(), tOther.leafSet.end())){
        set<Split> MappedSplits;
        for (Split spt : intSplits){
            MappedSplits.insert(spt.TDR(tOther.leafSet));
        }
        if (includes(MappedSplits.begin(),MappedSplits.end(),
                     tOther.intSplits.begin(),tOther.intSplits.end())){
            return true;
        }
    }
    return false;
}

bool pTree::covers(pTree tOther){
    if(over(tOther) && (rank == tOther.rank + 1)){
        return true;
    }
    return false;
}
