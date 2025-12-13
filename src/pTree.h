#ifndef PTREE_H
#define PTREE_H

#include <set>
#include <string>
#include <vector>

// =====================
// class Split
// =====================
class Split{//This object represents an edge in the tree, by having two sets of strings, each representing a side of the split
    public:
        std::set<std::string> side1;
        std::set<std::string> side2;
    
    Split();
    
    Split(std::set<std::string> inputSet1, std::set<std::string> inputSet2);
    
    bool operator<(const Split& other) const;
    
    bool operator==(const Split& other) const;
    
    void addLeaf(int side, std::string nLeaf);
    
    const void print();
    
    std::string printSt() const;
    
    Split TDR(std::set<std::string> L);
    
    bool isInternal();
    
    std::set<std::string> LeavesInSplit();
    
    bool contains(Split otherS);
};

// =====================
// class pTree
// =====================
class pTree{
    public:
        std::set<std::string> leafSet;
        std::set<Split> intSplits;
        float complexity;
        int rank;
    
    pTree();

    pTree(std::string newick);

    pTree(std::set<std::string> nwLeafSet, std::set<Split> nwIntSplits);
    
    pTree(std::set<std::string> nwLeafSet, std::set<Split> nwIntSplits, int nwComplexity);
    
    void setComplexity(float w);
    
    bool operator==(const pTree& other) const;

    void print();

    std::string printSt();

    pTree TDR(std::set<std::string> L);
    
    pTree Remove(Split s);
    
    pTree Remove(std::string a);
    
    pTree Insert(Split s);

    int returnRank();
    
    int returnComplexity();

    bool over(pTree tOther);
    
    bool covers(pTree tOther);
    
};

#endif