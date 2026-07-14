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
    
    // Trusted constructor: skips the disjointness check in the normal
    // constructor. Only use when you know side1 and side2 are already
    // disjoint (e.g. results coming out of TDR).
    Split(std::set<std::string> inputSet1, std::set<std::string> inputSet2,
          bool trusted);
    
    bool operator<(const Split& other) const;
    
    bool operator==(const Split& other) const;
    
    void addLeaf(int side, std::string nLeaf);
    
    const void print();
    
    std::string printSt() const;
    
    Split TDR(std::set<std::string> L);
    Split TDR(std::set<std::string> L) const;        // const overload for rho
    
    bool isInternal();
    bool isInternal() const;                         // const overload for rho
    
    std::set<std::string> LeavesInSplit();
    std::set<std::string> LeavesInSplit() const;     // const overload for rho
    
    bool contains(Split otherS);
    bool contains(Split otherS) const;               // const overload for rho
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
    pTree TDR(std::set<std::string> L) const;        // const overload for rho/commonLower

    pTree Remove(Split s);
    
    pTree Remove(std::string a);
    
    pTree Insert(Split s);
    pTree Insert(Split s) const;                     // const overload for rho

    pTree InsertCached(const Split& s, const std::set<std::string>& sLeaves);
    
    int returnRank();
    
    int returnComplexity();

    bool over(pTree tOther);
    bool over(pTree tOther) const;
    
    bool covers(pTree tOther);
    bool covers(pTree tOther) const;
    
};

pTree commonLower(const pTree& T1, const pTree& T2, const std::set<std::string>& sLeaves);

#endif