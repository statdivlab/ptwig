#ifndef SUBPOSET_H
#define SUBPOSET_H

#include "pTree.h"
#include <set>
#include <queue>
#include <vector>
#include <stack>
#include <string>

class spNode{
    public:
        pTree Tree;
        int kappa;
    
        int next;
        std::vector<int> over;
        std::vector<int> under;
        
    
    spNode(pTree eTree);
    
    void addChild(int newChild);
    
    void addParent(int newParent);
    
    void setNext(int newNext);
    
    void setKappa(int newKappa);
    
    void print();
    
    void printRd();
};


class subPoset{//Class describing a complex that is used to compute rho(T1,T2)
    public:
        std::vector<spNode> Poset;
    
        std::vector<int> firstRank;
        std::vector<int> lastRank;
    
        int Msize;
    
    subPoset(std::vector<pTree> initT, std::vector<pTree> Sample, std::set<std::string> compLeafSet, int rb);
    
    subPoset(std::vector<pTree> initT, std::vector<pTree> Sample, std::vector<int> nSample, std::set<std::string> compLeafSet, int rb);
    
    void print();
    
    void printRd();
};
    
#endif