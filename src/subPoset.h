#ifndef SUBPOSET_H
#define SUBPOSET_H

#include "pTree.h"
#include "mPhylo.h"
#include <set>
#include <queue>
#include <vector>
#include <stack>
#include <string>

class spNode{
    public:
        pTree Tree;
        int kappa;
        float zeta;
    
        int next;
        std::vector<int> over;
        std::vector<int> under;
        std::vector<int64_t> chainCountIe;
        std::vector<int> boundAntichain;
        
    
    spNode(pTree eTree);
    
    void addChild(int newChild);
    
    void addParent(int newParent);
    
    void setNext(int newNext);
    
    void setKappa(int newKappa);
    
    void setZeta (float newZ);
    
    void print();
    
    void printRd();
};


class subPoset{
    public:
        std::vector<spNode> Poset;
    
        std::vector<int> firstRank;
        std::vector<int> lastRank;
    
        int Msize;
    
    subPoset(std::vector<pTree> initT, std::vector<pTree> Sample, std::set<std::string> compLeafSet, int rb);
    
    subPoset(std::vector<pTree> initT, std::vector<pTree> Sample, std::vector<int> nSample, std::set<std::string> compLeafSet, int rb);
    
    subPoset(std::vector<pTree> Sample, std::vector<int> nSample, std::set<std::string> compLeafSet, int Mt, int rb);
    
    subPoset(std::vector<pTree> Sample, std::set<std::string> compLeafSet, int Mt, int rb);
    
    subPoset(std::vector<pTree> Sample, std::vector<int> nSample, std::set<std::string> compLeafSet, int Mt, int rb, bool Constructive);
    
    subPoset(std::vector<pTree> Sample, std::vector<int> nSample, std::set<std::string> compLeafSet, float q);

    subPoset(std::vector<pTree> Sample, std::vector<int> nSample, std::set<std::string> compLeafSet, int top_width, int bottom_width, std::string orientation);

    void print();
    
    void printRd();
};

std::vector<bool> computeDesc(const subPoset& SP, int y);
std::vector<bool> computeAnc(const subPoset& SP, int x);
int64_t countMaximalChains(const subPoset& SP, const std::vector<bool>& desc);
int maxLevelInIe(const subPoset& SP, const std::vector<bool>& desc, const std::vector<bool>& anc, const int& curRank);
void computeAllChainCounts(subPoset& SP);
void computeAllMaxLevelBounds(subPoset& SP);
    
#endif