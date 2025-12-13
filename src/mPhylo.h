#ifndef MPHYLO_H
#define MPHYLO_H

#include "pTree.h"
#include <set>
#include <queue>
#include <vector>
#include <stack>
#include <string>


class mPhylo{
    public:
        std::vector<std::vector<int>> edge;
        int Nnode;
        std::vector<std::string> tipLabel;
        std::vector<std::vector<int>> node2edge;
        std::unordered_map<std::string, int> label2node;
        
    mPhylo(pTree T);
    
    void print();
    
    std::string toNewick();
    
};

#endif