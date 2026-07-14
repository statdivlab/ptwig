#ifndef RHO_H
#define RHO_H

#include "pTree.h"

int rho(const pTree& T1, const pTree& T2);

class oRho{
    public:
        int rho;
        std::vector<std::set<std::string>> presLeaves;
    
    oRho();
    
    oRho(int nrho, std::vector<std::set<std::string>> newPresLeaves);
};
oRho rho(const pTree& U, const pTree& V, const oRho& baseORho, const pTree& Tl, const std::string& extral);

oRho rho(const pTree& U, const pTree& V, const oRho& baseORho, const pTree& Tl, const Split& extraS);

oRho rho(const pTree& V, const pTree& U, const oRho& baseORho, const pTree& Tl);

#endif