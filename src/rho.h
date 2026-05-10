#ifndef RHO_H
#define RHO_H

#include "pTree.h"

int rho(pTree T1, pTree T2);

class oRho{
    public:
        int rho;
        std::vector<std::set<std::string>> presLeaves;
    
    oRho(int nrho, std::vector<std::set<std::string>> newPresLeaves);
};

oRho rho(pTree U, pTree V, oRho baseORho, pTree Tl, std::string extral);

oRho rho(pTree U, pTree V, oRho baseORho, pTree Tl, Split extraS);

oRho rho(pTree V, pTree U, oRho baseORho, pTree Tl);

#endif