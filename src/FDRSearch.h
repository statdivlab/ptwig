#ifndef FDRSEARCH_H
#define FDRSEARCH_H

#include "pTree.h"
#include "subPoset.h"
#include "rho.h"

pTree FDRSearchGreedy(std::vector<pTree> treeSample, std::vector<int> nSample, std::vector<std::vector<oRho>> storedORho, subPoset SP, float q);

#endif