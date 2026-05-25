#ifndef FDRSEARCH_H
#define FDRSEARCH_H

#include "pTree.h"
#include "subPoset.h"

std::vector<pTree> FDRSearch(std::vector<pTree> treeSample, subPoset SP, float q);
std::vector<pTree> FDRSearch(std::vector<pTree> treeSample, std::vector<int> nSample, subPoset SP, float q);
std::vector<pTree> FDRSearch(std::vector<pTree> treeSample, subPoset SP, std::vector<float> lbEta, float q);
std::vector<pTree> FDRSearch(std::vector<pTree> treeSample, std::vector<int> nSample, subPoset SP, std::vector<float> lbEta, float q);

pTree FDRSearchGreedy(std::vector<pTree> treeSample, std::vector<int> nSample, subPoset SP, std::vector<float> lbEta, float q);

pTree FDRSearchGreedy(std::vector<pTree> treeSample, std::vector<int> nSample, std::vector<std::vector<oRho>> storedORho, subPoset SP, std::vector<float> lbEta, float q);
    
pTree FDRSearchGreedy(std::vector<pTree> treeSample, std::vector<std::vector<oRho>> storedORho, subPoset SP, std::vector<float> lbEta, float q);

pTree FDRSearchGreedy(std::vector<pTree> treeSample, std::vector<int> nSample, std::vector<std::vector<oRho>> storedORho, subPoset SP, float q);

#endif