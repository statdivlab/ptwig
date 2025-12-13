#ifndef FDRSEARCH_H
#define FDRSEARCH_H

#include "pTree.h"
#include "subPoset.h"

std::vector<pTree> FDRSearch(std::vector<pTree> treeSample, subPoset SP, float q);
std::vector<pTree> FDRSearch(std::vector<pTree> treeSample, std::vector<int> nSample, subPoset SP, float q);

#endif