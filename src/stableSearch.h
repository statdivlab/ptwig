#ifndef STABLESEARCH_H
#define STABLESEARCH_H

#include "pTree.h"

std::vector<pTree> stableSearch(std::vector<pTree> treeSample, std::set<std::string> compLeafSet, float alpha);
std::vector<pTree> stableSearch(std::vector<pTree> treeSample, std::vector<int> nSample, std::set<std::string> compLeafSet, float alpha);

#endif