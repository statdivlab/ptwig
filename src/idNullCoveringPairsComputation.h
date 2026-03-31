#ifndef IDNULLCOVERINGPAIRSCOMPUTATION_H
#define IDNULLCOVERINGPAIRSCOMPUTATION_H

#include "pTree.h"
#include "subPoset.h"

float idNullCoveringPairsComputation(pTree tStar, std::vector<pTree> treeSample, subPoset SP);
float idNullCoveringPairsComputation(pTree tStar, std::vector<pTree> treeSample, std::vector<int> nSample, subPoset SP);

#endif