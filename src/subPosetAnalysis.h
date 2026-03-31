#ifndef SUBPOSETANALYSIS_H
#define SUBPOSETANALYSIS_H

#include "pTree.h"
#include "subPoset.h"
#include <vector>
#include <string>

struct subPosetOutput {
    
    float nullCoveringProb; // Probability of null pairs being correctly identified.
    int minLower; // Member of the null covering pair that is being covered where minimum was achieved.
    int minUpper; // Member of the null covering pair that is covering where minimum was achieved.
    std::vector<std::pair<int,int>> edges;     // Covering pairs in subposet (under, over)
    std::vector<bool>        nullCovering;     // Indicator saying if the covering pair is null with respect to the target tree
    std::vector<float>       coveringMean;    // Sample mean of not-null indicator
    std::vector<float>       coveringVariance; // Sample variance of not-null indicator
};

subPosetOutput SPanalisys(pTree tStar, std::vector<pTree> treeSample, std::vector<pTree> bigTreeSample, subPoset SP);
subPosetOutput SPanalisys(pTree tStar, std::vector<pTree> treeSample, std::vector<int> nSample, std::vector<pTree> bigTreeSample, std::vector<int> nBSample, subPoset SP);

#endif