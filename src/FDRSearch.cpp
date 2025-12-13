#include "pTree.h"
#include "mPhylo.h"
#include "subPoset.h"
#include "rho.h"
#include "coverTrees.h"
#include <iostream>
#include <set>
#include <queue>
#include <vector>
#include <stack>
#include <fstream>
#include <iostream>
#include <string>
#include <chrono>
#include <random>
#include <variant>
#include <utility>
#include <filesystem> // C++17
#include <future>
#include <mutex>
#include <cmath>
#include <Rcpp.h>

using namespace Rcpp;
using namespace std;

float alphaT(int M,int kappa, float q, int rmax, int r){
    float internal = log(static_cast<float>(kappa)) - log(q) + log(static_cast<float>(rmax - r + 1)) - log(rmax);
    float firstTerm = 0;
    if (internal > 0){
        firstTerm = sqrt(internal/M);
    }
    return (firstTerm + 0.5);
}

vector<pTree> FDRSearch(vector<pTree> treeSample, subPoset SP, float q){
    
    vector<int> FDRtrees;
    
    int B = static_cast<int>(treeSample.size());
    int rmax = static_cast<int>(SP.firstRank.size()-1);
    
    stack<int> toProcess;
    stack<float> minScore;
    vector<int> nodeCheck(static_cast<int>(SP.Poset.size()), 0);
    
    int curIndx = SP.firstRank.at(1);
    vector<float> Scores1;
    vector<int> Indexes1;
    while (curIndx > -1) {
        Indexes1.push_back(curIndx);
        float sum = 0;
        for (pTree T : treeSample){
            if (rho(SP.Poset.at(curIndx).Tree, T) > 0){
                sum++;
            }
        }
        Scores1.push_back(sum/B);
        curIndx = SP.Poset.at(curIndx).next;
    }
    
    
    std::vector<int> order(Indexes1.size());
    std::iota(order.begin(), order.end(), 0);

    std::sort(order.begin(), order.end(),
        [&](int a, int b) { return Scores1[a] < Scores1[b]; });
    
    for (int p : order){
        toProcess.push(Indexes1.at(p));
        minScore.push(Scores1.at(p));
    }
    
    int bestRank = 0;
    
    while(!toProcess.empty()){
        curIndx = toProcess.top();
        
        if (nodeCheck.at(curIndx)<0){
            toProcess.pop();
            minScore.pop();
        } else if (nodeCheck.at(curIndx) > 0){
            toProcess.pop();
            
            if ((!toProcess.empty()) && (SP.Poset.at(toProcess.top()).Tree.returnRank() > SP.Poset.at(curIndx).Tree.returnRank())){
                minScore.pop();
            } else {
                
                toProcess.push(curIndx);
                
                vector<float> Scores;
                vector<int> Indexes;
                for(int upIndx : SP.Poset.at(curIndx).over){
                    if (nodeCheck.at(upIndx) == 0){
                        Indexes.push_back(upIndx);
                        float sum = 0;
                        for (pTree T : treeSample){
                            if ((rho(SP.Poset.at(upIndx).Tree, T) - rho(SP.Poset.at(curIndx).Tree, T)) > 0){
                                sum++;
                            }
                        }
                        Scores.push_back(sum/B);
                    }
                }
                
                if (Indexes.empty()){
                    toProcess.pop();
                    minScore.pop();
                    continue;
                }
                
                std::vector<int> order2(Indexes.size());
                std::iota(order2.begin(), order2.end(), 0);

                std::sort(order2.begin(), order2.end(),
                    [&](int a, int b) { return Scores[a] < Scores[b]; });
                
                for (int p : order2){
                    toProcess.push(Indexes.at(p));
                    minScore.push(Scores.at(p));
                }
            }
            
        } else {
            bool checkBelow = false;
            for (int downIndx : SP.Poset.at(curIndx).under){
                if (nodeCheck.at(downIndx) == 0){
                    toProcess.push(downIndx);
                    minScore.push(1.0);
                    checkBelow = true;
                    break;
                }
            }
            if (checkBelow){
                continue;
            }
            
            float curScore = minScore.top();
            
            if (SP.Poset.at(curIndx).under.empty()){
                float sum = 0;
                for (pTree T : treeSample){
                    if (rho(SP.Poset.at(curIndx).Tree, T) > 0){
                        sum++;
                    }
                }
                curScore = sum/B;
            } else {
                for(int downIndx : SP.Poset.at(curIndx).under){
                    float sum = 0;
                    for (pTree T : treeSample){
                        if ((rho(SP.Poset.at(curIndx).Tree, T) - rho(SP.Poset.at(downIndx).Tree, T)) > 0){
                            sum++;
                        }
                    }
                    if (sum/B < curScore){
                        curScore = sum/B;
                    }
                }
            }
            
            minScore.pop();
            minScore.push(curScore);
            
            float alphat = alphaT(B, SP.Poset.at(curIndx).kappa, q, rmax, SP.Poset.at(curIndx).Tree.returnRank());
            
            if (curScore >= alphat){
                if (bestRank < SP.Poset.at(curIndx).Tree.returnRank()){
                    FDRtrees.clear();
                    FDRtrees.push_back(curIndx);
                    bestRank = SP.Poset.at(curIndx).Tree.returnRank();
                } else if (bestRank == SP.Poset.at(curIndx).Tree.returnRank()){
                    FDRtrees.push_back(curIndx);
                }
                nodeCheck.at(curIndx) = 1;
            } else {
                nodeCheck.at(curIndx) = -1;
                stack<int> cleanUp;
                for (int p : SP.Poset.at(curIndx).over){
                    cleanUp.push(p);
                }
                while (!cleanUp.empty()){
                    int tIndx = cleanUp.top();
                    nodeCheck.at(tIndx) = -1;
                    cleanUp.pop();
                    for (int p : SP.Poset.at(tIndx).over){
                        cleanUp.push(p);
                    }
                }
            }
        }
    }
    
    vector<pTree> Results;
    for (int p : FDRtrees){
        Results.push_back(SP.Poset.at(p).Tree);
    }
    
    return Results;
}

vector<pTree> FDRSearch(vector<pTree> treeSample, vector<int> nSample, subPoset SP, float q){
    
    vector<int> FDRtrees;
    
    int B = std::accumulate(nSample.begin(), nSample.end(), 0);
    int rmax = static_cast<int>(SP.firstRank.size()-1);
    
    stack<int> toProcess;
    stack<float> minScore;
    vector<int> nodeCheck(static_cast<int>(SP.Poset.size()), 0);
    
    int curIndx = SP.firstRank.at(1);
    vector<float> Scores1;
    vector<int> Indexes1;
    while (curIndx > -1) {
        Indexes1.push_back(curIndx);
        float sum = 0;
        int treeIndx = 0;
        for (pTree T : treeSample){
            if (rho(SP.Poset.at(curIndx).Tree, T) > 0){
                sum += nSample.at(treeIndx);
            }
            treeIndx++;
        }
        Scores1.push_back(sum/B);
        curIndx = SP.Poset.at(curIndx).next;
    }
    
    std::vector<int> order(Indexes1.size());
    std::iota(order.begin(), order.end(), 0);

    std::sort(order.begin(), order.end(),
        [&](int a, int b) { return Scores1[a] < Scores1[b]; });
    
    for (int p : order){
        toProcess.push(Indexes1.at(p));
        minScore.push(Scores1.at(p));
    }
    
    int bestRank = 0;
    
    while(!toProcess.empty()){
        curIndx = toProcess.top();
        
        if (nodeCheck.at(curIndx)<0){
            toProcess.pop();
            minScore.pop();
        } else if (nodeCheck.at(curIndx) > 0){
            toProcess.pop();
            
            if ((!toProcess.empty()) && (SP.Poset.at(toProcess.top()).Tree.returnRank() > SP.Poset.at(curIndx).Tree.returnRank())){
                minScore.pop();
            } else {
                
                toProcess.push(curIndx);
                
                vector<float> Scores;
                vector<int> Indexes;
                for(int upIndx : SP.Poset.at(curIndx).over){
                    if (nodeCheck.at(upIndx) == 0){
                        Indexes.push_back(upIndx);
                        float sum = 0;
                        int treeIndx = 0;
                        for (pTree T : treeSample){
                            if ((rho(SP.Poset.at(upIndx).Tree, T) - rho(SP.Poset.at(curIndx).Tree, T)) > 0){
                                sum += nSample.at(treeIndx);
                            }
                            treeIndx++;
                        }
                        Scores.push_back(sum/B);
                    }
                }
                
                if (Indexes.empty()){
                    toProcess.pop();
                    minScore.pop();
                    continue;
                }
                
                std::vector<int> order2(Indexes.size());
                std::iota(order2.begin(), order2.end(), 0);

                std::sort(order2.begin(), order2.end(),
                    [&](int a, int b) { return Scores[a] < Scores[b]; });
                
                for (int p : order2){
                    toProcess.push(Indexes.at(p));
                    minScore.push(Scores.at(p));
                }
            }
            
        } else {
            bool checkBelow = false;
            for (int downIndx : SP.Poset.at(curIndx).under){
                if (nodeCheck.at(downIndx) == 0){
                    toProcess.push(downIndx);
                    minScore.push(1.0);
                    checkBelow = true;
                    break;
                }
            }
            if (checkBelow){
                continue;
            }
            
            float curScore = minScore.top();
            
            if (SP.Poset.at(curIndx).under.empty()){
                float sum = 0;
                int treeIndx = 0;
                for (pTree T : treeSample){
                    if (rho(SP.Poset.at(curIndx).Tree, T) > 0){
                        sum += nSample.at(treeIndx);
                    }
                    treeIndx++;
                }
                curScore = sum/B;
            } else {
                for(int downIndx : SP.Poset.at(curIndx).under){
                    float sum = 0;
                    int treeIndx = 0;
                    for (pTree T : treeSample){
                        if ((rho(SP.Poset.at(curIndx).Tree, T) - rho(SP.Poset.at(downIndx).Tree, T)) > 0){
                            sum += nSample.at(treeIndx);
                        }
                        treeIndx++;
                    }
                    if (sum/B < curScore){
                        curScore = sum/B;
                    }
                }
            }
            
            minScore.pop();
            minScore.push(curScore);
            
            float alphat = alphaT(B, SP.Poset.at(curIndx).kappa, q, rmax, SP.Poset.at(curIndx).Tree.returnRank());
            
            if (curScore >= alphat){
                if (bestRank < SP.Poset.at(curIndx).Tree.returnRank()){
                    FDRtrees.clear();
                    FDRtrees.push_back(curIndx);
                    bestRank = SP.Poset.at(curIndx).Tree.returnRank();
                } else if (bestRank == SP.Poset.at(curIndx).Tree.returnRank()){
                    FDRtrees.push_back(curIndx);
                }
                nodeCheck.at(curIndx) = 1;
            } else {
                nodeCheck.at(curIndx) = -1;
                stack<int> cleanUp;
                for (int p : SP.Poset.at(curIndx).over){
                    cleanUp.push(p);
                }
                while (!cleanUp.empty()){
                    int tIndx = cleanUp.top();
                    nodeCheck.at(tIndx) = -1;
                    cleanUp.pop();
                    for (int p : SP.Poset.at(tIndx).over){
                        cleanUp.push(p);
                    }
                }
            }
        }
    }
    
    vector<pTree> Results;
    for (int p : FDRtrees){
        Results.push_back(SP.Poset.at(p).Tree);
    }
    
    return Results;
}