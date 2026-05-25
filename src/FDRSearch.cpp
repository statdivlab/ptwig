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

float kapThreshold(float omegaT, float q, int N, int kappa){
    float prelim = log(kappa) + log(omegaT) - log(q);
    //cout << "Prelim is " << prelim << "\n"<< std::flush; 
    float thrs1 = std::max(prelim/(2*N), 0.0f);
    
    //cout << "Thrs1 " << thrs1 << "\n"<< std::flush; 
    
    return(sqrt(thrs1));
    
}

float nuThreshold(float omegaT, float q, int N, int nu, float zeta){
    float prelim = log(nu) + log(omegaT) - log(q);
    //cout << "Prelim is " << prelim << "\n"<< std::flush; 
    float thrs1 = std::max(prelim/(2*N), 0.0f);
    
    //cout << "Thrs1 " << thrs1 << "\n"<< std::flush; 
    
    return(sqrt(thrs1) + zeta);
    
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
                cout<< "Colored Green At: " << curIndx << "\n"<< std::flush;
            } else {
                nodeCheck.at(curIndx) = -1;
                cout<< "Colored Red At: " << curIndx << "\n"<< std::flush;
                stack<int> cleanUp;
                for (int p : SP.Poset.at(curIndx).over){
                    cleanUp.push(p);
                }
                while (!cleanUp.empty()){
                    int tIndx = cleanUp.top();
                    nodeCheck.at(tIndx) = -1;
                    cout<< "Colored Red At: " << tIndx << "\n"<< std::flush;
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
                cout<< "Colored Green At: " << curIndx << "\n"<< std::flush;
            } else {
                nodeCheck.at(curIndx) = -1;
                cout<< "Colored Red At: " << curIndx << "\n"<< std::flush;
                stack<int> cleanUp;
                for (int p : SP.Poset.at(curIndx).over){
                    cleanUp.push(p);
                }
                while (!cleanUp.empty()){
                    int tIndx = cleanUp.top();
                    nodeCheck.at(tIndx) = -1;
                    cout<< "Colored Red At: " << tIndx << "\n"<< std::flush;
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

vector<pTree> FDRSearch(vector<pTree> treeSample, subPoset SP, vector<float> lbEta, float q){
    
    vector<int> FDRtrees;
    
    int B = static_cast<int>(treeSample.size());
    int rmax = static_cast<int>(SP.firstRank.size()-1);
    
    stack<int> toProcess;
    stack<float> minScore;
    vector<int> nodeCheck(static_cast<int>(SP.Poset.size()), 0);
    
    std::vector<std::vector<int>> allRho(static_cast<int>(SP.Poset.size()), std::vector<int>(static_cast<int>(treeSample.size()), -1));
    
    int curIndx = SP.firstRank.at(1);
    vector<float> Scores1;
    vector<int> Indexes1;
    cout << "Computing rhos and scores for the first level \n"<< std::flush;
    while (curIndx > -1) {
        Indexes1.push_back(curIndx);
        float sum = 0;
        cout << "Processing for index" << curIndx << "\n"<< std::flush;
        for (int i = 0; i < treeSample.size(); ++i){
            pTree T = treeSample.at(i);
            allRho[curIndx][i] = rho(SP.Poset.at(curIndx).Tree, T);
            if (allRho[curIndx][i] > 0){
                sum++;
            }
        }
        Scores1.push_back((sum/B) + lbEta[lbEta.size() - 1] - 1);
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
                        for (int treeIndx = 0; treeIndx < treeSample.size(); treeIndx++){
                            pTree T = treeSample.at(treeIndx);
                            if (allRho[curIndx][treeIndx] < 0 ){
                                allRho[curIndx][treeIndx] = rho(SP.Poset.at(curIndx).Tree, T);
                            }
                            if (allRho[upIndx][treeIndx] < 0 ){
                                allRho[upIndx][treeIndx] = rho(SP.Poset.at(upIndx).Tree, T);
                            }
                            if ((allRho[upIndx][treeIndx] - allRho[curIndx][treeIndx]) > 0){
                                sum++;
                            }
                        }
                        Scores.push_back((sum/B) + lbEta[curIndx] - 1);
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
                    minScore.push(2.0);
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
                if (allRho[curIndx][0] < 0){
                    cout << "Computing rhos for the index " << curIndx << "\n"<< std::flush;
                    for (int treeIndx = 0; treeIndx < treeSample.size(); treeIndx++){
                        pTree T = treeSample.at(treeIndx);
                        allRho[curIndx][treeIndx] = rho(SP.Poset.at(curIndx).Tree, T);
                    }
                }
                for (int treeIndx = 0; treeIndx < treeSample.size(); treeIndx++){
                    if (allRho[curIndx][treeIndx] > 0){
                        sum++;
                    }
                }
                curScore = (sum/B) + lbEta[lbEta.size() - 1] - 1;
            } else {
                if (allRho[curIndx][0] < 0){
                    cout << "Computing rhos for the index " << curIndx << "\n"<< std::flush;
                    for (int treeIndx = 0; treeIndx < treeSample.size(); treeIndx++){
                        pTree T = treeSample.at(treeIndx);
                        allRho[curIndx][treeIndx] = rho(SP.Poset.at(curIndx).Tree, T);
                    }
                }
                for(int downIndx : SP.Poset.at(curIndx).under){
                    float sum = 0;
                    if (allRho[downIndx][0] < 0){
                        cout << "Computing rhos for the index " << downIndx << "\n"<< std::flush;
                        for (int treeIndx = 0; treeIndx < treeSample.size(); treeIndx++){
                            pTree T = treeSample.at(treeIndx);
                            allRho[downIndx][treeIndx] = rho(SP.Poset.at(downIndx).Tree, T);
                        }
                    }
                    for (int treeIndx = 0; treeIndx < treeSample.size(); treeIndx++){
                        pTree T = treeSample.at(treeIndx);
                        if ((allRho[curIndx][treeIndx] - allRho[downIndx][treeIndx]) > 0){
                            sum++;
                        }
                    }
                    if (((sum/B) + lbEta[downIndx] - 1) < curScore){
                        curScore = (sum/B) + lbEta[downIndx] - 1;
                    }
                }
            }
            
            minScore.pop();
            minScore.push(curScore);
            
            float omegaTemp =  static_cast<float>(rmax - SP.Poset.at(curIndx).Tree.rank + 1)/(static_cast<float>(rmax));
            
            float premThrs = kapThreshold(omegaTemp, q, B, SP.Poset.at(curIndx).kappa);
            
            if (curScore >= premThrs){
                if (bestRank < SP.Poset.at(curIndx).Tree.returnRank()){
                    FDRtrees.clear();
                    FDRtrees.push_back(curIndx);
                    bestRank = SP.Poset.at(curIndx).Tree.returnRank();
                } else if (bestRank == SP.Poset.at(curIndx).Tree.returnRank()){
                    FDRtrees.push_back(curIndx);
                }
                nodeCheck.at(curIndx) = 1;
                cout<< "Colored Green At: " << curIndx << "\n"<< std::flush;
            } else {
                nodeCheck.at(curIndx) = -1;
                cout<< "Colored Red At: " << curIndx << "\n"<< std::flush;
                stack<int> cleanUp;
                for (int p : SP.Poset.at(curIndx).over){
                    cleanUp.push(p);
                }
                while (!cleanUp.empty()){
                    int tIndx = cleanUp.top();
                    nodeCheck.at(tIndx) = -1;
                    cout<< "Colored Red At: " << tIndx << "\n"<< std::flush;
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

vector<pTree> FDRSearch(vector<pTree> treeSample, vector<int> nSample, subPoset SP, vector<float> lbEta, float q){
    
    vector<int> FDRtrees;
    
    int B = std::accumulate(nSample.begin(), nSample.end(), 0);
    int rmax = static_cast<int>(SP.firstRank.size()-1);
    
    stack<int> toProcess;
    stack<float> minScore;
    vector<int> nodeCheck(static_cast<int>(SP.Poset.size()), 0);
    
    std::vector<std::vector<int>> allRho(static_cast<int>(SP.Poset.size()), std::vector<int>(static_cast<int>(treeSample.size()), -1));
    
    int curIndx = SP.firstRank.at(1);
    vector<float> Scores1;
    vector<int> Indexes1;
    cout << "Computing rhos and scores for the first level \n"<< std::flush;
    while (curIndx > -1) {
        Indexes1.push_back(curIndx);
        float sum = 0;
        cout << "Processing for index" << curIndx << "\n"<< std::flush;
        for (int i = 0; i < treeSample.size(); ++i){
            pTree T = treeSample.at(i);
            allRho[curIndx][i] = rho(SP.Poset.at(curIndx).Tree, T);
            if (allRho[curIndx][i] > 0){
                sum += nSample.at(i);
            }
        }
        Scores1.push_back((sum/B) + lbEta[lbEta.size() - 1] - 1);
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
                        for (int treeIndx = 0; treeIndx < treeSample.size(); treeIndx++){
                            pTree T = treeSample.at(treeIndx);
                            if (allRho[curIndx][treeIndx] < 0 ){
                                allRho[curIndx][treeIndx] = rho(SP.Poset.at(curIndx).Tree, T);
                            }
                            if (allRho[upIndx][treeIndx] < 0 ){
                                allRho[upIndx][treeIndx] = rho(SP.Poset.at(upIndx).Tree, T);
                            }
                            if ((allRho[upIndx][treeIndx] - allRho[curIndx][treeIndx]) > 0){
                                sum += nSample.at(treeIndx);
                            }
                        }
                        Scores.push_back((sum/B) + lbEta[curIndx] - 1);
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
                    minScore.push(2.0);
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
                if (allRho[curIndx][0] < 0){
                    cout << "Computing rhos for the index " << curIndx << "\n"<< std::flush;
                    for (int treeIndx = 0; treeIndx < treeSample.size(); treeIndx++){
                        pTree T = treeSample.at(treeIndx);
                        allRho[curIndx][treeIndx] = rho(SP.Poset.at(curIndx).Tree, T);
                    }
                }
                for (int treeIndx = 0; treeIndx < treeSample.size(); treeIndx++){
                    if (allRho[curIndx][treeIndx] > 0){
                        sum += nSample.at(treeIndx);
                    }
                }
                curScore = (sum/B) + lbEta[lbEta.size() - 1] - 1;
            } else {
                if (allRho[curIndx][0] < 0){
                    cout << "Computing rhos for the index " << curIndx << "\n"<< std::flush;
                    for (int treeIndx = 0; treeIndx < treeSample.size(); treeIndx++){
                        pTree T = treeSample.at(treeIndx);
                        allRho[curIndx][treeIndx] = rho(SP.Poset.at(curIndx).Tree, T);
                    }
                }
                for(int downIndx : SP.Poset.at(curIndx).under){
                    float sum = 0;
                    if (allRho[downIndx][0] < 0){
                        cout << "Computing rhos for the index " << downIndx << "\n"<< std::flush;
                        for (int treeIndx = 0; treeIndx < treeSample.size(); treeIndx++){
                            pTree T = treeSample.at(treeIndx);
                            allRho[downIndx][treeIndx] = rho(SP.Poset.at(downIndx).Tree, T);
                        }
                    }
                    for (int treeIndx = 0; treeIndx < treeSample.size(); treeIndx++){
                        pTree T = treeSample.at(treeIndx);
                        if ((allRho[curIndx][treeIndx] - allRho[downIndx][treeIndx]) > 0){
                            sum += nSample.at(treeIndx);
                        }
                    }
                    if (((sum/B) + lbEta[downIndx] - 1) < curScore){
                        curScore = (sum/B) + lbEta[downIndx] - 1;
                    }
                }
            }
            
            minScore.pop();
            minScore.push(curScore);
            
            float omegaTemp =  static_cast<float>(rmax - SP.Poset.at(curIndx).Tree.rank + 1)/(static_cast<float>(rmax));
            
            float premThrs = kapThreshold(omegaTemp, q, B, SP.Poset.at(curIndx).kappa);
            
            if (curScore >= premThrs){
                if (bestRank < SP.Poset.at(curIndx).Tree.returnRank()){
                    FDRtrees.clear();
                    FDRtrees.push_back(curIndx);
                    bestRank = SP.Poset.at(curIndx).Tree.returnRank();
                } else if (bestRank == SP.Poset.at(curIndx).Tree.returnRank()){
                    FDRtrees.push_back(curIndx);
                }
                nodeCheck.at(curIndx) = 1;
                cout<< "Colored Green At: " << curIndx << "\n"<< std::flush;
            } else {
                nodeCheck.at(curIndx) = -1;
                cout<< "Colored Red At: " << curIndx << "\n"<< std::flush;
                stack<int> cleanUp;
                for (int p : SP.Poset.at(curIndx).over){
                    cleanUp.push(p);
                }
                while (!cleanUp.empty()){
                    int tIndx = cleanUp.top();
                    nodeCheck.at(tIndx) = -1;
                    cout<< "Colored Red At: " << tIndx << "\n"<< std::flush;
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

pTree FDRSearchGreedy(vector<pTree> treeSample, vector<int> nSample, subPoset SP, vector<float> lbEta, float q) {

    int B = std::accumulate(nSample.begin(), nSample.end(), 0);
    int rmax = static_cast<int>(SP.firstRank.size() - 1);
    int numTrees = static_cast<int>(treeSample.size());
    int numNodes = static_cast<int>(SP.Poset.size());

    // Lazily computed rho cache: -1 means not yet computed
    std::vector<std::vector<int>> allRho(numNodes, std::vector<int>(numTrees, -1));

    // Helper: ensure all rho values for a given node are computed
    auto ensureRho = [&](int nodeIdx) {
        if (allRho[nodeIdx][0] < 0) {
            for (int t = 0; t < numTrees; ++t)
                allRho[nodeIdx][t] = rho(SP.Poset.at(nodeIdx).Tree, treeSample.at(t));
        }
    };

    // Score for the bottom-level transition (rank-1 node standing alone):
    // fraction of trees for which rho(node, T) > 0, shifted by lbEta
    auto scoreBase = [&](int nodeIdx) -> float {
        ensureRho(nodeIdx);
        float sum = 0;
        for (int t = 0; t < numTrees; ++t)
            if (allRho[nodeIdx][t] > 0)
                sum += nSample.at(t);
        return (sum / B) + lbEta[lbEta.size() - 1] - 1;
    };

    // Score for the transition from parentIdx -> childIdx (moving up):
    // fraction of trees where rho increases, shifted by lbEta of the parent
    auto scoreTransition = [&](int parentIdx, int childIdx) -> float {
        ensureRho(parentIdx);
        ensureRho(childIdx);
        float sum = 0;
        for (int t = 0; t < numTrees; ++t)
            if ((allRho[childIdx][t] - allRho[parentIdx][t]) > 0)
                sum += nSample.at(t);
        return (sum / B) + lbEta[parentIdx] - 1;
    };

    // Threshold for a given node
    auto threshold = [&](int nodeIdx) -> float {
        float omega = static_cast<float>(rmax - SP.Poset.at(nodeIdx).Tree.rank + 1)
                    / static_cast<float>(rmax);
        return kapThreshold(omega, q, B, SP.Poset.at(nodeIdx).kappa);
    };

    // ----------------------------------------------------------------
    // Step 1: scan rank-1 nodes in random order, pick first that passes
    // ----------------------------------------------------------------
    int curIndx = SP.firstRank.at(1);

    // Collect all rank-1 node indices
    vector<int> rank1Nodes;
    while (curIndx > -1) {
        rank1Nodes.push_back(curIndx);
        curIndx = SP.Poset.at(curIndx).next;
    }

    // Shuffle for random order
    auto rd  = std::random_device{};
    auto rng = std::default_random_engine{ rd() };
    shuffle(rank1Nodes.begin(), rank1Nodes.end(), rng);

    int current = -1;
    for (int idx : rank1Nodes) {
        float s = scoreBase(idx);
        cout << "Rank-1 node " << idx << " score=" << s
             << " thresh=" << threshold(idx) << "\n" << std::flush;
        if (s >= threshold(idx)) {
            current = idx;
            cout << "Selected rank-1 node " << idx << "\n" << std::flush;
            break;
        }
    }

    if (current == -1) {
        cout << "No rank-1 node passes threshold. Returning empty.\n" << std::flush;
        return {pTree("();")};
    }

    // ----------------------------------------------------------------
    // Step 2: greedily climb upward
    // ----------------------------------------------------------------
    while (!SP.Poset.at(current).over.empty()) {
        const vector<int>& candidates = SP.Poset.at(current).over;

        // Shuffle candidates for random order
        vector<int> shuffled(candidates.begin(), candidates.end());
        auto rd2  = std::random_device{};
        auto rng2 = std::default_random_engine{ rd2() };
        shuffle(shuffled.begin(), shuffled.end(),rng2);

        int next = -1;
        for (int upIdx : shuffled) {
            float s = scoreTransition(current, upIdx);
            cout << "  Transition " << current << " -> " << upIdx
                 << " score=" << s << " thresh=" << threshold(upIdx) << "\n" << std::flush;
            if (s >= threshold(upIdx)) {
                next = upIdx;
                cout << "  Moving up to " << upIdx << "\n" << std::flush;
                break;
            }
        }

        if (next == -1) {
            cout << "No upward transition passes from node " << current
                 << ". Stopping.\n" << std::flush;
            break;  // stuck — report current as result
        }

        current = next;
    }

    cout << "Final node: " << current << "\n" << std::flush;
    return { SP.Poset.at(current).Tree };
}

pTree FDRSearchGreedy(vector<pTree> treeSample, vector<int> nSample, vector<vector<oRho>> storedORho, subPoset SP, vector<float> lbEta, float q) {

    int B = std::accumulate(nSample.begin(), nSample.end(), 0);
    int rmax = static_cast<int>(SP.firstRank.size() - 1);
    int numTrees = static_cast<int>(treeSample.size());
    int numNodes = static_cast<int>(SP.Poset.size());

    // Score for the bottom-level transition (rank-1 node standing alone):
    // fraction of trees for which rho(node, T) > 0, shifted by lbEta
    auto scoreBase = [&](int nodeIdx) -> float {
        float sum = 0;
        for (int t = 0; t < numTrees; ++t)
            if (storedORho.at(nodeIdx).at(t).rho > 0)
                sum += nSample.at(t);
        return (sum / B) + lbEta[lbEta.size() - 1] - 1;
    };

    // Score for the transition from parentIdx -> childIdx (moving up):
    // fraction of trees where rho increases, shifted by lbEta of the parent
    auto scoreTransition = [&](int TaIdx, int TbIdx) -> float {
        float sum = 0;
        for (int t = 0; t < numTrees; ++t)
            if ((storedORho.at(TbIdx).at(t).rho - storedORho.at(TaIdx).at(t).rho) > 0)
                sum += nSample.at(t);
        return (sum / B) + lbEta[TaIdx] - 1;
    };

    // Threshold for a given node
    auto threshold = [&](int nodeIdx) -> float {
        float omega = static_cast<float>(rmax - SP.Poset.at(nodeIdx).Tree.rank + 1)
                    / static_cast<float>(rmax);
        return kapThreshold(omega, q, B, SP.Poset.at(nodeIdx).kappa);
    };

    // ----------------------------------------------------------------
    // Step 1: scan rank-1 nodes in random order, pick first that passes
    // ----------------------------------------------------------------
    int curIndx = SP.firstRank.at(1);

    // Collect all rank-1 node indices
    vector<int> rank1Nodes;
    while (curIndx > -1) {
        rank1Nodes.push_back(curIndx);
        curIndx = SP.Poset.at(curIndx).next;
    }

    // Shuffle for random order
    auto rd  = std::random_device{};
    auto rng = std::default_random_engine{ rd() };
    shuffle(rank1Nodes.begin(), rank1Nodes.end(), rng);

    int current = -1;
    for (int idx : rank1Nodes) {
        float s = scoreBase(idx);
        cout << "Rank-1 node " << idx << " score=" << s
             << " thresh=" << threshold(idx) << "\n" << std::flush;
        if (s >= threshold(idx)) {
            current = idx;
            cout << "Selected rank-1 node " << idx << "\n" << std::flush;
            break;
        }
    }

    if (current == -1) {
        cout << "No rank-1 node passes threshold. Returning empty.\n" << std::flush;
        return {pTree("();")};
    }

    // ----------------------------------------------------------------
    // Step 2: greedily climb upward
    // ----------------------------------------------------------------
    while (!SP.Poset.at(current).over.empty()) {
        const vector<int>& candidates = SP.Poset.at(current).over;

        // Shuffle candidates for random order
        vector<int> shuffled(candidates.begin(), candidates.end());
        auto rd2  = std::random_device{};
        auto rng2 = std::default_random_engine{ rd2() };
        shuffle(shuffled.begin(), shuffled.end(),rng2);

        int next = -1;
        for (int upIdx : shuffled) {
            float s = scoreTransition(current, upIdx);
            cout << "  Transition " << current << " -> " << upIdx
                 << " score=" << s << " thresh=" << threshold(upIdx) << "\n" << std::flush;
            if (s >= threshold(upIdx)) {
                next = upIdx;
                cout << "  Moving up to " << upIdx << "\n" << std::flush;
                break;
            }
        }

        if (next == -1) {
            cout << "No upward transition passes from node " << current
                 << ". Stopping.\n" << std::flush;
            break;  // stuck — report current as result
        }

        current = next;
    }

    cout << "Final node: " << current << "\n" << std::flush;
    return { SP.Poset.at(current).Tree };
}

pTree FDRSearchGreedy(vector<pTree> treeSample, vector<vector<oRho>> storedORho, subPoset SP, vector<float> lbEta, float q) {

    int rmax = static_cast<int>(SP.firstRank.size() - 1);
    int numTrees = static_cast<int>(treeSample.size());
    int numNodes = static_cast<int>(SP.Poset.size());

    // Score for the bottom-level transition (rank-1 node standing alone):
    // fraction of trees for which rho(node, T) > 0, shifted by lbEta
    auto scoreBase = [&](int nodeIdx) -> float {
        float sum = 0;
        for (int t = 0; t < numTrees; ++t)
            if (storedORho.at(nodeIdx).at(t).rho > 0)
                sum++;
        return (sum / numTrees) + lbEta[lbEta.size() - 1] - 1;
    };

    // Score for the transition from parentIdx -> childIdx (moving up):
    // fraction of trees where rho increases, shifted by lbEta of the parent
    auto scoreTransition = [&](int TaIdx, int TbIdx) -> float {
        float sum = 0;
        for (int t = 0; t < numTrees; ++t)
            if ((storedORho.at(TbIdx).at(t).rho - storedORho.at(TaIdx).at(t).rho) > 0)
                sum++;
        return (sum / numTrees) + lbEta[TaIdx] - 1;
    };

    // Threshold for a given node
    auto threshold = [&](int nodeIdx) -> float {
        float omega = static_cast<float>(rmax - SP.Poset.at(nodeIdx).Tree.rank + 1)
                    / static_cast<float>(rmax);
        return kapThreshold(omega, q, numTrees, SP.Poset.at(nodeIdx).kappa);
    };

    // ----------------------------------------------------------------
    // Step 1: scan rank-1 nodes in random order, pick first that passes
    // ----------------------------------------------------------------
    int curIndx = SP.firstRank.at(1);

    // Collect all rank-1 node indices
    vector<int> rank1Nodes;
    while (curIndx > -1) {
        rank1Nodes.push_back(curIndx);
        curIndx = SP.Poset.at(curIndx).next;
    }

    // Shuffle for random order
    auto rd  = std::random_device{};
    auto rng = std::default_random_engine{ rd() };
    shuffle(rank1Nodes.begin(), rank1Nodes.end(), rng);

    int current = -1;
    for (int idx : rank1Nodes) {
        float s = scoreBase(idx);
        cout << "Rank-1 node " << idx << " score=" << s
             << " thresh=" << threshold(idx) << "\n" << std::flush;
        if (s >= threshold(idx)) {
            current = idx;
            cout << "Selected rank-1 node " << idx << "\n" << std::flush;
            break;
        }
    }

    if (current == -1) {
        cout << "No rank-1 node passes threshold. Returning empty.\n" << std::flush;
        return {pTree("();")};
    }

    // ----------------------------------------------------------------
    // Step 2: greedily climb upward
    // ----------------------------------------------------------------
    while (!SP.Poset.at(current).over.empty()) {
        const vector<int>& candidates = SP.Poset.at(current).over;

        // Shuffle candidates for random order
        vector<int> shuffled(candidates.begin(), candidates.end());
        auto rd2  = std::random_device{};
        auto rng2 = std::default_random_engine{ rd2() };
        shuffle(shuffled.begin(), shuffled.end(),rng2);

        int next = -1;
        for (int upIdx : shuffled) {
            float s = scoreTransition(current, upIdx);
            cout << "  Transition " << current << " -> " << upIdx
                 << " score=" << s << " thresh=" << threshold(upIdx) << "\n" << std::flush;
            if (s >= threshold(upIdx)) {
                next = upIdx;
                cout << "  Moving up to " << upIdx << "\n" << std::flush;
                break;
            }
        }

        if (next == -1) {
            cout << "No upward transition passes from node " << current
                 << ". Stopping.\n" << std::flush;
            break;  // stuck — report current as result
        }

        current = next;
    }

    cout << "Final node: " << current << "\n" << std::flush;
    return { SP.Poset.at(current).Tree };
}

pTree FDRSearchGreedy(vector<pTree> treeSample, vector<int> nSample, vector<vector<oRho>> storedORho, subPoset SP, float q){

    int B = std::accumulate(nSample.begin(), nSample.end(), 0);
    int rmax = static_cast<int>(SP.firstRank.size() - 1);
    int numTrees = static_cast<int>(treeSample.size());
    int numNodes = static_cast<int>(SP.Poset.size());

    // Score for the bottom-level transition (rank-1 node standing alone):
    // fraction of trees for which rho(node, T) > 0, shifted by lbEta
    auto scoreBase = [&](int nodeIdx) -> float {
        float sum = 0;
        for (int t = 0; t < numTrees; ++t)
            if (storedORho.at(nodeIdx).at(t).rho > 0)
                sum += nSample.at(t);
        return (sum / B);
    };

    // Score for the transition from parentIdx -> childIdx (moving up):
    // fraction of trees where rho increases, shifted by lbEta of the parent
    auto scoreTransition = [&](int TaIdx, int TbIdx) -> float {
        float sum = 0;
        for (int t = 0; t < numTrees; ++t)
            if ((storedORho.at(TbIdx).at(t).rho - storedORho.at(TaIdx).at(t).rho) > 0)
                sum += nSample.at(t);
        return (sum / B);
    };

    // Threshold for a given node
    auto threshold = [&](int nodeIdx, int childIdx) -> float {
        float omega = static_cast<float>(rmax - SP.Poset.at(nodeIdx).Tree.rank + 1)
                    / static_cast<float>(rmax);
        float taZeta = 1.0f/(3.0f);
        if (SP.Poset.at(nodeIdx).Tree.rank > 1){
            taZeta = min(SP.Poset.at(SP.Poset.at(nodeIdx).under[childIdx]).zeta, 0.5f);
        } 
        
        return kapThreshold(omega, q, B, SP.Poset.at(nodeIdx).boundAntichain[childIdx]);
    };

    // ----------------------------------------------------------------
    // Step 1: scan rank-1 nodes in random order, pick first that passes
    // ----------------------------------------------------------------
    int curIndx = SP.firstRank.at(1);

    // Collect all rank-1 node indices
    vector<int> rank1Nodes;
    while (curIndx > -1) {
        rank1Nodes.push_back(curIndx);
        curIndx = SP.Poset.at(curIndx).next;
    }

    // Shuffle for random order
    auto rd  = std::random_device{};
    auto rng = std::default_random_engine{ rd() };
    shuffle(rank1Nodes.begin(), rank1Nodes.end(), rng);

    int current = -1;
    for (int idx : rank1Nodes) {
        float s = scoreBase(idx);
        cout << "Rank-1 node " << idx << " score=" << s
             << " thresh=" << threshold(idx,0) << "\n" << std::flush;
        if (s >= threshold(idx, 0)) {
            current = idx;
            cout << "Selected rank-1 node " << idx << "\n" << std::flush;
            break;
        }
    }

    if (current == -1) {
        cout << "No rank-1 node passes threshold. Returning empty.\n" << std::flush;
        return {pTree("();")};
    }

    // ----------------------------------------------------------------
    // Step 2: greedily climb upward
    // ----------------------------------------------------------------
    while (!SP.Poset.at(current).over.empty()) {
        const vector<int>& candidates = SP.Poset.at(current).over;

        // Shuffle candidates for random order
        vector<int> shuffled(candidates.begin(), candidates.end());
        auto rd2  = std::random_device{};
        auto rng2 = std::default_random_engine{ rd2() };
        shuffle(shuffled.begin(), shuffled.end(),rng2);

        int next = -1;
        for (int upIdx : shuffled) {
            float s = scoreTransition(current, upIdx);
            int chldIdx = static_cast<int>(find(SP.Poset.at(upIdx).under.begin(), SP.Poset.at(upIdx).under.end(), current) - SP.Poset.at(upIdx).under.begin());
            cout << "  Transition " << current << " -> " << upIdx
                 << " score=" << s << " thresh=" << threshold(upIdx, chldIdx) << "\n" << std::flush;
            if (s >= threshold(upIdx, chldIdx)) {
                next = upIdx;
                cout << "  Moving up to " << upIdx << "\n" << std::flush;
                break;
            }
        }

        if (next == -1) {
            cout << "No upward transition passes from node " << current
                 << ". Stopping.\n" << std::flush;
            break;  // stuck — report current as result
        }

        current = next;
    }

    cout << "Final node: " << current << "\n" << std::flush;
    return { SP.Poset.at(current).Tree };
}