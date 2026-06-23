#include "pTree.h"
#include "mPhylo.h"
#include "rho.h"
#include "coverTrees.h"
#include <Rcpp.h>
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

using namespace Rcpp;
using namespace std;

// [[Rcpp::export]]
bool VerifyingRho(CharacterVector treeR, CharacterVector compLeafSetR) {
    
    pTree Tl = pTree(as<std::string>(treeR));
    
    std::set<std::string> compLeafSet;
    for (int i = 0; i < compLeafSetR.size(); i++) {
        if (compLeafSetR[i] == NA_STRING)
            stop("compLeafSet cannot contain NA.");
        compLeafSet.insert(as<std::string>(compLeafSetR[i]));
    }
    
    pTree curTree = pTree("();");
    vector<set<string>> emptyLeaves;
    oRho curORho = oRho(0,emptyLeaves);
    int curRho = 0;
    
    auto rd = std::random_device {};
    auto rng = std::default_random_engine { rd() };
    
    int rmax = 2*static_cast<int> (compLeafSet.size()) - 7;
    
    int mCount = 0;
    while(curTree.rank < rmax){
        cout << "Entered in "<< mCount << " cycle " << std::endl;
        mCount++;
        vector<pTree> AllV = coverTrees(curTree, compLeafSet);
        shuffle(begin(AllV), std::end(AllV), rng);
        int nRho = -1;
        oRho nORho = oRho(0,emptyLeaves);
        
        cout << "About to start the counter" << std::endl;
        int counter = -1;
        while (nRho <= curRho){
            counter++;
            auto start1 = std::chrono::high_resolution_clock::now();
            nRho = rho(AllV.at(counter),Tl);
            auto end1 = std::chrono::high_resolution_clock::now();
            
            auto start2 = std::chrono::high_resolution_clock::now();
            nORho = rho(AllV.at(counter), curTree, curORho, Tl);
            auto end2 = std::chrono::high_resolution_clock::now();
            
            if (nRho != nORho.rho){
                cout << "Rho's are different for " << mPhylo(AllV.at(counter)).toNewick() << std::endl;
                cout << "Original rho was" << nRho << std::endl;
                cout << "The new rho was" << nORho.rho << std::endl;
                cout << "The U tree is " << mPhylo(curTree).toNewick() << std::endl;
                cout << "The listed leaves where"<< std::endl;
                for (set<string> setLeaves : curORho.presLeaves){
                    cout << "     ";
                    for (string l : setLeaves) cout<< l << " ";
                    cout << " " << std::endl;
                }
                return false;
            }
            
            if (nRho > curRho){
                double ms1 = std::chrono::duration<double, std::milli>(end1 - start1).count();
                double ms2 = std::chrono::duration<double, std::milli>(end2 - start2).count();
                cout << "Computing the first rho took " << ms1 << " ms" << std::endl;
                cout << "Computing the second rho took " << ms2 << " ms" << std::endl;
                
                if (mCount>1){
                cout << "BTW... the U is "<< mPhylo(curTree).toNewick() << std::endl;
                cout << "The listed leaves where"<< std::endl;
                for (set<string> setLeaves : curORho.presLeaves){
                    cout << "     ";
                    for (string l : setLeaves) cout<< l << " ";
                    cout << " " << std::endl;
                } }
            }
            
        }
        
        curTree = AllV.at(counter);
        curRho = nRho;
        curORho = nORho;
        // Save curTree so it remains unchanged after the loop
        pTree loopTree = curTree;
        //int curRho2 = curRho;
        oRho curORho2 = curORho;
        int nRho2 = -1;
        oRho nORho2 = curORho;

        // Build a list of things we can remove: leaves (as strings) and splits
        // We use a variant to hold either a string or a Split
        
        while (loopTree.rank > 2) {
    
            // Collect removal candidates
            std::vector<std::variant<std::string, Split>> candidates;
            for (const std::string& leaf : loopTree.leafSet) {
                candidates.push_back(leaf);
            }
            for (const Split& sp : loopTree.intSplits) {
                candidates.push_back(sp);
            }

            // Shuffle and try candidates until we find one that reduces rank by exactly 1
            std::shuffle(candidates.begin(), candidates.end(),
                         std::mt19937{std::random_device{}()});

            pTree curTree2;
            std::variant<std::string, Split> chosen;
            bool found = false;

            for (const auto& candidate : candidates) {
                pTree attempt;
                if (std::holds_alternative<std::string>(candidate)) {
                    attempt = loopTree.Remove(std::get<std::string>(candidate));
                } else {
                    attempt = loopTree.Remove(std::get<Split>(candidate));
                }
                if (attempt.rank == loopTree.rank - 1) {
                    curTree2 = attempt;
                    chosen = candidate;
                    found = true;
                    break;
                }
            }

            if (!found) {
                cout << "Could not find a valid removal that reduces rank by exactly 1." << std::endl;
                break;
            }

            // Compute rho two ways and time them
            auto start1 = std::chrono::high_resolution_clock::now();
            nRho2 = rho(curTree2, Tl);
            auto end1 = std::chrono::high_resolution_clock::now();

            auto start2 = std::chrono::high_resolution_clock::now();
            if (std::holds_alternative<std::string>(chosen)) {
                nORho2 = rho(curTree2, loopTree, curORho2, Tl, std::get<std::string>(chosen));
            } else {
                nORho2 = rho(curTree2, loopTree, curORho2, Tl, std::get<Split>(chosen));
            }
            auto end2 = std::chrono::high_resolution_clock::now();

            // Verify they agree
            if (nRho2 != nORho2.rho) {
                cout << "Rho's are different for " << mPhylo(curTree2).toNewick() << std::endl;
                cout << "Original rho was " << nRho2 << std::endl;
                cout << "The new rho was " << nORho2.rho << std::endl;
                cout << "The V tree is " << mPhylo(loopTree).toNewick() << std::endl;
                cout << "The listed leaves where:" << std::endl;
                for (const std::set<std::string>& setLeaves : curORho2.presLeaves) {
                    cout << "     ";
                    for (const std::string& l : setLeaves) cout << l << " ";
                    cout << std::endl;
                }
                return false;
            }

            // Report timing
            cout << "Going down from tree V = " << mPhylo(loopTree).toNewick() << std::endl;
            cout << "Towards the tree U = " << mPhylo(curTree2).toNewick() << std::endl;
            double ms1 = std::chrono::duration<double, std::milli>(end1 - start1).count();
            double ms2 = std::chrono::duration<double, std::milli>(end2 - start2).count();
            cout << "Computing the first rho took " << ms1 << " ms" << std::endl;
            cout << "Computing the second rho took " << ms2 << " ms" << std::endl;

            // Step down: the tree we just computed becomes the new starting point
            loopTree = curTree2;
            curORho2 = nORho2;
            //curRho2 = nRho2;
        }
        
    }
    
    return true;
}