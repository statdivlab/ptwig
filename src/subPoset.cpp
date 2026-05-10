#include "pTree.h"
#include "mPhylo.h"
#include "subPoset.h"
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
    
spNode::spNode(pTree eTree){
    Tree = eTree;
    next = -1;
    kappa = -1;
}

void spNode::addChild(int newChild){
    under.push_back(newChild);
}

void spNode::addParent(int newParent){
    over.push_back(newParent);
}

void spNode::setNext(int newNext){
    next = newNext;
}

void spNode::setKappa(int newKappa){
    kappa = newKappa;
}

void spNode::setLBE (float newLBE){
    lbEta = newLBE;
}

void spNode::print(){
    cout << "This node with kappa " + to_string(kappa) + " has tree \n"<< std::flush;
    Tree.print();
    cout << "\n parents: [";
    for (int i : over){
        cout << to_string(i) << ", ";
    }
    cout << "] \n"<< std::flush;
    cout << "children: [";
    for (int i : under){
        cout << to_string(i) << ", ";
    }
    cout << "] \n"<< std::flush;
    cout << "Next: "<< to_string(next) << " \n"<< std::flush;

}

void spNode::printRd(){
    cout << "This node with kappa " + to_string(kappa) + " \n"<< std::flush;
    cout << "\n parents: [";
    for (int i : over){
        cout << to_string(i) << ", ";
    }
    cout << "] \n"<< std::flush;
    cout << "children: [";
    for (int i : under){
        cout << to_string(i) << ", ";
    }
    cout << "] \n"<< std::flush;
    cout << "Next: "<< to_string(next) << " \n"<< std::flush;
}
    
subPoset::subPoset(vector<pTree> initT, vector<pTree> Sample, set<string> compLeafSet, int rb){
    int rmax = 2*static_cast<int> (compLeafSet.size()) - 7;
    firstRank = std::vector<int>(rmax+1, -1);
    lastRank = std::vector<int>(rmax+1, -1);

    Poset.push_back(spNode(initT.at(0)));
    cout << "Subposet node " << 0 << " added \n"<< std::flush;
    firstRank.at(initT.at(0).rank) = 0;
    lastRank.at(initT.at(0).rank) = 0;

    int initRank = initT.at(0).rank;
    
    Msize = 0;
    
    if (initRank == rmax){
        Msize = static_cast<int>(initT.size());
    }

    for (int i = 1; i < initT.size(); i++){
        Poset.at(lastRank.at(initRank)).setNext(i);
        lastRank.at(initRank) = i;
        Poset.push_back(spNode(initT.at(i)));
        cout << "Subposet node " << i << " added \n"<< std::flush;
    }

    //Creating things above the initial trees.
    cout << "Creating things above initial trees in SubPoset \n"<< std::flush;
    int curRank = initRank;
    int B = static_cast<int>(Sample.size());
    int curIndx = 0;


    if (initRank < rmax){
        firstRank.at(curRank + 1) = -1;
        lastRank.at(curRank + 1) = -1;
    }
    
    int nodesCount = (int)initT.size();
    
    while(curRank < rmax){
        cout << "Entering the cycle for potential node " << nodesCount << "\n"<< std::flush;
        pTree U = Poset.at(curIndx).Tree;
        vector<pTree> AllV = coverTrees(U, compLeafSet);

        auto rd = std::random_device {};
        auto rng = std::default_random_engine { rd() };
        shuffle(begin(AllV), std::end(AllV), rng);
        
        // Precompute rho(U, Z) for all sample trees once
        cout<< "Computing rho for U\n"<< std::flush;
        vector<double> rhoU(Sample.size());
        for (int i = 0; i < Sample.size(); i++){
            rhoU[i] = rho(U, Sample[i]);
        }
        
        //Version where maximal value is found.

        int IndexMaxS = -1;
        int IndexMaxL = -1;
        int IndexMax  = -1;
        float MaxValueS = 0;
        float MaxValueL = 0;
        //float MaxValue  = 0;
        int indx = 0;
        cout << "About to search through candidates for V \n"<< std::flush;
        for (pTree V : AllV){
            double tempSum = 0;
            for (int i = 0; i < Sample.size(); i++){
                // rho(U, Z) is already cached; only rho(V, Z) is computed fresh
                if ((rho(V, Sample[i]) - rhoU[i]) > 0){
                    tempSum++;
                }
            }
            if ((V.leafSet.size() > U.leafSet.size()) && (V.leafSet.size()>4)){
                if ((tempSum/B) > MaxValueL){
                    MaxValueL = tempSum/B;
                    IndexMaxL = indx;
                }
            } else {
                if ((tempSum/B) > MaxValueS){
                    MaxValueS = tempSum/B;
                    IndexMaxS = indx;
                }
            }
            indx++;
        }
        if (MaxValueL > MaxValueS){
            //MaxValue = MaxValueL;
            IndexMax = IndexMaxL;
        } else {
            //MaxValue = MaxValueS;
            IndexMax = IndexMaxS;
        }

        pTree V = AllV.at(IndexMax);

        Poset.push_back(spNode(V));
        cout << "Subposet node " << nodesCount << " added \n"<< std::flush;
        nodesCount++;

        if (V.rank == rmax){
            Msize++;
        }

        int runIndx = firstRank.at(curRank);

        while (runIndx > -1){
            if (V.covers(Poset.at(runIndx).Tree)){
                Poset.back().addChild(runIndx);
                Poset.at(runIndx).addParent(static_cast<int> (Poset.size()) - 1);
            }
            runIndx = Poset.at(runIndx).next;
        }

        if (firstRank.at(curRank+1) == -1){
            firstRank.at(curRank+1) = static_cast<int> (Poset.size()) - 1;
        }

        if (lastRank.at(curRank+1) > -1){
            Poset.at(lastRank.at(curRank+1)).next = static_cast<int> (Poset.size()) - 1;
        }

        lastRank.at(curRank+1) = static_cast<int> (Poset.size()) - 1;

        curIndx = Poset.at(curIndx).next;
        while(curIndx > -1){
            if (!Poset.at(curIndx).over.empty()){
                curIndx = Poset.at(curIndx).next;
            } else {
                break;
            }
        }

        if (curIndx == -1){
            curRank++;
            curIndx = firstRank.at(curRank);
        }
    }

    //Constructing below
    cout << "We are constructing subposet below initial trees \n"<< std::flush;
    curRank = initRank;
    curIndx = firstRank.at(curRank);

    if(curRank > 1){
        firstRank.at(curRank - 1) = -1;
        lastRank.at(curRank - 1) = -1;
    }
    int Counter2 = 0;
    while (curRank > 1) {
        pTree V = Poset.at(curIndx).Tree;

        int toAdd = 0;

        if (V.rank > rb){
            toAdd = 1 - static_cast<int> (Poset.at(curIndx).under.size());
        } else {
            toAdd = 2 - static_cast<int> (Poset.at(curIndx).under.size());
        }
        
        cout << "We are adding extra " << toAdd << " nodes \n"<< std::flush;
        if (toAdd > 0){
            // Precompute rho(V, T) for all sample trees once per node
            cout << "Precomputing rhos for V \n"<< std::flush;
            vector<double> rhoV(Sample.size());
            for (int i = 0; i < Sample.size(); i++){
                rhoV[i] = rho(V, Sample[i]);
            }

            // Unified candidate scoring: collect (stb, U) pairs for all candidates
            // then pick the 1 or 2 with lowest stb not already in under.
            float min1 = 1.2, min2 = 1.2;
            bool addU1 = false, addU2 = false;
            pTree U1, U2;

            auto considerCandidate = [&](pTree U){
                if (U.rank < V.rank - 1) return;

                for (int k1 : Poset.at(curIndx).under){
                    if (U == Poset.at(k1).Tree) return;
                }

                float Sum = 0;
                for (int i = 0; i < Sample.size(); i++){
                    // rho(V, T) is cached; only rho(U, T) is computed fresh
                    if (rhoV[i] > rho(U, Sample[i])){
                        Sum++;
                    }
                }
                float stb = Sum / B;

                if (stb < min1){
                    min2 = min1; U2 = U1; if (addU1) addU2 = true;
                    min1 = stb;  U1 = U;  addU1 = true;
                } else if (stb < min2){
                    min2 = stb; U2 = U; addU2 = true;
                }
            };
            
            cout << "Evaluating candidates to add below V \n"<< std::flush;
            for (string a : V.leafSet)    considerCandidate(V.Remove(a));
            for (Split s : V.intSplits)   considerCandidate(V.Remove(s));

            // Helper lambda to insert a new node at curRank-1
            auto insertNode = [&](pTree Unew){
                Poset.push_back(spNode(Unew));
                cout << "Subposet node " << nodesCount << " added \n"<< std::flush;
                nodesCount++;

                int runIndx = firstRank.at(curRank);
                while (runIndx > -1){
                    if (Poset.at(runIndx).Tree.covers(Unew)){
                        Poset.back().addParent(runIndx);
                        Poset.at(runIndx).addChild(static_cast<int>(Poset.size()) - 1);
                    }
                    runIndx = Poset.at(runIndx).next;
                }

                if (firstRank.at(curRank-1) == -1){
                    firstRank.at(curRank-1) = static_cast<int>(Poset.size()) - 1;
                }
                if (lastRank.at(curRank-1) > -1){
                    Poset.at(lastRank.at(curRank-1)).next = static_cast<int>(Poset.size()) - 1;
                }
                lastRank.at(curRank-1) = static_cast<int>(Poset.size()) - 1;
            };

            if (addU1) insertNode(U1);
            if (toAdd == 2 && addU2) insertNode(U2);
        }
        
        Counter2++;
        cout << "End of cycle below number " << Counter2 << "\n"<< std::flush;

        curIndx = Poset.at(curIndx).next;

        if (curIndx == -1){
            curRank--;
            curIndx = firstRank.at(curRank);
        }
    }
    
    //Assigning kappa to upper trees.
    cout << "Assigning kappas \n"<< std::flush;
    curIndx = firstRank.at(rmax);
    while(curIndx > -1){
        Poset.at(curIndx).setKappa(Msize*(static_cast<int>(Poset.at(curIndx).under.size())));
        curIndx = Poset.at(curIndx).next;
    }

    for (int j = rmax - 1; j > 0; j--){
        curIndx = firstRank.at(j);
        while(curIndx > -1){
            int kap = -1;
            for (int k : Poset.at(curIndx).over){
                if (kap < Poset.at(k).kappa){
                    kap = Poset.at(k).kappa;
                }
            }
            Poset.at(curIndx).setKappa(kap*(static_cast<int>(Poset.at(curIndx).under.size())));
            curIndx = Poset.at(curIndx).next;
        }
    }
    
    curIndx = firstRank.at(1);

    while(curIndx > -1){
        int kap = -1;
        for (int k : Poset.at(curIndx).over){
            if (kap < Poset.at(k).kappa){
                kap = Poset.at(k).kappa;
            }
        }
        Poset.at(curIndx).setKappa(kap);
        curIndx = Poset.at(curIndx).next;
    }
}

subPoset::subPoset(vector<pTree> initT, vector<pTree> Sample, vector<int> nSample, set<string> compLeafSet, int rb){
    int rmax = 2*static_cast<int> (compLeafSet.size()) - 7;
    firstRank = std::vector<int>(rmax+1, -1);
    lastRank = std::vector<int>(rmax+1, -1);
    
    Poset.push_back(spNode(initT.at(0)));
    cout << "Subposet node " << 0 << " added \n"<< std::flush;
    
    firstRank.at(initT.at(0).rank) = 0;
    lastRank.at(initT.at(0).rank) = 0;

    int initRank = initT.at(0).rank;
    
    Msize = 0;
    
    if (initRank == rmax){
        Msize = static_cast<int>(initT.size());
    }


    for (int i = 1; i < initT.size(); i++){
        Poset.at(lastRank.at(initRank)).setNext(i);
        lastRank.at(initRank) = i;
        Poset.push_back(spNode(initT.at(i)));
        cout << "Subposet node " << i  << " added \n"<< std::flush;
    }

    //Creating things above the initial trees.
    cout << "Creating things above initial trees in SubPoset \n"<< std::flush;

    int curRank = initRank;
    int B = std::accumulate(nSample.begin(), nSample.end(), 0);
    int curIndx = 0;


    if (initRank < rmax){
        firstRank.at(curRank + 1) = -1;
        lastRank.at(curRank + 1) = -1;
    }
    
    int nodesCount = (int)initT.size();
    while(curRank < rmax){
        cout << "Entering the cycle for potential node " << nodesCount << "\n"<< std::flush;
        pTree U = Poset.at(curIndx).Tree;
        vector<pTree> AllV = coverTrees(U, compLeafSet);

        auto rd = std::random_device {};
        auto rng = std::default_random_engine { rd() };
        shuffle(begin(AllV), std::end(AllV), rng);
        
        // Precompute rho(U, Z) for all sample trees once
        cout<< "Computing rho for U\n"<< std::flush;
        vector<double> rhoU(Sample.size());
        for (int i = 0; i < Sample.size(); i++){
            rhoU[i] = rho(U, Sample[i]);
        }
        
        //Version where maximal value is found.

        int IndexMaxS = -1;
        int IndexMaxL = -1;
        int IndexMax  = -1;
        float MaxValueS = 0;
        float MaxValueL = 0;
        //float MaxValue  = 0;
        int indx = 0;
        cout << "About to search through candidates for V \n"<< std::flush;
        for (pTree V : AllV){
            double tempSum = 0;
            for (int i = 0; i < Sample.size(); i++){
                // rho(U, Z) is already cached; only rho(V, Z) is computed fresh
                if ((rho(V, Sample[i]) - rhoU[i]) > 0){
                    tempSum += nSample[i];
                }
            }
            if ((V.leafSet.size() > U.leafSet.size()) && (V.leafSet.size()>4)){
                if ((tempSum/B) > MaxValueL){
                    MaxValueL = tempSum/B;
                    IndexMaxL = indx;
                }
            } else {
                if ((tempSum/B) > MaxValueS){
                    MaxValueS = tempSum/B;
                    IndexMaxS = indx;
                }
            }
            indx++;
        }
        if (MaxValueL > MaxValueS){
            //MaxValue = MaxValueL;
            IndexMax = IndexMaxL;
        } else {
            //MaxValue = MaxValueS;
            IndexMax = IndexMaxS;
        }

        pTree V = AllV.at(IndexMax);

        Poset.push_back(spNode(V));
        cout << "Subposet node " << nodesCount << " added \n"<< std::flush;
        nodesCount++;

        if (V.rank == rmax){
            Msize++;
        }

        int runIndx = firstRank.at(curRank);

        while (runIndx > -1){
            if (V.covers(Poset.at(runIndx).Tree)){
                Poset.back().addChild(runIndx);
                Poset.at(runIndx).addParent(static_cast<int> (Poset.size()) - 1);
            }
            runIndx = Poset.at(runIndx).next;
        }

        if (firstRank.at(curRank+1) == -1){
            firstRank.at(curRank+1) = static_cast<int> (Poset.size()) - 1;
        }

        if (lastRank.at(curRank+1) > -1){
            Poset.at(lastRank.at(curRank+1)).next = static_cast<int> (Poset.size()) - 1;
        }

        lastRank.at(curRank+1) = static_cast<int> (Poset.size()) - 1;

        curIndx = Poset.at(curIndx).next;
        while(curIndx > -1){
            if (!Poset.at(curIndx).over.empty()){
                curIndx = Poset.at(curIndx).next;
            } else {
                break;
            }
        }

        if (curIndx == -1){
            curRank++;
            curIndx = firstRank.at(curRank);
        }
    }

    //Constructing below
    cout << "We are constructing subposet below initial trees \n"<< std::flush;
    curRank = initRank;
    curIndx = firstRank.at(curRank);

    if(curRank > 1){
        firstRank.at(curRank - 1) = -1;
        lastRank.at(curRank - 1) = -1;

    }
    int Counter2 = 0;
    while (curRank > 1) {
        pTree V = Poset.at(curIndx).Tree;

        int toAdd = 0;

        if (V.rank > rb){
            toAdd = 1 - static_cast<int> (Poset.at(curIndx).under.size());
        } else {
            toAdd = 2 - static_cast<int> (Poset.at(curIndx).under.size());
        }
        
        cout << "We are adding extra " << toAdd << " nodes \n"<< std::flush;
        if (toAdd > 0){
            // Precompute rho(V, T) for all sample trees once per node
            cout << "Precomputing rhos for V \n"<< std::flush;
            vector<double> rhoV(Sample.size());
            for (int i = 0; i < Sample.size(); i++){
                rhoV[i] = rho(V, Sample[i]);
            }
            
            float min1 = 1.2, min2 = 1.2;
            bool addU1 = false, addU2 = false;
            pTree U1, U2;
            
            auto considerCandidate = [&](pTree U){
                if (U.rank < V.rank - 1) return;

                for (int k1 : Poset.at(curIndx).under){
                    if (U == Poset.at(k1).Tree) return;
                }

                float Sum = 0;
                for (int i = 0; i < Sample.size(); i++){
                    // rho(V, T) is cached; only rho(U, T) is computed fresh
                    if (rhoV[i] > rho(U, Sample[i])){
                        Sum += nSample[i];
                    }
                }
                float stb = Sum / B;

                if (stb < min1){
                    min2 = min1; U2 = U1; if (addU1) addU2 = true;
                    min1 = stb;  U1 = U;  addU1 = true;
                } else if (stb < min2){
                    min2 = stb; U2 = U; addU2 = true;
                }
            };
            
            cout << "Evaluating candidates to add below V \n"<< std::flush;
            for (string a : V.leafSet)    considerCandidate(V.Remove(a));
            for (Split s : V.intSplits)   considerCandidate(V.Remove(s));
            
            // Helper lambda to insert a new node at curRank-1
            auto insertNode = [&](pTree Unew){
                Poset.push_back(spNode(Unew));
                cout << "Subposet node " << nodesCount << " added \n"<< std::flush;
                nodesCount++;

                int runIndx = firstRank.at(curRank);
                while (runIndx > -1){
                    if (Poset.at(runIndx).Tree.covers(Unew)){
                        Poset.back().addParent(runIndx);
                        Poset.at(runIndx).addChild(static_cast<int>(Poset.size()) - 1);
                    }
                    runIndx = Poset.at(runIndx).next;
                }

                if (firstRank.at(curRank-1) == -1){
                    firstRank.at(curRank-1) = static_cast<int>(Poset.size()) - 1;
                }
                if (lastRank.at(curRank-1) > -1){
                    Poset.at(lastRank.at(curRank-1)).next = static_cast<int>(Poset.size()) - 1;
                }
                lastRank.at(curRank-1) = static_cast<int>(Poset.size()) - 1;
            };

            if (addU1) insertNode(U1);
            if (toAdd == 2 && addU2) insertNode(U2);

        }
        
        Counter2++;
        cout << "End of cycle below number " << Counter2 << "\n"<< std::flush;
        
        curIndx = Poset.at(curIndx).next;

        if (curIndx == -1){
            curRank--;
            curIndx = firstRank.at(curRank);
        }
    }
    
    //Assigning kappa to upper trees.
    cout << "Assigning kappas \n"<< std::flush;
    curIndx = firstRank.at(rmax);
    while(curIndx > -1){
        Poset.at(curIndx).setKappa(Msize*(static_cast<int>(Poset.at(curIndx).under.size())));
        curIndx = Poset.at(curIndx).next;
    }

    for (int j = rmax - 1; j > 0; j--){
        curIndx = firstRank.at(j);
        while(curIndx > -1){
            int kap = -1;
            for (int k : Poset.at(curIndx).over){
                if (kap < Poset.at(k).kappa){
                    kap = Poset.at(k).kappa;
                }
            }
            Poset.at(curIndx).setKappa(kap*(static_cast<int>(Poset.at(curIndx).under.size())));
            curIndx = Poset.at(curIndx).next;
        }
    }

    curIndx = firstRank.at(1);

    while(curIndx > -1){
        int kap = -1;
        for (int k : Poset.at(curIndx).over){
            if (kap < Poset.at(k).kappa){
                kap = Poset.at(k).kappa;
            }
        }
        Poset.at(curIndx).setKappa(kap);
        curIndx = Poset.at(curIndx).next;
    }
}

subPoset::subPoset(vector<pTree> Sample, set<string> compLeafSet, int Mt, int rb) {

    int B      = static_cast<int>(Sample.size());
    int maxRnk = 2 * static_cast<int>(compLeafSet.size()) - 7;

    firstRank.assign(2 * static_cast<int>(compLeafSet.size()) - 6, -1);
    lastRank.assign(2 * static_cast<int>(compLeafSet.size()) - 6, -1);

    // ---------------------------------------------------------------
    // Shared helpers
    // ---------------------------------------------------------------

    struct BeamEntry {
        pTree          tree;
        vector<float>  rhoVec;
    };

    
    struct Candidate { 
        pTree tree; 
        vector<float> rhoVec; 
        float score; };

    // Score of V relative to a cached rho baseline
    auto computeScore = [&](const vector<float>& rhoBase,
                            const vector<float>& rhoV) -> float {
        float sum = 0;
        for (int i = 0; i < B; i++)
            if (rhoV[i] - rhoBase[i] > 0)
                sum ++;
        return sum / B;
    };

    // Compute and return rho(T, Sample[i]) for all i
    auto buildRhoVec = [&](const pTree& T) -> vector<float> {
        vector<float> rv(B);
        for (int i = 0; i < B; i++)
            rv[i] = rho(T, Sample[i]);
        return rv;
    };

    auto toNwk = [&](const pTree& T) -> string {
        mPhylo mp = mPhylo(T);
        return mp.toNewick();
    };
 
    // ---------------------------------------------------------------
    // PHASE 1: Beam search upward with beam width Mt
    // ---------------------------------------------------------------

    vector<BeamEntry> currentLevel;

    // --- Seed: find best Mt trees directly above the empty tree ---
    {
        pTree emptyTree  = pTree("();");
        vector<float> rhoEmpty(B,0.0f);
        vector<pTree> above    = coverTrees(emptyTree, compLeafSet);

        vector<Candidate> candidates;
        candidates.reserve(above.size());
        
        cout << "Computing rhos for the first level in the poset \n"<< std::flush;
        for (const pTree& V : above) {
            vector<float> rv = buildRhoVec(V);
            float sc = computeScore(rhoEmpty, rv);
            candidates.push_back({V, rv, sc});
        }

        sort(candidates.begin(), candidates.end(),
             [](const Candidate& a, const Candidate& b){
                 return a.score > b.score; });

        int take = min(Mt, (int)candidates.size());
        for (int i = 0; i < take; i++)
            currentLevel.push_back({candidates[i].tree, candidates[i].rhoVec});
    }

    // --- Beam search: ascend rank by rank until maxRnk ---
    cout << "We start searching upwards for best trees at max rank \n"<< std::flush;
    int counter1 = 0;
    while ((int)currentLevel[0].tree.rank < maxRnk) {
        counter1++;
        cout << "Going up for a " << counter1 << "time \n"<< std::flush; 
        set<string>    seenNewick;
        vector<BeamEntry> futureLevel;
        futureLevel.reserve(Mt);
        
        vector<vector<Candidate>> aboveCurrents;
        
        // Step A: from each beam entry, compute aboves and order by scores
        int counterEntry = 0;
        for (const BeamEntry& entry : currentLevel) {
            vector<pTree> above = coverTrees(entry.tree, compLeafSet);
            
            vector<Candidate> aboveOne;
            cout << "Computing rhos for above trees of entry "<< counterEntry << " \n"<< std::flush;
            counterEntry++;
            for (const pTree& V : above) {
                vector<float> rv = buildRhoVec(V);
                float sc = computeScore(entry.rhoVec, rv);
                
                aboveOne.push_back({V, rv, sc});
            }
            
            sort(aboveOne.begin(), aboveOne.end(),
             [](const Candidate& a, const Candidate& b){
                 return a.score > b.score; });
            
            aboveCurrents.push_back(aboveOne);

            
        }
        
        // Step B: Fill up with aboves, avoiding repetitions. We have a running indexes for each
        cout<< "Fillin the next level with best candidates \n"<< std::flush;
        vector<int> tryIndx;
        
        tryIndx.assign(Mt, 0);
        
        int currentIndex = 0;

        while(((int)futureLevel.size() < Mt) && (accumulate(tryIndx.begin(), tryIndx.end(), 0) > -(int)currentLevel.size())){
            
            if ((tryIndx.at(currentIndex) == -1) || (tryIndx.at(currentIndex) >= aboveCurrents[currentIndex].size())){
                tryIndx.at(currentIndex) = -1;
                continue;
            }
            Candidate canU = aboveCurrents[currentIndex][tryIndx.at(currentIndex)];
            
            string nwk = toNwk(canU.tree);
            if (!seenNewick.count(nwk)){
                futureLevel.push_back({canU.tree, canU.rhoVec});
                seenNewick.insert(nwk);
            }
            tryIndx.at(currentIndex)++;
            currentIndex = (currentIndex + 1) % ((int)currentLevel.size());
        }
        /*
        if ((int)futureLevel.size() < Mt) {

            // Find best currentLevel entry by MinimumStability
            int   bestIdx = 0;
            float bestMS  = -1.0f;
            for (int ci = 0; ci < (int)currentLevel.size(); ci++) {
                float ms = MinimumStability(Sample, nSample,
                                            currentLevel[ci].rhoVec,
                                            currentLevel[ci].tree, B);
                if (ms > bestMS) {
                    bestMS  = ms;
                    bestIdx = ci;
                }
            }

            const BeamEntry& filler = currentLevel[bestIdx];
            vector<pTree> above = coverTrees(filler.tree, compLeafSet);

            struct Candidate { pTree tree; vector<float> rhoVec; float score; };
            vector<Candidate> candidates;
            candidates.reserve(above.size());

            for (const pTree& V : above) {
                string nwk = toNwk(V);
                if (seenNewick.count(nwk)) continue;

                vector<float> rv = buildRhoVec(V);
                float sc = computeScore(filler.rhoVec, rv);
                candidates.push_back({V, rv, sc});
            }

            sort(candidates.begin(), candidates.end(),
                 [](const Candidate& a, const Candidate& b){
                     return a.score > b.score; });

            for (auto& cand : candidates) {
                if ((int)futureLevel.size() >= Mt) break;
                string nwk = toNwk(cand.tree);
                if (seenNewick.count(nwk)) continue;
                seenNewick.insert(nwk);
                futureLevel.push_back({cand.tree, cand.rhoVec});
            }
        } */
        
        currentLevel = std::move(futureLevel);
    }

    /*// --- Rank currentLevel by MinimumStability, keep top Mt ---
    vector<ScoredTree> topTrees;
    topTrees.reserve(currentLevel.size());
    for (auto& entry : currentLevel) {
        float ms = MinimumStability(Sample, nSample, entry.rhoVec, entry.tree, B);
        topTrees.push_back({entry.tree, entry.rhoVec, ms});
    }
    sort(topTrees.begin(), topTrees.end(),
         [](const ScoredTree& a, const ScoredTree& b){
             return a.minStab > b.minStab; });
    if ((int)topTrees.size() > Mt)
        topTrees.resize(Mt);*/

    // ---------------------------------------------------------------
    // PHASE 2: Insert top Mt trees into the poset at maxRnk,
    //          seeding rhoCache directly from Phase 1 results
    // ---------------------------------------------------------------

    vector<vector<float>> rhoCache;
    rhoCache.reserve(Mt);
    int nodesCount = 0;
    for (int m = 0; m < (int)currentLevel.size(); m++) {
        Poset.push_back(spNode(currentLevel[m].tree));
        cout << "Subposet node " << nodesCount << " added \n"<< std::flush;
        nodesCount++;
        rhoCache.push_back(currentLevel[m].rhoVec);   // no recomputation

        if (m < (int)currentLevel.size() - 1)
            Poset[m].setNext(m + 1);
        // last node: next stays -1 from spNode constructor

        // Maximal nodes: kappa placeholder — corrected in Phase 3 kappa pass
        Poset[m].setKappa(Mt);
    }

    firstRank[maxRnk] = 0;
    lastRank[maxRnk] = (int)currentLevel.size() - 1;
    

    // ---------------------------------------------------------------
    // PHASE 3: Build downward iteratively, rank by rank
    // ---------------------------------------------------------------

    // Helper: insert a child U into the poset at curRank-1,
    //         wiring all parent edges from curRank
    auto insertChild = [&](const pTree& U, const vector<float>& rhoU,
                           int curRank) {

        // Check if U already exists at curRank-1
        int existIdx = -1;
        int runIndx  = firstRank[curRank - 1];
        while (runIndx > -1) {
            if (Poset[runIndx].Tree == U) {
                existIdx = runIndx;
                break;
            }
            runIndx = Poset[runIndx].next;
        }

        int childIdx;
        if (existIdx > -1) {
            childIdx = existIdx;
        } else {
            Poset.push_back(spNode(U));
            cout << "Subposet node " << nodesCount << " added \n"<< std::flush;
            nodesCount++;
            childIdx = (int)Poset.size() - 1;
            rhoCache.push_back(rhoU);

            if (firstRank[curRank - 1] == -1)
                firstRank[curRank - 1] = childIdx;
            if (lastRank[curRank - 1] > -1)
                Poset[lastRank[curRank - 1]].setNext(childIdx);
            lastRank[curRank - 1] = childIdx;
        }

        // Wire parent edges: scan all nodes at curRank
        runIndx = firstRank[curRank];
        while (runIndx > -1) {
            if (Poset[runIndx].Tree.covers(U)) {
                bool alreadyLinked = false;
                for (int k : Poset[childIdx].over)
                    if (k == runIndx) { alreadyLinked = true; break; }
                if (!alreadyLinked) {
                    Poset[childIdx].addParent(runIndx);
                    Poset[runIndx].addChild(childIdx);
                }
            }
            runIndx = Poset[runIndx].next;
        }
    };
    int curRank = maxRnk;
    int curIndx = firstRank[curRank];

    if (curRank > 1) {
        firstRank[curRank - 1] = -1;
        lastRank[curRank - 1] = -1;
    }
    while (curRank > 1) {

        pTree          V    = Poset[curIndx].Tree;
        vector<float>  rhoV = rhoCache[curIndx];

        int toAdd = (V.rank > rb ? 1 : 2)
                    - (int)Poset[curIndx].under.size();

        if (toAdd > 0) {

            // Evaluate all candidate children, using cached rhoV
            vector<Candidate> candidates;

            auto evalFeature = [&](pTree U) {
                if (U.rank < V.rank - 1) return;
                for (int k : Poset[curIndx].under)
                    if (U == Poset[k].Tree) return;

                vector<float> ru = buildRhoVec(U);
                float sum = 0;
                for (int i = 0; i < B; i++)
                    if (rhoV[i] - ru[i] > 0)
                        sum++;
                candidates.push_back({U, ru, sum / B});
            };
            
            cout << "Evaluating candidates below \n"<< std::flush;
            for (string a : V.leafSet)   evalFeature(V.Remove(a));
            for (Split  s : V.intSplits) evalFeature(V.Remove(s));

            // Sort by stability ascending: lowest stability first
            sort(candidates.begin(), candidates.end(),
                 [](const Candidate& a, const Candidate& b){
                     return a.score < b.score; });

            for (int t = 0; t < min(toAdd, (int)candidates.size()); t++)
                insertChild(candidates[t].tree, candidates[t].rhoVec, curRank);
        }

        curIndx = Poset[curIndx].next;

        // Finished all nodes at curRank: run kappa pass then step down
        if (curIndx == -1) {

            int scanIndx = firstRank[curRank];
            while (scanIndx > -1) {
                spNode& node = Poset[scanIndx];
                if (node.over.empty()) {
                    node.setKappa(Mt * (int)node.under.size());
                } else {
                    int maxKappa = 0;
                    for (int k : node.over)
                        maxKappa = max(maxKappa, Poset[k].kappa);
                    node.setKappa(maxKappa * (int)node.under.size());
                }
                scanIndx = node.next;
            }

            curRank--;
            //if (curRank > 1) {
            //    firstRank[curRank - 1] = -1;
            //    lastRank [curRank - 1] = -1;
            //}
            curIndx = firstRank[curRank];
        }
    }
    // Kappa pass for rank 1
    cout << "Assigning kappa values \n"<< std::flush;
    if (curRank == 1 && firstRank[1] > -1) {
        int scanIndx = firstRank[1];
        while (scanIndx > -1) {
            spNode& node = Poset[scanIndx];
            int maxKappa = 0;
            for (int k : node.over)
                maxKappa = max(maxKappa, Poset[k].kappa);
            node.setKappa(maxKappa);
            scanIndx = node.next;
        }
    }
}


subPoset::subPoset(vector<pTree> Sample, vector<int> nSample,
                   set<string> compLeafSet, int Mt, int rb) {
    
    int N      = static_cast<int>(Sample.size());
    int B      = accumulate(nSample.begin(), nSample.end(), 0);
    int maxRnk = 2 * static_cast<int>(compLeafSet.size()) - 7;

    firstRank.assign(2 * static_cast<int>(compLeafSet.size()) - 6, -1);
    lastRank.assign(2 * static_cast<int>(compLeafSet.size()) - 6, -1);
    
    // ---------------------------------------------------------------
    // Shared helpers
    // ---------------------------------------------------------------

    struct BeamEntry {
        pTree          tree;
        vector<float>  rhoVec;
    };

    
    struct Candidate { 
        pTree tree; 
        vector<float> rhoVec; 
        float score; };

    // Score of V relative to a cached rho baseline
    auto computeScore = [&](const vector<float>& rhoBase,
                            const vector<float>& rhoV) -> float {
        float sum = 0;
        for (int i = 0; i < N; i++)
            if (rhoV[i] - rhoBase[i] > 0)
                sum += nSample[i];
        return sum / B;
    };

    // Compute and return rho(T, Sample[i]) for all i
    auto buildRhoVec = [&](const pTree& T) -> vector<float> {
        vector<float> rv(N);
        for (int i = 0; i < N; i++)
            rv[i] = rho(T, Sample[i]);
        return rv;
    };

    auto toNwk = [&](const pTree& T) -> string {
        mPhylo mp = mPhylo(T);
        return mp.toNewick();
    };

    // ---------------------------------------------------------------
    // PHASE 1: Beam search upward with beam width Mt
    // ---------------------------------------------------------------

    vector<BeamEntry> currentLevel;

    // --- Seed: find best Mt trees directly above the empty tree ---
    {
        pTree emptyTree  = pTree("();");
        vector<float> rhoEmpty = buildRhoVec(emptyTree);
        vector<pTree> above    = coverTrees(emptyTree, compLeafSet);

        vector<Candidate> candidates;
        candidates.reserve(above.size());
        
        cout << "Computing rhos for the first level in the poset \n"<< std::flush;
        for (const pTree& V : above) {
            vector<float> rv = buildRhoVec(V);
            float sc = computeScore(rhoEmpty, rv);
            candidates.push_back({V, rv, sc});
        }

        sort(candidates.begin(), candidates.end(),
             [](const Candidate& a, const Candidate& b){
                 return a.score > b.score; });

        int take = min(Mt, (int)candidates.size());
        for (int i = 0; i < take; i++)
            currentLevel.push_back({candidates[i].tree, candidates[i].rhoVec});
    }

    // --- Beam search: ascend rank by rank until maxRnk ---
    cout << "We start searching upwards for best trees at max rank \n"<< std::flush;
    int counter1 = 0;
    while ((int)currentLevel[0].tree.rank < maxRnk) {
        counter1++;
        cout << "Going up for a " << counter1 << "time \n"<< std::flush; 
        set<string>    seenNewick;
        vector<BeamEntry> futureLevel;
        futureLevel.reserve(Mt);
        
        vector<vector<Candidate>> aboveCurrents;
        
        // Step A: from each beam entry, compute aboves and order by scores
        int counterEntry = 0;
        for (const BeamEntry& entry : currentLevel) {
            vector<pTree> above = coverTrees(entry.tree, compLeafSet);
            
            vector<Candidate> aboveOne;
            
            cout << "Computing rhos for above trees of entry "<< counterEntry << " \n"<< std::flush;
            counterEntry++;
            for (const pTree& V : above) {
                vector<float> rv = buildRhoVec(V);
                float sc = computeScore(entry.rhoVec, rv);
                
                aboveOne.push_back({V, rv, sc});
            }
            
            sort(aboveOne.begin(), aboveOne.end(),
             [](const Candidate& a, const Candidate& b){
                 return a.score > b.score; });
            
            aboveCurrents.push_back(aboveOne);

            
        }
        
        // Step B: Fill up with aboves, avoiding repetitions. We have a running indexes for each
        cout<< "Fillin the next level with best candidates \n"<< std::flush;
        vector<int> tryIndx;
        
        tryIndx.assign(Mt, 0);
        
        int currentIndex = 0;
        
        while(((int)futureLevel.size() < Mt) && (accumulate(tryIndx.begin(), tryIndx.end(), 0) > -(int)currentLevel.size())){
            
            if ((tryIndx.at(currentIndex) == -1) || (tryIndx.at(currentIndex) >= aboveCurrents[currentIndex].size())){
                tryIndx.at(currentIndex) = -1;
                continue;
            }
            Candidate canU = aboveCurrents[currentIndex][tryIndx.at(currentIndex)];
            
            string nwk = toNwk(canU.tree);
            if (!seenNewick.count(nwk)){
                futureLevel.push_back({canU.tree, canU.rhoVec});
                seenNewick.insert(nwk);
            }
            tryIndx.at(currentIndex)++;
            currentIndex = (currentIndex + 1) % ((int)currentLevel.size());
        }
        /*
        if ((int)futureLevel.size() < Mt) {

            // Find best currentLevel entry by MinimumStability
            int   bestIdx = 0;
            float bestMS  = -1.0f;
            for (int ci = 0; ci < (int)currentLevel.size(); ci++) {
                float ms = MinimumStability(Sample, nSample,
                                            currentLevel[ci].rhoVec,
                                            currentLevel[ci].tree, B);
                if (ms > bestMS) {
                    bestMS  = ms;
                    bestIdx = ci;
                }
            }

            const BeamEntry& filler = currentLevel[bestIdx];
            vector<pTree> above = coverTrees(filler.tree, compLeafSet);

            struct Candidate { pTree tree; vector<float> rhoVec; float score; };
            vector<Candidate> candidates;
            candidates.reserve(above.size());

            for (const pTree& V : above) {
                string nwk = toNwk(V);
                if (seenNewick.count(nwk)) continue;

                vector<float> rv = buildRhoVec(V);
                float sc = computeScore(filler.rhoVec, rv);
                candidates.push_back({V, rv, sc});
            }

            sort(candidates.begin(), candidates.end(),
                 [](const Candidate& a, const Candidate& b){
                     return a.score > b.score; });

            for (auto& cand : candidates) {
                if ((int)futureLevel.size() >= Mt) break;
                string nwk = toNwk(cand.tree);
                if (seenNewick.count(nwk)) continue;
                seenNewick.insert(nwk);
                futureLevel.push_back({cand.tree, cand.rhoVec});
            }
        } */
        
        currentLevel = std::move(futureLevel);
    }

    /*// --- Rank currentLevel by MinimumStability, keep top Mt ---
    vector<ScoredTree> topTrees;
    topTrees.reserve(currentLevel.size());
    for (auto& entry : currentLevel) {
        float ms = MinimumStability(Sample, nSample, entry.rhoVec, entry.tree, B);
        topTrees.push_back({entry.tree, entry.rhoVec, ms});
    }
    sort(topTrees.begin(), topTrees.end(),
         [](const ScoredTree& a, const ScoredTree& b){
             return a.minStab > b.minStab; });
    if ((int)topTrees.size() > Mt)
        topTrees.resize(Mt);*/

    // ---------------------------------------------------------------
    // PHASE 2: Insert top Mt trees into the poset at maxRnk,
    //          seeding rhoCache directly from Phase 1 results
    // ---------------------------------------------------------------

    vector<vector<float>> rhoCache;
    rhoCache.reserve(Mt);
    
    int nodesCount = 0;
    for (int m = 0; m < (int)currentLevel.size(); m++) {
        Poset.push_back(spNode(currentLevel[m].tree));
        cout << "Subposet node " << nodesCount << " added \n"<< std::flush;
        nodesCount++;
        rhoCache.push_back(currentLevel[m].rhoVec);   // no recomputation

        if (m < (int)currentLevel.size() - 1)
            Poset[m].setNext(m + 1);
        // last node: next stays -1 from spNode constructor

        // Maximal nodes: kappa placeholder — corrected in Phase 3 kappa pass
        Poset[m].setKappa(Mt);
    }

    firstRank[maxRnk] = 0;
    lastRank[maxRnk] = (int)currentLevel.size() - 1;

    // ---------------------------------------------------------------
    // PHASE 3: Build downward iteratively, rank by rank
    // ---------------------------------------------------------------

    // Helper: insert a child U into the poset at curRank-1,
    //         wiring all parent edges from curRank
    auto insertChild = [&](const pTree& U, const vector<float>& rhoU,
                           int curRank) {

        // Check if U already exists at curRank-1
        int existIdx = -1;
        int runIndx  = firstRank[curRank - 1];
        while (runIndx > -1) {
            if (Poset[runIndx].Tree == U) {
                existIdx = runIndx;
                break;
            }
            runIndx = Poset[runIndx].next;
        }

        int childIdx;
        if (existIdx > -1) {
            childIdx = existIdx;
        } else {
            Poset.push_back(spNode(U));
            cout << "Subposet node " << nodesCount << " added \n"<< std::flush;
            nodesCount++;
            childIdx = (int)Poset.size() - 1;
            rhoCache.push_back(rhoU);

            if (firstRank[curRank - 1] == -1)
                firstRank[curRank - 1] = childIdx;
            if (lastRank[curRank - 1] > -1)
                Poset[lastRank[curRank - 1]].setNext(childIdx);
            lastRank[curRank - 1] = childIdx;
        }

        // Wire parent edges: scan all nodes at curRank
        runIndx = firstRank[curRank];
        while (runIndx > -1) {
            if (Poset[runIndx].Tree.covers(U)) {
                bool alreadyLinked = false;
                for (int k : Poset[childIdx].over)
                    if (k == runIndx) { alreadyLinked = true; break; }
                if (!alreadyLinked) {
                    Poset[childIdx].addParent(runIndx);
                    Poset[runIndx].addChild(childIdx);
                }
            }
            runIndx = Poset[runIndx].next;
        }
    };
    int curRank = maxRnk;
    int curIndx = firstRank[curRank];

    if (curRank > 1) {
        firstRank[curRank - 1] = -1;
        lastRank[curRank - 1] = -1;
    }
    while (curRank > 1) {

        pTree          V    = Poset[curIndx].Tree;
        vector<float>  rhoV = rhoCache[curIndx];

        int toAdd = (V.rank > rb ? 1 : 2)
                    - (int)Poset[curIndx].under.size();

        if (toAdd > 0) {

            // Evaluate all candidate children, using cached rhoV
            vector<Candidate> candidates;

            auto evalFeature = [&](pTree U) {
                if (U.rank < V.rank - 1) return;
                for (int k : Poset[curIndx].under)
                    if (U == Poset[k].Tree) return;

                vector<float> ru = buildRhoVec(U);
                float sum = 0;
                for (int i = 0; i < N; i++)
                    if (rhoV[i] - ru[i] > 0)
                        sum += nSample[i];
                candidates.push_back({U, ru, sum / B});
            };
            
            cout << "Evaluating candidates below \n"<< std::flush;
            for (string a : V.leafSet)   evalFeature(V.Remove(a));
            for (Split  s : V.intSplits) evalFeature(V.Remove(s));

            // Sort by stability ascending: lowest stability first
            sort(candidates.begin(), candidates.end(),
                 [](const Candidate& a, const Candidate& b){
                     return a.score < b.score; });

            for (int t = 0; t < min(toAdd, (int)candidates.size()); t++)
                insertChild(candidates[t].tree, candidates[t].rhoVec, curRank);
        }

        curIndx = Poset[curIndx].next;

        // Finished all nodes at curRank: run kappa pass then step down
        if (curIndx == -1) {

            int scanIndx = firstRank[curRank];
            while (scanIndx > -1) {
                spNode& node = Poset[scanIndx];
                if (node.over.empty()) {
                    node.setKappa(Mt * (int)node.under.size());
                } else {
                    int maxKappa = 0;
                    for (int k : node.over)
                        maxKappa = max(maxKappa, Poset[k].kappa);
                    node.setKappa(maxKappa * (int)node.under.size());
                }
                scanIndx = node.next;
            }

            curRank--;
            //if (curRank > 1) {
            //    firstRank[curRank - 1] = -1;
            //    lastRank [curRank - 1] = -1;
            //}
            curIndx = firstRank[curRank];
        }
    }
    // Kappa pass for rank 1
    cout << "Assigning kappa values \n"<< std::flush;
    if (curRank == 1 && firstRank[1] > -1) {
        int scanIndx = firstRank[1];
        while (scanIndx > -1) {
            spNode& node = Poset[scanIndx];
            int maxKappa = 0;
            for (int k : node.over)
                maxKappa = max(maxKappa, Poset[k].kappa);
            node.setKappa(maxKappa);
            scanIndx = node.next;
        }
    }
}


void subPoset::print(){
    for (int j = 1; j < firstRank.size(); j++){
        int curIndx = firstRank.at(j);
        while (curIndx > -1){
            cout << "\n For node " << to_string(curIndx) << "\n"<< std::flush;
            Poset.at(curIndx).print();
            curIndx = Poset.at(curIndx).next;
        }
    }

    cout<<"Additional information: \n"<< std::flush;
    cout<<"    First in ranks: ";
    for (int j : firstRank){
        cout<< to_string(j) << " ";
    }

    cout<<"\n    Last in ranks: ";
    for (int j : lastRank){
        cout<< to_string(j) << " ";
    }
    cout << "\n"<< std::flush;
}

void subPoset::printRd(){
    for (int j = 1; j < firstRank.size(); j++){
        int curIndx = firstRank.at(j);
        while (curIndx > -1){
            cout << "\n For node " << to_string(curIndx) << "\n"<< std::flush;
            Poset.at(curIndx).printRd();
            curIndx = Poset.at(curIndx).next;
        }
    }

    cout<<"Additional information: \n"<< std::flush;
    cout<<"    First in ranks: ";
    for (int j : firstRank){
        cout<< to_string(j) << " ";
    }

    cout<<"\n    Last in ranks: ";
    for (int j : lastRank){
        cout<< to_string(j) << " ";
    }
    cout << "\n"<< std::flush;
}

