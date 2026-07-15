#include "pTree.h"
#include "mPhylo.h"
#include "subPoset.h"
#include "rho.h"
#include "coverTrees.h"
#include <rcpptimer.h>
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
#include <unordered_map>
#include <algorithm>
#include <filesystem> // C++17
#include <future>
#include <mutex>
#include <cmath>

using namespace Rcpp;
using namespace std;
    
spNode::spNode(pTree eTree){
    Tree = eTree;
    next = -1;
    zeta = 0;
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

void spNode::setZeta(float newZ){
    zeta = newZ;
}

void spNode::print(){
    cout << "This node has tree \n"<< std::flush;
    cout << "This node has zeta " << zeta << "\n" << std::flush;
    cout << mPhylo(Tree).toNewick() << "\n" << std::flush;
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
    cout << "Nu for those children: [";
    for (int nu : boundAntichain){
        cout << to_string(nu) << ", ";
    }
    cout << "] \n"<< std::flush;
    cout << "Next: "<< to_string(next) << " \n"<< std::flush;

}

void spNode::printRd(){
    cout << "This node \n"<< std::flush;
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
    

subPoset::subPoset(vector<pTree> initT, vector<pTree> Sample, vector<int> nSample, set<string> compLeafSet, int rb){
    int rmax = 2*static_cast<int> (compLeafSet.size()) - 7;
    firstRank = std::vector<int>(rmax+1, -1);
    lastRank = std::vector<int>(rmax+1, -1);
    
    Poset.push_back(spNode(initT.at(0)));
    
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
    }

    //Creating things above the initial trees.

    int curRank = initRank;
    int B = std::accumulate(nSample.begin(), nSample.end(), 0);
    int curIndx = 0;


    if (initRank < rmax){
        firstRank.at(curRank + 1) = -1;
        lastRank.at(curRank + 1) = -1;
    }
    
    int nodesCount = (int)initT.size();
    while(curRank < rmax){
        pTree U = Poset.at(curIndx).Tree;
        vector<pTree> AllV = coverTrees(U, compLeafSet);

        auto rd = std::random_device {};
        auto rng = std::default_random_engine { rd() };
        shuffle(begin(AllV), std::end(AllV), rng);
        
        // Precompute rho(U, Z) for all sample trees once
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
    curRank = initRank;
    curIndx = firstRank.at(curRank);

    if(curRank > 1){
        firstRank.at(curRank - 1) = -1;
        lastRank.at(curRank - 1) = -1;

    }
    //int Counter2 = 0;
    while (curRank > 1) {
        pTree V = Poset.at(curIndx).Tree;

        int toAdd = 0;

        if (V.rank > rb){
            toAdd = 1 - static_cast<int> (Poset.at(curIndx).under.size());
        } else {
            toAdd = 2 - static_cast<int> (Poset.at(curIndx).under.size());
        }
        
        if (toAdd > 0){
            // Precompute rho(V, T) for all sample trees once per node
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
            
            for (string a : V.leafSet)    considerCandidate(V.Remove(a));
            for (Split s : V.intSplits)   considerCandidate(V.Remove(s));
            
            // Helper lambda to insert a new node at curRank-1
            auto insertNode = [&](pTree Unew){
                Poset.push_back(spNode(Unew));
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
        
        //Counter2++;
        
        curIndx = Poset.at(curIndx).next;

        if (curIndx == -1){
            curRank--;
            curIndx = firstRank.at(curRank);
        }
    }
    
    // Zeta pass: for each node, zeta is the mean over its covers of the fraction
    // of (weighted) treeSample1 trees whose rho strictly increases. This matches
    // the definition the fixed-width / basic_score constructor uses.
    for (int idx = 0; idx < static_cast<int>(Poset.size()); idx++){
        pTree Ta = Poset.at(idx).Tree;

        vector<double> rhoTa(Sample.size());
        for (int i = 0; i < static_cast<int>(Sample.size()); i++){
            rhoTa[i] = rho(Ta, Sample[i]);
        }

        vector<pTree> covers = coverTrees(Ta, compLeafSet);
        float zSum = 0.0f;
        int zCnt = 0;
        for (const pTree& Tb : covers){
            double improving = 0;
            for (int i = 0; i < static_cast<int>(Sample.size()); i++){
                if (rho(Tb, Sample[i]) > rhoTa[i]){
                    improving += nSample[i];
                }
            }
            zSum += static_cast<float>(improving / B);
            zCnt++;
        }
        Poset.at(idx).setZeta(zCnt > 0 ? zSum / zCnt : 0.0f);
    }
}



subPoset::subPoset(vector<pTree> Sample, vector<int> nSample,
                   set<string> compLeafSet, int Mt, int rb) {
    
    Rcpp::Timer timer("times_basic");
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
        vector<oRho>  rhoVec;
    };

    
    struct Candidate { 
        pTree tree; 
        vector<oRho> rhoVec; 
        float score; };

    // Score of V relative to a cached rho baseline
    auto computeScore = [&](const vector<oRho>& rhoBase,
                            const vector<oRho>& rhoV) -> float {
        Rcpp::Timer::ScopedTimer st(timer, "Basic_computeScore");
        float sum = 0;
        for (int i = 0; i < N; i++)
            if (rhoV[i].rho - rhoBase[i].rho > 0)
                sum += nSample[i];
        return sum / B;
    };

    // Compute and return rho(T, Sample[i]) for all i
    auto buildRhoVec = [&](const pTree& T, const pTree& TB, const vector<oRho>& rhoBase) -> vector<oRho> {
        Rcpp::Timer::ScopedTimer st(timer, "Basic_buildRhoVec");
        vector<oRho> rv(N);
        for (int i = 0; i < N; i++)
            rv[i] = rho(T, TB, rhoBase[i], Sample[i]);
        return rv;
    };
    
    auto buildRhoVecL = [&](const pTree& T, const pTree& TB, const vector<oRho>& rhoBase, const string& a) -> vector<oRho> {
        Rcpp::Timer::ScopedTimer st(timer, "Basic_buildRhoVecL");
        vector<oRho> rv(N);
        for (int i = 0; i < N; i++)
            rv[i] = rho(T, TB, rhoBase[i], Sample[i], a);
        return rv;
    };
    
    auto buildRhoVecS = [&](const pTree& T, const pTree& TB, const vector<oRho>& rhoBase, const Split& s) -> vector<oRho> {
        Rcpp::Timer::ScopedTimer st(timer, "Basic_buildRhoVecS");
        vector<oRho> rv(N);
        for (int i = 0; i < N; i++)
            rv[i] = rho(T, TB, rhoBase[i], Sample[i], s);
        return rv;
    };


    auto toNwk = [&](const pTree& T) -> string {
        Rcpp::Timer::ScopedTimer st(timer, "Basic_toNwk");
        mPhylo mp = mPhylo(T);
        return mp.toNewick();
    };

    // ---------------------------------------------------------------
    // PHASE 1: Beam search upward with beam width Mt
    // ---------------------------------------------------------------

    vector<BeamEntry> currentLevel;

    // --- Seed: find best Mt trees directly above the empty tree ---
    timer.tic("Basic_firstlevel");
    {
        pTree emptyTree  = pTree("();");
        //vector<set<string>> emptyLeaves;
        oRho curORho = oRho();
        vector<oRho> rhoEmpty(N, curORho);
        vector<pTree> above    = coverTrees(emptyTree, compLeafSet);

        vector<Candidate> candidates;
        candidates.reserve(above.size());
        
        for (const pTree& V : above) {
            vector<oRho> rv = buildRhoVec(V, emptyTree, rhoEmpty);
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
    timer.toc("Basic_firstlevel");

    // --- Beam search: ascend rank by rank until maxRnk ---
    timer.tic("Basic_beamUpwards");
    while ((int)currentLevel[0].tree.rank < maxRnk) {
        
        set<string>    seenNewick;
        vector<BeamEntry> futureLevel;
        futureLevel.reserve(Mt);
        
        vector<vector<Candidate>> aboveCurrents;
        
        // Step A: from each beam entry, compute aboves and order by scores
        for (const BeamEntry& entry : currentLevel) {
            vector<pTree> above = coverTrees(entry.tree, compLeafSet);
            
            vector<Candidate> aboveOne;
            
            for (const pTree& V : above) {
                vector<oRho> rv = buildRhoVec(V, entry.tree, entry.rhoVec);
                float sc = computeScore(entry.rhoVec, rv);
                
                aboveOne.push_back({V, rv, sc});
            }
            
            sort(aboveOne.begin(), aboveOne.end(),
             [](const Candidate& a, const Candidate& b){
                 return a.score > b.score; });
            
            aboveCurrents.push_back(aboveOne);

            
        }
        
        // Step B: Fill up with aboves, avoiding repetitions. We have a running indexes for each
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
        
        currentLevel = std::move(futureLevel);
    }

    timer.toc("Basic_beamUpwards");
    // ---------------------------------------------------------------
    // PHASE 2: Insert top Mt trees into the poset at maxRnk,
    //          seeding rhoCache directly from Phase 1 results
    // ---------------------------------------------------------------
    
    timer.tic("Basic_InsertingTop");
    vector<vector<oRho>> rhoCache;
    rhoCache.reserve(Mt);
    
    for (int m = 0; m < (int)currentLevel.size(); m++) {
        Poset.push_back(spNode(currentLevel[m].tree));
        rhoCache.push_back(currentLevel[m].rhoVec);   // no recomputation

        if (m < (int)currentLevel.size() - 1)
            Poset[m].setNext(m + 1);
        // last node: next stays -1 from spNode constructor

    }

    firstRank[maxRnk] = 0;
    lastRank[maxRnk] = (int)currentLevel.size() - 1;
    timer.toc("Basic_InsertingTop");

    // ---------------------------------------------------------------
    // PHASE 3: Build downward iteratively, rank by rank
    // ---------------------------------------------------------------

    // Helper: insert a child U into the poset at curRank-1,
    //         wiring all parent edges from curRank
    auto insertChild = [&](const pTree& U, const vector<oRho>& rhoU,
                           int curRank) {
        Rcpp::Timer::ScopedTimer st(timer, "Basic_insertChild");
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
    
    timer.tic("Basic_buildingDownwards");
    int curRank = maxRnk;
    int curIndx = firstRank[curRank];

    if (curRank > 1) {
        firstRank[curRank - 1] = -1;
        lastRank[curRank - 1] = -1;
    }
    while (curRank > 1) {

        pTree          V    = Poset[curIndx].Tree;
        vector<oRho>  rhoV = rhoCache[curIndx];

        int toAdd = (V.rank > rb ? 1 : 2)
                    - (int)Poset[curIndx].under.size();

        if (toAdd > 0) {

            // Evaluate all candidate children, using cached rhoV
            vector<Candidate> candidates;

            auto evalFeatureL = [&](pTree U, string a) {
                Rcpp::Timer::ScopedTimer st(timer, "Basic_evalFeatureL");
                if (U.rank < V.rank - 1) return;
                for (int k : Poset[curIndx].under)
                    if (U == Poset[k].Tree) return;

                vector<oRho> ru = buildRhoVecL(U, V, rhoV, a);
                float sum = 0;
                for (int i = 0; i < N; i++)
                    if (rhoV[i].rho - ru[i].rho > 0)
                        sum += nSample[i];
                candidates.push_back({U, ru, sum / B});
            };
            
            auto evalFeatureS = [&](pTree U, Split s) {
                Rcpp::Timer::ScopedTimer st(timer, "Basic_evalFeatureS");
                if (U.rank < V.rank - 1) return;
                for (int k : Poset[curIndx].under)
                    if (U == Poset[k].Tree) return;

                vector<oRho> ru = buildRhoVecS(U, V, rhoV, s);
                float sum = 0;
                for (int i = 0; i < N; i++)
                    if (rhoV[i].rho - ru[i].rho > 0)
                        sum += nSample[i];
                candidates.push_back({U, ru, sum / B});
            };
            
            for (string a : V.leafSet)   evalFeatureL(V.Remove(a), a);
            for (Split  s : V.intSplits) evalFeatureS(V.Remove(s), s);

            // Sort by stability ascending: lowest stability first
            sort(candidates.begin(), candidates.end(),
                 [](const Candidate& a, const Candidate& b){
                     return a.score < b.score; });

            for (int t = 0; t < min(toAdd, (int)candidates.size()); t++)
                insertChild(candidates[t].tree, candidates[t].rhoVec, curRank);
        }

        curIndx = Poset[curIndx].next;

        // Finished all nodes at curRank: step down
        if (curIndx == -1) {
            curRank--;
            curIndx = firstRank[curRank];
        }
    }
    timer.toc("Basic_buildingDownwards");

    // Zeta pass: for each node, zeta is the mean over its covers of the fraction
    // of (weighted) treeSample1 trees whose rho strictly increases - the same
    // definition the fixed-width / basic_score constructor uses. rhoCache[idx]
    // already holds rho(node, Sample[i]), so we reuse it.
    for (int idx = 0; idx < static_cast<int>(Poset.size()); idx++){
        vector<pTree> covers = coverTrees(Poset[idx].Tree, compLeafSet);
        float zSum = 0.0f;
        int zCnt = 0;
        for (const pTree& Tb : covers){
            vector<oRho> rvb = buildRhoVec(Tb, Poset[idx].Tree, rhoCache[idx]);
            zSum += computeScore(rhoCache[idx], rvb);
            zCnt++;
        }
        Poset[idx].setZeta(zCnt > 0 ? zSum / zCnt : 0.0f);
    }
}

// New builder:
// subPoset::subPoset(vector<pTree> Sample, vector<int> nSample,
//                    set<string> compLeafSet, int Mt, int rb, float q)
//
// Add to subPoset.h:
//   subPoset(std::vector<pTree> Sample, std::vector<int> nSample,
//            std::set<std::string> compLeafSet, int Mt, int rb, float q);
//
// Also add to pTree.h and pTree.cpp:
//   bool covers(pTree tOther) const;

subPoset::subPoset(vector<pTree> Sample, vector<int> nSample,
                   set<string> compLeafSet, float q)
{   //static std::ofstream //l ogFile("/tmp/rho_timing.log", std::ios::app);
    Rcpp::Timer timer("times_upwards");
    int N         = static_cast<int>(Sample.size());
    int B         = accumulate(nSample.begin(), nSample.end(), 0);
    int maxRnkAbs = 2 * static_cast<int>(compLeafSet.size()) - 7; // absolute max possible rank

    firstRank.assign(maxRnkAbs + 1, -1);
    lastRank.assign(maxRnkAbs + 1, -1);
    Msize = 0;

    // current max rank in the subposet being built — updated as nodes are added
    int curMaxRnk = 0;

    // ---------------------------------------------------------------
    // Shared types
    // ---------------------------------------------------------------

    struct NodeEntry {
        pTree        tree;
        vector<oRho> rhoVec;   // rho(tree, Sample[i]) for all i
        int          posetIdx; // index in Poset[]
        float        zeta;     // zeta(tree), estimated at insertion time
    };

    // rhoCache and zetaCache mirror Poset indices
    vector<vector<oRho>> rhoCache;
    vector<float>        zetaCache;

    // RNG shared across all shuffles
    auto rd  = random_device{};
    auto rng = default_random_engine{rd()};

    static constexpr int ZETA_SAMPLE = 30; // max trees sampled to estimate zeta

    // ---------------------------------------------------------------
    // Shared helpers
    // ---------------------------------------------------------------

    auto buildRhoVec = [&](const pTree& T, const pTree& TB,
                           const vector<oRho>& rhoBase) -> vector<oRho> {
        Rcpp::Timer::ScopedTimer st(timer, "Upwards_buildRhoVec");
        vector<oRho> rv(N);
        for (int i = 0; i < N; i++){
            //auto t0 = chrono::high_resolution_clock::now();
            rv[i] = rho(T, TB, rhoBase[i], Sample[i]);
            //auto t1 = chrono::high_resolution_clock::now();
            //double totalBuildRho = chrono::duration<double, milli>(t1 - t0).count();
            //l ogFile << "          One takes" << totalBuildRho << " ms \n" << std::flush;
            
        }
        return rv;
    };

    auto computeScore = [&](const vector<oRho>& rhoBase,
                            const vector<oRho>& rhoV) -> float {
        Rcpp::Timer::ScopedTimer st(timer, "Upwards_computeScore");
        float sum = 0;
        for (int i = 0; i < N; i++)
            if (rhoV[i].rho - rhoBase[i].rho > 0)
                sum += nSample[i];
        return sum / B;
    };

    // Estimate zeta for Ta using a random sample of at most ZETA_SAMPLE
    // cover trees. allCovers is the full (shuffled) list of cover trees —
    // we take the first min(ZETA_SAMPLE, size) entries.
    // Returns zeta and fills zetaCandRho with the rhoVecs for the sampled trees.
    auto estimateZeta = [&](const NodeEntry& ta,
                            const vector<pTree>& allCovers,
                            vector<pair<int,vector<oRho>>>& zetaCandRho) -> float {
        // zetaCandRho: pairs of (index into allCovers, rhoVec)
        Rcpp::Timer::ScopedTimer st(timer, "Upwards_estimateZeta");
        int nsamp = min(ZETA_SAMPLE, static_cast<int>(allCovers.size()));
        zetaCandRho.clear();
        zetaCandRho.reserve(nsamp);
        float zetaSum = 0.0f;
        //l ogFile << "    The nsamp = " << nsamp << "\n" << std::flush;
        for (int i = 0; i < nsamp; i++) {
            //l ogFile << "      For i = " << i << " the construction of rhoVec \n"<< std::flush;
            //auto t0 = chrono::high_resolution_clock::now();
            vector<oRho> rv = buildRhoVec(allCovers[i], ta.tree, ta.rhoVec);
            //auto t1 = chrono::high_resolution_clock::now();
            //double totalBuildRho = chrono::duration<double, milli>(t1 - t0).count();
            //l ogFile << "      it took " << totalBuildRho << " ms \n\n" << std::flush;
            float cnt = 0;
            //auto t2 = chrono::high_resolution_clock::now();
            for (int j = 0; j < N; j++)
                if (rv[j].rho - ta.rhoVec[j].rho > 0)
                    cnt += nSample[j];
            zetaSum += cnt / B;
            //auto t3 = chrono::high_resolution_clock::now();
            //double innerLoopT = chrono::duration<double, milli>(t3 - t2).count();
            zetaCandRho.push_back({i, rv});
            //l ogFile << "      And the final inner loop " << innerLoopT << " ms \n" << std::flush;
        }
        return (nsamp > 0) ? zetaSum / nsamp : 0.0f;
    };

    auto computeNu = [&](int xIdx, int yIdx) -> int {
        Rcpp::Timer::ScopedTimer st(timer, "Upwards_computeNu");
        vector<bool> desc = computeDesc(*this, yIdx);
        vector<bool> anc  = computeAnc (*this, xIdx);
        return maxLevelInIe(*this, desc, anc, Poset[yIdx].Tree.rank);
    };

    auto computeNuEmpty = [&](int yIdx) -> int {
        Rcpp::Timer::ScopedTimer st(timer, "Upwards_computeNuEmpty");
        vector<bool> desc = computeDesc(*this, yIdx);
        vector<bool> anc(static_cast<int>(Poset.size()), false);
        return maxLevelInIe(*this, desc, anc, Poset[yIdx].Tree.rank);
    };

    // Threshold functions use curMaxRnk (current subposet max rank)
    auto thresholdFn = [&](float zeta, int nu, int rankTb) -> float {
        Rcpp::Timer::ScopedTimer st(timer, "Upwards_thresholdFn");
        float base  = min(0.5f, zeta);
        float inner = static_cast<float>(nu * max((curMaxRnk - rankTb + 1),1))
                      / static_cast<float>(q * max(curMaxRnk,1));
        float logVal = logf(inner);
        if (logVal <= 0.0f) return base;     // inner < 1 → log negative → no bonus
        return base + sqrtf(logVal / (2.0f * static_cast<float>(B)));
    };

    auto laxThresholdFn = [&](float zeta, int rankTb) -> float {
        Rcpp::Timer::ScopedTimer st(timer, "Upwards_laxThresholdFn");
        float base  = min(0.5f, zeta);
        float inner = static_cast<float>(max((curMaxRnk - rankTb + 1),1))
                      / static_cast<float>(q * max(curMaxRnk,1));
        float logVal = logf(inner);
        if (logVal <= 0.0f) return base;     // inner < 1 → log negative → no bonus
        return base + sqrtf(logVal / (2.0f * static_cast<float>(B)));
        
    };

    // ---------------------------------------------------------------
    // insertNode
    // ---------------------------------------------------------------
    auto insertNode = [&](const pTree& T, const vector<oRho>& rhoT,
                          float zetaVal) -> int {
        Rcpp::Timer::ScopedTimer st(timer, "Upwards_insertNode");
        int rnk = T.rank;
        int runIdx = firstRank[rnk];
        while (runIdx > -1) {
            if (Poset[runIdx].Tree == T) return runIdx;
            runIdx = Poset[runIdx].next;
        }
        int newIdx = static_cast<int>(Poset.size());
        spNode spT = spNode(T);
        Poset.push_back(spT);
        rhoCache.push_back(rhoT);
        zetaCache.push_back(zetaVal);
        

        if (firstRank[rnk] == -1) firstRank[rnk] = newIdx;
        if (lastRank[rnk]  > -1)  Poset[lastRank[rnk]].setNext(newIdx);
        lastRank[rnk] = newIdx;

        if ((rnk + 1 <= maxRnkAbs) && (firstRank[rnk + 1] > -1)) {
            int scanIdx = firstRank[rnk + 1];
            while (scanIdx > -1) {
                if (Poset[scanIdx].Tree.covers(T)) {
                    Poset[newIdx].addParent(scanIdx);
                    Poset[scanIdx].addChild(newIdx);
                }
                scanIdx = Poset[scanIdx].next;
            }
        }
        if ((rnk - 1 >= 1) && (firstRank[rnk - 1] > -1)) {
            int scanIdx = firstRank[rnk - 1];
            while (scanIdx > -1) {
                if (T.covers(Poset[scanIdx].Tree)) {
                    Poset[newIdx].addChild(scanIdx);
                    Poset[scanIdx].addParent(newIdx);
                }
                scanIdx = Poset[scanIdx].next;
            }
        }
        if (rnk == maxRnkAbs) Msize++;
        if (rnk > curMaxRnk) curMaxRnk = rnk;
        return newIdx;
    };

    // ---------------------------------------------------------------
    // removeLastNode
    // ---------------------------------------------------------------
    auto removeLastNode = [&](int idx) {
        Rcpp::Timer::ScopedTimer st(timer, "Upwards_removeLastNode");
        int rnk = Poset[idx].Tree.rank;

        for (int p : Poset[idx].over) {
            auto& ch = Poset[p].under;
            ch.erase(remove(ch.begin(), ch.end(), idx), ch.end());
        }
        for (int c : Poset[idx].under) {
            auto& pr = Poset[c].over;
            pr.erase(remove(pr.begin(), pr.end(), idx), pr.end());
        }
        if (firstRank[rnk] == idx) {
            firstRank[rnk] = -1;
            lastRank[rnk]  = -1;
        } else {
            int prev = firstRank[rnk];
            while (prev > -1 && Poset[prev].next != idx)
                prev = Poset[prev].next;
            if (prev > -1) {
                Poset[prev].setNext(-1);
                lastRank[rnk] = prev;
            }
        }
        if (rnk == maxRnkAbs) Msize--;
        // Recompute curMaxRnk if needed
        if (rnk == curMaxRnk) {
            curMaxRnk = 0;
            for (int r = rnk - 1; r >= 1; r--) {
                if (firstRank[r] > -1) { curMaxRnk = r; break; }
            }
        }
        Poset.pop_back();
        rhoCache.pop_back();
        zetaCache.pop_back();
    };

    // ---------------------------------------------------------------
    // caveatBlocked: true if adding Tb would create an unreachable gap.
    // For each node h at rank >= Tb.rank+2 that is above Tb, check that
    // a downward path through the existing poset connects h to a node
    // at Tb.rank+1 that covers Tb.
    // ---------------------------------------------------------------
    auto caveatBlocked = [&](const pTree& Tb) -> bool {
        Rcpp::Timer::ScopedTimer st(timer, "Upwards_caveatBlocked");
        int tbRnk = Tb.rank;
        if (tbRnk + 2 > curMaxRnk) return false;

        vector<int> higherNodes;
        for (int r = tbRnk + 2; r <= curMaxRnk; r++) {
            int scanIdx = firstRank[r];
            while (scanIdx > -1) {
                if (Poset[scanIdx].Tree.over(Tb))
                    higherNodes.push_back(scanIdx);
                scanIdx = Poset[scanIdx].next;
            }
        }
        if (higherNodes.empty()) return false;

        int n = static_cast<int>(Poset.size());
        for (int h : higherNodes) {
            bool connected = false;
            queue<int> bfsQ;
            vector<bool> visited(n, false);
            bfsQ.push(h);
            visited[h] = true;
            while (!bfsQ.empty() && !connected) {
                int v = bfsQ.front(); bfsQ.pop();
                if (Poset[v].Tree.rank == tbRnk + 1) {
                    if (Poset[v].Tree.covers(Tb)) connected = true;
                    continue;
                }
                for (int c : Poset[v].under) {
                    if (!visited[c] && Poset[c].Tree.rank > tbRnk) {
                        visited[c] = true;
                        bfsQ.push(c);
                    }
                }
            }
            if (!connected) return true;
        }
        return false;
    };

    // ---------------------------------------------------------------
    // hasFullThresholdChain: forward DP using cached rho and zeta.
    // ---------------------------------------------------------------
    auto hasFullThresholdChain = [&]() -> bool {
        Rcpp::Timer::ScopedTimer st(timer, "Upwards_hasFullThresholdChain");
        int n = static_cast<int>(Poset.size());
        vector<bool> strictReach(n, false);

        int scanIdx = firstRank[1];
        while (scanIdx > -1) {
            strictReach[scanIdx] = true;
            scanIdx = Poset[scanIdx].next;
        }

        for (int r = 1; r < curMaxRnk; r++) {
            int u = firstRank[r];
            while (u > -1) {
                if (strictReach[u]) {
                    float zeta = zetaCache[u];
                    for (int v : Poset[u].over) {
                        float score = computeScore(rhoCache[u], rhoCache[v]);
                        int   nu    = computeNu(u, v);
                        float thr   = thresholdFn(zeta, nu, Poset[v].Tree.rank);
                        if (score >= thr)
                            strictReach[v] = true;
                    }
                }
                u = Poset[u].next;
            }
        }

        int topIdx = firstRank[curMaxRnk];
        while (topIdx > -1) {
            if (strictReach[topIdx]) return true;
            topIdx = Poset[topIdx].next;
        }
        return false;
    };

    // ---------------------------------------------------------------
    // tryGrowUp: search through allCovers of ta one at a time (in the
    // shuffled order) until one passes the strict threshold.
    // Zeta is estimated from the first ZETA_SAMPLE covers (already
    // shuffled). We reuse the rhoVecs computed during zeta estimation
    // for those candidates, and compute on-demand for the rest.
    // Returns true and fills tbEntry on success.
    // ---------------------------------------------------------------
    auto tryGrowUp = [&](const NodeEntry& ta,
                         vector<pTree>& allCovers, // pre-shuffled
                         float zeta,
                         const vector<pair<int,vector<oRho>>>& zetaCandRho,
                         // zetaCandRho[i] = {index into allCovers, rhoVec}
                         // for the first ZETA_SAMPLE entries
                         NodeEntry& tbEntry) -> bool {
        Rcpp::Timer::ScopedTimer st(timer, "Upwards_tryGrowUp");
        // Build a lookup from cover index -> rhoVec for the zeta sample
        // so we don't recompute rho for those candidates
        unordered_map<int,int> zetaRhoIdx; // allCovers index -> zetaCandRho index
        zetaRhoIdx.reserve(zetaCandRho.size());
        for (int i = 0; i < static_cast<int>(zetaCandRho.size()); i++)
            zetaRhoIdx[zetaCandRho[i].first] = i;
        
        pTree taCopy = ta.tree;

        for (int ci = 0; ci < static_cast<int>(allCovers.size()); ci++) {
            const pTree& Tb = allCovers[ci];
            pTree tbCopy = Tb;
            if (caveatBlocked(Tb)){
                //if (taCopy.rank == 0){
                    //l ogFile << "For ta = empty and tb = " << mPhylo(tbCopy).toNewick() << " the tree was blocked by caveat \n" << flush;
                //} else {
                    //l ogFile << "For ta = " << mPhylo(taCopy).toNewick() << " and tb = " << mPhylo(tbCopy).toNewick() << " the tree was blocked by caveat \n" << flush;
                //}
                
                continue;}

            // Get or compute rhoVec for Tb
            vector<oRho> rhoTb;
            auto it = zetaRhoIdx.find(ci);
            if (it != zetaRhoIdx.end()) {
                rhoTb = zetaCandRho[it->second].second; // reuse from zeta sample
            } else {
                rhoTb = buildRhoVec(Tb, ta.tree, ta.rhoVec);
            }

            float score = computeScore(ta.rhoVec, rhoTb);

            int  tbIdx  = insertNode(Tb, rhoTb, 0.0f); // placeholder zeta
            bool wasNew = (tbIdx == static_cast<int>(Poset.size()) - 1);
            int   nu  = (ta.posetIdx >= 0)
                        ? computeNu(ta.posetIdx, tbIdx)
                        : computeNuEmpty(tbIdx);
            float thr = thresholdFn(zeta, nu, Tb.rank);
            
            //if (taCopy.rank == 0){
                //l ogFile << "For ta = empty and tb = "<< mPhylo(tbCopy).toNewick() << " the score was" << score << " and the threshold was "<< thr << "\n" << flush;
            //}else {
                //l ogFile << "For ta = " << mPhylo(taCopy).toNewick() << " and tb = "<< mPhylo(tbCopy).toNewick() << " the score was" << score << " and the threshold was "<< thr << "\n" << flush;
            //}
            
            //l ogFile << "   zeta = " << zeta << "\n" << flush;
            //l ogFile << "   nu = " << nu << "\n" << flush;
            //l ogFile << "   Tb.rank = " << Tb.rank << "\n" << flush;
            //l ogFile << "   curMaxRnk = " << curMaxRnk << "\n" << flush;
            //l ogFile << "   B = " << B << "\n" << flush;

            if (score >= thr) {
                // Accepted — estimate zeta for Tb and cache it
                //l ogFile << "Computing cover trees for tree with rank " << Tb.rank << "\n" << flush;
                vector<pTree> tbCovers = coverTrees(Tb, compLeafSet);
                //l ogFile << "It finished computing the cover trees, with a total of " << tbCovers.size() << "\n" << flush;
                
                shuffle(tbCovers.begin(), tbCovers.end(), rng);
                vector<pair<int,vector<oRho>>> tbZetaCandRho;
                //l ogFile << "Estimating zeta \n" << flush;
                float zetaTb = estimateZeta({Tb, rhoTb, tbIdx, 0.0f},
                                            tbCovers, tbZetaCandRho);
                //l ogFile << "Finished estimating zeta \n" << flush;
                zetaCache[tbIdx] = zetaTb;
                tbEntry = {Tb, rhoTb, tbIdx, zetaTb};
                return true;
            } else {
                if (wasNew) removeLastNode(tbIdx);
            }
        }
        return false;
    };

    // ---------------------------------------------------------------
    // PHASE 1: Grow initial spine from empty tree upward
    // ---------------------------------------------------------------
    
    timer.tic("Upwards_SpineFirstTree");
    vector<NodeEntry> chain;

    {
        pTree emptyTree = pTree("();");
        vector<set<string>> emptyLeaves;
        oRho   baseORho = oRho(0, emptyLeaves);
        vector<oRho> rhoEmpty(N, baseORho);
        
        //l ogFile << "Computing cover trees for tree with rank 0 \n" << flush;
        vector<pTree> rank1Trees = coverTrees(emptyTree, compLeafSet);
        //l ogFile << "It finished computing the cover trees, with a total of " << rank1Trees.size() << "\n" << flush;
        shuffle(rank1Trees.begin(), rank1Trees.end(), rng);

        // Estimate zeta for empty tree using sample of rank-1 trees
        vector<pair<int,vector<oRho>>> zetaCandRho;
        NodeEntry emptyEntry = {emptyTree, rhoEmpty, -1, 0.0f};
        //l ogFile << "Estimating zeta \n" << flush;
        float zeta0 = estimateZeta(emptyEntry, rank1Trees, zetaCandRho);
        //l ogFile << "Finished estimating zeta \n" << flush;
        

        // Set curMaxRnk to rank of first accepted tree (rank 1)
        // but we need something > 0 for thresholdFn — set temporarily
        curMaxRnk = 1;

        NodeEntry tbEntry;
        if (!tryGrowUp(emptyEntry, rank1Trees, zeta0, zetaCandRho, tbEntry)) {
            return;
        }
        chain.push_back(tbEntry);
    }
    timer.toc("Upwards_SpineFirstTree");

    // Grow spine upward one rank at a time
    timer.tic("Upwards_SpineBuilding");
    while (chain.back().tree.rank < maxRnkAbs) {
        NodeEntry& ta = chain.back();
        
        //l ogFile << "Computing cover trees for tree with rank " << ta.tree.rank << "\n" << flush;
        vector<pTree> covers = coverTrees(ta.tree, compLeafSet);
        //l ogFile << "It finished computing the cover trees, with a total of " << covers.size() << "\n" << flush;
        
        pTree taCopy  = ta.tree;
        if (covers.empty()){ 
            break;}
        shuffle(covers.begin(), covers.end(), rng);

        // Estimate zeta using sample (reuse rhoVecs for first ZETA_SAMPLE)
        vector<pair<int,vector<oRho>>> zetaCandRho;
        //l ogFile << "Estimating zeta \n" << flush;
        float zeta = estimateZeta(ta, covers, zetaCandRho);
        //l ogFile << "Finished estimating zeta \n" << flush;
        
        // Update cached zeta now that we have a better estimate
        zetaCache[ta.posetIdx] = zeta;
        ta.zeta = zeta;

        NodeEntry tbEntry;
        if (tryGrowUp(ta, covers, zeta, zetaCandRho, tbEntry)) {
            chain.push_back(tbEntry);
        } else {
            break;
        }
    }
    timer.toc("Upwards_SpineBuilding");
    //l ogFile << "Spine built. Length: " << chain.size() << "  curMaxRnk: " << curMaxRnk << "\n" << flush;

    // ---------------------------------------------------------------
    // PHASE 2: Backward bifurcation
    // ---------------------------------------------------------------
    
    //l ogFile << "BIFURCATION! \n" << flush; 
    timer.tic("Upwards_Phase2");
    for (int k = static_cast<int>(chain.size()) - 2; k >= 0; k--) {
        NodeEntry& ta = chain[k];
        const pTree* spineSucc = (k + 1 < static_cast<int>(chain.size()))
                                 ? &chain[k + 1].tree : nullptr;

        //l ogFile << "Bifurcating from spine node " << k << " (rank " << ta.tree.rank << ")\n" << flush;
        
        //l ogFile << "Computing cover trees for tree with rank " << ta.tree.rank << "\n" << flush;
        vector<pTree> covers = coverTrees(ta.tree, compLeafSet);
        //l ogFile << "It finished computing the cover trees, with a total of " << covers.size() << "\n" << flush;
        
        if (covers.empty()) continue;

        // Remove spine successor from the candidate list
        if (spineSucc) {
            covers.erase(remove_if(covers.begin(), covers.end(),
                                   [&](const pTree& T){ return T == *spineSucc; }),
                         covers.end());
        }
        if (covers.empty()) continue;
        shuffle(covers.begin(), covers.end(), rng);

        // Reuse cached zeta (already estimated during spine growth)
        float zeta = ta.zeta;

        // Estimate zeta over the reduced candidate set for bifurcation
        // (spine successor excluded). Small correction but avoids bias.
        vector<pair<int,vector<oRho>>> zetaCandRho;
        
        //l ogFile << "Estimating zeta \n" << flush;
        zeta = estimateZeta(ta, covers, zetaCandRho);
        //l ogFile << "Finished estimating zeta \n" << flush;
        

        // Search one at a time for a candidate passing threshold + full chain
        unordered_map<int,int> zetaRhoIdx;
        zetaRhoIdx.reserve(zetaCandRho.size());
        for (int i = 0; i < static_cast<int>(zetaCandRho.size()); i++)
            zetaRhoIdx[zetaCandRho[i].first] = i;

        for (int ci = 0; ci < static_cast<int>(covers.size()); ci++) {
            const pTree& Tb = covers[ci];
            if (caveatBlocked(Tb)) continue;

            vector<oRho> rhoTb;
            auto it = zetaRhoIdx.find(ci);
            if (it != zetaRhoIdx.end())
                rhoTb = zetaCandRho[it->second].second;
            else
                rhoTb = buildRhoVec(Tb, ta.tree, ta.rhoVec);

            float score = computeScore(ta.rhoVec, rhoTb);

            int  tbIdx  = insertNode(Tb, rhoTb, 0.0f);
            bool wasNew = (tbIdx == static_cast<int>(Poset.size()) - 1);

            int   nu  = computeNu(ta.posetIdx, tbIdx);
            float thr = thresholdFn(zeta, nu, Tb.rank);

            if (score >= thr && hasFullThresholdChain()) {
                // Accepted — estimate zeta for Tb
                //l ogFile << "Computing cover trees for tree with rank " << Tb.rank << "\n" << flush;
                vector<pTree> tbCovers = coverTrees(Tb, compLeafSet);
                //l ogFile << "It finished computing the cover trees, with a total of " << tbCovers.size() << "\n" << flush;
        
                
                shuffle(tbCovers.begin(), tbCovers.end(), rng);
                vector<pair<int,vector<oRho>>> tbZetaCandRho;
                //l ogFile << "Estimating zeta \n" << flush;
                float zetaTb = estimateZeta({Tb, rhoTb, tbIdx, 0.0f},
                                            tbCovers, tbZetaCandRho);
                //l ogFile << "Finished estimating zeta \n" << flush;
                
                zetaCache[tbIdx] = zetaTb;

                //l ogFile << "  Bifurcation accepted at rank " << Tb.rank << "\n" << flush;

                // Grow upward from Tb
                NodeEntry curEntry = {Tb, rhoTb, tbIdx, zetaTb};
                vector<pTree>* curCovers = &tbCovers;
                vector<pair<int,vector<oRho>>>* curZetaCandRho = &tbZetaCandRho;

                while (curEntry.tree.rank < curMaxRnk) {
                    NodeEntry nextEntry;
                    if (tryGrowUp(curEntry, *curCovers,
                                  curEntry.zeta, *curZetaCandRho, nextEntry)) {
                        curEntry = nextEntry;
                        // Rebuild covers for next iteration
                        //l ogFile << "Computing cover trees for tree with rank " << curEntry.tree.rank << "\n" << flush;
                        tbCovers = coverTrees(curEntry.tree, compLeafSet);
                        //l ogFile << "It finished computing the cover trees, with a total of " << tbCovers.size() << "\n" << flush; 
                        shuffle(tbCovers.begin(), tbCovers.end(), rng);
                        tbZetaCandRho.clear();
                        //l ogFile << "Estimating zeta \n" << flush;
                        float zetaNext = estimateZeta(curEntry, tbCovers, tbZetaCandRho);
                        //l ogFile << "Finished estimating zeta \n" << flush;
                        zetaCache[curEntry.posetIdx] = zetaNext;
                        curEntry.zeta = zetaNext;
                    } else {
                        break;
                    }
                }
                break; // move to next bifurcation point
            } else {
                if (wasNew) removeLastNode(tbIdx);
            }
        }
    }
    timer.toc("Upwards_Phase2");

    //l ogFile << "Bifurcation done. Poset size: " << Poset.size() << "\n" << flush;

    // ---------------------------------------------------------------
    // PHASE 3: Zeta assignment
    // ---------------------------------------------------------------
    timer.tic("Upwards_SavingZeta");
    {
        int scanIdx = firstRank[curMaxRnk];
        while (scanIdx > -1) {
            Poset[scanIdx].setZeta(zetaCache[scanIdx]);
            scanIdx = Poset[scanIdx].next;
        }
    }
    for (int r = curMaxRnk - 1; r > 1; r--) {
        int scanIdx = firstRank[r];
        while (scanIdx > -1) {
            Poset[scanIdx].setZeta(zetaCache[scanIdx]);
            scanIdx = Poset[scanIdx].next;
        }
    }
    if (firstRank[1] > -1) {
        int scanIdx = firstRank[1];
        while (scanIdx > -1) {
            Poset[scanIdx].setZeta(zetaCache[scanIdx]);
            scanIdx = Poset[scanIdx].next;
        }
    }
    timer.toc("Upwards_SavingZeta");
    //l ogFile << "Constructive builder done. Poset size: " << Poset.size() << "\n" << flush;
}


// ---------------------------------------------------------------------
// Fixed-width builder: bottom-up, with an explicit per-rank width derived
// from top_width / bottom_width and an orientation ("upwards"/"downwards").
//
//   upwards   : width plateaus at bottom_width near the bottom and ramps
//               up by +1 per rank to top_width at the maximal rank
//               (so the subposet is widest at the top); expects
//               top_width >= bottom_width.
//   downwards : width starts at bottom_width at rank 1 and ramps down by
//               -1 per rank to top_width, then plateaus at top_width
//               (so the subposet is widest at the bottom); expects
//               bottom_width >= top_width.
//
// Each rank is filled (bottom-up) with up to width(rank) trees: the
// best-scoring covers of the trees already admitted one rank below,
// deduped by Newick string, taking the best score when a tree covers
// several admitted parents. A candidate is admitted only if it honors the
// connectivity caveat: every already-admitted tree that sits below it must
// be reachable from it by descending one covering step at a time without
// leaving the subposet. If a rank runs out of caveat-respecting candidates
// before reaching its target width, it is left narrower and a warning that
// names the rank is emitted.
//
// The empty tree is the implicit bottom and is NOT stored as a node (as in
// the other builders); rank-1 nodes are the minimal nodes and carry no
// `under` edges. The empty tree's zeta is always 1/3, so it is never
// computed or stored here.
//
// Each node stores zeta = mean, over ALL of its covers, of
// score(node -> cover). Because every cover of a node is enumerated while
// the rank above it is filled, this is computed exactly there (no sampling);
// nodes at the maximal rank have no covers and keep zeta = 0.
// ---------------------------------------------------------------------
subPoset::subPoset(vector<pTree> Sample, vector<int> nSample,
                   set<string> compLeafSet,
                   int top_width, int bottom_width, string orientation) {

    Rcpp::Timer timer("times_fixedwidth");
    int N      = static_cast<int>(Sample.size());
    int B      = accumulate(nSample.begin(), nSample.end(), 0);
    int maxRnk = 2 * static_cast<int>(compLeafSet.size()) - 7;

    firstRank.assign(2 * static_cast<int>(compLeafSet.size()) - 6, -1);
    lastRank.assign(2 * static_cast<int>(compLeafSet.size()) - 6, -1);

    bool upwards = (orientation != "downwards");   // default to upwards
    if (orientation != "upwards" && orientation != "downwards")
        Rcpp::warning("subPoset(fixed width): unrecognized orientation '" +
                      orientation + "', defaulting to 'upwards'.");

    // Target width at a given rank (always >= 1).
    auto widthAt = [&](int r) -> int {
        int w = upwards ? max(bottom_width, top_width - (maxRnk - r))
                        : max(bottom_width - (r - 1), top_width);
        return max(1, w);
    };

    struct Candidate { pTree tree; vector<oRho> rhoVec; float score; };

    // Weighted fraction of sample trees whose rho strictly improves.
    auto computeScore = [&](const vector<oRho>& rhoBase,
                            const vector<oRho>& rhoV) -> float {
        Rcpp::Timer::ScopedTimer st(timer, "FixedW_computeScore");
        float sum = 0;
        for (int i = 0; i < N; i++)
            if (rhoV[i].rho - rhoBase[i].rho > 0)
                sum += nSample[i];
        return sum / B;
    };

    // Incremental rho of cover T (one step above TB) for every sample tree.
    auto buildRhoVec = [&](const pTree& T, const pTree& TB,
                           const vector<oRho>& rhoBase) -> vector<oRho> {
        Rcpp::Timer::ScopedTimer st(timer, "FixedW_buildRhoVec");
        vector<oRho> rv(N);
        for (int i = 0; i < N; i++)
            rv[i] = rho(T, TB, rhoBase[i], Sample[i]);
        return rv;
    };

    auto toNwk = [&](const pTree& T) -> string {
        Rcpp::Timer::ScopedTimer st(timer, "FixedW_toNwk");
        mPhylo mp = mPhylo(T);
        return mp.toNewick();
    };

    // rhoCache mirrors Poset indices.
    vector<vector<oRho>> rhoCache;

    // Connectivity caveat (used for rank >= 2). A candidate Y (rank r) may be
    // admitted only if every already-admitted node X with X < Y is reachable
    // from Y's in-poset lower covers by walking down `under` edges. lowerCovers
    // are the rank r-1 nodes that Y covers. The implicit empty bottom sits
    // below everything, so rank-1 candidates need no check and skip this.
    auto caveatOk = [&](const pTree& Y, const vector<int>& lowerCovers,
                        int r) -> bool {
        Rcpp::Timer::ScopedTimer st(timer, "FixedW_caveat");
        if (lowerCovers.empty()) return false;          // nothing to attach to
        vector<char> reach(Poset.size(), 0);
        queue<int> bfs;
        for (int w : lowerCovers)
            if (!reach[w]) { reach[w] = 1; bfs.push(w); }
        while (!bfs.empty()) {
            int u = bfs.front(); bfs.pop();
            for (int d : Poset[u].under)
                if (!reach[d]) { reach[d] = 1; bfs.push(d); }
        }
        for (int x = 0; x < (int)Poset.size(); x++) {
            if (reach[x]) continue;
            if (Poset[x].Tree.rank >= r) continue;      // only nodes below Y
            if (Y.over(Poset[x].Tree)) return false;    // X < Y but disconnected
        }
        return true;
    };

    // -----------------------------------------------------------------
    // Ranks 1..maxRnk: fill each level bottom-up. The empty tree is the
    // implicit bottom (not stored); rank-1 nodes seed from its covers.
    // -----------------------------------------------------------------
    pTree emptyTree = pTree("();");
    vector<oRho> rhoEmpty(N, oRho());

    for (int r = 1; r <= maxRnk; r++) {
        int wTarget = widthAt(r);

        // Gather candidate covers for this rank. For r == 1 they are the covers
        // of the (implicit) empty tree; for r >= 2 they are the covers of every
        // admitted node at rank r-1. While expanding a real parent we finalize
        // its zeta as the mean score over ALL of its covers (admitted or not).
        vector<Candidate>          candList;
        unordered_map<string,int>  candIndex;

        auto collectFrom = [&](const pTree& parent, const vector<oRho>& rhoParent,
                               float& zSum, int& zCnt) {
            vector<pTree> above = coverTrees(parent, compLeafSet);
            for (const pTree& Y : above) {
                vector<oRho> rv = buildRhoVec(Y, parent, rhoParent);
                float sc = computeScore(rhoParent, rv);
                zSum += sc;
                zCnt++;

                string nwk = toNwk(Y);
                auto it = candIndex.find(nwk);
                if (it == candIndex.end()) {
                    candIndex[nwk] = (int)candList.size();
                    candList.push_back({Y, rv, sc});
                } else if (sc > candList[it->second].score) {
                    candList[it->second].score = sc;    // keep best parent score
                }
            }
        };

        if (r == 1) {
            // Covers of the empty tree. Its zeta is the known constant 1/3 and
            // is not stored (the empty tree is not a node).
            float zSum = 0.0f; int zCnt = 0;
            collectFrom(emptyTree, rhoEmpty, zSum, zCnt);
        } else {
            int p = firstRank[r - 1];
            while (p > -1) {
                float zSum = 0.0f; int zCnt = 0;
                collectFrom(Poset[p].Tree, rhoCache[p], zSum, zCnt);
                Poset[p].setZeta(zCnt > 0 ? zSum / zCnt : 0.0f);
                p = Poset[p].next;
            }
        }

        // Best score first.
        sort(candList.begin(), candList.end(),
             [](const Candidate& a, const Candidate& b){
                 return a.score > b.score; });

        // Admit up to wTarget candidates that honor the caveat.
        int admitted = 0;
        for (const Candidate& cand : candList) {
            if (admitted >= wTarget) break;

            // In-poset lower covers: rank r-1 nodes covered by this tree.
            // Rank-1 nodes attach to the implicit empty bottom — no lower
            // covers, and the caveat is trivially satisfied.
            vector<int> lowerCovers;
            if (r >= 2) {
                int q = firstRank[r - 1];
                while (q > -1) {
                    if (cand.tree.covers(Poset[q].Tree))
                        lowerCovers.push_back(q);
                    q = Poset[q].next;
                }
                if (!caveatOk(cand.tree, lowerCovers, r)) continue;
            }

            // Insert and wire cover edges to every rank r-1 node it covers.
            Poset.push_back(spNode(cand.tree));
            int idx = (int)Poset.size() - 1;
            rhoCache.push_back(cand.rhoVec);

            if (firstRank[r] == -1) firstRank[r] = idx;
            if (lastRank[r]  > -1)  Poset[lastRank[r]].setNext(idx);
            lastRank[r] = idx;

            for (int w : lowerCovers) {
                Poset[idx].addChild(w);     // w sits below idx
                Poset[w].addParent(idx);    // idx sits above w
            }
            admitted++;
        }

        if (admitted < wTarget)
            Rcpp::warning("subPoset(fixed width): rank " + std::to_string(r) +
                          " filled " + std::to_string(admitted) + " of " +
                          std::to_string(wTarget) +
                          " trees (caveat-respecting candidates exhausted).");
    }
}


void subPoset::print(){
    for (int j = 1; j < firstRank.size(); j++){
        int curIndx = firstRank.at(j);
        while (curIndx > -1){
            Poset.at(curIndx).print();
            curIndx = Poset.at(curIndx).next;
        }
    }

    cout <<"Additional information: \n"<< std::flush;
    cout <<"    First in ranks: ";
    for (int j : firstRank){
        cout << to_string(j) << " ";
    }

    cout <<"\n    Last in ranks: ";
    for (int j : lastRank){
        cout << to_string(j) << " ";
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

    cout <<"Additional information: \n"<< std::flush;
    cout <<"    First in ranks: ";
    for (int j : firstRank){
        cout << to_string(j) << " ";
    }

    cout <<"\n    Last in ranks: ";
    for (int j : lastRank){
        cout << to_string(j) << " ";
    }
    cout << "\n"<< std::flush;
}

std::vector<bool> computeDesc(const subPoset& SP, int y) {
    int n = SP.Poset.size();
    std::vector<bool> desc(n, false);
    std::queue<int> q;
    q.push(y);
    desc[y] = true;
    while (!q.empty()) {
        int v = q.front(); q.pop();
        for (int w : SP.Poset[v].over) {
            if (!desc[w]) {
                desc[w] = true;
                q.push(w);
            }
        }
    }
    return desc;
}

// BFS upward from `start` following "under" links
// Returns a boolean mask of size nodes.size()
std::vector<bool> computeAnc(const subPoset& SP, int x) {
    int n = SP.Poset.size();
    std::vector<bool> anc(n, false);
    std::queue<int> q;
    q.push(x);
    anc[x] = true;
    while (!q.empty()) {
        int u = q.front(); q.pop();
        for (int w : SP.Poset[u].under) {
            if (!anc[w]) {
                anc[w] = true;
                q.push(w);
            }
        }
    }
    return anc;
}

int64_t countMaximalChains(const subPoset& SP, const std::vector<bool>& desc,  // Desc(y)
                              const std::vector<bool>& anc)   // Anc(x)
{
    int n = SP.Poset.size();
    std::vector<int64_t> paths(n, 0);

    // Rank-1 nodes: one path from implicit bottom, unless they are in Anc(x)
    // (meaning the edge bottom->v is removed because v in Anc(x))
    int v = SP.firstRank[1];
    while (v != -1) {
        if (!anc[v]) {
            paths[v] = 1;
        }
        v = SP.Poset[v].next;
    }

    int64_t total = 0;

    for (int r = 1; r < (int)SP.firstRank.size(); ++r) {
        int u = SP.firstRank[r];
        while (u != -1) {
                // Check if u has any surviving outgoing edges
                bool hasSurvivingEdge = false;
                if (!desc[u]) {
                    for (int w : SP.Poset[u].over) {
                        if (!anc[w]) {
                            hasSurvivingEdge = true;
                            break;
                        }
                    }
                }

                if (!hasSurvivingEdge) {
                    // u is a sink in I(e): chain ends here
                    total += paths[u];
                } else {
                    // Propagate to surviving over-neighbors
                    for (int w : SP.Poset[u].over) {
                        if (!anc[w]) {
                            if (paths[u] == 0){
                                paths[w]++;
                            } else {
                                paths[w] += paths[u];}
                        }
                    }
                }
            u = SP.Poset[u].next;
        }
    }
    
    return total;
}

// For a fixed e = (x, y), count surviving covering pairs at each rank level.
// A pair (x', y') survives in I(e) iff:
//   x' not in Desc(y)  →  !desc[x']
//   y' not in Anc(x)   →  !anc[y']
// Returns the size of the largest level in I(e).
int maxLevelInIe(const subPoset& SP,
                 const std::vector<bool>& desc,  // Desc(y)
                 const std::vector<bool>& anc,
                 const int& curRank)   // Anc(x)
{
    int maxCount = 0;

    for (int r = 1; r < (int)SP.firstRank.size(); ++r) {
        int count = 0;
        int v = SP.firstRank[r];
        while (v != -1) {
            if (!anc[v]) {
                if (r == 1){
                    count++;
                }
                // x' = v survives; now check each covering pair (v, w)
                for (int u : SP.Poset[v].under) {
                    if (!desc[u]) {
                        ++count;
                    }
                }
            }
            v = SP.Poset[v].next;
        }
        if (r != curRank){
            count++;
        }
        maxCount = std::max(maxCount, count);
    }

    return maxCount;
}

void computeAllChainCounts(subPoset& SP) {
    int n = SP.Poset.size();
    for (int y = 0; y < n; ++y) {
        std::vector<bool> desc = computeDesc(SP, y);
        
        if (SP.Poset[y].Tree.rank == 1) {
            // Implicit pair (bottom, y): Anc(bottom) = empty
            std::vector<bool> anc(n, false);
            SP.Poset[y].chainCountIe.resize(1, countMaximalChains(SP, desc, anc));
        } else {
            SP.Poset[y].chainCountIe.resize(SP.Poset[y].under.size());
            for (int i = 0; i < (int)SP.Poset[y].under.size(); ++i) {
                int x = SP.Poset[y].under[i];
                std::vector<bool> anc = computeAnc(SP, x);
                SP.Poset[y].chainCountIe[i] = countMaximalChains(SP, desc, anc);
            }
        }
    }
}

void computeAllMaxLevelBounds(subPoset& SP) {
    int n = SP.Poset.size();

    for (int y = 0; y < n; ++y) {
        // Compute Desc(y) once for all pairs (x, y)
        std::vector<bool> desc = computeDesc(SP, y);

        if (SP.Poset[y].Tree.rank == 1) {
            // Implicit pair (bottom, y): Anc(bottom) = empty
            std::vector<bool> anc(n, false);
            SP.Poset[y].boundAntichain.resize(1, maxLevelInIe(SP, desc, anc, 1));
        } else {
            SP.Poset[y].boundAntichain.resize(SP.Poset[y].under.size());
            for (int i = 0; i < (int)SP.Poset[y].under.size(); ++i) {
                int x = SP.Poset[y].under[i];
                std::vector<bool> anc = computeAnc(SP, x);
                SP.Poset[y].boundAntichain[i] = maxLevelInIe(SP, desc, anc, SP.Poset[y].Tree.rank);
            }
        }
    }
}
