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

void spNode::print(){
    cout << "This node with kappa " + to_string(kappa) + " has tree \n";
    Tree.print();
    cout << "\n parents: [";
    for (int i : over){
        cout << to_string(i) << ", ";
    }
    cout << "] \n";
    cout << "children: [";
    for (int i : under){
        cout << to_string(i) << ", ";
    }
    cout << "] \n";
    cout << "Next: "<< to_string(next) << " \n";

}

void spNode::printRd(){
    cout << "This node with kappa " + to_string(kappa) + " \n";
    cout << "\n parents: [";
    for (int i : over){
        cout << to_string(i) << ", ";
    }
    cout << "] \n";
    cout << "children: [";
    for (int i : under){
        cout << to_string(i) << ", ";
    }
    cout << "] \n";
    cout << "Next: "<< to_string(next) << " \n";
}
    
subPoset::subPoset(vector<pTree> initT, vector<pTree> Sample, set<string> compLeafSet, int rb){
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
    int B = static_cast<int>(Sample.size());
    int curIndx = 0;


    if (initRank < rmax){
        firstRank.at(curRank + 1) = -1;
        lastRank.at(curRank + 1) = -1;
    }
    
    cout << "The rmax = " << rmax << "\n";

    while(curRank < rmax){

        pTree U = Poset.at(curIndx).Tree;
        vector<pTree> AllV = coverTrees(U, compLeafSet);

        auto rd = std::random_device {};
        auto rng = std::default_random_engine { rd() };
        shuffle(begin(AllV), std::end(AllV), rng);

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
            for (pTree Z : Sample){
                if ((rho(V, Z) - rho(U, Z)) > 0){
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
    while (curRank > 1) {
        pTree V = Poset.at(curIndx).Tree;

        int toAdd = 0;

        if (V.rank > rb){
            toAdd = 1 - static_cast<int> (Poset.at(curIndx).under.size());
        } else {
            toAdd = 2 - static_cast<int> (Poset.at(curIndx).under.size());
        }

        if (toAdd == 1){
            float min1 = 1.2;
            bool addU1 = false;
            pTree U1;

            for (string a : V.leafSet){
                pTree U = V.Remove(a);

                if (U.rank < V.rank - 1){
                    continue;
                }

                bool continueOuter = false;
                for (int k1 : Poset.at(curIndx).under){
                    if (U == Poset.at(k1).Tree){
                        continueOuter = true;
                        break;
                    }
                }
                if (continueOuter){
                    continue;
                }

                float Sum = 0;

                for (pTree T : Sample){
                    if ((rho(V, T) > rho(U, T)) > 0){
                        Sum++;
                    }
                }
                float stb = Sum/B;
                if (stb < min1){
                    min1 = stb;
                    U1 = U;
                    addU1 = true;
                }
            }

            for (Split s : V.intSplits){
                pTree U = V.Remove(s);

                if (U.rank < V.rank - 1){
                    continue;
                }

                bool continueOuter = false;
                for (int k1 : Poset.at(curIndx).under){
                    if (U == Poset.at(k1).Tree){
                        continueOuter = true;
                        break;
                    }
                }
                if (continueOuter){
                    continue;
                }

                float Sum = 0;

                for (pTree T : Sample){
                    if ((rho(V, T) > rho(U, T)) > 0){
                        Sum++;
                    }
                }
                float stb = Sum/B;
                if (stb < min1){
                    min1 = stb;
                    U1 = U;
                    addU1 = true;
                }
            }

            if(addU1){
                Poset.push_back(spNode(U1));

                int runIndx = firstRank.at(curRank);

                while (runIndx > -1){
                    if (Poset.at(runIndx).Tree.covers(U1)){
                        Poset.back().addParent(runIndx);
                        Poset.at(runIndx).addChild(static_cast<int> (Poset.size()) - 1);
                    }
                    runIndx = Poset.at(runIndx).next;
                }

                if (firstRank.at(curRank-1) == -1){
                    firstRank.at(curRank-1) = static_cast<int> (Poset.size()) - 1;
                }

                if (lastRank.at(curRank - 1) > -1){
                    Poset.at(lastRank.at(curRank-1)).next = static_cast<int> (Poset.size()) - 1;
                }

                lastRank.at(curRank-1) = static_cast<int> (Poset.size()) - 1;
            }

        } else if (toAdd == 2){
            float min1 = 1.2;
            float min2 = 1.2;

            bool addU1 = false;
            bool addU2 = false;

            pTree U1;
            pTree U2;

            for (string a : V.leafSet){
                pTree U = V.Remove(a);

                if (U.rank < V.rank - 1){
                    continue;
                }

                bool continueOuter = false;
                for (int k1 : Poset.at(curIndx).under){
                    if (U == Poset.at(k1).Tree){
                        continueOuter = true;
                        break;
                    }
                }
                if (continueOuter){
                    continue;
                }


                float Sum = 0;

                for (pTree T : Sample){
                    if ((rho(V, T) > rho(U, T)) > 0){
                        Sum++;
                    }
                }
                float stb = Sum/B;

                if (stb < min1){
                    min2 = min1;
                    U2 = U1;
                    min1 = stb;
                    U1 = U;
                    if (addU1){
                        addU2 = true;
                    }
                    addU1 = true;
                } else if (stb < min2){
                    min2 = stb;
                    U2 = U;
                    addU2 = true;
                }
            }

            for (Split s : V.intSplits){
                pTree U = V.Remove(s);

                if (U.rank < V.rank - 1){
                    continue;
                }

                bool continueOuter = false;
                for (int k1 : Poset.at(curIndx).under){
                    if (U == Poset.at(k1).Tree){
                        continueOuter = true;
                        break;
                    }
                }
                if (continueOuter){
                    continue;
                }

                float Sum = 0;

                for (pTree T : Sample){
                    if ((rho(V, T) > rho(U, T)) > 0){
                        Sum++;
                    }
                }
                float stb = Sum/B;
                if (stb < min1){
                    min2 = min1;
                    U2 = U1;
                    min1 = stb;
                    U1 = U;
                    if (addU1){
                        addU2 = true;
                    }
                    addU1 = true;
                } else if (stb < min2){
                    min2 = stb;
                    U2 = U;
                    addU2 = true;
                }
            }

            if (addU1){
                Poset.push_back(spNode(U1));

                int runIndx = firstRank.at(curRank);

                while (runIndx > -1){
                    if (Poset.at(runIndx).Tree.covers(U1)){
                        Poset.back().addParent(runIndx);
                        Poset.at(runIndx).addChild(static_cast<int> (Poset.size()) - 1);
                    }
                    runIndx = Poset.at(runIndx).next;
                }

                if (firstRank.at(curRank-1) == -1){
                    firstRank.at(curRank-1) = static_cast<int> (Poset.size()) - 1;
                }

                if (lastRank.at(curRank - 1) > -1){
                    Poset.at(lastRank.at(curRank-1)).next = static_cast<int> (Poset.size()) - 1;
                }

                lastRank.at(curRank-1) = static_cast<int> (Poset.size()) - 1;
            }

            if (addU2){
                Poset.push_back(spNode(U2));

                int runIndx = firstRank.at(curRank);

                while (runIndx > -1){
                    if (Poset.at(runIndx).Tree.covers(U2)){
                        Poset.back().addParent(runIndx);
                        Poset.at(runIndx).addChild(static_cast<int> (Poset.size()) - 1);
                    }
                    runIndx = Poset.at(runIndx).next;
                }

                Poset.at(lastRank.at(curRank-1)).next = static_cast<int> (Poset.size()) - 1;
                lastRank.at(curRank-1) = static_cast<int> (Poset.size()) - 1;
            }

        }

        curIndx = Poset.at(curIndx).next;

        if (curIndx == -1){
            curRank--;
            curIndx = firstRank.at(curRank);
        }
    }
    
    cout << "The Msize = " << Msize << "\n";
    //Assigning kappa to upper trees.
    curIndx = firstRank.at(rmax);
    while(curIndx > -1){
        cout << "For curIndx = " << curIndx << "\n";
        cout << "And because of this, we have " << static_cast<int>(Poset.at(curIndx).under.size()) << "\n\n";
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
            cout << "For curIndx = " << curIndx << "\n";
            cout << "Second the kappa is " << kap << "\n";
            cout << "Second under size " << static_cast<int>(Poset.at(curIndx).under.size()) << "\n\n";
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
        cout << "For curIndx = " << curIndx << "\n";
        cout << "Fourth the kappa is " << kap << "\n\n";
        Poset.at(curIndx).setKappa(kap);
        curIndx = Poset.at(curIndx).next;
    }
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
    
     cout << "The rmax = " << rmax << "\n";

    while(curRank < rmax){

        pTree U = Poset.at(curIndx).Tree;
        vector<pTree> AllV = coverTrees(U, compLeafSet);

        auto rd = std::random_device {};
        auto rng = std::default_random_engine { rd() };
        shuffle(begin(AllV), std::end(AllV), rng);

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
            int currentIndx = 0;
            for (pTree Z : Sample){
                if ((rho(V, Z) - rho(U, Z)) > 0){
                    tempSum += nSample.at(currentIndx);
                }
                currentIndx++;
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
    while (curRank > 1) {
        pTree V = Poset.at(curIndx).Tree;

        int toAdd = 0;

        if (V.rank > rb){
            toAdd = 1 - static_cast<int> (Poset.at(curIndx).under.size());
        } else {
            toAdd = 2 - static_cast<int> (Poset.at(curIndx).under.size());
        }

        if (toAdd == 1){
            float min1 = 1.2;
            bool addU1 = false;
            pTree U1;

            for (string a : V.leafSet){
                pTree U = V.Remove(a);

                if (U.rank < V.rank - 1){
                    continue;
                }

                bool continueOuter = false;
                for (int k1 : Poset.at(curIndx).under){
                    if (U == Poset.at(k1).Tree){
                        continueOuter = true;
                        break;
                    }
                }
                if (continueOuter){
                    continue;
                }

                float Sum = 0;
                int currentIndx = 0;
                for (pTree T : Sample){
                    if ((rho(V, T) > rho(U, T)) > 0){
                        Sum += nSample.at(currentIndx);
                    }
                    currentIndx++;
                }
                float stb = Sum/B;
                if (stb < min1){
                    min1 = stb;
                    U1 = U;
                    addU1 = true;
                }
            }

            for (Split s : V.intSplits){
                pTree U = V.Remove(s);

                if (U.rank < V.rank - 1){
                    continue;
                }

                bool continueOuter = false;
                for (int k1 : Poset.at(curIndx).under){
                    if (U == Poset.at(k1).Tree){
                        continueOuter = true;
                        break;
                    }
                }
                if (continueOuter){
                    continue;
                }

                float Sum = 0;
                int currentIndx = 0;
                
                for (pTree T : Sample){
                    if ((rho(V, T) > rho(U, T)) > 0){
                        Sum += nSample.at(currentIndx);
                    }
                    currentIndx++;
                }
                float stb = Sum/B;
                if (stb < min1){
                    min1 = stb;
                    U1 = U;
                    addU1 = true;
                }
            }

            if(addU1){
                Poset.push_back(spNode(U1));

                int runIndx = firstRank.at(curRank);

                while (runIndx > -1){
                    if (Poset.at(runIndx).Tree.covers(U1)){
                        Poset.back().addParent(runIndx);
                        Poset.at(runIndx).addChild(static_cast<int> (Poset.size()) - 1);
                    }
                    runIndx = Poset.at(runIndx).next;
                }

                if (firstRank.at(curRank-1) == -1){
                    firstRank.at(curRank-1) = static_cast<int> (Poset.size()) - 1;
                }

                if (lastRank.at(curRank - 1) > -1){
                    Poset.at(lastRank.at(curRank-1)).next = static_cast<int> (Poset.size()) - 1;
                }

                lastRank.at(curRank-1) = static_cast<int> (Poset.size()) - 1;
            }

        } else if (toAdd == 2){
            float min1 = 1.2;
            float min2 = 1.2;

            bool addU1 = false;
            bool addU2 = false;

            pTree U1;
            pTree U2;

            for (string a : V.leafSet){
                pTree U = V.Remove(a);

                if (U.rank < V.rank - 1){
                    continue;
                }

                bool continueOuter = false;
                for (int k1 : Poset.at(curIndx).under){
                    if (U == Poset.at(k1).Tree){
                        continueOuter = true;
                        break;
                    }
                }
                if (continueOuter){
                    continue;
                }


                float Sum = 0;
                int currentIndx = 0;
                for (pTree T : Sample){
                    if ((rho(V, T) > rho(U, T)) > 0){
                        Sum += nSample.at(currentIndx);
                    }
                    currentIndx++;
                }
                float stb = Sum/B;

                if (stb < min1){
                    min2 = min1;
                    U2 = U1;
                    min1 = stb;
                    U1 = U;
                    if (addU1){
                        addU2 = true;
                    }
                    addU1 = true;
                } else if (stb < min2){
                    min2 = stb;
                    U2 = U;
                    addU2 = true;
                }
            }

            for (Split s : V.intSplits){
                pTree U = V.Remove(s);

                if (U.rank < V.rank - 1){
                    continue;
                }

                bool continueOuter = false;
                for (int k1 : Poset.at(curIndx).under){
                    if (U == Poset.at(k1).Tree){
                        continueOuter = true;
                        break;
                    }
                }
                if (continueOuter){
                    continue;
                }

                float Sum = 0;
                int currentIndx = 0;
                
                for (pTree T : Sample){
                    if ((rho(V, T) > rho(U, T)) > 0){
                        Sum += nSample.at(currentIndx);
                    }
                    currentIndx++;
                }
                float stb = Sum/B;
                if (stb < min1){
                    min2 = min1;
                    U2 = U1;
                    min1 = stb;
                    U1 = U;
                    if (addU1){
                        addU2 = true;
                    }
                    addU1 = true;
                } else if (stb < min2){
                    min2 = stb;
                    U2 = U;
                    addU2 = true;
                }
            }

            if (addU1){
                Poset.push_back(spNode(U1));

                int runIndx = firstRank.at(curRank);

                while (runIndx > -1){
                    if (Poset.at(runIndx).Tree.covers(U1)){
                        Poset.back().addParent(runIndx);
                        Poset.at(runIndx).addChild(static_cast<int> (Poset.size()) - 1);
                    }
                    runIndx = Poset.at(runIndx).next;
                }

                if (firstRank.at(curRank-1) == -1){
                    firstRank.at(curRank-1) = static_cast<int> (Poset.size()) - 1;
                }

                if (lastRank.at(curRank - 1) > -1){
                    Poset.at(lastRank.at(curRank-1)).next = static_cast<int> (Poset.size()) - 1;
                }

                lastRank.at(curRank-1) = static_cast<int> (Poset.size()) - 1;
            }

            if (addU2){
                Poset.push_back(spNode(U2));

                int runIndx = firstRank.at(curRank);

                while (runIndx > -1){
                    if (Poset.at(runIndx).Tree.covers(U2)){
                        Poset.back().addParent(runIndx);
                        Poset.at(runIndx).addChild(static_cast<int> (Poset.size()) - 1);
                    }
                    runIndx = Poset.at(runIndx).next;
                }

                Poset.at(lastRank.at(curRank-1)).next = static_cast<int> (Poset.size()) - 1;
                lastRank.at(curRank-1) = static_cast<int> (Poset.size()) - 1;
            }

        }

        curIndx = Poset.at(curIndx).next;

        if (curIndx == -1){
            curRank--;
            curIndx = firstRank.at(curRank);
        }
    }
    
    cout << "The Msize = " << Msize << "\n";
    //Assigning kappa to upper trees.
    curIndx = firstRank.at(rmax);
    while(curIndx > -1){
        cout << "For curIndx = " << curIndx << "\n";
        cout << "And because of this, we have " << static_cast<int>(Poset.at(curIndx).under.size()) << "\n\n";
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
            cout << "For curIndx = " << curIndx << "\n";
            cout << "Second the kappa is " << kap << "\n";
            cout << "Second under size " << static_cast<int>(Poset.at(curIndx).under.size()) << "\n\n";
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
        cout << "For curIndx = " << curIndx << "\n";
        cout << "Fourth the kappa is " << kap << "\n\n";
        Poset.at(curIndx).setKappa(kap);
        curIndx = Poset.at(curIndx).next;
    }
}

subPoset::subPoset(vector<pTree> Sample, set<string> compLeafSet, int Mt, int rb) {
    
    cout << "It entered the correct builder \n";

    int B      = static_cast<int>(Sample.size());
    int maxRnk = 2 * static_cast<int>(compLeafSet.size()) - 7;

    firstRank.assign(2 * static_cast<int>(compLeafSet.size()) - 6, -1);
    lastRank.assign(2 * static_cast<int>(compLeafSet.size()) - 6, -1);
    
    cout<< "Reached 1 \n" ;

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

    cout<< "Reached 2 \n"; 
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

    cout<< "Reached 3 \n" ;
    int counterTemp = 0;
    // --- Beam search: ascend rank by rank until maxRnk ---
    while ((int)currentLevel[0].tree.rank < maxRnk) {
        
        counterTemp++;
        cout << "Cycle number" << counterTemp << "\n";

        set<string>    seenNewick;
        vector<BeamEntry> futureLevel;
        futureLevel.reserve(Mt);
        
        vector<vector<Candidate>> aboveCurrents;
        
        cout << "Before step A \n";
        // Step A: from each beam entry, compute aboves and order by scores
        for (const BeamEntry& entry : currentLevel) {
            vector<pTree> above = coverTrees(entry.tree, compLeafSet);
            
            vector<Candidate> aboveOne;
            
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
        
        cout << "Before step B \n";
        // Step B: Fill up with aboves, avoiding repetitions. We have a running indexes for each
        vector<int> tryIndx;
        cout << "Step B 1 \n";
        
        tryIndx.assign(Mt, 0);
        cout << "Step B 2 \n";
        
        int currentIndex = 0;
        
        cout << "BTW the size of currentLevel is " << (int)currentLevel.size() << "\n";
        
        int counterA = 0;
        while(((int)futureLevel.size() < Mt) && (accumulate(tryIndx.begin(), tryIndx.end(), 0) > -(int)currentLevel.size())){
            counterA++;
            cout << "Inside loop " << counterA << "\n";
            cout << "currentIndex = " << currentIndex << "\n";
            cout << "trying index = " << tryIndx.at(currentIndex) << "\n";
            
            if ((tryIndx.at(currentIndex) == -1) || (tryIndx.at(currentIndex) >= aboveCurrents[currentIndex].size())){
                tryIndx.at(currentIndex) = -1;
                continue;
            }
            Candidate canU = aboveCurrents[currentIndex][tryIndx.at(currentIndex)];
            
            cout << "Step B.1 \n";
            string nwk = toNwk(canU.tree);
            cout << "Step B.2 \n";
            if (!seenNewick.count(nwk)){
                cout << "Step B.3 \n";
                futureLevel.push_back({canU.tree, canU.rhoVec});
                cout << "Step B.4 \n";
                seenNewick.insert(nwk);
                cout << "Step B.5 \n";
            }
            cout << "Step B.6 \n";
            tryIndx.at(currentIndex)++;
            cout << "Step B.7 \n";
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
        
        cout << "Before potential problem \n";
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

    cout<< "Reached 4 \n" ;
    vector<vector<float>> rhoCache;
    rhoCache.reserve(Mt);

    for (int m = 0; m < (int)currentLevel.size(); m++) {
        Poset.push_back(spNode(currentLevel[m].tree));
        rhoCache.push_back(currentLevel[m].rhoVec);   // no recomputation

        if (m < (int)currentLevel.size() - 1)
            Poset[m].setNext(m + 1);
        // last node: next stays -1 from spNode constructor

        // Maximal nodes: kappa placeholder — corrected in Phase 3 kappa pass
        Poset[m].setKappa(Mt);
    }

    firstRank[maxRnk] = 0;
    lastRank[maxRnk] = (int)currentLevel.size() - 1;
    
    cout<< "Reached 5 \n" ;

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
    cout<< "Reached 6 \n" ;
    int curRank = maxRnk;
    int curIndx = firstRank[curRank];

    if (curRank > 1) {
        firstRank[curRank - 1] = -1;
        lastRank[curRank - 1] = -1;
    }
    cout<< "Reached 7 \n" ;
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
    cout<< "Reached 8 \n" ;
    // Kappa pass for rank 1
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
    
    cout << "It entered the correct builder \n";

    int N      = static_cast<int>(Sample.size());
    int B      = accumulate(nSample.begin(), nSample.end(), 0);
    int maxRnk = 2 * static_cast<int>(compLeafSet.size()) - 7;

    firstRank.assign(2 * static_cast<int>(compLeafSet.size()) - 6, -1);
    lastRank.assign(2 * static_cast<int>(compLeafSet.size()) - 6, -1);
    
    cout<< "Reached 1 \n" ;

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

    cout<< "Reached 2 \n"; 
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

    cout<< "Reached 3 \n" ;
    int counterTemp = 0;
    // --- Beam search: ascend rank by rank until maxRnk ---
    while ((int)currentLevel[0].tree.rank < maxRnk) {
        
        counterTemp++;
        cout << "Cycle number" << counterTemp << "\n";

        set<string>    seenNewick;
        vector<BeamEntry> futureLevel;
        futureLevel.reserve(Mt);
        
        vector<vector<Candidate>> aboveCurrents;
        
        cout << "Before step A \n";
        // Step A: from each beam entry, compute aboves and order by scores
        for (const BeamEntry& entry : currentLevel) {
            vector<pTree> above = coverTrees(entry.tree, compLeafSet);
            
            vector<Candidate> aboveOne;
            
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
        
        cout << "Before step B \n";
        // Step B: Fill up with aboves, avoiding repetitions. We have a running indexes for each
        vector<int> tryIndx;
        cout << "Step B 1 \n";
        
        tryIndx.assign(Mt, 0);
        cout << "Step B 2 \n";
        
        int currentIndex = 0;
        
        cout << "BTW the size of currentLevel is " << (int)currentLevel.size() << "\n";
        
        int counterA = 0;
        while(((int)futureLevel.size() < Mt) && (accumulate(tryIndx.begin(), tryIndx.end(), 0) > -(int)currentLevel.size())){
            counterA++;
            cout << "Inside loop " << counterA << "\n";
            cout << "currentIndex = " << currentIndex << "\n";
            cout << "trying index = " << tryIndx.at(currentIndex) << "\n";
            
            if ((tryIndx.at(currentIndex) == -1) || (tryIndx.at(currentIndex) >= aboveCurrents[currentIndex].size())){
                tryIndx.at(currentIndex) = -1;
                continue;
            }
            Candidate canU = aboveCurrents[currentIndex][tryIndx.at(currentIndex)];
            
            cout << "Step B.1 \n";
            string nwk = toNwk(canU.tree);
            cout << "Step B.2 \n";
            if (!seenNewick.count(nwk)){
                cout << "Step B.3 \n";
                futureLevel.push_back({canU.tree, canU.rhoVec});
                cout << "Step B.4 \n";
                seenNewick.insert(nwk);
                cout << "Step B.5 \n";
            }
            cout << "Step B.6 \n";
            tryIndx.at(currentIndex)++;
            cout << "Step B.7 \n";
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
        
        cout << "Before potential problem \n";
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

    cout<< "Reached 4 \n" ;
    vector<vector<float>> rhoCache;
    rhoCache.reserve(Mt);

    for (int m = 0; m < (int)currentLevel.size(); m++) {
        Poset.push_back(spNode(currentLevel[m].tree));
        rhoCache.push_back(currentLevel[m].rhoVec);   // no recomputation

        if (m < (int)currentLevel.size() - 1)
            Poset[m].setNext(m + 1);
        // last node: next stays -1 from spNode constructor

        // Maximal nodes: kappa placeholder — corrected in Phase 3 kappa pass
        Poset[m].setKappa(Mt);
    }

    firstRank[maxRnk] = 0;
    lastRank[maxRnk] = (int)currentLevel.size() - 1;
    
    cout<< "Reached 5 \n" ;

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
    cout<< "Reached 6 \n" ;
    int curRank = maxRnk;
    int curIndx = firstRank[curRank];

    if (curRank > 1) {
        firstRank[curRank - 1] = -1;
        lastRank[curRank - 1] = -1;
    }
    cout<< "Reached 7 \n" ;
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
    cout<< "Reached 8 \n" ;
    // Kappa pass for rank 1
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
            cout << "\n For node " << to_string(curIndx) << "\n";
            Poset.at(curIndx).print();
            curIndx = Poset.at(curIndx).next;
        }
    }

    cout<<"Additional information: \n";
    cout<<"    First in ranks: ";
    for (int j : firstRank){
        cout<< to_string(j) << " ";
    }

    cout<<"\n    Last in ranks: ";
    for (int j : lastRank){
        cout<< to_string(j) << " ";
    }
    cout << "\n";
}

void subPoset::printRd(){
    for (int j = 1; j < firstRank.size(); j++){
        int curIndx = firstRank.at(j);
        while (curIndx > -1){
            cout << "\n For node " << to_string(curIndx) << "\n";
            Poset.at(curIndx).printRd();
            curIndx = Poset.at(curIndx).next;
        }
    }

    cout<<"Additional information: \n";
    cout<<"    First in ranks: ";
    for (int j : firstRank){
        cout<< to_string(j) << " ";
    }

    cout<<"\n    Last in ranks: ";
    for (int j : lastRank){
        cout<< to_string(j) << " ";
    }
    cout << "\n";
}

