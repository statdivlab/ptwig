#include "pTree.h"
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

    int initRank = initT.at(0).returnRank();


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

        Msize = 0;
    }

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
        float MaxValue  = 0;
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
            MaxValue = MaxValueL;
            IndexMax = IndexMaxL;
        } else {
            MaxValue = MaxValueS;
            IndexMax = IndexMaxS;
        }

        pTree V = AllV.at(IndexMax);

        Poset.push_back(spNode(V));

        if (V.returnRank() == rmax){
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

    //Assigning kappa to upper trees.
    curIndx = firstRank.at(rmax);
    while(curIndx > -1){
        Poset.at(curIndx).setKappa(Msize*(static_cast<int>(Poset.at(curIndx).under.size())));
        curIndx = Poset.at(curIndx).next;
    }

    for (int j = rmax - 1; j > initRank; j--){
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

        if (V.returnRank() > rb){
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

                if (U.returnRank() < V.returnRank() - 1){
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

                if (U.returnRank() < V.returnRank() - 1){
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

                if (U.returnRank() < V.returnRank() - 1){
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

                if (U.returnRank() < V.returnRank() - 1){
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
            curIndx = firstRank.at(curRank);
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

            curRank--;
            curIndx = firstRank.at(curRank);
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
    firstRank.at(initT.at(0).rank) = 0;
    lastRank.at(initT.at(0).rank) = 0;

    int initRank = initT.at(0).returnRank();


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

        Msize = 0;
    }

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
        float MaxValue  = 0;
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
            MaxValue = MaxValueL;
            IndexMax = IndexMaxL;
        } else {
            MaxValue = MaxValueS;
            IndexMax = IndexMaxS;
        }

        pTree V = AllV.at(IndexMax);

        Poset.push_back(spNode(V));

        if (V.returnRank() == rmax){
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

    //Assigning kappa to upper trees.
    curIndx = firstRank.at(rmax);
    while(curIndx > -1){
        Poset.at(curIndx).setKappa(Msize*(static_cast<int>(Poset.at(curIndx).under.size())));
        curIndx = Poset.at(curIndx).next;
    }

    for (int j = rmax - 1; j > initRank; j--){
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

        if (V.returnRank() > rb){
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

                if (U.returnRank() < V.returnRank() - 1){
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

                if (U.returnRank() < V.returnRank() - 1){
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

                if (U.returnRank() < V.returnRank() - 1){
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

                if (U.returnRank() < V.returnRank() - 1){
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
            curIndx = firstRank.at(curRank);
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

            curRank--;
            curIndx = firstRank.at(curRank);
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