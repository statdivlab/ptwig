#include "pTree.h"
#include "mPhylo.h"
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


        
mPhylo::mPhylo(pTree T){
    Nnode = static_cast<int>(T.intSplits.size()) + 1;

    vector<set<int>> independentgroups;

    set<int> LastIndependentGroup;

    int count = 0;
    for(string l : T.leafSet){
        tipLabel.push_back(l);
        label2node.insert({l,count});
        LastIndependentGroup.insert(count);
        count++;
        node2edge.push_back(vector<int> {});
    }

    int nCur = static_cast<int>(tipLabel.size());

    for (Split s : T.intSplits){
        set<int> tSide;

        if (s.side1.size() <= s.side2.size()){
            for (string a : s.side1){
                tSide.insert(label2node.at(a));
            }
        } else {
            for (string a : s.side2){
                tSide.insert(label2node.at(a));
            }
        }

        independentgroups.push_back(tSide);
    }

    std::sort(independentgroups.begin(), independentgroups.end(), [](const set<int> & a, const set<int> & b){ return a.size() < b.size(); });


    independentgroups.push_back(LastIndependentGroup);

    int p = 0;

    while (p < independentgroups.size()){
        set<int> curGroup = independentgroups.at(p);

        node2edge.push_back(vector<int>{});

        for (int nnode : curGroup){
            node2edge.at(nnode).push_back(static_cast<int>(edge.size()));
            node2edge.at(nCur).push_back(static_cast<int>(edge.size()));

            edge.push_back(vector<int> ({nnode,nCur}));
        }

        for (int q = 0; q < independentgroups.size(); q++){
            if (std::includes(independentgroups.at(q).begin(), independentgroups.at(q).end(),
                              curGroup.begin(), curGroup.end())){
                for (int nnode : curGroup){
                    independentgroups.at(q).erase(nnode);
                }
                independentgroups.at(q).insert(nCur);
            }
        }

        std::sort(independentgroups.begin(), independentgroups.end(), [](const set<int> & a, const set<int> & b){ return a.size() < b.size(); });

        nCur++;
        p++;
    }
}

void mPhylo::print(){
    cout << "The tree has " << tipLabel.size() << " leaves and " << Nnode << " internal nodes \n Its edges are: \n";
    for (vector<int> ed : edge){
        cout << "(" << ed.at(0) << ", " << ed.at(1) << ")\n";
    }
    cout << "\n The labels of the leaves are: ";

    for (string st : tipLabel){
        cout << st << " ";
    }
    cout << "\n \n These are associated to the nodes by: \n";

    for (auto i : label2node){
        cout << i.first << " : " << i.second << endl;
    }

    cout<< "\n \n And finally, here are the associations between the nodes and the edges: \n";

    for (int p = 0; p < node2edge.size(); p++){
        cout << p << ": {";
        for (int en : node2edge.at(p)){
            cout<< en << ", ";
        }
        cout << "} \n";
    }

}

string mPhylo::toNewick(){
    vector<bool> nodeVisited(tipLabel.size() + Nnode, false);
    stack<int> curNode;
    stack<int> curIndx;
    stack<string> curNewick;

    int stNode = edge.at(node2edge.at(0).at(0)).at(0);
    if (stNode == 0){
        stNode = edge.at(node2edge.at(0).at(0)).at(1);
    }

    string Result = "";

    curNode.push(stNode);
    curIndx.push(0);
    curNewick.push("");

    nodeVisited.at(stNode) = true;

    while(!curNode.empty()){
        int cNode = curNode.top();
        int cIndx = curIndx.top();


        if (cIndx >= node2edge.at(cNode).size()){
            string cString = curNewick.top();

            curNode.pop();
            curIndx.pop();
            curNewick.pop();

            if (!curNewick.empty()){
                string tempNewick = curNewick.top();
                if (tempNewick.empty()){
                    tempNewick =  "(" + cString + ")";
                } else {
                    tempNewick =  tempNewick + ",(" + cString + ")";
                }
                curNewick.pop();
                curNewick.push(tempNewick);
            } else {
                Result = cString;
            }

        } else {
            int nextNode = edge.at(node2edge.at(cNode).at(cIndx)).at(0);
            if (nextNode == cNode){
                nextNode = edge.at(node2edge.at(cNode).at(cIndx)).at(1);
            }

            curIndx.pop();
            curIndx.push(cIndx+1);

            if (!nodeVisited.at(nextNode)){
                nodeVisited.at(nextNode) = true;

                if (node2edge.at(nextNode).size() == 1){
                    string cString = curNewick.top();
                    curNewick.pop();
                    if (cString.empty()){
                        curNewick.push(tipLabel.at(nextNode));
                    } else {
                        curNewick.push(cString + "," + tipLabel.at(nextNode));
                    }
                } else {
                    curNode.push(nextNode);
                    curIndx.push(0);
                    curNewick.push("");
                }
            }
        }
    }

    return("(" + Result + ");");
}
