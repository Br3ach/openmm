/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2019-2024 Stanford University and the Authors.      *
 * Authors: Peter Eastman                                                     *
 * Contributors:                                                              *
 *                                                                            *
 * This program is free software: you can redistribute it and/or modify       *
 * it under the terms of the GNU Lesser General Public License as published   *
 * by the Free Software Foundation, either version 3 of the License, or       *
 * (at your option) any later version.                                        *
 *                                                                            *
 * This program is distributed in the hope that it will be useful,            *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of             *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the              *
 * GNU Lesser General Public License for more details.                        *
 *                                                                            *
 * You should have received a copy of the GNU Lesser General Public License   *
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.      *
 * -------------------------------------------------------------------------- */

#include "../platforms/common/src/MoleculeConnectivity.h"
#include "openmm/internal/ContextImpl.h"
#include <algorithm>
#include <iostream>
#include <random>
#include <stdexcept>
#include <unordered_set>
using namespace OpenMM;
using namespace std;
struct Graph {
    int particles;
    vector<pair<int,int>> constraints;
    vector<vector<vector<int>>> forces;
};
static vector<vector<int>> legacy(const Graph& graph) {
    vector<unordered_set<int>> sets(graph.particles);
    for(auto edge:graph.constraints){sets[edge.first].insert(edge.second);sets[edge.second].insert(edge.first);}
    for(const auto& force:graph.forces)
        for(const auto& group:force)
            for(int i=1;i<(int)group.size();i++){sets[group[i]].insert(group[i-1]);sets[group[i-1]].insert(group[i]);}
    vector<vector<int>> adjacent(graph.particles);
    for(int i=0;i<graph.particles;i++)adjacent[i].assign(sets[i].begin(),sets[i].end());
    return ContextImpl::findMolecules(graph.particles,adjacent);
}
static vector<vector<int>> compact(const Graph& graph) {
    detail::MoleculeConnectivity connectivity(graph.particles);
    for(auto edge:graph.constraints)connectivity.connect(edge.first,edge.second);
    for(const auto& force:graph.forces)
        for(const auto& group:force)
            for(int i=1;i<(int)group.size();i++)connectivity.connect(group[i],group[i-1]);
    auto result=connectivity.createMolecules();
    if(result!=connectivity.createMolecules())throw runtime_error("Repeated enumeration changed components");
    return result;
}
static void check(const Graph& graph) {
    auto expected=legacy(graph),actual=compact(graph);
    if(expected!=actual)throw runtime_error("Ordered components differ from legacy algorithm");
    int previous=-1,total=0;
    for(const auto& molecule:actual){
        if(molecule.empty()||molecule.front()<=previous||!is_sorted(molecule.begin(),molecule.end()))throw runtime_error("Wrong ordering");
        previous=molecule.front();total+=(int)molecule.size();
    }
    if(total!=graph.particles)throw runtime_error("Incomplete particle coverage");
}
int main() {
    try {
        check({0,{},{}});check({1,{{0,0}},{{{}, {0}, {0,0}}}});
        check({12,{{0,1},{1,0},{2,3},{3,3}},{{{}, {1,2,2,0}, {8}, {11,9,10}, {}},{{4,5},{6,7},{5,6}}}});
        mt19937 random(20261024);
        for(int trial=0;trial<500;trial++){
            Graph graph;graph.particles=1+random()%200;
            for(int i=0;i<trial%170;i++)graph.constraints.push_back({random()%graph.particles,random()%graph.particles});
            graph.forces.resize(trial%6);
            for(auto& force:graph.forces){
                force.resize(random()%30);
                for(auto& group:force){group.resize(random()%8);for(auto& atom:group)atom=random()%graph.particles;}
            }
            check(graph);
            reverse(graph.constraints.begin(),graph.constraints.end());
            for(auto& edge:graph.constraints)swap(edge.first,edge.second);
            check(graph);
        }
        Graph chain;chain.particles=200000;
        for(int i=chain.particles-1;i>0;i--)chain.constraints.push_back({i,i-1});check(chain);
        Graph sites;sites.particles=20;sites.forces={{{3,0,1,2},{4,3,1},{7,5,6},{8,7,6},{12,9,10,11},{13,12,4},{18,14,15,16,17},{19,18,13}}};check(sites);
        detail::MoleculeConnectivity mutableGraph(5);auto separate=mutableGraph.createMolecules();mutableGraph.connect(4,1);mutableGraph.connect(1,4);mutableGraph.connect(4,4);
        if(separate.size()!=5||mutableGraph.createMolecules()!=vector<vector<int>>{{0},{1,4},{2},{3}})throw runtime_error("Union or snapshot semantics changed");
        cout<<"Ordered connectivity passed:1000 randomized orientations, crafted groups,200k chain, nested virtual-site chains"<<endl;return 0;
    }catch(const exception& error){cerr<<error.what()<<endl;return 1;}
}
