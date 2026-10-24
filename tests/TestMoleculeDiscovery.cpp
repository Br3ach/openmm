/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2026 the Authors.                                *
 * Authors: OpenMM contributors                                                     *
 * Contributors:                                                              *
 *                                                                            *
 * Permission is hereby granted, free of charge, to any person obtaining a    *
 * copy of this software and associated documentation files (the "Software"), *
 * to deal in the Software without restriction, including without limitation  *
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,   *
 * and/or sell copies of the Software, and to permit persons to whom the      *
 * Software is furnished to do so, subject to the following conditions:       *
 *                                                                            *
 * The above copyright notice and this permission notice shall be included in *
 * all copies or substantial portions of the Software.                        *
 *                                                                            *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR *
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,   *
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL    *
 * THE AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,    *
 * DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR      *
 * OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE  *
 * USE OR OTHER DEALINGS IN THE SOFTWARE.                                     *
 * -------------------------------------------------------------------------- */

#include "OpenMM.h"
#include "openmm/internal/ContextImpl.h"
#include "openmm/internal/ForceImpl.h"
#include "openmm/internal/AssertionUtilities.h"
#include <algorithm>
#include <iostream>
#include <random>
using namespace OpenMM;
using namespace std;

class GraphForce;
class GraphForceImpl : public ForceImpl {
public:
    explicit GraphForceImpl(const GraphForce& owner) : owner(owner) {}
    void initialize(ContextImpl&) override {}
    const Force& getOwner() const override;
    double calcForcesAndEnergy(ContextImpl&, bool, bool, int) override { return 0; }
    map<string,double> getDefaultParameters() override { return {}; }
    vector<string> getKernelNames() override { return {}; }
    vector<pair<int,int> > getBondedParticles() const override;
private:
    const GraphForce& owner;
};
class GraphForce : public Force {
public:
    vector<pair<int,int> > edges;
    vector<int>& calls;
    int id;
    mutable bool fail=false;
    GraphForce(vector<int>& calls, int id) : calls(calls), id(id) {}
    bool usesPeriodicBoundaryConditions() const override { return false; }
protected:
    ForceImpl* createImpl() const override { return new GraphForceImpl(*this); }
};
const Force& GraphForceImpl::getOwner() const { return owner; }
vector<pair<int,int> > GraphForceImpl::getBondedParticles() const {
    owner.calls.push_back(owner.id);
    if(owner.fail)throw OpenMMException("Intentional bonded-particle callback failure");
    return owner.edges;
}
static vector<vector<int> > expected(int n,const vector<pair<int,int> >& edges) {
    vector<vector<int> > adjacent(n);
    for(auto e:edges){adjacent[e.first].push_back(e.second);adjacent[e.second].push_back(e.first);}
    return ContextImpl::findMolecules(n,adjacent);
}
static void check(const vector<vector<int> >& actual,const vector<vector<int> >& reference) {
    ASSERT_EQUAL_CONTAINERS(reference,actual);
    int previous=-1;
    for(auto& molecule:actual){
        ASSERT(!molecule.empty());ASSERT(molecule.front()>previous);previous=molecule.front();
        ASSERT(is_sorted(molecule.begin(),molecule.end()));
    }
}
static void testGraphs() {
    mt19937 random(20261024);
    for(int repeat=0;repeat<50;repeat++) {
        int n=1+random()%200;vector<int> calls;System system;
        for(int i=0;i<n;i++)system.addParticle(1);
        auto* a=new GraphForce(calls,1);auto* b=new GraphForce(calls,2);
        system.addForce(a);system.addForce(b);
        vector<pair<int,int> > all;
        for(int i=0;i<repeat*2;i++){
            pair<int,int> edge(random()%n,random()%n);(i%2?a:b)->edges.push_back(edge);all.push_back(edge);
            if(i%4==0){a->edges.push_back(edge);all.push_back(edge);}
        }
        VerletIntegrator integrator(0.001);Context context(system,integrator,Platform::getPlatformByName("Reference"));
        ASSERT(calls.empty());
        auto& result=context.getMolecules();check(result,expected(n,all));
        ASSERT_EQUAL_CONTAINERS((vector<int>{1,2}),calls);
        ASSERT(&result==&context.getMolecules());ASSERT_EQUAL(2,calls.size());
        if(n>1){b->edges.push_back({0,n-1});all.push_back({0,n-1});}
        context.reinitialize();ASSERT_EQUAL(2,calls.size());
        check(context.getMolecules(),expected(n,all));
        ASSERT_EQUAL_CONTAINERS((vector<int>{1,2,1,2}),calls);
    }
    // A long component must not introduce recursive stack use.
    vector<int> calls;System system;const int n=20000;
    auto* force=new GraphForce(calls,1);system.addForce(force);
    for(int i=0;i<n;i++){system.addParticle(1);if(i>0)force->edges.push_back({i,i-1});}
    reverse(force->edges.begin(),force->edges.end());
    VerletIntegrator integrator(0.001);Context context(system,integrator,Platform::getPlatformByName("Reference"));
    check(context.getMolecules(),expected(n,force->edges));
}
static void testConstraintsAndSites() {
    vector<int> calls;System system;const int n=15;
    for(int i=0;i<n;i++)system.addParticle(i>=11?0:1);
    system.addConstraint(8,9,1);
    system.setVirtualSite(11,new TwoParticleAverageSite(0,1,0.5,0.5));
    system.setVirtualSite(12,new ThreeParticleAverageSite(2,3,4,0.2,0.3,0.5));
    system.setVirtualSite(13,new OutOfPlaneSite(11,5,6,0.2,0.3,0.1));
    system.setVirtualSite(14,new LocalCoordinatesSite(12,13,7,Vec3(1,0,0),Vec3(-1,1,0),Vec3(-1,0,1),Vec3(0.1,0.2,0.3)));
    auto* a=new GraphForce(calls,1);auto* empty=new GraphForce(calls,2);
    a->edges={{9,10},{10,9},{10,10}};system.addForce(a);system.addForce(empty);
    vector<pair<int,int> > all=a->edges;all.push_back({8,9});
    for(int i=0;i<n;i++)if(system.isVirtualSite(i)){
        auto& site=system.getVirtualSite(i);for(int j=0;j<site.getNumParticles();j++)all.push_back({i,site.getParticle(j)});
    }
    VerletIntegrator integrator(0.001);Context context(system,integrator,Platform::getPlatformByName("Reference"));
    check(context.getMolecules(),expected(n,all));ASSERT_EQUAL_CONTAINERS((vector<int>{1,2}),calls);
}
static void testCallbackFailure() {
    vector<int> calls;System system;for(int i=0;i<6;i++)system.addParticle(1);
    auto* a=new GraphForce(calls,1);auto* b=new GraphForce(calls,2);
    a->edges={{4,5},{3,4}};b->edges={{0,2},{1,3},{2,5}};b->fail=true;
    system.addForce(a);system.addForce(b);
    VerletIntegrator integrator(0.001);Context context(system,integrator,Platform::getPlatformByName("Reference"));
    bool threw=false;try{context.getMolecules();}catch(const OpenMMException&){threw=true;}ASSERT(threw);
    b->fail=false;vector<pair<int,int> > all=a->edges;all.insert(all.end(),b->edges.begin(),b->edges.end());
    check(context.getMolecules(),expected(6,all));
    ASSERT_EQUAL_CONTAINERS((vector<vector<int> >{{0,1,2,3,4,5}}),context.getMolecules());
    ASSERT_EQUAL_CONTAINERS((vector<int>{1,2,1,2}),calls);
}
int main(){try{testGraphs();testConstraintsAndSites();testCallbackFailure();cout<<"Ordered molecule discovery, callbacks, cache, reinitialization and virtual sites passed"<<endl;return 0;}
catch(const exception& e){cerr<<e.what()<<endl;return 1;}}
