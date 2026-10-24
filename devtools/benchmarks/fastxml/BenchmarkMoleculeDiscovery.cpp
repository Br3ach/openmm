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

#define NOMINMAX
#include <windows.h>
#include <psapi.h>
#include "OpenMM.h"
#include "openmm/internal/ContextImpl.h"
#include "openmm/internal/ForceImpl.h"
#include <chrono>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <iomanip>
using namespace OpenMM;
using namespace std;
class GraphForce;
class GraphImpl : public ForceImpl {
    const GraphForce& owner;
public:
    explicit GraphImpl(const GraphForce& owner):owner(owner){}
    void initialize(ContextImpl&) override {}
    const Force& getOwner() const override;
    double calcForcesAndEnergy(ContextImpl&,bool,bool,int) override{return 0;}
    map<string,double> getDefaultParameters() override{return {};}
    vector<string> getKernelNames() override{return {};}
    vector<pair<int,int> > getBondedParticles() const override;
};
class GraphForce : public Force {
public:
    int count;string kind;mutable int calls=0;
    GraphForce(int count,string kind):count(count),kind(kind){}
    bool usesPeriodicBoundaryConditions() const override{return false;}
protected:
    ForceImpl* createImpl() const override{return new GraphImpl(*this);}
};
const Force& GraphImpl::getOwner() const{return owner;}
vector<pair<int,int> > GraphImpl::getBondedParticles() const{
    ++owner.calls;vector<pair<int,int> > edges;
    if(owner.kind=="water"){
        edges.reserve(owner.count);
        for(int i=0;i+2<owner.count;i+=3){edges.push_back({i,i+1});edges.push_back({i+1,i+2});edges.push_back({i+2,i});}
    }else if(owner.kind=="chain"){
        edges.reserve(owner.count-1);for(int i=1;i<owner.count;i++)edges.push_back({i,i-1});
    }
    return edges;
}
static PROCESS_MEMORY_COUNTERS_EX memory(){PROCESS_MEMORY_COUNTERS_EX m={};m.cb=sizeof(m);if(!GetProcessMemoryInfo(GetCurrentProcess(),(PROCESS_MEMORY_COUNTERS*)&m,sizeof(m)))throw runtime_error("GetProcessMemoryInfo failed");return m;}
static uint64_t hashInt(uint64_t h,int value){for(int i=0;i<4;i++)h=(h^((value>>(8*i))&255))*1099511628211ULL;return h;}
int main(int argc,char**argv){try{
    if(argc!=3)throw runtime_error("Usage: BenchmarkMoleculeDiscovery particles water|chain|isolated");
    const int count=atoi(argv[1]);string kind=argv[2];if(count<1||count>10000000||(kind!="water"&&kind!="chain"&&kind!="isolated"))throw runtime_error("Invalid arguments");
    HANDLE job=CreateJobObjectW(nullptr,nullptr);JOBOBJECT_EXTENDED_LIMIT_INFORMATION limit={};limit.BasicLimitInformation.LimitFlags=JOB_OBJECT_LIMIT_PROCESS_MEMORY;limit.ProcessMemoryLimit=6ULL*1024*1024*1024;
    if(!job||!SetInformationJobObject(job,JobObjectExtendedLimitInformation,&limit,sizeof(limit))||!AssignProcessToJobObject(job,GetCurrentProcess()))throw runtime_error("Could not establish 6 GiB process memory limit");
    System system;for(int i=0;i<count;i++)system.addParticle(1);
    auto* force=new GraphForce(count,kind);system.addForce(force);
    VerletIntegrator integrator(0.001);Context context(system,integrator,Platform::getPlatformByName("Reference"));
    if(force->calls!=0)throw runtime_error("Molecules were discovered during setup");
    auto before=memory();auto start=chrono::steady_clock::now();
    const auto& molecules=context.getMolecules();
    double seconds=chrono::duration<double>(chrono::steady_clock::now()-start).count();auto after=memory();
    if(force->calls!=1||&molecules!=&context.getMolecules()||force->calls!=1)throw runtime_error("Unexpected callback/cache behavior");
    int expected=kind=="water"?count/3+count%3:kind=="chain"?1:count;
    if(molecules.size()!=expected)throw runtime_error("Wrong molecule count");
    uint64_t h=hashInt(14695981039346656037ULL,(int)molecules.size());int total=0,previous=-1;
    for(auto& molecule:molecules){
        if(molecule.empty()||molecule.front()<=previous)throw runtime_error("Wrong component order");previous=molecule.front();
        h=hashInt(h,(int)molecule.size());int previousParticle=-1;
        for(int p:molecule){if(p<=previousParticle||p>=count)throw runtime_error("Wrong particle order");previousParticle=p;h=hashInt(h,p);++total;}
    }
    if(total!=count)throw runtime_error("Wrong particle total");
    cout<<setprecision(10)<<"{\"particles\":"<<count<<",\"kind\":\""<<kind<<"\",\"molecules\":"<<molecules.size()<<",\"seconds\":"<<seconds<<",\"baseline_working_set_kib\":"<<before.WorkingSetSize/1024<<",\"baseline_peak_working_set_kib\":"<<before.PeakWorkingSetSize/1024<<",\"peak_working_set_kib\":"<<after.PeakWorkingSetSize/1024<<",\"retained_working_set_kib\":"<<after.WorkingSetSize/1024<<",\"private_kib\":"<<after.PrivateUsage/1024<<",\"fingerprint\":\""<<hex<<h<<"\",\"peak_scope\":\"process_lifetime\",\"platform\":\"Reference\"}"<<endl;
    CloseHandle(job);return 0;
}catch(const exception&e){cerr<<e.what()<<endl;return 1;}}
