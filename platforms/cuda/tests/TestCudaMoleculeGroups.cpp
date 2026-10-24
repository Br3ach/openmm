/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2026 the Authors.           *
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
#include "CudaPlatform.h"
#include "CudaContext.h"
#include "openmm/common/ContextSelector.h"
#include "openmm/internal/ForceImpl.h"
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <iostream>
#include <numeric>
#include <sstream>
using namespace OpenMM;
using namespace std;
static uint64_t hashInt(uint64_t h,uint64_t value){for(int i=0;i<8;i++)h=(h^((value>>(8*i))&255))*1099511628211ULL;return h;}
static uint64_t hashDouble(uint64_t h,double value){if(!std::isfinite(value))throw runtime_error("Nonfinite state value");uint64_t bits;memcpy(&bits,&value,8);return hashInt(h,bits);}
static vector<int> groupCallbackTrace;
struct InfoState {
    vector<vector<int>> groups;
    vector<int> labels,groupLabels,trace;
    ComputeForceInfo* info=nullptr;
    CudaContext* context=nullptr;
    int id;
    explicit InfoState(int id):id(id){}
};
class RecordedInfo : public ComputeForceInfo {
    InfoState& state;
public:
    explicit RecordedInfo(InfoState& state):state(state){}
    bool areParticlesIdentical(int a,int b) override{return state.labels[a]==state.labels[b];}
    int getNumParticleGroups() override{groupCallbackTrace.push_back(state.id);groupCallbackTrace.push_back(-1);state.trace.push_back(-100-state.id);return (int)state.groups.size();}
    void getParticlesInGroup(int index,vector<int>& particles) override{groupCallbackTrace.push_back(state.id);groupCallbackTrace.push_back(index);state.trace.push_back(index);particles=state.groups[index];}
    bool areGroupsIdentical(int a,int b) override{return state.groupLabels[a]==state.groupLabels[b]&&state.groups[a].size()==state.groups[b].size();}
};
class RecordedForce;
class RecordedImpl : public ForceImpl {
    const RecordedForce& owner;
public:
    explicit RecordedImpl(const RecordedForce& owner):owner(owner){}
    void initialize(ContextImpl& context) override;
    const Force& getOwner() const override;
    double calcForcesAndEnergy(ContextImpl&,bool,bool,int) override{return 0;}
    map<string,double> getDefaultParameters() override{return {};}
    vector<string> getKernelNames() override{return {};}
};
class RecordedForce : public Force {
public:
    InfoState& first;InfoState& second;
    RecordedForce(InfoState& first,InfoState& second):first(first),second(second){}
    bool usesPeriodicBoundaryConditions() const override{return false;}
protected:
    ForceImpl* createImpl() const override{return new RecordedImpl(*this);}
};
void RecordedImpl::initialize(ContextImpl& context){
    auto* data=static_cast<CudaPlatform::PlatformData*>(context.getPlatformData());
    for(auto* cc:data->contexts){owner.first.context=cc;owner.second.context=cc;owner.first.info=new RecordedInfo(owner.first);owner.second.info=new RecordedInfo(owner.second);cc->addForce(owner.first.info);cc->addForce(owner.second.info);}
}
const Force& RecordedImpl::getOwner() const{return owner;}
static NonbondedForce* makeSystem(System& system,vector<Vec3>& positions,InfoState& first,InfoState& second){
    auto* nb=new NonbondedForce();nb->setNonbondedMethod(NonbondedForce::CutoffPeriodic);nb->setCutoffDistance(1.0);nb->setUseDispersionCorrection(false);
    system.setDefaultPeriodicBoxVectors(Vec3(6,0,0),Vec3(0,6,0),Vec3(0,0,6));
    first.groups.push_back({});first.groupLabels.push_back(0);second.groups.push_back({});second.groupLabels.push_back(0);
    for(int block=0;block<12;block++){
        int base=system.getNumParticles();double d=block==11?0.12:0.1;Vec3 origin(0.5+((11-block)%4)*0.65,0.5+((11-block)/4)*0.65,0.5);
        vector<Vec3> local={origin,origin+Vec3(d,0,0),origin+Vec3(-0.28*d,0.96*d,0),origin};
        for(int j=0;j<4;j++){
            system.addParticle(j==0?16.0:j==3?0.0:(block==10?2.0:1.0));positions.push_back(local[j]);nb->addParticle(0.0,0.1,0.001);first.labels.push_back(j);second.labels.push_back(j);
        }
        system.addConstraint(base,base+1,d);system.addConstraint(base+2,base, d);system.addConstraint(base+1,base+2,1.6*d);
        if(block==9)system.setVirtualSite(base+3,new TwoParticleAverageSite(base,base+1,0.3,0.7));
        else system.setVirtualSite(base+3,new ThreeParticleAverageSite(base,base+1,base+2,0.6,0.2,0.2));
        for(int a=0;a<4;a++)for(int b=a+1;b<4;b++)nb->addException(base+a,base+b,0,1,0);
        first.groups.push_back({base});first.groupLabels.push_back(1);
        first.groups.push_back({base+2,base+1,base+2,base});first.groupLabels.push_back(2);
        second.groups.push_back({base+3,base,base+1,base+2});second.groupLabels.push_back(3);
    }
    // Nested site and a ForceInfo-only cross-water connection, absent from ForceImpl bonds.
    int site=system.getNumParticles();system.addParticle(0);positions.push_back(Vec3());nb->addParticle(0,0.1,0);first.labels.push_back(4);second.labels.push_back(4);
    system.setVirtualSite(site,new TwoParticleAverageSite(3,0,0.5,0.5));
    for(int j=0;j<4;j++)nb->addException(site,j,0,1,0);
    first.groups.push_back({3,7,3});first.groupLabels.push_back(4);first.groups.push_back({});first.groupLabels.push_back(0);
    system.addForce(nb);system.addForce(new RecordedForce(first,second));return nb;
}
static void reorder(CudaContext& cc){ContextSelector select(cc);cc.forceReorder();cc.reorderAtoms();}
static bool invalidate(CudaContext& cc,ComputeForceInfo* info){ContextSelector select(cc);return cc.invalidateMolecules(info,true,true);}
static uint64_t stateHash(Context& context){
    State s=context.getState(State::Positions|State::Velocities|State::Forces|State::Energy|State::Parameters);uint64_t h=14695981039346656037ULL;
    for(const auto* values:{&s.getPositions(),&s.getVelocities(),&s.getForces()})for(auto v:*values)for(int d=0;d<3;d++)h=hashDouble(h,v[d]);
    h=hashDouble(hashDouble(hashDouble(h,s.getTime()),s.getPotentialEnergy()),s.getKineticEnergy());return hashInt(h,s.getStepCount());
}
static uint64_t indexHash(CudaContext& cc){auto order=cc.getAtomIndex();auto sorted=order;sort(sorted.begin(),sorted.end());for(int i=0;i<(int)sorted.size();i++)if(sorted[i]!=i)throw runtime_error("Atom index is not a permutation");uint64_t h=14695981039346656037ULL;for(int i:order)h=hashInt(h,i);return h;}
static uint64_t traceHash(const InfoState& first,const InfoState& second){uint64_t h=14695981039346656037ULL;for(int x:first.trace)h=hashInt(h,(uint32_t)x);h=hashInt(h,123456);for(int x:second.trace)h=hashInt(h,(uint32_t)x);for(int x:groupCallbackTrace)h=hashInt(h,(uint32_t)x);return h;}
int main(int argc,char**argv){try{
    string precision=argc>1?argv[1]:"mixed",mode=argc>2?argv[2]:"fresh";if(mode!="fresh"&&mode!="save"&&mode!="load")throw runtime_error("Invalid mode");if(mode!="fresh"&&argc!=4)throw runtime_error("Checkpoint path required");
    InfoState first(1),second(2);System system;vector<Vec3> positions;NonbondedForce* nb=makeSystem(system,positions,first,second);
    VerletIntegrator integrator(0.0005);CudaPlatform platform;map<string,string> properties;properties["Precision"]=precision;properties["UseBlockingSync"]="true";properties["DeterministicForces"]="true";
    Context context(system,integrator,platform,properties);if(!first.context)throw runtime_error("Missing CUDA test context");CudaContext& cc=*first.context;
    if(!cc.getNonbondedUtilities().getUseCutoff()||!first.info||!second.info)throw runtime_error("Reorder/ForceInfo setup missing");
    size_t forceCount=cc.getForceInfos().size();context.setPositions(positions);context.computeVirtualSites();context.setVelocitiesToTemperature(2.0,20261024);context.applyVelocityConstraints(1e-7);
    auto before=context.getState(State::Positions).getPositions();reorder(cc);auto after=context.getState(State::Positions).getPositions();
    for(int i=0;i<(int)before.size();i++)for(int d=0;d<3;d++)if(!std::isfinite(before[i][d])||!std::isfinite(after[i][d])||abs(before[i][d]-after[i][d])>1e-6)throw runtime_error("Reorder changed positions");
    bool changed=false;for(int i=0;i<system.getNumParticles();i++)changed|=cc.getAtomIndex()[i]!=i;if(!changed)throw runtime_error("No nonidentity reordering occurred");
    uint64_t sequence=hashInt(stateHash(context),indexHash(cc));integrator.step(260);reorder(cc);sequence=hashInt(hashInt(sequence,stateHash(context)),indexHash(cc));int regroups=0;
    auto record=[&](){if(cc.getForceInfos().size()!=forceCount)throw runtime_error("Regroup reconstructed ForceInfo cache");sequence=hashInt(hashInt(sequence,stateHash(context)),indexHash(cc));};
    first.labels[16]=99;if(!invalidate(cc,first.info))throw runtime_error("Particle split was not detected");regroups++;record();
    first.labels[16]=0;first.labels[20]=99;if(!invalidate(cc,first.info))throw runtime_error("Second particle split was not detected");regroups++;record();
    first.labels[20]=0;second.groupLabels[7]=99;if(!invalidate(cc,second.info))throw runtime_error("Group split was not detected");regroups++;record();
    second.groupLabels[7]=3;double q,sigma,epsilon;nb->getParticleParameters(32,q,sigma,epsilon);nb->setParticleParameters(32,q,sigma,2*epsilon);nb->updateParametersInContext(context);record();
    nb->setParticleParameters(32,q,sigma,epsilon);nb->updateParametersInContext(context);first.labels[28]=99;if(!invalidate(cc,first.info))throw runtime_error("Restored-group merge trigger was not detected");regroups++;record();
    if(mode=="load"){ifstream file(argv[3],ios::binary);if(!file)throw runtime_error("Checkpoint missing");context.loadCheckpoint(file);}
    stringstream checkpoint(ios::in|ios::out|ios::binary);context.createCheckpoint(checkpoint);
    if(mode=="save"){ofstream file(argv[3],ios::binary);file<<checkpoint.rdbuf();if(!file)throw runtime_error("Checkpoint write failed");}
    integrator.step(32);reorder(cc);uint64_t finalState=stateHash(context),finalIndex=indexHash(cc);
    if(mode!="load"){checkpoint.clear();checkpoint.seekg(0);context.loadCheckpoint(checkpoint);integrator.step(32);reorder(cc);if(finalState!=stateHash(context)||finalIndex!=indexHash(cc))throw runtime_error("Reordered checkpoint replay differs");}
    cout<<"{\"precision\":\""<<precision<<"\",\"sequence\":\""<<hex<<sequence<<"\",\"state\":\""<<finalState<<"\",\"atom_index\":\""<<finalIndex<<"\",\"trace\":\""<<traceHash(first,second)<<dec<<"\",\"trace_entries\":"<<first.trace.size()+second.trace.size()<<",\"regroups\":"<<regroups<<",\"force_infos\":"<<forceCount<<",\"nonidentity_reorder\":true,\"self_replay_verified\":"<<(mode=="load"?"false":"true")<<",\"mode\":\""<<mode<<"\"}"<<endl;return 0;
}catch(const exception&error){cerr<<error.what()<<endl;return 1;}}
