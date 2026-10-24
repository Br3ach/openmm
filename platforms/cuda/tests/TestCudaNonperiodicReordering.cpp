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

// Reuse the constrained multi-site fixture and finite full-state hash helpers.
#define main originalMoleculeGroupsMain
#include "TestCudaMoleculeGroups.cpp"
#undef main
int main(int argc,char**argv){try{
    string precision=argc>1?argv[1]:"mixed";InfoState first(1),second(2);System system;vector<Vec3> positions;
    auto* nb=makeSystem(system,positions,first,second);nb->setNonbondedMethod(NonbondedForce::CutoffNonPeriodic);
    VerletIntegrator integrator(.0005);CudaPlatform platform;map<string,string> props;props["Precision"]=precision;props["UseBlockingSync"]="true";props["DeterministicForces"]="true";
    Context context(system,integrator,platform,props);if(!first.context)throw runtime_error("Missing context");auto& cc=*first.context;
    context.setPositions(positions);context.computeVirtualSites();context.setVelocitiesToTemperature(2,20261024);context.applyVelocityConstraints(1e-7);
    reorder(cc);bool changed=false;for(int i=0;i<system.getNumParticles();i++)changed|=cc.getAtomIndex()[i]!=i;if(!changed)throw runtime_error("No nonidentity reordering");
    integrator.step(260);reorder(cc);stringstream checkpoint(ios::in|ios::out|ios::binary);context.createCheckpoint(checkpoint);
    integrator.step(32);reorder(cc);auto state=stateHash(context),index=indexHash(cc);checkpoint.seekg(0);context.loadCheckpoint(checkpoint);integrator.step(32);reorder(cc);
    if(state!=stateHash(context)||index!=indexHash(cc))throw runtime_error("Nonperiodic checkpoint continuation differs");
    cout<<"{\"precision\":\""<<precision<<"\",\"state\":\""<<hex<<state<<"\",\"atom_index\":\""<<index<<"\",\"trace\":\""<<traceHash(first,second)<<"\",\"self_replay_verified\":true}"<<endl;return 0;
}catch(const exception&e){cerr<<e.what()<<endl;return 1;}}
