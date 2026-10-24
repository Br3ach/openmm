/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2015 Stanford University and the Authors.           *
 * Authors: Peter Eastman                                                     *
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
#include "CudaContext.h"
#include "CudaPlatform.h"
#include "openmm/common/ContextSelector.h"
#include "openmm/common/IntegrationUtilities.h"
#include "openmm/internal/ContextImpl.h"
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <streambuf>
using namespace OpenMM;using namespace std;
class AccessPlatform : public CudaPlatform { public:CudaContext& cuda(Context& c){return *static_cast<PlatformData*>(getContextImpl(c).getPlatformData())->contexts[0];} };
static PROCESS_MEMORY_COUNTERS_EX memory(){PROCESS_MEMORY_COUNTERS_EX m={};m.cb=sizeof(m);if(!GetProcessMemoryInfo(GetCurrentProcess(),(PROCESS_MEMORY_COUNTERS*)&m,sizeof(m)))throw runtime_error("Memory query");return m;}
static void printMemory(const char* phase,const PROCESS_MEMORY_COUNTERS_EX&m){cout<<",\""<<phase<<"_working_set_kib\":"<<m.WorkingSetSize/1024<<",\""<<phase<<"_peak_working_set_kib\":"<<m.PeakWorkingSetSize/1024<<",\""<<phase<<"_private_kib\":"<<m.PrivateUsage/1024;}
class ZeroInput : public streambuf {public:size_t consumed=0;protected:streamsize xsgetn(char* dest,streamsize n) override{memset(dest,0,(size_t)n);consumed+=(size_t)n;return n;}};
int main(int argc,char**argv){try{
 if(argc!=3)throw runtime_error("Usage: BenchmarkCudaHostMemory particles startup|validate|rngsave|rngload|reorder");int n=atoi(argv[1]);string mode=argv[2];
 if(n<4||n>2000000||n%4||!(mode=="startup"||mode=="validate"||mode=="rngsave"||mode=="rngload"||mode=="reorder"))throw runtime_error("Invalid bounded case");
 HANDLE job=CreateJobObjectW(NULL,NULL);JOBOBJECT_EXTENDED_LIMIT_INFORMATION limit={};limit.BasicLimitInformation.LimitFlags=JOB_OBJECT_LIMIT_PROCESS_MEMORY;limit.ProcessMemoryLimit=6ULL*1024*1024*1024;if(!job||!SetInformationJobObject(job,JobObjectExtendedLimitInformation,&limit,sizeof(limit))||!AssignProcessToJobObject(job,GetCurrentProcess()))throw runtime_error("Cannot establish6GiB cap");
 const bool water=mode=="startup"||mode=="validate"||mode=="reorder";System system;vector<Vec3> positions;NonbondedForce* nb=NULL;int side=(int)ceil(cbrt(n/4.0));double box=(side+2)*.5;
 if(mode=="reorder"){nb=new NonbondedForce();nb->setNonbondedMethod(NonbondedForce::CutoffPeriodic);nb->setCutoffDistance(.3);nb->setUseDispersionCorrection(false);system.setDefaultPeriodicBoxVectors(Vec3(box,0,0),Vec3(0,box,0),Vec3(0,0,box));positions.resize(n);}
 for(int i=0;i<n;i++){system.addParticle(water?(i%4==0?16.:i%4==3?0.:1.):1.);if(nb)nb->addParticle(0,.1,0);}
 if(water)for(int i=0;i<n;i+=4){system.addConstraint(i,i+1,.1);system.addConstraint(i+1,i+2,.16);system.addConstraint(i+2,i,.1);system.setVirtualSite(i+3,new ThreeParticleAverageSite(i,i+1,i+2,.6,.2,.2));
  if(nb){int j=n/4-1-i/4;Vec3 origin(.5+(j%side)*.5,.5+((j/side)%side)*.5,.5+(j/(side*side))*.5);positions[i]=origin;positions[i+1]=origin+Vec3(.1,0,0);positions[i+2]=origin+Vec3(-.028,.096,0);positions[i+3]=origin;for(int a=0;a<4;a++)for(int b=a+1;b<4;b++)nb->addException(i+a,i+b,0,1,0);}
 }
 if(nb)system.addForce(nb);auto before=memory();VerletIntegrator integrator(.001);AccessPlatform platform;map<string,string> props;props["Precision"]="mixed";props["UseBlockingSync"]="true";
 auto start=chrono::steady_clock::now();Context context(system,integrator,platform,props);double constructorSeconds=chrono::duration<double>(chrono::steady_clock::now()-start).count();auto constructed=memory();CudaContext& cc=platform.cuda(context);ContextSelector select(cc);size_t randomBytes=0,ioBytes=0;int repeats=mode=="validate"?5:mode=="reorder"?3:1;
 if(mode=="rngsave"||mode=="rngload"){auto& utilities=cc.getIntegrationUtilities();utilities.initRandomNumberGenerator(20261024);utilities.prepareRandomNumbers(4*cc.getPaddedNumAtoms());if(cuStreamSynchronize(cc.getCurrentStream())!=CUDA_SUCCESS)throw runtime_error("RNG sync failed");randomBytes=utilities.getRandom().getSize()*16;}
 auto phaseBefore=memory();start=chrono::steady_clock::now();
 if(mode=="validate")for(int i=0;i<repeats;i++)cc.validateAtomOrder();
 else if(mode=="rngsave"){ofstream out("NUL",ios::binary);cc.getIntegrationUtilities().createCheckpoint(out);if(!out)throw runtime_error("NUL checkpoint write failed");}
 else if(mode=="rngload"){ZeroInput buffer;istream in(&buffer);cc.getIntegrationUtilities().loadCheckpoint(in);ioBytes=buffer.consumed;size_t expected=sizeof(int)+randomBytes+cc.getNumThreadBlocks()*64*16; if(ioBytes!=expected||!in)throw runtime_error("Checkpoint input size mismatch");}
 else if(mode=="reorder"){context.setPositions(positions);context.computeVirtualSites();for(int i=0;i<repeats;i++){cc.forceReorder();cc.reorderAtoms();}}
 double phaseSeconds=chrono::duration<double>(chrono::steady_clock::now()-start).count();auto after=memory();
 // Verify after recording the measured highwater, so verification buffers are excluded.
 if(nb){auto state=context.getState(State::Positions);for(int i=0;i<n;i++)for(int d=0;d<3;d++)if(!isfinite(state.getPositions()[i][d]))throw runtime_error("Nonfinite reordered position");}
 cout<<setprecision(10)<<"{\"particles\":"<<n<<",\"mode\":\""<<mode<<"\",\"constraints\":"<<system.getNumConstraints()<<",\"precision\":\"mixed\",\"construction_seconds\":"<<constructorSeconds<<",\"phase_seconds\":"<<phaseSeconds<<",\"repeats\":"<<repeats<<",\"random_bytes\":"<<randomBytes<<",\"input_bytes\":"<<ioBytes;
 printMemory("before",before);printMemory("constructed",constructed);printMemory("phase_before",phaseBefore);printMemory("phase_after",after);cout<<",\"peak_scope\":\"process_lifetime_through_measured_phase_before_verification\"}"<<endl;CloseHandle(job);return 0;
}catch(const exception&e){cerr<<e.what()<<endl;return 1;}}
