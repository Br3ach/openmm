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
#include "CudaContext.h"
#include "CudaPlatform.h"
#include "../src/CudaCheckpointUtilities.h"
#include "openmm/common/IntegrationUtilities.h"
#include "openmm/internal/ContextImpl.h"
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <iostream>
#include <sstream>
using namespace OpenMM;using namespace std;
class AccessPlatform : public CudaPlatform { public:
    CudaContext& cuda(Context& c){return *static_cast<PlatformData*>(getContextImpl(c).getPlatformData())->contexts[0];}
};
static uint64_t hashBytes(const string& s){uint64_t h=14695981039346656037ULL;for(unsigned char c:s)h=(h^c)*1099511628211ULL;return h;}
static uint64_t stateHash(Context& context){auto s=context.getState(State::Positions|State::Velocities|State::Forces|State::Energy);string bytes;
 for(const auto* values:{&s.getPositions(),&s.getVelocities(),&s.getForces()})for(auto v:*values)for(int d=0;d<3;d++){if(!isfinite(v[d]))throw runtime_error("Nonfinite state");double x=v[d];bytes.append((char*)&x,8);}
 for(double x:{s.getTime(),s.getPotentialEnergy(),s.getKineticEnergy()}){if(!isfinite(x))throw runtime_error("Nonfinite energy/time");bytes.append((char*)&x,8);}return hashBytes(bytes);
}
class LimitedOutput : public streambuf { public:size_t remaining,written=0;explicit LimitedOutput(size_t n):remaining(n){}protected:
 streamsize xsputn(const char*,streamsize n) override{size_t count=min(remaining,(size_t)n);remaining-=count;written+=count;return count;}
};
static void require(bool ok,const char* what){if(!ok)throw runtime_error(what);}
static void cudaCheck(CUresult r){if(r!=CUDA_SUCCESS)throw runtime_error(CudaContext::getErrorString(r));}
static void testErrors(const string& precision){
 System system;system.addParticle(1);VerletIntegrator integrator(.001);AccessPlatform platform;map<string,string> props;props["Precision"]=precision;props["UseBlockingSync"]="true";Context context(system,integrator,platform,props);CudaContext& cc=platform.cuda(context);ContextSelector selector(cc);
 stringstream none(ios::in|ios::out|ios::binary);cc.getIntegrationUtilities().createCheckpoint(none);require(none.str().empty(),"Uninitialized RNG wrote data");cc.getIntegrationUtilities().loadCheckpoint(none);
 ComputeQueue previous=cc.getCurrentQueue();ComputeQueue queue=cc.createQueue();cc.setCurrentQueue(queue);
 try {
  const size_t chunk=detail::CudaCheckpointChunkBytes;
  for(size_t count:{size_t(1),chunk-1,chunk,chunk+1,2*chunk+3}){
   CudaArray array(cc,count,1,"checkpoint-boundary");cudaCheck(cuMemsetD8Async(array.getDevicePointer(),0x5a,count,cc.getCurrentStream()));
   ostringstream out(ios::out|ios::binary);detail::writeCudaCheckpointArray(array,out);require(out.str()==string(count,char(0x5a)),"Chunk output/race mismatch");
   istringstream in(string(count,char(0xa5)),ios::in|ios::binary);detail::readCudaCheckpointArray(array,in);vector<char> data;array.download(data);require(data==vector<char>(count,char(0xa5)),"Chunk input mismatch");
  }
  const size_t count=2*chunk+3;CudaArray array(cc,count,1,"checkpoint-failure");
  for(bool exceptions:{false,true}){
   cudaCheck(cuMemsetD8Async(array.getDevicePointer(),0x39,count,cc.getCurrentStream()));
   istringstream shortInput(string(chunk+17,char(0xa5)),ios::in|ios::binary);if(exceptions)shortInput.exceptions(ios::badbit|ios::failbit);bool failed=false;
   try{detail::readCudaCheckpointArray(array,shortInput);}catch(const exception&){failed=true;}require(failed,"Short read accepted");
   vector<char> data;array.download(data);for(size_t i=0;i<count;i++)require(data[i]==char(i<chunk?0xa5:0x39),"Incomplete chunk was uploaded");
   LimitedOutput sink(chunk+17);ostream shortOutput(&sink);if(exceptions)shortOutput.exceptions(ios::badbit|ios::failbit);failed=false;
   try{detail::writeCudaCheckpointArray(array,shortOutput);}catch(const exception&){failed=true;}require(failed&&!shortOutput.good()&&sink.written==chunk+17,"Short write accepted");
  }
 }catch(...){cc.setCurrentQueue(previous);throw;}
 cc.setCurrentQueue(previous);queue.reset();
 cout<<"{\"precision\":\""<<precision<<"\",\"chunk_boundaries\":5,\"short_reads\":2,\"short_writes\":2,\"nondefault_stream\":true,\"uninitialized_rng\":true}"<<endl;
}
int main(int argc,char**argv){try{
 string precision=argc>1?argv[1]:"mixed",mode=argc>2?argv[2]:"fresh";if(mode=="errors"){testErrors(precision);return 0;}
 if(mode!="fresh"&&mode!="save"&&mode!="load")throw runtime_error("Invalid mode");if(mode!="fresh"&&argc!=4)throw runtime_error("Checkpoint path required");
 System system;vector<Vec3> positions;for(int i=0;i<257;i++){system.addParticle(1);positions.push_back(Vec3(.1*(i%7),.1*((i/7)%7),.1*(i/49)));}
 LangevinMiddleIntegrator integrator(2,1,.0005);integrator.setRandomNumberSeed(20261024);AccessPlatform platform;map<string,string> props;props["Precision"]=precision;props["UseBlockingSync"]="true";props["DeterministicForces"]="true";
 Context context(system,integrator,platform,props);context.setPositions(positions);context.setVelocitiesToTemperature(2,20261024);integrator.step(3);CudaContext& cc=platform.cuda(context);
 {ContextSelector select(cc);cc.getIntegrationUtilities().prepareRandomNumbers(2*(detail::CudaCheckpointChunkBytes/16)+3);cudaCheck(cuStreamSynchronize(cc.getCurrentStream()));}
 if(mode=="load"){ifstream file(argv[3],ios::binary);if(!file)throw runtime_error("Missing checkpoint");context.loadCheckpoint(file);}
 stringstream checkpoint(ios::in|ios::out|ios::binary);context.createCheckpoint(checkpoint);string original=checkpoint.str();
 if(mode=="save"){ofstream file(argv[3],ios::binary);file.write(original.data(),original.size());require(bool(file),"Write failed");}
 integrator.step(7);uint64_t state=stateHash(context);
 if(mode!="load"){checkpoint.clear();checkpoint.seekg(0);context.loadCheckpoint(checkpoint);integrator.step(7);require(state==stateHash(context),"RNG checkpoint continuation differs");}
 cout<<"{\"precision\":\""<<precision<<"\",\"state\":\""<<hex<<state<<"\",\"checkpoint\":\""<<hashBytes(original)<<dec<<"\",\"checkpoint_bytes\":"<<original.size()<<",\"grown_rng\":true,\"self_replay_verified\":"<<(mode=="load"?"false":"true")<<"}"<<endl;return 0;
}catch(const exception&e){cerr<<e.what()<<endl;return 1;}}
