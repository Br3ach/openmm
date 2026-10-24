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
#include "CudaPlatform.h"
#include <chrono>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <memory>
using namespace OpenMM;
using namespace std;
static PROCESS_MEMORY_COUNTERS_EX memory() {
    PROCESS_MEMORY_COUNTERS_EX result={};result.cb=sizeof(result);
    if(!GetProcessMemoryInfo(GetCurrentProcess(),(PROCESS_MEMORY_COUNTERS*)&result,sizeof(result)))throw runtime_error("Memory query failed");
    return result;
}
class MeasuredCudaPlatform : public CudaPlatform {
public:
    mutable PROCESS_MEMORY_COUNTERS_EX created={};
    void contextCreated(ContextImpl& context,const map<string,string>& properties) const override {
        CudaPlatform::contextCreated(context,properties);created=memory();
    }
};
static void printMemory(const char* phase,const PROCESS_MEMORY_COUNTERS_EX& m) {
    cout<<",\""<<phase<<"_working_set_kib\":"<<m.WorkingSetSize/1024
        <<",\""<<phase<<"_peak_working_set_kib\":"<<m.PeakWorkingSetSize/1024
        <<",\""<<phase<<"_private_kib\":"<<m.PrivateUsage/1024;
}
int main(int argc,char** argv) {
    try {
        if(argc!=3)throw runtime_error("Usage: BenchmarkCudaSettleStartup particles water|unconstrained");
        int count=atoi(argv[1]);string kind=argv[2];
        if(count<4||count>500000||count%4!=0||(kind!="water"&&kind!="unconstrained"))throw runtime_error("Invalid bounded case");
        HANDLE job=CreateJobObjectW(nullptr,nullptr);JOBOBJECT_EXTENDED_LIMIT_INFORMATION limit={};limit.BasicLimitInformation.LimitFlags=JOB_OBJECT_LIMIT_PROCESS_MEMORY;limit.ProcessMemoryLimit=6ULL*1024*1024*1024;
        if(!job||!SetInformationJobObject(job,JobObjectExtendedLimitInformation,&limit,sizeof(limit))||!AssignProcessToJobObject(job,GetCurrentProcess()))throw runtime_error("Could not establish6GiB process limit");
        System system;
        for(int i=0;i<count;i++)system.addParticle(kind=="water"?(i%4==0?16.0:(i%4==3?0.0:1.0)):1.0);
        if(kind=="water")for(int i=0;i<count;i+=4) {
            system.addConstraint(i,i+1,0.1);system.addConstraint(i+1,i+2,0.16);system.addConstraint(i+2,i,0.1);
            system.setVirtualSite(i+3,new ThreeParticleAverageSite(i,i+1,i+2,0.6,0.2,0.2));
        }
        VerletIntegrator integrator(0.001);MeasuredCudaPlatform platform;
        map<string,string> properties;properties["Precision"]="mixed";properties["UseBlockingSync"]="true";
        auto before=memory();auto start=chrono::steady_clock::now();
        unique_ptr<Context> context(new Context(system,integrator,platform,properties));
        double seconds=chrono::duration<double>(chrono::steady_clock::now()-start).count();auto constructed=memory();
        if(context->getTime()!=0.0||platform.created.cb==0)throw runtime_error("Unexpected Context");
        context.reset();auto destroyed=memory();
        cout<<setprecision(10)<<"{\"particles\":"<<count<<",\"constraints\":"<<system.getNumConstraints()<<",\"kind\":\""<<kind<<"\",\"precision\":\"mixed\",\"construction_seconds\":"<<seconds;
        printMemory("before",before);printMemory("platform_created",platform.created);printMemory("constructed",constructed);printMemory("destroyed",destroyed);
        cout<<",\"peak_scope\":\"process_lifetime\",\"platform\":\"CUDA\"}"<<endl;
        CloseHandle(job);return 0;
    } catch(const exception& error) {cerr<<error.what()<<endl;return 1;}
}
