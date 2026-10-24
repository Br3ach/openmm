/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2012-2015 Stanford University and the Authors.      *
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
#include <chrono>
#include <cstdint>
#include <cstring>
#include <cstdlib>
#include <iostream>
#include <iomanip>
using namespace OpenMM;
using namespace std;
static PROCESS_MEMORY_COUNTERS_EX memory(){PROCESS_MEMORY_COUNTERS_EX m={};m.cb=sizeof(m);if(!GetProcessMemoryInfo(GetCurrentProcess(),(PROCESS_MEMORY_COUNTERS*)&m,sizeof(m)))throw runtime_error("GetProcessMemoryInfo failed");return m;}
static uint64_t hashBytes(uint64_t h,const void* p,size_t n){const unsigned char* b=(const unsigned char*)p;for(size_t i=0;i<n;i++)h=(h^b[i])*1099511628211ULL;return h;}
int main(int argc,char**argv){try{
 if(argc!=3)throw runtime_error("Usage: BenchmarkContextState particles mask");
 HANDLE job=CreateJobObjectW(nullptr,nullptr);JOBOBJECT_EXTENDED_LIMIT_INFORMATION lim={};
 lim.BasicLimitInformation.LimitFlags=JOB_OBJECT_LIMIT_PROCESS_MEMORY;lim.ProcessMemoryLimit=6ULL*1024*1024*1024;
 if(!job||!SetInformationJobObject(job,JobObjectExtendedLimitInformation,&lim,sizeof(lim))||!AssignProcessToJobObject(job,GetCurrentProcess()))throw runtime_error("Could not establish 6 GiB process memory limit");
 const int count=atoi(argv[1]),mask=atoi(argv[2]);if(count<1||count>10000000||mask<1||mask>31)throw runtime_error("Invalid arguments");
 System system;for(int i=0;i<count;i++)system.addParticle(1.0);
 VerletIntegrator integ(0.001);Context context(system,integ,Platform::getPlatformByName("Reference"));
 {vector<Vec3> a(count);for(int i=0;i<count;i++)a[i]=Vec3((i%997)*0.001,(i%991)*0.002,(i%983)*0.003);context.setPositions(a);for(int i=0;i<count;i++)a[i]=Vec3((i%17)*0.0001,(i%19)*-0.0002,(i%23)*0.0003);context.setVelocities(a);}
 auto before=memory();auto start=chrono::steady_clock::now();
 State state=context.getState(mask,false,-1);
 double seconds=chrono::duration<double>(chrono::steady_clock::now()-start).count();auto after=memory();
 if(state.getDataTypes()!=mask)throw runtime_error("Wrong State mask");
 uint64_t h=14695981039346656037ULL;
 if(mask&State::Positions){auto&a=state.getPositions();if(a.size()!=count)throw runtime_error("positions count");h=hashBytes(h,a.data(),a.size()*sizeof(Vec3));}
 if(mask&State::Velocities){auto&a=state.getVelocities();if(a.size()!=count)throw runtime_error("velocities count");h=hashBytes(h,a.data(),a.size()*sizeof(Vec3));}
 if(mask&State::Forces){auto&a=state.getForces();if(a.size()!=count)throw runtime_error("forces count");h=hashBytes(h,a.data(),a.size()*sizeof(Vec3));}
 if(mask&State::Energy){double e[2]={state.getPotentialEnergy(),state.getKineticEnergy()};h=hashBytes(h,e,sizeof(e));}
 cout<<setprecision(10)<<"{\"particles\":"<<count<<",\"mask\":"<<mask<<",\"seconds\":"<<seconds<<",\"baseline_working_set_kib\":"<<before.WorkingSetSize/1024<<",\"baseline_peak_working_set_kib\":"<<before.PeakWorkingSetSize/1024<<",\"peak_working_set_kib\":"<<after.PeakWorkingSetSize/1024<<",\"retained_working_set_kib\":"<<after.WorkingSetSize/1024<<",\"private_kib\":"<<after.PrivateUsage/1024<<",\"fingerprint\":\""<<hex<<h<<"\",\"platform\":\"Reference\",\"peak_scope\":\"process_lifetime\"}"<<endl;
 CloseHandle(job);return 0;
}catch(const exception&e){cerr<<e.what()<<endl;return 1;}}
