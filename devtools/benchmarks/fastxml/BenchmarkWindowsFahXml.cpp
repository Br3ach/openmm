// Native Windows measurements for external, uncompressed System/State XML.
// Reads a stable copy of an input; never open an actively changing checkpoint.
// PeakWorkingSetSize is a process lifetime metric, not an operation-reset peak.
#define NOMINMAX
#include <windows.h>
#include <psapi.h>
#include "openmm/System.h"
#include "openmm/State.h"
#include "openmm/serialization/XmlSerializer.h"
#include <chrono>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
using namespace OpenMM;
using namespace std;
PROCESS_MEMORY_COUNTERS_EX memory() {
    PROCESS_MEMORY_COUNTERS_EX m={};
    if (!GetProcessMemoryInfo(GetCurrentProcess(),reinterpret_cast<PROCESS_MEMORY_COUNTERS*>(&m),sizeof(m)))
        throw runtime_error("GetProcessMemoryInfo failed");
    return m;
}
int main(int argc,char** argv) {
    try {
        if(argc!=6) throw runtime_error("usage: BenchmarkWindowsFahXml kind input output-or-dash fast memory-limit-MiB");
        if(string(argv[1])!="system" && string(argv[1])!="state")
            throw runtime_error("kind must be system or state");
        if(string(argv[4])!="0" && string(argv[4])!="1")
            throw runtime_error("fast must be 0 or 1");
        const auto limitMiB=stoull(argv[5]);
        if(limitMiB<1 || limitMiB>65536)
            throw runtime_error("memory limit must be between 1 and 65536 MiB");
        HANDLE job=CreateJobObjectA(NULL,NULL);
        if(!job) throw runtime_error("CreateJobObject failed");
        JOBOBJECT_EXTENDED_LIMIT_INFORMATION limits={};
        limits.BasicLimitInformation.LimitFlags=JOB_OBJECT_LIMIT_PROCESS_MEMORY;
        limits.ProcessMemoryLimit=limitMiB*1024*1024;
        if(!SetInformationJobObject(job,JobObjectExtendedLimitInformation,&limits,sizeof(limits)) ||
           !AssignProcessToJobObject(job,GetCurrentProcess())) throw runtime_error("Cannot enforce benchmark memory limit");
        _putenv_s("OPENMM_FAST_XML",argv[4]);
        _putenv_s("OPENMM_FAST_SYSTEM_WRITE",argv[4]);
        _putenv_s("OPENMM_FAST_STATE_WRITE",argv[4]);
        bool state=string(argv[1])=="state";
        unique_ptr<System> sy; unique_ptr<State> st;
        ifstream input(argv[2],ios::binary);
        if(!input) throw runtime_error("Cannot open input");
        auto before=memory();
        auto start=chrono::steady_clock::now();
        if(state) st.reset(XmlSerializer::deserialize<State>(input)); else sy.reset(XmlSerializer::deserialize<System>(input));
        double readSeconds=chrono::duration<double>(chrono::steady_clock::now()-start).count();
        auto readMemory=memory();
        size_t count=state ? st->getPositions().size() : sy->getNumParticles();
        double writeSeconds=0;
        if(string(argv[3])!="-") {
            ofstream output(argv[3],ios::binary);
            if(!output) throw runtime_error("Cannot open output");
            start=chrono::steady_clock::now();
            if(state) XmlSerializer::serialize(st.get(),"State",output); else XmlSerializer::serialize(sy.get(),"System",output);
            output.flush(); if(!output) throw runtime_error("XML output failed");
            writeSeconds=chrono::duration<double>(chrono::steady_clock::now()-start).count();
        }
        auto after=memory();
        cout<<setprecision(9)<<"{\"kind\":\""<<argv[1]<<"\",\"fast\":"<<argv[4]<<",\"particles\":"<<count
            <<",\"read_seconds\":"<<readSeconds<<",\"write_seconds\":"<<writeSeconds
            <<",\"baseline_working_set_kib\":"<<before.WorkingSetSize/1024
            <<",\"read_peak_working_set_kib\":"<<readMemory.PeakWorkingSetSize/1024
            <<",\"read_private_kib\":"<<readMemory.PrivateUsage/1024
            <<",\"lifetime_peak_working_set_kib\":"<<after.PeakWorkingSetSize/1024
            <<",\"peak_scope\":\"process_lifetime\"";
        if(!state) cout<<",\"constraints\":"<<sy->getNumConstraints()<<",\"forces\":"<<sy->getNumForces();
        else cout<<",\"data_types\":"<<st->getDataTypes()<<",\"velocities\":"<<st->getVelocities().size()<<",\"step_count\":"<<st->getStepCount();
        cout<<"}"<<endl;
        CloseHandle(job);
    } catch(const exception& e) { cerr<<e.what()<<endl;return 1; }
}