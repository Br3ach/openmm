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
#include <chrono>
#include <cstdint>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <memory>
#include <stdexcept>
using namespace OpenMM;
using namespace std;

static PROCESS_MEMORY_COUNTERS_EX memory() {
    PROCESS_MEMORY_COUNTERS_EX result={};
    result.cb=sizeof(result);
    if (!GetProcessMemoryInfo(GetCurrentProcess(), (PROCESS_MEMORY_COUNTERS*) &result, sizeof(result)))
        throw runtime_error("GetProcessMemoryInfo failed");
    return result;
}
static void printMemory(const char* phase, const PROCESS_MEMORY_COUNTERS_EX& m) {
    cout << ",\"" << phase << "_working_set_kib\":" << m.WorkingSetSize/1024
         << ",\"" << phase << "_peak_working_set_kib\":" << m.PeakWorkingSetSize/1024
         << ",\"" << phase << "_private_kib\":" << m.PrivateUsage/1024;
}
static uint64_t hashInt(uint64_t h, int value) {
    for (int i=0; i<4; i++)
        h=(h^((value>>(8*i))&255))*1099511628211ULL;
    return h;
}
int main(int argc, char** argv) {
    try {
        if (argc!=3)
            throw runtime_error("Usage: BenchmarkContextConstruction particles constraints");
        const int count=atoi(argv[1]), constraints=atoi(argv[2]);
        if (count<3 || count>10000000 || constraints<0 || constraints>3*(count/3))
            throw runtime_error("Invalid particle or constraint count");
        HANDLE job=CreateJobObjectW(nullptr, nullptr);
        JOBOBJECT_EXTENDED_LIMIT_INFORMATION limit={};
        limit.BasicLimitInformation.LimitFlags=JOB_OBJECT_LIMIT_PROCESS_MEMORY;
        limit.ProcessMemoryLimit=6ULL*1024*1024*1024;
        if (!job || !SetInformationJobObject(job, JobObjectExtendedLimitInformation, &limit, sizeof(limit)) ||
                !AssignProcessToJobObject(job, GetCurrentProcess()))
            throw runtime_error("Could not establish 6 GiB process memory limit");
        System system;
        for (int i=0; i<count; i++)
            system.addParticle(i%3==0 ? 16.0 : 1.0);
        uint64_t fingerprint=hashInt(hashInt(14695981039346656037ULL, count), constraints);
        for (int i=0; i<constraints; i++) {
            int base=3*(i/3), edge=i%3;
            int a=base+(edge==2 ? 1 : 0), b=base+(edge==0 ? 1 : 2);
            system.addConstraint(a, b, edge==2 ? 0.16 : 0.1);
            fingerprint=hashInt(hashInt(fingerprint, a), b);
        }
        VerletIntegrator integrator(0.001);
        Platform& platform=Platform::getPlatformByName("Reference");
        auto before=memory();
        auto start=chrono::steady_clock::now();
        unique_ptr<Context> context(new Context(system, integrator, platform));
        double constructionSeconds=chrono::duration<double>(chrono::steady_clock::now()-start).count();
        auto constructed=memory();
        if (context->getTime()!=0.0 || context->getSystem().getNumConstraints()!=constraints)
            throw runtime_error("Unexpected constructed Context");
        start=chrono::steady_clock::now();
        context.reset();
        double destructionSeconds=chrono::duration<double>(chrono::steady_clock::now()-start).count();
        auto destroyed=memory();
        cout << setprecision(10) << "{\"particles\":" << count << ",\"constraints\":" << constraints
             << ",\"construction_seconds\":" << constructionSeconds << ",\"destruction_seconds\":" << destructionSeconds;
        printMemory("before", before);
        printMemory("constructed", constructed);
        printMemory("destroyed", destroyed);
        cout << ",\"input_fingerprint\":\"" << hex << fingerprint
             << "\",\"platform\":\"Reference\",\"peak_scope\":\"process_lifetime\"}" << endl;
        CloseHandle(job);
        return 0;
    }
    catch (const exception& error) {
        cerr << error.what() << endl;
        return 1;
    }
}
