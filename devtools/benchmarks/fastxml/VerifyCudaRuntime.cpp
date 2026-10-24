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
#include "OpenMM.h"
#include "CudaPlatform.h"
#include <cmath>
#include <iostream>
using namespace OpenMM;
using namespace std;
static string quoted(const string& text) {
    string out="\"";
    for(char c:text){if(c=='\\'||c=='\"')out+='\\';out+=c;}
    return out+"\"";
}
int main() {
    try {
        System system;system.addParticle(16);system.addParticle(1);system.addParticle(1);system.addParticle(0);
        system.addConstraint(0,1,0.1);system.addConstraint(0,2,0.1);system.addConstraint(1,2,0.16);
        system.setVirtualSite(3,new ThreeParticleAverageSite(0,1,2,0.6,0.2,0.2));
        system.setDefaultPeriodicBoxVectors(Vec3(2,0,0),Vec3(0,2,0),Vec3(0,0,2));
        auto* force=new NonbondedForce();force->setNonbondedMethod(NonbondedForce::PME);force->setCutoffDistance(0.7);
        for(double charge:{-0.8,0.4,0.4,0.0})force->addParticle(charge,0.3,0.0);
        force->addException(0,1,0,1,0);force->addException(0,2,0,1,0);force->addException(1,2,0,1,0);system.addForce(force);
        VerletIntegrator integrator(0.001);CudaPlatform platform;map<string,string> properties;
        properties["Precision"]="mixed";properties["UseBlockingSync"]="true";
        Context context(system,integrator,platform,properties);
        context.setPositions({Vec3(0.5,0.5,0.5),Vec3(0.6,0.5,0.5),Vec3(0.472,0.596,0.5),Vec3()});context.computeVirtualSites();
        double energy=context.getState(State::Energy).getPotentialEnergy();if(!isfinite(energy))throw runtime_error("PME energy not finite");
        cout<<"{\"pme_energy_finite\":true,\"loaded_modules\":{";int count=0;
        for(const char* name:{"OpenMM.dll","OpenMMCUDA.dll","cufft64_12.dll","nvrtc64_130_0.dll","nvrtc-builtins64_134.dll"}) {
            char path[32768]={};HMODULE module=GetModuleHandleA(name);
            if(!module||!GetModuleFileNameA(module,path,sizeof(path)))throw runtime_error(string("Module missing: ")+name);
            if(count++)cout<<',';cout<<quoted(name)<<':'<<quoted(path);
        }
        cout<<"}}"<<endl;return 0;
    }catch(const exception& error){cerr<<error.what()<<endl;return 1;}
}
