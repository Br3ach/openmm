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

#include "OpenMM.h"
#include "CudaPlatform.h"
#include <cmath>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
using namespace OpenMM;
using namespace std;

static uint64_t hashWord(uint64_t h, uint64_t word) {
    for (int i=0; i<8; i++)
        h=(h^((word>>(8*i))&255))*1099511628211ULL;
    return h;
}
static uint64_t hashDouble(uint64_t h, double value) {
    uint64_t word;memcpy(&word,&value,sizeof(word));return hashWord(h,word);
}
static uint64_t snapshot(Context& context) {
    State state=context.getState(State::Positions|State::Velocities|State::Forces|State::Energy|State::Parameters);
    uint64_t h=14695981039346656037ULL;
    for (const auto* values : {&state.getPositions(), &state.getVelocities(), &state.getForces()})
        for (const auto& value:*values)
            for (int d=0; d<3; d++)
                h=hashDouble(h,value[d]);
    h=hashDouble(hashDouble(hashDouble(h,state.getTime()),state.getPotentialEnergy()),state.getKineticEnergy());
    h=hashWord(h,state.getStepCount());
    for (const auto& parameter:state.getParameters()) {
        for (unsigned char c:parameter.first)h=hashWord(h,c);
        h=hashDouble(h,parameter.second);
    }
    const System& system=context.getSystem();
    for (int i=0; i<system.getNumConstraints(); i++) {
        int a,b;double distance;system.getConstraintParameters(i,a,b,distance);
        Vec3 delta=state.getPositions()[a]-state.getPositions()[b];
        if (abs(sqrt(delta.dot(delta))-distance)>1e-5)
            throw runtime_error("Constraint distance changed");
    }
    return h;
}
static void makeSystem(System& system, vector<Vec3>& positions) {
    // Sorted neighbor order must choose the same central atom at all three IDs.
    for (int group=0; group<6; group++) {
        int first=system.getNumParticles(), central=group%3;
        Vec3 origin(group*0.6,0,0);
        for (int j=0; j<3; j++) {
            system.addParticle(j==central ? 16.0 : 1.0);
            positions.push_back(origin);
        }
        int outer1=(central+1)%3, outer2=(central+2)%3;
        positions[first+outer1]=origin+Vec3(0.1,0,0);
        positions[first+outer2]=origin+Vec3(-0.028,0.096,0);
        for (int edge=2; edge>=0; edge--) {
            int a=edge,b=(edge+1)%3;
            double distance=(a==central || b==central ? 0.1 : 0.16);
            if (group==3 && distance==0.1)distance+=1e-11; // Equal after float rounding.
            system.addConstraint(first+b,first+a,distance);
        }
        system.addParticle(0.0);positions.push_back(origin);
        system.setVirtualSite(first+3,new ThreeParticleAverageSite(first,first+1,first+2,0.6,0.2,0.2));
    }
    // A nonisosceles triangle must continue down the CCMA path.
    int first=system.getNumParticles();
    for (int i=0;i<3;i++)system.addParticle(12.0);
    double x=0.0935,y=sqrt(0.16*0.16-x*x);
    positions.push_back(Vec3(0,1,0));positions.push_back(Vec3(0.1,1,0));positions.push_back(Vec3(x,1+y,0));
    system.addConstraint(first,first+1,0.1);system.addConstraint(first+1,first+2,0.13);system.addConstraint(first+2,first,0.16);
    // A degree2 cycle is not a SETTLE triangle.
    first=system.getNumParticles();
    for (int i=0;i<4;i++)system.addParticle(12.0);
    for (Vec3 p:{Vec3(0,2,0),Vec3(0.12,2,0),Vec3(0.12,2.12,0),Vec3(0,2.12,0)})positions.push_back(p);
    for (int i=0;i<4;i++)system.addConstraint(first+i,first+(i+1)%4,0.12);
    // Valid SHAKE star, plus an open chain requiring CCMA.
    first=system.getNumParticles();system.addParticle(12.0);positions.push_back(Vec3(0,3,0));
    for (int i=0;i<3;i++) {
        system.addParticle(1.0);Vec3 p(0,3,0);p[i]+=0.11;positions.push_back(p);system.addConstraint(first,first+i+1,0.11);
    }
    first=system.getNumParticles();
    for (int i=0;i<4;i++) {system.addParticle(12.0);positions.push_back(Vec3(0.15*i,4,0));if(i>0)system.addConstraint(first+i-1,first+i,0.15);}
    auto* force=new CustomExternalForce("0.005*(x*x+y*y+z*z)");
    for (int i=0;i<system.getNumParticles();i++)force->addParticle(i,vector<double>());
    system.addForce(force);
}
int main(int argc,char** argv) {
    try {
        string precision=argc>1?argv[1]:"single",mode=argc>2?argv[2]:"fresh";
        if(mode!="fresh"&&mode!="save"&&mode!="load")throw runtime_error("Invalid mode");
        if(mode!="fresh"&&argc!=4)throw runtime_error("Checkpoint path required");
        System system;vector<Vec3> positions;makeSystem(system,positions);
        LangevinMiddleIntegrator integrator(300.0,1.0,0.0005);integrator.setRandomNumberSeed(12345);integrator.setConstraintTolerance(1e-6);
        CudaPlatform platform;map<string,string> properties;
        properties["Precision"]=precision;properties["DeterministicForces"]="true";properties["UseBlockingSync"]="true";
        Context context(system,integrator,platform,properties);
        if(mode=="load") {
            ifstream checkpoint(argv[3],ios::binary);if(!checkpoint)throw runtime_error("Cannot open checkpoint");context.loadCheckpoint(checkpoint);
        } else {
            context.setPositions(positions);context.computeVirtualSites();context.setVelocitiesToTemperature(300.0,20261024);context.applyVelocityConstraints(1e-6);integrator.step(16);
        }
        stringstream saved(ios::in|ios::out|ios::binary);context.createCheckpoint(saved);
        if(mode=="save") {ofstream file(argv[3],ios::binary);file<<saved.rdbuf();if(!file)throw runtime_error("Cannot write checkpoint");}
        integrator.step(32);uint64_t first=snapshot(context);
        if(mode!="load") {
            saved.clear();saved.seekg(0);context.loadCheckpoint(saved);integrator.step(32);
            if(snapshot(context)!=first)throw runtime_error("Checkpoint replay differs");
        }
        cout<<"{\"precision\":\""<<precision<<"\",\"particles\":"<<system.getNumParticles()<<",\"constraints\":"<<system.getNumConstraints()<<",\"fingerprint\":\""<<hex<<first<<"\",\"mode\":\""<<mode<<"\",\"self_replay_verified\":"<<(mode=="load"?"false":"true")<<"}"<<endl;
        return 0;
    } catch(const exception& error) {cerr<<error.what()<<endl;return 1;}
}
