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

#include "OpenMM.h"
#include "openmm/internal/AssertionUtilities.h"
#include <sstream>
#include <iostream>
using namespace OpenMM;
using namespace std;

// A State must remain an independent snapshot after later getState() calls,
// Context mutation/destruction, and integrator-parameter mutation.
static void testSnapshots() {
    State saved;
    string originalXml;
    {
        System system;
        system.addParticle(1); system.addParticle(2);
        system.setDefaultPeriodicBoxVectors(Vec3(3,0,0),Vec3(0,4,0),Vec3(0,0,5));
        auto* bond=new CustomBondForce("k*(r-r0)^2");
        bond->addGlobalParameter("k",10);bond->addGlobalParameter("r0",0.2);
        bond->addEnergyParameterDerivative("k");bond->addBond(0,1);
        bond->setForceGroup(1);system.addForce(bond);
        CustomIntegrator integrator(0.001);
        integrator.addGlobalVariable("counter",3);
        integrator.addPerDofVariable("history",1);
        Context context(system,integrator,Platform::getPlatformByName("Reference"));
        vector<Vec3> positions={Vec3(4.1,1,1),Vec3(4.2,1,1)};
        vector<Vec3> velocities={Vec3(0.1,0.2,0.3),Vec3(-0.1,-0.2,-0.3)};
        context.setPositions(positions);context.setVelocities(velocities);
        context.setTime(1.25);context.setStepCount(123);
        const int all=State::Positions|State::Velocities|State::Forces|State::Energy|
                      State::Parameters|State::ParameterDerivatives|State::IntegratorParameters;
        saved=context.getState(all,true,2);
        ASSERT_EQUAL(all,saved.getDataTypes());
        ASSERT_EQUAL_VEC(Vec3(1.1,1,1),saved.getPositions()[0],1e-12);
        ASSERT_EQUAL_VEC(velocities[1],saved.getVelocities()[1],1e-12);
        ASSERT_EQUAL_TOL(0.1,saved.getPotentialEnergy(),1e-12);
        ASSERT_EQUAL_TOL(0.01,saved.getEnergyParameterDerivatives().at("k"),1e-12);
        ASSERT_EQUAL(10,saved.getParameters().at("k"));
        ASSERT_EQUAL(123,saved.getStepCount());ASSERT_EQUAL(1.25,saved.getTime());
        stringstream first;XmlSerializer::serialize(&saved,"State",first);originalXml=first.str();
        for(int mask=0;mask<=all;mask++) {
            State next=context.getState(mask,false,2);
            ASSERT_EQUAL(mask,next.getDataTypes());
            if(mask&State::Positions)ASSERT_EQUAL_VEC(positions[0],next.getPositions()[0],1e-12);
            if(mask&State::Velocities)ASSERT_EQUAL_VEC(velocities[0],next.getVelocities()[0],1e-12);
            if(mask&State::Forces)ASSERT_EQUAL_VEC(saved.getForces()[0],next.getForces()[0],1e-12);
            if(mask&State::Energy)ASSERT_EQUAL_TOL(saved.getPotentialEnergy(),next.getPotentialEnergy(),1e-12);
            if(mask&State::Parameters)ASSERT_EQUAL_CONTAINERS(saved.getParameters(),next.getParameters());
            if(mask&State::ParameterDerivatives)ASSERT_EQUAL_CONTAINERS(saved.getEnergyParameterDerivatives(),next.getEnergyParameterDerivatives());
            if(!(mask&State::Positions)) {
                bool threw=false;try{next.getPositions();}catch(const OpenMMException&){threw=true;}
                ASSERT(threw);
            }
        }
        // Wrapped output must not change the Context's coordinates.
        ASSERT_EQUAL_VEC(positions[0],context.getState(State::Positions).getPositions()[0],1e-12);
        State excluded=context.getState(State::Forces|State::Energy,false,1);
        ASSERT_EQUAL_TOL(0,excluded.getPotentialEnergy(),1e-12);
        ASSERT_EQUAL_VEC(Vec3(),excluded.getForces()[0],1e-12);
        context.setPositions({Vec3(1,2,3),Vec3(2,3,4)});
        context.setVelocities({Vec3(3,2,1),Vec3(4,3,2)});
        context.setParameter("k",20);context.setTime(2.5);context.setStepCount(456);
        integrator.setGlobalVariable(0,9);
        integrator.setPerDofVariable(0,{Vec3(9,8,7),Vec3(6,5,4)});
        context.setPeriodicBoxVectors(Vec3(6,0,0),Vec3(0,7,0),Vec3(0,0,8));
        State later=context.getState(all);
        ASSERT_EQUAL(20,later.getParameters().at("k"));
        stringstream second;XmlSerializer::serialize(&saved,"State",second);
        ASSERT_EQUAL(originalXml,second.str());
    }
    stringstream afterDestruction;XmlSerializer::serialize(&saved,"State",afterDestruction);
    ASSERT_EQUAL(originalXml,afterDestruction.str());
    State copy=saved;
    stringstream copyXml;XmlSerializer::serialize(&copy,"State",copyXml);
    ASSERT_EQUAL(originalXml,copyXml.str());
}
int main(){try{testSnapshots();cout<<"Snapshot ownership and all 128 data masks passed"<<endl;return 0;}
catch(const exception& e){cerr<<e.what()<<endl;return 1;}}
