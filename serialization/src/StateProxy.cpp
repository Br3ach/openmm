/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2010-2021 Stanford University and the Authors.      *
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

#include "openmm/serialization/StateProxy.h"
#include "openmm/Platform.h"
#include "openmm/State.h"
#include "openmm/Vec3.h"
#include <cstddef>
#include <cstdlib>
#include <cstring>
#include <map>

using namespace std;
using namespace OpenMM;

// Implemented in XmlSerializer.cpp.  The normal XmlSerializer<T> template in
// the already-built caller invokes StateProxy::serialize() first and then calls
// XmlSerializer::serialize(node, stream).  For very large States we defer the
// array payload here and let XmlSerializer.cpp stream it directly to XML,
// avoiding millions of SerializationNode objects without changing the caller.
namespace OpenMM {
void clearFastStateXmlSerialization();
void registerFastStateXmlSerialization(const State* state, SerializationNode* node);
}

static bool fastXmlEnabledByEnvironment() {
    const char* value = std::getenv("OPENMM_FAST_XML");
    if (value == NULL)
        return true;
    return !(std::strcmp(value, "0") == 0 ||
             std::strcmp(value, "false") == 0 ||
             std::strcmp(value, "FALSE") == 0 ||
             std::strcmp(value, "off") == 0 ||
             std::strcmp(value, "OFF") == 0);
}

static bool fastStateWriteEnabledByEnvironment() {
    if (!fastXmlEnabledByEnvironment())
        return false;
    const char* value = std::getenv("OPENMM_FAST_STATE_WRITE");
    if (value == NULL)
        return true; // v7: fast writer is enabled by default.
    return !(std::strcmp(value, "0") == 0 ||
             std::strcmp(value, "false") == 0 ||
             std::strcmp(value, "FALSE") == 0 ||
             std::strcmp(value, "off") == 0 ||
             std::strcmp(value, "OFF") == 0);
}

static bool shouldUseFastStateXml(const State& state, const SerializationNode& node) {
    if (!fastStateWriteEnabledByEnvironment())
        return false;
    // XmlSerializer::clone() invokes the proxy with an unnamed node and then
    // immediately deserializes that node without passing through the XML
    // writer.  Likewise, callers may choose an arbitrary XML root name.
    // Restrict this DLL-only handoff to the canonical State XML serialization
    // used by existing Core 28 callers.  In particular, clone() and custom-root
    // serializations retain the normal complete SerializationNode path.
    if (node.getName() != "State")
        return false;

    const int types = state.getDataTypes();

    // IntegratorParameters are themselves a SerializationNode tree.  Keep the
    // stock path for those States so all existing semantics are preserved.
    if ((types & State::IntegratorParameters) != 0)
        return false;

    // Restrict the workaround to genuinely large canonical State XML documents.
    // This keeps ordinary OpenMM serializations on the untouched stock path.
    const size_t minVectors = 100000;
    if ((types & State::Positions) != 0 && state.getPositions().size() >= minVectors)
        return true;
    if ((types & State::Velocities) != 0 && state.getVelocities().size() >= minVectors)
        return true;
    if ((types & State::Forces) != 0 && state.getForces().size() >= minVectors)
        return true;
    return false;
}

StateProxy::StateProxy() : SerializationProxy("State") {

}

void StateProxy::serialize(const void* object, SerializationNode& node) const {
    // Discard any abandoned handoff from an earlier State serialization on
    // this thread before registering a new one.
    clearFastStateXmlSerialization();

    node.setIntProperty("version", 1);
    node.setStringProperty("openmmVersion", Platform::getOpenMMVersion());
    const State& s = *reinterpret_cast<const State*>(object);
    node.setDoubleProperty("time", s.getTime());
    node.setLongProperty("stepCount", s.getStepCount());

    if (shouldUseFastStateXml(s, node)) {
        // Keep only the small root metadata in the node.  The caller will add
        // type="State" exactly as before, then XmlSerializer.cpp consumes this
        // deferred State and emits the normal XML schema directly.
        registerFastStateXmlSerialization(&s, &node);
        return;
    }

    Vec3 a,b,c;
    s.getPeriodicBoxVectors(a,b,c);
    SerializationNode& boxVectorsNode = node.createChildNode("PeriodicBoxVectors");
    boxVectorsNode.createChildNode("A").setDoubleProperty("x", a[0]).setDoubleProperty("y", a[1]).setDoubleProperty("z", a[2]);
    boxVectorsNode.createChildNode("B").setDoubleProperty("x", b[0]).setDoubleProperty("y", b[1]).setDoubleProperty("z", b[2]);
    boxVectorsNode.createChildNode("C").setDoubleProperty("x", c[0]).setDoubleProperty("y", c[1]).setDoubleProperty("z", c[2]);
    if ((s.getDataTypes()&State::Parameters) != 0) {
        s.getParameters();
        SerializationNode& parametersNode = node.createChildNode("Parameters");
        for (auto& param : s.getParameters())
            parametersNode.setDoubleProperty(param.first, param.second);
    }
    if ((s.getDataTypes()&State::Energy) != 0) {
        s.getPotentialEnergy();
        SerializationNode& energiesNode = node.createChildNode("Energies");
        energiesNode.setDoubleProperty("PotentialEnergy", s.getPotentialEnergy());
        energiesNode.setDoubleProperty("KineticEnergy", s.getKineticEnergy());
    }
    if ((s.getDataTypes()&State::Positions) != 0) {
        s.getPositions();
        SerializationNode& positionsNode = node.createChildNode("Positions");
        vector<Vec3> statePositions = s.getPositions();
        for (int i=0; i<statePositions.size();i++) {
           positionsNode.createChildNode("Position").setDoubleProperty("x", statePositions[i][0]).setDoubleProperty("y", statePositions[i][1]).setDoubleProperty("z", statePositions[i][2]);
        }
    }
    if ((s.getDataTypes()&State::Velocities) != 0) {
        s.getVelocities();
        SerializationNode& velocitiesNode = node.createChildNode("Velocities");
        vector<Vec3> stateVelocities = s.getVelocities();
        for (int i=0; i<stateVelocities.size();i++) {
           velocitiesNode.createChildNode("Velocity").setDoubleProperty("x", stateVelocities[i][0]).setDoubleProperty("y", stateVelocities[i][1]).setDoubleProperty("z", stateVelocities[i][2]);
        }
    }
    if ((s.getDataTypes()&State::Forces) != 0) {
        s.getForces();
        SerializationNode& forcesNode = node.createChildNode("Forces");
        vector<Vec3> stateForces = s.getForces();
        for (int i=0; i<stateForces.size();i++) {
            forcesNode.createChildNode("Force").setDoubleProperty("x", stateForces[i][0]).setDoubleProperty("y", stateForces[i][1]).setDoubleProperty("z", stateForces[i][2]);
        }
    }
    if ((s.getDataTypes()&State::IntegratorParameters) != 0) {
        node.getChildren().push_back(s.getIntegratorParameters());
    }
}

void* StateProxy::deserialize(const SerializationNode& node) const {
    if (node.getIntProperty("version") != 1)
        throw OpenMMException("Unsupported version number");
    double outTime = node.getDoubleProperty("time");
    long long outStepCount = node.getLongProperty("stepCount", 0);
    const SerializationNode& boxVectorsNode = node.getChildNode("PeriodicBoxVectors");
    const SerializationNode& AVec = boxVectorsNode.getChildNode("A");
    Vec3 outAVec(AVec.getDoubleProperty("x"),AVec.getDoubleProperty("y"),AVec.getDoubleProperty("z"));
    const SerializationNode& BVec = boxVectorsNode.getChildNode("B");
    Vec3 outBVec(BVec.getDoubleProperty("x"),BVec.getDoubleProperty("y"),BVec.getDoubleProperty("z"));
    const SerializationNode& CVec = boxVectorsNode.getChildNode("C");
    Vec3 outCVec(CVec.getDoubleProperty("x"),CVec.getDoubleProperty("y"),CVec.getDoubleProperty("z"));
    int types = 0;
    vector<int> arraySizes;
    State::StateBuilder builder(outTime, outStepCount);
    for (auto& child : node.getChildren()) {
        if (child.getName() == "Parameters") {
            map<string, double> outStateParams;
            for (auto& param : child.getProperties())
                outStateParams[param.first] = child.getDoubleProperty(param.first);
            builder.setParameters(outStateParams);
        }
        else if (child.getName() == "Energies") {
            double potentialEnergy = child.getDoubleProperty("PotentialEnergy");
            double kineticEnergy = child.getDoubleProperty("KineticEnergy");
            builder.setEnergy(kineticEnergy, potentialEnergy);
        }
        else if (child.getName() == "Positions") {
            vector<Vec3> outPositions;
            for (auto& particle : child.getChildren())
                outPositions.push_back(Vec3(particle.getDoubleProperty("x"),particle.getDoubleProperty("y"),particle.getDoubleProperty("z")));
            builder.setPositions(outPositions);
            arraySizes.push_back(outPositions.size());
        }
        else if (child.getName() == "Velocities") {
            vector<Vec3> outVelocities;
            for (auto& particle : child.getChildren())
                outVelocities.push_back(Vec3(particle.getDoubleProperty("x"),particle.getDoubleProperty("y"),particle.getDoubleProperty("z")));
            builder.setVelocities(outVelocities);
            arraySizes.push_back(outVelocities.size());
        }
        else if (child.getName() == "Forces") {
            vector<Vec3> outForces;
            for (auto& particle : child.getChildren())
                outForces.push_back(Vec3(particle.getDoubleProperty("x"),particle.getDoubleProperty("y"),particle.getDoubleProperty("z")));
            builder.setForces(outForces);
            arraySizes.push_back(outForces.size());
        }
        else if (child.getName() == "IntegratorParameters") {
            builder.updateIntegratorParameters() = child;
        }
    }
    for (int i = 1; i < arraySizes.size(); i++) {
        if (arraySizes[i] != arraySizes[i-1]) {
            throw(OpenMMException("State Deserialization Particle Size Mismatch, check number of particles in Forces, Velocities, Positions!"));
        }
    }
    builder.setPeriodicBoxVectors(outAVec, outBVec, outCVec);
    State *s = new State();
    *s = builder.getState();
    return s;
}
