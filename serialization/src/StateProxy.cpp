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
#include "FastStateXml.h"
#include <cstddef>
#include <cstdlib>
#include <cstring>
#include <map>
#include <utility>
#include <vector>

using namespace std;
using namespace OpenMM;

namespace {

bool environmentFlagEnabled(const char* name, bool defaultValue=true) {
    const char* value = std::getenv(name);
    if (value == NULL)
        return defaultValue;
    return !(std::strcmp(value, "0") == 0 ||
             std::strcmp(value, "false") == 0 ||
             std::strcmp(value, "FALSE") == 0 ||
             std::strcmp(value, "off") == 0 ||
             std::strcmp(value, "OFF") == 0);
}

bool shouldUseFastStateXml(const State& state, const SerializationNode& node) {
    if (!environmentFlagEnabled("OPENMM_FAST_XML") ||
            !environmentFlagEnabled("OPENMM_FAST_STATE_WRITE"))
        return false;

    // XmlSerializer::clone() uses an unnamed node and custom-root callers may
    // choose another name.  Only intercept canonical State XML output.
    if (node.getName() != "State")
        return false;

    const int types = state.getDataTypes();
    if ((types & State::IntegratorParameters) != 0)
        return false;

    const size_t minVectors = 100000;
    return ((types & State::Positions) != 0 && state.getPositions().size() >= minVectors) ||
           ((types & State::Velocities) != 0 && state.getVelocities().size() >= minVectors) ||
           ((types & State::Forces) != 0 && state.getForces().size() >= minVectors);
}

void recordArraySize(size_t size, bool& haveArraySize, size_t& arraySize) {
    if (!haveArraySize) {
        arraySize = size;
        haveArraySize = true;
    }
    else if (arraySize != size) {
        throw OpenMMException("State Deserialization Particle Size Mismatch, check number of particles in Forces, Velocities, Positions!");
    }
}

} // anonymous namespace

StateProxy::StateProxy() : SerializationProxy("State") {
}

void StateProxy::serialize(const void* object, SerializationNode& node) const {
    // Discard an abandoned handoff from an earlier serialization on this thread.
    clearFastStateXmlSerialization();

    node.setIntProperty("version", 1);
    node.setStringProperty("openmmVersion", Platform::getOpenMMVersion());
    const State& state = *reinterpret_cast<const State*>(object);
    node.setDoubleProperty("time", state.getTime());
    node.setLongProperty("stepCount", state.getStepCount());

    if (shouldUseFastStateXml(state, node)) {
        registerFastStateXmlSerialization(&state, &node);
        return;
    }

    Vec3 a, b, c;
    state.getPeriodicBoxVectors(a, b, c);
    SerializationNode& boxVectorsNode = node.createChildNode("PeriodicBoxVectors");
    boxVectorsNode.createChildNode("A").setDoubleProperty("x", a[0]).setDoubleProperty("y", a[1]).setDoubleProperty("z", a[2]);
    boxVectorsNode.createChildNode("B").setDoubleProperty("x", b[0]).setDoubleProperty("y", b[1]).setDoubleProperty("z", b[2]);
    boxVectorsNode.createChildNode("C").setDoubleProperty("x", c[0]).setDoubleProperty("y", c[1]).setDoubleProperty("z", c[2]);

    if ((state.getDataTypes() & State::Parameters) != 0) {
        SerializationNode& parametersNode = node.createChildNode("Parameters");
        for (const auto& param : state.getParameters())
            parametersNode.setDoubleProperty(param.first, param.second);
    }
    if ((state.getDataTypes() & State::Energy) != 0) {
        SerializationNode& energiesNode = node.createChildNode("Energies");
        energiesNode.setDoubleProperty("PotentialEnergy", state.getPotentialEnergy());
        energiesNode.setDoubleProperty("KineticEnergy", state.getKineticEnergy());
    }
    if ((state.getDataTypes() & State::Positions) != 0) {
        SerializationNode& positionsNode = node.createChildNode("Positions");
        for (const Vec3& position : state.getPositions())
            positionsNode.createChildNode("Position").setDoubleProperty("x", position[0]).setDoubleProperty("y", position[1]).setDoubleProperty("z", position[2]);
    }
    if ((state.getDataTypes() & State::Velocities) != 0) {
        SerializationNode& velocitiesNode = node.createChildNode("Velocities");
        for (const Vec3& velocity : state.getVelocities())
            velocitiesNode.createChildNode("Velocity").setDoubleProperty("x", velocity[0]).setDoubleProperty("y", velocity[1]).setDoubleProperty("z", velocity[2]);
    }
    if ((state.getDataTypes() & State::Forces) != 0) {
        SerializationNode& forcesNode = node.createChildNode("Forces");
        for (const Vec3& force : state.getForces())
            forcesNode.createChildNode("Force").setDoubleProperty("x", force[0]).setDoubleProperty("y", force[1]).setDoubleProperty("z", force[2]);
    }
    if ((state.getDataTypes() & State::IntegratorParameters) != 0)
        node.getChildren().push_back(state.getIntegratorParameters());
}

void* StateProxy::deserialize(const SerializationNode& node) const {
    if (node.getIntProperty("version") != 1)
        throw OpenMMException("Unsupported version number");

    State::StateBuilder builder(node.getDoubleProperty("time"), node.getLongProperty("stepCount", 0));

    const SerializationNode& boxVectorsNode = node.getChildNode("PeriodicBoxVectors");
    const SerializationNode& aNode = boxVectorsNode.getChildNode("A");
    const SerializationNode& bNode = boxVectorsNode.getChildNode("B");
    const SerializationNode& cNode = boxVectorsNode.getChildNode("C");
    builder.setPeriodicBoxVectors(
        Vec3(aNode.getDoubleProperty("x"), aNode.getDoubleProperty("y"), aNode.getDoubleProperty("z")),
        Vec3(bNode.getDoubleProperty("x"), bNode.getDoubleProperty("y"), bNode.getDoubleProperty("z")),
        Vec3(cNode.getDoubleProperty("x"), cNode.getDoubleProperty("y"), cNode.getDoubleProperty("z")));

    bool haveArraySize = false;
    size_t arraySize = 0;
    for (const auto& child : node.getChildren()) {
        if (child.getName() == "Parameters") {
            map<string, double> parameters;
            for (const auto& param : child.getProperties())
                parameters[param.first] = child.getDoubleProperty(param.first);
            builder.setParameters(std::move(parameters));
        }
        else if (child.getName() == "Energies") {
            builder.setEnergy(child.getDoubleProperty("KineticEnergy"), child.getDoubleProperty("PotentialEnergy"));
        }
        else if (child.getName() == "Positions") {
            vector<Vec3> positions;
            positions.reserve(child.getChildren().size());
            for (const auto& particle : child.getChildren())
                positions.push_back(Vec3(particle.getDoubleProperty("x"), particle.getDoubleProperty("y"), particle.getDoubleProperty("z")));
            recordArraySize(positions.size(), haveArraySize, arraySize);
            builder.setPositions(std::move(positions));
        }
        else if (child.getName() == "Velocities") {
            vector<Vec3> velocities;
            velocities.reserve(child.getChildren().size());
            for (const auto& particle : child.getChildren())
                velocities.push_back(Vec3(particle.getDoubleProperty("x"), particle.getDoubleProperty("y"), particle.getDoubleProperty("z")));
            recordArraySize(velocities.size(), haveArraySize, arraySize);
            builder.setVelocities(std::move(velocities));
        }
        else if (child.getName() == "Forces") {
            vector<Vec3> forces;
            forces.reserve(child.getChildren().size());
            for (const auto& particle : child.getChildren())
                forces.push_back(Vec3(particle.getDoubleProperty("x"), particle.getDoubleProperty("y"), particle.getDoubleProperty("z")));
            recordArraySize(forces.size(), haveArraySize, arraySize);
            builder.setForces(std::move(forces));
        }
        else if (child.getName() == "IntegratorParameters") {
            builder.updateIntegratorParameters() = child;
        }
    }

    return new State(builder.takeState());
}
