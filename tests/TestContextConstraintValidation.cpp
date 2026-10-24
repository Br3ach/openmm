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

#include "OpenMM.h"
#include "openmm/internal/AssertionUtilities.h"
#include <iostream>
#include <string>
using namespace OpenMM;
using namespace std;

static const string selfError = "A constraint cannot connect a particle to itself";
static const string rangeError = "Illegal particle index in constraint";
static const string massError = "A constraint cannot involve a massless particle";
static const string duplicateError = "The System has two constraints between the same atoms.  This will produce a singular constraint matrix.";

static void expectFailure(System& system, const string& expected, bool badProperty=false) {
    VerletIntegrator integrator(0.001);
    map<string, string> properties;
    if (badProperty)
        properties["NotAReferenceProperty"] = "x";
    try {
        Context context(system, integrator, Platform::getPlatformByName("Reference"), properties);
    }
    catch (const OpenMMException& error) {
        ASSERT_EQUAL(expected, string(error.what()));
        return;
    }
    throw runtime_error("Expected Context construction to fail: "+expected);
}

static void testValid() {
    for (int kind=0; kind<3; kind++) {
        System system;
        for (int i=0; i<3; i++)
            system.addParticle(kind==2 ? 0.0 : 1.0);
        if (kind!=0) {
            system.addConstraint(0, 1, 0.1);
            system.addConstraint(0, 2, 0.1);
            system.addConstraint(1, 2, 0.16);
        }
        VerletIntegrator integrator(0.001);
        Context context(system, integrator, Platform::getPlatformByName("Reference"));
        ASSERT_EQUAL(0.0, context.getTime());
        context.reinitialize();
        ASSERT_EQUAL(0.0, context.getTime());
    }
}

static void testErrors() {
    System system;
    system.addParticle(1.0);
    system.addParticle(1.0);
    system.addParticle(0.0);
    system.addConstraint(0, 0, 0.1);
    expectFailure(system, selfError);
    system.setConstraintParameters(0, -1, -1, 0.1);
    expectFailure(system, selfError); // Self check precedes the range check.
    system.setConstraintParameters(0, -1, 0, 0.1);
    expectFailure(system, rangeError);
    system.setConstraintParameters(0, 0, 3, 0.1);
    expectFailure(system, rangeError);
    system.setConstraintParameters(0, 0, 2, 0.1);
    expectFailure(system, massError);
    system.setConstraintParameters(0, 2, 0, 0.1);
    expectFailure(system, massError);
    system.setConstraintParameters(0, 0, 1, 0.1);
    system.addConstraint(0, 1, 0.1);
    expectFailure(system, duplicateError);
    system.setConstraintParameters(1, 1, 0, 0.2);
    expectFailure(system, duplicateError);
    expectFailure(system, duplicateError, true); // Constraints precede properties.
    system.addConstraint(1, 1, 0.1);
    expectFailure(system, duplicateError); // First invalid constraint wins.
    system.setConstraintParameters(1, 1, 1, 0.1);
    expectFailure(system, selfError);
}

int main() {
    try {
        testValid();
        testErrors();
        cout << "Done" << endl;
        return 0;
    }
    catch (const exception& error) {
        cerr << error.what() << endl;
        return 1;
    }
}
