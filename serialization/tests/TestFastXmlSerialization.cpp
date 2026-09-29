/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * Regression tests for the streaming State/System XML fast paths.            *
 * -------------------------------------------------------------------------- */

#include "openmm/internal/AssertionUtilities.h"
#include "openmm/CMMotionRemover.h"
#include "openmm/CustomBondForce.h"
#include "openmm/HarmonicAngleForce.h"
#include "openmm/HarmonicBondForce.h"
#include "openmm/MonteCarloBarostat.h"
#include "openmm/NonbondedForce.h"
#include "openmm/OpenMMException.h"
#include "openmm/PeriodicTorsionForce.h"
#include "openmm/State.h"
#include "openmm/System.h"
#include "openmm/serialization/SerializationNode.h"
#include "openmm/serialization/XmlSerializer.h"
#include <cstdlib>
#include <iostream>
#include <map>
#include <memory>
#include <sstream>
#include <string>
#include <streambuf>
#include <vector>

using namespace OpenMM;
using namespace std;

namespace {

static void setEnvironment(const char* name, const char* value) {
#ifdef _WIN32
    _putenv_s(name, value == NULL ? "" : value);
#else
    if (value == NULL)
        unsetenv(name);
    else
        setenv(name, value, 1);
#endif
}

class EnvironmentGuard {
public:
    explicit EnvironmentGuard(const char* name) : name(name), hadValue(false) {
        const char* value = getenv(name);
        if (value != NULL) {
            hadValue = true;
            oldValue = value;
        }
    }
    ~EnvironmentGuard() {
        setEnvironment(name.c_str(), hadValue ? oldValue.c_str() : NULL);
    }
private:
    string name;
    string oldValue;
    bool hadValue;
};

// A seekable stream buffer that can be told to fail a specific seekg(pos)
// operation.  This makes fast-path fallback/rewind behavior observable without
// depending on the stock parser's treatment of deliberately non-canonical XML.
class FailOnSeekStringBuf : public stringbuf {
public:
    FailOnSeekStringBuf(const string& value, int failOnSeek) :
            stringbuf(value, ios_base::in), failOnSeek(failOnSeek), seekCount(0) {
    }

protected:
    pos_type seekpos(pos_type position, ios_base::openmode mode = ios_base::in | ios_base::out) override {
        ++seekCount;
        if (seekCount == failOnSeek)
            return pos_type(off_type(-1));
        return stringbuf::seekpos(position, mode);
    }

private:
    int failOnSeek;
    int seekCount;
};

static void setFastXml(bool enabled) {
    setEnvironment("OPENMM_FAST_XML", enabled ? "1" : "0");
}

static void setFastStateWrite(bool enabled) {
    setEnvironment("OPENMM_FAST_STATE_WRITE", enabled ? "1" : "0");
}

template <class T>
static string serializeObject(const T& object, const char* rootName, bool fastXml) {
    setFastXml(fastXml);
    stringstream stream;
    XmlSerializer::serialize<T>(&object, rootName, stream);
    return stream.str();
}

template <class T>
static unique_ptr<T> deserializeObject(const string& xml, bool fastXml) {
    setFastXml(fastXml);
    stringstream stream(xml);
    return unique_ptr<T>(XmlSerializer::deserialize<T>(stream));
}

template <class F>
static void assertThrowsOpenMM(F operation) {
    bool threw = false;
    try {
        operation();
    }
    catch (const OpenMMException&) {
        threw = true;
    }
    ASSERT(threw);
}

static void compareStates(const State& expected, const State& found) {
    ASSERT_EQUAL(expected.getDataTypes(), found.getDataTypes());
    ASSERT_EQUAL(expected.getTime(), found.getTime());
    ASSERT_EQUAL(expected.getStepCount(), found.getStepCount());

    Vec3 ea, eb, ec, fa, fb, fc;
    expected.getPeriodicBoxVectors(ea, eb, ec);
    found.getPeriodicBoxVectors(fa, fb, fc);
    ASSERT_EQUAL_VEC(ea, fa, 0);
    ASSERT_EQUAL_VEC(eb, fb, 0);
    ASSERT_EQUAL_VEC(ec, fc, 0);

    if (expected.getDataTypes() & State::Positions) {
        const vector<Vec3>& a = expected.getPositions();
        const vector<Vec3>& b = found.getPositions();
        ASSERT_EQUAL(a.size(), b.size());
        for (size_t i = 0; i < a.size(); ++i)
            ASSERT_EQUAL_VEC(a[i], b[i], 0);
    }
    if (expected.getDataTypes() & State::Velocities) {
        const vector<Vec3>& a = expected.getVelocities();
        const vector<Vec3>& b = found.getVelocities();
        ASSERT_EQUAL(a.size(), b.size());
        for (size_t i = 0; i < a.size(); ++i)
            ASSERT_EQUAL_VEC(a[i], b[i], 0);
    }
    if (expected.getDataTypes() & State::Forces) {
        const vector<Vec3>& a = expected.getForces();
        const vector<Vec3>& b = found.getForces();
        ASSERT_EQUAL(a.size(), b.size());
        for (size_t i = 0; i < a.size(); ++i)
            ASSERT_EQUAL_VEC(a[i], b[i], 0);
    }
    if (expected.getDataTypes() & State::Energy) {
        ASSERT_EQUAL(expected.getKineticEnergy(), found.getKineticEnergy());
        ASSERT_EQUAL(expected.getPotentialEnergy(), found.getPotentialEnergy());
    }
    if (expected.getDataTypes() & State::Parameters)
        ASSERT_EQUAL_CONTAINERS(expected.getParameters(), found.getParameters());
}

static State createLargeState() {
    const int count = 100000; // Activates the production direct-State writer threshold.
    vector<Vec3> positions;
    positions.reserve(count);
    for (int i = 0; i < count; ++i)
        positions.push_back(Vec3(0.001*i, -0.002*i, 0.003*i));

    map<string, double> parameters;
    parameters["alpha"] = 1.25;
    parameters["beta"] = -3.5;

    State::StateBuilder builder(2.5, 1234567890123LL);
    builder.setPeriodicBoxVectors(Vec3(4.0, 0.1, 0.2), Vec3(0.0, 5.0, 0.3), Vec3(0.0, 0.0, 6.0));
    builder.setPositions(std::move(positions));
    builder.setParameters(std::move(parameters));
    builder.setEnergy(7.25, -19.5);
    return builder.takeState();
}

static void testLargeStateFastWriterAndReaderParity() {
    State state = createLargeState();
    setFastStateWrite(true);

    // The direct writer is intended to be schema-compatible with the stock
    // SerializationNode writer, including deterministic attribute ordering.
    const string stockXml = serializeObject(state, "State", false);
    const string fastXml = serializeObject(state, "State", true);
    ASSERT_EQUAL(stockXml, fastXml);

    unique_ptr<State> stockCopy = deserializeObject<State>(stockXml, false);
    unique_ptr<State> fastCopy = deserializeObject<State>(stockXml, true);
    compareStates(*stockCopy, *fastCopy);
    compareStates(state, *fastCopy);
}

static unique_ptr<System> createSupportedSystem() {
    unique_ptr<System> system(new System());
    system->setDefaultPeriodicBoxVectors(Vec3(3.0, 0.0, 0.0), Vec3(0.1, 4.0, 0.0), Vec3(0.2, 0.3, 5.0));
    for (int i = 0; i < 4; ++i)
        system->addParticle(10.0+i);
    system->addConstraint(0, 1, 0.15);

    NonbondedForce* nonbonded = new NonbondedForce();
    nonbonded->setNonbondedMethod(NonbondedForce::CutoffPeriodic);
    nonbonded->setCutoffDistance(1.1);
    nonbonded->setUseSwitchingFunction(true);
    nonbonded->setSwitchingDistance(0.9);
    nonbonded->setEwaldErrorTolerance(1e-5);
    nonbonded->setReactionFieldDielectric(70.0);
    nonbonded->setUseDispersionCorrection(false);
    nonbonded->setIncludeDirectSpace(true);
    nonbonded->setReciprocalSpaceForceGroup(3);
    nonbonded->setExceptionsUsePeriodicBoundaryConditions(true);
    const int global = nonbonded->addGlobalParameter("lambda", 0.75);
    (void) global;
    for (int i = 0; i < 4; ++i)
        nonbonded->addParticle(0.1*i, 0.2+0.01*i, 0.3+0.02*i);
    const int exception = nonbonded->addException(0, 1, 0.0, 0.25, 0.0);
    nonbonded->addParticleParameterOffset("lambda", 2, 0.01, 0.02, 0.03);
    nonbonded->addExceptionParameterOffset("lambda", exception, 0.04, 0.05, 0.06);
    system->addForce(nonbonded);

    HarmonicBondForce* bonds = new HarmonicBondForce();
    bonds->setUsesPeriodicBoundaryConditions(true);
    bonds->addBond(0, 1, 0.15, 250.0);
    system->addForce(bonds);

    HarmonicAngleForce* angles = new HarmonicAngleForce();
    angles->setUsesPeriodicBoundaryConditions(true);
    angles->addAngle(0, 1, 2, 1.75, 80.0);
    system->addForce(angles);

    PeriodicTorsionForce* torsions = new PeriodicTorsionForce();
    torsions->setUsesPeriodicBoundaryConditions(true);
    torsions->addTorsion(0, 1, 2, 3, 3, 0.25, 5.0);
    system->addForce(torsions);

    system->addForce(new CMMotionRemover(7));
    system->addForce(new MonteCarloBarostat(1.5, 310.0, 19));
    return system;
}

static void testSupportedSystemFastReaderParity() {
    unique_ptr<System> system = createSupportedSystem();
    const string xml = serializeObject(*system, "System", false);

    unique_ptr<System> stockCopy = deserializeObject<System>(xml, false);
    unique_ptr<System> fastCopy = deserializeObject<System>(xml, true);

    // Re-serialize both parsed objects through the stock writer.  This compares
    // all fields represented in the canonical System/Force XML schemas.
    const string stockRoundTrip = serializeObject(*stockCopy, "System", false);
    const string fastRoundTrip = serializeObject(*fastCopy, "System", false);
    ASSERT_EQUAL(stockRoundTrip, fastRoundTrip);
    ASSERT_EQUAL(xml, fastRoundTrip);
}

static void testUnsupportedForceFallback() {
    System system;
    system.addParticle(12.0);
    system.addParticle(13.0);
    CustomBondForce* force = new CustomBondForce("0.5*k*(r-r0)^2");
    force->addPerBondParameter("r0");
    force->addPerBondParameter("k");
    vector<double> parameters;
    parameters.push_back(0.2);
    parameters.push_back(125.0);
    force->addBond(0, 1, parameters);
    system.addForce(force);

    const string xml = serializeObject(system, "System", false);
    unique_ptr<System> copy = deserializeObject<System>(xml, true);
    ASSERT_EQUAL(1, copy->getNumForces());
    ASSERT(dynamic_cast<CustomBondForce*>(&copy->getForce(0)) != NULL);
    ASSERT_EQUAL(xml, serializeObject(*copy, "System", false));
}

static void testEscapedStringFallback() {
    System system;
    system.addParticle(1.0);
    system.addParticle(1.0);
    HarmonicBondForce* force = new HarmonicBondForce();
    force->setName("bond & fallback");
    force->addBond(0, 1, 0.1, 10.0);
    system.addForce(force);

    const string xml = serializeObject(system, "System", false);
    ASSERT(xml.find("bond &amp; fallback") != string::npos);
    unique_ptr<System> copy = deserializeObject<System>(xml, true);
    ASSERT_EQUAL(string("bond & fallback"), copy->getForce(0).getName());
}

static void testIntegratorParametersFallback() {
    State::StateBuilder builder(1.0, 9);
    builder.setPeriodicBoxVectors(Vec3(1, 0, 0), Vec3(0, 1, 0), Vec3(0, 0, 1));
    SerializationNode& parameters = builder.updateIntegratorParameters();
    parameters.setIntProperty("version", 17);
    parameters.createChildNode("Global").setStringProperty("name", "x").setDoubleProperty("value", 2.0);
    State state = builder.takeState();

    const string xml = serializeObject(state, "State", true);
    unique_ptr<State> copy = deserializeObject<State>(xml, true);
    ASSERT(copy->getDataTypes() & State::IntegratorParameters);
    // getIntegratorParameters() is intentionally private on State.  Verify that
    // the fallback preserved it by round-tripping the reconstructed State
    // through the stock serializer and comparing canonical XML.
    const string roundTrip = serializeObject(*copy, "State", false);
    ASSERT_EQUAL(xml, roundTrip);
}

static string forceOpeningTag(const string& xml) {
    const size_t start = xml.find("<Force ");
    ASSERT(start != string::npos);
    const size_t end = xml.find('>', start);
    ASSERT(end != string::npos);
    return xml.substr(start, end-start+1);
}

static string replaceForceOpeningTag(const string& xml, const string& tag) {
    const size_t start = xml.find("<Force ");
    ASSERT(start != string::npos);
    const size_t end = xml.find('>', start);
    ASSERT(end != string::npos);
    return xml.substr(0, start)+tag+xml.substr(end+1);
}

static string removeTagAttribute(string tag, const string& name) {
    const string prefix = " "+name+"=\"";
    const size_t start = tag.find(prefix);
    if (start == string::npos)
        return tag;
    const size_t valueEnd = tag.find('"', start+prefix.size());
    ASSERT(valueEnd != string::npos);
    tag.erase(start, valueEnd-start+1);
    return tag;
}

static string setTagAttribute(string tag, const string& name, const string& value) {
    const string prefix = " "+name+"=\"";
    const size_t start = tag.find(prefix);
    ASSERT(start != string::npos);
    const size_t valueStart = start+prefix.size();
    const size_t valueEnd = tag.find('"', valueStart);
    ASSERT(valueEnd != string::npos);
    tag.replace(valueStart, valueEnd-valueStart, value);
    return tag;
}

template <class ForceType>
static string createSingleForceSystemXml(ForceType* force) {
    System system;
    system.addParticle(1.0);
    system.addParticle(1.0);
    system.addParticle(1.0);
    system.addParticle(1.0);
    system.addForce(force);
    return serializeObject(system, "System", false);
}

static void testHistoricalPeriodicForceVersions() {
    vector<string> version2Xml;

    HarmonicBondForce* bond = new HarmonicBondForce();
    bond->setUsesPeriodicBoundaryConditions(true);
    bond->addBond(0, 1, 0.2, 10.0);
    version2Xml.push_back(createSingleForceSystemXml(bond));

    HarmonicAngleForce* angle = new HarmonicAngleForce();
    angle->setUsesPeriodicBoundaryConditions(true);
    angle->addAngle(0, 1, 2, 1.5, 20.0);
    version2Xml.push_back(createSingleForceSystemXml(angle));

    PeriodicTorsionForce* torsion = new PeriodicTorsionForce();
    torsion->setUsesPeriodicBoundaryConditions(true);
    torsion->addTorsion(0, 1, 2, 3, 2, 0.5, 3.0);
    version2Xml.push_back(createSingleForceSystemXml(torsion));

    for (size_t i = 0; i < version2Xml.size(); ++i) {
        // Current version 2 requires usesPeriodic, matching the stock proxies.
        unique_ptr<System> v2 = deserializeObject<System>(version2Xml[i], true);
        ASSERT_EQUAL(1, v2->getNumForces());

        string missing = replaceForceOpeningTag(version2Xml[i],
                                                removeTagAttribute(forceOpeningTag(version2Xml[i]), "usesPeriodic"));
        assertThrowsOpenMM([&]() { deserializeObject<System>(missing, true); });
        assertThrowsOpenMM([&]() { deserializeObject<System>(missing, false); });

        // Historical version 1 predates usesPeriodic and must continue to load.
        string v1tag = forceOpeningTag(version2Xml[i]);
        v1tag = setTagAttribute(v1tag, "version", "1");
        v1tag = removeTagAttribute(v1tag, "usesPeriodic");
        const string v1xml = replaceForceOpeningTag(version2Xml[i], v1tag);
        unique_ptr<System> v1 = deserializeObject<System>(v1xml, true);
        ASSERT_EQUAL(1, v1->getNumForces());
    }
}

static string removeXmlSection(string xml, const string& name) {
    const size_t forceStart = xml.find("<Force ");
    ASSERT(forceStart != string::npos);
    const string open = "<"+name;
    size_t start = xml.find(open, forceStart);
    ASSERT(start != string::npos);
    const size_t openEnd = xml.find('>', start);
    ASSERT(openEnd != string::npos);
    size_t eraseEnd;
    if (openEnd > start && xml[openEnd-1] == '/')
        eraseEnd = openEnd+1;
    else {
        const string close = "</"+name+">";
        const size_t closeStart = xml.find(close, openEnd+1);
        ASSERT(closeStart != string::npos);
        eraseEnd = closeStart+close.size();
    }
    if (eraseEnd < xml.size() && xml[eraseEnd] == '\n')
        ++eraseEnd;
    xml.erase(start, eraseEnd-start);
    return xml;
}

static void testHistoricalNonbondedVersions() {
    System system;
    system.addParticle(1.0);
    system.addParticle(1.0);
    NonbondedForce* force = new NonbondedForce();
    force->addGlobalParameter("lambda", 0.5);
    force->addParticle(0.2, 0.3, 0.4);
    force->addParticle(-0.2, 0.3, 0.4);
    const int exception = force->addException(0, 1, 0.0, 0.3, 0.0);
    force->addParticleParameterOffset("lambda", 0, 0.1, 0.0, 0.0);
    force->addExceptionParameterOffset("lambda", exception, 0.1, 0.0, 0.0);
    system.addForce(force);
    const string v4xml = serializeObject(system, "System", false);

    for (int version = 1; version <= 4; ++version) {
        string xml = v4xml;
        string tag = forceOpeningTag(xml);
        stringstream versionText;
        versionText << version;
        tag = setTagAttribute(tag, "version", versionText.str());
        xml = replaceForceOpeningTag(xml, tag);
        if (version < 3) {
            xml = removeXmlSection(xml, "GlobalParameters");
            xml = removeXmlSection(xml, "ParticleOffsets");
            xml = removeXmlSection(xml, "ExceptionOffsets");
        }
        unique_ptr<System> copy = deserializeObject<System>(xml, true);
        ASSERT_EQUAL(1, copy->getNumForces());
        NonbondedForce* decoded = dynamic_cast<NonbondedForce*>(&copy->getForce(0));
        ASSERT(decoded != NULL);
        ASSERT_EQUAL(2, decoded->getNumParticles());
        ASSERT_EQUAL(1, decoded->getNumExceptions());
    }
}

static void testCanonicalBooleanAndCheckedRewind() {
    // Canonical OpenMM boolean properties are emitted as numeric 0/1.  Make a
    // supported v2 force non-canonical, then arrange for the fallback rewind to
    // fail.  If the fast parser incorrectly accepts textual "true", it returns
    // successfully and this test fails; rejecting it reaches the checked rewind
    // and produces the expected OpenMMException.
    System system;
    system.addParticle(1.0);
    system.addParticle(1.0);
    HarmonicBondForce* force = new HarmonicBondForce();
    force->setUsesPeriodicBoundaryConditions(true);
    force->addBond(0, 1, 0.2, 10.0);
    system.addForce(force);
    string xml = serializeObject(system, "System", false);
    const string canonical = "usesPeriodic=\"1\"";
    const size_t boolPos = xml.find(canonical);
    ASSERT(boolPos != string::npos);
    xml.replace(boolPos, canonical.size(), "usesPeriodic=\"true\"");

    setFastXml(true);
    FailOnSeekStringBuf boolBuffer(xml, 4); // probe x2, pre-fast rewind, fallback rewind
    istream boolStream(&boolBuffer);
    assertThrowsOpenMM([&]() {
        unique_ptr<System> copy(XmlSerializer::deserialize<System>(boolStream));
    });

    // Independently verify that a failed rewind is never ignored.  The third
    // positional seek is the checked rewind immediately before fast parsing.
    const string stateXml =
        "<?xml version=\"1.0\" ?>\n"
        "<State type=\"State\" version=\"1\" time=\"0\" stepCount=\"0\">\n"
        "\t<PeriodicBoxVectors>\n"
        "\t\t<A x=\"1\" y=\"0\" z=\"0\"/>\n"
        "\t\t<B x=\"0\" y=\"1\" z=\"0\"/>\n"
        "\t\t<C x=\"0\" y=\"0\" z=\"1\"/>\n"
        "\t</PeriodicBoxVectors>\n"
        "</State>";
    FailOnSeekStringBuf rewindBuffer(stateXml, 3);
    istream rewindStream(&rewindBuffer);
    assertThrowsOpenMM([&]() {
        unique_ptr<State> copy(XmlSerializer::deserialize<State>(rewindStream));
    });
}

static void testMalformedInputAndStreamPosition() {
    const string stateXml =
        "<?xml version=\"1.0\" ?>\n"
        "<State type=\"State\" version=\"1\" time=\"2.5\" stepCount=\"4\">\n"
        "\t<PeriodicBoxVectors>\n"
        "\t\t<A x=\"1\" y=\"0\" z=\"0\"/>\n"
        "\t\t<B x=\"0\" y=\"1\" z=\"0\"/>\n"
        "\t\t<C x=\"0\" y=\"0\" z=\"1\"/>\n"
        "\t</PeriodicBoxVectors>\n"
        "</State>";

    setFastXml(true);
    stringstream stream(stateXml+"TRAIL");
    stream.exceptions(ios_base::badbit | ios_base::failbit);
    unique_ptr<State> state(XmlSerializer::deserialize<State>(stream));
    ASSERT_EQUAL(2.5, state->getTime());
    char trailing[6] = {0, 0, 0, 0, 0, 0};
    stream.read(trailing, 5);
    ASSERT_EQUAL(string("TRAIL"), string(trailing, 5));

    string malformed = stateXml;
    const string needle = "<A x=\"1\" y=\"0\" z=\"0\"/>";
    const size_t pos = malformed.find(needle);
    ASSERT(pos != string::npos);
    malformed.replace(pos, needle.size(), "<A x=\"1\" y=\"0\"/>");
    assertThrowsOpenMM([&]() { deserializeObject<State>(malformed, true); });
    assertThrowsOpenMM([&]() { deserializeObject<State>(malformed, false); });
}

} // anonymous namespace

int main() {
    EnvironmentGuard fastXmlGuard("OPENMM_FAST_XML");
    EnvironmentGuard fastWriteGuard("OPENMM_FAST_STATE_WRITE");
    try {
        testLargeStateFastWriterAndReaderParity();
        testSupportedSystemFastReaderParity();
        testUnsupportedForceFallback();
        testEscapedStringFallback();
        testIntegratorParametersFallback();
        testHistoricalPeriodicForceVersions();
        testHistoricalNonbondedVersions();
        testCanonicalBooleanAndCheckedRewind();
        testMalformedInputAndStreamPosition();
    }
    catch (const exception& e) {
        cout << "exception: " << e.what() << endl;
        return 1;
    }
    cout << "Done" << endl;
    return 0;
}
