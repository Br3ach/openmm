/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * Regression tests for the streaming State/System XML fast paths.            *
 * -------------------------------------------------------------------------- */

#include "openmm/internal/AssertionUtilities.h"
#include "openmm/CMMotionRemover.h"
#include "openmm/Context.h"
#include "openmm/CustomBondForce.h"
#include "openmm/HarmonicAngleForce.h"
#include "openmm/HarmonicBondForce.h"
#include "openmm/MonteCarloBarostat.h"
#include "openmm/NonbondedForce.h"
#include "openmm/OpenMMException.h"
#include "openmm/PeriodicTorsionForce.h"
#include "openmm/Platform.h"
#include "openmm/State.h"
#include "openmm/System.h"
#include "openmm/VerletIntegrator.h"
#include "openmm/VirtualSite.h"
#include "openmm/serialization/SerializationNode.h"
#include "openmm/serialization/XmlSerializer.h"
#include <algorithm>
#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <fstream>
#include <iterator>
#include <map>
#include <memory>
#include <sstream>
#include <string>
#include <streambuf>
#include <vector>
#include <utility>

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
            stringbuf(value, ios_base::in), failOnSeek(failOnSeek), seekCount(0), legacyReader(false) {
    }
    bool usedLegacyReader() const { return legacyReader; }

protected:
    pos_type seekoff(off_type offset, ios_base::seekdir direction, ios_base::openmode mode) override {
        if (direction == ios_base::end)
            legacyReader = true;
        return stringbuf::seekoff(offset, direction, mode);
    }
    pos_type seekpos(pos_type position, ios_base::openmode mode = ios_base::in | ios_base::out) override {
        ++seekCount;
        if (seekCount == failOnSeek)
            return pos_type(off_type(-1));
        return stringbuf::seekpos(position, mode);
    }

private:
    int failOnSeek;
    int seekCount;
    bool legacyReader;
};

static void setFastXml(bool enabled) {
    setEnvironment("OPENMM_FAST_XML", enabled ? "1" : "0");
}

class ShortRefillStringBuf : public FailOnSeekStringBuf {
public:
    ShortRefillStringBuf(const string& xml, streamsize limit) :
            FailOnSeekStringBuf(xml, 0), limit(limit) {
    }
protected:
    streamsize xsgetn(char* target, streamsize count) override {
        // Limit only the fast reader's refills, leaving its root probe intact.
        if (limit != 0 && count == 1024*1024)
            count = min(count, limit);
        return stringbuf::xsgetn(target, count);
    }
private:
    streamsize limit;
};


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

static State createLargeState(bool integratorParameters=false) {
    const int count = 100000; // Exercise large arrays without a full reader tree.
    vector<Vec3> positions;
    vector<Vec3> velocities, forces;
    positions.reserve(count);
    velocities.reserve(count);
    forces.reserve(count);
    for (int i = 0; i < count; ++i) {
        positions.push_back(Vec3(0.001*i, -0.002*i, 0.003*i));
        velocities.push_back(Vec3(-0.004*i, 0.005*i, -0.006*i));
        forces.push_back(Vec3(0.007*i, -0.008*i, 0.009*i));
    }

    map<string, double> parameters;
    parameters["alpha"] = 1.25;
    parameters["beta"] = -3.5;

    State::StateBuilder builder(2.5, 1234567890123LL);
    builder.setPeriodicBoxVectors(Vec3(4.0, 0.1, 0.2), Vec3(0.0, 5.0, 0.3), Vec3(0.0, 0.0, 6.0));
    builder.setPositions(std::move(positions));
    builder.setVelocities(std::move(velocities));
    builder.setForces(std::move(forces));
    builder.setParameters(std::move(parameters));
    builder.setEnergy(7.25, -19.5);
    if (integratorParameters)
        builder.updateIntegratorParameters().setIntProperty("version", 17);
    return builder.takeState();
}

static void testLargeStateReaderParity() {
    State state = createLargeState();
    ASSERT_EQUAL(31, state.getDataTypes()); // Core28's full-State request: positions, velocities, forces, energy, parameters.

    // The Core28 build always retains the complete-node stock writer.
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

static void testExceptionIndexReconstruction() {
    unique_ptr<System> original = createSupportedSystem();
    const string xml = serializeObject(*original, "System", false);
    auto exercise = [&](NonbondedForce& force) {
        assertThrowsOpenMM([&]() { force.setExceptionParameters(-1, 0, 1, 0, 1, 0); });
        assertThrowsOpenMM([&]() { force.setExceptionParameters(force.getNumExceptions(), 0, 1, 0, 1, 0); });
        force.setExceptionParameters(0, 0, 1, .5, .25, .125); // Numerical edit before any rebuild.
        force.setExceptionParameters(0, 2, 3, .5, .25, .125);
        // Preserve the stock lookup's historical keys after a parent edit.
        ASSERT_EQUAL(1, force.addException(2, 3, .25, .5, .125));
        assertThrowsOpenMM([&]() { force.addException(1, 0, 0, 1, 0); });
        ASSERT_EQUAL(0, force.addException(1, 0, .25, .5, .125, true));
        ASSERT_EQUAL(2, force.addException(2, 2, 0, 1, 0));
        assertThrowsOpenMM([&]() { force.addException(2, 2, 0, 1, 0); });
        ASSERT_EQUAL(2, force.addException(2, 2, .5, 1, 0, true));
        return serializeObject(force, "Force", false);
    };
    vector<string> expected;
    for (bool fast : {false, true}) {
        unique_ptr<System> system = deserializeObject<System>(xml, fast);
        const NonbondedForce& source = dynamic_cast<const NonbondedForce&>(system->getForce(0));
        ASSERT_EQUAL(xml, serializeObject(*system, "System", false)); // Reads/serialization before mutation.
        NonbondedForce copied(source);
        NonbondedForce assigned(source);
        assigned.addException(1, 2, 0, 1, 0); // Assign over an indexed target.
        assigned = source;
        NonbondedForce toMove(source);
        NonbondedForce moved(std::move(toMove));
        NonbondedForce toMoveAssign(source), moveAssigned(source);
        moveAssigned.addException(1, 2, 0, 1, 0);
        moveAssigned = std::move(toMoveAssign);
        vector<string> found;
        for (NonbondedForce* force : {&copied, &assigned, &moved, &moveAssigned})
            found.push_back(exercise(*force));
        if (!fast)
            expected = found;
        else
            ASSERT_EQUAL_CONTAINERS(expected, found);
        NonbondedForce added(source);
        assertThrowsOpenMM([&]() { added.addException(1, 0, 0, 1, 0); });
        ASSERT_EQUAL(1, added.addException(1, 2, 0, 1, 0));
    }
}

static void testExceptionContextUpdate() {
    System original;
    original.addParticle(12);
    original.addParticle(12);
    NonbondedForce* force = new NonbondedForce();
    force->addParticle(0, .1, 0);
    force->addParticle(0, .1, 0);
    force->addException(0, 1, .5, .1, 0);
    original.addForce(force);
    const string xml = serializeObject(original, "System", false);
    unique_ptr<System> stock = deserializeObject<System>(xml, false);
    unique_ptr<System> fast = deserializeObject<System>(xml, true);
    VerletIntegrator stockIntegrator(.001), fastIntegrator(.001);
    Platform& platform = Platform::getPlatformByName("Reference");
    Context stockContext(*stock, stockIntegrator, platform);
    Context fastContext(*fast, fastIntegrator, platform);
    const vector<Vec3> positions = {Vec3(0, 0, 0), Vec3(1, 0, 0)};
    stockContext.setPositions(positions);
    fastContext.setPositions(positions);
    const int types = State::Energy | State::Forces;
    compareStates(stockContext.getState(types), fastContext.getState(types));
    const double before = fastContext.getState(types).getPotentialEnergy();
    for (System* system : {stock.get(), fast.get()})
        dynamic_cast<NonbondedForce&>(system->getForce(0)).setExceptionParameters(0, 0, 1, .75, .1, 0);
    dynamic_cast<NonbondedForce&>(stock->getForce(0)).updateParametersInContext(stockContext);
    dynamic_cast<NonbondedForce&>(fast->getForce(0)).updateParametersInContext(fastContext);
    compareStates(stockContext.getState(types), fastContext.getState(types));
    ASSERT(fastContext.getState(types).getPotentialEnergy() != before);
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

class CountingProxy : public SerializationProxy {
public:
    CountingProxy(const type_info& type) : SerializationProxy(getProxy(type).getTypeName()),
            calls(0), type(type), original(getProxy(type)) {
        registerProxy(type, this);
    }
    ~CountingProxy() { registerProxy(type, &original); }
    void serialize(const void* object, SerializationNode& node) const override {
        original.serialize(object, node);
    }
    void* deserialize(const SerializationNode& node) const override {
        ++calls;
        return original.deserialize(node);
    }
    mutable int calls;
private:
    const type_info& type;
    const SerializationProxy& original;
};

static void testRegisteredProxyOverrides() {
    unique_ptr<System> system = createSupportedSystem();
    const string xml = serializeObject(*system, "System", false);
    {
        CountingProxy proxy(typeid(System));
        for (bool fast : {false, true}) {
            unique_ptr<System> copy = deserializeObject<System>(xml, fast);
            ASSERT_EQUAL(xml, serializeObject(*copy, "System", false));
        }
        ASSERT_EQUAL(2, proxy.calls);
    }
    State state = createLargeState();
    const string stateXml = serializeObject(state, "State", false);
    {
        CountingProxy proxy(typeid(State));
        for (bool fast : {false, true}) {
            unique_ptr<State> copy = deserializeObject<State>(stateXml, fast);
            compareStates(state, *copy);
        }
        ASSERT_EQUAL(2, proxy.calls);
    }
    for (int i = 0; i < system->getNumForces(); ++i) {
        CountingProxy proxy(typeid(system->getForce(i)));
        for (bool fast : {false, true}) {
            unique_ptr<System> copy = deserializeObject<System>(xml, fast);
            ASSERT_EQUAL(xml, serializeObject(*copy, "System", false));
        }
        ASSERT_EQUAL(2, proxy.calls);
    }
    System empty;
    const string emptyXml = serializeObject(empty, "System", false);
    CountingProxy unused(typeid(HarmonicBondForce));
    setFastXml(true);
    FailOnSeekStringBuf buffer(emptyXml, 0);
    istream stream(&buffer);
    unique_ptr<System> copy(XmlSerializer::deserialize<System>(stream));
    ASSERT(!buffer.usedLegacyReader());
    ASSERT_EQUAL(0, unused.calls);
}

static void testFileStreamPosition() {
    State::StateBuilder builder(2.5, 7);
    builder.setPeriodicBoxVectors(Vec3(1, 0, 0), Vec3(0, 1, 0), Vec3(0, 0, 1));
    State state = builder.takeState();
    System system;
    system.addParticle(12);
    const string tail = "\r\nTAIL\r\nSECOND\r\n";
    const char* path = "TestFastXmlTextStream.xml";
    try {
        for (string xml : {serializeObject(state, "State", false), serializeObject(system, "System", false)}) {
            xml.erase(xml.size()-1); // Start the known tail immediately after the root.
            for (size_t padding : {size_t(0), size_t(1024*1024)}) {
                string padded = xml;
                padded.insert(padded.find("<PeriodicBoxVectors"), padding, ' ');
                for (bool crlf : {false, true}) {
                    string text;
                    for (char c : padded) {
                        if (crlf && c == '\n')
                            text += '\r';
                        text += c;
                    }
                    for (bool binary : {false, true}) {
                        { ofstream output(path, ios_base::binary); output << text << tail; ASSERT(output.good()); }
                        setFastXml(true);
                        ifstream input(path, binary ? ios_base::binary : ios_base::in);
                        if (xml.find("<State ") != string::npos) {
                            unique_ptr<State> copy(XmlSerializer::deserialize<State>(input));
                            compareStates(state, *copy);
                        }
                        else {
                            unique_ptr<System> copy(XmlSerializer::deserialize<System>(input));
                            ASSERT_EQUAL(12, copy->getParticleMass(0));
                        }
                        const string found((istreambuf_iterator<char>(input)), istreambuf_iterator<char>());
#ifdef _WIN32
                        ASSERT_EQUAL(binary ? tail : string("\nTAIL\nSECOND\n"), found);
#else
                        ASSERT_EQUAL(tail, found);
#endif
                    }
                }
            }
        }
    }
    catch (...) {
        remove(path);
        throw;
    }
    ASSERT_EQUAL(0, remove(path));
}

static void testReadAheadRecovery() {
    class ReplayStringBuf : public FailOnSeekStringBuf {
    public:
        ReplayStringBuf(const string& xml, int fault) : FailOnSeekStringBuf(xml, fault == 3 ? 4 : 0),
                fault(fault), bodyRead(false), replay(false) {}
    protected:
        pos_type seekpos(pos_type position, ios_base::openmode mode) override {
            if (bodyRead)
                replay = true;
            return FailOnSeekStringBuf::seekpos(position, mode);
        }
        streamsize xsgetn(char* target, streamsize count) override {
            if (replay && !usedLegacyReader()) {
                const int failure = fault;
                fault = 0;
                if (failure == 1)
                    return 0;
                if (failure == 2)
                    throw ios_base::failure("injected replay failure");
                count = min(count, streamsize(2));
            }
            if (count == 1024*1024)
                bodyRead = true;
            return stringbuf::xsgetn(target, count);
        }
    private:
        int fault;
        bool bodyRead, replay;
    };
    State::StateBuilder builder(2.5, 7);
    builder.setPeriodicBoxVectors(Vec3(1, 0, 0), Vec3(0, 1, 0), Vec3(0, 0, 1));
    State state = builder.takeState();
    const string xml = serializeObject(state, "State", false);
    // Short positive reads, zero read, read exception, and failed finish seek.
    for (int fault : {0, 1, 2, 3}) {
        setFastXml(true);
        ReplayStringBuf buffer(xml+"TRAIL", fault);
        istream stream(&buffer);
        unique_ptr<State> copy(XmlSerializer::deserialize<State>(stream));
        compareStates(state, *copy);
        ASSERT_EQUAL(fault != 0, buffer.usedLegacyReader());
        if (fault == 0) {
            const string tail((istreambuf_iterator<char>(stream)), istreambuf_iterator<char>());
            ASSERT_EQUAL(string("\nTRAIL"), tail);
        }
    }
}

static void testThreeParticleVirtualSiteAttributes() {
    System system;
    for (int i = 0; i < 4; ++i)
        system.addParticle(i == 3 ? 0 : 12);
    system.setVirtualSite(3, new ThreeParticleAverageSite(0, 1, 2, 0.125, 0.25, 0.625));
    const string xml = serializeObject(system, "System", false);
    const size_t start = xml.find("<ThreeParticleAverageSite");
    ASSERT(start != string::npos);
    const size_t length = xml.find("/>", start)+2-start;
    const string attributes[] = {"p1=\"0\"", "p2=\"1\"", "p3=\"2\"", "w1=\".125\"", "w2=\".25\"", "w3=\".625\""};
    const string canonical = "<ThreeParticleAverageSite p1=\"0\" p2=\"1\" p3=\"2\" w1=\".125\" w2=\".25\" w3=\".625\"";
    const vector<string> valid = {
        canonical+"/>",
        "<ThreeParticleAverageSite w3='.625' p3='2' w2='2.5e-1' p1='+0' w1='0.125' p2='1'/>",
        canonical+" ignored=\"not numeric\" ignored=\"still ignored\"/>",
        canonical+"></ThreeParticleAverageSite>",
        "<ThreeParticleAverageSite p1='-2147483648' p2='2147483647' p3='-1' w1='.125' w2='.25' w3='.625'/>"
    };
    for (const string& tag : valid) {
        string input = xml;
        input.replace(start, length, tag);
        unique_ptr<System> stock = deserializeObject<System>(input, false);
        const string expected = serializeObject(*stock, "System", false);
        setFastXml(true);
        FailOnSeekStringBuf buffer(input, 0);
        istream stream(&buffer);
        unique_ptr<System> fast(XmlSerializer::deserialize<System>(stream));
        ASSERT(!buffer.usedLegacyReader());
        ASSERT_EQUAL(expected, serializeObject(*fast, "System", false));
        setFastXml(true);
        ShortRefillStringBuf shortBuffer(input, 7);
        istream shortStream(&shortBuffer);
        unique_ptr<System> split(XmlSerializer::deserialize<System>(shortStream));
        ASSERT(!shortBuffer.usedLegacyReader());
        ASSERT_EQUAL(expected, serializeObject(*split, "System", false));
    }
    vector<string> invalid;
    for (const string& attribute : attributes) {
        string missing = canonical;
        missing.erase(missing.find(attribute), attribute.size());
        invalid.push_back(missing+"/>");
        invalid.push_back(canonical+" "+attribute+"/>");
    }
    invalid.push_back(canonical+" bad=unquoted/>");
    invalid.push_back(canonical+"/garbage>");
    for (const string& replacement : {"p1=\"2147483648\"", "p1=\"-2147483649\"", "p1=\"abc\"", "w1=\".125x\""}) {
        string bad = canonical;
        const string original = replacement[0] == 'p' ? attributes[0] : attributes[3];
        bad.replace(bad.find(original), original.size(), replacement);
        invalid.push_back(bad+"/>");
    }
    for (const string& tag : invalid) {
        string input = xml;
        input.replace(start, length, tag);
        setFastXml(true);
        FailOnSeekStringBuf buffer(input, 4);
        istream stream(&buffer);
        assertThrowsOpenMM([&]() { unique_ptr<System> copy(XmlSerializer::deserialize<System>(stream)); });
    }
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

static void testTagBufferBoundaries() {
    System system;
    system.addParticle(1.0);
    system.addParticle(1.0);
    HarmonicBondForce* force = new HarmonicBondForce();
    force->setName("quoted>name");
    force->addBond(0, 1, 0.2, 10.0);
    system.addForce(force);
    const string canonical = serializeObject(system, "System", false);
    string xml = canonical;
    const size_t escaped = xml.find("quoted&gt;name");
    ASSERT(escaped != string::npos);
    xml.replace(escaped, string("quoted&gt;name").size(), "quoted>name");
    const size_t start = xml.find("<Force ");
    const size_t end = xml.find(">\n", start);
    const size_t openQuote = xml.find("name=\"", start)+5;
    const size_t closeQuote = xml.find('"', openQuote+1);
    const size_t quotedEnd = xml.find('>', openQuote);
    const size_t block = 1024*1024;
    auto check = [&](const string& input, streamsize limit) {
        // The legacy reader seeks to the end to determine the document size.
        setFastXml(true);
        ShortRefillStringBuf buffer(input+"TRAIL", limit);
        istream stream(&buffer);
        unique_ptr<System> copy(XmlSerializer::deserialize<System>(stream));
        ASSERT(!buffer.usedLegacyReader());
        ASSERT_EQUAL(canonical, serializeObject(*copy, "System", false));
        if (input.back() == '\n')
            ASSERT_EQUAL('\n', stream.get());
        char trailing[5];
        stream.read(trailing, 5);
        ASSERT_EQUAL(string("TRAIL"), string(trailing, 5));
    };
    for (char quote : {'"', '\''}) {
        string quoted = xml;
        quoted[openQuote] = quoted[closeQuote] = quote;
        for (size_t marker : {openQuote, closeQuote, quotedEnd, end}) {
            for (int delta = -1; delta <= 1; ++delta) {
                string aligned = quoted;
                aligned.insert(start, block+delta-marker, ' ');
                check(aligned, 0);
            }
        }
    }
    string prefixed = xml;
    prefixed.insert(start, "<!-- comment --><?probe instruction?>");
    prefixed.insert(0, "\xEF\xBB\xBF");
    for (streamsize limit : {1, 2, 3, 7, 31})
        check(prefixed, limit);
    for (int delta = -1; delta <= 1; ++delta) {
        string sized = xml;
        sized.insert(start+6, block+delta-(end-start+1), ' ');
        if (delta <= 0)
            check(sized, 0);
        else
            assertThrowsOpenMM([&]() { check(sized, 0); });
    }
    for (size_t cut : {closeQuote, end})
        assertThrowsOpenMM([&]() { check(xml.substr(0, cut), 0); });
}

static void testAttributeEscaping() {
    System system;
    CustomBondForce* force = new CustomBondForce("ESCAPE_SENTINEL");
    system.addForce(force);
    const string pattern = serializeObject(system, "System", false);
    const size_t offset = pattern.find("ESCAPE_SENTINEL");
    ASSERT(offset != string::npos);
    auto check = [&](const string& value, const string& escaped) {
        force->setEnergyFunction(value);
        string expected = pattern;
        expected.replace(offset, string("ESCAPE_SENTINEL").size(), escaped);
        ASSERT_EQUAL(expected, serializeObject(system, "System", false));
    };
    const pair<string, string> cases[] = {
        {"", ""}, {"-1.234e-12", "-1.234e-12"},
        {"&<>\"'", "&amp;&lt;&gt;&quot;&apos;"},
        {"&#xA9;", "&#xA9;"}, {"&#x;", "&#x;"},
        {"&#xZZ;", "&#xZZ;"}, {"&#x", "&#x"},
        {"&#xA9", "&#xA9"}, {"&#xA9&", "&#xA9&amp;"},
        {"&#xA9<", "&#xA9&lt;"}, {"&#xA9\"", "&#xA9&quot;"},
        {"&#169;", "&amp;#169;"}, {"&#XA9;", "&amp;#XA9;"},
        {"&amp;", "&amp;amp;"}, {"&&#x41;", "&amp;&#x41;"}
    };
    for (const auto& value : cases)
        check(value.first, value.second);
    for (unsigned int c = 0; c < 32; ++c) {
        char reference[8];
        snprintf(reference, sizeof(reference), "&#x%02X;", c);
        check(string(1, static_cast<char>(c)), reference);
    }
    string highBytes;
    for (unsigned int c = 128; c < 256; ++c)
        highBytes += static_cast<char>(c);
    check(highBytes, highBytes);
    check("caf\xC3\xA9", "caf\xC3\xA9");
}

} // anonymous namespace

static void testDirectProxyAndMovedArrays() {
    State state = createLargeState();
    SerializationNode node;
    node.setName("State");
    const SerializationProxy& proxy = SerializationProxy::getProxy(typeid(state));
    proxy.serialize(&state, node);
    // Direct callers must receive a complete node, without a deferred handoff.
    ASSERT_EQUAL(state.getPositions().size(), node.getChildNode("Positions").getChildren().size());
    unique_ptr<State> copy(reinterpret_cast<State*>(proxy.deserialize(node)));
    compareStates(state, *copy);

    vector<Vec3> positions(4, Vec3(1, 2, 3));
    const Vec3* storage = positions.data();
    State::StateBuilder builder(2.0, 7);
    builder.setPositions(std::move(positions));
    State moved = builder.takeState();
    ASSERT(storage == moved.getPositions().data());
    ASSERT_EQUAL(7, moved.getStepCount());
}

static void testPrepopulatedStateProxy() {
    State state = createLargeState(true);
    SerializationNode node;
    node.setName("State");
    node.createChildNode("Existing").setStringProperty("value", "preserved");
    const SerializationProxy& proxy = SerializationProxy::getProxy(typeid(state));
    proxy.serialize(&state, node);
    ASSERT_EQUAL(8, node.getChildren().size());
    ASSERT_EQUAL(string("preserved"), node.getChildren()[0].getStringProperty("value"));
    ASSERT_EQUAL(17, node.getChildNode("IntegratorParameters").getIntProperty("version"));
    unique_ptr<State> copy(reinterpret_cast<State*>(proxy.deserialize(node)));
    compareStates(state, *copy);
    ASSERT_EQUAL(serializeObject(state, "State", false), serializeObject(*copy, "State", false));
    unique_ptr<State> again(reinterpret_cast<State*>(proxy.deserialize(node)));
    compareStates(state, *again);
    unique_ptr<State> fallback = deserializeObject<State>(serializeObject(state, "State", false), true);
    compareStates(state, *fallback);
    unique_ptr<State> clone(XmlSerializer::clone(state));
    compareStates(state, *clone);
}

static void testLegacyArraySizes() {
    State::StateBuilder builder(1.0, 4);
    builder.setPeriodicBoxVectors(Vec3(1, 0, 0), Vec3(0, 1, 0), Vec3(0, 0, 1));
    builder.setPositions(vector<Vec3>(2));
    builder.setVelocities(vector<Vec3>(3));
    State mismatched = builder.takeState();
    const SerializationProxy& proxy = SerializationProxy::getProxy(typeid(mismatched));
    SerializationNode node;
    proxy.serialize(&mismatched, node);
    assertThrowsOpenMM([&]() { unique_ptr<State> copy(reinterpret_cast<State*>(proxy.deserialize(node))); });
    assertThrowsOpenMM([&]() { deserializeObject<State>(serializeObject(mismatched, "State", false), false); });

    node.getChildNode("Velocities").getChildren().pop_back();
    SerializationNode duplicate = node.getChildNode("Positions");
    duplicate.getChildren().pop_back();
    node.getChildren().push_back(duplicate);
    assertThrowsOpenMM([&]() { unique_ptr<State> copy(reinterpret_cast<State*>(proxy.deserialize(node))); });

    builder = State::StateBuilder(1.0, 4);
    builder.setPeriodicBoxVectors(Vec3(1, 0, 0), Vec3(0, 1, 0), Vec3(0, 0, 1));
    builder.setPositions(vector<Vec3>());
    builder.setVelocities(vector<Vec3>());
    State empty = builder.takeState();
    unique_ptr<State> copy = deserializeObject<State>(serializeObject(empty, "State", false), false);
    compareStates(empty, *copy);
    ASSERT_EQUAL(State::Positions | State::Velocities, copy->getDataTypes());
}

int main() {
    EnvironmentGuard fastXmlGuard("OPENMM_FAST_XML");
    try {
        testDirectProxyAndMovedArrays();
        testPrepopulatedStateProxy();
        testLegacyArraySizes();
        testLargeStateReaderParity();
        testSupportedSystemFastReaderParity();
        testExceptionIndexReconstruction();
        testExceptionContextUpdate();
        testRegisteredProxyOverrides();
        testFileStreamPosition();
        testReadAheadRecovery();
        testThreeParticleVirtualSiteAttributes();
        testUnsupportedForceFallback();
        testEscapedStringFallback();
        testIntegratorParametersFallback();
        testHistoricalPeriodicForceVersions();
        testHistoricalNonbondedVersions();
        testCanonicalBooleanAndCheckedRewind();
        testMalformedInputAndStreamPosition();
        testTagBufferBoundaries();
        testAttributeEscaping();
    }
    catch (const exception& e) {
        cout << "exception: " << e.what() << endl;
        return 1;
    }
    cout << "Done" << endl;
    return 0;
}
