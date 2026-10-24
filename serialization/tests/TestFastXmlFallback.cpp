/* -------------------------------------------------------------------------- *
 * Force-scoped XML fallback regression tests.                                 *
 * -------------------------------------------------------------------------- */
#include "openmm/internal/AssertionUtilities.h"
#include "openmm/CMMotionRemover.h"
#include "openmm/CustomBondForce.h"
#include "openmm/HarmonicBondForce.h"
#include "openmm/NonbondedForce.h"
#include "openmm/System.h"
#include "openmm/serialization/NonbondedForceProxy.h"
#include "openmm/serialization/XmlSerializer.h"
#include <cstdlib>
#include <iostream>
#include <memory>
#include <sstream>
#include <string>
#include <vector>
using namespace OpenMM;
using namespace std;

namespace {
void fastXml(bool enabled) {
#ifdef _WIN32
    _putenv_s("OPENMM_FAST_XML", enabled ? "1" : "0");
#else
    setenv("OPENMM_FAST_XML", enabled ? "1" : "0", 1);
#endif
}
string serialize(const System& system) {
    fastXml(false);
    stringstream out;
    XmlSerializer::serialize(&system, "System", out);
    return out.str();
}
unique_ptr<System> read(const string& xml, bool fast, bool suffix = true) {
    fastXml(fast);
    stringstream in(xml+(suffix ? "TAIL" : ""));
    unique_ptr<System> result(XmlSerializer::deserialize<System>(in));
    if (fast && suffix) {
        string tail;
        in >> tail;
        ASSERT_EQUAL(string("TAIL"), tail);
    }
    return result;
}
void equivalentWithLocalFallback(const string& xml) {
    unique_ptr<System> stock = read(xml, false, false);
    unique_ptr<System> fast = read(xml, true);
    ASSERT_EQUAL(serialize(*stock), serialize(*fast));
}
void replaceOnce(string& text, const string& oldValue, const string& newValue) {
    size_t at = text.find(oldValue);
    ASSERT(at != string::npos);
    text.replace(at, oldValue.size(), newValue);
}
NonbondedForce* nonbonded() {
    NonbondedForce* force = new NonbondedForce();
    force->addParticle(0.1,0.2,0.3);
    force->addParticle(-0.1,0.2,0.3);
    force->addException(0,1,0,1,0);
    return force;
}
CustomBondForce* custom(int count = 1) {
    CustomBondForce* force = new CustomBondForce("k*(r-0.1)^2");
    force->setName("custom & <escaped> \"force\"");
    force->addPerBondParameter("k");
    for (int i=0; i<count; ++i)
        force->addBond(0,1,vector<double>(1,0.125*(i+1)));
    return force;
}
void testMixedForces() {
    for (int position=0; position<3; ++position) {
        System system;
        system.addParticle(1); system.addParticle(2);
        if (position == 0) system.addForce(custom());
        system.addForce(nonbonded());
        if (position == 1) system.addForce(custom());
        system.addForce(new CMMotionRemover(3));
        if (position == 2) system.addForce(custom());
        equivalentWithLocalFallback(serialize(system));
    }
}
void testKnownForceFallback() {
    System system;
    system.addParticle(1); system.addParticle(2);
    system.addForce(nonbonded());
    system.addForce(new CMMotionRemover(4));
    string original = serialize(system);
    // Direct parser rejects this child only after reading the force arrays.
    string late = original;
    replaceOnce(late,"</Force>","<Unrecognized/>\n</Force>");
    equivalentWithLocalFallback(late);
    // Escaped type names must still go through irrXML's attribute decoding.
    string escaped = original;
    replaceOnce(escaped,"type=\"NonbondedForce\"","type=\"NonbondedForc&#101;\"");
    equivalentWithLocalFallback(escaped);
    string escapedName = original;
    replaceOnce(escapedName,"name=\"NonbondedForce\"","name=\"A&amp;B&#x20;force\"");
    equivalentWithLocalFallback(escapedName);

    System bonds;
    bonds.addParticle(1); bonds.addParticle(2);
    HarmonicBondForce* force = new HarmonicBondForce();
    force->addBond(0,1,0.1,2);
    force->setUsesPeriodicBoundaryConditions(true);
    bonds.addForce(force);
    string noncanonical = serialize(bonds);
    replaceOnce(noncanonical,"usesPeriodic=\"1\"","usesPeriodic=\"2\"");
    equivalentWithLocalFallback(noncanonical);
}
void testLateFallbackAcrossInputBuffers() {
    System system;
    NonbondedForce* force = new NonbondedForce();
    for (int i=0; i<25000; ++i) {
        system.addParticle(1);
        force->addParticle(0.1,0.2,0.3);
    }
    system.addForce(force);
    system.addForce(new CMMotionRemover(4));
    string xml = serialize(system);
    ASSERT(xml.size() > 1024*1024);
    replaceOnce(xml,"</Force>","<Unrecognized/>\n</Force>");
    equivalentWithLocalFallback(xml);
}
void testNonzeroDocumentOffset() {
    System system;
    system.addParticle(1); system.addParticle(2);
    system.addForce(nonbonded()); system.addForce(custom());
    const string xml = serialize(system);
    const string prefix = "unrelated prefix bytes\n";
    stringstream input(prefix+xml+"TAIL");
    input.seekg(prefix.size());
    fastXml(true);
    unique_ptr<System> result(XmlSerializer::deserialize<System>(input));
    string tail;
    input >> tail;
    ASSERT_EQUAL(string("TAIL"),tail);
    ASSERT_EQUAL(xml,serialize(*result));
}
void testLargeSubtreeAndComments() {
    System system;
    system.addParticle(1); system.addParticle(2);
    system.addForce(nonbonded());
    system.addForce(custom(25000));
    system.addForce(new CMMotionRemover(4));
    string xml = serialize(system);
    ASSERT(xml.size() > 1024*1024);
    replaceOnce(xml,"<Bonds>","<Bonds><!-- comment with < and > --><?ignored value?>");
    equivalentWithLocalFallback(xml);
}
class ReplacementNonbondedProxy : public NonbondedForceProxy {
public:
    mutable int calls;
    ReplacementNonbondedProxy() : calls(0) {}
    void* deserialize(const SerializationNode& node) const override {
        ++calls;
        NonbondedForce* force = reinterpret_cast<NonbondedForce*>(NonbondedForceProxy::deserialize(node));
        force->setName("replacement proxy used");
        return force;
    }
};
class ProxyGuard {
public:
    ProxyGuard(const SerializationProxy* replacement) : original(&SerializationProxy::getProxy(typeid(NonbondedForce))) {
        SerializationProxy::registerProxy(typeid(NonbondedForce), replacement);
    }
    ~ProxyGuard() { SerializationProxy::registerProxy(typeid(NonbondedForce), original); }
private:
    const SerializationProxy* original;
};
void testReplacedProxy() {
    System system;
    system.addParticle(1); system.addParticle(2);
    system.addForce(nonbonded());
    system.addForce(new CMMotionRemover(2));
    string xml = serialize(system);
    ReplacementNonbondedProxy replacement;
    ProxyGuard guard(&replacement);
    unique_ptr<System> result = read(xml,true);
    ASSERT_EQUAL(1,replacement.calls);
    ASSERT_EQUAL(string("replacement proxy used"),result->getForce(0).getName());
}
class ThrowingNonbondedProxy : public NonbondedForceProxy {
public:
    mutable int calls;
    ThrowingNonbondedProxy() : calls(0) {}
    void* deserialize(const SerializationNode&) const override {
        ++calls;
        throw OpenMMException("intentional proxy failure");
    }
};
void testFailingProxyIsNotRetried() {
    System system;
    system.addParticle(1); system.addParticle(2);
    system.addForce(nonbonded());
    string xml = serialize(system);
    ThrowingNonbondedProxy replacement;
    ProxyGuard guard(&replacement);
    bool threw = false;
    try {
        read(xml,true);
    }
    catch (const OpenMMException& e) {
        threw = true;
        ASSERT_EQUAL(string("intentional proxy failure"),string(e.what()));
    }
    ASSERT(threw);
    ASSERT_EQUAL(1,replacement.calls);
}
class FailingSubtreeSeekBuf : public stringbuf {
public:
    explicit FailingSubtreeSeekBuf(const string& value) : stringbuf(value,ios::in), failed(false) {}
    bool failed;
protected:
    pos_type seekpos(pos_type position, ios_base::openmode mode) override {
        // Root rewinds remain available, but a force-local rewind is denied.
        if (position > pos_type(0)) {
            failed = true;
            return pos_type(off_type(-1));
        }
        return stringbuf::seekpos(position,mode);
    }
};
void testFailedLocalRewind() {
    System system;
    system.addParticle(1); system.addParticle(2);
    system.addForce(nonbonded()); system.addForce(custom());
    string xml = serialize(system);
    fastXml(true);
    FailingSubtreeSeekBuf buffer(xml);
    istream input(&buffer);
    unique_ptr<System> result(XmlSerializer::deserialize<System>(input));
    ASSERT(buffer.failed);
    ASSERT_EQUAL(xml,serialize(*result));
}
}
int main() {
    try {
        testMixedForces();
        testKnownForceFallback();
        testLargeSubtreeAndComments();
        testLateFallbackAcrossInputBuffers();
        testNonzeroDocumentOffset();
        testReplacedProxy();
        testFailingProxyIsNotRetried();
        testFailedLocalRewind();
    }
    catch (const exception& e) {
        cerr << "exception: " << e.what() << endl;
        return 1;
    }
    cout << "Done" << endl;
    return 0;
}
