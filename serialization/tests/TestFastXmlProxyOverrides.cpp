/* Regression coverage for user-replaced State and System serialization proxies. */
#include "openmm/State.h"
#include "openmm/System.h"
#include "openmm/internal/AssertionUtilities.h"
#include "openmm/serialization/StateProxy.h"
#include "openmm/serialization/SystemProxy.h"
#include "openmm/serialization/XmlSerializer.h"
#include <cstdlib>
#include <iostream>
#include <memory>
#include <sstream>
#include <string>
#include <typeinfo>

using namespace OpenMM;
using namespace std;

namespace {

class EnvironmentGuard {
public:
    EnvironmentGuard() : present(getenv("OPENMM_FAST_XML") != NULL), value(present ? getenv("OPENMM_FAST_XML") : "") {
        set("1");
    }
    ~EnvironmentGuard() { set(present ? value.c_str() : NULL); }
private:
    bool present;
    string value;
    static void set(const char* value) {
#ifdef _WIN32
        _putenv_s("OPENMM_FAST_XML", value == NULL ? "" : value);
#else
        if (value == NULL)
            unsetenv("OPENMM_FAST_XML");
        else
            setenv("OPENMM_FAST_XML", value, 1);
#endif
    }
};

template<class T> class ProxyGuard {
public:
    explicit ProxyGuard(const SerializationProxy* replacement) : original(&SerializationProxy::getProxy(typeid(T))) {
        SerializationProxy::registerProxy(typeid(T), replacement);
    }
    ~ProxyGuard() { SerializationProxy::registerProxy(typeid(T), original); }
private:
    const SerializationProxy* original;
};

class CountingStateProxy : public StateProxy {
public:
    mutable int calls = 0;
    void* deserialize(const SerializationNode& node) const override {
        ++calls;
        return StateProxy::deserialize(node);
    }
};

class MarkerSystemProxy : public SystemProxy {
public:
    mutable int calls = 0;
    void serialize(const void* object, SerializationNode& node) const override {
        SystemProxy::serialize(object, node);
        node.createChildNode("CustomMetadata").setStringProperty("marker", "replacement proxy");
    }
    void* deserialize(const SerializationNode& node) const override {
        ++calls;
        ASSERT_EQUAL(string("replacement proxy"), node.getChildNode("CustomMetadata").getStringProperty("marker"));
        return SystemProxy::deserialize(node);
    }
};

static void testCustomStateReader() {
    State::StateBuilder builder(2.0, 12);
    builder.setPeriodicBoxVectors(Vec3(1, 0, 0), Vec3(0, 1, 0), Vec3(0, 0, 1));
    vector<Vec3> positions(3, Vec3(0.25, 0.5, 0.75));
    builder.setPositions(positions);
    State state = builder.getState();
    stringstream stream;
    XmlSerializer::serialize(&state, "State", stream);
    CountingStateProxy proxy;
    ProxyGuard<State> guard(&proxy);
    unique_ptr<State> copy(XmlSerializer::deserialize<State>(stream));
    ASSERT_EQUAL(1, proxy.calls);
    ASSERT_EQUAL(state.getTime(), copy->getTime());
    ASSERT_EQUAL(positions.size(), copy->getPositions().size());
}

static void testCustomLargeSystemRoundTrip() {
    System system;
    for (int i = 0; i < 100000; ++i)
        system.addParticle(1.0+i%3);
    MarkerSystemProxy proxy;
    ProxyGuard<System> guard(&proxy);
    stringstream stream;
    XmlSerializer::serialize(&system, "System", stream);
    ASSERT(stream.str().find("CustomMetadata") != string::npos);
    ASSERT(stream.str().find("__OpenMMFastSystemCookie") == string::npos);
    unique_ptr<System> copy(XmlSerializer::deserialize<System>(stream));
    ASSERT_EQUAL(1, proxy.calls);
    ASSERT_EQUAL(system.getNumParticles(), copy->getNumParticles());
    ASSERT_EQUAL(system.getParticleMass(99999), copy->getParticleMass(99999));
    unique_ptr<System> clone(XmlSerializer::clone(system));
    ASSERT_EQUAL(2, proxy.calls);
    ASSERT_EQUAL(system.getNumParticles(), clone->getNumParticles());
}

} // namespace

int main() {
    EnvironmentGuard env;
    try {
        testCustomStateReader();
        testCustomLargeSystemRoundTrip();
    }
    catch (const exception& e) {
        cout << "exception: " << e.what() << endl;
        return 1;
    }
    cout << "Done" << endl;
    return 0;
}
