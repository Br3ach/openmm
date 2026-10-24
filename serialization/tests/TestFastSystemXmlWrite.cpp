/* Regression coverage for the bounded-memory System XML writer. */
#include "openmm/internal/AssertionUtilities.h"
#include "openmm/CMMotionRemover.h"
#include "openmm/CustomBondForce.h"
#include "openmm/HarmonicAngleForce.h"
#include "openmm/HarmonicBondForce.h"
#include "openmm/MonteCarloBarostat.h"
#include "openmm/NonbondedForce.h"
#include "openmm/PeriodicTorsionForce.h"
#include "openmm/System.h"
#include "openmm/VirtualSite.h"
#include "openmm/serialization/NonbondedForceProxy.h"
#include "openmm/serialization/SystemProxy.h"
#include "openmm/serialization/XmlSerializer.h"
#include <cstdlib>
#include <iostream>
#include <memory>
#include <locale>
#include <sstream>
#include <streambuf>
#include <vector>

using namespace OpenMM;
using namespace std;

static void setFlag(const char* name, const char* value) {
#ifdef _WIN32
    _putenv_s(name, value);
#else
    setenv(name, value, 1);
#endif
}

static string serializeSystem(const System& system, bool fast, const char* root="System") {
    setFlag("OPENMM_FAST_XML", fast ? "1" : "0");
    stringstream stream;
    XmlSerializer::serialize(&system, root, stream);
    return stream.str();
}

static unique_ptr<System> makeSystem(bool withForces) {
    unique_ptr<System> system(new System());
    for (int i = 0; i < 100000; ++i)
        system->addParticle(i%3 == 0 ? 16.0 : 1.0);
    system->setDefaultPeriodicBoxVectors(Vec3(4,0,0), Vec3(0.1,5,0), Vec3(0.2,0.3,6));
    if (!withForces)
        return system;
    system->setVirtualSite(99996, new TwoParticleAverageSite(0,1,0.25,0.75));
    system->setVirtualSite(99997, new ThreeParticleAverageSite(0,1,2,0.2,0.3,0.5));
    system->setVirtualSite(99998, new OutOfPlaneSite(0,1,2,0.2,0.3,0.4));
    vector<int> particles = {0,1,2};
    vector<double> wo = {0.2,0.3,0.5}, wx = {-1,1,0}, wy = {-1,0,1};
    system->setVirtualSite(99999, new LocalCoordinatesSite(particles,wo,wx,wy,Vec3(0.01,0.02,0.03)));
    system->addConstraint(0,1,0.125);
    auto* nb = new NonbondedForce();
    nb->setName("NB & <escaped> \" ' name");
    nb->setForceGroup(3);
    nb->setNonbondedMethod(NonbondedForce::LJPME);
    nb->setUseSwitchingFunction(true);
    nb->setSwitchingDistance(0.75);
    nb->setCutoffDistance(1.2);
    nb->setEwaldErrorTolerance(0.0002);
    nb->setReactionFieldDielectric(70.0);
    nb->setUseDispersionCorrection(false);
    nb->setExceptionsUsePeriodicBoundaryConditions(true);
    nb->setIncludeDirectSpace(false);
    nb->setPMEParameters(0.25,16,17,18);
    nb->setLJPMEParameters(0.5,19,20,21);
    nb->setReciprocalSpaceForceGroup(4);
    nb->addGlobalParameter("lambda&x",0.25);
    for (int i = 0; i < 100000; ++i)
        nb->addParticle((i%3-1)*0.25,0.125,0.0625);
    nb->addException(0,1,0.25,0.5,0.75);
    nb->addParticleParameterOffset("lambda&x",2,0.1,0.2,0.3);
    nb->addExceptionParameterOffset("lambda&x",0,0.4,0.5,0.6);
    system->addForce(nb);
    auto* bonds = new HarmonicBondForce();
    bonds->setUsesPeriodicBoundaryConditions(true); bonds->addBond(0,1,0.125,200);
    system->addForce(bonds);
    auto* angles = new HarmonicAngleForce();
    angles->setUsesPeriodicBoundaryConditions(true); angles->addAngle(0,1,2,1.25,100);
    system->addForce(angles);
    auto* torsions = new PeriodicTorsionForce();
    torsions->setUsesPeriodicBoundaryConditions(true); torsions->addTorsion(0,1,2,3,3,0.5,5);
    system->addForce(torsions);
    // Empty direct-force sections must use the stock self-closing form.
    system->addForce(new NonbondedForce());
    system->addForce(new HarmonicBondForce());
    system->addForce(new HarmonicAngleForce());
    system->addForce(new PeriodicTorsionForce());
    system->addForce(new CMMotionRemover(10));
    system->addForce(new MonteCarloBarostat(1.0,300.0,25));
    auto* custom = new CustomBondForce("k*(r-r0)^2");
    custom->addGlobalParameter("k",2); custom->addGlobalParameter("r0",0.1);
    custom->addBond(0,1); system->addForce(custom);
    return system;
}

static void testCanonicalOutputAndControls() {
    setFlag("OPENMM_FAST_SYSTEM_WRITE","1");
    for (bool withForces : {false,true}) {
        auto system=makeSystem(withForces);
        const string stock=serializeSystem(*system,false);
        ASSERT_EQUAL(stock,serializeSystem(*system,true));
        setFlag("OPENMM_FAST_SYSTEM_WRITE","0");
        ASSERT_EQUAL(stock,serializeSystem(*system,true));
        setFlag("OPENMM_FAST_SYSTEM_WRITE","1");
        ASSERT_EQUAL(serializeSystem(*system,false,"CustomSystem"),serializeSystem(*system,true,"CustomSystem"));
        setFlag("OPENMM_FAST_XML","1");
        unique_ptr<System> clone(XmlSerializer::clone(*system));
        ASSERT_EQUAL(stock,serializeSystem(*clone,true));
    }
}

class FailOutput : public streambuf {
protected:
    streamsize xsputn(const char*, streamsize) override {return 0;}
    int_type overflow(int_type) override {return traits_type::eof();}
};

static void testWriteFailure() {
    auto system=makeSystem(false);
    setFlag("OPENMM_FAST_XML","1");
    setFlag("OPENMM_FAST_SYSTEM_WRITE","1");
    FailOutput buffer;
    ostream stream(&buffer);
    bool threw=false;
    try {XmlSerializer::serialize(system.get(),"System",stream);} catch (const OpenMMException&) {threw=true;}
    ASSERT(threw); ASSERT(stream.bad());
    // A failed write must not affect a subsequent serialization.
    ASSERT_EQUAL(serializeSystem(*system,false),serializeSystem(*system,true));
    setFlag("OPENMM_FAST_XML","1");
    stringstream failed;
    failed.setstate(ios_base::failbit);
    threw=false;
    try {XmlSerializer::serialize(system.get(),"System",failed);} catch (const OpenMMException&) {threw=true;}
    ASSERT(threw);
}

class ModifiedNonbondedProxy : public NonbondedForceProxy {
public:
    void serialize(const void* object, SerializationNode& node) const override {
        NonbondedForceProxy::serialize(object,node);
        node.setStringProperty("custom-marker","preserved");
    }
};

class DelegatingSystemProxy : public SerializationProxy {
public:
    DelegatingSystemProxy() : SerializationProxy("System") {}
    void serialize(const void* object, SerializationNode& node) const override {
        SystemProxy original;
        original.serialize(object,node);
        node.setStringProperty("custom-marker","preserved");
    }
    void* deserialize(const SerializationNode& node) const override {
        SystemProxy original;
        return original.deserialize(node);
    }
};

template<class Object> class ProxyGuard {
public:
    ProxyGuard(const SerializationProxy* replacement) : previous(&SerializationProxy::getProxy(typeid(Object))) {
        SerializationProxy::registerProxy(typeid(Object),replacement);
    }
    ~ProxyGuard() {SerializationProxy::registerProxy(typeid(Object),previous);}
private:
    const SerializationProxy* previous;
};

static void testReplacementProxies() {
    auto system=makeSystem(true);
    ModifiedNonbondedProxy nonbonded;
    {
        ProxyGuard<NonbondedForce> guard(&nonbonded);
        const string stock=serializeSystem(*system,false);
        ASSERT(stock.find("custom-marker=\"preserved\"")!=string::npos);
        ASSERT_EQUAL(stock,serializeSystem(*system,true));
    }
    DelegatingSystemProxy customSystem;
    {
        ProxyGuard<System> guard(&customSystem);
        const string stock=serializeSystem(*system,false);
        ASSERT(stock.find("custom-marker=\"preserved\"")!=string::npos);
        ASSERT_EQUAL(stock,serializeSystem(*system,true));
    }
}

// A plugin can legitimately embed multiple canonical System child nodes.
// The direct writer must never change SystemProxy's generic tree contract.
class NestedSystemForce : public Force {
public:
    unique_ptr<System> first, second;
    NestedSystemForce() : first(new System()), second(new System()) {}
protected:
    ForceImpl* createImpl() const override {
        throw OpenMMException("Serialization-only test force");
    }
};

class NestedSystemProxy : public SerializationProxy {
public:
    NestedSystemProxy() : SerializationProxy("NestedSystemForce") {}
    void serialize(const void* object, SerializationNode& node) const override {
        const NestedSystemForce& force = *reinterpret_cast<const NestedSystemForce*>(object);
        node.createChildNode("System", force.first.get());
        node.createChildNode("System", force.second.get());
    }
    void* deserialize(const SerializationNode& node) const override {
        unique_ptr<NestedSystemForce> force(new NestedSystemForce());
        force->first.reset(node.getChildren()[0].decodeObject<System>());
        force->second.reset(node.getChildren()[1].decodeObject<System>());
        return force.release();
    }
};

static void testDirectProxySnapshot() {
    auto system = makeSystem(false);
    setFlag("OPENMM_FAST_XML", "1");
    SerializationNode node;
    node.setName("System");
    SystemProxy proxy;
    proxy.serialize(system.get(), node);
    system.reset();
    unique_ptr<System> decoded(reinterpret_cast<System*>(proxy.deserialize(node)));
    ASSERT_EQUAL(100000, decoded->getNumParticles());
    ASSERT_EQUAL(16.0, decoded->getParticleMass(0));
    ASSERT_EQUAL(1.0, decoded->getParticleMass(1));
}

static void testNestedSystems() {
    static NestedSystemProxy proxy;
    SerializationProxy::registerProxy(typeid(NestedSystemForce), &proxy);
    auto system = makeSystem(false);
    auto* nested = new NestedSystemForce();
    nested->first = makeSystem(false);
    nested->second = makeSystem(false);
    system->addForce(nested);
    const string stock = serializeSystem(*system, false);
    ASSERT_EQUAL(stock, serializeSystem(*system, true));
    stringstream input(stock);
    unique_ptr<System> decoded(XmlSerializer::deserialize<System>(input));
    const auto& force = dynamic_cast<const NestedSystemForce&>(decoded->getForce(0));
    ASSERT_EQUAL(100000, force.first->getNumParticles());
    ASSERT_EQUAL(100000, force.second->getNumParticles());
}

class GroupedIntegers : public numpunct<char> {
protected:
    string do_grouping() const override {return "\3";}
    char do_thousands_sep() const override {return ',';}
};

static void testGlobalLocale() {
    auto system = makeSystem(false);
    system->addConstraint(1000,2000,0.125);
    const locale previous;
    locale::global(locale(locale::classic(), new GroupedIntegers()));
    string stock, fast;
    try {
        stock = serializeSystem(*system, false);
        fast = serializeSystem(*system, true);
    }
    catch (...) {
        locale::global(previous);
        throw;
    }
    locale::global(previous);
    ASSERT(stock.find("p1=\"1,000\"") != string::npos);
    ASSERT_EQUAL(stock, fast);
}

int main() {
    try {
        testCanonicalOutputAndControls();
        testWriteFailure();
        testReplacementProxies();
        testNestedSystems();
        testDirectProxySnapshot();
        testGlobalLocale();
    } catch (const exception& e) {
        cerr<<e.what()<<endl;return 1;
    }
    cout<<"Done"<<endl;
    return 0;
}
