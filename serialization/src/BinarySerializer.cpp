/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * Experimental portable binary serialization.                               *
 * -------------------------------------------------------------------------- */

#include "openmm/serialization/BinarySerializer.h"
#include "openmm/serialization/SerializationNode.h"
#include "openmm/serialization/SerializationProxy.h"
#include "openmm/OpenMMException.h"
#include "openmm/Force.h"
#include "openmm/HarmonicAngleForce.h"
#include "openmm/HarmonicBondForce.h"
#include "openmm/NonbondedForce.h"
#include "openmm/PeriodicTorsionForce.h"
#include "openmm/State.h"
#include "openmm/System.h"
#include "openmm/Vec3.h"
#include "openmm/VirtualSite.h"

#include <algorithm>
#include <cstdint>
#include <cstring>
#include <istream>
#include <limits>
#include <map>
#include <memory>
#include <ostream>
#include <string>
#include <type_traits>
#include <typeinfo>
#include <vector>
#include "irrXML.h"
#include <cstdlib>
#include <sstream>
#include <streambuf>
#include <utility>
#include "openmm/CMMotionRemover.h"
#include "openmm/MonteCarloBarostat.h"

using namespace OpenMM;
using namespace std;
using namespace irr;
using namespace io;

namespace {

static const char MAGIC[8] = {'O','M','M','B','I','N','1','\0'};
static const uint32_t FORMAT_VERSION = 1;
static const uint64_t MAX_CONTAINER_ELEMENTS = UINT64_C(1000000000);
static const uint32_t MAX_STRING_BYTES = 256u*1024u*1024u;
static const int MAX_NODE_DEPTH = 1024;

enum ObjectEncoding {
    GENERIC_NODE = 0,
    NATIVE_BINARY = 1
};

static bool littleEndianHost() {
    const uint16_t x = 1;
    return *reinterpret_cast<const unsigned char*>(&x) == 1;
}

class Writer {
public:
    explicit Writer(ostream& stream) : sink(stream.rdbuf()), buffer(4*1024*1024), used(0) {
        static_assert(sizeof(double) == 8, "BinarySerializer requires 64-bit double");
        static_assert(numeric_limits<double>::is_iec559, "BinarySerializer requires IEEE-754 double");
        if (sink == NULL)
            throw OpenMMException("BinarySerializer: output stream has no stream buffer");
        if (!stream.good())
            throw OpenMMException("BinarySerializer: output stream is not writable");
    }

    void flush() {
        if (used != 0) {
            writeExact(reinterpret_cast<const char*>(&buffer[0]), used);
            used = 0;
        }
    }

    void bytes(const void* data, size_t size) {
        if (size == 0)
            return;
        const unsigned char* src = reinterpret_cast<const unsigned char*>(data);

        // Large payloads (positions/velocities/etc.) bypass the staging buffer.
        if (size >= buffer.size()) {
            flush();
            writeExact(reinterpret_cast<const char*>(src), size);
            return;
        }
        if (used+size > buffer.size())
            flush();
        memcpy(&buffer[used], src, size);
        used += size;
    }

    void u8(uint8_t value) {
        bytes(&value, 1);
    }

    void boolean(bool value) {
        u8(value ? 1 : 0);
    }

    void u32(uint32_t value) {
        unsigned char b[4];
        b[0] = static_cast<unsigned char>(value);
        b[1] = static_cast<unsigned char>(value >> 8);
        b[2] = static_cast<unsigned char>(value >> 16);
        b[3] = static_cast<unsigned char>(value >> 24);
        bytes(b, sizeof(b));
    }

    void i32(int32_t value) {
        u32(static_cast<uint32_t>(value));
    }

    void u64(uint64_t value) {
        unsigned char b[8];
        for (int i = 0; i < 8; i++)
            b[i] = static_cast<unsigned char>(value >> (8*i));
        bytes(b, sizeof(b));
    }

    void i64(int64_t value) {
        u64(static_cast<uint64_t>(value));
    }

    void real(double value) {
        uint64_t bits;
        memcpy(&bits, &value, sizeof(bits));
        u64(bits);
    }

    void stringValue(const string& value) {
        if (value.size() > numeric_limits<uint32_t>::max())
            throw OpenMMException("BinarySerializer: string is too large");
        u32(static_cast<uint32_t>(value.size()));
        bytes(value.data(), value.size());
    }

    void vec3Array(const vector<Vec3>& values) {
        u64(static_cast<uint64_t>(values.size()));
        if (values.empty())
            return;
        if (values.size() > MAX_CONTAINER_ELEMENTS)
            throw OpenMMException("BinarySerializer: Vec3 array is too large");
        static_assert(is_standard_layout<Vec3>::value, "Vec3 must be standard layout");
        static_assert(sizeof(Vec3) == 3*sizeof(double), "Unexpected Vec3 layout");
        if (littleEndianHost()) {
            bytes(&values[0], values.size()*sizeof(Vec3));
        }
        else {
            for (size_t i = 0; i < values.size(); i++) {
                real(values[i][0]);
                real(values[i][1]);
                real(values[i][2]);
            }
        }
    }

private:
    streambuf* sink;
    vector<unsigned char> buffer;
    size_t used;

    void writeExact(const char* data, size_t size) {
        size_t written = 0;
        while (written < size) {
            const size_t remaining = size-written;
            const streamsize request = static_cast<streamsize>(remaining);
            streamsize count;
            try {
                count = sink->sputn(data+written, request);
            }
            catch (const ios_base::failure&) {
                throw OpenMMException("BinarySerializer: error writing stream");
            }
            if (count <= 0)
                throw OpenMMException("BinarySerializer: error writing stream");
            written += static_cast<size_t>(count);
        }
    }
};

class Reader {
public:
    explicit Reader(istream& stream) :
            source(stream.rdbuf()), buffer(4*1024*1024), pos(0), end(0), buffered(canSeek(source)) {
        static_assert(sizeof(double) == 8, "BinarySerializer requires 64-bit double");
        static_assert(numeric_limits<double>::is_iec559, "BinarySerializer requires IEEE-754 double");
        if (source == NULL)
            throw OpenMMException("BinarySerializer: input stream has no stream buffer");
        if (!stream.good())
            throw OpenMMException("BinarySerializer: input stream is not readable");
    }

    void bytes(void* data, size_t size) {
        if (!buffered) {
            readExact(data, size);
            return;
        }

        unsigned char* dst = reinterpret_cast<unsigned char*>(data);
        while (size != 0) {
            size_t available = end-pos;
            if (available != 0) {
                size_t take = std::min(size, available);
                memcpy(dst, &buffer[pos], take);
                dst += take;
                pos += take;
                size -= take;
                continue;
            }

            // Large reads bypass the staging buffer after already-buffered bytes
            // have been consumed.  They read exactly the requested payload and
            // therefore cannot consume bytes belonging to a following object.
            if (size >= buffer.size()) {
                readExact(dst, size);
                return;
            }
            refill();
        }
    }

    uint8_t u8() {
        if (!buffered) {
            uint8_t value;
            readExact(&value, 1);
            return value;
        }
        if (pos == end)
            refill();
        return buffer[pos++];
    }

    bool boolean() {
        uint8_t value = u8();
        if (value > 1)
            throw OpenMMException("BinarySerializer: invalid boolean value");
        return value != 0;
    }

    uint32_t u32() {
        unsigned char b[4];
        bytes(b, sizeof(b));
        return static_cast<uint32_t>(b[0]) |
               (static_cast<uint32_t>(b[1]) << 8) |
               (static_cast<uint32_t>(b[2]) << 16) |
               (static_cast<uint32_t>(b[3]) << 24);
    }

    int32_t i32() {
        return static_cast<int32_t>(u32());
    }

    uint64_t u64() {
        unsigned char b[8];
        bytes(b, sizeof(b));
        uint64_t value = 0;
        for (int i = 0; i < 8; i++)
            value |= static_cast<uint64_t>(b[i]) << (8*i);
        return value;
    }

    int64_t i64() {
        return static_cast<int64_t>(u64());
    }

    double real() {
        uint64_t bits = u64();
        double value;
        memcpy(&value, &bits, sizeof(value));
        return value;
    }

    string stringValue() {
        uint32_t size = u32();
        if (size > MAX_STRING_BYTES)
            throw OpenMMException("BinarySerializer: unreasonable string size");
        string value(size, '\0');
        if (size != 0)
            bytes(&value[0], size);
        return value;
    }

    uint64_t count() {
        uint64_t value = u64();
        if (value > MAX_CONTAINER_ELEMENTS)
            throw OpenMMException("BinarySerializer: unreasonable element count");
        return value;
    }

    vector<Vec3> vec3Array() {
        uint64_t n = count();
        if (n > static_cast<uint64_t>(numeric_limits<size_t>::max()/sizeof(Vec3)))
            throw OpenMMException("BinarySerializer: Vec3 array size overflow");
        vector<Vec3> values(static_cast<size_t>(n));
        if (n == 0)
            return values;
        static_assert(is_standard_layout<Vec3>::value, "Vec3 must be standard layout");
        static_assert(sizeof(Vec3) == 3*sizeof(double), "Unexpected Vec3 layout");
        if (littleEndianHost()) {
            bytes(&values[0], static_cast<size_t>(n)*sizeof(Vec3));
        }
        else {
            for (uint64_t i = 0; i < n; i++) {
                double x = real();
                double y = real();
                double z = real();
                values[static_cast<size_t>(i)] = Vec3(x, y, z);
            }
        }
        return values;
    }

    // A buffered Reader may have prefetched bytes belonging to the next
    // serialized object.  Restore the underlying stream to the exact logical
    // end of the object before the Reader is discarded.  Non-seekable streams
    // use exact reads from the outset and therefore need no restoration.
    void finish() {
        if (!buffered || end == pos)
            return;

        const streamoff unread = static_cast<streamoff>(end-pos);
        streampos restored;
        try {
            restored = source->pubseekoff(-unread, ios_base::cur, ios_base::in);
        }
        catch (const ios_base::failure&) {
            throw OpenMMException("BinarySerializer: failed to restore input stream position");
        }
        if (restored == streampos(streamoff(-1)))
            throw OpenMMException("BinarySerializer: failed to restore input stream position");
        pos = end = 0;
    }

private:
    streambuf* source;
    vector<unsigned char> buffer;
    size_t pos, end;
    bool buffered;

    static bool canSeek(streambuf* source) {
        if (source == NULL)
            return false;
        try {
            return source->pubseekoff(0, ios_base::cur, ios_base::in) != streampos(streamoff(-1));
        }
        catch (const ios_base::failure&) {
            return false;
        }
    }

    void readExact(void* data, size_t size) {
        unsigned char* dst = reinterpret_cast<unsigned char*>(data);
        size_t received = 0;
        while (received < size) {
            const size_t remaining = size-received;
            streamsize got;
            try {
                got = source->sgetn(reinterpret_cast<char*>(dst+received), static_cast<streamsize>(remaining));
            }
            catch (const ios_base::failure&) {
                throw OpenMMException("BinarySerializer: truncated or unreadable stream");
            }
            if (got <= 0)
                throw OpenMMException("BinarySerializer: truncated or unreadable stream");
            received += static_cast<size_t>(got);
        }
    }

    void refill() {
        streamsize got;
        try {
            got = source->sgetn(reinterpret_cast<char*>(&buffer[0]), static_cast<streamsize>(buffer.size()));
        }
        catch (const ios_base::failure&) {
            throw OpenMMException("BinarySerializer: unexpected end of stream");
        }
        if (got <= 0)
            throw OpenMMException("BinarySerializer: unexpected end of stream");
        pos = 0;
        end = static_cast<size_t>(got);
    }
};

static void writeStringDoubleMap(Writer& out, const map<string, double>& values) {
    out.u64(static_cast<uint64_t>(values.size()));
    for (map<string, double>::const_iterator it = values.begin(); it != values.end(); ++it) {
        out.stringValue(it->first);
        out.real(it->second);
    }
}

static map<string, double> readStringDoubleMap(Reader& in) {
    uint64_t n = in.count();
    map<string, double> values;
    for (uint64_t i = 0; i < n; i++) {
        string name = in.stringValue();
        double value = in.real();
        values[name] = value;
    }
    return values;
}

static void writeNode(Writer& out, const SerializationNode& node, int depth) {
    if (depth > MAX_NODE_DEPTH)
        throw OpenMMException("BinarySerializer: SerializationNode tree is too deep");
    out.stringValue(node.getName());
    const map<string, string> properties = node.getProperties();
    out.u64(static_cast<uint64_t>(properties.size()));
    for (map<string, string>::const_iterator it = properties.begin(); it != properties.end(); ++it) {
        out.stringValue(it->first);
        out.stringValue(it->second);
    }
    const vector<SerializationNode>& children = node.getChildren();
    out.u64(static_cast<uint64_t>(children.size()));
    for (size_t i = 0; i < children.size(); i++)
        writeNode(out, children[i], depth+1);
}

static SerializationNode readNode(Reader& in, int depth) {
    if (depth > MAX_NODE_DEPTH)
        throw OpenMMException("BinarySerializer: SerializationNode tree is too deep");
    SerializationNode node;
    node.setName(in.stringValue());
    uint64_t numProperties = in.count();
    for (uint64_t i = 0; i < numProperties; i++) {
        string name = in.stringValue();
        string value = in.stringValue();
        node.setStringProperty(name, value);
    }
    uint64_t numChildren = in.count();
    vector<SerializationNode>& children = node.getChildren();
    if (numChildren > static_cast<uint64_t>(numeric_limits<size_t>::max()))
        throw OpenMMException("BinarySerializer: child count overflow");
    children.reserve(static_cast<size_t>(numChildren));
    for (uint64_t i = 0; i < numChildren; i++)
        children.push_back(readNode(in, depth+1));
    return node;
}

class FastProxy {
public:
    FastProxy(const string& stableName, const type_info& cppType) : stableName(stableName), cppTypeName(cppType.name()) {}
    virtual ~FastProxy() {}
    const string stableName;
    const string cppTypeName;
    virtual bool canSerialize(const void* object) const { return true; }
    virtual void serialize(const void* object, Writer& out) const = 0;
    virtual void* deserialize(Reader& in) const = 0;
};

static map<string, const FastProxy*>& proxiesByCppType() {
    static map<string, const FastProxy*> value;
    return value;
}

static map<string, const FastProxy*>& proxiesByStableName() {
    static map<string, const FastProxy*> value;
    return value;
}

static void addProxy(const FastProxy& proxy) {
    proxiesByCppType()[proxy.cppTypeName] = &proxy;
    proxiesByStableName()[proxy.stableName] = &proxy;
}

static void writeObject(Writer& out, const void* object, const type_info& type);
static void* readObject(Reader& in);

class StateFastProxy : public FastProxy {
public:
    StateFastProxy() : FastProxy("State", typeid(State)) {}

    bool canSerialize(const void* object) const {
        const State& state = *reinterpret_cast<const State*>(object);
        // IntegratorParameters are intentionally left to the generic fallback:
        // State does not expose them publicly, and duplicating private access in
        // a first binary-format patch would make the implementation more invasive.
        return (state.getDataTypes() & State::IntegratorParameters) == 0;
    }

    void serialize(const void* object, Writer& out) const {
        const State& state = *reinterpret_cast<const State*>(object);
        out.u32(1); // State binary schema version.
        int types = state.getDataTypes();
        out.i32(types);
        out.real(state.getTime());
        out.i64(state.getStepCount());

        Vec3 a, b, c;
        state.getPeriodicBoxVectors(a, b, c);
        out.real(a[0]); out.real(a[1]); out.real(a[2]);
        out.real(b[0]); out.real(b[1]); out.real(b[2]);
        out.real(c[0]); out.real(c[1]); out.real(c[2]);

        if (types & State::Parameters)
            writeStringDoubleMap(out, state.getParameters());
        if (types & State::ParameterDerivatives)
            writeStringDoubleMap(out, state.getEnergyParameterDerivatives());
        if (types & State::Energy) {
            out.real(state.getKineticEnergy());
            out.real(state.getPotentialEnergy());
        }
        if (types & State::Positions)
            out.vec3Array(state.getPositions());
        if (types & State::Velocities)
            out.vec3Array(state.getVelocities());
        if (types & State::Forces)
            out.vec3Array(state.getForces());
    }

    void* deserialize(Reader& in) const {
        uint32_t version = in.u32();
        if (version != 1)
            throw OpenMMException("BinarySerializer: unsupported State schema version");
        int types = in.i32();
        double time = in.real();
        long long stepCount = static_cast<long long>(in.i64());
        State::StateBuilder builder(time, stepCount);
        double ax = in.real(); double ay = in.real(); double az = in.real();
        double bx = in.real(); double by = in.real(); double bz = in.real();
        double cx = in.real(); double cy = in.real(); double cz = in.real();
        Vec3 a(ax, ay, az), b(bx, by, bz), c(cx, cy, cz);
        builder.setPeriodicBoxVectors(a, b, c);

        if (types & State::Parameters)
            builder.setParameters(readStringDoubleMap(in));
        if (types & State::ParameterDerivatives)
            builder.setEnergyParameterDerivatives(readStringDoubleMap(in));
        if (types & State::Energy) {
            double kinetic = in.real();
            double potential = in.real();
            builder.setEnergy(kinetic, potential);
        }
        if (types & State::Positions)
            builder.setPositions(in.vec3Array());
        if (types & State::Velocities)
            builder.setVelocities(in.vec3Array());
        if (types & State::Forces)
            builder.setForces(in.vec3Array());
        State* state = new State();\n        *state = builder.getState();\n        return state;
    }
};

enum VirtualSiteCode {
    VS_TWO_PARTICLE = 1,
    VS_THREE_PARTICLE = 2,
    VS_OUT_OF_PLANE = 3,
    VS_LOCAL_COORDINATES = 4
};

class SystemFastProxy : public FastProxy {
public:
    SystemFastProxy() : FastProxy("System", typeid(System)) {}

    bool canSerialize(const void* object) const {
        const System& system = *reinterpret_cast<const System*>(object);
        for (int i = 0; i < system.getNumParticles(); i++) {
            if (!system.isVirtualSite(i))
                continue;
            const VirtualSite& site = system.getVirtualSite(i);
            const type_info& t = typeid(site);
            if (t != typeid(TwoParticleAverageSite) &&
                t != typeid(ThreeParticleAverageSite) &&
                t != typeid(OutOfPlaneSite) &&
                t != typeid(LocalCoordinatesSite))
                return false;
        }
        return true;
    }

    void serialize(const void* object, Writer& out) const {
        const System& system = *reinterpret_cast<const System*>(object);
        out.u32(1); // System schema version.

        Vec3 a, b, c;
        system.getDefaultPeriodicBoxVectors(a, b, c);
        out.real(a[0]); out.real(a[1]); out.real(a[2]);
        out.real(b[0]); out.real(b[1]); out.real(b[2]);
        out.real(c[0]); out.real(c[1]); out.real(c[2]);

        const int numParticles = system.getNumParticles();
        out.u64(static_cast<uint64_t>(numParticles));
        // Masses are written as one canonical fixed-width array.  Avoid making
        // an 8M-node tree simply to represent 8M doubles.
        for (int i = 0; i < numParticles; i++)
            out.real(system.getParticleMass(i));

        uint64_t numVirtualSites = 0;
        for (int i = 0; i < numParticles; i++)
            if (system.isVirtualSite(i))
                numVirtualSites++;
        out.u64(numVirtualSites);
        for (int i = 0; i < numParticles; i++) {
            if (!system.isVirtualSite(i))
                continue;
            out.u32(static_cast<uint32_t>(i));
            const VirtualSite& site = system.getVirtualSite(i);
            if (typeid(site) == typeid(TwoParticleAverageSite)) {
                const TwoParticleAverageSite& s = dynamic_cast<const TwoParticleAverageSite&>(site);
                out.u8(VS_TWO_PARTICLE);
                out.i32(s.getParticle(0)); out.i32(s.getParticle(1));
                out.real(s.getWeight(0)); out.real(s.getWeight(1));
            }
            else if (typeid(site) == typeid(ThreeParticleAverageSite)) {
                const ThreeParticleAverageSite& s = dynamic_cast<const ThreeParticleAverageSite&>(site);
                out.u8(VS_THREE_PARTICLE);
                out.i32(s.getParticle(0)); out.i32(s.getParticle(1)); out.i32(s.getParticle(2));
                out.real(s.getWeight(0)); out.real(s.getWeight(1)); out.real(s.getWeight(2));
            }
            else if (typeid(site) == typeid(OutOfPlaneSite)) {
                const OutOfPlaneSite& s = dynamic_cast<const OutOfPlaneSite&>(site);
                out.u8(VS_OUT_OF_PLANE);
                out.i32(s.getParticle(0)); out.i32(s.getParticle(1)); out.i32(s.getParticle(2));
                out.real(s.getWeight12()); out.real(s.getWeight13()); out.real(s.getWeightCross());
            }
            else if (typeid(site) == typeid(LocalCoordinatesSite)) {
                const LocalCoordinatesSite& s = dynamic_cast<const LocalCoordinatesSite&>(site);
                out.u8(VS_LOCAL_COORDINATES);
                int n = s.getNumParticles();
                out.u32(static_cast<uint32_t>(n));
                vector<double> wo, wx, wy;
                s.getOriginWeights(wo);
                s.getXWeights(wx);
                s.getYWeights(wy);
                for (int j = 0; j < n; j++) {
                    out.i32(s.getParticle(j));
                    out.real(wo[j]); out.real(wx[j]); out.real(wy[j]);
                }
                Vec3 p = s.getLocalPosition();
                out.real(p[0]); out.real(p[1]); out.real(p[2]);
            }
            else {
                throw OpenMMException("BinarySerializer: unsupported virtual site type");
            }
        }

        const int numConstraints = system.getNumConstraints();
        out.u64(static_cast<uint64_t>(numConstraints));
        for (int i = 0; i < numConstraints; i++) {
            int p1, p2;
            double distance;
            system.getConstraintParameters(i, p1, p2, distance);
            out.i32(p1); out.i32(p2); out.real(distance);
        }

        const int numForces = system.getNumForces();
        out.u64(static_cast<uint64_t>(numForces));
        for (int i = 0; i < numForces; i++) {
            const Force& force = system.getForce(i);
            writeObject(out, &force, typeid(force));
        }
    }

    void* deserialize(Reader& in) const {
        uint32_t version = in.u32();
        if (version != 1)
            throw OpenMMException("BinarySerializer: unsupported System schema version");
        unique_ptr<System> system(new System());
        double ax = in.real(); double ay = in.real(); double az = in.real();
        double bx = in.real(); double by = in.real(); double bz = in.real();
        double cx = in.real(); double cy = in.real(); double cz = in.real();
        Vec3 a(ax, ay, az), b(bx, by, bz), c(cx, cy, cz);
        system->setDefaultPeriodicBoxVectors(a, b, c);

        uint64_t numParticles = in.count();
        for (uint64_t i = 0; i < numParticles; i++)
            system->addParticle(in.real());

        uint64_t numVirtualSites = in.count();
        for (uint64_t i = 0; i < numVirtualSites; i++) {
            uint32_t particle = in.u32();
            if (particle >= numParticles)
                throw OpenMMException("BinarySerializer: virtual site particle index out of range");
            uint8_t code = in.u8();
            if (code == VS_TWO_PARTICLE) {
                int p1 = in.i32(); int p2 = in.i32();
                double w1 = in.real(); double w2 = in.real();
                system->setVirtualSite(particle, new TwoParticleAverageSite(p1, p2, w1, w2));
            }
            else if (code == VS_THREE_PARTICLE) {
                int p1 = in.i32(); int p2 = in.i32(); int p3 = in.i32();
                double w1 = in.real(); double w2 = in.real(); double w3 = in.real();
                system->setVirtualSite(particle, new ThreeParticleAverageSite(p1, p2, p3, w1, w2, w3));
            }
            else if (code == VS_OUT_OF_PLANE) {
                int p1 = in.i32(); int p2 = in.i32(); int p3 = in.i32();
                double w12 = in.real(); double w13 = in.real(); double wc = in.real();
                system->setVirtualSite(particle, new OutOfPlaneSite(p1, p2, p3, w12, w13, wc));
            }
            else if (code == VS_LOCAL_COORDINATES) {
                uint32_t n = in.u32();
                vector<int> particles(n);
                vector<double> wo(n), wx(n), wy(n);
                for (uint32_t j = 0; j < n; j++) {
                    particles[j] = in.i32();
                    wo[j] = in.real(); wx[j] = in.real(); wy[j] = in.real();
                }
                double px = in.real(); double py = in.real(); double pz = in.real();
                Vec3 p(px, py, pz);
                system->setVirtualSite(particle, new LocalCoordinatesSite(particles, wo, wx, wy, p));
            }
            else
                throw OpenMMException("BinarySerializer: unknown virtual site code");
        }

        uint64_t numConstraints = in.count();
        for (uint64_t i = 0; i < numConstraints; i++) {
            int p1 = in.i32();
            int p2 = in.i32();
            double distance = in.real();
            system->addConstraint(p1, p2, distance);
        }

        uint64_t numForces = in.count();
        for (uint64_t i = 0; i < numForces; i++)
            system->addForce(reinterpret_cast<Force*>(readObject(in)));
        return system.release();
    }
};

class HarmonicBondFastProxy : public FastProxy {
public:
    HarmonicBondFastProxy() : FastProxy("HarmonicBondForce", typeid(HarmonicBondForce)) {}
    void serialize(const void* object, Writer& out) const {
        const HarmonicBondForce& force = *reinterpret_cast<const HarmonicBondForce*>(object);
        out.u32(1);
        out.i32(force.getForceGroup());
        out.stringValue(force.getName());
        out.boolean(force.usesPeriodicBoundaryConditions());
        out.u64(static_cast<uint64_t>(force.getNumBonds()));
        for (int i = 0; i < force.getNumBonds(); i++) {
            int p1, p2; double d, k;
            force.getBondParameters(i, p1, p2, d, k);
            out.i32(p1); out.i32(p2); out.real(d); out.real(k);
        }
    }
    void* deserialize(Reader& in) const {
        if (in.u32() != 1)
            throw OpenMMException("BinarySerializer: unsupported HarmonicBondForce schema version");
        unique_ptr<HarmonicBondForce> force(new HarmonicBondForce());
        force->setForceGroup(in.i32());
        force->setName(in.stringValue());
        force->setUsesPeriodicBoundaryConditions(in.boolean());
        uint64_t n = in.count();
        for (uint64_t i = 0; i < n; i++) {
            int p1 = in.i32();
            int p2 = in.i32();
            double d = in.real();
            double k = in.real();
            force->addBond(p1, p2, d, k);
        }
        return force.release();
    }
};

class HarmonicAngleFastProxy : public FastProxy {
public:
    HarmonicAngleFastProxy() : FastProxy("HarmonicAngleForce", typeid(HarmonicAngleForce)) {}
    void serialize(const void* object, Writer& out) const {
        const HarmonicAngleForce& force = *reinterpret_cast<const HarmonicAngleForce*>(object);
        out.u32(1);
        out.i32(force.getForceGroup());
        out.stringValue(force.getName());
        out.boolean(force.usesPeriodicBoundaryConditions());
        out.u64(static_cast<uint64_t>(force.getNumAngles()));
        for (int i = 0; i < force.getNumAngles(); i++) {
            int p1, p2, p3; double a, k;
            force.getAngleParameters(i, p1, p2, p3, a, k);
            out.i32(p1); out.i32(p2); out.i32(p3); out.real(a); out.real(k);
        }
    }
    void* deserialize(Reader& in) const {
        if (in.u32() != 1)
            throw OpenMMException("BinarySerializer: unsupported HarmonicAngleForce schema version");
        unique_ptr<HarmonicAngleForce> force(new HarmonicAngleForce());
        force->setForceGroup(in.i32());
        force->setName(in.stringValue());
        force->setUsesPeriodicBoundaryConditions(in.boolean());
        uint64_t n = in.count();
        for (uint64_t i = 0; i < n; i++) {
            int p1 = in.i32();
            int p2 = in.i32();
            int p3 = in.i32();
            double angle = in.real();
            double k = in.real();
            force->addAngle(p1, p2, p3, angle, k);
        }
        return force.release();
    }
};

class PeriodicTorsionFastProxy : public FastProxy {
public:
    PeriodicTorsionFastProxy() : FastProxy("PeriodicTorsionForce", typeid(PeriodicTorsionForce)) {}
    void serialize(const void* object, Writer& out) const {
        const PeriodicTorsionForce& force = *reinterpret_cast<const PeriodicTorsionForce*>(object);
        out.u32(1);
        out.i32(force.getForceGroup());
        out.stringValue(force.getName());
        out.boolean(force.usesPeriodicBoundaryConditions());
        out.u64(static_cast<uint64_t>(force.getNumTorsions()));
        for (int i = 0; i < force.getNumTorsions(); i++) {
            int p1, p2, p3, p4, periodicity; double phase, k;
            force.getTorsionParameters(i, p1, p2, p3, p4, periodicity, phase, k);
            out.i32(p1); out.i32(p2); out.i32(p3); out.i32(p4);
            out.i32(periodicity); out.real(phase); out.real(k);
        }
    }
    void* deserialize(Reader& in) const {
        if (in.u32() != 1)
            throw OpenMMException("BinarySerializer: unsupported PeriodicTorsionForce schema version");
        unique_ptr<PeriodicTorsionForce> force(new PeriodicTorsionForce());
        force->setForceGroup(in.i32());
        force->setName(in.stringValue());
        force->setUsesPeriodicBoundaryConditions(in.boolean());
        uint64_t n = in.count();
        for (uint64_t i = 0; i < n; i++) {
            int p1 = in.i32(); int p2 = in.i32(); int p3 = in.i32(); int p4 = in.i32();
            int periodicity = in.i32();
            double phase = in.real(); double k = in.real();
            force->addTorsion(p1, p2, p3, p4, periodicity, phase, k);
        }
        return force.release();
    }
};

class NonbondedFastProxy : public FastProxy {
public:
    NonbondedFastProxy() : FastProxy("NonbondedForce", typeid(NonbondedForce)) {}

    void serialize(const void* object, Writer& out) const {
        const NonbondedForce& force = *reinterpret_cast<const NonbondedForce*>(object);
        out.u32(1);
        out.i32(force.getForceGroup());
        out.stringValue(force.getName());
        out.i32(static_cast<int>(force.getNonbondedMethod()));
        out.real(force.getCutoffDistance());
        out.boolean(force.getUseSwitchingFunction());
        out.real(force.getSwitchingDistance());
        out.real(force.getEwaldErrorTolerance());
        out.real(force.getReactionFieldDielectric());
        out.boolean(force.getUseDispersionCorrection());
        out.boolean(force.getExceptionsUsePeriodicBoundaryConditions());
        out.boolean(force.getIncludeDirectSpace());

        double alpha; int nx, ny, nz;
        force.getPMEParameters(alpha, nx, ny, nz);
        out.real(alpha); out.i32(nx); out.i32(ny); out.i32(nz);
        force.getLJPMEParameters(alpha, nx, ny, nz);
        out.real(alpha); out.i32(nx); out.i32(ny); out.i32(nz);
        out.i32(force.getReciprocalSpaceForceGroup());

        out.u64(static_cast<uint64_t>(force.getNumGlobalParameters()));
        for (int i = 0; i < force.getNumGlobalParameters(); i++) {
            out.stringValue(force.getGlobalParameterName(i));
            out.real(force.getGlobalParameterDefaultValue(i));
        }

        out.u64(static_cast<uint64_t>(force.getNumParticleParameterOffsets()));
        for (int i = 0; i < force.getNumParticleParameterOffsets(); i++) {
            string parameter; int particle; double q, sig, eps;
            force.getParticleParameterOffset(i, parameter, particle, q, sig, eps);
            out.stringValue(parameter); out.i32(particle); out.real(q); out.real(sig); out.real(eps);
        }

        out.u64(static_cast<uint64_t>(force.getNumExceptionParameterOffsets()));
        for (int i = 0; i < force.getNumExceptionParameterOffsets(); i++) {
            string parameter; int exception; double q, sig, eps;
            force.getExceptionParameterOffset(i, parameter, exception, q, sig, eps);
            out.stringValue(parameter); out.i32(exception); out.real(q); out.real(sig); out.real(eps);
        }

        // Fixed-size records: no node allocation, no decimal conversion, no
        // repeated XML attribute names.  At 7.95M particles this is ~191 MB.
        out.u64(static_cast<uint64_t>(force.getNumParticles()));
        for (int i = 0; i < force.getNumParticles(); i++) {
            double q, sig, eps;
            force.getParticleParameters(i, q, sig, eps);
            out.real(q); out.real(sig); out.real(eps);
        }

        out.u64(static_cast<uint64_t>(force.getNumExceptions()));
        for (int i = 0; i < force.getNumExceptions(); i++) {
            int p1, p2; double q, sig, eps;
            force.getExceptionParameters(i, p1, p2, q, sig, eps);
            out.i32(p1); out.i32(p2); out.real(q); out.real(sig); out.real(eps);
        }
    }

    void* deserialize(Reader& in) const {
        if (in.u32() != 1)
            throw OpenMMException("BinarySerializer: unsupported NonbondedForce schema version");
        unique_ptr<NonbondedForce> force(new NonbondedForce());
        force->setForceGroup(in.i32());
        force->setName(in.stringValue());
        force->setNonbondedMethod(static_cast<NonbondedForce::NonbondedMethod>(in.i32()));
        force->setCutoffDistance(in.real());
        force->setUseSwitchingFunction(in.boolean());
        force->setSwitchingDistance(in.real());
        force->setEwaldErrorTolerance(in.real());
        force->setReactionFieldDielectric(in.real());
        force->setUseDispersionCorrection(in.boolean());
        force->setExceptionsUsePeriodicBoundaryConditions(in.boolean());
        force->setIncludeDirectSpace(in.boolean());
        double alpha = in.real(); int nx = in.i32(); int ny = in.i32(); int nz = in.i32();
        force->setPMEParameters(alpha, nx, ny, nz);
        alpha = in.real(); nx = in.i32(); ny = in.i32(); nz = in.i32();
        force->setLJPMEParameters(alpha, nx, ny, nz);
        force->setReciprocalSpaceForceGroup(in.i32());

        uint64_t n = in.count();
        for (uint64_t i = 0; i < n; i++) {
            string name = in.stringValue();
            double value = in.real();
            force->addGlobalParameter(name, value);
        }

        n = in.count();
        for (uint64_t i = 0; i < n; i++) {
            string parameter = in.stringValue();
            int particle = in.i32();
            double q = in.real(); double sig = in.real(); double eps = in.real();
            force->addParticleParameterOffset(parameter, particle, q, sig, eps);
        }

        n = in.count();
        for (uint64_t i = 0; i < n; i++) {
            string parameter = in.stringValue();
            int exception = in.i32();
            double q = in.real(); double sig = in.real(); double eps = in.real();
            force->addExceptionParameterOffset(parameter, exception, q, sig, eps);
        }

        n = in.count();
        for (uint64_t i = 0; i < n; i++) {
            double q = in.real();
            double sig = in.real();
            double eps = in.real();
            force->addParticle(q, sig, eps);
        }

        n = in.count();
        for (uint64_t i = 0; i < n; i++) {
            int p1 = in.i32(); int p2 = in.i32();
            double q = in.real(); double sig = in.real(); double eps = in.real();
            force->addException(p1, p2, q, sig, eps);
        }
        return force.release();
    }
};

struct BuiltinProxyRegistry {
    StateFastProxy stateProxy;
    SystemFastProxy systemProxy;
    NonbondedFastProxy nonbondedProxy;
    HarmonicBondFastProxy bondProxy;
    HarmonicAngleFastProxy angleProxy;
    PeriodicTorsionFastProxy torsionProxy;

    BuiltinProxyRegistry() {
        addProxy(stateProxy);
        addProxy(systemProxy);
        addProxy(nonbondedProxy);
        addProxy(bondProxy);
        addProxy(angleProxy);
        addProxy(torsionProxy);
    }
};

static void ensureBuiltinsRegistered() {
    // Function-local static initialization is thread-safe in C++11.
    static BuiltinProxyRegistry registry;
    (void) registry;
}

static void writeObject(Writer& out, const void* object, const type_info& type) {
    ensureBuiltinsRegistered();
    map<string, const FastProxy*>::const_iterator fast = proxiesByCppType().find(type.name());
    if (fast != proxiesByCppType().end() && fast->second->canSerialize(object)) {
        out.u8(NATIVE_BINARY);
        out.stringValue(fast->second->stableName);
        fast->second->serialize(object, out);
        return;
    }

    // Generic compatibility path.  It still removes XML parsing/formatting and
    // works for every currently registered OpenMM or plugin SerializationProxy.
    const SerializationProxy& proxy = SerializationProxy::getProxy(type);
    out.u8(GENERIC_NODE);
    out.stringValue(proxy.getTypeName());
    SerializationNode node;
    node.setName(proxy.getTypeName());
    proxy.serialize(object, node);
    writeNode(out, node, 0);
}

static void* readObject(Reader& in) {
    ensureBuiltinsRegistered();
    uint8_t encoding = in.u8();
    string typeName = in.stringValue();
    if (encoding == NATIVE_BINARY) {
        map<string, const FastProxy*>::const_iterator proxy = proxiesByStableName().find(typeName);
        if (proxy == proxiesByStableName().end())
            throw OpenMMException("BinarySerializer: no native proxy for type '"+typeName+"'");
        return proxy->second->deserialize(in);
    }
    if (encoding == GENERIC_NODE) {
        SerializationNode node = readNode(in, 0);
        return SerializationProxy::getProxy(typeName).deserialize(node);
    }
    throw OpenMMException("BinarySerializer: unknown object encoding");
}

} // anonymous namespace

void BinarySerializer::serializeImpl(const void* object, const type_info& type, ostream& stream) {
    Writer out(stream);
    out.bytes(MAGIC, sizeof(MAGIC));
    out.u32(FORMAT_VERSION);
    writeObject(out, object, type);
    out.flush();
}

void* BinarySerializer::deserializeImpl(istream& stream) {
    Reader in(stream);
    char magic[sizeof(MAGIC)];
    in.bytes(magic, sizeof(magic));
    if (memcmp(magic, MAGIC, sizeof(MAGIC)) != 0)
        throw OpenMMException("BinarySerializer: invalid file signature");
    uint32_t version = in.u32();
    if (version != FORMAT_VERSION)
        throw OpenMMException("BinarySerializer: unsupported file format version");

    void* object = readObject(in);
    in.finish();
    return object;
}

