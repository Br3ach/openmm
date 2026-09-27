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

// OpenMM's internal locale-independent floating-point parser.  Reusing it keeps
// numeric behavior aligned with SerializationNode/stock XML deserialization.
extern "C" double strtod2(const char* s00, char** se);

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

class FastXmlStreamReader : public IFileReadCallBack {
public:
    explicit FastXmlStreamReader(istream& stream) : stream(stream) {
        stream.clear();
        start = stream.tellg();
        if (start == streampos(-1))
            throw OpenMMException("BinarySerializer: XML stream is not seekable");

        stream.seekg(0, ios_base::end);
        streampos end = stream.tellg();
        if (end == streampos(-1) || end < start)
            throw OpenMMException("BinarySerializer: could not determine XML stream size");

        streamoff length = end-start;
        if (length > numeric_limits<int>::max())
            throw OpenMMException("BinarySerializer: XML stream is too large for irrXML");

        size = static_cast<int>(length);
        stream.clear();
        stream.seekg(start);
    }

    int read(void* buffer, int sizeToRead) override {
        stream.read(reinterpret_cast<char*>(buffer), sizeToRead);
        return static_cast<int>(stream.gcount());
    }

    int getSize() override {
        return size;
    }

private:
    istream& stream;
    streampos start;
    int size;
};

static const char* requiredXmlAttribute(IrrXMLReader& xml, const char* name) {
    const char* value = xml.getAttributeValue(name);
    if (value == NULL)
        throw OpenMMException(string("BinarySerializer: missing XML attribute '")+name+"'");
    return value;
}

static double xmlDouble(IrrXMLReader& xml, const char* name) {
    return strtod2(requiredXmlAttribute(xml, name), NULL);
}

static long long xmlLong(IrrXMLReader& xml, const char* name, long long defaultValue) {
    const char* value = xml.getAttributeValue(name);
    if (value == NULL)
        return defaultValue;

    stringstream s(value);
    long long result;
    s >> result;
    if (!s)
        throw OpenMMException(string("BinarySerializer: invalid integer XML attribute '")+name+"'");
    return result;
}

static int parseXmlIntValue(const char* value, const char* name) {
    if (value == NULL || *value == '\0')
        throw OpenMMException(
            string("BinarySerializer: invalid integer XML attribute '")+name+"'");

    const char* p = value;
    bool negative = false;

    if (*p == '-') {
        negative = true;
        ++p;
    }
    else if (*p == '+') {
        ++p;
    }

    if (*p < '0' || *p > '9')
        throw OpenMMException(
            string("BinarySerializer: invalid integer XML attribute '")+name+"'");

    const uint64_t limit = negative
        ? static_cast<uint64_t>(numeric_limits<int>::max()) + 1
        : static_cast<uint64_t>(numeric_limits<int>::max());

    uint64_t result = 0;

    while (*p >= '0' && *p <= '9') {
        unsigned int digit = static_cast<unsigned int>(*p-'0');

        if (result > (limit-digit)/10)
            throw OpenMMException(
                string("BinarySerializer: integer XML attribute out of range '")+name+"'");

        result = result*10 + digit;
        ++p;
    }

    if (*p != '\0')
        throw OpenMMException(
            string("BinarySerializer: invalid integer XML attribute '")+name+"'");

    if (negative) {
        if (result == static_cast<uint64_t>(numeric_limits<int>::max())+1)
            return numeric_limits<int>::min();
        return -static_cast<int>(result);
    }

    return static_cast<int>(result);
}

static int xmlInt(IrrXMLReader& xml, const char* name) {
    return parseXmlIntValue(requiredXmlAttribute(xml, name), name);
}

static int xmlInt(IrrXMLReader& xml, const char* name, int defaultValue) {
    const char* value = xml.getAttributeValue(name);
    if (value == NULL)
        return defaultValue;

    return parseXmlIntValue(value, name);
}

static bool xmlBool(IrrXMLReader& xml, const char* name, bool defaultValue) {
    const char* value = xml.getAttributeValue(name);
    if (value == NULL)
        return defaultValue;

    if (strcmp(value, "1") == 0 || strcmp(value, "true") == 0)
        return true;
    if (strcmp(value, "0") == 0 || strcmp(value, "false") == 0)
        return false;

    throw OpenMMException(
        string("BinarySerializer: invalid boolean XML attribute '")+name+"'");
}

static string xmlString(IrrXMLReader& xml, const char* name,
                        const string& defaultValue) {
    const char* value = xml.getAttributeValue(name);
    return (value == NULL ? defaultValue : string(value));
}

// Lightweight, genuinely streaming reader for canonical OpenMM XML.
// Unlike irrXML, this never materializes the full XML document.  Seekable
// streams are read in 1 MiB chunks; unread bytes are put back when parsing
// completes so the caller's logical stream position is exact.  Non-seekable
// streams use exact one-byte extraction to avoid consuming following objects.
class FastXmlInput {
public:
    explicit FastXmlInput(istream& stream) : source(stream.rdbuf()), pos(0), end(0), seekable(false) {
        if (source == NULL)
            throw OpenMMException("BinarySerializer: XML stream has no stream buffer");
        seekable = canSeek(source);
        // Keep the v4-style 1 MiB buffered hot path for seekable streams.
        // Non-seekable streams use a one-byte buffer, which preserves exact
        // post-object position without adding a seekability branch to every
        // byte consumed by the normal State/System path.
        buffer.resize(seekable ? 1024*1024 : 1);
    }

    bool get(char& value) {
        if (pos == end && !refill())
            return false;
        value = buffer[pos++];
        return true;
    }

    void finish() {
        if (!seekable || pos == end)
            return;
        const streamoff unread = static_cast<streamoff>(end-pos);
        streampos restored;
        try {
            restored = source->pubseekoff(-unread, ios_base::cur, ios_base::in);
        }
        catch (const ios_base::failure&) {
            throw OpenMMException("BinarySerializer: failed to restore XML stream position");
        }
        if (restored == streampos(streamoff(-1)))
            throw OpenMMException("BinarySerializer: failed to restore XML stream position");
        pos = end = 0;
    }

private:
    streambuf* source;
    vector<char> buffer;
    size_t pos, end;
    bool seekable;

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

    bool refill() {
        streamsize got;
        try {
            got = source->sgetn(&buffer[0], static_cast<streamsize>(buffer.size()));
        }
        catch (const ios_base::failure&) {
            throw OpenMMException("BinarySerializer: error reading XML stream");
        }
        if (got <= 0) {
            pos = end = 0;
            return false;
        }
        pos = 0;
        end = static_cast<size_t>(got);
        return true;
    }
};

class FastXmlTagReader {
public:
    explicit FastXmlTagReader(istream& stream) : input(stream), firstByte(true) {
        tag.reserve(256);
    }

    const string& next() {
        char c;
        while (input.get(c)) {
            if (firstByte) {
                firstByte = false;
                if (static_cast<unsigned char>(c) == 0xEF) {
                    char b1, b2;
                    if (!input.get(b1) || !input.get(b2) ||
                            static_cast<unsigned char>(b1) != 0xBB ||
                            static_cast<unsigned char>(b2) != 0xBF)
                        throw OpenMMException("BinarySerializer: unsupported XML encoding");
                    continue;
                }
            }

            if (isSpace(c))
                continue;
            if (c != '<')
                throw OpenMMException("BinarySerializer: unexpected text in XML");

            char nextChar;
            if (!input.get(nextChar))
                throw OpenMMException("BinarySerializer: truncated XML markup");

            if (nextChar == '?') {
                skipSequence("?>");
                continue;
            }
            if (nextChar == '!') {
                char a, b;
                if (!input.get(a) || !input.get(b))
                    throw OpenMMException("BinarySerializer: truncated XML declaration");
                if (a == '-' && b == '-') {
                    skipSequence("-->");
                    continue;
                }
                throw OpenMMException("BinarySerializer: unsupported XML declaration");
            }

            tag.clear();
            tag.push_back('<');
            tag.push_back(nextChar);
            char quote = 0;
            while (input.get(c)) {
                tag.push_back(c);
                if (tag.size() > MAX_TAG_BYTES)
                    throw OpenMMException("BinarySerializer: XML tag is unreasonably large");
                if (quote != 0) {
                    if (c == quote)
                        quote = 0;
                }
                else if (c == '\'' || c == '"')
                    quote = c;
                else if (c == '>')
                    return tag;
            }
            throw OpenMMException("BinarySerializer: truncated XML tag");
        }
        throw OpenMMException("BinarySerializer: unexpected end of XML");
    }

    void finish() {
        input.finish();
    }

private:
    static const size_t MAX_TAG_BYTES = 1024*1024;
    FastXmlInput input;
    string tag;
    bool firstByte;

    static bool isSpace(char c) {
        return c == ' ' || c == '\t' || c == '\r' || c == '\n';
    }

    void skipSequence(const char* terminator) {
        const size_t n = strlen(terminator);
        size_t matched = 0;
        char c;
        while (input.get(c)) {
            if (c == terminator[matched]) {
                matched++;
                if (matched == n)
                    return;
            }
            else
                matched = (c == terminator[0] ? 1 : 0);
        }
        throw OpenMMException("BinarySerializer: truncated XML declaration");
    }
};

struct FastXmlTag {
    const string* text;
    size_t nameStart;
    size_t nameLength;
    size_t attributesStart;
    bool closing;
    bool empty;
};

struct FastXmlSpan {
    const char* data;
    size_t size;
};

static bool fastXmlSpace(char c) {
    return c == ' ' || c == '\t' || c == '\r' || c == '\n';
}

static bool fastXmlNameChar(char c) {
    return (c >= 'A' && c <= 'Z') || (c >= 'a' && c <= 'z') ||
           (c >= '0' && c <= '9') || c == '_' || c == ':' ||
           c == '-' || c == '.';
}

static FastXmlTag parseFastXmlTag(const string& text) {
    if (text.size() < 3 || text[0] != '<' || text[text.size()-1] != '>')
        throw OpenMMException("BinarySerializer: malformed XML tag");

    size_t pos = 1;
    bool closing = false;
    if (text[pos] == '/') {
        closing = true;
        ++pos;
    }
    const size_t nameStart = pos;
    while (pos < text.size()-1 && fastXmlNameChar(text[pos]))
        ++pos;
    if (pos == nameStart)
        throw OpenMMException("BinarySerializer: malformed XML element name");

    size_t tail = text.size()-1;
    while (tail > pos && fastXmlSpace(text[tail-1]))
        --tail;
    bool empty = false;
    if (!closing && tail > pos && text[tail-1] == '/')
        empty = true;
    if (closing) {
        for (size_t i = pos; i < tail; ++i) {
            if (!fastXmlSpace(text[i]))
                throw OpenMMException("BinarySerializer: malformed XML closing tag");
        }
    }

    FastXmlTag tag;
    tag.text = &text;
    tag.nameStart = nameStart;
    tag.nameLength = pos-nameStart;
    tag.attributesStart = pos;
    tag.closing = closing;
    tag.empty = empty;
    return tag;
}

static bool fastXmlNameEquals(const FastXmlTag& tag, const char* expected) {
    const size_t length = strlen(expected);
    return tag.nameLength == length &&
           memcmp(tag.text->data()+tag.nameStart, expected, length) == 0;
}

class FastXmlAttributes {
public:
    explicit FastXmlAttributes(const FastXmlTag& tag) : text(*tag.text), pos(tag.attributesStart) {}

    bool next(FastXmlSpan& name, FastXmlSpan& value) {
        const size_t limit = text.size()-1;
        while (pos < limit && fastXmlSpace(text[pos]))
            ++pos;
        if (pos >= limit || text[pos] == '/' || text[pos] == '>')
            return false;

        const size_t nameStart = pos;
        while (pos < limit && fastXmlNameChar(text[pos]))
            ++pos;
        if (pos == nameStart)
            throw OpenMMException("BinarySerializer: malformed XML attribute name");
        name.data = text.data()+nameStart;
        name.size = pos-nameStart;

        while (pos < limit && fastXmlSpace(text[pos]))
            ++pos;
        if (pos >= limit || text[pos] != '=')
            throw OpenMMException("BinarySerializer: malformed XML attribute");
        ++pos;
        while (pos < limit && fastXmlSpace(text[pos]))
            ++pos;
        if (pos >= limit || (text[pos] != '"' && text[pos] != '\''))
            throw OpenMMException("BinarySerializer: malformed XML attribute value");
        const char quote = text[pos++];
        const size_t valueStart = pos;
        while (pos < limit && text[pos] != quote)
            ++pos;
        if (pos >= limit)
            throw OpenMMException("BinarySerializer: unterminated XML attribute value");
        value.data = text.data()+valueStart;
        value.size = pos-valueStart;
        ++pos;
        return true;
    }

private:
    const string& text;
    size_t pos;
};

static bool fastXmlSpanEquals(const FastXmlSpan& span, const char* value) {
    const size_t n = strlen(value);
    return span.size == n && memcmp(span.data, value, n) == 0;
}

static double fastXmlDoubleSpan(const FastXmlSpan& span, const char* attributeName) {
    char* end = NULL;
    double value = strtod2(span.data, &end);
    if (end != span.data+span.size)
        throw OpenMMException(string("BinarySerializer: invalid floating point XML attribute '")+attributeName+"'");
    return value;
}

static long long fastXmlLongSpan(const FastXmlSpan& span, const char* attributeName) {
    if (span.size == 0)
        throw OpenMMException(string("BinarySerializer: invalid integer XML attribute '")+attributeName+"'");
    size_t pos = 0;
    bool negative = false;
    if (span.data[pos] == '-' || span.data[pos] == '+') {
        negative = span.data[pos] == '-';
        if (++pos == span.size)
            throw OpenMMException(string("BinarySerializer: invalid integer XML attribute '")+attributeName+"'");
    }
    uint64_t value = 0;
    const uint64_t positiveLimit = static_cast<uint64_t>(numeric_limits<long long>::max());
    const uint64_t limit = negative ? positiveLimit+UINT64_C(1) : positiveLimit;
    for (; pos < span.size; ++pos) {
        const char c = span.data[pos];
        if (c < '0' || c > '9')
            throw OpenMMException(string("BinarySerializer: invalid integer XML attribute '")+attributeName+"'");
        const unsigned int digit = static_cast<unsigned int>(c-'0');
        if (value > (limit-digit)/10)
            throw OpenMMException(string("BinarySerializer: integer XML attribute out of range '")+attributeName+"'");
        value = value*10 + digit;
    }
    if (!negative)
        return static_cast<long long>(value);
    if (value == positiveLimit+UINT64_C(1))
        return numeric_limits<long long>::min();
    return -static_cast<long long>(value);
}

static string fastXmlPlainName(const FastXmlSpan& span) {
    // The canonical OpenMM State writer uses parameter names as XML attribute
    // names.  If escaping or unusual syntax is present, let the stock parser
    // handle it rather than implementing the full XML Names grammar here.
    for (size_t i = 0; i < span.size; ++i) {
        if (!fastXmlNameChar(span.data[i]))
            throw OpenMMException("BinarySerializer: unsupported escaped State parameter name");
    }
    return string(span.data, span.size);
}

static bool findFastXmlAttribute(const FastXmlTag& tag, const char* expected,
                                 FastXmlSpan& value) {
    FastXmlAttributes attrs(tag);
    FastXmlSpan name, candidate;
    while (attrs.next(name, candidate)) {
        if (fastXmlSpanEquals(name, expected)) {
            value = candidate;
            return true;
        }
    }
    return false;
}

static FastXmlSpan requireFastXmlAttribute(const FastXmlTag& tag, const char* name) {
    FastXmlSpan value;
    if (!findFastXmlAttribute(tag, name, value))
        throw OpenMMException(string("BinarySerializer: missing XML attribute '")+name+"'");
    return value;
}

static string fastXmlPlainString(const FastXmlSpan& span, const char* attributeName) {
    // The fast path intentionally handles only the canonical strings produced by
    // OpenMM that need no entity decoding.  Escaped strings fall back to the
    // stock XML parser rather than risking a semantic mismatch.
    for (size_t i = 0; i < span.size; ++i) {
        if (span.data[i] == '&')
            throw OpenMMException(string("BinarySerializer: escaped XML string requires stock parser for attribute '")+attributeName+"'");
    }
    return string(span.data, span.size);
}

static string fastXmlString(const FastXmlTag& tag, const char* name) {
    return fastXmlPlainString(requireFastXmlAttribute(tag, name), name);
}

static string fastXmlString(const FastXmlTag& tag, const char* name,
                            const string& defaultValue) {
    FastXmlSpan value;
    if (!findFastXmlAttribute(tag, name, value))
        return defaultValue;
    return fastXmlPlainString(value, name);
}

static double fastXmlDouble(const FastXmlTag& tag, const char* name) {
    return fastXmlDoubleSpan(requireFastXmlAttribute(tag, name), name);
}

static double fastXmlDouble(const FastXmlTag& tag, const char* name, double defaultValue) {
    FastXmlSpan value;
    if (!findFastXmlAttribute(tag, name, value))
        return defaultValue;
    return fastXmlDoubleSpan(value, name);
}

static int fastXmlInt(const FastXmlTag& tag, const char* name) {
    const long long value = fastXmlLongSpan(requireFastXmlAttribute(tag, name), name);
    if (value < numeric_limits<int>::min() || value > numeric_limits<int>::max())
        throw OpenMMException(string("BinarySerializer: integer XML attribute out of range '")+name+"'");
    return static_cast<int>(value);
}

static int fastXmlInt(const FastXmlTag& tag, const char* name, int defaultValue) {
    FastXmlSpan value;
    if (!findFastXmlAttribute(tag, name, value))
        return defaultValue;
    const long long parsed = fastXmlLongSpan(value, name);
    if (parsed < numeric_limits<int>::min() || parsed > numeric_limits<int>::max())
        throw OpenMMException(string("BinarySerializer: integer XML attribute out of range '")+name+"'");
    return static_cast<int>(parsed);
}

static bool fastXmlBool(const FastXmlTag& tag, const char* name, bool defaultValue) {
    FastXmlSpan value;
    if (!findFastXmlAttribute(tag, name, value))
        return defaultValue;
    if (fastXmlSpanEquals(value, "1") || fastXmlSpanEquals(value, "true"))
        return true;
    if (fastXmlSpanEquals(value, "0") || fastXmlSpanEquals(value, "false"))
        return false;
    throw OpenMMException(string("BinarySerializer: invalid boolean XML attribute '")+name+"'");
}

static void consumeLeafClose(FastXmlTagReader& reader, const FastXmlTag& opening,
                             const char* elementName) {
    if (opening.empty)
        return;
    FastXmlTag close = parseFastXmlTag(reader.next());
    if (!close.closing || !fastXmlNameEquals(close, elementName))
        throw OpenMMException(string("BinarySerializer: unexpected child in XML element '")+elementName+"'");
}

static Vec3 parseStateVec3Tag(const FastXmlTag& tag, const char* expectedName) {
    if (tag.closing || !tag.empty || !fastXmlNameEquals(tag, expectedName))
        throw OpenMMException(string("BinarySerializer: unexpected element in State XML array; expected ")+expectedName);

    bool haveX = false, haveY = false, haveZ = false;
    double x = 0, y = 0, z = 0;
    FastXmlAttributes attrs(tag);
    FastXmlSpan name, value;
    while (attrs.next(name, value)) {
        if (fastXmlSpanEquals(name, "x")) {
            x = fastXmlDoubleSpan(value, "x");
            haveX = true;
        }
        else if (fastXmlSpanEquals(name, "y")) {
            y = fastXmlDoubleSpan(value, "y");
            haveY = true;
        }
        else if (fastXmlSpanEquals(name, "z")) {
            z = fastXmlDoubleSpan(value, "z");
            haveZ = true;
        }
    }
    if (!haveX || !haveY || !haveZ)
        throw OpenMMException("BinarySerializer: incomplete Vec3 in State XML");
    return Vec3(x, y, z);
}

static void parseStateBox(FastXmlTagReader& reader, const FastXmlTag& opening,
                          Vec3& a, Vec3& b, Vec3& c) {
    if (opening.closing || opening.empty || !fastXmlNameEquals(opening, "PeriodicBoxVectors"))
        throw OpenMMException("BinarySerializer: malformed PeriodicBoxVectors section");

    FastXmlTag tagA = parseFastXmlTag(reader.next());
    a = parseStateVec3Tag(tagA, "A");
    FastXmlTag tagB = parseFastXmlTag(reader.next());
    b = parseStateVec3Tag(tagB, "B");
    FastXmlTag tagC = parseFastXmlTag(reader.next());
    c = parseStateVec3Tag(tagC, "C");
    FastXmlTag close = parseFastXmlTag(reader.next());
    if (!close.closing || !fastXmlNameEquals(close, "PeriodicBoxVectors"))
        throw OpenMMException("BinarySerializer: malformed PeriodicBoxVectors closing tag");
}

static map<string, double> parseStateParameters(const FastXmlTag& tag) {
    if (tag.closing || !tag.empty || !fastXmlNameEquals(tag, "Parameters"))
        throw OpenMMException("BinarySerializer: malformed Parameters section");
    map<string, double> result;
    FastXmlAttributes attrs(tag);
    FastXmlSpan name, value;
    while (attrs.next(name, value)) {
        string parameter = fastXmlPlainName(name);
        result[parameter] = fastXmlDoubleSpan(value, parameter.c_str());
    }
    return result;
}

static void parseStateEnergies(const FastXmlTag& tag, double& kinetic, double& potential) {
    if (tag.closing || !tag.empty || !fastXmlNameEquals(tag, "Energies"))
        throw OpenMMException("BinarySerializer: malformed Energies section");
    bool haveKinetic = false, havePotential = false;
    FastXmlAttributes attrs(tag);
    FastXmlSpan name, value;
    while (attrs.next(name, value)) {
        if (fastXmlSpanEquals(name, "KineticEnergy")) {
            kinetic = fastXmlDoubleSpan(value, "KineticEnergy");
            haveKinetic = true;
        }
        else if (fastXmlSpanEquals(name, "PotentialEnergy")) {
            potential = fastXmlDoubleSpan(value, "PotentialEnergy");
            havePotential = true;
        }
    }
    if (!haveKinetic || !havePotential)
        throw OpenMMException("BinarySerializer: incomplete Energies section");
}

static void parseStateVec3Section(FastXmlTagReader& reader, const FastXmlTag& opening,
                                  const char* sectionName, const char* elementName,
                                  vector<Vec3>& values, size_t reserveHint) {
    if (opening.closing || !fastXmlNameEquals(opening, sectionName))
        throw OpenMMException(string("BinarySerializer: malformed ")+sectionName+" section");
    if (reserveHint != 0)
        values.reserve(reserveHint);
    if (opening.empty)
        return;

    while (true) {
        FastXmlTag tag = parseFastXmlTag(reader.next());
        if (tag.closing) {
            if (!fastXmlNameEquals(tag, sectionName))
                throw OpenMMException(string("BinarySerializer: unexpected closing tag in ")+sectionName);
            return;
        }
        if (values.size() >= static_cast<size_t>(MAX_CONTAINER_ELEMENTS))
            throw OpenMMException("BinarySerializer: State XML array is too large");
        values.push_back(parseStateVec3Tag(tag, elementName));
    }
}

static State* parseStateXmlStreaming(istream& stream) {
    FastXmlTagReader reader(stream);
    FastXmlTag root = parseFastXmlTag(reader.next());
    if (root.closing || root.empty || !fastXmlNameEquals(root, "State"))
        throw OpenMMException("BinarySerializer: XML object is not a State");

    bool haveType = false, haveVersion = false, haveTime = false;
    string type;
    long long stepCount = 0;
    double time = 0;
    FastXmlAttributes rootAttrs(root);
    FastXmlSpan attrName, attrValue;
    while (rootAttrs.next(attrName, attrValue)) {
        if (fastXmlSpanEquals(attrName, "type")) {
            type.assign(attrValue.data, attrValue.size);
            haveType = true;
        }
        else if (fastXmlSpanEquals(attrName, "version")) {
            if (!fastXmlSpanEquals(attrValue, "1"))
                throw OpenMMException("BinarySerializer: unsupported State XML version");
            haveVersion = true;
        }
        else if (fastXmlSpanEquals(attrName, "time")) {
            time = fastXmlDoubleSpan(attrValue, "time");
            haveTime = true;
        }
        else if (fastXmlSpanEquals(attrName, "stepCount"))
            stepCount = fastXmlLongSpan(attrValue, "stepCount");
    }
    if (!haveType || type != "State")
        throw OpenMMException("BinarySerializer: XML type is not State");
    if (!haveVersion || !haveTime)
        throw OpenMMException("BinarySerializer: incomplete State root attributes");

    Vec3 a, b, c;
    bool haveBox = false;
    int types = 0;
    int stage = 0;
    map<string, double> parameters;
    double kinetic = 0, potential = 0;
    vector<Vec3> positions, velocities, forces;

    while (true) {
        FastXmlTag tag = parseFastXmlTag(reader.next());
        if (tag.closing) {
            if (!fastXmlNameEquals(tag, "State"))
                throw OpenMMException("BinarySerializer: unexpected closing tag in State XML");
            break;
        }

        if (fastXmlNameEquals(tag, "PeriodicBoxVectors")) {
            if (stage != 0)
                throw OpenMMException("BinarySerializer: unexpected State XML ordering");
            parseStateBox(reader, tag, a, b, c);
            haveBox = true;
            stage = 1;
        }
        else if (fastXmlNameEquals(tag, "Parameters")) {
            if (!haveBox || stage >= 2)
                throw OpenMMException("BinarySerializer: unexpected Parameters section");
            parameters = parseStateParameters(tag);
            types |= State::Parameters;
            stage = 2;
        }
        else if (fastXmlNameEquals(tag, "Energies")) {
            if (!haveBox || stage >= 3)
                throw OpenMMException("BinarySerializer: unexpected Energies section");
            parseStateEnergies(tag, kinetic, potential);
            types |= State::Energy;
            stage = 3;
        }
        else if (fastXmlNameEquals(tag, "Positions")) {
            if (!haveBox || stage >= 4)
                throw OpenMMException("BinarySerializer: unexpected Positions section");
            parseStateVec3Section(reader, tag, "Positions", "Position", positions, 0);
            types |= State::Positions;
            stage = 4;
        }
        else if (fastXmlNameEquals(tag, "Velocities")) {
            if (!haveBox || stage >= 5)
                throw OpenMMException("BinarySerializer: unexpected Velocities section");
            parseStateVec3Section(reader, tag, "Velocities", "Velocity", velocities,
                                  (types & State::Positions) ? positions.size() : 0);
            if ((types & State::Positions) && velocities.size() != positions.size())
                throw OpenMMException("BinarySerializer: State particle array size mismatch");
            types |= State::Velocities;
            stage = 5;
        }
        else if (fastXmlNameEquals(tag, "Forces")) {
            if (!haveBox || stage >= 6)
                throw OpenMMException("BinarySerializer: unexpected Forces section");
            const size_t hint = (types & State::Positions) ? positions.size() :
                                ((types & State::Velocities) ? velocities.size() : 0);
            parseStateVec3Section(reader, tag, "Forces", "Force", forces, hint);
            if (((types & State::Positions) && forces.size() != positions.size()) ||
                    (!(types & State::Positions) && (types & State::Velocities) &&
                     forces.size() != velocities.size()))
                throw OpenMMException("BinarySerializer: State particle array size mismatch");
            types |= State::Forces;
            stage = 6;
        }
        else if (fastXmlNameEquals(tag, "IntegratorParameters"))
            throw OpenMMException("BinarySerializer: IntegratorParameters are not supported by the fast State path");
        else
            throw OpenMMException("BinarySerializer: unsupported State XML element");
    }

    if (!haveBox)
        throw OpenMMException("BinarySerializer: State has no PeriodicBoxVectors");

    // If Positions are absent, still enforce consistency between any other two
    // particle arrays that are present.
    if (!(types & State::Positions) && (types & State::Velocities) &&
            (types & State::Forces) && velocities.size() != forces.size())
        throw OpenMMException("BinarySerializer: State particle array size mismatch");

    // Restore any read-ahead so a successful XML read leaves a seekable stream
    // positioned immediately after </State>.
    reader.finish();

    State::StateBuilder builder(time, stepCount);
    builder.setPeriodicBoxVectors(a, b, c);
    if (types & State::Parameters)
        builder.setParameters(std::move(parameters));
    if (types & State::Energy)
        builder.setEnergy(kinetic, potential);
    if (types & State::Positions)
        builder.setPositions(std::move(positions));
    if (types & State::Velocities)
        builder.setVelocities(std::move(velocities));
    if (types & State::Forces)
        builder.setForces(std::move(forces));
    return new State(builder.takeState());
}

// Streaming System XML reader for the canonical OpenMM schema used by Core 28.
// It builds System and the supported Force objects directly while keeping only
// one tag and a bounded input buffer in memory.  Unsupported schema features
// throw and are handled by XmlSerializer's stock-parser fallback.
static void parseSystemBoxStreaming(FastXmlTagReader& reader, const FastXmlTag& opening,
                                    System& system) {
    if (opening.closing || opening.empty || !fastXmlNameEquals(opening, "PeriodicBoxVectors"))
        throw OpenMMException("BinarySerializer: malformed System PeriodicBoxVectors section");

    FastXmlTag aTag = parseFastXmlTag(reader.next());
    Vec3 a = parseStateVec3Tag(aTag, "A");
    FastXmlTag bTag = parseFastXmlTag(reader.next());
    Vec3 b = parseStateVec3Tag(bTag, "B");
    FastXmlTag cTag = parseFastXmlTag(reader.next());
    Vec3 c = parseStateVec3Tag(cTag, "C");
    FastXmlTag close = parseFastXmlTag(reader.next());
    if (!close.closing || !fastXmlNameEquals(close, "PeriodicBoxVectors"))
        throw OpenMMException("BinarySerializer: malformed System PeriodicBoxVectors closing tag");
    system.setDefaultPeriodicBoxVectors(a, b, c);
}

static void parseSystemVirtualSiteStreaming(FastXmlTagReader& reader, const FastXmlTag& tag,
                                            System& system, int particle) {
    if (tag.closing)
        throw OpenMMException("BinarySerializer: malformed virtual site element");

    if (fastXmlNameEquals(tag, "TwoParticleAverageSite")) {
        system.setVirtualSite(particle, new TwoParticleAverageSite(
            fastXmlInt(tag, "p1"), fastXmlInt(tag, "p2"),
            fastXmlDouble(tag, "w1"), fastXmlDouble(tag, "w2")));
        consumeLeafClose(reader, tag, "TwoParticleAverageSite");
    }
    else if (fastXmlNameEquals(tag, "ThreeParticleAverageSite")) {
        system.setVirtualSite(particle, new ThreeParticleAverageSite(
            fastXmlInt(tag, "p1"), fastXmlInt(tag, "p2"), fastXmlInt(tag, "p3"),
            fastXmlDouble(tag, "w1"), fastXmlDouble(tag, "w2"), fastXmlDouble(tag, "w3")));
        consumeLeafClose(reader, tag, "ThreeParticleAverageSite");
    }
    else if (fastXmlNameEquals(tag, "OutOfPlaneSite")) {
        system.setVirtualSite(particle, new OutOfPlaneSite(
            fastXmlInt(tag, "p1"), fastXmlInt(tag, "p2"), fastXmlInt(tag, "p3"),
            fastXmlDouble(tag, "w12"), fastXmlDouble(tag, "w13"), fastXmlDouble(tag, "wc")));
        consumeLeafClose(reader, tag, "OutOfPlaneSite");
    }
    else if (fastXmlNameEquals(tag, "LocalCoordinatesSite")) {
        vector<int> particles;
        vector<double> wo, wx, wy;
        for (int j = 1; ; ++j) {
            stringstream ss;
            ss << j;
            const string index = ss.str();
            const string pName = "p"+index;
            FastXmlSpan pValue;
            if (!findFastXmlAttribute(tag, pName.c_str(), pValue))
                break;
            const long long p = fastXmlLongSpan(pValue, pName.c_str());
            if (p < numeric_limits<int>::min() || p > numeric_limits<int>::max())
                throw OpenMMException("BinarySerializer: virtual-site particle index out of range");
            particles.push_back(static_cast<int>(p));
            wo.push_back(fastXmlDouble(tag, ("wo"+index).c_str()));
            wx.push_back(fastXmlDouble(tag, ("wx"+index).c_str()));
            wy.push_back(fastXmlDouble(tag, ("wy"+index).c_str()));
        }
        if (particles.empty())
            throw OpenMMException("BinarySerializer: LocalCoordinatesSite has no particles");
        Vec3 position(fastXmlDouble(tag, "pos1"), fastXmlDouble(tag, "pos2"),
                      fastXmlDouble(tag, "pos3"));
        system.setVirtualSite(particle, new LocalCoordinatesSite(particles, wo, wx, wy, position));
        consumeLeafClose(reader, tag, "LocalCoordinatesSite");
    }
    else {
        throw OpenMMException("BinarySerializer: unsupported virtual site type in fast System path");
    }
}

static void parseSystemParticlesStreaming(FastXmlTagReader& reader, const FastXmlTag& opening,
                                          System& system) {
    if (opening.closing || !fastXmlNameEquals(opening, "Particles"))
        throw OpenMMException("BinarySerializer: malformed System Particles section");
    if (opening.empty)
        return;

    while (true) {
        FastXmlTag tag = parseFastXmlTag(reader.next());
        if (tag.closing) {
            if (!fastXmlNameEquals(tag, "Particles"))
                throw OpenMMException("BinarySerializer: unexpected closing tag in System Particles");
            return;
        }
        if (!fastXmlNameEquals(tag, "Particle"))
            throw OpenMMException("BinarySerializer: unexpected element in System Particles");
        if (system.getNumParticles() >= numeric_limits<int>::max())
            throw OpenMMException("BinarySerializer: too many System particles");

        const int particle = system.getNumParticles();
        system.addParticle(fastXmlDouble(tag, "mass"));
        if (tag.empty)
            continue;

        FastXmlTag child = parseFastXmlTag(reader.next());
        if (child.closing) {
            if (!fastXmlNameEquals(child, "Particle"))
                throw OpenMMException("BinarySerializer: malformed Particle closing tag");
            continue;
        }
        parseSystemVirtualSiteStreaming(reader, child, system, particle);
        FastXmlTag close = parseFastXmlTag(reader.next());
        if (!close.closing || !fastXmlNameEquals(close, "Particle"))
            throw OpenMMException("BinarySerializer: Particle has multiple or malformed virtual sites");
    }
}

static void parseSystemConstraintsStreaming(FastXmlTagReader& reader, const FastXmlTag& opening,
                                            System& system) {
    if (opening.closing || !fastXmlNameEquals(opening, "Constraints"))
        throw OpenMMException("BinarySerializer: malformed System Constraints section");
    if (opening.empty)
        return;

    while (true) {
        FastXmlTag tag = parseFastXmlTag(reader.next());
        if (tag.closing) {
            if (!fastXmlNameEquals(tag, "Constraints"))
                throw OpenMMException("BinarySerializer: unexpected closing tag in Constraints");
            return;
        }
        if (!fastXmlNameEquals(tag, "Constraint"))
            throw OpenMMException("BinarySerializer: unexpected element in System Constraints");
        system.addConstraint(fastXmlInt(tag, "p1"), fastXmlInt(tag, "p2"),
                             fastXmlDouble(tag, "d"));
        consumeLeafClose(reader, tag, "Constraint");
    }
}

enum FastNonbondedSectionKind {
    FAST_NB_GLOBAL_PARAMETERS,
    FAST_NB_PARTICLE_OFFSETS,
    FAST_NB_EXCEPTION_OFFSETS,
    FAST_NB_PARTICLES,
    FAST_NB_EXCEPTIONS
};

static void parseNonbondedSectionStreaming(FastXmlTagReader& reader, const FastXmlTag& opening,
                                           FastNonbondedSectionKind kind,
                                           NonbondedForce& force) {
    const char* sectionName = NULL;
    const char* elementName = NULL;
    switch (kind) {
        case FAST_NB_GLOBAL_PARAMETERS: sectionName = "GlobalParameters"; elementName = "Parameter"; break;
        case FAST_NB_PARTICLE_OFFSETS: sectionName = "ParticleOffsets"; elementName = "Offset"; break;
        case FAST_NB_EXCEPTION_OFFSETS: sectionName = "ExceptionOffsets"; elementName = "Offset"; break;
        case FAST_NB_PARTICLES: sectionName = "Particles"; elementName = "Particle"; break;
        case FAST_NB_EXCEPTIONS: sectionName = "Exceptions"; elementName = "Exception"; break;
    }
    if (opening.closing || !fastXmlNameEquals(opening, sectionName))
        throw OpenMMException("BinarySerializer: malformed NonbondedForce section");
    if (opening.empty)
        return;

    while (true) {
        FastXmlTag tag = parseFastXmlTag(reader.next());
        if (tag.closing) {
            if (!fastXmlNameEquals(tag, sectionName))
                throw OpenMMException("BinarySerializer: unexpected closing tag in NonbondedForce section");
            return;
        }
        if (!fastXmlNameEquals(tag, elementName))
            throw OpenMMException("BinarySerializer: unexpected element in NonbondedForce section");

        switch (kind) {
            case FAST_NB_GLOBAL_PARAMETERS:
                force.addGlobalParameter(fastXmlString(tag, "name"), fastXmlDouble(tag, "default"));
                break;
            case FAST_NB_PARTICLE_OFFSETS:
                force.addParticleParameterOffset(fastXmlString(tag, "parameter"),
                    fastXmlInt(tag, "particle"), fastXmlDouble(tag, "q"),
                    fastXmlDouble(tag, "sig"), fastXmlDouble(tag, "eps"));
                break;
            case FAST_NB_EXCEPTION_OFFSETS:
                force.addExceptionParameterOffset(fastXmlString(tag, "parameter"),
                    fastXmlInt(tag, "exception"), fastXmlDouble(tag, "q"),
                    fastXmlDouble(tag, "sig"), fastXmlDouble(tag, "eps"));
                break;
            case FAST_NB_PARTICLES:
                force.addParticle(fastXmlDouble(tag, "q"), fastXmlDouble(tag, "sig"),
                                  fastXmlDouble(tag, "eps"));
                break;
            case FAST_NB_EXCEPTIONS:
                force.addException(fastXmlInt(tag, "p1"), fastXmlInt(tag, "p2"),
                    fastXmlDouble(tag, "q"), fastXmlDouble(tag, "sig"),
                    fastXmlDouble(tag, "eps"));
                break;
        }
        consumeLeafClose(reader, tag, elementName);
    }
}

static Force* parseNonbondedForceStreaming(FastXmlTagReader& reader, const FastXmlTag& opening) {
    const int version = fastXmlInt(opening, "version");
    if (version < 1 || version > 4)
        throw OpenMMException("BinarySerializer: unsupported NonbondedForce XML version");

    unique_ptr<NonbondedForce> force(new NonbondedForce());
    force->setForceGroup(fastXmlInt(opening, "forceGroup", 0));
    force->setName(fastXmlString(opening, "name", force->getName()));
    force->setNonbondedMethod(static_cast<NonbondedForce::NonbondedMethod>(fastXmlInt(opening, "method")));
    force->setCutoffDistance(fastXmlDouble(opening, "cutoff"));
    force->setUseSwitchingFunction(fastXmlBool(opening, "useSwitchingFunction", false));
    force->setSwitchingDistance(fastXmlDouble(opening, "switchingDistance", -1.0));
    force->setEwaldErrorTolerance(fastXmlDouble(opening, "ewaldTolerance"));
    force->setReactionFieldDielectric(fastXmlDouble(opening, "rfDielectric"));
    force->setUseDispersionCorrection(fastXmlInt(opening, "dispersionCorrection") != 0);
    force->setIncludeDirectSpace(fastXmlBool(opening, "includeDirectSpace", force->getIncludeDirectSpace()));
    force->setPMEParameters(fastXmlDouble(opening, "alpha", 0.0),
                            fastXmlInt(opening, "nx", 0), fastXmlInt(opening, "ny", 0),
                            fastXmlInt(opening, "nz", 0));
    if (version >= 2)
        force->setLJPMEParameters(fastXmlDouble(opening, "ljAlpha", 0.0),
                                  fastXmlInt(opening, "ljnx", 0), fastXmlInt(opening, "ljny", 0),
                                  fastXmlInt(opening, "ljnz", 0));
    force->setReciprocalSpaceForceGroup(fastXmlInt(opening, "recipForceGroup", -1));
    if (version >= 4)
        force->setExceptionsUsePeriodicBoundaryConditions(fastXmlInt(opening, "exceptionsUsePeriodic") != 0);

    bool haveGlobal = false, haveParticleOffsets = false, haveExceptionOffsets = false;
    bool haveParticles = false, haveExceptions = false;
    if (opening.empty)
        throw OpenMMException("BinarySerializer: incomplete NonbondedForce");

    while (true) {
        FastXmlTag tag = parseFastXmlTag(reader.next());
        if (tag.closing) {
            if (!fastXmlNameEquals(tag, "Force"))
                throw OpenMMException("BinarySerializer: unexpected closing tag in NonbondedForce");
            if (!haveParticles || !haveExceptions)
                throw OpenMMException("BinarySerializer: incomplete NonbondedForce");
            if (version >= 3 && (!haveGlobal || !haveParticleOffsets || !haveExceptionOffsets))
                throw OpenMMException("BinarySerializer: incomplete NonbondedForce parameter sections");
            return force.release();
        }
        if (fastXmlNameEquals(tag, "GlobalParameters")) {
            if (version < 3 || haveGlobal)
                throw OpenMMException("BinarySerializer: unexpected GlobalParameters");
            parseNonbondedSectionStreaming(reader, tag, FAST_NB_GLOBAL_PARAMETERS, *force);
            haveGlobal = true;
        }
        else if (fastXmlNameEquals(tag, "ParticleOffsets")) {
            if (version < 3 || haveParticleOffsets)
                throw OpenMMException("BinarySerializer: unexpected ParticleOffsets");
            parseNonbondedSectionStreaming(reader, tag, FAST_NB_PARTICLE_OFFSETS, *force);
            haveParticleOffsets = true;
        }
        else if (fastXmlNameEquals(tag, "ExceptionOffsets")) {
            if (version < 3 || haveExceptionOffsets)
                throw OpenMMException("BinarySerializer: unexpected ExceptionOffsets");
            parseNonbondedSectionStreaming(reader, tag, FAST_NB_EXCEPTION_OFFSETS, *force);
            haveExceptionOffsets = true;
        }
        else if (fastXmlNameEquals(tag, "Particles")) {
            if (haveParticles)
                throw OpenMMException("BinarySerializer: duplicate NonbondedForce Particles");
            parseNonbondedSectionStreaming(reader, tag, FAST_NB_PARTICLES, *force);
            haveParticles = true;
        }
        else if (fastXmlNameEquals(tag, "Exceptions")) {
            if (haveExceptions)
                throw OpenMMException("BinarySerializer: duplicate NonbondedForce Exceptions");
            parseNonbondedSectionStreaming(reader, tag, FAST_NB_EXCEPTIONS, *force);
            haveExceptions = true;
        }
        else
            throw OpenMMException("BinarySerializer: unsupported NonbondedForce child in fast System path");
    }
}

static Force* parseHarmonicBondForceStreaming(FastXmlTagReader& reader, const FastXmlTag& opening) {
    const int version = fastXmlInt(opening, "version");
    if (version < 1 || version > 2)
        throw OpenMMException("BinarySerializer: unsupported HarmonicBondForce XML version");
    unique_ptr<HarmonicBondForce> force(new HarmonicBondForce());
    force->setForceGroup(fastXmlInt(opening, "forceGroup", 0));
    force->setName(fastXmlString(opening, "name", force->getName()));
    if (version > 1)
        force->setUsesPeriodicBoundaryConditions(fastXmlBool(opening, "usesPeriodic", false));
    bool haveBonds = false;
    if (opening.empty)
        throw OpenMMException("BinarySerializer: HarmonicBondForce has no Bonds section");
    while (true) {
        FastXmlTag tag = parseFastXmlTag(reader.next());
        if (tag.closing) {
            if (!fastXmlNameEquals(tag, "Force"))
                throw OpenMMException("BinarySerializer: unexpected closing tag in HarmonicBondForce");
            if (!haveBonds)
                throw OpenMMException("BinarySerializer: HarmonicBondForce has no Bonds section");
            return force.release();
        }
        if (!fastXmlNameEquals(tag, "Bonds") || haveBonds)
            throw OpenMMException("BinarySerializer: unexpected HarmonicBondForce element");
        haveBonds = true;
        if (tag.empty)
            continue;
        while (true) {
            FastXmlTag bond = parseFastXmlTag(reader.next());
            if (bond.closing) {
                if (!fastXmlNameEquals(bond, "Bonds"))
                    throw OpenMMException("BinarySerializer: unexpected closing tag in Bonds");
                break;
            }
            if (!fastXmlNameEquals(bond, "Bond"))
                throw OpenMMException("BinarySerializer: unexpected element in HarmonicBondForce Bonds");
            force->addBond(fastXmlInt(bond, "p1"), fastXmlInt(bond, "p2"),
                           fastXmlDouble(bond, "d"), fastXmlDouble(bond, "k"));
            consumeLeafClose(reader, bond, "Bond");
        }
    }
}

static Force* parseHarmonicAngleForceStreaming(FastXmlTagReader& reader, const FastXmlTag& opening) {
    const int version = fastXmlInt(opening, "version");
    if (version < 1 || version > 2)
        throw OpenMMException("BinarySerializer: unsupported HarmonicAngleForce XML version");
    unique_ptr<HarmonicAngleForce> force(new HarmonicAngleForce());
    force->setForceGroup(fastXmlInt(opening, "forceGroup", 0));
    force->setName(fastXmlString(opening, "name", force->getName()));
    if (version > 1)
        force->setUsesPeriodicBoundaryConditions(fastXmlBool(opening, "usesPeriodic", false));
    bool haveAngles = false;
    if (opening.empty)
        throw OpenMMException("BinarySerializer: HarmonicAngleForce has no Angles section");
    while (true) {
        FastXmlTag tag = parseFastXmlTag(reader.next());
        if (tag.closing) {
            if (!fastXmlNameEquals(tag, "Force"))
                throw OpenMMException("BinarySerializer: unexpected closing tag in HarmonicAngleForce");
            if (!haveAngles)
                throw OpenMMException("BinarySerializer: HarmonicAngleForce has no Angles section");
            return force.release();
        }
        if (!fastXmlNameEquals(tag, "Angles") || haveAngles)
            throw OpenMMException("BinarySerializer: unexpected HarmonicAngleForce element");
        haveAngles = true;
        if (tag.empty)
            continue;
        while (true) {
            FastXmlTag angle = parseFastXmlTag(reader.next());
            if (angle.closing) {
                if (!fastXmlNameEquals(angle, "Angles"))
                    throw OpenMMException("BinarySerializer: unexpected closing tag in Angles");
                break;
            }
            if (!fastXmlNameEquals(angle, "Angle"))
                throw OpenMMException("BinarySerializer: unexpected element in HarmonicAngleForce Angles");
            force->addAngle(fastXmlInt(angle, "p1"), fastXmlInt(angle, "p2"), fastXmlInt(angle, "p3"),
                            fastXmlDouble(angle, "a"), fastXmlDouble(angle, "k"));
            consumeLeafClose(reader, angle, "Angle");
        }
    }
}

static Force* parsePeriodicTorsionForceStreaming(FastXmlTagReader& reader, const FastXmlTag& opening) {
    const int version = fastXmlInt(opening, "version");
    if (version < 1 || version > 2)
        throw OpenMMException("BinarySerializer: unsupported PeriodicTorsionForce XML version");
    unique_ptr<PeriodicTorsionForce> force(new PeriodicTorsionForce());
    force->setForceGroup(fastXmlInt(opening, "forceGroup", 0));
    force->setName(fastXmlString(opening, "name", force->getName()));
    if (version > 1)
        force->setUsesPeriodicBoundaryConditions(fastXmlBool(opening, "usesPeriodic", false));
    bool haveTorsions = false;
    if (opening.empty)
        throw OpenMMException("BinarySerializer: PeriodicTorsionForce has no Torsions section");
    while (true) {
        FastXmlTag tag = parseFastXmlTag(reader.next());
        if (tag.closing) {
            if (!fastXmlNameEquals(tag, "Force"))
                throw OpenMMException("BinarySerializer: unexpected closing tag in PeriodicTorsionForce");
            if (!haveTorsions)
                throw OpenMMException("BinarySerializer: PeriodicTorsionForce has no Torsions section");
            return force.release();
        }
        if (!fastXmlNameEquals(tag, "Torsions") || haveTorsions)
            throw OpenMMException("BinarySerializer: unexpected PeriodicTorsionForce element");
        haveTorsions = true;
        if (tag.empty)
            continue;
        while (true) {
            FastXmlTag torsion = parseFastXmlTag(reader.next());
            if (torsion.closing) {
                if (!fastXmlNameEquals(torsion, "Torsions"))
                    throw OpenMMException("BinarySerializer: unexpected closing tag in Torsions");
                break;
            }
            if (!fastXmlNameEquals(torsion, "Torsion"))
                throw OpenMMException("BinarySerializer: unexpected element in PeriodicTorsionForce Torsions");
            force->addTorsion(fastXmlInt(torsion, "p1"), fastXmlInt(torsion, "p2"),
                fastXmlInt(torsion, "p3"), fastXmlInt(torsion, "p4"),
                fastXmlInt(torsion, "periodicity"), fastXmlDouble(torsion, "phase"),
                fastXmlDouble(torsion, "k"));
            consumeLeafClose(reader, torsion, "Torsion");
        }
    }
}

static Force* parseCMMotionRemoverStreaming(FastXmlTagReader& reader, const FastXmlTag& opening) {
    if (fastXmlInt(opening, "version") != 1)
        throw OpenMMException("BinarySerializer: unsupported CMMotionRemover XML version");
    unique_ptr<CMMotionRemover> force(new CMMotionRemover(fastXmlInt(opening, "frequency")));
    force->setForceGroup(fastXmlInt(opening, "forceGroup", 0));
    force->setName(fastXmlString(opening, "name", force->getName()));
    consumeLeafClose(reader, opening, "Force");
    return force.release();
}

static Force* parseMonteCarloBarostatStreaming(FastXmlTagReader& reader, const FastXmlTag& opening) {
    if (fastXmlInt(opening, "version") != 1)
        throw OpenMMException("BinarySerializer: unsupported MonteCarloBarostat XML version");
    unique_ptr<MonteCarloBarostat> force(new MonteCarloBarostat(
        fastXmlDouble(opening, "pressure"), fastXmlDouble(opening, "temperature"),
        fastXmlInt(opening, "frequency")));
    force->setForceGroup(fastXmlInt(opening, "forceGroup", 0));
    force->setName(fastXmlString(opening, "name", force->getName()));
    force->setRandomNumberSeed(fastXmlInt(opening, "randomSeed"));
    consumeLeafClose(reader, opening, "Force");
    return force.release();
}

static Force* parseForceStreaming(FastXmlTagReader& reader, const FastXmlTag& opening) {
    if (opening.closing || !fastXmlNameEquals(opening, "Force"))
        throw OpenMMException("BinarySerializer: malformed Force element");
    const string type = fastXmlString(opening, "type");
    if (type == "NonbondedForce")
        return parseNonbondedForceStreaming(reader, opening);
    if (type == "HarmonicBondForce")
        return parseHarmonicBondForceStreaming(reader, opening);
    if (type == "HarmonicAngleForce")
        return parseHarmonicAngleForceStreaming(reader, opening);
    if (type == "PeriodicTorsionForce")
        return parsePeriodicTorsionForceStreaming(reader, opening);
    if (type == "CMMotionRemover")
        return parseCMMotionRemoverStreaming(reader, opening);
    if (type == "MonteCarloBarostat")
        return parseMonteCarloBarostatStreaming(reader, opening);
    throw OpenMMException(string("BinarySerializer: unsupported Force type '")+type+"'");
}

static void parseSystemForcesStreaming(FastXmlTagReader& reader, const FastXmlTag& opening,
                                       System& system) {
    if (opening.closing || !fastXmlNameEquals(opening, "Forces"))
        throw OpenMMException("BinarySerializer: malformed System Forces section");
    if (opening.empty)
        return;
    while (true) {
        FastXmlTag tag = parseFastXmlTag(reader.next());
        if (tag.closing) {
            if (!fastXmlNameEquals(tag, "Forces"))
                throw OpenMMException("BinarySerializer: unexpected closing tag in Forces");
            return;
        }
        unique_ptr<Force> force(parseForceStreaming(reader, tag));
        system.addForce(force.release());
    }
}

static System* parseSystemXmlStreaming(istream& stream) {
    FastXmlTagReader reader(stream);
    FastXmlTag root = parseFastXmlTag(reader.next());
    if (root.closing || root.empty || !fastXmlNameEquals(root, "System"))
        throw OpenMMException("BinarySerializer: XML object is not a System");
    if (fastXmlString(root, "type") != "System")
        throw OpenMMException("BinarySerializer: XML type is not System");
    if (fastXmlInt(root, "version") != 1)
        throw OpenMMException("BinarySerializer: unsupported System XML version");

    unique_ptr<System> system(new System());
    bool haveBox = false, haveParticles = false, haveConstraints = false, haveForces = false;
    int stage = 0;
    while (true) {
        FastXmlTag tag = parseFastXmlTag(reader.next());
        if (tag.closing) {
            if (!fastXmlNameEquals(tag, "System"))
                throw OpenMMException("BinarySerializer: unexpected closing tag in System XML");
            break;
        }
        if (fastXmlNameEquals(tag, "PeriodicBoxVectors")) {
            if (stage != 0 || haveBox)
                throw OpenMMException("BinarySerializer: unexpected System PeriodicBoxVectors");
            parseSystemBoxStreaming(reader, tag, *system);
            haveBox = true;
            stage = 1;
        }
        else if (fastXmlNameEquals(tag, "Particles")) {
            if (!haveBox || stage != 1 || haveParticles)
                throw OpenMMException("BinarySerializer: unexpected System Particles section");
            parseSystemParticlesStreaming(reader, tag, *system);
            haveParticles = true;
            stage = 2;
        }
        else if (fastXmlNameEquals(tag, "Constraints")) {
            if (!haveParticles || stage != 2 || haveConstraints)
                throw OpenMMException("BinarySerializer: unexpected System Constraints section");
            parseSystemConstraintsStreaming(reader, tag, *system);
            haveConstraints = true;
            stage = 3;
        }
        else if (fastXmlNameEquals(tag, "Forces")) {
            if (!haveConstraints || stage != 3 || haveForces)
                throw OpenMMException("BinarySerializer: unexpected System Forces section");
            parseSystemForcesStreaming(reader, tag, *system);
            haveForces = true;
            stage = 4;
        }
        else
            throw OpenMMException("BinarySerializer: unsupported System XML element in fast path");
    }
    if (!haveBox || !haveParticles || !haveConstraints || !haveForces)
        throw OpenMMException("BinarySerializer: incomplete System XML");
    reader.finish();
    return system.release();
}

static void parseNonbondedSection(IrrXMLReader& xml,
                                  const char* sectionName,
                                  NonbondedForce& force) {
    // irrXML owns getNodeName() storage and may reuse it on the next read().
    // Copy the section name before advancing the parser.
    const string section(sectionName);

    enum SectionKind {
        GLOBAL_PARAMETERS,
        PARTICLE_OFFSETS,
        EXCEPTION_OFFSETS,
        PARTICLES,
        EXCEPTIONS
    };

    SectionKind kind;
    if (section == "GlobalParameters")
        kind = GLOBAL_PARAMETERS;
    else if (section == "ParticleOffsets")
        kind = PARTICLE_OFFSETS;
    else if (section == "ExceptionOffsets")
        kind = EXCEPTION_OFFSETS;
    else if (section == "Particles")
        kind = PARTICLES;
    else if (section == "Exceptions")
        kind = EXCEPTIONS;
    else
        throw OpenMMException(
            string("BinarySerializer: unknown NonbondedForce section ")+section);

    if (xml.isEmptyElement())
        return;

    while (xml.read()) {
        if (xml.getNodeType() == EXN_ELEMENT) {
            const char* name = xml.getNodeName();

            switch (kind) {
            case GLOBAL_PARAMETERS:
                if (strcmp(name, "Parameter") != 0)
                    throw OpenMMException(
                        "BinarySerializer: unexpected element in GlobalParameters");
                force.addGlobalParameter(
                    requiredXmlAttribute(xml, "name"),
                    xmlDouble(xml, "default"));
                break;

            case PARTICLE_OFFSETS:
                if (strcmp(name, "Offset") != 0)
                    throw OpenMMException(
                        "BinarySerializer: unexpected element in ParticleOffsets");
                force.addParticleParameterOffset(
                    requiredXmlAttribute(xml, "parameter"),
                    xmlInt(xml, "particle"),
                    xmlDouble(xml, "q"),
                    xmlDouble(xml, "sig"),
                    xmlDouble(xml, "eps"));
                break;

            case EXCEPTION_OFFSETS:
                if (strcmp(name, "Offset") != 0)
                    throw OpenMMException(
                        "BinarySerializer: unexpected element in ExceptionOffsets");
                force.addExceptionParameterOffset(
                    requiredXmlAttribute(xml, "parameter"),
                    xmlInt(xml, "exception"),
                    xmlDouble(xml, "q"),
                    xmlDouble(xml, "sig"),
                    xmlDouble(xml, "eps"));
                break;

            case PARTICLES:
                if (strcmp(name, "Particle") != 0)
                    throw OpenMMException(
                        "BinarySerializer: unexpected element in NonbondedForce Particles");
                force.addParticle(
                    xmlDouble(xml, "q"),
                    xmlDouble(xml, "sig"),
                    xmlDouble(xml, "eps"));
                break;

            case EXCEPTIONS:
                if (strcmp(name, "Exception") != 0)
                    throw OpenMMException(
                        "BinarySerializer: unexpected element in NonbondedForce Exceptions");
                force.addException(
                    xmlInt(xml, "p1"),
                    xmlInt(xml, "p2"),
                    xmlDouble(xml, "q"),
                    xmlDouble(xml, "sig"),
                    xmlDouble(xml, "eps"));
                break;
            }
        }
        else if (xml.getNodeType() == EXN_ELEMENT_END &&
                 strcmp(xml.getNodeName(), section.c_str()) == 0) {
            return;
        }
    }

    throw OpenMMException(
        string("BinarySerializer: truncated NonbondedForce section ")+section);
}

static Force* parseNonbondedForceXml(IrrXMLReader& xml) {
    int version = xmlInt(xml, "version");
    if (version < 1 || version > 4)
        throw OpenMMException(
            "BinarySerializer: unsupported NonbondedForce XML version");

    unique_ptr<NonbondedForce> force(new NonbondedForce());

    force->setForceGroup(xmlInt(xml, "forceGroup", 0));
    force->setName(xmlString(xml, "name", force->getName()));

    force->setNonbondedMethod(
        static_cast<NonbondedForce::NonbondedMethod>(
            xmlInt(xml, "method")));

    force->setCutoffDistance(xmlDouble(xml, "cutoff"));
    force->setUseSwitchingFunction(
        xmlBool(xml, "useSwitchingFunction", false));

    const char* switching = xml.getAttributeValue("switchingDistance");
    force->setSwitchingDistance(
        switching == NULL ? -1.0 : strtod2(switching, NULL));

    force->setEwaldErrorTolerance(xmlDouble(xml, "ewaldTolerance"));
    force->setReactionFieldDielectric(xmlDouble(xml, "rfDielectric"));
    force->setUseDispersionCorrection(
        xmlInt(xml, "dispersionCorrection") != 0);

    const char* includeDirect = xml.getAttributeValue("includeDirectSpace");
    if (includeDirect != NULL)
        force->setIncludeDirectSpace(
            strcmp(includeDirect, "1") == 0 ||
            strcmp(includeDirect, "true") == 0);

    double alpha = 0.0;
    int nx = 0, ny = 0, nz = 0;

    const char* alphaText = xml.getAttributeValue("alpha");
    if (alphaText != NULL)
        alpha = strtod2(alphaText, NULL);
    nx = xmlInt(xml, "nx", 0);
    ny = xmlInt(xml, "ny", 0);
    nz = xmlInt(xml, "nz", 0);
    force->setPMEParameters(alpha, nx, ny, nz);

    if (version >= 2) {
        alpha = 0.0;
        const char* ljAlpha = xml.getAttributeValue("ljAlpha");
        if (ljAlpha != NULL)
            alpha = strtod2(ljAlpha, NULL);

        nx = xmlInt(xml, "ljnx", 0);
        ny = xmlInt(xml, "ljny", 0);
        nz = xmlInt(xml, "ljnz", 0);
        force->setLJPMEParameters(alpha, nx, ny, nz);
    }

    force->setReciprocalSpaceForceGroup(
        xmlInt(xml, "recipForceGroup", -1));

    if (version >= 4)
        force->setExceptionsUsePeriodicBoundaryConditions(
            xmlInt(xml, "exceptionsUsePeriodic") != 0);

    bool haveGlobalParameters = false;
    bool haveParticleOffsets = false;
    bool haveExceptionOffsets = false;
    bool haveParticles = false;
    bool haveExceptions = false;

    while (xml.read()) {
        if (xml.getNodeType() == EXN_ELEMENT) {
            const char* name = xml.getNodeName();

            if (strcmp(name, "GlobalParameters") == 0) {
                if (version < 3)
                    throw OpenMMException(
                        "BinarySerializer: unexpected GlobalParameters");
                parseNonbondedSection(xml, name, *force);
                haveGlobalParameters = true;
            }
            else if (strcmp(name, "ParticleOffsets") == 0) {
                if (version < 3)
                    throw OpenMMException(
                        "BinarySerializer: unexpected ParticleOffsets");
                parseNonbondedSection(xml, name, *force);
                haveParticleOffsets = true;
            }
            else if (strcmp(name, "ExceptionOffsets") == 0) {
                if (version < 3)
                    throw OpenMMException(
                        "BinarySerializer: unexpected ExceptionOffsets");
                parseNonbondedSection(xml, name, *force);
                haveExceptionOffsets = true;
            }
            else if (strcmp(name, "Particles") == 0) {
                parseNonbondedSection(xml, name, *force);
                haveParticles = true;
            }
            else if (strcmp(name, "Exceptions") == 0) {
                parseNonbondedSection(xml, name, *force);
                haveExceptions = true;
            }
            else {
                throw OpenMMException(
                    string("BinarySerializer: unexpected NonbondedForce element '")+
                    name+"'");
            }
        }
        else if (xml.getNodeType() == EXN_ELEMENT_END &&
                 strcmp(xml.getNodeName(), "Force") == 0) {

            if (!haveParticles || !haveExceptions)
                throw OpenMMException(
                    "BinarySerializer: incomplete NonbondedForce");

            if (version >= 3 &&
                (!haveGlobalParameters ||
                 !haveParticleOffsets ||
                 !haveExceptionOffsets))
                throw OpenMMException(
                    "BinarySerializer: incomplete NonbondedForce parameter sections");

            return force.release();
        }
    }

    throw OpenMMException("BinarySerializer: truncated NonbondedForce");
}
static Force* parseHarmonicBondForceXml(IrrXMLReader& xml) {
    int version = xmlInt(xml, "version");
    if (version < 1 || version > 2)
        throw OpenMMException(
            "BinarySerializer: unsupported HarmonicBondForce XML version");

    unique_ptr<HarmonicBondForce> force(new HarmonicBondForce());

    force->setForceGroup(xmlInt(xml, "forceGroup", 0));
    force->setName(xmlString(xml, "name", force->getName()));

    if (version > 1)
        force->setUsesPeriodicBoundaryConditions(
            xmlBool(xml, "usesPeriodic", false));

    bool haveBonds = false;

    while (xml.read()) {
        if (xml.getNodeType() == EXN_ELEMENT) {
            const char* name = xml.getNodeName();

            if (strcmp(name, "Bonds") == 0) {
                haveBonds = true;

                if (xml.isEmptyElement())
                    continue;

                while (xml.read()) {
                    if (xml.getNodeType() == EXN_ELEMENT) {
                        if (strcmp(xml.getNodeName(), "Bond") != 0)
                            throw OpenMMException(
                                "BinarySerializer: unexpected element in HarmonicBondForce Bonds");

                        force->addBond(
                            xmlInt(xml, "p1"),
                            xmlInt(xml, "p2"),
                            xmlDouble(xml, "d"),
                            xmlDouble(xml, "k"));
                    }
                    else if (xml.getNodeType() == EXN_ELEMENT_END &&
                             strcmp(xml.getNodeName(), "Bonds") == 0) {
                        break;
                    }
                }
            }
            else {
                throw OpenMMException(
                    string("BinarySerializer: unexpected HarmonicBondForce element '")+
                    name+"'");
            }
        }
        else if (xml.getNodeType() == EXN_ELEMENT_END &&
                 strcmp(xml.getNodeName(), "Force") == 0) {
            if (!haveBonds)
                throw OpenMMException(
                    "BinarySerializer: HarmonicBondForce has no Bonds section");
            return force.release();
        }
    }

    throw OpenMMException("BinarySerializer: truncated HarmonicBondForce");
}

static Force* parseHarmonicAngleForceXml(IrrXMLReader& xml) {
    int version = xmlInt(xml, "version");
    if (version < 1 || version > 2)
        throw OpenMMException(
            "BinarySerializer: unsupported HarmonicAngleForce XML version");

    unique_ptr<HarmonicAngleForce> force(new HarmonicAngleForce());

    force->setForceGroup(xmlInt(xml, "forceGroup", 0));
    force->setName(xmlString(xml, "name", force->getName()));

    if (version > 1)
        force->setUsesPeriodicBoundaryConditions(
            xmlBool(xml, "usesPeriodic", false));

    bool haveAngles = false;

    while (xml.read()) {
        if (xml.getNodeType() == EXN_ELEMENT) {
            const char* name = xml.getNodeName();

            if (strcmp(name, "Angles") == 0) {
                haveAngles = true;

                if (xml.isEmptyElement())
                    continue;

                while (xml.read()) {
                    if (xml.getNodeType() == EXN_ELEMENT) {
                        if (strcmp(xml.getNodeName(), "Angle") != 0)
                            throw OpenMMException(
                                "BinarySerializer: unexpected element in HarmonicAngleForce Angles");

                        force->addAngle(
                            xmlInt(xml, "p1"),
                            xmlInt(xml, "p2"),
                            xmlInt(xml, "p3"),
                            xmlDouble(xml, "a"),
                            xmlDouble(xml, "k"));
                    }
                    else if (xml.getNodeType() == EXN_ELEMENT_END &&
                             strcmp(xml.getNodeName(), "Angles") == 0) {
                        break;
                    }
                }
            }
            else {
                throw OpenMMException(
                    string("BinarySerializer: unexpected HarmonicAngleForce element '")+
                    name+"'");
            }
        }
        else if (xml.getNodeType() == EXN_ELEMENT_END &&
                 strcmp(xml.getNodeName(), "Force") == 0) {
            if (!haveAngles)
                throw OpenMMException(
                    "BinarySerializer: HarmonicAngleForce has no Angles section");
            return force.release();
        }
    }

    throw OpenMMException("BinarySerializer: truncated HarmonicAngleForce");
}

static Force* parsePeriodicTorsionForceXml(IrrXMLReader& xml) {
    int version = xmlInt(xml, "version");
    if (version < 1 || version > 2)
        throw OpenMMException(
            "BinarySerializer: unsupported PeriodicTorsionForce XML version");

    unique_ptr<PeriodicTorsionForce> force(new PeriodicTorsionForce());

    force->setForceGroup(xmlInt(xml, "forceGroup", 0));
    force->setName(xmlString(xml, "name", force->getName()));

    if (version > 1)
        force->setUsesPeriodicBoundaryConditions(
            xmlBool(xml, "usesPeriodic", false));

    bool haveTorsions = false;

    while (xml.read()) {
        if (xml.getNodeType() == EXN_ELEMENT) {
            const char* name = xml.getNodeName();

            if (strcmp(name, "Torsions") == 0) {
                haveTorsions = true;

                if (xml.isEmptyElement())
                    continue;

                while (xml.read()) {
                    if (xml.getNodeType() == EXN_ELEMENT) {
                        if (strcmp(xml.getNodeName(), "Torsion") != 0)
                            throw OpenMMException(
                                "BinarySerializer: unexpected element in PeriodicTorsionForce Torsions");

                        force->addTorsion(
                            xmlInt(xml, "p1"),
                            xmlInt(xml, "p2"),
                            xmlInt(xml, "p3"),
                            xmlInt(xml, "p4"),
                            xmlInt(xml, "periodicity"),
                            xmlDouble(xml, "phase"),
                            xmlDouble(xml, "k"));
                    }
                    else if (xml.getNodeType() == EXN_ELEMENT_END &&
                             strcmp(xml.getNodeName(), "Torsions") == 0) {
                        break;
                    }
                }
            }
            else {
                throw OpenMMException(
                    string("BinarySerializer: unexpected PeriodicTorsionForce element '")+
                    name+"'");
            }
        }
        else if (xml.getNodeType() == EXN_ELEMENT_END &&
                 strcmp(xml.getNodeName(), "Force") == 0) {
            if (!haveTorsions)
                throw OpenMMException(
                    "BinarySerializer: PeriodicTorsionForce has no Torsions section");
            return force.release();
        }
    }

    throw OpenMMException("BinarySerializer: truncated PeriodicTorsionForce");
}
static Force* parseCMMotionRemoverXml(IrrXMLReader& xml) {
    int version = xmlInt(xml, "version");
    if (version != 1)
        throw OpenMMException(
            "BinarySerializer: unsupported CMMotionRemover XML version");

    unique_ptr<CMMotionRemover> force(
        new CMMotionRemover(xmlInt(xml, "frequency")));

    force->setForceGroup(xmlInt(xml, "forceGroup", 0));
    force->setName(xmlString(xml, "name", force->getName()));

    if (xml.isEmptyElement())
        return force.release();

    while (xml.read()) {
        if (xml.getNodeType() == EXN_ELEMENT)
            throw OpenMMException(
                "BinarySerializer: unexpected child in CMMotionRemover");

        if (xml.getNodeType() == EXN_ELEMENT_END &&
            strcmp(xml.getNodeName(), "Force") == 0)
            return force.release();
    }

    throw OpenMMException("BinarySerializer: truncated CMMotionRemover");
}

static Force* parseMonteCarloBarostatXml(IrrXMLReader& xml) {
    int version = xmlInt(xml, "version");
    if (version != 1)
        throw OpenMMException(
            "BinarySerializer: unsupported MonteCarloBarostat XML version");

    unique_ptr<MonteCarloBarostat> force(
        new MonteCarloBarostat(
            xmlDouble(xml, "pressure"),
            xmlDouble(xml, "temperature"),
            xmlInt(xml, "frequency")));

    force->setForceGroup(xmlInt(xml, "forceGroup", 0));
    force->setName(xmlString(xml, "name", force->getName()));
    force->setRandomNumberSeed(xmlInt(xml, "randomSeed"));

    if (xml.isEmptyElement())
        return force.release();

    while (xml.read()) {
        if (xml.getNodeType() == EXN_ELEMENT)
            throw OpenMMException(
                "BinarySerializer: unexpected child in MonteCarloBarostat");

        if (xml.getNodeType() == EXN_ELEMENT_END &&
            strcmp(xml.getNodeName(), "Force") == 0)
            return force.release();
    }

    throw OpenMMException("BinarySerializer: truncated MonteCarloBarostat");
}
static Force* parseForceXml(IrrXMLReader& xml) {
    const char* type = requiredXmlAttribute(xml, "type");

    if (strcmp(type, "NonbondedForce") == 0)
        return parseNonbondedForceXml(xml);

    if (strcmp(type, "HarmonicBondForce") == 0)
        return parseHarmonicBondForceXml(xml);

    if (strcmp(type, "HarmonicAngleForce") == 0)
        return parseHarmonicAngleForceXml(xml);

    if (strcmp(type, "PeriodicTorsionForce") == 0)
        return parsePeriodicTorsionForceXml(xml);

    if (strcmp(type, "CMMotionRemover") == 0)
        return parseCMMotionRemoverXml(xml);

    if (strcmp(type, "MonteCarloBarostat") == 0)
        return parseMonteCarloBarostatXml(xml);

    throw OpenMMException(
        string("BinarySerializer: unsupported Force type '")+type+"'");
}
static void parseSystemBoxVectors(IrrXMLReader& xml, System& system) {
    Vec3 values[3];
    bool found[3] = {false, false, false};

    while (xml.read()) {
        if (xml.getNodeType() == EXN_ELEMENT) {
            const char* name = xml.getNodeName();

            int index = -1;
            if (strcmp(name, "A") == 0)
                index = 0;
            else if (strcmp(name, "B") == 0)
                index = 1;
            else if (strcmp(name, "C") == 0)
                index = 2;
            else
                throw OpenMMException(
                    "BinarySerializer: unexpected element in System PeriodicBoxVectors");

            values[index] = Vec3(
                xmlDouble(xml, "x"),
                xmlDouble(xml, "y"),
                xmlDouble(xml, "z"));
            found[index] = true;
        }
        else if (xml.getNodeType() == EXN_ELEMENT_END &&
                 strcmp(xml.getNodeName(), "PeriodicBoxVectors") == 0) {
            if (!found[0] || !found[1] || !found[2])
                throw OpenMMException(
                    "BinarySerializer: incomplete System PeriodicBoxVectors");

            system.setDefaultPeriodicBoxVectors(
                values[0], values[1], values[2]);
            return;
        }
    }

    throw OpenMMException(
        "BinarySerializer: truncated System PeriodicBoxVectors");
}

static void parseSystemVirtualSite(IrrXMLReader& xml,
                                   System& system,
                                   int particle) {
    string type = xml.getNodeName();

    if (type == "TwoParticleAverageSite") {
        system.setVirtualSite(
            particle,
            new TwoParticleAverageSite(
                xmlInt(xml, "p1"),
                xmlInt(xml, "p2"),
                xmlDouble(xml, "w1"),
                xmlDouble(xml, "w2")));
    }
    else if (type == "ThreeParticleAverageSite") {
        system.setVirtualSite(
            particle,
            new ThreeParticleAverageSite(
                xmlInt(xml, "p1"),
                xmlInt(xml, "p2"),
                xmlInt(xml, "p3"),
                xmlDouble(xml, "w1"),
                xmlDouble(xml, "w2"),
                xmlDouble(xml, "w3")));
    }
    else if (type == "OutOfPlaneSite") {
        system.setVirtualSite(
            particle,
            new OutOfPlaneSite(
                xmlInt(xml, "p1"),
                xmlInt(xml, "p2"),
                xmlInt(xml, "p3"),
                xmlDouble(xml, "w12"),
                xmlDouble(xml, "w13"),
                xmlDouble(xml, "wc")));
    }
    else if (type == "LocalCoordinatesSite") {
        vector<int> particles;
        vector<double> wo, wx, wy;

        for (int j = 1; ; j++) {
            stringstream indexStream;
            indexStream << j;
            string index = indexStream.str();

            string pName = "p"+index;
            if (xml.getAttributeValue(pName.c_str()) == NULL)
                break;

            string woName = "wo"+index;
            string wxName = "wx"+index;
            string wyName = "wy"+index;

            particles.push_back(xmlInt(xml, pName.c_str()));
            wo.push_back(xmlDouble(xml, woName.c_str()));
            wx.push_back(xmlDouble(xml, wxName.c_str()));
            wy.push_back(xmlDouble(xml, wyName.c_str()));
        }

        if (particles.empty())
            throw OpenMMException(
                "BinarySerializer: LocalCoordinatesSite has no particles");

        Vec3 position(
            xmlDouble(xml, "pos1"),
            xmlDouble(xml, "pos2"),
            xmlDouble(xml, "pos3"));

        system.setVirtualSite(
            particle,
            new LocalCoordinatesSite(
                particles, wo, wx, wy, position));
    }
    else {
        throw OpenMMException(
            string("BinarySerializer: unsupported virtual site type '")+
            type+"'");
    }
}

static void parseSystemParticles(IrrXMLReader& xml, System& system) {
    if (xml.isEmptyElement())
        return;

    while (xml.read()) {
        if (xml.getNodeType() == EXN_ELEMENT) {
            if (strcmp(xml.getNodeName(), "Particle") != 0)
                throw OpenMMException(
                    "BinarySerializer: unexpected element in System Particles");

            int particle = system.getNumParticles();
            system.addParticle(xmlDouble(xml, "mass"));

            if (xml.isEmptyElement())
                continue;

            bool haveVirtualSite = false;

            while (xml.read()) {
                if (xml.getNodeType() == EXN_ELEMENT) {
                    if (haveVirtualSite)
                        throw OpenMMException(
                            "BinarySerializer: Particle has multiple virtual sites");

                    parseSystemVirtualSite(xml, system, particle);
                    haveVirtualSite = true;
                }
                else if (xml.getNodeType() == EXN_ELEMENT_END &&
                         strcmp(xml.getNodeName(), "Particle") == 0) {
                    break;
                }
            }
        }
        else if (xml.getNodeType() == EXN_ELEMENT_END &&
                 strcmp(xml.getNodeName(), "Particles") == 0) {
            return;
        }
    }

    throw OpenMMException(
        "BinarySerializer: truncated System Particles section");
}

static void parseSystemConstraints(IrrXMLReader& xml, System& system) {
    if (xml.isEmptyElement())
        return;

    while (xml.read()) {
        if (xml.getNodeType() == EXN_ELEMENT) {
            if (strcmp(xml.getNodeName(), "Constraint") != 0)
                throw OpenMMException(
                    "BinarySerializer: unexpected element in System Constraints");

            system.addConstraint(
                xmlInt(xml, "p1"),
                xmlInt(xml, "p2"),
                xmlDouble(xml, "d"));
        }
        else if (xml.getNodeType() == EXN_ELEMENT_END &&
                 strcmp(xml.getNodeName(), "Constraints") == 0) {
            return;
        }
    }

    throw OpenMMException(
        "BinarySerializer: truncated System Constraints section");
}

static void parseSystemForces(IrrXMLReader& xml, System& system) {
    if (xml.isEmptyElement())
        return;

    while (xml.read()) {
        if (xml.getNodeType() == EXN_ELEMENT) {
            if (strcmp(xml.getNodeName(), "Force") != 0)
                throw OpenMMException(
                    "BinarySerializer: unexpected element in System Forces");

            unique_ptr<Force> force(parseForceXml(xml));
            system.addForce(force.release());
        }
        else if (xml.getNodeType() == EXN_ELEMENT_END &&
                 strcmp(xml.getNodeName(), "Forces") == 0) {
            return;
        }
    }

    throw OpenMMException(
        "BinarySerializer: truncated System Forces section");
}

static System* parseSystemXml(istream& stream) {
    FastXmlStreamReader callback(stream);
    unique_ptr<IrrXMLReader> xml(createIrrXMLReader(&callback));

    if (!xml)
        throw OpenMMException(
            "BinarySerializer: could not create XML reader");

    while (xml->read() &&
           xml->getNodeType() != EXN_ELEMENT)
        ;

    if (xml->getNodeType() != EXN_ELEMENT ||
        strcmp(xml->getNodeName(), "System") != 0)
        throw OpenMMException(
            "BinarySerializer: XML object is not a System");

    const char* type = xml->getAttributeValue("type");
    if (type == NULL || strcmp(type, "System") != 0)
        throw OpenMMException(
            "BinarySerializer: XML type is not System");

    if (xmlInt(*xml, "version") != 1)
        throw OpenMMException(
            "BinarySerializer: unsupported System XML version");

    unique_ptr<System> system(new System());

    bool haveBox = false;
    bool haveParticles = false;
    bool haveConstraints = false;
    bool haveForces = false;

    while (xml->read()) {
        if (xml->getNodeType() == EXN_ELEMENT) {
            string name = xml->getNodeName();

            if (name == "PeriodicBoxVectors") {
                if (haveBox)
                    throw OpenMMException(
                        "BinarySerializer: duplicate PeriodicBoxVectors");
                parseSystemBoxVectors(*xml, *system);
                haveBox = true;
            }
            else if (name == "Particles") {
                if (haveParticles)
                    throw OpenMMException(
                        "BinarySerializer: duplicate Particles section");
                parseSystemParticles(*xml, *system);
                haveParticles = true;
            }
            else if (name == "Constraints") {
                if (haveConstraints)
                    throw OpenMMException(
                        "BinarySerializer: duplicate Constraints section");
                parseSystemConstraints(*xml, *system);
                haveConstraints = true;
            }
            else if (name == "Forces") {
                if (haveForces)
                    throw OpenMMException(
                        "BinarySerializer: duplicate Forces section");
                parseSystemForces(*xml, *system);
                haveForces = true;
            }
            else {
                throw OpenMMException(
                    string("BinarySerializer: unexpected System element '")+
                    name+"'");
            }
        }
        else if (xml->getNodeType() == EXN_ELEMENT_END &&
                 strcmp(xml->getNodeName(), "System") == 0) {
            if (!haveBox ||
                !haveParticles ||
                !haveConstraints ||
                !haveForces)
                throw OpenMMException(
                    "BinarySerializer: incomplete System XML");

            return system.release();
        }
    }

    throw OpenMMException(
        "BinarySerializer: truncated System XML");
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
        return new State(builder.takeState());
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

void* BinarySerializer::deserializeStateXml(istream& stream) {
    return parseStateXmlStreaming(stream);
}
void* BinarySerializer::deserializeSystemXml(istream& stream) {
    return parseSystemXmlStreaming(stream);
}