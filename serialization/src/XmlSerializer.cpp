/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2010-2015 Stanford University and the Authors.      *
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

#include "openmm/serialization/XmlSerializer.h"
#include "openmm/serialization/StateProxy.h"
#include "openmm/serialization/SystemProxy.h"
#include "openmm/State.h"
#include "irrXML.h"
#include "FastStateXml.h"
#include <cstdlib>
#include <vector>
#include <cstring>
#include <iostream>
#include <map>
#include <cstdint>
#include <cstdio>
#include <atomic>
#if defined(_WIN32) && defined(_MSC_VER) && defined(_M_X64)
#define NOMINMAX
#define WIN32_LEAN_AND_MEAN
#include <windows.h>
#include <bcrypt.h>
#endif

using namespace OpenMM;
using namespace std;
using namespace irr;
using namespace io;

extern "C" char* g_fmt(char*, double);
static bool fastXmlEnabledByEnvironment();

namespace {
const char* core28Cookie = "__OpenMMCore28State";
struct DeferredState {
    const State* state;
    const SerializationNode* node;
    uint64_t generation;
};
thread_local DeferredState deferredState = {};
thread_local uint64_t stateGeneration = 0;
atomic<bool> stateProxyUnchanged(true);

#if defined(_WIN32) && defined(_MSC_VER) && defined(_M_X64)
// This is a fixed-binary adapter, not an alternative SerializationNode format.
struct Core28Host {
    const unsigned char* base;
    unsigned char wrapper[0xafb9-0xad80];
};

static bool core28Hash(HANDLE file, const unsigned char* data, ULONG size, unsigned char digest[32]) {
    BCRYPT_ALG_HANDLE algorithm = NULL;
    BCRYPT_HASH_HANDLE hash = NULL;
    unsigned char object[1024], buffer[8192];
    ULONG length = 0, returned = 0;
    bool ok = false;
    do {
        if (BCryptOpenAlgorithmProvider(&algorithm, BCRYPT_SHA256_ALGORITHM, NULL, 0) < 0 ||
                BCryptGetProperty(algorithm, BCRYPT_OBJECT_LENGTH, reinterpret_cast<PUCHAR>(&length),
                                  sizeof(length), &returned, 0) < 0 || length > sizeof(object) ||
                BCryptCreateHash(algorithm, &hash, object, length, NULL, 0, 0) < 0)
            break;
        if (file != INVALID_HANDLE_VALUE) {
            DWORD count;
            while ((ok = ReadFile(file, buffer, sizeof(buffer), &count, NULL) != 0) && count != 0)
                if (BCryptHashData(hash, buffer, count, 0) < 0) { ok = false; break; }
        }
        else
            ok = BCryptHashData(hash, const_cast<PUCHAR>(data), size, 0) >= 0;
        if (ok)
            ok = BCryptFinishHash(hash, digest, 32, 0) >= 0;
    } while (false);
    if (hash != NULL) BCryptDestroyHash(hash);
    if (algorithm != NULL) BCryptCloseAlgorithmProvider(algorithm, 0);
    return ok;
}

static bool core28DigestMatches(const unsigned char digest[32], const char* hex) {
    const char* digits = "0123456789abcdef";
    for (int i = 0; i < 32; i++)
        if (digits[digest[i]>>4] != hex[2*i] || digits[digest[i]&15] != hex[2*i+1])
            return false;
    return true;
}

static Core28Host inspectCore28Host() {
    Core28Host host = {};
    wchar_t path[MAX_PATH];
    DWORD length = GetModuleFileNameW(NULL, path, sizeof(path)/sizeof(path[0]));
    if (length == 0 || length >= sizeof(path)/sizeof(path[0]))
        return host;
    HANDLE file = CreateFileW(path, GENERIC_READ, FILE_SHARE_READ | FILE_SHARE_DELETE,
                              NULL, OPEN_EXISTING, FILE_ATTRIBUTE_NORMAL, NULL);
    if (file == INVALID_HANDLE_VALUE)
        return host;
    LARGE_INTEGER size;
    unsigned char digest[32];
    bool matched = GetFileSizeEx(file, &size) && size.QuadPart == 2507776 &&
        core28Hash(file, NULL, 0, digest) && core28DigestMatches(digest,
        "9e27f4493c83a31e8fb0403248629909a67f8fc9c3a7ac6a41afa6dc93133436");
    CloseHandle(file);
    if (!matched)
        return host;
    const unsigned char* base = reinterpret_cast<const unsigned char*>(GetModuleHandleW(NULL));
    const IMAGE_DOS_HEADER* dos = reinterpret_cast<const IMAGE_DOS_HEADER*>(base);
    if (dos->e_magic != IMAGE_DOS_SIGNATURE || dos->e_lfanew < 0 || dos->e_lfanew > 4096)
        return host;
    const IMAGE_NT_HEADERS64* pe = reinterpret_cast<const IMAGE_NT_HEADERS64*>(base+dos->e_lfanew);
    if (pe->Signature != IMAGE_NT_SIGNATURE || pe->FileHeader.Machine != IMAGE_FILE_MACHINE_AMD64 ||
            pe->OptionalHeader.SizeOfImage < 0xafb9 ||
            !core28Hash(INVALID_HANDLE_VALUE, base+0xad80, sizeof(host.wrapper), digest) ||
            !core28DigestMatches(digest,
            "b3d29e3afd1715820203cee78d004b3ce8a48712d6e810dbe96ef454b2a226f9"))
        return host;
    memcpy(host.wrapper, base+0xad80, sizeof(host.wrapper));
    host.base = base;
    return host;
}
#endif

static bool core28Caller(const void* caller, size_t rva) {
#if defined(_WIN32) && defined(_MSC_VER) && defined(_M_X64)
    const unsigned char* base = reinterpret_cast<const unsigned char*>(GetModuleHandleW(NULL));
    if (caller != base+rva)
        return false;
    static const Core28Host host = inspectCore28Host();
    return host.base == base &&
        memcmp(host.wrapper, host.base+0xad80, sizeof(host.wrapper)) == 0;
#else
    return false;
#endif
}
} // namespace

bool OpenMM::canDeferCore28StateXml(const State& state, const StateProxy& proxy,
        const SerializationNode& node, const void* caller) {
    return stateProxyUnchanged.load(memory_order_relaxed) && fastXmlEnabledByEnvironment() &&
        deferredState.node == NULL && core28Caller(caller, 0xadee) &&
        typeid(proxy) == typeid(StateProxy) && state.getDataTypes() == 0x1f &&
        node.getName() == "State" && node.getChildren().empty() && node.getProperties().empty() &&
        !state.getPositions().empty() && state.getPositions().size() == state.getVelocities().size() &&
        state.getPositions().size() == state.getForces().size();
}

bool OpenMM::deferCore28StateXml(const State& state, SerializationNode& node) {
    if (deferredState.node != NULL)
        return false;
    if (++stateGeneration == 0) ++stateGeneration;
    const uint64_t generation = stateGeneration;
    char cookie[32];
    snprintf(cookie, sizeof(cookie), "%llu", static_cast<unsigned long long>(generation));
    node.setStringProperty(core28Cookie, cookie);
    // Publish last: exceptions above leave no borrowed pointer behind.
    deferredState = {&state, &node, generation};
    return true;
}

void OpenMM::cancelCore28StateXml(const SerializationNode* node) noexcept {
    if (deferredState.node == node)
        deferredState = {};
}

void OpenMM::disableCore28StateXml() noexcept {
    stateProxyUnchanged.store(false, memory_order_relaxed);
}

/**
 * Apply XML encoding to a string.  This is adapted from TinyXML (written by Lee Thomason).
 */
static void encodeString(const string& str, string* outString) {
    int i=0;

    while (i<(int)str.length()) {
        unsigned char c = (unsigned char) str[i];
        const char* entity = NULL;
        switch (c) {
            case '&': entity = "&amp;"; break;
            case '<': entity = "&lt;"; break;
            case '>': entity = "&gt;"; break;
            case '"': entity = "&quot;"; break;
            case '\'': entity = "&apos;"; break;
        }

        if (c == '&' 
             && i < ((int)str.length() - 2)
             && str[i+1] == '#'
             && str[i+2] == 'x') {
            // Hexadecimal character reference.
            // Pass through unchanged.
            // &#xA9;	-- copyright symbol, for example.
            //
            // The -1 is a bug fix from Rob Laveaux. It keeps
            // an overflow from happening if there is no ';'.
            // There are actually 2 ways to exit this loop -
            // while fails (error case) and break (semicolon found).
            // However, there is no mechanism (currently) for
            // this function to return an error.
            while (i<(int)str.length()-1) {
                outString->append(str.c_str() + i, 1);
                ++i;
                if (str[i] == ';')
                    break;
            }
        }
        else if (entity != NULL) {
            outString->append(entity);
            ++i;
        }
        else if (c < 32) {
            // Easy pass at non-alpha/numeric/symbol
            // Below 32 is symbolic.
            char buf[ 32 ];

            sprintf(buf, "&#x%02X;", (unsigned) (c & 0xff));

            //*ME:	warning C4267: convert 'size_t' to 'int'
            //*ME:	Int-Cast to make compiler happy ...
            outString->append(buf, (int)strlen(buf));
            ++i;
        }
        else {
            //char realc = (char) c;
            //outString->append(&realc, 1);
            *outString += (char) c;	// somewhat more efficient function call.
            ++i;
        }
    }
}

namespace {
// Reuse the bounded writer from the fast-XML work, with ordinary ostream
// writes so stream state, exceptions and the caller's compression stay in charge.
class FastXmlOutput {
public:
    explicit FastXmlOutput(ostream& stream) : stream(stream) {
        if (stream.good())
            buffer.reserve(limit);
        else
            stream.write("", 0); // Apply the ordinary ostream sentry and exception policy.
    }
    bool good() const { return stream.good(); }
    void text(const char* value) { append(value, strlen(value)); }
    void text(const string& value) { append(value.data(), value.size()); }
    void real(double value) { char formatted[32]; g_fmt(formatted, value); text(formatted); }
    void encoded(const string& value) { string encoded; encodeString(value, &encoded); text(encoded); }
    void flush() {
        if (!buffer.empty()) {
            stream.write(buffer.data(), buffer.size());
            buffer.clear();
        }
    }
private:
    enum { limit = 4*1024*1024 };
    ostream& stream;
    string buffer;
    void append(const char* value, size_t size) {
        if (!stream.good()) return;
        while (size != 0) {
            if (buffer.size() == limit) { flush(); if (!stream.good()) return; }
            size_t count = min(size, limit-buffer.size());
            buffer.append(value, count); value += count; size -= count;
        }
    }
};

static void writeVec3Section(FastXmlOutput& out, const char* name, const char* element, const vector<Vec3>& values) {
    out.text("\t<"); out.text(name);
    if (values.empty()) { out.text("/>\n"); return; }
    out.text(">\n");
    for (const Vec3& value : values) {
        if (!out.good()) break;
        out.text("\t\t<"); out.text(element);
        out.text(" x=\""); out.real(value[0]); out.text("\" y=\""); out.real(value[1]);
        out.text("\" z=\""); out.real(value[2]); out.text("\"/>\n");
    }
    out.text("\t</"); out.text(name); out.text(">\n");
}

static void serializeStateDirect(const SerializationNode& node, const State& state, ostream& stream) {
    FastXmlOutput out(stream);
    if (!out.good()) return;
    out.text("<?xml version=\"1.0\" ?>\n<"); out.text(node.getName());
    for (const auto& prop : node.getProperties()) {
        if (prop.first == core28Cookie) continue;
        out.text(" "); out.encoded(prop.first); out.text("=\""); out.encoded(prop.second); out.text("\"");
    }
    out.text(">\n\t<PeriodicBoxVectors>\n");
    Vec3 box[3]; state.getPeriodicBoxVectors(box[0], box[1], box[2]);
    const char* names[] = {"A", "B", "C"};
    for (int i = 0; i < 3; i++) {
        out.text("\t\t<"); out.text(names[i]); out.text(" x=\""); out.real(box[i][0]);
        out.text("\" y=\""); out.real(box[i][1]); out.text("\" z=\""); out.real(box[i][2]); out.text("\"/>\n");
    }
    out.text("\t</PeriodicBoxVectors>\n");
    const int types = state.getDataTypes();
    if (types & State::Parameters) {
        out.text("\t<Parameters");
        for (const auto& param : state.getParameters()) {
            out.text(" "); out.encoded(param.first); out.text("=\""); out.real(param.second); out.text("\"");
        }
        out.text("/>\n");
    }
    if (types & State::Energy) {
        out.text("\t<Energies KineticEnergy=\""); out.real(state.getKineticEnergy());
        out.text("\" PotentialEnergy=\""); out.real(state.getPotentialEnergy()); out.text("\"/>\n");
    }
    if (types & State::Positions) writeVec3Section(out, "Positions", "Position", state.getPositions());
    if (types & State::Velocities) writeVec3Section(out, "Velocities", "Velocity", state.getVelocities());
    if (types & State::Forces) writeVec3Section(out, "Forces", "Force", state.getForces());
    out.text("</"); out.text(node.getName()); out.text(">\n"); out.flush();
}
} // namespace

static const State* takeCore28StateXml(const SerializationNode& node) {
    if (deferredState.node != &node)
        return NULL; // Unrelated/nested nodes, including user cookie properties, are ordinary trees.
    DeferredState pending = deferredState;
    deferredState = {}; // Clear before any validation, allocation or output can throw.
    char cookie[32];
    snprintf(cookie, sizeof(cookie), "%llu", static_cast<unsigned long long>(pending.generation));
    if (node.getName() != "State" ||
            node.getStringProperty("type", "") != "State" || node.getStringProperty(core28Cookie, "") != cookie)
        throw OpenMMException("XmlSerializer: invalid deferred Core28 State");
    return pending.state;
}

bool OpenMM::writeCore28StateXml(const SerializationNode& node, ostream& stream, const void* caller) {
    const State* state = takeCore28StateXml(node);
    if (state == NULL)
        return false;
    if (!core28Caller(caller, 0xaf1f))
        throw OpenMMException("XmlSerializer: invalid deferred Core28 caller");
    serializeStateDirect(node, *state, stream);
    return true;
}

void XmlSerializer::serialize(const SerializationNode& node, std::ostream& stream) {
    if (writeCore28StateXml(node, stream, OPENMM_XML_CALLER()))
        return;
    stream << "<?xml version=\"1.0\" ?>\n";
    encodeNode(node, stream, 0);
}

void XmlSerializer::encodeNode(const SerializationNode& node, std::ostream& stream, int depth) {
    for (int i = 0; i < depth; i++)
        stream << '\t';
    stream << '<' << node.getName();
    for (auto& prop : node.getProperties()) {
        string name, value;
        encodeString(prop.first, &name);
        encodeString(prop.second, &value);
        stream << ' ' << name << "=\"" << value << '\"';
    }
    const vector<SerializationNode>& children = node.getChildren();
    if (children.size() == 0)
        stream << "/>\n";
    else {
        stream << ">\n";
        for (auto& child : children)
            encodeNode(child, stream, depth+1);
        for (int i = 0; i < depth; i++)
            stream << '\t';
        stream << "</" << node.getName() << ">\n";
    }
}

/**
 * Adapter class to let irrXML read a C++ stream.
 */
class XmlSerializer::StreamReader : public IFileReadCallBack {
public:
    StreamReader(std::istream& stream) : stream(stream) {
        stream.seekg(0, ios_base::end);
        size = stream.tellg();
        stream.seekg(0);
    }
    int read(void* buffer, int sizeToRead) {
        stream.read((char*) buffer, sizeToRead);
        return stream.gcount();
    }
    int getSize() {
        return size;
    }
private:
    std::istream& stream;
    int size;
};

/**
 * Process an XML node, storing its content into a SerializationNode.
 */
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

/**
 * Safely inspect a bounded prefix of a seekable XML stream.
 *
 * irrXML reads the complete size advertised by IFileReadCallBack into one
 * buffer before parsing it.  Feeding it a truncated prefix is therefore not
 * safe: a construct such as a comment that begins inside the prefix and ends
 * after it can make the bundled parser scan beyond that buffer.
 *
 * This scanner deliberately understands only enough XML syntax to identify a
 * simple canonical OpenMM root element and its type attribute.  It bounds
 * every access.  Any declaration it does not explicitly understand, any
 * incomplete token at the prefix boundary, or any unusual encoding simply
 * returns false so the caller uses OpenMM's stock XML parser.
 */
static bool isProbeNameChar(char c) {
    return (c >= 'a' && c <= 'z') ||
           (c >= 'A' && c <= 'Z') ||
           (c >= '0' && c <= '9') ||
           c == '_' || c == ':' || c == '-' || c == '.';
}

static bool isProbeSpace(char c) {
    return c == ' ' || c == '\t' || c == '\r' || c == '\n';
}

static void rewindXmlStreamOrThrow(std::istream& stream, std::streampos start) {
    // Rewinding is part of the compatibility contract for the opportunistic
    // fast path.  Do it with exceptions temporarily disabled, then normalize
    // any seek failure to OpenMMException instead of entering the stock parser
    // with a failed stream.
    const std::ios_base::iostate exceptionMask = stream.exceptions();
    stream.exceptions(std::ios_base::goodbit);

    bool rewound = false;
    try {
        stream.clear();
        stream.seekg(start);
        rewound = static_cast<bool>(stream);
    }
    catch (...) {
        rewound = false;
    }

    // Restore the caller's exception mask only from a good state so
    // exceptions(mask) itself cannot throw because of the failed rewind.
    if (!rewound)
        stream.clear();
    stream.exceptions(exceptionMask);

    if (!rewound)
        throw OpenMMException("XmlSerializer: failed to rewind XML stream");
}

static bool startsWithAt(const vector<char>& data, size_t pos, const char* text) {
    const size_t length = strlen(text);
    return pos <= data.size() && length <= data.size()-pos &&
           memcmp(&data[pos], text, length) == 0;
}

static bool skipUntil(const vector<char>& data, size_t& pos, const char* terminator) {
    const size_t length = strlen(terminator);
    while (pos <= data.size()) {
        if (length <= data.size()-pos &&
                memcmp(&data[pos], terminator, length) == 0) {
            pos += length;
            return true;
        }
        if (pos == data.size())
            break;
        ++pos;
    }
    return false;
}

static bool probeXmlRoot(std::istream& stream, std::streampos start,
                         string& rootName, string& type) {
    const size_t probeLimit = 64*1024;
    vector<char> data(probeLimit);

    // A normal short read sets eofbit/failbit.  Temporarily disable stream
    // exceptions so probing a valid document shorter than probeLimit cannot
    // throw when the caller has enabled failbit/badbit exceptions.  Restore
    // both the original position and exception mask before returning.
    const std::ios_base::iostate exceptionMask = stream.exceptions();
    stream.exceptions(std::ios_base::goodbit);

    stream.clear();
    stream.seekg(start);
    if (!stream) {
        stream.clear();
        stream.exceptions(exceptionMask);
        return false;
    }

    stream.read(&data[0], static_cast<streamsize>(data.size()));
    const streamsize count = stream.gcount();

    // Reaching EOF during the bounded probe is normal.  Clear it before
    // rewinding, then restore the caller's exception mask only after the stream
    // is back in a good state so exceptions(mask) itself cannot throw.
    stream.clear();
    stream.seekg(start);
    const bool rewound = static_cast<bool>(stream);
    stream.clear();
    stream.exceptions(exceptionMask);
    if (!rewound || count <= 0)
        return false;

    data.resize(static_cast<size_t>(count));
    size_t pos = 0;

    // UTF-8 BOM is safe to skip.  Let the stock parser handle UTF-16/UTF-32 or
    // any other encoding whose byte representation this small scanner does not
    // understand.
    if (data.size() >= 3 &&
            static_cast<unsigned char>(data[0]) == 0xEF &&
            static_cast<unsigned char>(data[1]) == 0xBB &&
            static_cast<unsigned char>(data[2]) == 0xBF)
        pos = 3;
    else if ((data.size() >= 2 &&
              ((static_cast<unsigned char>(data[0]) == 0xFF &&
                static_cast<unsigned char>(data[1]) == 0xFE) ||
               (static_cast<unsigned char>(data[0]) == 0xFE &&
                static_cast<unsigned char>(data[1]) == 0xFF))) ||
             (data.size() >= 4 &&
              ((static_cast<unsigned char>(data[0]) == 0x00 &&
                static_cast<unsigned char>(data[1]) == 0x00 &&
                static_cast<unsigned char>(data[2]) == 0xFE &&
                static_cast<unsigned char>(data[3]) == 0xFF) ||
               (static_cast<unsigned char>(data[0]) == 0xFF &&
                static_cast<unsigned char>(data[1]) == 0xFE &&
                static_cast<unsigned char>(data[2]) == 0x00 &&
                static_cast<unsigned char>(data[3]) == 0x00))))
        return false;

    // Skip whitespace, XML declarations/processing instructions, and complete
    // comments.  For DOCTYPE or any other <! ...> declaration, conservatively
    // fall back rather than implementing more XML grammar here.
    while (true) {
        while (pos < data.size() && isProbeSpace(data[pos]))
            ++pos;
        if (pos == data.size())
            return false;

        if (startsWithAt(data, pos, "<?")) {
            pos += 2;
            if (!skipUntil(data, pos, "?>"))
                return false;
            continue;
        }

        if (startsWithAt(data, pos, "<!--")) {
            pos += 4;
            if (!skipUntil(data, pos, "-->"))
                return false;
            continue;
        }

        if (startsWithAt(data, pos, "<!"))
            return false;

        break;
    }

    if (pos >= data.size() || data[pos] != '<')
        return false;
    ++pos;

    if (pos >= data.size() || data[pos] == '/' || data[pos] == '?' ||
            data[pos] == '!')
        return false;

    const size_t rootStart = pos;
    while (pos < data.size() && isProbeNameChar(data[pos]))
        ++pos;
    if (pos == rootStart)
        return false;
    rootName.assign(&data[rootStart], pos-rootStart);

    // We only need the fast path for the two canonical OpenMM root names.  A
    // custom root can immediately use the stock path without scanning its
    // attributes.
    if (rootName != "State" && rootName != "System")
        return true;

    while (true) {
        while (pos < data.size() && isProbeSpace(data[pos]))
            ++pos;
        if (pos == data.size())
            return false;  // root start tag crosses the probe boundary

        if (data[pos] == '>') {
            ++pos;
            return true;
        }
        if (data[pos] == '/') {
            if (pos+1 >= data.size())
                return false;
            if (data[pos+1] != '>')
                return false;
            pos += 2;
            return true;
        }

        const size_t nameStart = pos;
        while (pos < data.size() && isProbeNameChar(data[pos]))
            ++pos;
        if (pos == nameStart)
            return false;
        const string attributeName(&data[nameStart], pos-nameStart);

        while (pos < data.size() && isProbeSpace(data[pos]))
            ++pos;
        if (pos == data.size() || data[pos] != '=')
            return false;
        ++pos;

        while (pos < data.size() && isProbeSpace(data[pos]))
            ++pos;
        if (pos == data.size() || (data[pos] != '"' && data[pos] != '\''))
            return false;

        const char quote = data[pos++];
        const size_t valueStart = pos;
        while (pos < data.size() && data[pos] != quote)
            ++pos;
        if (pos == data.size())
            return false;  // attribute value crosses the probe boundary

        if (attributeName == "type")
            type.assign(&data[valueStart], pos-valueStart);
        ++pos;
    }
}

/**
 * Process an XML node, storing its content into a SerializationNode.
 */
static void decodeNode(SerializationNode& node, IrrXMLReader& xml) {
    for (int i = 0; i < xml.getAttributeCount(); i++)
        node.setStringProperty(xml.getAttributeName(i), xml.getAttributeValue(i));
    if (xml.isEmptyElement())
        return;
    while (xml.read()) {
        switch (xml.getNodeType()) {
            case EXN_ELEMENT:
            {
                SerializationNode& childNode = node.createChildNode(xml.getNodeName());
                decodeNode(childNode, xml);
                break;
            }
            case EXN_ELEMENT_END:
                return;
        }
    }
}

void* XmlSerializer::deserializeStream(std::istream& stream) {
    // Opportunistic fast path for canonical State/System XML.  It must never
    // narrow XmlSerializer's accepted input: custom root names, unsupported
    // force/plugin types, IntegratorParameters, or any other schema feature
    // understood by the registered SerializationProxy path fall back to the
    // original reader below.
    std::streampos start = std::streampos(-1);
    if (fastXmlEnabledByEnvironment() && stream.rdbuf() != NULL) {
        try {
            start = stream.rdbuf()->pubseekoff(0, ios_base::cur, ios_base::in);
        }
        catch (const ios_base::failure&) {
            start = std::streampos(-1);
        }
    }
    if (start != std::streampos(-1)) {
        // Inspect at most 64 KiB with a bounds-checked scanner.  Do not
        // hand a truncated prefix to irrXML: the bundled parser assumes its
        // advertised buffer contains complete constructs and can scan beyond a
        // truncated comment or processing instruction.
        string rootName;
        string type;
        const bool haveRoot = probeXmlRoot(stream, start, rootName, type);

        rewindXmlStreamOrThrow(stream, start);

        try {
            if (haveRoot && rootName == "State" && type == "State" &&
                    typeid(SerializationProxy::getProxy(type)) == typeid(StateProxy))
                return deserializeStateXmlFast(stream);

            if (haveRoot && rootName == "System" && type == "System" &&
                    typeid(SerializationProxy::getProxy(type)) == typeid(SystemProxy))
                return deserializeSystemXmlFast(stream);
        }
        catch (const OpenMMException&) {
            // The fast readers intentionally support only the schemas that are
            // profitable to specialize.  Fall through to the checked rewind
            // below and preserve the behavior of the existing generic
            // SerializationProxy reader.  If the document is genuinely
            // invalid, that reader will report the normal deserialization error.
        }

        // A custom root name, unsupported schema, or fast-path parse failure
        // belongs on the original reader.  Never enter it unless rewind is
        // known to have succeeded.
        rewindXmlStreamOrThrow(stream, start);
    }


    SerializationNode root;
    StreamReader reader(stream);
    IrrXMLReader* xml = createIrrXMLReader(&reader);
    
    // Find the root node in the file.
    
    while (xml->read() && xml->getNodeType() != EXN_ELEMENT)
        ;
    decodeNode(root, *xml);
    delete xml;
    
    // Process the SerializationNodes.
    
    const SerializationProxy& proxy = SerializationProxy::getProxy(root.getStringProperty("type"));
    return proxy.deserialize(root);
}
