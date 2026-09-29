/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * Streaming XML fast paths for very large canonical State/System documents. *
 * -------------------------------------------------------------------------- */

#include "FastStateXml.h"
#include "openmm/OpenMMException.h"
#include "openmm/CMMotionRemover.h"
#include "openmm/Force.h"
#include "openmm/HarmonicAngleForce.h"
#include "openmm/HarmonicBondForce.h"
#include "openmm/MonteCarloBarostat.h"
#include "openmm/NonbondedForce.h"
#include "openmm/PeriodicTorsionForce.h"
#include "openmm/State.h"
#include "openmm/System.h"
#include "openmm/Vec3.h"
#include "openmm/VirtualSite.h"
#include <cstdint>
#include <cstring>
#include <istream>
#include <limits>
#include <map>
#include <memory>
#include <sstream>
#include <streambuf>
#include <string>
#include <utility>
#include <vector>

using namespace OpenMM;
using namespace std;

extern "C" double strtod2(const char* s00, char** se);

namespace {

class FastXmlInput {
public:
    explicit FastXmlInput(istream& stream) : source(stream.rdbuf()), buffer(1024*1024), pos(0), end(0) {
        if (source == NULL)
            throw OpenMMException("XmlSerializer: XML stream has no stream buffer");
        try {
            if (source->pubseekoff(0, ios_base::cur, ios_base::in) == streampos(streamoff(-1)))
                throw OpenMMException("XmlSerializer: fast XML path requires a seekable stream");
        }
        catch (const ios_base::failure&) {
            throw OpenMMException("XmlSerializer: fast XML path requires a seekable stream");
        }
    }

    bool get(char& value) {
        if (pos == end && !refill())
            return false;
        value = buffer[pos++];
        return true;
    }

    void finish() {
        if (pos == end)
            return;
        const streamoff unread = static_cast<streamoff>(end-pos);
        streampos restored;
        try {
            restored = source->pubseekoff(-unread, ios_base::cur, ios_base::in);
        }
        catch (const ios_base::failure&) {
            throw OpenMMException("XmlSerializer: failed to restore XML stream position");
        }
        if (restored == streampos(streamoff(-1)))
            throw OpenMMException("XmlSerializer: failed to restore XML stream position");
        pos = end = 0;
    }

private:
    streambuf* source;
    vector<char> buffer;
    size_t pos, end;

    bool refill() {
        streamsize got;
        try {
            got = source->sgetn(&buffer[0], static_cast<streamsize>(buffer.size()));
        }
        catch (const ios_base::failure&) {
            throw OpenMMException("XmlSerializer: error reading XML stream");
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
                        throw OpenMMException("XmlSerializer: unsupported XML encoding");
                    continue;
                }
            }

            if (isSpace(c))
                continue;
            if (c != '<')
                throw OpenMMException("XmlSerializer: unexpected text in XML");

            char nextChar;
            if (!input.get(nextChar))
                throw OpenMMException("XmlSerializer: truncated XML markup");

            if (nextChar == '?') {
                skipSequence("?>");
                continue;
            }
            if (nextChar == '!') {
                char a, b;
                if (!input.get(a) || !input.get(b))
                    throw OpenMMException("XmlSerializer: truncated XML declaration");
                if (a == '-' && b == '-') {
                    skipSequence("-->");
                    continue;
                }
                throw OpenMMException("XmlSerializer: unsupported XML declaration");
            }

            tag.clear();
            tag.push_back('<');
            tag.push_back(nextChar);
            char quote = 0;
            while (input.get(c)) {
                tag.push_back(c);
                if (tag.size() > MAX_TAG_BYTES)
                    throw OpenMMException("XmlSerializer: XML tag is unreasonably large");
                if (quote != 0) {
                    if (c == quote)
                        quote = 0;
                }
                else if (c == '\'' || c == '"')
                    quote = c;
                else if (c == '>')
                    return tag;
            }
            throw OpenMMException("XmlSerializer: truncated XML tag");
        }
        throw OpenMMException("XmlSerializer: unexpected end of XML");
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
        throw OpenMMException("XmlSerializer: truncated XML declaration");
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
        throw OpenMMException("XmlSerializer: malformed XML tag");

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
        throw OpenMMException("XmlSerializer: malformed XML element name");

    size_t tail = text.size()-1;
    while (tail > pos && fastXmlSpace(text[tail-1]))
        --tail;
    bool empty = false;
    if (!closing && tail > pos && text[tail-1] == '/')
        empty = true;
    if (closing) {
        for (size_t i = pos; i < tail; ++i) {
            if (!fastXmlSpace(text[i]))
                throw OpenMMException("XmlSerializer: malformed XML closing tag");
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
            throw OpenMMException("XmlSerializer: malformed XML attribute name");
        name.data = text.data()+nameStart;
        name.size = pos-nameStart;

        while (pos < limit && fastXmlSpace(text[pos]))
            ++pos;
        if (pos >= limit || text[pos] != '=')
            throw OpenMMException("XmlSerializer: malformed XML attribute");
        ++pos;
        while (pos < limit && fastXmlSpace(text[pos]))
            ++pos;
        if (pos >= limit || (text[pos] != '"' && text[pos] != '\''))
            throw OpenMMException("XmlSerializer: malformed XML attribute value");
        const char quote = text[pos++];
        const size_t valueStart = pos;
        while (pos < limit && text[pos] != quote)
            ++pos;
        if (pos >= limit)
            throw OpenMMException("XmlSerializer: unterminated XML attribute value");
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
        throw OpenMMException(string("XmlSerializer: invalid floating point XML attribute '")+attributeName+"'");
    return value;
}

static long long fastXmlLongSpan(const FastXmlSpan& span, const char* attributeName) {
    if (span.size == 0)
        throw OpenMMException(string("XmlSerializer: invalid integer XML attribute '")+attributeName+"'");
    size_t pos = 0;
    bool negative = false;
    if (span.data[pos] == '-' || span.data[pos] == '+') {
        negative = span.data[pos] == '-';
        if (++pos == span.size)
            throw OpenMMException(string("XmlSerializer: invalid integer XML attribute '")+attributeName+"'");
    }
    uint64_t value = 0;
    const uint64_t positiveLimit = static_cast<uint64_t>(numeric_limits<long long>::max());
    const uint64_t limit = negative ? positiveLimit+UINT64_C(1) : positiveLimit;
    for (; pos < span.size; ++pos) {
        const char c = span.data[pos];
        if (c < '0' || c > '9')
            throw OpenMMException(string("XmlSerializer: invalid integer XML attribute '")+attributeName+"'");
        const unsigned int digit = static_cast<unsigned int>(c-'0');
        if (value > (limit-digit)/10)
            throw OpenMMException(string("XmlSerializer: integer XML attribute out of range '")+attributeName+"'");
        value = value*10 + digit;
    }
    if (!negative)
        return static_cast<long long>(value);
    if (value == positiveLimit+UINT64_C(1))
        return numeric_limits<long long>::min();
    return -static_cast<long long>(value);
}

static int fastXmlIntSpan(const FastXmlSpan& span, const char* attributeName) {
    const long long value = fastXmlLongSpan(span, attributeName);
    if (value < numeric_limits<int>::min() || value > numeric_limits<int>::max())
        throw OpenMMException(string("XmlSerializer: integer XML attribute out of range '")+attributeName+"'");
    return static_cast<int>(value);
}

static string fastXmlPlainName(const FastXmlSpan& span) {
    // The canonical OpenMM State writer uses parameter names as XML attribute
    // names.  If escaping or unusual syntax is present, let the stock parser
    // handle it rather than implementing the full XML Names grammar here.
    for (size_t i = 0; i < span.size; ++i) {
        if (!fastXmlNameChar(span.data[i]))
            throw OpenMMException("XmlSerializer: unsupported escaped State parameter name");
    }
    return string(span.data, span.size);
}

static bool findFastXmlAttribute(const FastXmlTag& tag, const char* expected,
                                 FastXmlSpan& value) {
    FastXmlAttributes attrs(tag);
    FastXmlSpan name, candidate;
    bool found = false;
    while (attrs.next(name, candidate)) {
        if (fastXmlSpanEquals(name, expected)) {
            if (found)
                throw OpenMMException(string("XmlSerializer: duplicate XML attribute '")+expected+"'");
            value = candidate;
            found = true;
        }
    }
    return found;
}

static void markFastXmlAttributeSeen(bool& seen, const char* name) {
    if (seen)
        throw OpenMMException(string("XmlSerializer: duplicate XML attribute '")+name+"'");
    seen = true;
}

static FastXmlSpan requireFastXmlAttribute(const FastXmlTag& tag, const char* name) {
    FastXmlSpan value;
    if (!findFastXmlAttribute(tag, name, value))
        throw OpenMMException(string("XmlSerializer: missing XML attribute '")+name+"'");
    return value;
}

static string fastXmlPlainString(const FastXmlSpan& span, const char* attributeName) {
    // The fast path intentionally handles only the canonical strings produced by
    // OpenMM that need no entity decoding.  Escaped strings fall back to the
    // stock XML parser rather than risking a semantic mismatch.
    for (size_t i = 0; i < span.size; ++i) {
        if (span.data[i] == '&')
            throw OpenMMException(string("XmlSerializer: escaped XML string requires stock parser for attribute '")+attributeName+"'");
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
    return fastXmlIntSpan(requireFastXmlAttribute(tag, name), name);
}

static int fastXmlInt(const FastXmlTag& tag, const char* name, int defaultValue) {
    FastXmlSpan value;
    if (!findFastXmlAttribute(tag, name, value))
        return defaultValue;
    return fastXmlIntSpan(value, name);
}

static bool fastXmlBoolSpan(const FastXmlSpan& value, const char* name) {
    // SerializationNode::setBoolProperty() writes canonical OpenMM booleans as
    // numeric 0/1.  Keep the fast path deliberately no broader than that
    // representation: non-canonical spellings fall back to the stock parser.
    if (fastXmlSpanEquals(value, "1"))
        return true;
    if (fastXmlSpanEquals(value, "0"))
        return false;
    throw OpenMMException(string("XmlSerializer: invalid boolean XML attribute '")+name+"'");
}

static bool fastXmlBool(const FastXmlTag& tag, const char* name) {
    return fastXmlBoolSpan(requireFastXmlAttribute(tag, name), name);
}

static bool fastXmlBool(const FastXmlTag& tag, const char* name, bool defaultValue) {
    FastXmlSpan value;
    if (!findFastXmlAttribute(tag, name, value))
        return defaultValue;
    return fastXmlBoolSpan(value, name);
}

static void consumeLeafClose(FastXmlTagReader& reader, const FastXmlTag& opening,
                             const char* elementName) {
    if (opening.empty)
        return;
    FastXmlTag close = parseFastXmlTag(reader.next());
    if (!close.closing || !fastXmlNameEquals(close, elementName))
        throw OpenMMException(string("XmlSerializer: unexpected child in XML element '")+elementName+"'");
}

static Vec3 parseStateVec3Tag(const FastXmlTag& tag, const char* expectedName) {
    if (tag.closing || !tag.empty || !fastXmlNameEquals(tag, expectedName))
        throw OpenMMException(string("XmlSerializer: unexpected element in State XML array; expected ")+expectedName);

    bool haveX = false, haveY = false, haveZ = false;
    double x = 0, y = 0, z = 0;
    FastXmlAttributes attrs(tag);
    FastXmlSpan name, value;
    while (attrs.next(name, value)) {
        if (fastXmlSpanEquals(name, "x")) {
            markFastXmlAttributeSeen(haveX, "x");
            x = fastXmlDoubleSpan(value, "x");
        }
        else if (fastXmlSpanEquals(name, "y")) {
            markFastXmlAttributeSeen(haveY, "y");
            y = fastXmlDoubleSpan(value, "y");
        }
        else if (fastXmlSpanEquals(name, "z")) {
            markFastXmlAttributeSeen(haveZ, "z");
            z = fastXmlDoubleSpan(value, "z");
        }
    }
    if (!haveX || !haveY || !haveZ)
        throw OpenMMException("XmlSerializer: incomplete Vec3 in State XML");
    return Vec3(x, y, z);
}

static void parseStateBox(FastXmlTagReader& reader, const FastXmlTag& opening,
                          Vec3& a, Vec3& b, Vec3& c) {
    if (opening.closing || opening.empty || !fastXmlNameEquals(opening, "PeriodicBoxVectors"))
        throw OpenMMException("XmlSerializer: malformed PeriodicBoxVectors section");

    FastXmlTag tagA = parseFastXmlTag(reader.next());
    a = parseStateVec3Tag(tagA, "A");
    FastXmlTag tagB = parseFastXmlTag(reader.next());
    b = parseStateVec3Tag(tagB, "B");
    FastXmlTag tagC = parseFastXmlTag(reader.next());
    c = parseStateVec3Tag(tagC, "C");
    FastXmlTag close = parseFastXmlTag(reader.next());
    if (!close.closing || !fastXmlNameEquals(close, "PeriodicBoxVectors"))
        throw OpenMMException("XmlSerializer: malformed PeriodicBoxVectors closing tag");
}

static map<string, double> parseStateParameters(const FastXmlTag& tag) {
    if (tag.closing || !tag.empty || !fastXmlNameEquals(tag, "Parameters"))
        throw OpenMMException("XmlSerializer: malformed Parameters section");
    map<string, double> result;
    FastXmlAttributes attrs(tag);
    FastXmlSpan name, value;
    while (attrs.next(name, value)) {
        string parameter = fastXmlPlainName(name);
        if (result.find(parameter) != result.end())
            throw OpenMMException(string("XmlSerializer: duplicate XML attribute '")+parameter+"'");
        result[parameter] = fastXmlDoubleSpan(value, parameter.c_str());
    }
    return result;
}

static void parseStateEnergies(const FastXmlTag& tag, double& kinetic, double& potential) {
    if (tag.closing || !tag.empty || !fastXmlNameEquals(tag, "Energies"))
        throw OpenMMException("XmlSerializer: malformed Energies section");
    bool haveKinetic = false, havePotential = false;
    FastXmlAttributes attrs(tag);
    FastXmlSpan name, value;
    while (attrs.next(name, value)) {
        if (fastXmlSpanEquals(name, "KineticEnergy")) {
            markFastXmlAttributeSeen(haveKinetic, "KineticEnergy");
            kinetic = fastXmlDoubleSpan(value, "KineticEnergy");
        }
        else if (fastXmlSpanEquals(name, "PotentialEnergy")) {
            markFastXmlAttributeSeen(havePotential, "PotentialEnergy");
            potential = fastXmlDoubleSpan(value, "PotentialEnergy");
        }
    }
    if (!haveKinetic || !havePotential)
        throw OpenMMException("XmlSerializer: incomplete Energies section");
}

static void parseStateVec3Section(FastXmlTagReader& reader, const FastXmlTag& opening,
                                  const char* sectionName, const char* elementName,
                                  vector<Vec3>& values, size_t reserveHint) {
    if (opening.closing || !fastXmlNameEquals(opening, sectionName))
        throw OpenMMException(string("XmlSerializer: malformed ")+sectionName+" section");
    if (reserveHint != 0)
        values.reserve(reserveHint);
    if (opening.empty)
        return;

    while (true) {
        FastXmlTag tag = parseFastXmlTag(reader.next());
        if (tag.closing) {
            if (!fastXmlNameEquals(tag, sectionName))
                throw OpenMMException(string("XmlSerializer: unexpected closing tag in ")+sectionName);
            return;
        }
        values.push_back(parseStateVec3Tag(tag, elementName));
    }
}

static State* parseStateXmlStreaming(istream& stream) {
    FastXmlTagReader reader(stream);
    FastXmlTag root = parseFastXmlTag(reader.next());
    if (root.closing || root.empty || !fastXmlNameEquals(root, "State"))
        throw OpenMMException("XmlSerializer: XML object is not a State");

    bool haveType = false, haveVersion = false, haveTime = false, haveStepCount = false;
    string type;
    long long stepCount = 0;
    double time = 0;
    FastXmlAttributes rootAttrs(root);
    FastXmlSpan attrName, attrValue;
    while (rootAttrs.next(attrName, attrValue)) {
        if (fastXmlSpanEquals(attrName, "type")) {
            markFastXmlAttributeSeen(haveType, "type");
            type.assign(attrValue.data, attrValue.size);
        }
        else if (fastXmlSpanEquals(attrName, "version")) {
            markFastXmlAttributeSeen(haveVersion, "version");
            if (!fastXmlSpanEquals(attrValue, "1"))
                throw OpenMMException("XmlSerializer: unsupported State XML version");
        }
        else if (fastXmlSpanEquals(attrName, "time")) {
            markFastXmlAttributeSeen(haveTime, "time");
            time = fastXmlDoubleSpan(attrValue, "time");
        }
        else if (fastXmlSpanEquals(attrName, "stepCount")) {
            markFastXmlAttributeSeen(haveStepCount, "stepCount");
            stepCount = fastXmlLongSpan(attrValue, "stepCount");
        }
    }
    if (!haveType || type != "State")
        throw OpenMMException("XmlSerializer: XML type is not State");
    if (!haveVersion || !haveTime)
        throw OpenMMException("XmlSerializer: incomplete State root attributes");

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
                throw OpenMMException("XmlSerializer: unexpected closing tag in State XML");
            break;
        }

        if (fastXmlNameEquals(tag, "PeriodicBoxVectors")) {
            if (stage != 0)
                throw OpenMMException("XmlSerializer: unexpected State XML ordering");
            parseStateBox(reader, tag, a, b, c);
            haveBox = true;
            stage = 1;
        }
        else if (fastXmlNameEquals(tag, "Parameters")) {
            if (!haveBox || stage >= 2)
                throw OpenMMException("XmlSerializer: unexpected Parameters section");
            parameters = parseStateParameters(tag);
            types |= State::Parameters;
            stage = 2;
        }
        else if (fastXmlNameEquals(tag, "Energies")) {
            if (!haveBox || stage >= 3)
                throw OpenMMException("XmlSerializer: unexpected Energies section");
            parseStateEnergies(tag, kinetic, potential);
            types |= State::Energy;
            stage = 3;
        }
        else if (fastXmlNameEquals(tag, "Positions")) {
            if (!haveBox || stage >= 4)
                throw OpenMMException("XmlSerializer: unexpected Positions section");
            parseStateVec3Section(reader, tag, "Positions", "Position", positions, 0);
            types |= State::Positions;
            stage = 4;
        }
        else if (fastXmlNameEquals(tag, "Velocities")) {
            if (!haveBox || stage >= 5)
                throw OpenMMException("XmlSerializer: unexpected Velocities section");
            parseStateVec3Section(reader, tag, "Velocities", "Velocity", velocities,
                                  (types & State::Positions) ? positions.size() : 0);
            if ((types & State::Positions) && velocities.size() != positions.size())
                throw OpenMMException("XmlSerializer: State particle array size mismatch");
            types |= State::Velocities;
            stage = 5;
        }
        else if (fastXmlNameEquals(tag, "Forces")) {
            if (!haveBox || stage >= 6)
                throw OpenMMException("XmlSerializer: unexpected Forces section");
            const size_t hint = (types & State::Positions) ? positions.size() :
                                ((types & State::Velocities) ? velocities.size() : 0);
            parseStateVec3Section(reader, tag, "Forces", "Force", forces, hint);
            if (((types & State::Positions) && forces.size() != positions.size()) ||
                    (!(types & State::Positions) && (types & State::Velocities) &&
                     forces.size() != velocities.size()))
                throw OpenMMException("XmlSerializer: State particle array size mismatch");
            types |= State::Forces;
            stage = 6;
        }
        else if (fastXmlNameEquals(tag, "IntegratorParameters"))
            throw OpenMMException("XmlSerializer: IntegratorParameters are not supported by the fast State path");
        else
            throw OpenMMException("XmlSerializer: unsupported State XML element");
    }

    if (!haveBox)
        throw OpenMMException("XmlSerializer: State has no PeriodicBoxVectors");

    // If Positions are absent, still enforce consistency between any other two
    // particle arrays that are present.
    if (!(types & State::Positions) && (types & State::Velocities) &&
            (types & State::Forces) && velocities.size() != forces.size())
        throw OpenMMException("XmlSerializer: State particle array size mismatch");

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
        throw OpenMMException("XmlSerializer: malformed System PeriodicBoxVectors section");

    FastXmlTag aTag = parseFastXmlTag(reader.next());
    Vec3 a = parseStateVec3Tag(aTag, "A");
    FastXmlTag bTag = parseFastXmlTag(reader.next());
    Vec3 b = parseStateVec3Tag(bTag, "B");
    FastXmlTag cTag = parseFastXmlTag(reader.next());
    Vec3 c = parseStateVec3Tag(cTag, "C");
    FastXmlTag close = parseFastXmlTag(reader.next());
    if (!close.closing || !fastXmlNameEquals(close, "PeriodicBoxVectors"))
        throw OpenMMException("XmlSerializer: malformed System PeriodicBoxVectors closing tag");
    system.setDefaultPeriodicBoxVectors(a, b, c);
}

static void parseSystemVirtualSiteStreaming(FastXmlTagReader& reader, const FastXmlTag& tag,
                                            System& system, int particle) {
    if (tag.closing)
        throw OpenMMException("XmlSerializer: malformed virtual site element");

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
                throw OpenMMException("XmlSerializer: virtual-site particle index out of range");
            particles.push_back(static_cast<int>(p));
            wo.push_back(fastXmlDouble(tag, ("wo"+index).c_str()));
            wx.push_back(fastXmlDouble(tag, ("wx"+index).c_str()));
            wy.push_back(fastXmlDouble(tag, ("wy"+index).c_str()));
        }
        if (particles.empty())
            throw OpenMMException("XmlSerializer: LocalCoordinatesSite has no particles");
        Vec3 position(fastXmlDouble(tag, "pos1"), fastXmlDouble(tag, "pos2"),
                      fastXmlDouble(tag, "pos3"));
        system.setVirtualSite(particle, new LocalCoordinatesSite(particles, wo, wx, wy, position));
        consumeLeafClose(reader, tag, "LocalCoordinatesSite");
    }
    else {
        throw OpenMMException("XmlSerializer: unsupported virtual site type in fast System path");
    }
}

static void parseSystemParticlesStreaming(FastXmlTagReader& reader, const FastXmlTag& opening,
                                          System& system) {
    if (opening.closing || !fastXmlNameEquals(opening, "Particles"))
        throw OpenMMException("XmlSerializer: malformed System Particles section");
    if (opening.empty)
        return;

    while (true) {
        FastXmlTag tag = parseFastXmlTag(reader.next());
        if (tag.closing) {
            if (!fastXmlNameEquals(tag, "Particles"))
                throw OpenMMException("XmlSerializer: unexpected closing tag in System Particles");
            return;
        }
        if (!fastXmlNameEquals(tag, "Particle"))
            throw OpenMMException("XmlSerializer: unexpected element in System Particles");
        if (system.getNumParticles() >= numeric_limits<int>::max())
            throw OpenMMException("XmlSerializer: too many System particles");

        const int particle = system.getNumParticles();
        system.addParticle(fastXmlDouble(tag, "mass"));
        if (tag.empty)
            continue;

        FastXmlTag child = parseFastXmlTag(reader.next());
        if (child.closing) {
            if (!fastXmlNameEquals(child, "Particle"))
                throw OpenMMException("XmlSerializer: malformed Particle closing tag");
            continue;
        }
        parseSystemVirtualSiteStreaming(reader, child, system, particle);
        FastXmlTag close = parseFastXmlTag(reader.next());
        if (!close.closing || !fastXmlNameEquals(close, "Particle"))
            throw OpenMMException("XmlSerializer: Particle has multiple or malformed virtual sites");
    }
}

static void parseSystemConstraintsStreaming(FastXmlTagReader& reader, const FastXmlTag& opening,
                                            System& system) {
    if (opening.closing || !fastXmlNameEquals(opening, "Constraints"))
        throw OpenMMException("XmlSerializer: malformed System Constraints section");
    if (opening.empty)
        return;

    while (true) {
        FastXmlTag tag = parseFastXmlTag(reader.next());
        if (tag.closing) {
            if (!fastXmlNameEquals(tag, "Constraints"))
                throw OpenMMException("XmlSerializer: unexpected closing tag in Constraints");
            return;
        }
        if (!fastXmlNameEquals(tag, "Constraint"))
            throw OpenMMException("XmlSerializer: unexpected element in System Constraints");

        bool haveP1 = false, haveP2 = false, haveDistance = false;
        int p1 = 0, p2 = 0;
        double distance = 0;
        FastXmlAttributes attrs(tag);
        FastXmlSpan name, value;
        while (attrs.next(name, value)) {
            if (fastXmlSpanEquals(name, "p1")) {
                markFastXmlAttributeSeen(haveP1, "p1");
                p1 = fastXmlIntSpan(value, "p1");
            }
            else if (fastXmlSpanEquals(name, "p2")) {
                markFastXmlAttributeSeen(haveP2, "p2");
                p2 = fastXmlIntSpan(value, "p2");
            }
            else if (fastXmlSpanEquals(name, "d")) {
                markFastXmlAttributeSeen(haveDistance, "d");
                distance = fastXmlDoubleSpan(value, "d");
            }
        }
        if (!haveP1 || !haveP2 || !haveDistance)
            throw OpenMMException("XmlSerializer: incomplete Constraint element");
        system.addConstraint(p1, p2, distance);
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
        throw OpenMMException("XmlSerializer: malformed NonbondedForce section");
    if (opening.empty)
        return;

    while (true) {
        FastXmlTag tag = parseFastXmlTag(reader.next());
        if (tag.closing) {
            if (!fastXmlNameEquals(tag, sectionName))
                throw OpenMMException("XmlSerializer: unexpected closing tag in NonbondedForce section");
            return;
        }
        if (!fastXmlNameEquals(tag, elementName))
            throw OpenMMException("XmlSerializer: unexpected element in NonbondedForce section");

        switch (kind) {
            case FAST_NB_GLOBAL_PARAMETERS: {
                bool haveName = false, haveDefault = false;
                string parameterName;
                double defaultValue = 0;
                FastXmlAttributes attrs(tag);
                FastXmlSpan name, value;
                while (attrs.next(name, value)) {
                    if (fastXmlSpanEquals(name, "name")) {
                        markFastXmlAttributeSeen(haveName, "name");
                        parameterName = fastXmlPlainString(value, "name");
                    }
                    else if (fastXmlSpanEquals(name, "default")) {
                        markFastXmlAttributeSeen(haveDefault, "default");
                        defaultValue = fastXmlDoubleSpan(value, "default");
                    }
                }
                if (!haveName || !haveDefault)
                    throw OpenMMException("XmlSerializer: incomplete NonbondedForce Parameter element");
                force.addGlobalParameter(parameterName, defaultValue);
                break;
            }
            case FAST_NB_PARTICLE_OFFSETS: {
                bool haveParameter = false, haveParticle = false;
                bool haveQ = false, haveSig = false, haveEps = false;
                string parameter;
                int particle = 0;
                double q = 0, sig = 0, eps = 0;
                FastXmlAttributes attrs(tag);
                FastXmlSpan name, value;
                while (attrs.next(name, value)) {
                    if (fastXmlSpanEquals(name, "parameter")) { markFastXmlAttributeSeen(haveParameter, "parameter"); parameter = fastXmlPlainString(value, "parameter"); }
                    else if (fastXmlSpanEquals(name, "particle")) { markFastXmlAttributeSeen(haveParticle, "particle"); particle = fastXmlIntSpan(value, "particle"); }
                    else if (fastXmlSpanEquals(name, "q")) { markFastXmlAttributeSeen(haveQ, "q"); q = fastXmlDoubleSpan(value, "q"); }
                    else if (fastXmlSpanEquals(name, "sig")) { markFastXmlAttributeSeen(haveSig, "sig"); sig = fastXmlDoubleSpan(value, "sig"); }
                    else if (fastXmlSpanEquals(name, "eps")) { markFastXmlAttributeSeen(haveEps, "eps"); eps = fastXmlDoubleSpan(value, "eps"); }
                }
                if (!haveParameter || !haveParticle || !haveQ || !haveSig || !haveEps)
                    throw OpenMMException("XmlSerializer: incomplete NonbondedForce particle Offset element");
                force.addParticleParameterOffset(parameter, particle, q, sig, eps);
                break;
            }
            case FAST_NB_EXCEPTION_OFFSETS: {
                bool haveParameter = false, haveException = false;
                bool haveQ = false, haveSig = false, haveEps = false;
                string parameter;
                int exception = 0;
                double q = 0, sig = 0, eps = 0;
                FastXmlAttributes attrs(tag);
                FastXmlSpan name, value;
                while (attrs.next(name, value)) {
                    if (fastXmlSpanEquals(name, "parameter")) { markFastXmlAttributeSeen(haveParameter, "parameter"); parameter = fastXmlPlainString(value, "parameter"); }
                    else if (fastXmlSpanEquals(name, "exception")) { markFastXmlAttributeSeen(haveException, "exception"); exception = fastXmlIntSpan(value, "exception"); }
                    else if (fastXmlSpanEquals(name, "q")) { markFastXmlAttributeSeen(haveQ, "q"); q = fastXmlDoubleSpan(value, "q"); }
                    else if (fastXmlSpanEquals(name, "sig")) { markFastXmlAttributeSeen(haveSig, "sig"); sig = fastXmlDoubleSpan(value, "sig"); }
                    else if (fastXmlSpanEquals(name, "eps")) { markFastXmlAttributeSeen(haveEps, "eps"); eps = fastXmlDoubleSpan(value, "eps"); }
                }
                if (!haveParameter || !haveException || !haveQ || !haveSig || !haveEps)
                    throw OpenMMException("XmlSerializer: incomplete NonbondedForce exception Offset element");
                force.addExceptionParameterOffset(parameter, exception, q, sig, eps);
                break;
            }
            case FAST_NB_PARTICLES: {
                bool haveQ = false, haveSig = false, haveEps = false;
                double q = 0, sig = 0, eps = 0;
                FastXmlAttributes attrs(tag);
                FastXmlSpan name, value;
                while (attrs.next(name, value)) {
                    if (fastXmlSpanEquals(name, "q")) { markFastXmlAttributeSeen(haveQ, "q"); q = fastXmlDoubleSpan(value, "q"); }
                    else if (fastXmlSpanEquals(name, "sig")) { markFastXmlAttributeSeen(haveSig, "sig"); sig = fastXmlDoubleSpan(value, "sig"); }
                    else if (fastXmlSpanEquals(name, "eps")) { markFastXmlAttributeSeen(haveEps, "eps"); eps = fastXmlDoubleSpan(value, "eps"); }
                }
                if (!haveQ || !haveSig || !haveEps)
                    throw OpenMMException("XmlSerializer: incomplete NonbondedForce Particle element");
                force.addParticle(q, sig, eps);
                break;
            }
            case FAST_NB_EXCEPTIONS: {
                bool haveP1 = false, haveP2 = false, haveQ = false, haveSig = false, haveEps = false;
                int p1 = 0, p2 = 0;
                double q = 0, sig = 0, eps = 0;
                FastXmlAttributes attrs(tag);
                FastXmlSpan name, value;
                while (attrs.next(name, value)) {
                    if (fastXmlSpanEquals(name, "p1")) { markFastXmlAttributeSeen(haveP1, "p1"); p1 = fastXmlIntSpan(value, "p1"); }
                    else if (fastXmlSpanEquals(name, "p2")) { markFastXmlAttributeSeen(haveP2, "p2"); p2 = fastXmlIntSpan(value, "p2"); }
                    else if (fastXmlSpanEquals(name, "q")) { markFastXmlAttributeSeen(haveQ, "q"); q = fastXmlDoubleSpan(value, "q"); }
                    else if (fastXmlSpanEquals(name, "sig")) { markFastXmlAttributeSeen(haveSig, "sig"); sig = fastXmlDoubleSpan(value, "sig"); }
                    else if (fastXmlSpanEquals(name, "eps")) { markFastXmlAttributeSeen(haveEps, "eps"); eps = fastXmlDoubleSpan(value, "eps"); }
                }
                if (!haveP1 || !haveP2 || !haveQ || !haveSig || !haveEps)
                    throw OpenMMException("XmlSerializer: incomplete NonbondedForce Exception element");
                force.addException(p1, p2, q, sig, eps);
                break;
            }
        }
        consumeLeafClose(reader, tag, elementName);
    }
}

static Force* parseNonbondedForceStreaming(FastXmlTagReader& reader, const FastXmlTag& opening) {
    const int version = fastXmlInt(opening, "version");
    if (version < 1 || version > 4)
        throw OpenMMException("XmlSerializer: unsupported NonbondedForce XML version");

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
        throw OpenMMException("XmlSerializer: incomplete NonbondedForce");

    while (true) {
        FastXmlTag tag = parseFastXmlTag(reader.next());
        if (tag.closing) {
            if (!fastXmlNameEquals(tag, "Force"))
                throw OpenMMException("XmlSerializer: unexpected closing tag in NonbondedForce");
            if (!haveParticles || !haveExceptions)
                throw OpenMMException("XmlSerializer: incomplete NonbondedForce");
            if (version >= 3 && (!haveGlobal || !haveParticleOffsets || !haveExceptionOffsets))
                throw OpenMMException("XmlSerializer: incomplete NonbondedForce parameter sections");
            return force.release();
        }
        if (fastXmlNameEquals(tag, "GlobalParameters")) {
            if (version < 3 || haveGlobal)
                throw OpenMMException("XmlSerializer: unexpected GlobalParameters");
            parseNonbondedSectionStreaming(reader, tag, FAST_NB_GLOBAL_PARAMETERS, *force);
            haveGlobal = true;
        }
        else if (fastXmlNameEquals(tag, "ParticleOffsets")) {
            if (version < 3 || haveParticleOffsets)
                throw OpenMMException("XmlSerializer: unexpected ParticleOffsets");
            parseNonbondedSectionStreaming(reader, tag, FAST_NB_PARTICLE_OFFSETS, *force);
            haveParticleOffsets = true;
        }
        else if (fastXmlNameEquals(tag, "ExceptionOffsets")) {
            if (version < 3 || haveExceptionOffsets)
                throw OpenMMException("XmlSerializer: unexpected ExceptionOffsets");
            parseNonbondedSectionStreaming(reader, tag, FAST_NB_EXCEPTION_OFFSETS, *force);
            haveExceptionOffsets = true;
        }
        else if (fastXmlNameEquals(tag, "Particles")) {
            if (haveParticles)
                throw OpenMMException("XmlSerializer: duplicate NonbondedForce Particles");
            parseNonbondedSectionStreaming(reader, tag, FAST_NB_PARTICLES, *force);
            haveParticles = true;
        }
        else if (fastXmlNameEquals(tag, "Exceptions")) {
            if (haveExceptions)
                throw OpenMMException("XmlSerializer: duplicate NonbondedForce Exceptions");
            parseNonbondedSectionStreaming(reader, tag, FAST_NB_EXCEPTIONS, *force);
            haveExceptions = true;
        }
        else
            throw OpenMMException("XmlSerializer: unsupported NonbondedForce child in fast System path");
    }
}

static Force* parseHarmonicBondForceStreaming(FastXmlTagReader& reader, const FastXmlTag& opening) {
    const int version = fastXmlInt(opening, "version");
    if (version < 1 || version > 2)
        throw OpenMMException("XmlSerializer: unsupported HarmonicBondForce XML version");
    unique_ptr<HarmonicBondForce> force(new HarmonicBondForce());
    force->setForceGroup(fastXmlInt(opening, "forceGroup", 0));
    force->setName(fastXmlString(opening, "name", force->getName()));
    if (version > 1)
        force->setUsesPeriodicBoundaryConditions(fastXmlBool(opening, "usesPeriodic"));
    bool haveBonds = false;
    if (opening.empty)
        throw OpenMMException("XmlSerializer: HarmonicBondForce has no Bonds section");
    while (true) {
        FastXmlTag tag = parseFastXmlTag(reader.next());
        if (tag.closing) {
            if (!fastXmlNameEquals(tag, "Force"))
                throw OpenMMException("XmlSerializer: unexpected closing tag in HarmonicBondForce");
            if (!haveBonds)
                throw OpenMMException("XmlSerializer: HarmonicBondForce has no Bonds section");
            return force.release();
        }
        if (!fastXmlNameEquals(tag, "Bonds") || haveBonds)
            throw OpenMMException("XmlSerializer: unexpected HarmonicBondForce element");
        haveBonds = true;
        if (tag.empty)
            continue;
        while (true) {
            FastXmlTag bond = parseFastXmlTag(reader.next());
            if (bond.closing) {
                if (!fastXmlNameEquals(bond, "Bonds"))
                    throw OpenMMException("XmlSerializer: unexpected closing tag in Bonds");
                break;
            }
            if (!fastXmlNameEquals(bond, "Bond"))
                throw OpenMMException("XmlSerializer: unexpected element in HarmonicBondForce Bonds");
            bool haveP1 = false, haveP2 = false, haveD = false, haveK = false;
            int p1 = 0, p2 = 0;
            double d = 0, k = 0;
            FastXmlAttributes attrs(bond);
            FastXmlSpan name, value;
            while (attrs.next(name, value)) {
                if (fastXmlSpanEquals(name, "p1")) { markFastXmlAttributeSeen(haveP1, "p1"); p1 = fastXmlIntSpan(value, "p1"); }
                else if (fastXmlSpanEquals(name, "p2")) { markFastXmlAttributeSeen(haveP2, "p2"); p2 = fastXmlIntSpan(value, "p2"); }
                else if (fastXmlSpanEquals(name, "d")) { markFastXmlAttributeSeen(haveD, "d"); d = fastXmlDoubleSpan(value, "d"); }
                else if (fastXmlSpanEquals(name, "k")) { markFastXmlAttributeSeen(haveK, "k"); k = fastXmlDoubleSpan(value, "k"); }
            }
            if (!haveP1 || !haveP2 || !haveD || !haveK)
                throw OpenMMException("XmlSerializer: incomplete HarmonicBondForce Bond element");
            force->addBond(p1, p2, d, k);
            consumeLeafClose(reader, bond, "Bond");
        }
    }
}

static Force* parseHarmonicAngleForceStreaming(FastXmlTagReader& reader, const FastXmlTag& opening) {
    const int version = fastXmlInt(opening, "version");
    if (version < 1 || version > 2)
        throw OpenMMException("XmlSerializer: unsupported HarmonicAngleForce XML version");
    unique_ptr<HarmonicAngleForce> force(new HarmonicAngleForce());
    force->setForceGroup(fastXmlInt(opening, "forceGroup", 0));
    force->setName(fastXmlString(opening, "name", force->getName()));
    if (version > 1)
        force->setUsesPeriodicBoundaryConditions(fastXmlBool(opening, "usesPeriodic"));
    bool haveAngles = false;
    if (opening.empty)
        throw OpenMMException("XmlSerializer: HarmonicAngleForce has no Angles section");
    while (true) {
        FastXmlTag tag = parseFastXmlTag(reader.next());
        if (tag.closing) {
            if (!fastXmlNameEquals(tag, "Force"))
                throw OpenMMException("XmlSerializer: unexpected closing tag in HarmonicAngleForce");
            if (!haveAngles)
                throw OpenMMException("XmlSerializer: HarmonicAngleForce has no Angles section");
            return force.release();
        }
        if (!fastXmlNameEquals(tag, "Angles") || haveAngles)
            throw OpenMMException("XmlSerializer: unexpected HarmonicAngleForce element");
        haveAngles = true;
        if (tag.empty)
            continue;
        while (true) {
            FastXmlTag angle = parseFastXmlTag(reader.next());
            if (angle.closing) {
                if (!fastXmlNameEquals(angle, "Angles"))
                    throw OpenMMException("XmlSerializer: unexpected closing tag in Angles");
                break;
            }
            if (!fastXmlNameEquals(angle, "Angle"))
                throw OpenMMException("XmlSerializer: unexpected element in HarmonicAngleForce Angles");
            bool haveP1 = false, haveP2 = false, haveP3 = false, haveA = false, haveK = false;
            int p1 = 0, p2 = 0, p3 = 0;
            double a = 0, k = 0;
            FastXmlAttributes attrs(angle);
            FastXmlSpan name, value;
            while (attrs.next(name, value)) {
                if (fastXmlSpanEquals(name, "p1")) { markFastXmlAttributeSeen(haveP1, "p1"); p1 = fastXmlIntSpan(value, "p1"); }
                else if (fastXmlSpanEquals(name, "p2")) { markFastXmlAttributeSeen(haveP2, "p2"); p2 = fastXmlIntSpan(value, "p2"); }
                else if (fastXmlSpanEquals(name, "p3")) { markFastXmlAttributeSeen(haveP3, "p3"); p3 = fastXmlIntSpan(value, "p3"); }
                else if (fastXmlSpanEquals(name, "a")) { markFastXmlAttributeSeen(haveA, "a"); a = fastXmlDoubleSpan(value, "a"); }
                else if (fastXmlSpanEquals(name, "k")) { markFastXmlAttributeSeen(haveK, "k"); k = fastXmlDoubleSpan(value, "k"); }
            }
            if (!haveP1 || !haveP2 || !haveP3 || !haveA || !haveK)
                throw OpenMMException("XmlSerializer: incomplete HarmonicAngleForce Angle element");
            force->addAngle(p1, p2, p3, a, k);
            consumeLeafClose(reader, angle, "Angle");
        }
    }
}

static Force* parsePeriodicTorsionForceStreaming(FastXmlTagReader& reader, const FastXmlTag& opening) {
    const int version = fastXmlInt(opening, "version");
    if (version < 1 || version > 2)
        throw OpenMMException("XmlSerializer: unsupported PeriodicTorsionForce XML version");
    unique_ptr<PeriodicTorsionForce> force(new PeriodicTorsionForce());
    force->setForceGroup(fastXmlInt(opening, "forceGroup", 0));
    force->setName(fastXmlString(opening, "name", force->getName()));
    if (version > 1)
        force->setUsesPeriodicBoundaryConditions(fastXmlBool(opening, "usesPeriodic"));
    bool haveTorsions = false;
    if (opening.empty)
        throw OpenMMException("XmlSerializer: PeriodicTorsionForce has no Torsions section");
    while (true) {
        FastXmlTag tag = parseFastXmlTag(reader.next());
        if (tag.closing) {
            if (!fastXmlNameEquals(tag, "Force"))
                throw OpenMMException("XmlSerializer: unexpected closing tag in PeriodicTorsionForce");
            if (!haveTorsions)
                throw OpenMMException("XmlSerializer: PeriodicTorsionForce has no Torsions section");
            return force.release();
        }
        if (!fastXmlNameEquals(tag, "Torsions") || haveTorsions)
            throw OpenMMException("XmlSerializer: unexpected PeriodicTorsionForce element");
        haveTorsions = true;
        if (tag.empty)
            continue;
        while (true) {
            FastXmlTag torsion = parseFastXmlTag(reader.next());
            if (torsion.closing) {
                if (!fastXmlNameEquals(torsion, "Torsions"))
                    throw OpenMMException("XmlSerializer: unexpected closing tag in Torsions");
                break;
            }
            if (!fastXmlNameEquals(torsion, "Torsion"))
                throw OpenMMException("XmlSerializer: unexpected element in PeriodicTorsionForce Torsions");
            bool haveP1 = false, haveP2 = false, haveP3 = false, haveP4 = false;
            bool havePeriodicity = false, havePhase = false, haveK = false;
            int p1 = 0, p2 = 0, p3 = 0, p4 = 0, periodicity = 0;
            double phase = 0, k = 0;
            FastXmlAttributes attrs(torsion);
            FastXmlSpan name, value;
            while (attrs.next(name, value)) {
                if (fastXmlSpanEquals(name, "p1")) { markFastXmlAttributeSeen(haveP1, "p1"); p1 = fastXmlIntSpan(value, "p1"); }
                else if (fastXmlSpanEquals(name, "p2")) { markFastXmlAttributeSeen(haveP2, "p2"); p2 = fastXmlIntSpan(value, "p2"); }
                else if (fastXmlSpanEquals(name, "p3")) { markFastXmlAttributeSeen(haveP3, "p3"); p3 = fastXmlIntSpan(value, "p3"); }
                else if (fastXmlSpanEquals(name, "p4")) { markFastXmlAttributeSeen(haveP4, "p4"); p4 = fastXmlIntSpan(value, "p4"); }
                else if (fastXmlSpanEquals(name, "periodicity")) { markFastXmlAttributeSeen(havePeriodicity, "periodicity"); periodicity = fastXmlIntSpan(value, "periodicity"); }
                else if (fastXmlSpanEquals(name, "phase")) { markFastXmlAttributeSeen(havePhase, "phase"); phase = fastXmlDoubleSpan(value, "phase"); }
                else if (fastXmlSpanEquals(name, "k")) { markFastXmlAttributeSeen(haveK, "k"); k = fastXmlDoubleSpan(value, "k"); }
            }
            if (!haveP1 || !haveP2 || !haveP3 || !haveP4 || !havePeriodicity || !havePhase || !haveK)
                throw OpenMMException("XmlSerializer: incomplete PeriodicTorsionForce Torsion element");
            force->addTorsion(p1, p2, p3, p4, periodicity, phase, k);
            consumeLeafClose(reader, torsion, "Torsion");
        }
    }
}

static Force* parseCMMotionRemoverStreaming(FastXmlTagReader& reader, const FastXmlTag& opening) {
    if (fastXmlInt(opening, "version") != 1)
        throw OpenMMException("XmlSerializer: unsupported CMMotionRemover XML version");
    unique_ptr<CMMotionRemover> force(new CMMotionRemover(fastXmlInt(opening, "frequency")));
    force->setForceGroup(fastXmlInt(opening, "forceGroup", 0));
    force->setName(fastXmlString(opening, "name", force->getName()));
    consumeLeafClose(reader, opening, "Force");
    return force.release();
}

static Force* parseMonteCarloBarostatStreaming(FastXmlTagReader& reader, const FastXmlTag& opening) {
    if (fastXmlInt(opening, "version") != 1)
        throw OpenMMException("XmlSerializer: unsupported MonteCarloBarostat XML version");
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
        throw OpenMMException("XmlSerializer: malformed Force element");
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
    throw OpenMMException(string("XmlSerializer: unsupported Force type '")+type+"'");
}

static void parseSystemForcesStreaming(FastXmlTagReader& reader, const FastXmlTag& opening,
                                       System& system) {
    if (opening.closing || !fastXmlNameEquals(opening, "Forces"))
        throw OpenMMException("XmlSerializer: malformed System Forces section");
    if (opening.empty)
        return;
    while (true) {
        FastXmlTag tag = parseFastXmlTag(reader.next());
        if (tag.closing) {
            if (!fastXmlNameEquals(tag, "Forces"))
                throw OpenMMException("XmlSerializer: unexpected closing tag in Forces");
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
        throw OpenMMException("XmlSerializer: XML object is not a System");
    if (fastXmlString(root, "type") != "System")
        throw OpenMMException("XmlSerializer: XML type is not System");
    if (fastXmlInt(root, "version") != 1)
        throw OpenMMException("XmlSerializer: unsupported System XML version");

    unique_ptr<System> system(new System());
    bool haveBox = false, haveParticles = false, haveConstraints = false, haveForces = false;
    int stage = 0;
    while (true) {
        FastXmlTag tag = parseFastXmlTag(reader.next());
        if (tag.closing) {
            if (!fastXmlNameEquals(tag, "System"))
                throw OpenMMException("XmlSerializer: unexpected closing tag in System XML");
            break;
        }
        if (fastXmlNameEquals(tag, "PeriodicBoxVectors")) {
            if (stage != 0 || haveBox)
                throw OpenMMException("XmlSerializer: unexpected System PeriodicBoxVectors");
            parseSystemBoxStreaming(reader, tag, *system);
            haveBox = true;
            stage = 1;
        }
        else if (fastXmlNameEquals(tag, "Particles")) {
            if (!haveBox || stage != 1 || haveParticles)
                throw OpenMMException("XmlSerializer: unexpected System Particles section");
            parseSystemParticlesStreaming(reader, tag, *system);
            haveParticles = true;
            stage = 2;
        }
        else if (fastXmlNameEquals(tag, "Constraints")) {
            if (!haveParticles || stage != 2 || haveConstraints)
                throw OpenMMException("XmlSerializer: unexpected System Constraints section");
            parseSystemConstraintsStreaming(reader, tag, *system);
            haveConstraints = true;
            stage = 3;
        }
        else if (fastXmlNameEquals(tag, "Forces")) {
            if (!haveConstraints || stage != 3 || haveForces)
                throw OpenMMException("XmlSerializer: unexpected System Forces section");
            parseSystemForcesStreaming(reader, tag, *system);
            haveForces = true;
            stage = 4;
        }
        else
            throw OpenMMException("XmlSerializer: unsupported System XML element in fast path");
    }
    if (!haveBox || !haveParticles || !haveConstraints || !haveForces)
        throw OpenMMException("XmlSerializer: incomplete System XML");
    reader.finish();
    return system.release();
}

} // anonymous namespace

namespace OpenMM {

State* deserializeStateXmlFast(istream& stream) {
    return parseStateXmlStreaming(stream);
}

System* deserializeSystemXmlFast(istream& stream) {
    return parseSystemXmlStreaming(stream);
}

} // namespace OpenMM
