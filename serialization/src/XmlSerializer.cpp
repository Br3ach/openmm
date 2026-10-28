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
#include "irrXML.h"
#include "FastStateXml.h"
#include <cstdlib>
#include <vector>
#include <cstring>
#include <iostream>
#include <map>

using namespace OpenMM;
using namespace std;
using namespace irr;
using namespace io;

/**
 * Apply XML encoding to a string.  This is adapted from TinyXML (written by Lee Thomason).
 */
static void encodeString(const string& str, string* outString) {
    static map<char, string> entities;
    static bool hasInitialized = false;
    if (!hasInitialized) {
        hasInitialized = true;
	entities['&'] = "&amp;";
	entities['<'] = "&lt;";
	entities['>'] = "&gt;";
	entities['\"'] = "&quot;";
	entities['\''] = "&apos;";
    }

    int i=0;

    while (i<(int)str.length()) {
        unsigned char c = (unsigned char) str[i];

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
        else if (entities.find(c) != entities.end()) {
            outString->append(entities[c]);
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

void XmlSerializer::serialize(const SerializationNode& node, std::ostream& stream) {
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
            if (haveRoot && rootName == "State" && type == "State")
                return deserializeStateXmlFast(stream);

            if (haveRoot && rootName == "System" && type == "System")
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
