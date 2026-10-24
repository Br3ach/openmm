#ifndef OPENMM_FAST_STATE_XML_H_
#define OPENMM_FAST_STATE_XML_H_

#include <iosfwd>
#include <exception>

namespace OpenMM {

class SerializationNode;
class State;
class System;

// A generic proxy has already accepted this subtree.  Preserve its exception
// without retrying the entire document and invoking the proxy a second time.
struct FastXmlProxyException {
    explicit FastXmlProxyException(std::exception_ptr exception) : exception(exception) {}
    std::exception_ptr exception;
};

// StateProxy uses these to defer large canonical State XML documents to the
// direct writer in XmlSerializer.cpp, avoiding a SerializationNode per vector.
void clearFastStateXmlSerialization();
void registerFastStateXmlSerialization(const State* state, SerializationNode* node);

// Streaming readers used by XmlSerializer for canonical State/System XML.
State* deserializeStateXmlFast(std::istream& stream);
System* deserializeSystemXmlFast(std::istream& stream);

} // namespace OpenMM

#endif /* OPENMM_FAST_STATE_XML_H_ */
