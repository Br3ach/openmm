#ifndef OPENMM_FAST_STATE_XML_H_
#define OPENMM_FAST_STATE_XML_H_

#include <iosfwd>

namespace OpenMM {

class SerializationNode;
class State;
class System;

// StateProxy uses these to defer large canonical State XML documents to the
// direct writer in XmlSerializer.cpp, avoiding a SerializationNode per vector.
void clearFastStateXmlSerialization();
void registerFastStateXmlSerialization(const State* state, SerializationNode* node);

// Streaming readers used by XmlSerializer for canonical State/System XML.
State* deserializeStateXmlFast(std::istream& stream);
System* deserializeSystemXmlFast(std::istream& stream);

} // namespace OpenMM

#endif /* OPENMM_FAST_STATE_XML_H_ */
