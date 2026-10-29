#ifndef OPENMM_FAST_STATE_XML_H_
#define OPENMM_FAST_STATE_XML_H_
#include <iosfwd>
#if defined(_MSC_VER) && defined(_M_X64)
#include <intrin.h>
#define OPENMM_XML_CALLER() _ReturnAddress()
#else
#define OPENMM_XML_CALLER() nullptr
#endif
namespace OpenMM {
class State;
class System;
class StateProxy;
class SerializationNode;
State* deserializeStateXmlFast(std::istream& stream);
System* deserializeSystemXmlFast(std::istream& stream);
bool canDeferCore28StateXml(const State&, const StateProxy&, const SerializationNode&, const void* caller);
bool deferCore28StateXml(const State&, SerializationNode&);
void cancelCore28StateXml(const SerializationNode*) noexcept;
void disableCore28StateXml() noexcept;
bool writeCore28StateXml(const SerializationNode&, std::ostream&, const void* caller);
}
#endif
