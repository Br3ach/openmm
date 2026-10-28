#ifndef OPENMM_FAST_STATE_XML_H_
#define OPENMM_FAST_STATE_XML_H_
#include <iosfwd>
namespace OpenMM {
class State;
class System;
State* deserializeStateXmlFast(std::istream& stream);
System* deserializeSystemXmlFast(std::istream& stream);
}
#endif
