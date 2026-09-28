#ifndef OPENMM_BINARY_SERIALIZER_H_
#define OPENMM_BINARY_SERIALIZER_H_

/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * Experimental portable binary serialization for OpenMM objects.            *
 *                                                                            *
 * This API deliberately coexists with XmlSerializer.  XML remains the        *
 * human-readable/interchange format; BinarySerializer is intended for large  *
 * objects where constructing/parsing XML is prohibitively expensive.         *
 * -------------------------------------------------------------------------- */

#include "openmm/internal/windowsExport.h"
#include <iosfwd>
#include <sstream>
#include <typeinfo>

namespace OpenMM {

/**
 * Serialize OpenMM objects in a portable, versioned binary format.
 *
 * The implementation has native streaming fast paths for State, System,
 * NonbondedForce, HarmonicBondForce, HarmonicAngleForce, and
 * PeriodicTorsionForce.  Any other object registered with OpenMM's existing
 * SerializationProxy registry is encoded through a compact binary
 * SerializationNode fallback, so plugins and uncommon force types continue to
 * work without being modified.
 *
 * The wire representation is canonical little-endian and does not contain
 * platform pointers, CUDA state, or other hardware-specific data.  It is
 * therefore distinct from Context::createCheckpoint().
 */
class OPENMM_EXPORT BinarySerializer {
public:
    /** Serialize an object to a binary stream. */
    template <class T>
    static void serialize(const T* object, std::ostream& stream) {
        serializeImpl(object, typeid(*object), stream);
    }

    /** Deserialize an object from a binary stream.  The caller owns it. */
    template <class T>
    static T* deserialize(std::istream& stream) {
        return reinterpret_cast<T*>(deserializeImpl(stream));
    }

    /**
     * Clone an object through the binary representation.  This convenience
     * method materializes the complete binary representation in memory before
     * deserializing it.
     */
    template <class T>
    static T* clone(const T& object) {
        std::stringstream stream(std::ios_base::in | std::ios_base::out | std::ios_base::binary);
        serialize<T>(&object, stream);
        stream.seekg(0);
        return deserialize<T>(stream);
    }

private:
    static void serializeImpl(const void* object, const std::type_info& type, std::ostream& stream);
    static void* deserializeImpl(std::istream& stream);
};

} // namespace OpenMM

#endif /* OPENMM_BINARY_SERIALIZER_H_ */
