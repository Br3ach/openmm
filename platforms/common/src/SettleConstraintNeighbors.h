#ifndef OPENMM_SETTLECONSTRAINTNEIGHBORS_H_
#define OPENMM_SETTLECONSTRAINTNEIGHBORS_H_

/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2026 the Authors.      *
 * Authors: OpenMM contributors                                                     *
 * Contributors:                                                              *
 *                                                                            *
 * This program is free software: you can redistribute it and/or modify       *
 * it under the terms of the GNU Lesser General Public License as published   *
 * by the Free Software Foundation, either version 3 of the License, or       *
 * (at your option) any later version.                                        *
 *                                                                            *
 * This program is distributed in the hope that it will be useful,            *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of             *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the              *
 * GNU Lesser General Public License for more details.                        *
 *                                                                            *
 * You should have received a copy of the GNU Lesser General Public License   *
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.      *
 * -------------------------------------------------------------------------- */

#include "openmm/OpenMMException.h"
#include <cstddef>
#include <utility>

namespace OpenMM {
namespace detail {

// Private temporary connectivity for particles whose active constraint degree is two.
class SettleConstraintNeighbors {
public:
    typedef std::pair<int, float>* iterator;
    typedef const std::pair<int, float>* const_iterator;

    SettleConstraintNeighbors() : count(0) {
    }
    std::size_t size() const {
        return count;
    }
    iterator begin() {
        return entries;
    }
    const_iterator begin() const {
        return entries;
    }
    iterator end() {
        return entries+count;
    }
    const_iterator end() const {
        return entries+count;
    }
    iterator find(int atom) {
        for (int i = 0; i < count; i++)
            if (entries[i].first == atom)
                return entries+i;
        return end();
    }
    const_iterator find(int atom) const {
        for (int i = 0; i < count; i++)
            if (entries[i].first == atom)
                return entries+i;
        return end();
    }
    float& operator[](int atom) {
        iterator entry = find(atom);
        if (entry != end())
            return entry->second;
        if (count == 2)
            throw OpenMMException("Internal error: too many SETTLE constraint neighbors");
        int position = (count == 1 && entries[0].first < atom ? 1 : 0);
        if (count == 1 && position == 0)
            entries[1] = entries[0];
        entries[position] = std::make_pair(atom, 0.0f);
        count++;
        return entries[position].second;
    }
    void clear() {
        count = 0;
    }
private:
    std::pair<int, float> entries[2];
    unsigned char count;
};

} // namespace detail
} // namespace OpenMM

#endif
