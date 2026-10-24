#ifndef OPENMM_VIRTUALSITESNAPSHOT_H_
#define OPENMM_VIRTUALSITESNAPSHOT_H_
/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2019-2024 Stanford University and the Authors.      *
 * Authors: Peter Eastman                                                     *
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


#include "openmm/System.h"
#include "openmm/VirtualSite.h"
#include <array>
#include <typeinfo>
#include <vector>

namespace OpenMM { namespace detail {
// Immutable reorder metadata. Unused legacy weight slots stay zero.
class VirtualSiteSnapshot {
public:
    explicit VirtualSiteSnapshot(const System& system) {
        size_t count = 0, particleCount = 0;
        for (int i = 0; i < system.getNumParticles(); i++)
            if (system.isVirtualSite(i)) {
                count++;
                particleCount += size_t(1)+size_t(system.getVirtualSite(i).getNumParticles());
            }
        types.reserve(count);
        offsets.reserve(count+1);
        weights.reserve(count);
        particles.reserve(particleCount);
        for (int i = 0; i < system.getNumParticles(); i++) {
            if (!system.isVirtualSite(i))
                continue;
            const VirtualSite& site = system.getVirtualSite(i);
            types.push_back(&typeid(site));
            offsets.push_back(particles.size());
            particles.push_back(i);
            for (int j = 0; j < site.getNumParticles(); j++)
                particles.push_back(site.getParticle(j));
            std::array<double, 3> value = {{0, 0, 0}};
            if (const auto* average = dynamic_cast<const TwoParticleAverageSite*>(&site)) {
                value[0] = average->getWeight(0);
                value[1] = average->getWeight(1);
            }
            else if (const auto* average = dynamic_cast<const ThreeParticleAverageSite*>(&site)) {
                for (int j = 0; j < 3; j++)
                    value[j] = average->getWeight(j);
            }
            else if (const auto* plane = dynamic_cast<const OutOfPlaneSite*>(&site)) {
                value[0] = plane->getWeight12();
                value[1] = plane->getWeight13();
                value[2] = plane->getWeightCross();
            }
            weights.push_back(value);
        }
        offsets.push_back(particles.size());
    }
    int getNumParticleGroups() const {
        return types.size();
    }
    void getParticlesInGroup(int index, std::vector<int>& result) const {
        result.assign(particles.begin()+offsets[index], particles.begin()+offsets[index+1]);
    }
    bool areGroupsIdentical(int group1, int group2) const {
        if (types[group1] != types[group2])
            return false;
        // Equal RTTI means the same legacy cast branch and weight count.
        for (int i = 0; i < 3; i++)
            if (weights[group1][i] != weights[group2][i])
                return false;
        return true;
    }
private:
    std::vector<const std::type_info*> types;
    std::vector<size_t> offsets;
    std::vector<int> particles;
    std::vector<std::array<double, 3> > weights;
};
} }
#endif
