#ifndef OPENMM_MOLECULECONNECTIVITY_H_
#define OPENMM_MOLECULECONNECTIVITY_H_

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

#include <utility>
#include <vector>

namespace OpenMM {
namespace detail {

// Private temporary connectivity. Negative roots store component sizes.
class MoleculeConnectivity {
public:
    explicit MoleculeConnectivity(int numParticles) : parent(numParticles, -1), numMolecules(numParticles) {
    }
    void connect(int particle1, int particle2) {
        int root1 = findRoot(particle1);
        int root2 = findRoot(particle2);
        if (root1 == root2)
            return;
        if (parent[root1] > parent[root2])
            std::swap(root1, root2);
        parent[root1] += parent[root2];
        parent[root2] = root1;
        numMolecules--;
    }
    std::vector<std::vector<int> > createMolecules() {
        std::vector<int> moleculeForRoot(parent.size(), -1);
        std::vector<std::vector<int> > result(numMolecules);
        int nextMolecule = 0;
        // First encountered members establish minimum-atom component order.
        for (int i = 0; i < (int) parent.size(); i++) {
            int root = findRoot(i);
            int& molecule = moleculeForRoot[root];
            if (molecule == -1) {
                molecule = nextMolecule++;
                result[molecule].reserve(-parent[root]);
            }
            result[molecule].push_back(i);
        }
        return result;
    }
private:
    int findRoot(int particle) {
        int root = particle;
        while (parent[root] >= 0)
            root = parent[root];
        while (particle != root) {
            int next = parent[particle];
            parent[particle] = root;
            particle = next;
        }
        return root;
    }
    std::vector<int> parent;
    int numMolecules;
};

} // namespace detail
} // namespace OpenMM

#endif
