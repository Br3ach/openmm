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

#include "../platforms/common/src/VirtualSiteSnapshot.h"
#include <iostream>
#include <limits>
#include <memory>
#include <stdexcept>
using namespace OpenMM;
using namespace std;
class LegacySnapshot {
public:
    LegacySnapshot(const System& system) {
        for (int i = 0; i < system.getNumParticles(); i++) {
            if (system.isVirtualSite(i)) {
                const VirtualSite& vsite = system.getVirtualSite(i);
                siteTypes.push_back(&typeid(vsite));
                vector<int> particles;
                particles.push_back(i);
                for (int j = 0; j < vsite.getNumParticles(); j++)
                    particles.push_back(vsite.getParticle(j));
                siteParticles.push_back(particles);
                vector<double> weights;
                if (dynamic_cast<const TwoParticleAverageSite*>(&vsite) != NULL) {
                    // A two particle average.

                    const TwoParticleAverageSite& site = dynamic_cast<const TwoParticleAverageSite&>(vsite);
                    weights.push_back(site.getWeight(0));
                    weights.push_back(site.getWeight(1));
                }
                else if (dynamic_cast<const ThreeParticleAverageSite*>(&vsite) != NULL) {
                    // A three particle average.

                    const ThreeParticleAverageSite& site = dynamic_cast<const ThreeParticleAverageSite&>(vsite);
                    weights.push_back(site.getWeight(0));
                    weights.push_back(site.getWeight(1));
                    weights.push_back(site.getWeight(2));
                }
                else if (dynamic_cast<const OutOfPlaneSite*>(&vsite) != NULL) {
                    // An out of plane site.

                    const OutOfPlaneSite& site = dynamic_cast<const OutOfPlaneSite&>(vsite);
                    weights.push_back(site.getWeight12());
                    weights.push_back(site.getWeight13());
                    weights.push_back(site.getWeightCross());
                }
                siteWeights.push_back(weights);
            }
        }
    }
    int getNumParticleGroups() {
        return siteTypes.size();
    }
    void getParticlesInGroup(int index, std::vector<int>& particles) {
        particles = siteParticles[index];
    }
    bool areGroupsIdentical(int group1, int group2) {
        if (siteTypes[group1] != siteTypes[group2])
            return false;
        int numParticles = siteWeights[group1].size();
        if (siteWeights[group2].size() != numParticles)
            return false;
        for (int i = 0; i < numParticles; i++)
            if (siteWeights[group1][i] != siteWeights[group2][i])
                return false;
        return true;
    }
private:
    vector<const type_info*> siteTypes;
    vector<vector<int> > siteParticles;
    vector<vector<double> > siteWeights;
};

class UnknownSite : public VirtualSite { public: UnknownSite(const vector<int>& p){setParticles(p);} };
class DerivedAverage : public ThreeParticleAverageSite {public: DerivedAverage():ThreeParticleAverageSite(0,1,2,.5,.25,.25){} };
static void compare(LegacySnapshot& legacy, detail::VirtualSiteSnapshot& compact) {
    if(legacy.getNumParticleGroups()!=compact.getNumParticleGroups())throw runtime_error("Group count differs");
    vector<int> a={99,98},b=a;
    for(int i=0;i<legacy.getNumParticleGroups();i++){
        legacy.getParticlesInGroup(i,a);compact.getParticlesInGroup(i,b);if(a!=b)throw runtime_error("Snapshot particle order differs");
        for(int j=0;j<legacy.getNumParticleGroups();j++)
            if(legacy.areGroupsIdentical(i,j)!=compact.areGroupsIdentical(i,j))throw runtime_error("Snapshot identity differs");
    }
}
int main(){try{
    System empty;LegacySnapshot e(empty);detail::VirtualSiteSnapshot c(empty);compare(e,c);
    unique_ptr<LegacySnapshot> before;unique_ptr<detail::VirtualSiteSnapshot> after;
    {
        System system;for(int i=0;i<50;i++)system.addParticle(1);
        for(int i=4;i<10;i++)system.setVirtualSite(i,new TwoParticleAverageSite(0,1,i==5?-0.0:0.0,1));
        for(int i=10;i<15;i++)system.setVirtualSite(i,new ThreeParticleAverageSite(0,1,2,.5,.25,.25));
        system.setVirtualSite(15,new DerivedAverage());
        for(int i=16;i<20;i++)system.setVirtualSite(i,new OutOfPlaneSite(0,1,2,.2,.3,.4));
        system.setVirtualSite(20,new UnknownSite({}));system.setVirtualSite(21,new UnknownSite({2,1,2,0,3}));
        system.setVirtualSite(22,new LocalCoordinatesSite({0,1,2,3},{1,0,0,0},{-1,1,0,0},{-1,0,1,0},Vec3(.1,.2,.3)));
        system.setVirtualSite(23,new LocalCoordinatesSite({0,1,2,3},{1,0,0,0},{-1,1,0,0},{-1,0,1,0},Vec3(.2,.5,.6)));
        system.setVirtualSite(24,new TwoParticleAverageSite(0,1,numeric_limits<double>::quiet_NaN(),1));
        system.setVirtualSite(25,new TwoParticleAverageSite(0,1,numeric_limits<double>::quiet_NaN(),1));
        system.setVirtualSite(26,new ThreeParticleAverageSite(0,1,2,numeric_limits<double>::infinity(),-.0,1));
        system.setVirtualSite(27,new ThreeParticleAverageSite(0,1,2,numeric_limits<double>::infinity(),.0,1));
        vector<int> many;for(int i=0;i<1000;i++)many.push_back(i%4);system.setVirtualSite(28,new UnknownSite(many));
        before.reset(new LegacySnapshot(system));after.reset(new detail::VirtualSiteSnapshot(system));compare(*before,*after);
        system.setVirtualSite(10,new TwoParticleAverageSite(3,2,.2,.8));compare(*before,*after);
    }
    compare(*before,*after);
    cout<<"Virtual-site snapshots match legacy identity/order, RTTI, arbitrary parents, ignored local weights, nonfinite weights and immutable ownership"<<endl;return 0;
}catch(const exception& e){cerr<<e.what()<<endl;return 1;}}
