#include "OpenMM.h"
#include "SimTKOpenMMRealType.h"
#include <cstdint>
#include <cstring>
#include <iomanip>
#include <iostream>
using namespace OpenMM;
using namespace std;
int main() { try {
    const int count=50000;
    System system;
    for (int i=0; i<count; i++) system.addParticle(1.0);
    CustomIntegrator integrator(0.001);
    Context context(system,integrator,Platform::getPlatformByName("Reference"));
    context.setPositions(vector<Vec3>(count));
    for(int trial=0;trial<101;trial++) {
        int seed=trial==100?20261024:trial;
        context.setVelocitiesToTemperature(300.0,seed);
        auto velocities=context.getState(State::Velocities).getVelocities();
        uint64_t hash=14695981039346656037ULL;double kinetic=0;
        for(const auto& v:velocities) {
            kinetic+=0.5*v.dot(v);
            for(int d=0;d<3;d++) {
                double value=v[d];uint64_t bits;memcpy(&bits,&value,sizeof(bits));
                for(int byte=0;byte<8;byte++)hash=(hash^((bits>>(byte*8))&255))*1099511628211ULL;
            }
        }
        double temperature=2*kinetic/(3*count*BOLTZ);
        cout<<setprecision(17)<<"{\"seed\":"<<seed<<",\"temperature\":"<<temperature<<",\"fingerprint\":\""<<hex<<hash<<dec<<"\"}"<<endl;
    }
    return 0;
} catch(const exception& e) {cerr<<e.what()<<endl;return 1;}}
