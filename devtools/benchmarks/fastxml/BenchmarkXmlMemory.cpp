#include "openmm/State.h"
#include "openmm/System.h"
#include "openmm/NonbondedForce.h"
#include "openmm/CustomBondForce.h"
#include "openmm/serialization/XmlSerializer.h"
#include <chrono>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <malloc.h>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <sys/resource.h>
#include <vector>
using namespace OpenMM;
using namespace std;

static long statusKB(const string& field) {
    ifstream f("/proc/self/status");
    string line;
    while (getline(f,line)) if (line.rfind(field+":",0)==0) {
        istringstream s(line.substr(field.size()+1)); long value; s>>value; return value;
    }
    throw runtime_error("Cannot read memory metric");
}
static void resetPeak() {
    ofstream f("/proc/self/clear_refs"); f<<"5\n"; f.close();
    if (!f) throw runtime_error("Cannot reset Linux RSS high-water mark");
}
static uint64_t digest=1469598103934665603ULL;
template<class T> static void hashValue(const T& v) {
    const unsigned char* p=reinterpret_cast<const unsigned char*>(&v);
    for(size_t j=0;j<sizeof(v);++j) {digest^=p[j]; digest*=1099511628211ULL;}
}
static void hashVec(const Vec3& v) {for(int j=0;j<3;j++) hashValue(v[j]);}
static unique_ptr<State> makeState(int n,bool forces) {
    State::StateBuilder b(2.5,1234567);
    b.setPeriodicBoxVectors(Vec3(40,0,0),Vec3(0,40,0),Vec3(0,0,40));
    vector<Vec3> v(n);
    for(int i=0;i<n;i++) v[i]=Vec3((i%1009)*0.03125,(i%1013)*0.015625,(i%1019)*0.0078125);
    b.setPositions(v);
    for(int i=0;i<n;i++) v[i]=Vec3(((i%97)-48)*0.0078125,((i%89)-44)*0.015625,((i%83)-41)*0.03125);
    b.setVelocities(v);
    if(forces) {
        for(int i=0;i<n;i++) v[i]=Vec3(((i%71)-35)*0.0625,((i%67)-33)*0.125,((i%61)-30)*0.25);
        b.setForces(v);
    }
    return unique_ptr<State>(new State(b.getState()));
}
static unique_ptr<System> makeSystem(int n,bool fallback) {
    unique_ptr<System> s(new System());
    s->setDefaultPeriodicBoxVectors(Vec3(40,0,0),Vec3(0,40,0),Vec3(0,0,40));
    auto* nb=new NonbondedForce();
    nb->setNonbondedMethod(NonbondedForce::CutoffPeriodic);
    nb->setCutoffDistance(1.0);
    for(int i=0;i<n;i++) {
        bool oxygen=i%3==0;
        s->addParticle(oxygen?15.999:1.008);
        nb->addParticle(oxygen?-0.834:0.417,oxygen?0.315075:1.0,oxygen?0.635968:0.0);
    }
    for(int i=0;i+2<n;i+=3) {
        s->addConstraint(i,i+1,0.09572); s->addConstraint(i,i+2,0.09572); s->addConstraint(i+1,i+2,0.15139);
        nb->addException(i,i+1,0,1,0); nb->addException(i,i+2,0,1,0); nb->addException(i+1,i+2,0,1,0);
    }
    s->addForce(nb);
    if(fallback) { // Last force deliberately triggers late whole-document fallback.
        auto* custom=new CustomBondForce("0.5*k*(r-r0)^2");
        custom->addGlobalParameter("k",100); custom->addGlobalParameter("r0",0.15);
        custom->addBond(0,1); s->addForce(custom);
    }
    return s;
}
static void checkState(const State& s,int n,bool forces) {
    if(s.getPositions().size()!=n || s.getVelocities().size()!=n || (forces && s.getForces().size()!=n)) throw runtime_error("State count mismatch");
    hashValue(s.getDataTypes()); hashValue(s.getTime()); hashValue(s.getStepCount());
    Vec3 a,b,c; s.getPeriodicBoxVectors(a,b,c);hashVec(a);hashVec(b);hashVec(c);
    for(auto& v:s.getPositions()) hashVec(v);
    for(auto& v:s.getVelocities()) hashVec(v);
    if(forces) for(auto& v:s.getForces()) hashVec(v);
}
static void checkSystem(const System& s,int n,bool fallback) {
    if(s.getNumParticles()!=n || s.getNumConstraints()!=3*(n/3) || s.getNumForces()!=(fallback?2:1)) throw runtime_error("System count mismatch");
    hashValue(s.getNumParticles()); hashValue(s.getNumForces());
    Vec3 a,b,c;s.getDefaultPeriodicBoxVectors(a,b,c);hashVec(a);hashVec(b);hashVec(c);
    for(int i=0;i<n;i++) hashValue(s.getParticleMass(i));
    for(int i=0;i<s.getNumConstraints();i++) {int a,b;double d;s.getConstraintParameters(i,a,b,d);hashValue(a);hashValue(b);hashValue(d);}
    auto& nb=dynamic_cast<const NonbondedForce&>(s.getForce(0));
    if(nb.getNumParticles()!=n || nb.getNumExceptions()!=3*(n/3)) throw runtime_error("Force count mismatch");
    hashValue(nb.getNonbondedMethod());hashValue(nb.getCutoffDistance());
    for(int i=0;i<n;i++) {double q,s,e;nb.getParticleParameters(i,q,s,e);hashValue(q);hashValue(s);hashValue(e);}
    for(int i=0;i<nb.getNumExceptions();i++) {int a,b;double q,s,e;nb.getExceptionParameters(i,a,b,q,s,e);hashValue(a);hashValue(b);hashValue(q);hashValue(s);hashValue(e);}
    if(fallback) {
        auto& f=dynamic_cast<const CustomBondForce&>(s.getForce(1));
        for(char c:f.getEnergyFunction()) hashValue(c);
        hashValue(f.getNumBonds());
        for(int i=0;i<f.getNumGlobalParameters();i++) hashValue(f.getGlobalParameterDefaultValue(i));
        int a,b;vector<double> v;f.getBondParameters(0,a,b,v);hashValue(a);hashValue(b);
    }
}
int main(int argc,char** argv) {
    try {
        if(argc!=6) throw runtime_error("usage: bench operation kind count path input_mode");
        string op=argv[1],kind=argv[2],path=argv[4],inputMode=argv[5];
        int n=stoi(argv[3]);bool state=kind.rfind("state",0)==0, forces=kind=="state-pvf", fallback=kind=="system-fallback";
        unique_ptr<State> st;unique_ptr<System> sy;
        if(op=="write" || op=="generate") {
            if(state) st=makeState(n,forces);else sy=makeSystem(n,fallback);
            malloc_trim(0); // Discard released setup allocations before measuring serialization.
        }
        long baseline=statusKB("VmRSS");resetPeak();
        auto start=chrono::steady_clock::now();
        ifstream input;unique_ptr<stringstream> memoryInput;string retainedXml;
        if(op=="read") {
            input.open(path,ios::binary);if(!input) throw runtime_error("Cannot open XML input");
            istream* in=&input;
            if(inputMode=="copied-string") {
                retainedXml.assign(istreambuf_iterator<char>(input),istreambuf_iterator<char>());
                memoryInput.reset(new stringstream(retainedXml)); in=memoryInput.get();
            }
            if(state) st.reset(XmlSerializer::deserialize<State>(*in));
            else sy.reset(XmlSerializer::deserialize<System>(*in));
        } else {
            ofstream out(path,ios::binary);if(!out) throw runtime_error("Cannot open XML output");
            if(state) XmlSerializer::serialize(st.get(),"State",out);
            else XmlSerializer::serialize(sy.get(),"System",out);
            out.flush();if(!out) throw runtime_error("XML output failed");
        }
        double elapsed=chrono::duration<double>(chrono::steady_clock::now()-start).count();
        long peak=statusKB("VmHWM"),retained=statusKB("VmRSS");
        rusage ru;getrusage(RUSAGE_SELF,&ru);
        // Validation is deliberately outside the timed and memory-measured interval.
        if(state) checkState(*st,n,forces);else checkSystem(*sy,n,fallback);
        cout<<setprecision(9)<<"{\"seconds\":"<<elapsed<<",\"baseline_kib\":"<<baseline<<",\"peak_kib\":"<<peak<<",\"retained_kib\":"<<retained<<",\"lifetime_peak_kib\":"<<ru.ru_maxrss<<",\"digest\":\""<<hex<<digest<<"\"}"<<endl;
    } catch(const exception& e) {cerr<<e.what()<<endl;return 1;}
}
