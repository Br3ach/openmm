/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2009-2025 Stanford University and the Authors.      *
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

#include "../platforms/common/src/SettleConstraintNeighbors.h"
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <iostream>
#include <iterator>
#include <limits>
#include <map>
#include <random>
#include <stdexcept>
#include <utility>
#include <vector>
using namespace std;
struct Constraint { double distance; int first,second; };
static_assert(sizeof(Constraint)==16,"Input record layout");
struct Input { vector<double> masses; vector<Constraint> constraints; };
using CompactRow=OpenMM::detail::SettleConstraintNeighbors;

struct Settle {int a,b,c;float d1,d2;};
static_assert(sizeof(Settle)==20,"Output record layout");
struct Result { vector<int> clusters; vector<Settle> atoms; vector<bool> isShake; uint64_t neighbors=14695981039346656037ULL; };
static uint64_t hashInt(uint64_t h,uint32_t x) {for(int i=0;i<4;i++)h=(h^((x>>(8*i))&255))*1099511628211ULL;return h;}
static uint32_t bits(float x) {uint32_t b;memcpy(&b,&x,4);return b;}
template<class Row> Result discover(const Input& input) {
    const int n=(int)input.masses.size();
    vector<int> atom1,atom2;
    vector<double> distance;
    vector<int> constraintCount(n,0);
    for(const auto& c:input.constraints) {
        if(input.masses[c.first]!=0||input.masses[c.second]!=0) {
            atom1.push_back(c.first);atom2.push_back(c.second);distance.push_back(c.distance);
            constraintCount[c.first]++;constraintCount[c.second]++;
        }
    }
    Result result;result.isShake.resize(n,false);
    {
        vector<Row> settleConstraints(n);
        for(int i=0;i<(int)atom1.size();i++) {
            if(constraintCount[atom1[i]]==2&&constraintCount[atom2[i]]==2) {
                settleConstraints[atom1[i]][atom2[i]]=(float)distance[i];
                settleConstraints[atom2[i]][atom1[i]]=(float)distance[i];
            }
        }
        // Same ascending clear cascade and short-circuit topology checks as IntegrationUtilities.
        for(int i=0;i<n;i++) {
            if(settleConstraints[i].size()==2) {
                int partner1=settleConstraints[i].begin()->first;
                int partner2=next(settleConstraints[i].begin())->first;
                if(settleConstraints[partner1].size()!=2||settleConstraints[partner2].size()!=2||
                   settleConstraints[partner1].find(partner2)==settleConstraints[partner1].end())
                    settleConstraints[i].clear();
                else if(i<partner1&&i<partner2)result.clusters.push_back(i);
            } else settleConstraints[i].clear();
        }
        for(int a:result.clusters) {
            int b=settleConstraints[a].begin()->first;
            int c=next(settleConstraints[a].begin())->first;
            float d12=settleConstraints[a].find(b)->second;
            float d13=settleConstraints[a].find(c)->second;
            float d23=settleConstraints[b].find(c)->second;
            if(d12==d13)result.atoms.push_back({a,b,c,d12,d23});
            else if(d12==d23)result.atoms.push_back({b,a,c,d12,d13});
            else if(d13==d23)result.atoms.push_back({c,a,b,d13,d12});
            else continue;
            result.isShake[a]=true;result.isShake[b]=true;result.isShake[c]=true;
        }
        // Hash after the measured stages. Include every remaining row in order.
        for(auto& row:settleConstraints) {
            result.neighbors=hashInt(result.neighbors,(uint32_t)row.size());
            for(auto entry:row) {result.neighbors=hashInt(result.neighbors,entry.first);result.neighbors=hashInt(result.neighbors,bits(entry.second));}
        }
    }
    return result;
}
static uint64_t fingerprint(const Result& result) {
    uint64_t h=14695981039346656037ULL;
    for(int c:result.clusters)h=hashInt(h,c);
    for(const auto& s:result.atoms) {h=hashInt(hashInt(hashInt(h,s.a),s.b),s.c);h=hashInt(hashInt(h,bits(s.d1)),bits(s.d2));}
    for(bool value:result.isShake)h=hashInt(h,value);
    return h;
}
static void compare(const Input& input) {
    auto a=discover<map<int,float>>(input),b=discover<CompactRow>(input);
    if(a.clusters!=b.clusters||a.atoms.size()!=b.atoms.size()||a.isShake!=b.isShake||a.neighbors!=b.neighbors)
        throw runtime_error("Classification differs");
    for(size_t i=0;i<a.atoms.size();i++)if(memcmp(&a.atoms[i],&b.atoms[i],sizeof(Settle))!=0)throw runtime_error("SETTLE output bytes differ");
    if(fingerprint(a)!=fingerprint(b))throw runtime_error("Fingerprint differs");
}
static void selftest() {
    mt19937 rng(20261024);
    for(int repeat=0;repeat<500;repeat++) {
        Input in;int n=1+rng()%90;in.masses.resize(n);
        for(auto& m:in.masses)m=rng()%4==0?0.0:1.0;
        for(int j=0;j<repeat%150;j++) {
            int a=rng()%n,b=rng()%n;
            in.constraints.push_back({0.1+0.01*(rng()%4),a,b});
        }
        compare(in);
    }
    const float nan=numeric_limits<float>::quiet_NaN(),inf=numeric_limits<float>::infinity();
    for(double a:{0.1,0.10000000001,0.2,-0.0,(double)nan,(double)inf})
        for(double b:{0.1,0.2,0.3,0.0,(double)nan,(double)inf})
            for(double c:{0.1,0.2,0.3}) {
                Input in;in.masses={16,1,1,0};in.constraints={{a,2,0},{b,1,2},{c,0,1}};compare(in);
            }
    Input cycle;cycle.masses.resize(20000,1.0);
    for(int i=0;i<20000;i++)cycle.constraints.push_back({0.1,i,(i+1)%20000});compare(cycle);
    for(int n:{4,5,12}) {
        Input path;path.masses.resize(n,1.0);
        for(int i=1;i<n;i++)path.constraints.push_back({0.1,i-1,i});
        compare(path);path.constraints.push_back({0.1,n-1,0});compare(path);
        reverse(path.constraints.begin(),path.constraints.end());
        for(auto& edge:path.constraints)swap(edge.first,edge.second);
        compare(path);
    }
    Input empty;compare(empty);
    Input massless;massless.masses={0,0,0};massless.constraints={{0.1,0,1},{0.1,0,2},{0.16,1,2}};compare(massless);
    Input rowCase;rowCase.masses={1,1,1};rowCase.constraints={{0.1,0,1},{0.2,0,1}};compare(rowCase);
    CompactRow row;row[8]=1;row[3]=2;row[8]=4;
    if(row.size()!=2||row.begin()->first!=3||next(row.begin())->first!=8||row.find(8)->second!=4)throw runtime_error("Row ordering/update");
    row.clear();row[2]=5;if(row.size()!=1||row.find(8)!=row.end())throw runtime_error("Row clear");
    cout<<"{\"selftest\":\"passed\",\"random_graphs\":500,\"float_triangles\":108,\"long_cycle_particles\":20000}"<<endl;
}
int main() {
    try {
        selftest();
        CompactRow row;row[1]=1;row[2]=2;
        try {row[3]=3;throw runtime_error("Overflow guard did not throw");}
        catch(const OpenMM::OpenMMException&) {}
        const CompactRow& fixed=row;
        if(fixed.size()!=2||fixed.find(1)==fixed.end()||fixed.find(3)!=fixed.end())throw runtime_error("Const lookup failed");
        return 0;
    } catch(const exception& error) {cerr<<error.what()<<endl;return 1;}
}
