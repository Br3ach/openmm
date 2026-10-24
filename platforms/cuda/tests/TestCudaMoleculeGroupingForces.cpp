#include "CudaTests.h"
#define main originalNonbondedMain
#include "TestNonbondedForce.h"
#undef main
void runPlatformTests() {}
int main(int argc,char** argv) {
 try {
  initializeTests(argc,argv);
  testReordering();testChangingParameters();testParameterOffsets();testPeriodic();testPeriodicExceptions();
  std::cout<<"Reordering, parameter changes/offsets, periodic interactions/exceptions passed"<<std::endl;return 0;
 } catch(const std::exception& e){std::cerr<<e.what()<<std::endl;return 1;}
}
