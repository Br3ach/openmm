@echo off
call "C:\Program Files\Microsoft Visual Studio\18\Professional\VC\Auxiliary\Build\vcvars64.bat"
if errorlevel 1 exit /b %errorlevel%
"C:\Program Files\Microsoft Visual Studio\18\Professional\Common7\IDE\CommonExtensions\Microsoft\CMake\CMake\bin\cmake.exe" -S source -B settle-opt\cmake -G "Visual Studio 18 2026" -A x64 -DOPENMM_BUILD_CUDA_LIB=ON -DOPENMM_BUILD_OPENCL_LIB=OFF -DOPENMM_BUILD_HIP_LIB=OFF -DOPENMM_BUILD_CPU_LIB=OFF -DOPENMM_BUILD_AMOEBA_PLUGIN=OFF -DOPENMM_BUILD_RPMD_PLUGIN=OFF -DOPENMM_BUILD_DRUDE_PLUGIN=OFF -DOPENMM_BUILD_PME_PLUGIN=OFF -DOPENMM_BUILD_C_AND_FORTRAN_WRAPPERS=OFF -DOPENMM_BUILD_PYTHON_WRAPPERS=OFF -DOPENMM_BUILD_REFERENCE_TESTS=OFF -DBUILD_TESTING=ON "-DCUDAToolkit_ROOT=C:/Program Files/NVIDIA GPU Computing Toolkit/CUDA/v13.4"
