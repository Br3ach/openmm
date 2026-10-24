from pathlib import Path
import shutil,json
base=Path.cwd();work=base/'settle-opt';sdk=Path(r'C:\Program Files\NVIDIA GPU Computing Toolkit\CUDA\v13.4');runtime=['cufft64_12.dll','nvrtc64_130_0.dll','nvrtc-builtins64_134.dll']
for variant in ['baseline','candidate']:
 stage=work/('runtime-'+variant);stage.mkdir(exist_ok=True)
 for name in runtime:shutil.copy2(sdk/'bin/x64'/name,stage/name)
 for name in ['OpenMM.dll','OpenMMCUDA.dll','VerifyCudaRuntime.exe','TestCudaFFT3D.exe','TestCudaSettleSnapshots.exe','BenchmarkCudaSettleStartup.exe']:shutil.copy2(work/variant/name,stage/name)
shutil.copy2(sdk/'EULA.txt',work/'NVIDIA-CUDA-EULA.txt')
print('Staged isolated runtime folders with all three CUDA dependencies',flush=True)
