from pathlib import Path
import shutil
base=Path.cwd();work=base/'settle-opt';scratch=work/'fixed-seed';scratch.mkdir(exist_ok=True)
header=(base/'source/tests/TestCustomIntegrator.h').read_text();needle='context.setVelocitiesToTemperature(targetTemperature);';assert header.count(needle)==1
(scratch/'TestCustomIntegrator.h').write_text(header.replace(needle,'context.setVelocitiesToTemperature(targetTemperature, 20261024);'))
shutil.copy2(base/'source/platforms/cuda/tests/TestCudaCustomIntegrator.cpp',scratch/'TestCudaCustomIntegratorFixedSeed.cpp')
s=(work/'build_snapshots.py').read_text().replace("root/'platforms/cuda/tests/TestCudaSettleSnapshots.cpp'","work/'fixed-seed/TestCudaCustomIntegratorFixedSeed.cpp'")
(work/'build_fixed.py').write_text(s)
s=(work/'build_snapshots.cmd').read_text().replace('build_snapshots.py','build_fixed.py');(work/'build_fixed.cmd').write_text(s)
