from pathlib import Path
import subprocess,json,concurrent.futures,shutil
base=Path.cwd();work=base/'settle-opt';root=base/'source';cfg=json.loads((work/'build-config.json').read_text());flags=[x for x in cfg['flags'] if not x.startswith('/D')]+['/DWIN32','/D_WINDOWS','/DNDEBUG','/I'+str(root/'platforms/cuda/tests')]
def run(args):
 p=subprocess.run(args,capture_output=True,text=True)
 if p.returncode:raise RuntimeError(p.stdout+p.stderr)
 return p
sources=[root/'platforms/cuda/tests'/('TestCuda'+name+'.cpp') for name in ['Settle','VirtualSites','Checkpoints','VerletIntegrator','LangevinMiddleIntegrator','CustomIntegrator','VariableVerletIntegrator','NoseHooverIntegrator','MonteCarloBarostat']]
sources += [root/'tests/TestSettleConnectivity.cpp',root/'devtools/benchmarks/fastxml/BenchmarkSettleStorage.cpp']
def build(src):
 out=work/'candidate';exe=out/(src.stem+'.exe');run(['cl',*flags,'/Fo'+str(out/(src.stem+'.obj')),'/Fe'+str(exe),str(src),str(work/'baseline/OpenMMCUDA.lib'),str(work/'baseline/OpenMM.lib'),'Psapi.lib']);shutil.copy2(exe,work/'baseline'/exe.name);print('Built',exe.name,flush=True)
with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:list(pool.map(build,sources))
print('Test callers built once and copied byte-identically',flush=True)
