from pathlib import Path
import os,json,subprocess,shutil
base=Path.cwd();work=base/'settle-opt';root=base/'source';cfg=json.loads((work/'build-config.json').read_text());flags=[x for x in cfg['flags'] if not x.startswith('/D')]+['/DWIN32','/D_WINDOWS','/DNDEBUG','/I'+str(root/'platforms/cuda/tests')]
src=root/'platforms/cuda/tests/TestCudaSettleSnapshots.cpp';out=work/'candidate';exe=out/(src.stem+'.exe')
p=subprocess.run(['cl',*flags,'/Fo'+str(out/(src.stem+'.obj')),'/Fe'+str(exe),str(src),str(work/'baseline/OpenMMCUDA.lib'),str(work/'baseline/OpenMM.lib'),'Psapi.lib'],capture_output=True,text=True)
if p.returncode:raise RuntimeError(p.stdout+p.stderr)
shutil.copy2(exe,work/'baseline'/exe.name);print('Snapshot executable built once/copied',flush=True)
