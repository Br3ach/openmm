from pathlib import Path
import os,subprocess,json,time
base=Path.cwd();work=base/'settle-opt';env=os.environ.copy();env['PATH']=r'C:\Program Files\NVIDIA GPU Computing Toolkit\CUDA\v13.4\bin\x64'+os.pathsep+env['PATH'];env['OPENMM_CACHE_DIR']=str(work/'cuda-cache')
rows=[]
for precision in ['single','mixed','double']:
 for name in ['Settle','VirtualSites','Checkpoints','VerletIntegrator','LangevinMiddleIntegrator','CustomIntegrator','VariableVerletIntegrator','NoseHooverIntegrator','MonteCarloBarostat']:
  exe=work/'candidate'/('TestCuda'+name+'.exe');start=time.perf_counter()
  p=subprocess.run([str(exe),precision],cwd=exe.parent,env=env,capture_output=True,text=True,timeout=600)
  row=dict(test=exe.stem,precision=precision,variant='candidate',returncode=p.returncode,seconds=time.perf_counter()-start,stdout=p.stdout,stderr=p.stderr);rows.append(row)
  (work/'results/gpu-tests-initial.json').write_text(json.dumps(rows,indent=2));print(row,flush=True)
print('Finished GPU suites',len(rows),'passes',sum(r['returncode']==0 for r in rows),flush=True)
