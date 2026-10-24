from pathlib import Path
import json
from run_monitored import run_case,work,dest
rows=[]
cases=[('candidate','TestCudaMonteCarloBarostat','mixed')]
cases += [('candidate','TestCuda'+name,'double') for name in ['Settle','VirtualSites','Checkpoints','VerletIntegrator','LangevinMiddleIntegrator','CustomIntegrator','VariableVerletIntegrator','NoseHooverIntegrator','MonteCarloBarostat']]
cases += [('baseline','TestCuda'+name,precision) for precision in ['single','mixed','double'] for name in ['Settle','VirtualSites','Checkpoints'] if not(name=='Settle' and precision=='single')]
cases += [(variant,'TestCudaCustomIntegratorFixedSeed','single') for variant in ['baseline','candidate']]
cases += [(variant,'TestCudaFFT3D',precision) for precision in ['single','mixed','double'] for variant in ['baseline','candidate']]
for variant,test,precision in cases:
 rows.append(run_case(variant,test,precision,600));(dest/'additional-tests.json').write_text(json.dumps(rows,indent=2))
print(json.dumps(dict(event='additional_complete',cases=len(rows),passes=sum(not r['returncode'] and not r['timeout'] for r in rows))),flush=True)
