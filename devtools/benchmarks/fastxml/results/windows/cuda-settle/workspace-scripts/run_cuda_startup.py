from pathlib import Path
import os,json,subprocess,time,statistics,hashlib,ctypes
base=Path.cwd();work=base/'settle-opt';dest=work/'results/cuda-startup';dest.mkdir(exist_ok=True)
env=os.environ.copy();env['PATH']=r'C:\Program Files\NVIDIA GPU Computing Toolkit\CUDA\v13.4\bin\x64'+os.pathsep+env['PATH'];env['OPENMM_CACHE_DIR']=str(work/'cuda-cache')
exes={v:work/v/'BenchmarkCudaSettleStartup.exe' for v in ['baseline','candidate']};assert len({hashlib.sha256(p.read_bytes()).hexdigest() for p in exes.values()})==1
rows=[];warmups=[];cases=[(250000,'unconstrained'),(250000,'water'),(500000,'water')]
def run(v,n,kind):
 gpu=subprocess.check_output(['nvidia-smi','--query-gpu=memory.free,utilization.gpu','--format=csv,noheader,nounits'],text=True).strip()
 if int(gpu.split(',')[0])<2048:raise RuntimeError('Insufficient free VRAM for scaled case: '+gpu)
 p=subprocess.run([str(exes[v]),str(n),kind],cwd=exes[v].parent,env=env,capture_output=True,text=True,timeout=180)
 if p.returncode:raise RuntimeError(v+': '+p.stdout+p.stderr)
 row=json.loads(p.stdout);row.update(variant=v,gpu_status_before=gpu);return row
for n,kind in cases:
 for v in ['baseline','candidate']:warmups.append(run(v,n,kind));(dest/'warmup.json').write_text(json.dumps(warmups,indent=2))
 for repeat in range(3):
  for v in (['baseline','candidate'] if repeat%2==0 else ['candidate','baseline']):
   row=run(v,n,kind);row['repeat']=repeat+1;rows.append(row);(dest/'raw.json').write_text(json.dumps(rows,indent=2));print(json.dumps(row),flush=True)
summary=[]
for n,kind in cases:
 selected=[r for r in rows if r['particles']==n and r['kind']==kind];assert len({r['constraints'] for r in selected})==1
 item=dict(particles=n,kind=kind,constraints=selected[0]['constraints'])
 for variant in ['baseline','candidate']:
  chosen=[r for r in selected if r['variant']==variant];keys=['construction_seconds']+[k for k in chosen[0] if k.endswith('_kib')]
  item[variant]={k:statistics.median(r[k] for r in chosen) for k in keys}
 summary.append(item)
(dest/'summary.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary,indent=2))
gpu=subprocess.check_output(['nvidia-smi','--query-gpu=memory.free,utilization.gpu','--format=csv,noheader,nounits'],text=True).strip()
(dest/'full-size-status.json').write_text(json.dumps(dict(performed=False,reason='Scaled benchmark is bounded to500000particles; active GPU workload, full-topology evidence retained in host probe',gpu_status_after=gpu),indent=2))
