from pathlib import Path
import os,json,subprocess,uuid,hashlib
base=Path.cwd();work=base/'settle-opt';dest=work/'results/resumed';rows=[]
for variant in ['baseline','candidate']:
 folder=work/('runtime-'+variant);cache=folder/('cache-'+uuid.uuid4().hex);cache.mkdir()
 env=os.environ.copy();env['PATH']=r'C:\Windows\System32;C:\Windows';env['OPENMM_CACHE_DIR']=str(cache)
 for name in list(env):
  if name.upper().startswith('CUDA_PATH') or name.upper() in ['CUDA_HOME','CUDATOOLKIT_ROOT']:env.pop(name,None)
 p=subprocess.run([str(folder/'VerifyCudaRuntime.exe')],cwd=folder,env=env,capture_output=True,text=True,timeout=180)
 row=dict(variant=variant,returncode=p.returncode,stdout=p.stdout,stderr=p.stderr,isolated_path=env['PATH'],fresh_cache=True)
 if not p.returncode:
  data=json.loads(p.stdout);row['data']=data
  for name,path in data['loaded_modules'].items():
   assert Path(path).resolve().parent==folder.resolve(),(name,path)
   row.setdefault('dll_sha256',{})[name]=hashlib.sha256(Path(path).read_bytes()).hexdigest()
 rows.append(row);(dest/'runtime-isolation.json').write_text(json.dumps(rows,indent=2));print(json.dumps(row),flush=True)
 if p.returncode:raise RuntimeError(p.stdout+p.stderr)
print('Both isolated PME/runtime checks passed',flush=True)
