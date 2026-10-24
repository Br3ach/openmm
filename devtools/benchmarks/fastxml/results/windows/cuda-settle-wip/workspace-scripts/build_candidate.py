from pathlib import Path
import json,subprocess
base=Path.cwd();work=base/'settle-opt';c=json.loads((work/'build-config.json').read_text());out=work/'candidate'
def run(args):
 p=subprocess.run(args,capture_output=True,text=True)
 if p.returncode:raise RuntimeError(p.stdout+p.stderr)
 return p
objects=[]
for i,src in enumerate(c['sources']):
 if Path(src).name=='IntegrationUtilities.cpp':
  obj=out/'IntegrationUtilities.obj';run(['cl','/c',*c['flags'],'/Fo'+str(obj),src])
 else:obj=work/'baseline'/f'unit-{i:03d}.obj'
 objects.append(str(obj))
rsp=out/'link.rsp';rsp.write_text('\n'.join(['"/DLL"','"/MACHINE:X64"','"/INCREMENTAL:NO"','"/OUT:'+str(out/'OpenMMCUDA.dll')+'"','"/IMPLIB:'+str(out/'OpenMMCUDA.lib')+'"']+['"'+p+'"' for p in objects+c['libs']]))
run(['link','/nologo','@'+str(rsp)]);print('Candidate plugin rebuilt from one changed translation unit',flush=True)
