from pathlib import Path
import xml.etree.ElementTree as E,subprocess,concurrent.futures,json,shutil,hashlib
base=Path.cwd();work=base/'settle-opt';root=base/'source';out=work/'baseline';ns={'m':'http://schemas.microsoft.com/developer/msbuild/2003'}
t=E.parse(work/'cmake/platforms/cuda/sharedTarget/OpenMMCUDA.vcxproj');g=next(x for x in t.findall('m:ItemDefinitionGroup',ns) if "'Release|x64'" in x.get('Condition',''));c=g.find('m:ClCompile',ns)
inc=[x for x in c.find('m:AdditionalIncludeDirectories',ns).text.split(';') if not x.startswith('%(')];inc += [r'C:\Program Files\NVIDIA GPU Computing Toolkit\CUDA\v13.4\include',r'C:\Program Files\NVIDIA GPU Computing Toolkit\CUDA\v13.4\include\cccl']
defs=[x for x in c.find('m:PreprocessorDefinitions',ns).text.split(';') if not x.startswith('%(') and not x.startswith('CMAKE_INTDIR')]
flags=['/nologo','/O2','/Ob2','/MD','/EHsc','/std:c++14','/W2','/wd4305','/wd4244']+['/I'+x for x in inc]+['/D'+x for x in defs]
libs=[str(base/'constraint-opt/candidate/OpenMM.lib')]+g.find('m:Link/m:AdditionalDependencies',ns).text.split(';')[1:]
sources=[str(Path(x.get('Include')).resolve()) for x in t.findall('m:ItemGroup/m:ClCompile',ns) if x.get('Include')]
config=dict(flags=flags,libs=libs,sources=sources);(work/'build-config.json').write_text(json.dumps(config,indent=2))
def run(args):
 p=subprocess.run(args,capture_output=True,text=True)
 if p.returncode:raise RuntimeError(p.stdout+p.stderr)
 return p

def compile_one(pair):
 i,src=pair;obj=out/f'unit-{i:03d}.obj'
 if obj.exists():return str(obj)
 run(['cl','/c',*flags,'/Fo'+str(obj),src]);print('Compiled',Path(src).name,flush=True);return str(obj)
with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:objects=list(pool.map(compile_one,enumerate(sources)))
rsp=out/'link.rsp';rsp.write_text('\n'.join(['"/DLL"','"/MACHINE:X64"','"/INCREMENTAL:NO"','"/OUT:'+str(out/'OpenMMCUDA.dll')+'"','"/IMPLIB:'+str(out/'OpenMMCUDA.lib')+'"']+['"'+p+'"' for p in objects+libs]))
run(['link','/nologo','@'+str(rsp)])
for variant in ['baseline','candidate']:
 for name in ['OpenMM.dll','OpenMM.lib']:shutil.copy2(base/'constraint-opt/candidate'/name,work/variant/name)
print('Baseline built',hashlib.sha256((out/'OpenMMCUDA.dll').read_bytes()).hexdigest(),flush=True)
