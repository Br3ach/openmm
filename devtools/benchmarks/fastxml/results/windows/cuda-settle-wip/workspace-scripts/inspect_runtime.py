from pathlib import Path
import struct,json
base=Path.cwd();sdk=Path(r'C:\Program Files\NVIDIA GPU Computing Toolkit\CUDA\v13.4\bin\x64')
def deps(p):
 b=p.read_bytes();pe=struct.unpack_from('<I',b,0x3c)[0];n=struct.unpack_from('<H',b,pe+6)[0];optsize=struct.unpack_from('<H',b,pe+20)[0];opt=pe+24;sections=[]
 for i in range(n):
  o=opt+optsize+40*i;size,rva,raw,ptr=struct.unpack_from('<IIII',b,o+8);sections.append((rva,max(size,raw),ptr))
 def at(rva):
  for start,size,ptr in sections:
   if start<=rva<start+size:return ptr+rva-start
  raise ValueError(rva)
 rva,size=struct.unpack_from('<II',b,opt+120);result=[]
 if rva:
  o=at(rva)
  while any(b[o:o+20]):
   nrva=struct.unpack_from('<I',b,o+12)[0];s=at(nrva);result.append(b[s:b.index(b'\0',s)].decode());o+=20
 return result
files=[base/'settle-opt/candidate/OpenMMCUDA.dll',base/'settle-opt/candidate/OpenMM.dll',sdk/'nvrtc-builtins64_134.dll'];seen={}
for p in files:
 if p.name in seen:continue
 seen[p.name]=dict(path=str(p),bytes=p.stat().st_size,imports=deps(p))
 for name in seen[p.name]['imports']:
  if (sdk/name).exists() and name not in seen:files.append(sdk/name)
Path('settle-opt/runtime-dependencies.json').write_text(json.dumps(seen,indent=2));print(json.dumps(seen,indent=2))
