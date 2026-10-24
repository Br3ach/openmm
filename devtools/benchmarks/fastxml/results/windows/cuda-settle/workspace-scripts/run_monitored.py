from pathlib import Path
import ctypes,json,os,subprocess,time,threading
from ctypes import wintypes
base=Path.cwd();work=base/'settle-opt';dest=work/'results/resumed';dest.mkdir(exist_ok=True)
env=os.environ.copy();env['PATH']=r'C:\Program Files\NVIDIA GPU Computing Toolkit\CUDA\v13.4\bin\x64'+os.pathsep+env['PATH'];env['OPENMM_CACHE_DIR']=str(work/'cuda-cache')
class Times(ctypes.Structure):_fields_=[('low',wintypes.DWORD),('high',wintypes.DWORD)]
class Memory(ctypes.Structure):_fields_=[('cb',wintypes.DWORD),('faults',wintypes.DWORD)]+[(n,ctypes.c_size_t) for n in ['peakWS','ws','peakPaged','paged','peakNonPaged','nonPaged','pagefile','peakPagefile','private']]
getTimes=ctypes.windll.kernel32.GetProcessTimes;getTimes.argtypes=[wintypes.HANDLE]+[ctypes.POINTER(Times)]*4
getMem=ctypes.windll.psapi.GetProcessMemoryInfo;getMem.argtypes=[wintypes.HANDLE,ctypes.POINTER(Memory),wintypes.DWORD]
def run_case(variant,test,precision,timeout=600,extra=()):
 folder=work/variant;exe=folder/(test+'.exe');start=time.perf_counter();lines=[];samples=[]
 p=subprocess.Popen([str(exe),precision,*extra],cwd=folder,env=env,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True)
 print(json.dumps(dict(event='start',variant=variant,test=test,precision=precision,pid=p.pid)),flush=True)
 def reader():
  for line in p.stdout:lines.append(line);print(json.dumps(dict(event='output',variant=variant,test=test,line=line.rstrip())),flush=True)
 thread=threading.Thread(target=reader,daemon=True);thread.start();last=-10
 while p.poll() is None:
  elapsed=time.perf_counter()-start
  if elapsed>=timeout:p.kill();p.wait();timedout=True;break
  if elapsed-last>=10:
   last=elapsed;t=[Times() for _ in range(4)];m=Memory();m.cb=ctypes.sizeof(m)
   ok1=getTimes(int(p._handle),*[ctypes.byref(x) for x in t]);ok2=getMem(int(p._handle),ctypes.byref(m),ctypes.sizeof(m))
   sample=dict(elapsed=round(elapsed,2),cpu_seconds=((t[2].high<<32)+t[2].low+(t[3].high<<32)+t[3].low)/1e7 if ok1 else None,working_set_mib=round(m.ws/2**20,2) if ok2 else None)
   samples.append(sample);print(json.dumps(dict(event='sample',variant=variant,test=test,**sample)),flush=True)
  time.sleep(0.25)
 else:timedout=False
 thread.join(timeout=5)
 row=dict(test=test,precision=precision,variant=variant,pid=p.pid,returncode=p.returncode,timeout=timedout,timeout_seconds=timeout,seconds=time.perf_counter()-start,output=''.join(lines),samples=samples)
 target=dest/(test+'-'+precision+'-'+variant+'.json');target.write_text(json.dumps(row,indent=2));print(json.dumps({k:v for k,v in row.items() if k!='samples'}),flush=True);return row
if __name__=='__main__':
 run_case('baseline','TestCudaNoseHooverIntegrator','mixed',600)
 run_case('candidate','TestCudaNoseHooverIntegrator','mixed',600)
