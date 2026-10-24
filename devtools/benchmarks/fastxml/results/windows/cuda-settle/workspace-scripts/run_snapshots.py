from pathlib import Path
import os,json,time,subprocess,uuid
from run_monitored import run_case,work,dest,env
rows=[]
def record(row):
 rows.append(row);(dest/'snapshot-comparisons-raw.json').write_text(json.dumps(rows,indent=2))
def get_fingerprint(row):
 if row['returncode'] or row['timeout']:raise RuntimeError('Snapshot execution failed: '+row['output'])
 return json.loads(row['output'].strip().splitlines()[-1])['fingerprint']
for precision in ['single','mixed','double']:
 paths={variant:dest/(precision+'-'+variant+'.chk') for variant in ['baseline','candidate']}
 saved={}
 for variant in ['baseline','candidate']:
  r=run_case(variant,'TestCudaSettleSnapshots',precision,180,('save',str(paths[variant])));r['mode']='save';record(r);saved[variant]=get_fingerprint(r)
 if len(set(saved.values()))!=1:raise RuntimeError('Fresh seeded trajectories differ: '+str(saved))
 for writer,reader in [('baseline','candidate'),('candidate','baseline')]:
  r=run_case(reader,'TestCudaSettleSnapshots',precision,180,('load',str(paths[writer])));r['mode']='load';r['checkpoint_writer']=writer;record(r)
  if get_fingerprint(r)!=saved[writer]:raise RuntimeError('Cross-plugin checkpoint continuation differs')
summary=dict(precisions=3,processes=len(rows),seeded_trajectories_identical=True,cross_plugin_checkpoint_directions=2,all_replays_identical=True)
(dest/'snapshot-comparisons-summary.json').write_text(json.dumps(summary,indent=2));print(summary,flush=True)
