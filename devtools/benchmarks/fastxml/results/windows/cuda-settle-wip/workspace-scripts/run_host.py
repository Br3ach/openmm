from pathlib import Path
import subprocess,json,statistics,hashlib
base=Path(__file__).resolve().parent;exe=base/'candidate/BenchmarkSettleStorage.exe';fixture=base.parent/'settle-assessment/fixture-constraints.bin';dest=base/'results/host';dest.mkdir(exist_ok=True);rows=[]
for repeat in range(3):
 order=['baseline','compact'] if repeat%2==0 else ['compact','baseline']
 for variant in order:
  p=subprocess.run([str(exe),variant,str(fixture)],cwd=base,capture_output=True,text=True,timeout=180)
  if p.returncode:raise RuntimeError(p.stdout+p.stderr)
  row=json.loads(p.stdout);row['repeat']=repeat+1;rows.append(row)
  (dest/'raw.json').write_text(json.dumps(rows,indent=2));print(json.dumps(row),flush=True)
assert len({(r['output_fingerprint'],r['neighbors_fingerprint'],r['topological_triangles'],r['settle_clusters'],r['is_shake_true']) for r in rows})==1
layout=json.loads(subprocess.check_output([str(exe),'layout'],text=True));summary=dict(layout=layout,particles=rows[0]['particles'],constraints=rows[0]['constraints'],active_constraints=rows[0]['active_constraints'],eligible_neighbor_entries=rows[0]['eligible_neighbor_entries'],topological_triangles=rows[0]['topological_triangles'],settle_clusters=rows[0]['settle_clusters'],is_shake_true=rows[0]['is_shake_true'],outputs_match_all_runs=True,exe_sha256=hashlib.sha256(exe.read_bytes()).hexdigest(),fixture_binary_sha256=hashlib.sha256(fixture.read_bytes()).hexdigest())
for variant in ['baseline','compact']:
 selected=[r for r in rows if r['variant']==variant]
 summary[variant]={key:statistics.median(r[key] for r in selected) for key in ['build_classify_record_seconds','row_release_seconds']}
 summary[variant]['phases']={}
 for phase in [p['phase'] for p in selected[0]['phases']]:
  points=[next(p for p in r['phases'] if p['phase']==phase) for r in selected]
  summary[variant]['phases'][phase]={key:statistics.median(p[key] for p in points) for key in ['working_set_kib','peak_working_set_kib','private_kib']}
summary['baseline_requested_container_bytes']=summary['particles']*(layout['map_object_bytes']+layout['map_sentinel_bytes'])+summary['eligible_neighbor_entries']*layout['map_entry_bytes']
summary['compact_requested_container_bytes']=summary['particles']*layout['compact_row_bytes']
(dest/'summary.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary,indent=2),flush=True)
