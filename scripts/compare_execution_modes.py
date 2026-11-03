"""Paired case study: fixed policy vs live Meta planning on identical public inputs."""
import datetime,json,pathlib,sys,uuid
ROOT=pathlib.Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from genetrace import demo
from genetrace.io import sha256
OUT=ROOT/'outputs/execution_comparison_20261103';OUT.mkdir(exist_ok=True)
reports=[]
for mode in ['baseline','meta']:
 if '--reuse' in sys.argv:
  reports.append(json.loads((OUT/(mode+'.json')).read_text()));continue
 r=dict(id=uuid.uuid4().hex,case_id='nimt',mode=mode,pipeline='execute',created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),model_calls=0,evidence=[],hypotheses=[],trace=[])
 demo.LOCK.acquire();demo.investigate(r);reports.append(r)
 print(mode,r['status'],'model_calls',r['model_calls'],'seconds',r.get('elapsed_seconds'),flush=True)
 if r['status']!='complete':raise RuntimeError('Investigation failed: '+str(r.get('error_type')))
 (OUT/(mode+'.json')).write_text(json.dumps(r,indent=2)+'\n')

def measures(r):
 c=next(e['result']['candidate'] for e in r['evidence'] if e['tool']=='align_candidate_reads')
 n=next(e['result']['neighbor_regions'] for e in r['evidence'] if e['tool']=='check_neighbor_regions')
 n=[dict(x,name=x.get('name',name)) for x,name in zip(n,['nimR','yeaO'])]
 return dict(candidate_metrics=c['metrics'],candidate_signals=c['signals'],neighbor_metrics={x['name']:x['metrics'] for x in n},neighbor_signals={x['name']:x['signals'] for x in n})
rows=[]
for r in reports:
 known={e['id'] for e in r['evidence']}
 rows.append(dict(mode=r['mode'],report_id=r['id'],status=r['status'],elapsed_seconds=r['elapsed_seconds'],model_calls=r['model_calls'],fresh_tools=[e['tool'] for e in r['evidence'] if e['acquisition']=='fresh_local_execution'],tool_order=[e['tool'] for e in r['evidence']],cited_hypotheses=len(r['hypotheses']),citation_membership_valid=all(h['evidence_ids'] and set(h['evidence_ids'])<=known for h in r['hypotheses']),guard_interventions=sum(t['kind']=='guard' for t in r['trace']),judgment=r['judgment']))
result=dict(input_manifest_sha256=sha256(ROOT/'demo/execution/nimt/manifest.json'),runs=rows,deterministic_measurements_equal=measures(reports[0])==measures(reports[1]),measurements=measures(reports[1]),scope='One curated case, one run per mode; fixed minimum-check policy shared by both modes. Not an accuracy, productivity, independence or causal diagnosis benchmark.',interpretation='Same tools generate the evidence. Meta supplies evidence-dependent selection rationales and cited synthesis. Different tool order is not proof of better decisions; equal guarded outcomes are expected under this shared policy.')
(OUT/'comparison.json').write_text(json.dumps(result,indent=2)+'\n')
(ROOT/'demo/execution_comparison.json').write_text(json.dumps(result,indent=2)+'\n')
(ROOT/'demo/verified_execution_report.json').write_text(json.dumps(reports[1],indent=2)+'\n')
lines=['# Fixed policy vs Meta-guided execution','','Same 756 reference-selected read pairs, assembly and annotation. Population context is a saved analysis.','', '| Mode | New tools executed | Meta responses | Elapsed (s) | Guard interventions | Cited interpretations |','|---|---:|---:|---:|---:|---:|']
for r in rows:lines.append(f"| {r['mode']} | {len(r['fresh_tools'])} | {r['model_calls']} | {r['elapsed_seconds']} | {r['guard_interventions']} | {r['cited_hypotheses']} |")
lines+=['',f"Deterministic observations equal: **{result['deterministic_measurements_equal']}**.",'',result['scope'],'',result['interpretation'],'','The structured judgment is a deterministic guard, not solely model-generated. Both modes can identify the neighbor discrepancy under the shared minimum-check policy. Claims of improved accuracy, faster expert judgment, or unique AI discovery are unsupported.']
(ROOT/'docs/submission/EXECUTION_COMPARISON.md').write_text('\n'.join(lines)+'\n')
print('comparison complete; equal measurements:',result['deterministic_measurements_equal'],flush=True)
