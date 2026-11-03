"""Small reproducible evidence-review demo; snapshots + optional live Meta planning.
Run: python3 -m genetrace.demo --port 8766. No original genomes/tools required.
"""
import argparse,datetime,json,re,threading,uuid
from http.server import ThreadingHTTPServer
from urllib.parse import urlsplit
from .agent import ROOT,model_config,request_model
from .io import sha256
from .web import Handler
DEMO=ROOT/'demo';RUNS=ROOT/'runs/submission_demo';LOCK=threading.Lock()
def load_case(key):
 manifest=json.loads((DEMO/'manifest.json').read_text());entry=next((c for c in manifest['cases'] if c['id']==key),None)
 if entry is None:raise ValueError('Unknown case')
 path=DEMO/'evidence'/(key+'.json')
 if sha256(path)!=entry['sha256']:raise ValueError('Evidence snapshot hash mismatch')
 return json.loads(path.read_text())
def persist(report):
 RUNS.mkdir(parents=True,exist_ok=True);path=RUNS/(report['id']+'.json');temp=path.with_suffix('.tmp');temp.write_text(json.dumps(report,ensure_ascii=False,indent=2)+'\n');temp.replace(path)
def investigate(report):
 try:
  case=load_case(report['case_id'])
  if report.get('pipeline')=='execute':
   from .execution import investigate as execute_investigation
   execute_investigation(report,case,persist);return
  cfg=model_config();cfg['response_language']='en';report.update(status='running',case_title=case['title'],model_id=cfg['model'] if report['mode']=='meta' else None);persist(report)
  remaining={t['tool']:t for t in case['tools']}
  def inspect(name):
   item=remaining.pop(name);data=json.loads(json.dumps(item['data']))
   if isinstance(data.get('read_provenance'),dict):data['read_provenance'].pop('note',None)
   report['evidence'].append(dict(id=f"E{len(report['evidence'])+1:03d}",tool=name,title=item['title'],result=data,sources=item['sources'],acquisition='lookup_completed_analysis_snapshot'));report['trace'].append(dict(kind='evidence_lookup',tool=name));persist(report)
  inspect('lookup_cohort_comparison')
  while remaining:
   if report['mode']=='meta':
    step,usage=request_model(dict(case=case['title'],question=case['question'],evidence_mode=case['evidence_mode'],instruction='Read the available evidence; provenance concerns and neighboring controls must be disclosed. Do not claim new genome analysis or a confirmed causal diagnosis.'),list(remaining),report['evidence'],cfg);report['model_calls']+=1;report['trace'].append(dict(kind='live_meta_planner',step=step,usage=usage));choice=step['next_tool']
    if choice is None:choice=next(iter(remaining));report['trace'].append(dict(kind='guard',reason='Read all case evidence before concluding',tool=choice))
   else:choice=next(iter(remaining));report['trace'].append(dict(kind='baseline_fixed_order',tool=choice))
   inspect(choice)
  if report['mode']=='meta':
   step,usage=request_model(dict(question=case['question'],instruction='All available case evidence has been inspected. Summarize provisional interpretations and limitations, with evidence IDs. next_tool must be null.'),[],report['evidence'],cfg);report['model_calls']+=1;report['hypotheses']=step['hypotheses'];report['trace'].append(dict(kind='live_meta_conclusion',step=step,usage=usage))
  report.update(status='complete',next_check=case['next_check'],analysis_execution='Saved analysis lookup only; live model inference only when mode=meta and model_calls>0',causal_diagnosis='unresolved',biological_loss_inferred=False)
 except Exception as e:report.update(status='failed',error='Investigation failed during evidence acquisition, analysis execution or Meta response validation. No baseline fallback was performed.',error_type=type(e).__name__)
 finally:
  persist(report);LOCK.release()
class DemoHandler(Handler):
 def do_GET(self):
  if not self.valid_host():return self.send(403,{'error':'Local host required'})
  path=urlsplit(self.path).path
  try:
   if path=='/api/demo/status':
    from .execution import availability
    c=model_config();return self.send(200,dict(execution_available=availability(),meta_configured=bool(c['key'] and c['model']),model=c['model'],busy=LOCK.locked(),snapshot_mode=True))
   if path=='/api/demo/cases':return self.send(200,json.loads((DEMO/'manifest.json').read_text()))
   m=re.fullmatch(r'/api/demo/jobs/([a-f0-9]{32})',path)
   if m:
    p=RUNS/(m[1]+'.json');return self.send(200,json.loads(p.read_text())) if p.exists() else self.send(404,{'error':'Unknown job'})
   files={'/':'index.html','/demo.js':'demo.js','/style.css':'style.css'}
   if path in files:
    p=DEMO/'web'/files[path];kind='text/html' if path=='/' else 'text/javascript' if path.endswith('.js') else 'text/css';return self.send(200,p.read_bytes(),kind+'; charset=utf-8')
   return self.send(404,{'error':'Not found'})
  except Exception:return self.send(500,{'error':'Demo evidence unavailable or invalid'})
 def do_POST(self):
  if not self.valid_host():return self.send(403,{'error':'Local host required'})
  origin=self.headers.get('Origin');allowed={f'http://127.0.0.1:{self.server.server_port}',f'http://localhost:{self.server.server_port}'}
  if origin and origin not in allowed:return self.send(403,{'error':'Origin not allowed'})
  if self.path!='/api/demo/investigate':return self.send(404,{'error':'Not found'})
  if self.headers.get('Content-Type','').split(';')[0]!='application/json':return self.send(415,{'error':'JSON required'})
  try:
   n=int(self.headers.get('Content-Length','0'))
   if not 0<n<=4096:return self.send(413,{'error':'Body too large'})
   body=json.loads(self.rfile.read(n));case=load_case(body.get('case_id'));mode=body.get('mode','baseline')
   if mode not in ['baseline','meta']:raise ValueError('Invalid mode')
   pipeline=body.get('pipeline','snapshot')
   if pipeline not in ['snapshot','execute']:raise ValueError('Invalid pipeline')
   if pipeline=='execute':
    from .execution import availability
    if case['id']!='nimt':return self.send(400,{'error':'Fresh execution is available for nimT only.'})
    if not availability():return self.send(409,{'error':'Install minimap2 and samtools and verify the execution case inputs.'})
   cfg=model_config()
   if mode=='meta' and not(cfg['key'] and cfg['model']):return self.send(409,{'error':'Set MODEL_API_KEY and META_MODEL in local .env; never paste keys into chat.'})
   if not LOCK.acquire(blocking=False):return self.send(409,{'error':'One investigation is already running'})
   report=dict(id=uuid.uuid4().hex,case_id=case['id'],mode=mode,pipeline=pipeline,status='queued',created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),model_calls=0,evidence=[],hypotheses=[],trace=[])
   try:persist(report);threading.Thread(target=investigate,args=(report,),daemon=True).start()
   except Exception:LOCK.release();raise
   return self.send(202,{'id':report['id']})
  except (ValueError,TypeError,AttributeError):return self.send(400,{'error':'Invalid case or mode'})
  except Exception:return self.send(500,{'error':'Could not start investigation'})
def main():
 p=argparse.ArgumentParser();p.add_argument('--port',type=int,default=8766);a=p.parse_args();server=ThreadingHTTPServer(('127.0.0.1',a.port),DemoHandler);print(f'GeneTrace submission demo: http://127.0.0.1:{a.port}/',flush=True)
 try:server.serve_forever()
 except KeyboardInterrupt:server.server_close()
if __name__=='__main__':main()
