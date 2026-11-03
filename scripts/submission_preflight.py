#!/usr/bin/env python3
"""No-network readiness report: credentials are never printed; mocks do not satisfy live usage."""
import json,pathlib,sys
ROOT=pathlib.Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from genetrace.agent import model_config
from genetrace.demo import load_case

def main():
 settings=json.loads((ROOT/'docs/submission/submission.json').read_text());config=model_config();checks=[]
 for key in ['nimt','provenance','contamination']:load_case(key)
 checks.append(dict(check='bundled_real_evidence_snapshots',ok=True))
 checks.append(dict(check='Meta_credentials_configured',ok=bool(config['key'] and config['model']),note='Presence only; does not prove live API success.'))
 live=[]
 for p in (ROOT/'runs/submission_demo').glob('*.json'):
  r=json.loads(p.read_text())
  if r.get('scientific_review',{}).get('status')!='needs_correction' and r.get('mode')=='meta' and r.get('status')=='complete' and r.get('model_calls',0)>0 and any(t.get('kind')=='live_meta_conclusion' for t in r.get('trace',[])):live.append(dict(report=str(p.relative_to(ROOT)),model=r.get('model_id'),calls=r['model_calls']))
 checks.append(dict(check='successful_live_Meta_demo',ok=bool(live),reports=live))
 for field in ['repository_url','video_url']:checks.append(dict(check=field+'_provided',ok=settings[field].startswith('https://'),note='URL presence only; public accessibility must be checked separately.'))
 checks.append(dict(check='owner_license_review_complete',ok=settings['owner_license_review_complete']))
 out=ROOT/'dist';out.mkdir(exist_ok=True);report=dict(ready_for_manual_submission_review=all(c['ok'] for c in checks),checks=checks);(out/'preflight.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2));return 0 if report['ready_for_manual_submission_review'] else 2
if __name__=='__main__':raise SystemExit(main())
