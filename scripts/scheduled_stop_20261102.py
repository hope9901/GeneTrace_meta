#!/usr/bin/env python3
"""One-shot user-authorized stop at 17:20 KST; no automatic analysis restart."""
import datetime,json,os,pathlib,subprocess,sys,time
from zoneinfo import ZoneInfo
ROOT=pathlib.Path(__file__).resolve().parents[1];STATUS=ROOT/'outputs/scheduled_stop_20261102_1720.json';TARGET=datetime.datetime(2026,11,2,17,20,tzinfo=ZoneInfo('Asia/Seoul'))
def save(stage,**kw):
 p=STATUS.with_suffix('.tmp');p.write_text(json.dumps(dict(stage=stage,target=TARGET.isoformat(),pid=os.getpid(),updated_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),**kw),indent=2)+'\n');p.replace(STATUS)
def main():
 save('scheduled')
 while time.time()<TARGET.timestamp():
  if (ROOT/'outputs/cancel_scheduled_stop_20261102_1720').exists():save('cancelled');return
  time.sleep(min(5,max(.1,TARGET.timestamp()-time.time())))
 save('stopping')
 p=subprocess.run([sys.executable,str(ROOT/'scripts/stop_project_jobs.py')],cwd=ROOT,capture_output=True,text=True)
 (ROOT/'outputs/scheduled_stop_20261102_1720.log').write_text(p.stdout+p.stderr)
 save('complete' if p.returncode==0 else 'failed',exit_code=p.returncode,result=p.stdout[-6000:],error=p.stderr[-2000:])
if __name__=='__main__':main()
