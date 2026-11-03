"""Read-only, bounded monitoring for the two known OrthoFinder workflows."""
import datetime,json,pathlib,re,shutil,subprocess,threading,time
from .engine import ROOT
LOCK=threading.Lock();CACHE=None;CACHED_AT=0
JOBS=[('yeast','효모 전체 유전자 OrthoFinder','outputs/orthology_20261030/orthofinder_status.json','outputs/orthology_20261030/core_inference.log','runs/orthofinder32_core_20261030'),('ecoli','대장균 전체 유전자 OrthoFinder','outputs/ecoli100_orthology_20261102/status.json','outputs/ecoli100_orthology_20261102/core_inference.log','runs/ecoli100_orthofinder_core_20261102')]
def read_json(path):
 try:return json.loads(path.read_text()),None
 except FileNotFoundError:return {},'not_available'
 except (json.JSONDecodeError,OSError):return {},'temporarily_unreadable'
def age_seconds(stamp,now):
 try:return max(0,(now-datetime.datetime.fromisoformat(stamp.replace('Z','+00:00'))).total_seconds())
 except (ValueError,TypeError,AttributeError):return None

def job_snapshot(config,processes,process_error,now):
 key,title,status_path,log_path,core_path=config;s,error=read_json(ROOT/status_path);pid=s.get('pid');process=processes.get(pid);expected=str(ROOT/core_path)
 valid=bool(process and expected in process['args'] and 'orthofinder' in process['args'] and not process['stat'].startswith('Z')) if not process_error else None
 children=set([pid]) if valid else set()
 for _ in range(8):children.update(p for p,r in processes.items() if r['ppid'] in children)
 active=[processes[p] for p in children if p in processes and not processes[p]['stat'].startswith('Z')];searching=[r for r in active if pathlib.Path(r['args'].split()[0]).name=='diamond'];age=age_seconds(s.get('updated_utc'),now)
 state=s.get('status','not_available')
 if state=='running':observed='running' if valid else 'unknown' if process_error else 'not_running'
 elif state in ['paused_by_user','resource_limit']:observed='running_despite_status' if valid else state
 else:observed=state
 progress=None
 if s.get('stage')=='core_inference':
  directories=sorted((ROOT/core_path).glob('Results_*/WorkingDirectory'));directory=directories[-1] if directories else None
  if directory:
   files=set(directory.glob('Blast*.txt.gz'));inflight=set()
   for row in searching:
    args=row['args'].split()
    if '-o' in args:
     out=pathlib.Path(args[args.index('-o')+1]);inflight.add(pathlib.Path(str(out)+'.gz'))
   total=s.get('core_size',0)**2
   if total:progress=dict(result_files=len(files-inflight),total=total,active_searches=len(searching),percent=round(min(100,len(files-inflight)/total*100),1),scope='core_search_result_files',validated_complete=False,note='실행 중인 출력은 제외한 검색 결과 파일 수입니다. 전체 분석 완료율이나 최종 검증 결과가 아닙니다.')
 log=ROOT/log_path
 reported_log=pathlib.Path(s.get('log') or log)
 if reported_log.parent.resolve()==(ROOT/status_path).parent.resolve() and reported_log.name in ['core_inference.log','all_genome_assignment.log']:log=reported_log
 try:
  with log.open('rb') as f:
   f.seek(0,2);end=f.tell();f.seek(max(0,end-12000));lines=f.read().decode(errors='replace').splitlines()
  # Only selected progress lines, never arbitrary commands or environment values.
  log_lines=[line.strip() for line in lines if re.search(r'(\d{4}-\d\d-\d\d \d\d:\d\d|Using \d+ thread|Running \d+ of|ERROR|Error:)',line)][-4:]
 except OSError:log_lines=[]
 return dict(id=key,title=title,status=state,observed_status=observed,status_read_error=error,stage=s.get('stage'),pid=pid,process_alive=valid,process_check_error=process_error,active_processes=len(active) if valid else None,cpu_percent=round(sum(r['cpu'] for r in active),1) if valid else None,configured_search_threads=s.get('threads'),assemblies=s.get('retained_assemblies'),coding_units=s.get('coding_units'),updated_utc=s.get('updated_utc'),heartbeat_age_seconds=round(age,1) if age is not None else None,heartbeat_stale=age is None or age>60,progress=progress,log_lines=log_lines,overall_percent=None)

def snapshot():
 global CACHE,CACHED_AT
 with LOCK:
  if CACHE is not None and time.monotonic()-CACHED_AT<5:return CACHE
  now=datetime.datetime.now(datetime.timezone.utc);processes={};process_error=None
  try:
   output=subprocess.check_output(['ps','-axo','pid,ppid,pcpu,stat,args'],text=True,timeout=3,stderr=subprocess.DEVNULL)
   for line in output.splitlines()[1:]:
    f=line.strip().split(None,4)
    if len(f)==5:processes[int(f[0])]=dict(ppid=int(f[1]),cpu=float(f[2]),stat=f[3],args=f[4])
  except (OSError,subprocess.SubprocessError,ValueError):process_error='process_check_unavailable'
  completed=[]
  for title,path in [('대장균 100개 Pfam 패널','outputs/ecoli100_domains_20261102/status.json'),('대장균 후보 raw read 검사','outputs/ecoli100_raw_validation_20261102/status.json'),('NCBI 오염 양성 3개 대조','outputs/contamination_validation_20261102/status.json')]:
   state,error=read_json(ROOT/path);completed.append(dict(title=title,stage=state.get('stage','not_available'),read_error=error))
  CACHE=dict(checked_utc=now.isoformat(),refresh_seconds=15,jobs=[job_snapshot(c,processes,process_error,now) for c in JOBS],free_disk_gib=round(shutil.disk_usage(ROOT).free/1024**3,1),reserve_disk_gib=6,other_workflows=completed,note='로컬 Mac에서 실행됩니다. 컴퓨터가 꺼지면 계산도 중단됩니다. 남은 전체 시간은 아직 산출하지 않습니다.');CACHED_AT=time.monotonic();return CACHE
