#!/usr/bin/env python3
"""Stop only GeneTrace project jobs, retaining resumable results."""
import datetime,json,os,pathlib,signal,subprocess,time
ROOT=pathlib.Path(__file__).resolve().parents[1]
NAMES={'run_scalable_orthology.py','run_ecoli100_orthology.py','finalize_ecoli_followup.py','validate_ecoli_candidate_reads.py','benchmark_fcs_ecoli_controls.py','scan_ecoli100_domains.py','collect_ecoli100_reads.py','cluster_ecoli100_proteins.py','build_ecoli100_context.py','audit_ecoli100_candidates.py','prepare_ecoli_quality100.py','finish_hn1_read_qc.py','download_verified_read.py','test_domain_bacteria.py','compare_locus_cohort.py'}
def processes():
 result={}
 for line in subprocess.check_output(['ps','-axo','pid,ppid,stat,args'],text=True).splitlines()[1:]:
  f=line.strip().split(None,3)
  if len(f)==4:result[int(f[0])]=dict(ppid=int(f[1]),stat=f[2],args=f[3])
 return result
rows=processes();own=os.getpid();roots=set()
for pid,row in rows.items():
 args=row['args'].split();exe=args[0] if args else ''
 if pid==own:continue
 if 'python' in pathlib.Path(exe).name.lower() and (any(pathlib.Path(a).name in NAMES for a in args[1:3]) or args[1:5] in [['-m','genetrace.web','--port','8765'],['-m','genetrace.demo','--port','8766']]):roots.add(pid)
 if str(ROOT/'tools') in exe and any(x in pathlib.Path(exe).name.lower() for x in ['orthofinder','diamond','mafft','fasttree','mmseqs','minimap2','samtools','hmmsearch']):roots.add(pid)
selected=set(roots)
while True:
 extra={pid for pid,r in rows.items() if r['ppid'] in selected and pid!=own}
 if extra<=selected:break
 selected|=extra
for pid,row in rows.items():
 args=row['args'].split()
 if args and pathlib.Path(args[0]).name=='caffeinate' and any(str(p) in args[1:] for p in selected):selected.add(pid)
stamp=datetime.datetime.now(datetime.timezone.utc).isoformat();tag=datetime.datetime.now().strftime('%Y%m%d_%H%M%S');snapshot=dict(request='user_requested_stop',stopped_utc=stamp,processes=[dict(pid=p,**rows[p]) for p in sorted(selected)],statuses={},interrupted_outputs=[])
status_paths=[ROOT/'outputs/orthology_20261030/orthofinder_status.json',ROOT/'outputs/ecoli100_orthology_20261102/status.json']
for path in status_paths:
 if path.exists():snapshot['statuses'][str(path.relative_to(ROOT))]=json.loads(path.read_text())
# Freeze supervisors before terminating workers so no next stage can be started.
for pid in selected:
 try:os.kill(pid,signal.SIGSTOP)
 except ProcessLookupError:pass
for pid in selected:
 try:os.kill(pid,signal.SIGTERM)
 except ProcessLookupError:pass
for pid in selected:
 try:os.kill(pid,signal.SIGCONT)
 except ProcessLookupError:pass
time.sleep(1)
remaining=processes()
for pid in selected:
 if pid in remaining and not remaining[pid]['stat'].startswith('Z'):
  try:os.kill(pid,signal.SIGKILL)
  except ProcessLookupError:pass
time.sleep(.5)
remaining=processes();live=[pid for pid in selected if pid in remaining and not remaining[pid]['stat'].startswith('Z')]
# Retain unfinished outputs separately; they must not masquerade as completed searches.
if not live:
 for pid in selected:
  args=rows[pid]['args'].split()
  if args and pathlib.Path(args[0]).name=='diamond' and '-o' in args:
   output=pathlib.Path(args[args.index('-o')+1]);output=pathlib.Path(str(output)+'.gz') if '--compress' in args else output
   if output.is_relative_to(ROOT) and output.exists():
    dest=pathlib.Path(str(output)+'.interrupted_'+tag);output.rename(dest);snapshot['interrupted_outputs'].append(dict(original=str(output),preserved=str(dest)))
 for path in status_paths:
  if path.exists():
   status=json.loads(path.read_text())
   if status.get('status') not in ['complete']:
    status.update(status='paused_by_user',stopped_utc=stamp,resume_required=True,note='Running workers stopped; completed files retained; interrupted outputs preserved separately.');path.write_text(json.dumps(status,indent=2)+'\n')
snapshot['remaining_live_processes']=live;snapshot['all_target_jobs_stopped']=not live
checkpoint=ROOT/'outputs'/('shutdown_checkpoint_'+tag+'.json');checkpoint.write_text(json.dumps(snapshot,indent=2)+'\n');print(json.dumps(dict(stopped_processes=len(selected),remaining_live_processes=live,checkpoint=str(checkpoint),interrupted_outputs_preserved=len(snapshot['interrupted_outputs'])),indent=2))
if live:raise SystemExit(1)
