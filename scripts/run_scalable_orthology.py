"""Run official core+assign OrthoFinder, retaining every genome and coding unit.

The supervisor continues from core inference to full-cohort assignment and final
accounting without another user action. Status is updated during long stages.
"""
import argparse
import csv
import json
import os
import shutil
import signal
import subprocess
import sys
import time
from datetime import datetime,timezone
from pathlib import Path

sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from genetrace.io import sha256,write_json

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'outputs/orthology_20261030'
SOURCE=ROOT/'data/orthology'
ENV=ROOT/'tools/orthofinder32-env'
DEPS=ROOT/'tools/orthofinder-env'
STATUS=OUT/'orthofinder_status.json'
STATE=OUT/'scalable_run_state.json'
CORE_BASE=ROOT/'runs/orthofinder32_core_20261030'


def stamp(): return datetime.now(timezone.utc).isoformat()


def run_stage(stage,command,env,state,result_parent,threads):
    with (OUT/(stage+'.log')).open('a',buffering=1) as log:
        log.write('\nCOMMAND '+json.dumps(command)+'\n')
        process=subprocess.Popen(command,stdout=log,stderr=subprocess.STDOUT,cwd=ROOT,env=env,start_new_session=True)
        resource_stop=False
        print(stage+' started; PID '+str(process.pid),flush=True)
        while process.poll() is None:
            free=shutil.disk_usage(ROOT).free
            status=dict(status='running',stage=stage,pid=process.pid,supervisor_pid=os.getpid(),
                        updated_utc=stamp(),free_disk_gib=free/1024**3,threads=threads,
                        retained_assemblies=state['retained_assemblies'],
                        coding_units=state['coding_units'],core_size=state['core_size'],
                        additional_size=state['additional_size'],log=str(OUT/(stage+'.log')),
                        all_retained_genomes_included=True,gene_subsampling=False)
            result_dirs=sorted(result_parent.glob('Results_*')) if result_parent.exists() else []
            status['result_directories']=[str(p) for p in result_dirs if p.is_dir()]
            write_json(STATUS,status)
            if free<6*1024**3:
                resource_stop=True
                os.killpg(process.pid,signal.SIGINT)
                break
            time.sleep(10)
        code=process.wait()
    if code!=0 or resource_stop:
        write_json(STATUS,dict(status='resource_limit' if resource_stop else 'failed',stage=stage,
                              exit_code=code,updated_utc=stamp(),log=str(OUT/(stage+'.log')),
                              checkpoints_retained=True))
        raise RuntimeError(stage+' failed; see stage log and resume checkpoints')


def final_result(directory,expected_columns,require_tree=False):
    candidates=[]
    for result in directory.glob('Results_*'):
        table=result/'Orthogroups/Orthogroups.tsv'
        if not table.exists(): continue
        with table.open() as handle:
            header=next(csv.reader(handle,delimiter='\t'))
        if set(header[1:])!=set(expected_columns): continue
        if require_tree and not (result/'Species_Tree/SpeciesTree_rooted.txt').exists(): continue
        candidates.append(result)
    if not candidates:
        raise RuntimeError('No complete result has the expected full cohort and required tree')
    return max(candidates,key=lambda p:p.stat().st_mtime)


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--threads',type=int,default=6)
    parser.add_argument('--analysis-threads',type=int,default=1)
    parser.add_argument('--resume',action='store_true')
    args=parser.parse_args()
    split=json.loads((OUT/'core_selection.json').read_text())
    input_summary=json.loads((OUT/'all_gene_input_summary.json').read_text())
    selected=json.loads((OUT/'selected_assemblies.json').read_text())
    all_ids={r['assembly_accession'] for r in selected}
    if set(split['core'])&set(split['additional']) or set(split['core'])|set(split['additional'])!=all_ids:
        raise ValueError('Core+additional split does not cover the full selected cohort exactly once')
    if {r['assembly_accession'] for r in input_summary}!=all_ids:
        raise ValueError('Missing prepared coding-locus inputs')
    for row in input_summary:
        if sha256(SOURCE/'proteomes'/(row['assembly_accession']+'.faa'))!=row['proteome_sha256']:
            raise ValueError('Prepared proteome modified after manifest creation')
    if STATE.exists():
        if not args.resume: raise FileExistsError('Run state already exists; use --resume')
        state=json.loads(STATE.read_text())
    else:
        if CORE_BASE.exists(): raise FileExistsError('Unexpected existing core run directory')
        state=dict(started_utc=stamp(),orthofinder_version='3.2',retained_assemblies=len(selected),
                   coding_units=sum(r['orthology_coding_units'] for r in input_summary),
                   source_protein_records=sum(r['source_protein_records'] for r in input_summary),
                   core_size=len(split['core']),additional_size=len(split['additional']),
                   all_genes_included=True,all_retained_genomes_included=True,commands=[])
        write_json(STATE,state)
    env=dict(os.environ)
    env.update(PATH=str(ENV/'bin')+os.pathsep+str(DEPS/'bin')+os.pathsep+env.get('PATH',''),
               OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1')
    base=[str(ENV/'bin/python'),str(ENV/'bin/orthofinder')]
    options=['-X','-S','diamond','-t',str(args.threads),'-a',str(args.analysis_threads)]
    if 'core_results' not in state:
        if CORE_BASE.exists():
            existing=sorted(CORE_BASE.glob('Results_*'))
            if not existing: raise RuntimeError('Core checkpoint directory missing')
            command=base+['-b',str(existing[-1])]+options
        else:
            command=base+['-f',str(SOURCE/'core_proteomes'),'-o',str(CORE_BASE)]+options
        state['commands'].append(command);write_json(STATE,state)
        run_stage('core_inference',command,env,state,CORE_BASE,args.threads)
        core=final_result(CORE_BASE,split['core'],require_tree=True)
        state['core_results']=str(core);write_json(STATE,state)
    core=Path(state['core_results'])
    if 'full_cohort_results' not in state:
        previous=[Path(p) for p in state.get('assignment_result_directories',[]) if Path(p).exists()]
        if previous:
            command=base+['-b',str(previous[-1])]+options
        else:
            command=base+['--assign',str(SOURCE/'additional_proteomes'),'--core',str(core)]+options
        # Completed core folder is retained. --assign creates a sibling result.
        state['commands'].append(command);write_json(STATE,state)
        try:
            run_stage('all_genome_assignment',command,env,state,core.parent,args.threads)
        finally:
            state['assignment_result_directories']=[str(p) for p in core.parent.glob('Results_*') if p!=core]
            write_json(STATE,state)
        complete=final_result(core.parent,all_ids)
        state['full_cohort_results']=str(complete);write_json(STATE,state)
    write_json(STATUS,dict(status='validating_all_gene_accounting',updated_utc=stamp(),retained_assemblies=len(all_ids)))
    subprocess.run([sys.executable,str(ROOT/'scripts/summarize_orthogroups.py'),
                    '--results',state['full_cohort_results']],check=True)
    state['completed_utc']=stamp();write_json(STATE,state)
    write_json(STATUS,dict(status='complete',updated_utc=stamp(),retained_assemblies=len(all_ids),
                          coding_units=state['coding_units'],results=state['full_cohort_results'],
                          output=str(OUT/'orthogroups')))
    print('All retained genomes and coding loci accounted for in final orthogroups.',flush=True)


if __name__=='__main__':
    main()
