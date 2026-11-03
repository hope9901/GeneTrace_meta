"""Loopback-only UI for local GeneTrace evidence investigations."""
import argparse,json,mimetypes,re,threading
from http.server import BaseHTTPRequestHandler,ThreadingHTTPServer
from urllib.parse import urlsplit,parse_qs
from .agent import ROOT,cohort,genes,load,model_config,Investigation,reference_for

WEB=ROOT/'web';BUSY=threading.Lock()

class Handler(BaseHTTPRequestHandler):
    def log_message(self,fmt,*args):pass
    def send(self,code,body,kind='application/json; charset=utf-8',filename=None):
        if isinstance(body,(dict,list)):body=json.dumps(body,ensure_ascii=False).encode()
        elif isinstance(body,str):body=body.encode()
        self.send_response(code);self.send_header('Content-Type',kind);self.send_header('Content-Length',str(len(body)))
        self.send_header('Cache-Control','no-store');self.send_header('X-Content-Type-Options','nosniff')
        self.send_header('Content-Security-Policy',"default-src 'self'; style-src 'self'; script-src 'self'; img-src 'self' data:; frame-ancestors 'none'; base-uri 'none'; form-action 'self'")
        if filename:self.send_header('Content-Disposition','attachment; filename="'+filename+'"')
        self.end_headers();self.wfile.write(body)
    def valid_host(self):
        return self.headers.get('Host','') in [f'127.0.0.1:{self.server.server_port}',f'localhost:{self.server.server_port}']
    def do_GET(self):
        if not self.valid_host():return self.send(403,dict(error='Local host required'))
        url=urlsplit(self.path);path=url.path
        try:
            if path=='/api/status':
                cfg=model_config()
                return self.send(200,dict(assemblies=len(cohort()),reference_genes=len(genes()),reference_counts={r:len(genes(r)) for r in sorted({reference_for(acc) for acc in cohort()})},busy=BUSY.locked(),meta_configured=bool(cfg['key'] and cfg['model']),meta_model=cfg['model'] or None,orthology=load(ROOT/'outputs/orthology_20261030/orthofinder_status.json',{}),read_pipeline=load(ROOT/'outputs/read_qc_HN1_20261030/pipeline_status.json',{})))
            if path=='/api/background-jobs':
                from .background import snapshot
                return self.send(200,snapshot())
            if path=='/api/assemblies':return self.send(200,[dict(**{k:r.get(k) for k in ['assembly_accession','resolved_strain','assembly_level','sequencing_group','number_of_proteins','raw_read_status']},species='Escherichia coli' if r.get('reference_accession') else 'Saccharomyces cerevisiae',reference_accession=reference_for(r['assembly_accession'])) for r in cohort().values()])
            if path=='/api/domain-results':
                from .domains import domain_snapshot
                query=parse_qs(url.query);acc=query.get('assembly',[''])[0];profile=query.get('profile',[None])[0]
                if acc not in cohort():return self.send(400,dict(error='Unknown assembly'))
                if profile is not None and not re.fullmatch(r'PF[0-9]{5}',profile):return self.send(400,dict(error='Invalid Pfam accession'))
                return self.send(200,domain_snapshot(acc,profile))
            if path=='/api/validation':
                directory=ROOT/'outputs/regional_validation_20261102'
                return self.send(200,dict(benchmark=load(directory/'benchmark_summary.json',{}),real_regions=load(directory/'real_region_summary.json',{}),causal_ambiguity=load(directory/'causal_ambiguity_demo.json',{})))
            if path=='/api/cohort-comparison':
                from .agent import resolve_gene
                query=parse_qs(url.query);acc=query.get('assembly',['GCA_903819125.2'])[0]
                if acc not in cohort():return self.send(400,dict(error='Unknown assembly'))
                gene=resolve_gene(query.get('gene',['YJU3'])[0],reference_for(acc));key=gene['locus_tag'] or gene['gene_name']
                directory=ROOT/'outputs/cohort_comparison_20261102'/key
                state=load(directory/'status.json',{'stage':'not_computed'});summary=load(directory/'summary.json',{})
                return self.send(200,dict(state=state,available_for_target=state.get('stage')=='complete' and summary.get('focal_assembly')==acc,summary=summary))
            if path=='/api/genes':
                query=parse_qs(url.query);q=query.get('q',[''])[0].upper()[:100];acc=query.get('assembly',['GCA_000146045.2'])[0]
                if acc not in cohort():return self.send(400,dict(error='Unknown assembly'))
                matches=[r for r in genes(reference_for(acc)) if q in (' '.join(r.get(k,'') for k in ['gene_name','locus_tag','gene_id','ncbi_protein_id'])).upper()]
                matches.sort(key=lambda r:(r['gene_name'].upper()!=q,r['gene_name']))
                return self.send(200,[{k:r[k] for k in ['gene_name','locus_tag','gene_id','ncbi_protein_id']} for r in matches[:40]])
            m=re.fullmatch(r'/api/jobs/([a-f0-9]{32})(/export)?',path)
            if m:
                report=load(ROOT/'runs/agent'/m[1]/'report.json')
                if report is None:return self.send(404,dict(error='Unknown investigation'))
                return self.send(200,report,filename='genetrace-'+m[1]+'.json' if m[2] else None)
            files={'/':'index.html','/app.js':'app.js','/style.css':'style.css'}
            if path in files:
                p=WEB/files[path];return self.send(200,p.read_bytes(),(mimetypes.guess_type(p.name)[0] or 'text/plain')+'; charset=utf-8')
            return self.send(404,dict(error='Not found'))
        except Exception:return self.send(500,dict(error='Could not read analysis state; retry shortly'))
    def do_POST(self):
        if not self.valid_host():return self.send(403,dict(error='Local host required'))
        allowed={f'http://127.0.0.1:{self.server.server_port}',f'http://localhost:{self.server.server_port}'}
        if self.headers.get('Origin') and self.headers['Origin'] not in allowed:return self.send(403,dict(error='Origin not allowed'))
        if self.path!='/api/investigate':return self.send(404,dict(error='Not found'))
        if self.headers.get('Content-Type','').split(';')[0]!='application/json':return self.send(415,dict(error='JSON required'))
        try:
            length=int(self.headers.get('Content-Length','0'))
            if not 0<length<=4096:return self.send(413,dict(error='Request too large'))
            body=json.loads(self.rfile.read(length));mode=body.get('mode','baseline')
            if mode=='meta':
                cfg=model_config()
                if not cfg['key'] or not cfg['model']:return self.send(409,dict(error='로컬 .env에 MODEL_API_KEY와 META_MODEL을 설정해야 실제 Meta 호출을 할 수 있습니다.'))
            if not BUSY.acquire(blocking=False):return self.send(409,dict(error='조사 한 건이 실행 중입니다. 완료 후 다시 시작해 주세요.'))
            try:job=Investigation(body.get('assembly',''),body.get('gene',''),mode)
            except Exception:BUSY.release();raise
            job.created=__import__('genetrace.agent',fromlist=['now']).now();job.persist('queued')
            def work():
                try:job.run()
                finally:BUSY.release()
            threading.Thread(target=work,daemon=True).start()
            return self.send(202,dict(id=job.id))
        except (ValueError,TypeError,AttributeError):return self.send(400,dict(error='유효한 assembly와 정확한 기준 유전자 이름을 선택하세요.'))
        except Exception:return self.send(500,dict(error='Investigation could not start'))

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--port',type=int,default=8765);args=ap.parse_args()
    server=ThreadingHTTPServer(('127.0.0.1',args.port),Handler)
    print(f'GeneTrace local workspace: http://127.0.0.1:{args.port}',flush=True)
    try:server.serve_forever()
    except KeyboardInterrupt:server.server_close()
if __name__=='__main__':main()
