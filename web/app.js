'use strict';
const $=id=>document.getElementById(id);let assemblies=[],job=null,timer=null,currentReport=null,referenceCounts={};
const escapeHTML=s=>String(s??'').replace(/[&<>"']/g,c=>({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]));
const toolNames={targeted_homology:'단백질·유전체 상동성 검색',inspect_annotation:'구조 주석과 주변 유전자 검사',inspect_sequence:'후보 서열과 contig 경계 검사',lookup_read_evidence:'Raw read 근거와 적용 범위 확인',lookup_orthology_status:'전체 orthogroup 계산 상태 확인',lookup_domain_evidence:'Pfam 도메인 근거와 패널 범위 확인',lookup_regional_qc:'구간 mismatch·indel·clipping 검사',lookup_cohort_comparison:'동일 종 1:N 비교 근거 확인'};
async function api(path,options){const res=await fetch(path,options);const data=await res.json();if(!res.ok)throw new Error(data.error||'요청에 실패했습니다.');return data;}
function setAssemblyList(){const species=$('species').value;$('assembly').innerHTML=assemblies.filter(r=>r.species===species).map(r=>`<option value="${r.assembly_accession}">${escapeHTML(r.resolved_strain)} · ${r.assembly_accession}</option>`).join('');}
function assemblyDetail(){const r=assemblies.find(r=>r.assembly_accession===$('assembly').value);$('assembly-detail').textContent=r?`${r.assembly_level} · 단백질 ${Number(r.number_of_proteins).toLocaleString()}개 · ${r.sequencing_group}`:'';const bacterial=r?.species==='Escherichia coli';$('species-title').textContent=' / '+(r?.species||'');$('reference-label').textContent=bacterial?'MG1655':'S288C';$('gene-count').textContent=`기준 코딩 유전자 ${(referenceCounts[r?.reference_accession]||0).toLocaleString()}개 · 유전자명/좌위 ID 검색`;$('gene-options').innerHTML='';domainView().catch(e=>{$('domain-state').textContent='결과 조회 재시도 필요';});cohortView().catch(()=>{});}
function modeNote(){$('mode-note').textContent=$('mode').value==='meta'?'Meta 모델이 허용된 검사 중 다음 도구를 선택합니다. API 키는 로컬 .env에만 설정하세요.':'실제 BLAST 검색과 주석·서열 검사를 실행합니다. 모델 호출은 없습니다.';}
async function health(){try{const r=await api('/api/status');$('cohort-count').textContent=r.assemblies;referenceCounts=r.reference_counts||{};const selected=assemblies.find(a=>a.assembly_accession===$('assembly').value);$('gene-count').textContent=`기준 코딩 유전자 ${(referenceCounts[selected?.reference_accession]||0).toLocaleString()}개 · 유전자명/좌위 ID 검색`;$('orthology-state').textContent=r.orthology.status==='running'?'계산 중':r.orthology.status||'상태 확인 필요';$('read-state').textContent=r.read_pipeline.stage==='complete'?'분석 완료':r.read_pipeline.stage==='full_read_mapping'?'전체 read 분석 중':r.read_pipeline.stage||'확인 필요';$('model-state').textContent=r.meta_configured?`Meta 준비됨 · ${r.meta_model}`:'Meta 미설정 · 근거 검사 사용 가능';}catch(e){$('model-state').textContent='로컬 서버 연결 확인 필요';}}
function metric(value,label){return `<div class="metric"><strong>${escapeHTML(value)}</strong><span>${escapeHTML(label)}</span></div>`;}
function evidenceNote(e){const r=e.result;if(e.status==='error')return '검사 실패 · 생물학적 결론에 사용하지 않습니다.';if(e.status==='not_applicable')return '후보 위치 미확정 · 검사 적용 불가';switch(e.tool){case'targeted_homology':return r.best_genome_hit?`DNA identity ${r.best_genome_hit.pident}% · query coverage ${(r.best_genome_hit.query_coverage*100).toFixed(1)}% · 단백질 hit ${r.protein_hit_count}개`:'설정 기준의 DNA hit 없음 · 소실 판정 불가';case'inspect_annotation':return `겹치는 CDS ${r.overlapping_cds?.length||0}개 · ${r.recorded_disruption?'파괴된 모델 기록 있음':'frameshift/pseudogene 기록 미검출'}`;case'inspect_sequence':return `N ${(r.n_fraction*100).toFixed(2)}% · 가까운 contig 경계까지 ${r.distance_to_contig_edge.toLocaleString()}bp`;case'lookup_read_evidence':return typeof r.regional_read_analysis==='object'?`${r.regional_read_analysis.partial_input?'부분 read 검사':'전체 read 검사'} · assembly 사용 run 연결은 미확정`:'이 대상 구간의 read 분석 없음 · 등록 정보만 조회';case'lookup_cohort_comparison':return r.available_for_target?`${r.comparison.total_assemblies}개 집단 · 위치/주변 유전자 지지 ${r.comparison.unambiguous_synteny_supported}개 · 오류 판정 아님`:`집단 비교 ${r.state.stage} · 비교 전에는 population 해석을 보류합니다.`;case'lookup_regional_qc':return r.available_for_target?r.regional_audits.map(a=>`mismatch ${(a.metrics.mismatch_fraction*100).toFixed(2)}% · spanning ${a.metrics.spanning_plus100_records} records · 원인 미해결`).join('; '):'이 대상 구간의 실제 read 검사는 연결되지 않았습니다.';case'lookup_domain_evidence':return '선택한 Pfam 6개 패널의 도메인 근거 · orthology/기능 판정 아님';case'lookup_orthology_status':return `${r.status.status||'확인 필요'} · 완료 전 orthogroup 소속 추정 안 함`;default:return '';}}
function drawTrack(report){const c=report.region;if(!c){$('region-label').textContent='Search pending';$('genome-track').className='track-empty';$('genome-track').textContent='상동성 검사 후 후보 위치를 표시합니다.';return;}
 $('region-label').textContent=`${c.sseqid}:${c.start.toLocaleString()}–${c.end.toLocaleString()}`;
 const ann=report.evidence.find(e=>e.tool==='inspect_annotation')?.result;const gs=ann?.neighboring_genes||[];const lo=Math.max(1,Math.min(c.start-1500,...gs.map(g=>g.start)));const hi=Math.max(c.end+1500,...gs.map(g=>g.end));const x=n=>35+(n-lo)/(hi-lo)*730;
 const rect=(start,end,y,color)=>`<rect x="${x(start)}" y="${y}" width="${Math.max(2,x(end)-x(start))}" height="10" rx="2" fill="${color}"/>`;
 $('genome-track').className='';$('genome-track').innerHTML=`<svg viewBox="0 0 800 125" role="img" aria-label="상동 후보와 주변 구조 주석의 위치"><line x1="35" x2="765" y1="67" y2="67" stroke="#dfe5de"/>${[0,1,2,3,4].map(i=>{const n=lo+(hi-lo)*i/4;return `<line x1="${x(n)}" x2="${x(n)}" y1="20" y2="90" stroke="#edf0e9"/><text x="${x(n)}" y="111" text-anchor="middle">${Math.round(n).toLocaleString()}</text>`;}).join('')}${gs.map(g=>rect(g.start,g.end,62,'#6c8a78')).join('')}${rect(c.start,c.end,30,'#b77c43')}<text x="${x(c.start)}" y="20">${escapeHTML(report.gene)} 상동 후보</text></svg>`;
}
function render(report){currentReport=report;const done=report.status==='complete',failed=report.status==='failed';$('execution-badge').textContent=report.execution_mode==='meta'?`Meta · ${report.model_calls}회 호출`:'근거 검사 · 모델 미사용';
 $('verdict').innerHTML=`<div class="eyebrow green">${escapeHTML(report.strain)} / ${escapeHTML(report.gene)} · ${failed?'FAILED':done?'INVESTIGATION COMPLETE':'INVESTIGATING'}</div><h3>${escapeHTML(failed?'조사를 완료하지 못했습니다':report.evidence.length?report.assessment.label:'근거를 모으고 있습니다')}</h3><p class="muted">${escapeHTML(report.error||report.assessment.explanation)}</p>`;
 const h=report.evidence.find(e=>e.tool==='targeted_homology')?.result;const read=report.evidence.find(e=>e.tool==='lookup_read_evidence')?.result.regional_read_analysis;
 $('metrics').innerHTML=h?metric(h.best_genome_hit?`${h.best_genome_hit.pident}%`:'—','DNA identity')+metric(h.best_genome_hit?`${(h.best_genome_hit.query_coverage*100).toFixed(1)}%`:'—','Query coverage')+metric(h.protein_hit_count??'—','단백질 hit')+metric(report.evidence.length,'수집된 검사 근거'):'';
 $('limits').textContent=read&&typeof read==='object'&&read.partial_input?'부분 read 근거입니다. coverage 부족·소실 판정에 사용할 수 없습니다.':report.limitations[0];
 if(read&&typeof read==='object'&&read.fastq?.all_quality_values_zero)$('limits').textContent+=' 관찰한 품질 문자가 모두 !입니다. BQ 필터로 전부 제거될 수 있으며, 이것만으로 저품질 sequencing을 뜻하지 않습니다.';
 $('evidence-list').className='';$('evidence-list').innerHTML=report.evidence.map(e=>`<div class="evidence-row"><span class="evidence-id">${e.id}</span><div><strong>${escapeHTML(toolNames[e.tool]||e.tool)}</strong><p>${escapeHTML(evidenceNote(e))}</p></div><button class="text-button inspect-evidence" data-id="${e.id}">원시 근거 ↗</button></div>`).join('')||'<p class="muted">실제 서열 검색을 실행하는 중입니다.</p>';
 document.querySelectorAll('.inspect-evidence').forEach(b=>b.addEventListener('click',()=>{const e=report.evidence.find(e=>e.id===b.dataset.id);$('dialog-title').textContent=`${e.id} · ${toolNames[e.tool]}`;$('dialog-content').textContent=JSON.stringify(e,null,2);$('evidence-dialog').showModal();}));
 let next=[];const code=report.assessment.code;
 const population=report.evidence.find(e=>e.tool==='lookup_cohort_comparison')?.result;
 if(!population?.available_for_target)next.push('먼저 동일 종의 1:N 비교에서 같은 위치·주변 유전자 배치·자연 변이·기술/프로젝트 분포를 확인합니다. 현재 집단 비교 맥락은 미확정입니다.');
 if(read&&typeof read==='object'&&read.fastq?.all_quality_values_zero)next.push('원본 품질값과 변환 이력을 확인합니다. ! 값만으로 sequencing 실패를 추론하거나 BQ 필터 결과를 진짜 결손으로 해석하지 않습니다.');
 if(code==='annotation_scope_mismatch')next.push('제공된 핵/소기관 주석의 범위를 맞춘 뒤 유전자군 비교를 다시 수행합니다.');
 if(code==='annotation_gap_candidate')next.push('후보 CDS의 시작·종결 코돈과 reading frame을 확인하고 동일한 기준으로 구조 주석을 비교합니다.');
 if(code==='disrupted_model_unresolved')next.push('기록된 frameshift 위치에서 read의 indel 지지와 haplotype을 검사합니다.');
 if(!h?.best_genome_hit)next.push('검색 기준·분절 hit·주변 synteny를 추가 확인합니다. 미검출을 소실로 해석하지 않습니다.');
 if(typeof read!=='object')next.push('해당 assembly에 실제 사용한 raw run을 확인하고 후보 구간의 read를 정렬합니다.');
 else if(read.partial_input)next.push('전체 파일 검증 후 진행 중인 full-read 분석을 검토합니다. 부분 coverage를 전체 coverage로 해석하지 않습니다.');
 else next.push('전체 read의 per-base consensus, 작은 indel과 구조적 연결을 추가 점검합니다.');
 next.push('전체 orthogroup 결과에서 paralog와 유전자군 소속을 확인하고, 반복·배수성·기능 해석은 별도로 검토합니다.');
 $('next-actions').innerHTML='<ol>'+next.map(n=>`<li>${escapeHTML(n)}</li>`).join('')+'</ol>';
 $('model-hypotheses').innerHTML=report.model_hypotheses.length?'<strong>Meta가 제안한 가설 · 연구자 검토 필요</strong>'+report.model_hypotheses.map(h=>`<p>${escapeHTML(h.text)} <span class="evidence-id">${h.evidence_ids.join(', ')}</span></p>`).join(''):'';
 $('export').hidden=false;drawTrack(report);renderRegional(report);
 if(done||failed){clearInterval(timer);timer=null;$('run-button').disabled=false;$('run-button').innerHTML='새 조사 시작 <span>↗</span>';if(failed)$('error').textContent=report.error;}
}
async function poll(){try{render(await api(`/api/jobs/${job}`));}catch(e){$('error').textContent=e.message;}}
$('investigate-form').addEventListener('submit',async e=>{e.preventDefault();$('error').textContent='';$('run-button').disabled=true;$('run-button').textContent='실제 근거 검사 중…';try{const r=await api('/api/investigate',{method:'POST',headers:{'Content-Type':'application/json'},body:JSON.stringify({assembly:$('assembly').value,gene:$('gene').value.trim(),mode:$('mode').value})});job=r.id;localStorage.setItem('genetrace-job',job);if(timer)clearInterval(timer);timer=setInterval(poll,1500);await poll();}catch(e){$('error').textContent=e.message;$('run-button').disabled=false;$('run-button').textContent='조사 시작 ↗';}});
$('assembly').addEventListener('change',assemblyDetail);$('mode').addEventListener('change',modeNote);let geneTimer;
$('gene').addEventListener('input',()=>{clearTimeout(geneTimer);geneTimer=setTimeout(async()=>{try{const rows=await api('/api/genes?q='+encodeURIComponent($('gene').value)+'&assembly='+encodeURIComponent($('assembly').value));$('gene-options').innerHTML=rows.map(g=>`<option value="${escapeHTML(g.gene_name||g.locus_tag)}">${escapeHTML(g.locus_tag+' · '+g.ncbi_protein_id)}</option>`).join('');}catch(e){$('error').textContent=e.message;}},200);});
$('close-dialog').addEventListener('click',()=>$('evidence-dialog').close());$('export').addEventListener('click',()=>{if(job)window.location.href=`/api/jobs/${job}/export`;});
document.querySelectorAll('.example').forEach(b=>b.addEventListener('click',()=>{$('species').value='Saccharomyces cerevisiae';setAssemblyList();$('assembly').value=b.dataset.assembly;$('gene').value=b.dataset.gene;assemblyDetail();$('run-button').focus();}));
(async()=>{try{assemblies=await api('/api/assemblies');setAssemblyList();$('assembly').value='GCA_903819125.2';assemblyDetail();await health();setInterval(()=>{health();domainView().catch(()=>{});cohortView().catch(()=>{});},15000);const old=localStorage.getItem('genetrace-job');if(old&&/^[a-f0-9]{32}$/.test(old)){try{const r=await api('/api/jobs/'+old);job=old;const prior=assemblies.find(a=>a.assembly_accession===r.assembly);if(prior){$('species').value=prior.species;setAssemblyList();$('assembly').value=r.assembly;$('gene').value=r.gene;$('mode').value=r.execution_mode;assemblyDetail();modeNote();}render(r);if(!['complete','failed'].includes(r.status)){timer=setInterval(poll,1500);$('run-button').disabled=true;}}catch(e){localStorage.removeItem('genetrace-job');}}}catch(e){$('error').textContent=e.message;}finally{if(!timer)$('run-button').disabled=false;}})();

async function domainView(){
 const acc=$('assembly').value;if(!acc)return;document.querySelectorAll('.scope-warning').forEach(e=>e.remove());const profile=$('domain-profile').value;
 const [all,detail]=await Promise.all([api('/api/domain-results?assembly='+acc),api('/api/domain-results?assembly='+acc+'&profile='+profile)]);
 if($('assembly').value!==acc||$('domain-profile').value!==profile)return;
 $('domain-state').textContent=`${all.status.stage==='complete'?'완료':'검사 중'} · ${all.status.completed||0}/${all.status.total||'—'} 유전체`;
 $('domain-counts').innerHTML=all.assembly_scanned?all.counts.map(r=>`<div class="domain-cell"><span>${escapeHTML(r.profile)}</span><strong>${r.proteins_with_domain} / ${r.domain_instances}</strong><span>단백질 수 / 도메인 반복 수</span></div>`).join(''):'<p class="muted small">이 assembly는 아직 검색되지 않았습니다. 0으로 처리하지 않습니다.</p>';
 if(all.counts.some(r=>r.annotation_scope_warning==='True'))$('domain-counts').insertAdjacentHTML('beforebegin','<p class="notice scope-warning">Organelle-only annotation: nuclear domain presence/absence is unresolved.</p>');
 $('domain-members').innerHTML=detail.members.length?`<p>${detail.matching_domain_instances} instances / first ${Math.min(20,detail.members.length)} shown</p>`+detail.members.slice(0,20).map(h=>`<div class="domain-member"><b>${escapeHTML(h.protein_id)}</b> · ${escapeHTML(h.profile_version)}<br>AA ${h.aa_start}–${h.aa_end} · HMM coverage ${(Number(h.profile_coverage)*100).toFixed(1)}% · ${escapeHTML(h.coverage_class)}</div>`).join(''):'선택 패널에서 보고된 hit 없음. 유전자/도메인 소실 판정 불가.';
}
$('domain-profile').addEventListener('change',()=>domainView().catch(e=>{$('error').textContent=e.message;}));
$('species').addEventListener('change',()=>{setAssemblyList();$('gene').value=$('species').value==='Escherichia coli'?'rpoB':'YJU3';assemblyDetail();});

function renderRegional(report){
 const e=report.evidence.find(e=>e.tool==='lookup_regional_qc');const rows=e?.result?.regional_audits||[];
 $('regional-qc-result').innerHTML=rows.length?rows.map(a=>`<div class="notice"><strong>${escapeHTML(e.id)} / ${escapeHTML(report.strain)} · ${escapeHTML(report.gene)}</strong><p>mismatch ${(a.metrics.mismatch_fraction*100).toFixed(3)}% · 평균 depth ${a.metrics.mean_depth.toFixed(1)} · 양쪽 100bp 지지 ${a.metrics.spanning_plus100_records} records</p><p>일관된 substitution ${a.coherent_substitution_count} / indel ${a.coherent_indel_count} · boundary clipping ${(a.metrics.boundary_clipped_fraction*100).toFixed(2)}%</p><p>${escapeHTML(a.signals.join(', '))} · 원인 미해결</p>${a.quality_information_unresolved?'<p>품질 정보 미확정: base quality로 신뢰도를 검증하지 못했습니다.</p>':''}</div>`).join(''):'이 대상 구간의 실제 read 검사가 아직 연결되지 않았습니다. 정상으로 처리하지 않습니다.';
}
async function validationHealth(){try{const data=await api('/api/validation');const b=data.benchmark;const test=(b.groups||[]).filter(g=>g.split==='holdout');if(!test.length){$('validation-summary').textContent='검증 작업 진행 중 또는 결과 없음';return;}const recovered=test.reduce((n,g)=>n+g.expected_signal_detected,0);const total=test.reduce((n,g)=>n+g.cases,0);$('validation-summary').textContent=`시뮬레이션 ${b.loci}개 구간 / ${b.cases}건 · holdout 기대 패턴 ${recovered}/${total}건 재현. 실제 오류 진단 성능은 아직 미확립입니다.`;}catch(e){$('validation-summary').textContent='검증 결과 조회 필요';}}
validationHealth();

async function cohortView(){
 const acc=$('assembly').value,gene=$('gene').value.trim();if(!acc||!gene)return;
 const data=await api('/api/cohort-comparison?assembly='+encodeURIComponent(acc)+'&gene='+encodeURIComponent(gene));
 if(acc!==$('assembly').value||gene!==$('gene').value.trim())return;
 $('cohort-comparison-state').textContent=data.state.stage==='complete'?'비교 완료':`${data.state.stage} · ${data.state.completed||0}/${data.state.total||'—'}`;
 if(!data.available_for_target){$('cohort-comparison-result').textContent='현재 대상의 1:N 비교가 계산 중이거나 아직 없습니다. 단일 비교만으로 데이터 오류를 판정하지 않습니다.';return;}
 const s=data.summary;$('cohort-comparison-result').innerHTML=`<p>전체 ${s.total_assemblies}개 / 동일 위치·양쪽 주변 유전자 지지 ${s.unambiguous_synteny_supported}개. 비교 불가·모호한 위치는 빈도 분모에서 제외합니다.</p>`+s.markers.map(m=>`<div class="notice"><strong>대상 좌표 ${m.focal_reference_position} · ${escapeHTML(m.repeat_base)} × ${m.focal_repeat_length}</strong><p>다른 균주 ${m.population.callable_assemblies}개에서 비교 가능 / ${m.population.uncallable_assemblies}개 비교 불가</p><p>반복 allele: ${Object.entries(m.population.allele_counts).map(([a,n])=>`${escapeHTML(a||'(deletion)')}: ${n}`).join(' · ')}</p><p>${s.focal_technology==='N/A'?'대상 sequencing 기술 N/A · 기술별 비교 불가':`같은 기술 조합 ${m.same_technology.callable_assemblies}개 / 다른 알려진 기술 ${m.different_known_technology.callable_assemblies}개.`} 빈도는 기술적 오류의 증명이 아닙니다.</p></div>`).join('');
}
$('gene').addEventListener('change',()=>cohortView().catch(()=>{}));

let backgroundLoading=false;
async function backgroundView(){
 if(backgroundLoading)return;backgroundLoading=true;$('refresh-background').disabled=true;
 try{
  const data=await api('/api/background-jobs');
  const time=new Date(data.checked_utc).toLocaleTimeString('ko-KR',{timeZone:'Asia/Seoul',hour12:false});
  $('background-summary').textContent=`최근 확인 ${time} KST · 디스크 여유 ${data.free_disk_gib} GiB · 작업 보호 기준 ${data.reserve_disk_gib} GiB`;
  const names={running:'실행 중',not_running:'프로세스 중단 · 확인 필요',unknown:'실행 여부 확인 불가',paused_by_user:'사용자 중지',resource_limit:'자원 부족으로 중지',complete:'완료',failed:'실패 · 확인 필요',failed_resumable:'실패 · 재개 가능',not_available:'등록된 상태 없음',running_despite_status:'상태와 실행 정보 불일치'};
  const stages={core_inference:'공통 그룹·계통 추론',all_genome_assignment:'나머지 유전체 배정',validating_all_gene_accounting:'전체 유전자 검증'};
  $('background-jobs').innerHTML=data.jobs.map(j=>{
   const p=j.progress;const percent=p?Math.min(100,Math.max(0,Number(p.percent))):null;
   const stamp=j.updated_utc?new Date(j.updated_utc).toLocaleTimeString('ko-KR',{timeZone:'Asia/Seoul',hour12:false}):'—';
   const stale=j.status==='running'&&j.heartbeat_stale;
   return `<article class="background-job"><div class="background-job-head"><strong>${escapeHTML(j.title)}</strong><span class="tag">${escapeHTML(names[j.observed_status]||j.observed_status)}</span></div><p class="small">${escapeHTML(stages[j.stage]||j.stage||'단계 확인 중')} · 유전체 ${Number(j.assemblies||0).toLocaleString()}개 · 코딩 유전자 ${Number(j.coding_units||0).toLocaleString()}개</p>${p?`<div class="background-progress-label"><span>현재 core 검색 결과 ${Number(p.result_files).toLocaleString()} / ${Number(p.total).toLocaleString()}</span><b>${percent}%</b></div><progress max="100" value="${percent}" aria-label="${escapeHTML(j.title)} core 검색 결과 파일 비율"></progress><p class="field-note">${escapeHTML(p.note)} ${p.active_searches}개 검색 실행 중.</p>`:'<p class="field-note">현재 단계의 수치 진행률은 제공되지 않습니다.</p>'}<p class="small muted">${j.process_alive===true?`실행 확인 · PID ${Number(j.pid)} · 작업 프로세스 ${Number(j.active_processes)}개 · CPU 합계 ${Number(j.cpu_percent)}%`:j.process_alive===false?'현재 PID에서 해당 작업을 확인하지 못했습니다.':'프로세스 직접 확인 불가'}<br>검색 스레드 ${Number(j.configured_search_threads||0)}개 · 상태 갱신 ${escapeHTML(stamp)} KST${stale?' · 60초 이상 갱신되지 않음':''}</p>${j.log_lines.length?`<details><summary>최근 진행 로그</summary><pre>${escapeHTML(j.log_lines.join('\n'))}</pre></details>`:''}</article>`;
  }).join('');
  $('background-completed').innerHTML=data.other_workflows.map(w=>`<p>${escapeHTML(w.title)} · ${escapeHTML(names[w.stage]||w.stage)}</p>`).join('');
 }catch(e){$('background-summary').textContent='상태 조회 실패 · 서버 연결을 확인하세요. 아래 정보는 마지막 조회 결과입니다.';}
 finally{backgroundLoading=false;$('refresh-background').disabled=false;}
}
$('refresh-background').addEventListener('click',backgroundView);
backgroundView();setInterval(backgroundView,15000);
