"""Submission-cut video. Actual UI close-ups, accessible language, checked biological context.
Video outputs stay local; do not publish the MP4 to GitHub.
"""
from pathlib import Path
import json,re,subprocess,wave,shutil,hashlib
from PIL import Image,ImageDraw,ImageFont,ImageOps
import imageio_ffmpeg
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'outputs/_video_build';OUT.mkdir(parents=True,exist_ok=True)
FF=imageio_ffmpeg.get_ffmpeg_exe();CAP=ROOT/'assets/video_sources'
SCENES=[
 dict(title='Can this genome region be trusted?',tag='GENETRACE · BEFORE THE BIOLOGICAL CLAIM',layout='plain',points=['A gene looks broken.','Is it biology, or a problem in the data?','Investigate before planning the next experiment.'],speech='A gene looks broken. Is that a real biological difference, or a problem in the data? GeneTrace investigates genome reliability before researchers plan their next experiment.'),
 dict(title='A missing label does not prove gene loss.',tag='THREE SIMPLE STEPS',layout='plain',points=['Reads: pieces of DNA sequence',"Assembly: reconstruction of a genome (the organism's DNA)",'Annotation: labels for genes (regions carrying biological instructions)'],speech='Sequencing reads DNA pieces. Assembly reconstructs the genome. Annotation labels the genes. Like a book, missing a chapter label does not prove the text is gone.'),
 dict(title='Meta plans. GeneTrace executes and checks.',tag='META MODEL API · muse-spark-1.3',layout='plain',points=['Meta receives the question and available evidence.','Code defines required checks and allowed choices.','Meta selects a check; code executes it.','New results return to Meta for a cited interpretation.'],speech="Meta's Muse Spark one point three receives the question and evidence. Code defines required checks and allowed choices. Meta selects a check, code executes it, and new results return for interpretation."),
 dict(title='Public data, with a clearly defined scope.',tag='E. COLI · A WELL-STUDIED BACTERIUM',layout='plain',points=['100 complete NCBI RefSeq genomes: comparison context','JEONG-1266 · GCF_001558995.2: the demonstration case','756 selected Illumina read pairs · SRR2762258','Regional reanalysis, not whole-run quality control'],speech='We use one hundred complete E. coli genomes for comparison context. This demonstration reanalyzes seven hundred fifty six selected public read pairs from one strain. It is not whole run quality control.'),
 dict(title='Why these three genes?',tag='BIOLOGICAL CONTEXT · NOT A FUNCTIONAL ASSAY',layout='genes',speech="Nim T is associated with transport of an antibacterial compound. Nim R regulates gene activity. Yea O's function remains unclear. A false broken gene label could misdirect biological research. Here, the neighbors are nearby comparison sites, not assumed error free controls."),
 dict(title='The candidate check really executes.',tag='RECORDED EXECUTION · ACTUAL UI CLOSE-UP',layout='capture',shot='candidate.png',crop=[30,345,450,545],note=['Fresh alignment of selected reads','3.09% mismatch at nimT','Mismatch = read / genome base disagreement'],speech='Nim T is labeled as disrupted. GeneTrace runs a fresh read alignment and measures disagreement with the genome. About three percent of compared bases differ. That observation needs context.'),
 dict(title='Meta requests the neighboring check.',tag='ACTUAL RECORDED DECISION AND EXECUTION',layout='capture',shot='neighbor_decision.png',crop=[35,345,448,628],note=['Candidate disagreement observed','Meta asks for nearby comparison','Code executes the selected check'],speech='Meta requests nearby genes to test whether the discrepancy is confined to the candidate. The selected check actually runs, and the new evidence updates the assessment.'),
 dict(title='Neighbors disagree too. The cause stays open.',tag='NEWLY CALCULATED REGIONAL RESULTS',layout='capture',shot='neighbor_results.png',crop=[30,345,450,550],note=['nimT: 3.09% mismatch','nimR: 2.07% · yeaO: 1.96%','A local-only cause is not established'],speech='Nim R and Yea O also show about two percent disagreement. These data cannot isolate a Nim T specific assembly error. Biological variation, technical effects, or mismatched datasets remain possible.'),
 dict(title='This example demonstrates the working loop.',tag='PROTOTYPE BEHAVIOR · NOT AN ACCURACY BENCHMARK',layout='capture',shot='judgment.png',crop=[43,330,444,610],note=['Model choice → real tool execution','New evidence → provisional reassessment','Sources and next verification remain visible'],speech='This example demonstrates model choice, real analysis, and reassessment. A fixed policy produced the same measurements; Meta adds selection reasons and cited interpretation. Exact read to assembly version linkage remains unverified.'),
 dict(title='Next: follow up BUSCO completeness flags.',tag='FUTURE DIRECTION · NOT YET IMPLEMENTED',layout='plain',points=['BUSCO: conserved genes usually present once in related species','Investigate apparently missing or fragmented markers','Connect genome, gene labels and read evidence','A missing marker is a question, not proof of sequencing error'],speech='Next, BUSCO could prioritize apparently missing or fragmented conserved genes for this investigation. These markers are usually present once in related species; they are not all essential genes. This extension is planned, not implemented.'),
 dict(title='Evidence before interpretation.',tag='GENETRACE · META HACKATHON',layout='plain',points=['Bring biological context and computational checks together.','Show what was found, what remains uncertain, and what to verify.','github.com/hope9901/GeneTrace_meta'],speech='GeneTrace connects biological context, real analysis, and Meta guided interpretation, helping researchers decide what to trust and what to verify next.')
]
REG='/System/Library/Fonts/Supplemental/Arial.ttf';BOLD='/System/Library/Fonts/Supplemental/Arial Bold.ttf'
def font(n,b=False):return ImageFont.truetype(BOLD if b else REG,n)
def lines(d,text,ft,width):
 out=[];line=''
 for word in text.split():
  candidate=(line+' '+word).strip()
  if line and d.textlength(candidate,font=ft)>width:out.append(line);line=word
  else:line=candidate
 return out+[line] if line else out

def text(d,x,y,value,size=38,width=1700,bold=False,color='#294b3a',gap=1.3):
 ft=font(size,bold);rows=lines(d,value,ft,width)
 for row in rows:d.text((x,y),row,font=ft,fill=color);y+=round(size*gap)
 return y

def render(s,i,caption):
 im=Image.new('RGB',(1920,1080),'#f4f6ef');d=ImageDraw.Draw(im);d.rectangle((0,0,18,1080),fill='#204b35')
 d.text((65,35),'GeneTrace',font=font(30,True),fill='#204b35');d.text((1740,42),f'{i+1:02}/{len(SCENES):02}',font=font(23),fill='#5d7766')
 d.text((65,100),s['tag'],font=font(23,True),fill='#52715f')
 title_end=text(d,65,151,s['title'],size=52,width=1780,bold=True)
 if s['layout']=='capture':
  # Crop an actual screenshot, preserving its pixels and the relevant surrounding context.
  shot=Image.open(CAP/s['shot']).convert('RGB').crop(tuple(s['crop']));shot=ImageOps.contain(shot,(1190,610),Image.Resampling.LANCZOS)
  x=65+(1190-shot.width)//2;y=245+(610-shot.height)//2
  d.rounded_rectangle((55,232,1267,871),15,fill='#dae4d4');im.paste(shot,(x,y))
  yy=285
  for point in s['note']:
   yy=text(d,1320,yy,point,size=35,width=520,bold=True)+36
  text(d,1320,yy+15,'Recorded run · enlarged UI excerpt',size=23,width=520,color='#647a6a')
 elif s['layout']=='genes':
  data=[('nimT','TARGET','Compound transporter','Associated with 2-nitroimidazole transport.'),('nimR','NEIGHBOR','Gene-activity regulator','Reference E. coli regulator of this locus.'),('yeaO','NEIGHBOR','Function not established','Precise biological role remains unclear.')]
  for j,(gene,role,title,body) in enumerate(data):
   x=65+j*605;d.rounded_rectangle((x,270,x+565,675),18,fill='#e2ebdb')
   text(d,x+28,300,role,size=22,color='#5c765e');text(d,x+28,355,gene,size=62,bold=True);text(d,x+28,450,title,size=33,width=500,bold=True);text(d,x+28,550,body,size=27,width=500)
  text(d,65,725,'Why check the neighbors? Ask whether the same disagreement extends nearby.',size=35,width=1780,bold=True)
  text(d,65,815,'Reference biology + strain annotation; no phenotype or gene-essentiality test in this demo.',size=24,width=1780)
 else:
  y=max(275,title_end+55)
  for point in s['points']:
   d.ellipse((70,y+17,83,y+30),fill='#428052');y=text(d,110,y,point,size=39,width=1730)+33
  if i==9:text(d,65,830,'BUSCO user guide · busco.ezlab.org/busco_userguide',size=23,color='#647a6a')
 caplines=lines(d,caption,font(32),1760);assert len(caplines)<=2, 'Caption exceeds two lines'
 d.rectangle((40,922,1880,1032),fill='#1e4833');y=937+(2-len(caplines))*19
 for row in caplines:d.text(((1920-d.textlength(row,font=font(32)))/2,y),row,font=font(32),fill='white');y+=41
 d.rectangle((40,1051,40+1840*(i+1)/len(SCENES),1059),fill='#4b885b')
 return im

def stamp(t):
 m=round(t*1000);return f'{m//3600000:02}:{m//60000%60:02}:{m//1000%60:02},{m%1000:03}'
cache={} # Voice files are regenerated; no dependency on deleted drafts.
silence={}
for duration in [2,3]:
 p=OUT/f'silence_{duration}.aiff'
 subprocess.run([FF,'-y','-v','error','-f','lavfi','-i','anullsrc=r=22050:cl=mono','-t',str(duration),'-c:a','pcm_s16be',str(p)],check=True);silence[duration]=p
pictures=[];audios=[];srt=[];timeline=[];cursor=0;n=0
for i,s in enumerate(SCENES):
 start=cursor
 for sentence in re.split(r'(?<!E\.)(?<=[.!?])\s+',s['speech']):
  stem=f'{n:03}';txt=OUT/(stem+'.txt');aif=OUT/(stem+'.aiff');wav=OUT/(stem+'.wav');png=OUT/(stem+'.png');txt.write_text(sentence)
  if sentence in cache:shutil.copyfile(cache[sentence].with_suffix('.aiff'),aif);shutil.copyfile(cache[sentence].with_suffix('.wav'),wav)
  else:
   subprocess.run(['say','-v','Samantha','-r','153','-f',str(txt),'-o',str(aif)],check=True)
   subprocess.run([FF,'-y','-v','error','-i',str(aif),'-c:a','pcm_s16le',str(wav)],check=True)
  with wave.open(str(wav),'rb') as w:duration=w.getnframes()/w.getframerate()
  assert duration>.1
  render(s,i,sentence).save(png);pictures +=[f"file '{png}'",f'duration {duration:.6f}'];audios.append(f"file '{aif}'");srt.append(f'{n+1}\n{stamp(cursor)} --> {stamp(cursor+duration)}\n{sentence}\n');cursor+=duration;n+=1
 hold=3 if i==len(SCENES)-1 else 2
 pictures += [f"file '{png}'",f'duration {hold:.6f}'];audios.append(f"file '{silence[hold]}'");cursor+=hold
 srt[-1]=re.sub(r' --> \d{2}:\d{2}:\d{2},\d{3}', ' --> '+stamp(cursor), srt[-1],count=1)
 timeline.append(dict(scene=i+1,title=s['title'],start_seconds=start,end_seconds=cursor));print('Scene',i+1,'ends',round(cursor,2),flush=True)
assert 100<cursor<180,cursor
assert abs(sum(float(x.split()[1]) for x in pictures if x.startswith('duration'))-cursor)<.001
pictures.append(f"file '{png}'");(OUT/'frames.ffconcat').write_text('\n'.join(pictures)+'\n');(OUT/'audio.ffconcat').write_text('\n'.join(audios)+'\n');(OUT/'GeneTrace_EN.srt').write_text('\n'.join(srt));(OUT/'timeline.json').write_text(json.dumps(timeline,indent=2));(OUT/'script.json').write_text(json.dumps(SCENES,indent=2));(OUT/'NARRATION_EN.md').write_text('# Final submission narration\n\n'+'\n\n'.join('## '+s['title']+'\n\n'+s['speech'] for s in SCENES))
subprocess.run([FF,'-y','-v','warning','-f','concat','-safe','0','-i',str(OUT/'frames.ffconcat'),'-f','concat','-safe','0','-i',str(OUT/'audio.ffconcat'),'-t',str(cursor),'-vf','fps=24','-c:v','libx264','-crf','21','-preset','veryfast','-threads','2','-pix_fmt','yuv420p','-c:a','aac','-ar','48000','-b:a','160k','-movflags','+faststart',str(OUT/'GeneTrace_EN.mp4')],check=True)
(OUT/'manifest.json').write_text(json.dumps(dict(duration_seconds=cursor,resolution='1920x1080',execution_report='f64d1b2f2a5c4761b1d023d3bd125796',narration='macOS Samantha synthesized voice',presentation='Actual recorded investigation UI close-ups plus explanatory cards; cropped for legibility; not a new live API run or continuous screen recording',gene_function_scope='Reference E. coli biology and distributed JEONG-1266 annotation; no functional phenotype or essentiality assay',BUSCO='Future extension; not implemented; conserved markers are not all essential genes',video_distribution='Local only. Do not upload MP4 to GitHub.'),indent=2))
print('COMPLETE',cursor,flush=True)
