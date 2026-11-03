"""Build an English narrated review video from verified project evidence and UI captures."""
from pathlib import Path
import subprocess, json, re, wave
from PIL import Image, ImageDraw, ImageFont, ImageOps
import imageio_ffmpeg
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'outputs/video_20261103_final'; OUT.mkdir(parents=True,exist_ok=True)
FF=imageio_ffmpeg.get_ffmpeg_exe()
SCENES=[('Before the biological claim', 'GeneTrace', ['Investigate genome reliability', 'Automate evidence gathering and analysis', 'A provisional judgment with traceable sources'], None, 'A gene appears to be missing. Is that biology, or a data problem? GeneTrace investigates genome reliability before researchers build experiments on an uncertain result.'), ('Why this matters', 'A missing label is not a missing gene.', ['Reads: DNA sequence fragments', 'Assembly: a reconstructed genome', 'Annotation: gene location and structure labels'], 'demo_why_verified_20261103.png', 'Sequencing produces DNA fragments called reads. Assembly reconstructs a genome, and annotation marks gene locations. Think of reconstructing a book and marking its chapters. A missing chapter label does not prove the text is gone.'), ('The research goal', 'Automate the investigation.', ['Collect evidence and run checks', 'Choose a follow-up based on the result', 'Reassess before drawing a biological conclusion'], None, 'Experts need data collection, computational analysis, and biological interpretation together. Our goal is to automate that investigation. A rare difference can be natural, so the agent must also know when the cause is unresolved.'), ('The demonstration dataset', 'Public genomes. Selected read checks.', ['100 quality-selected complete E. coli genomes', 'NCBI RefSeq: sequences and structural annotation', 'MG1655 comparison anchor · selected ENA reads', 'Not reads from all 100 genomes'], None, 'We use one hundred quality selected, complete Escherichia coli genomes from NCBI RefSeq, with gene annotations. MG1655 anchors the comparison panel. For this demonstration, seven hundred fifty six read pairs were selected from a public Illumina run. This reference selected subset is not whole run quality control or independent validation.'), ('Actual execution captures · time condensed', 'Meta selects a check. Code executes it.', ['JEONG-1266 · nimT · GCF_001558995.2', 'SRR2762258: selected original Illumina records', 'Fresh minimap2 alignment + regional analysis'], 'live_capture_20261103/input.png', 'We select nimT, a gene with frameshift related annotation, and enable fresh analysis. Meta requests structural annotation inspection, then read alignment. The application actually executes minimap2 and regional checks on fixed public inputs, and returns newly computed evidence.'), ('Actual execution captures · time condensed', 'A new result opens a follow-up.', ['Candidate shows coherent read disagreement', 'Neighboring-control analysis becomes available', 'Meta requests local controls and reassesses'], 'live_capture_20261103/frame_012.png', 'The candidate disagrees with selected reads. That observation exposes a neighboring control tool. Meta requests it to investigate whether the discrepancy is local. The application calculates nimR and yeaO from the new alignment, then updates the judgment. These are executed checks, not just saved result lookups.'), ('Fresh result · guarded judgment', 'A discrepancy needs context.', ['nimT: 3.09% mismatch', 'nimR: 2.07% · yeaO: 1.96%', 'Possible causes, limitations and next verification'], 'live_capture_20261103/judgment.png', 'The neighboring genes also disagree. A strictly local assembly error is therefore insufficient as an explanation. The structured result lists possible biological, assembly, read, or dataset linkage causes, with limitations and next verification. Exact run to assembly version linkage remains unverified; the cause is unresolved.'), ('Paired case study · same inputs and policy', 'What does Meta add?', ['Fixed policy: four fresh checks, no model', 'Meta: four fresh checks, five model responses', 'Equal measurements; selection reasons + cited synthesis'], None, 'We ran the same case with a fixed policy and with Meta. Both executed four fresh checks and produced equal measurements. Meta added evidence dependent selection reasons and cited synthesis. Both share minimum check safeguards, so this single case does not establish improved accuracy or unique AI discovery.'), ('Scope and reproducibility', 'Show what the evidence can support.', ['Population context: saved analysis', 'Candidate and controls: fresh execution', 'Commands, versions, input hashes and output hashes'], None, 'Population context is still a saved analysis. Fresh execution is enabled for one bounded case. Commands, versions, hashes, and evidence exports make the investigation reviewable. This is a working prototype, not a validated diagnostic system.'), ('GeneTrace · Meta hackathon', 'Evidence before interpretation.', ['Real tool execution + Meta-guided investigation', 'Explicit uncertainty and next verification', 'github.com/hope9901/GeneTrace_meta'], None, 'GeneTrace connects biological context, real analysis execution, and Meta guided interpretation, helping researchers investigate uncertain genome regions before deciding what to trust and what to verify next.')]
FONT='/System/Library/Fonts/Supplemental/Arial.ttf'; BOLD='/System/Library/Fonts/Supplemental/Arial Bold.ttf'
def f(n,b=False): return ImageFont.truetype(BOLD if b else FONT,n)
def wrap(draw,text,font,width):
    lines=[]; line=''
    for word in text.split():
        trial=(line+' '+word).strip()
        if draw.textlength(trial,font=font)>width and line: lines.append(line); line=word
        else: line=trial
    return lines+[line] if line else lines

def frame(scene,index,caption):
    eyebrow,title,points,screen,_=scene
    im=Image.new('RGB',(1920,1080),'#f3f5ed'); d=ImageDraw.Draw(im)
    d.rectangle((0,0,18,1080),fill='#24553d')
    d.text((75,46),'GeneTrace',font=f(35,True),fill='#183b2b')
    d.text((1530,54),f'{index+1:02d} / {len(SCENES):02d}',font=f(25),fill='#587063')
    d.text((75,139),eyebrow.upper(),font=f(24,True),fill='#52705e')
    width=1110 if screen else 1720
    y=198
    for line in wrap(d,title,f(68,True),width):
        d.text((72,y),line,font=f(68,True),fill='#173c2b'); y+=82
    y+=40
    for p in points:
        d.ellipse((80,y+16,93,y+29),fill='#34794d')
        lines=wrap(d,p,f(35),width-65)
        for line in lines:
            d.text((119,y),line,font=f(35),fill='#314d3d'); y+=49
        y+=25
    if screen:
        shot=Image.open(ROOT/'outputs'/screen).convert('RGB')
        shot.thumbnail((555,700),Image.Resampling.LANCZOS)
        x=1305+(555-shot.width)//2; sy=150
        d.rounded_rectangle((x-12,sy-12,x+shot.width+12,sy+shot.height+12),18,fill='#d9e2d5')
        im.paste(shot,(x,sy))
        d.text((1305,875),'ACTUAL EXECUTION CAPTURES' if index in [4,5,6] else 'ACTUAL UI · SAVED CONTEXT',font=f(18,True),fill='#52705e')
    d.rectangle((40,919,1880,1030),fill='#183b2b')
    caplines=wrap(d,caption,f(31),1740)
    cy=933+(2-len(caplines))*17
    for line in caplines:
        d.text(((1920-d.textlength(line,font=f(31)))/2,cy),line,font=f(31),fill='white');cy+=40
    d.rectangle((40,1050,40+1840*(index+1)/len(SCENES),1057),fill='#4c8b58')
    d.text((70,870),'PUBLIC-DATA PROTOTYPE  /  NO CAUSAL DIAGNOSIS OR ACCURACY CLAIM',font=f(18),fill='#617266')
    return im

def stamp(t):
    ms=round(t*1000); return f'{ms//3600000:02d}:{ms//60000%60:02d}:{ms//1000%60:02d},{ms%1000:03d}'
images=[]; audios=[]; subs=[]; cursor=0; n=0
for i,s in enumerate(SCENES):
    sentences=re.split(r'(?<!E\.)(?<=[.!?])\s+',s[4])
    for sentence in sentences:
        stem=f'{n:03d}'; txt=OUT/(stem+'.txt');audio=OUT/(stem+'.aiff');png=OUT/(stem+'.png')
        txt.write_text(sentence)
        subprocess.run(['say','-v','Samantha','-r','153','-f',str(txt),'-o',str(audio)],check=True)
        wav=OUT/(stem+'.wav')
        subprocess.run([FF,'-y','-v','error','-i',str(audio),'-c:a','pcm_s16le',str(wav)],check=True)
        with wave.open(str(wav),'rb') as a: dur=a.getnframes()/a.getframerate()
        sequence=([f'live_capture_20261103/frame_{j:03d}.png' for j in (range(0,12) if i==4 else range(12,24))] if i in [4,5] else [s[3]])
        sentence_index=sentences.index(sentence)
        picks=sequence[round(sentence_index*len(sequence)/len(sentences)):max(round((sentence_index+1)*len(sequence)/len(sentences)),round(sentence_index*len(sequence)/len(sentences))+1)]
        for k,shot in enumerate(picks):
            selected=(*s[:3],shot,s[4]);png=OUT/(stem+f'_{k:02d}.png')
            frame(selected,i,sentence).save(png)
            images.extend([f"file '{png}'",f'duration {dur/len(picks):.6f}'])
        audios.append(f"file '{audio}'")
        subs.append(f'{n+1}\n{stamp(cursor)} --> {stamp(cursor+dur)}\n{sentence}\n')
        cursor+=dur;n+=1
    print(f'scene {i+1} rendered; elapsed {cursor:.1f}s',flush=True)
images.append(f"file '{png}'")
(OUT/'frames.ffconcat').write_text('\n'.join(images)+'\n')
(OUT/'audio.ffconcat').write_text('\n'.join(audios)+'\n')
(OUT/'GeneTrace_EN.srt').write_text('\n'.join(subs))
(OUT/'script.json').write_text(json.dumps(SCENES,indent=2))
(OUT/'NARRATION_EN.md').write_text('# GeneTrace video narration\n\n'+'\n\n'.join('## '+s[1]+'\n\n'+s[4] for s in SCENES))
assert 90<cursor<180, f'Unexpected video duration: {cursor}'
result=OUT/'GeneTrace_EN.mp4'
subprocess.run([FF,'-y','-v','warning','-f','concat','-safe','0','-i',str(OUT/'frames.ffconcat'),'-f','concat','-safe','0','-i',str(OUT/'audio.ffconcat'),'-map','0:v','-map','1:a','-t',str(cursor),'-vf','fps=24','-c:v','libx264','-preset','veryfast','-crf','22','-threads','2','-pix_fmt','yuv420p','-c:a','aac','-b:a','160k','-movflags','+faststart',str(result)],check=True)
(OUT/'manifest.json').write_text(json.dumps({'duration_seconds':cursor,'resolution':'1920x1080','narration':'macOS Samantha synthesized English voice','presentation':'Edited actual execution capture sequence and explanatory cards; elapsed time condensed, not a continuous recording','execution_report':'f64d1b2f2a5c4761b1d023d3bd125796','new_alignment_performed':True,'fresh_tool_execution_loop':'bounded nimT case verified; other cases retain saved review','source_evidence':['demo/evidence/nimt.json','demo/evidence/provenance.json','demo/evidence/contamination.json']},indent=2))
print(f'COMPLETE {cursor:.2f}s {result}',flush=True)
