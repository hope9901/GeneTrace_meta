"""Normalize only explicit NCBI provenance and compare sequencing groups."""
import csv
import json
import os
import re
from collections import Counter
from pathlib import Path

os.environ.setdefault('MPLCONFIGDIR', '/tmp/genetrace-matplotlib')
os.environ.setdefault('XDG_CACHE_HOME', '/tmp/genetrace-cache')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import StrMethodFormatter
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / 'data/inventory/s_cerevisiae_methods_2026-10-30'
BASE = ROOT / 'outputs/s_cerevisiae_20261030/analysis.json'
OUT = ROOT / 'outputs/s_cerevisiae_methods_20261030'
NA = 'N/A'

ASSEMBLERS = {
    'Canu': r'\bcanu\b', 'SMARTdenovo': r'\bsmartdenovo\b', 'MIRA': r'\bmira\b',
    'Newbler': r'\bnewbler\b', 'SPAdes': r'\b(?:meta)?spades\b', 'ABySS': r'\babyss\b',
    'HGAP': r'\bhgap\w*', 'CLC Genomics Workbench': r'\bclc\b', 'Phrap': r'\bphrap\b',
    'RaGOO': r'\bragoo\b', 'RagTag': r'\bragtag\b', 'DBG2OLC': r'\bdbg2olc\b',
    'Velvet': r'\bvelvet\b', 'Maia': r'\bmaia\b', 'Sprai': r'\bsprai\b',
    'Flye': r'\bflye\b', 'Hifiasm': r'\bhifiasm\b', 'FALCON': r'\bfalcon(?:-unzip)?\b',
    'wtdbg2': r'\b(?:wtdbg2|redbean)\b', 'SOAPdenovo': r'\bsoapdenovo\w*',
    'ALLPATHS-LG': r'\ballpaths(?:-lg)?\b', 'Celera/CABOG': r'\b(?:celera|cabog)\b',
    'Unicycler': r'\bunicycler\b', 'MaSuRCA': r'\bmasurca\b', 'Raven': r'\braven\b',
    'IDBA': r'\bidba(?:-ud)?\b', 'Minia': r'\bminia\b', 'Edena': r'\bedena\b',
    'A5': r'\ba5(?:-miseq)?\b', 'Shasta': r'\bshasta\b', '3D-DNA': r'\b3d-dna\b',
    'SALSA': r'\bsalsa\w*', 'Pilon': r'\bpilon\b', 'Racon': r'\bracon\b',
    'Medaka': r'\bmedaka\b', 'Quickmerge': r'\bquickmerge\b', 'PBcR': r'\bpbcr\b',
}
STRUCTURAL = {
    'AUGUSTUS': r'\baugustus\b', 'GeneMark': r'\bgenemark(?:[.-]\w+)*',
    'BRAKER': r'\bbraker\d*\b', 'MAKER': r'\bmaker\d*\b', 'SNAP': r'\bsnap\b',
    'GlimmerHMM': r'\bglimmerhmm\b', 'Gnomon': r'\bgnomon\b',
    'EGAP': r'\b(?:egap|ncbi eukaryotic genome annotation pipeline)\b',
    'PGAP': r'\b(?:pgap|ncbi prokaryotic genome annotation pipeline)\b',
    'Funannotate': r'\bfunannotate\b', 'PASA': r'\bpasa\b',
    'EVidenceModeler': r'\b(?:evidencemodeler|evidence ?modeler|evm)\b',
    'Liftoff': r'\bliftoff\b', 'Liftover': r'\bliftover\b', 'RATT': r'\bratt\b',
    'YGAP': r'\b(?:ygap|yeast genome annotation pipeline)\b', 'AGAPE': r'\bagape\b',
    'Exonerate': r'\bexonerate\b', 'GeneWise': r'\bgenewise\b',
    'Prodigal': r'\bprodigal\b', 'MITOS': r'\bmitos\d*\b', 'MFannot': r'\bmfannot\b',
    'TransDecoder': r'\btransdecoder\b', 'GeMoMa': r'\bgemoma\b',
}


def clean(value):
    value = re.sub(r'\s+', ' ', value or '').strip()
    return NA if value.lower() in ('', 'na', 'n/a', 'unknown', 'not provided', 'not available') else value


def structured_header_fields(header):
    result = {}
    for line in header.splitlines():
        match = re.match(r'^\s+([^:]+?)\s*::\s*(.*?)\s*$', line)
        if match:
            result[match[1].strip()] = match[2].strip()
    return result


def technology(raw, assembly_raw):
    text = raw.lower() if raw != NA else ''
    platforms = []
    # Solexa and the literal llumina typo are explicit aliases, not guesses.
    if re.search(r'illumina|\bllumina\b|solexa', text):
        platforms.append('Illumina')
    if re.search(r'nanopore|\bont\b|mini?ion|promethion|gridion', text):
        platforms.append('Nanopore')
    pacbio = bool(re.search(r'pacbio|pacific biosciences|\bsmrt\b|\bhifi\b|\bccs\b', text))
    hifi = bool(re.search(r'\bhifi\b|hi-fi|\bccs\b|circular consensus', text))
    clr = bool(re.search(r'\bclr\b|continuous long reads?', text))
    if pacbio:
        platforms.append('PacBio HiFi/CCS' if hifi else ('PacBio CLR' if clr else 'PacBio (mode unspecified)'))
    if re.search(r'\b454\b|roche', text):
        platforms.append('454')
    if re.search(r'sanger|\babi\s*37[03]0|\b3730xl\b', text):
        platforms.append('Sanger/ABI capillary')
    if re.search(r'dnbseq|\bmgi\b', text):
        platforms.append('MGI/DNBSEQ')
    if re.search(r'ion torrent|ion proton|ion pgm|ion s5', text):
        platforms.append('Ion Torrent')
    if re.search(r'solid', text):
        platforms.append('SOLiD')
    hic = bool(re.search(r'hi[- ]?c\b', text + ' ' + (assembly_raw if assembly_raw != NA else '').lower()))
    if not platforms:
        platforms = [NA if raw == NA else 'Other: ' + raw]
    group = ' + '.join(platforms)
    if hic:
        group += ' + Hi-C'
    return dict(sequencing_platforms='; '.join(platforms), sequencing_group=group,
                hic_reported='Yes' if hic else NA,
                pacbio_hifi_reported='Yes' if hifi else NA,
                pacbio_mode=('HiFi/CCS' if hifi else ('CLR' if clr else 'Unspecified')) if pacbio else NA)


def extract(row):
    fields = row['report_fields']
    gb = structured_header_fields(row['genbank_first_record_header'])
    raw_technology = clean(fields.get('Sequencing technology') or gb.get('Sequencing Technology'))
    raw_assembly = clean(fields.get('Assembly method') or gb.get('Assembly Method'))
    assembly_tools = [name for name, pattern in ASSEMBLERS.items() if re.search(pattern, raw_assembly, re.I)]
    if not assembly_tools and raw_assembly != NA:
        # Keep unknown reported method available rather than silently dropping it.
        assembly_tools = ['Other reported method: ' + raw_assembly]
    annotation_parts = []
    for key in ('Annotation Pipeline', 'Annotation Method', 'Annotation Software Version'):
        if key in gb:
            annotation_parts.append(key + ': ' + gb[key])
    evidence = list(annotation_parts)
    evidence += [s for s in row['gff_structural_feature_sources'] if s.lower() not in
                 ('genbank', 'embl', 'ddbj', 'tpg', '.', 'refseq', 'bestrefseq', 'sgd', 'ncbi')]
    evidence += row.get('gff_structural_method_evidence', [])
    evidence += [line for line in row['gff_header'] if re.search(r'annotation|gene.?predict|pipeline', line, re.I)]
    comment = re.search(r'^COMMENT\b.*', row['genbank_first_record_header'], re.M | re.S)
    if comment:
        for line in comment[0].splitlines():
            if any(re.search(pattern, line, re.I) for pattern in STRUCTURAL.values()):
                evidence.append(line.strip())
    combined = '\n'.join(evidence)
    annotation_tools = [name for name, pattern in STRUCTURAL.items() if re.search(pattern, combined, re.I)]
    annotation_versions = []
    for name in annotation_tools:
        for match in re.finditer(STRUCTURAL[name] + r'\s*(?:ver\.?|version|v\.?)\s*([0-9][\w.\-]*)', combined, re.I):
            annotation_versions.append(name + ' ' + match[1].rstrip('.'))
    # Only structural metadata is inspected. No product names, GO, Pfam, BLAST,
    # InterPro or other functional annotation attributes are used.
    return dict(sequencing_technology_raw=raw_technology, **technology(raw_technology, raw_assembly),
                assembly_method_raw=raw_assembly, assembly_tools='; '.join(assembly_tools) or NA,
                genome_coverage_raw=clean(fields.get('Genome coverage') or gb.get('Genome Coverage')),
                structural_annotation_tools='; '.join(annotation_tools) or NA,
                structural_annotation_tool_versions='; '.join(sorted(set(annotation_versions))) or NA,
                structural_annotation_method_raw=clean(gb.get('Annotation Method')),
                structural_annotation_pipeline_raw=clean(gb.get('Annotation Pipeline')),
                structural_annotation_software_version_raw=clean(gb.get('Annotation Software Version')),
                structural_annotation_evidence=clean(' | '.join(evidence)) if annotation_tools else NA,
                gff_structural_feature_sources='; '.join(row['gff_structural_feature_sources']) or NA,
                assembly_report_url=row['report_url'], annotation_gff_url=row['gff_url'],
                annotation_genbank_url=row['gbff_url'],
                method_source_scope='NCBI assembly report; full GFF structural feature sources/inference; first GenBank record header',
                missing_method_meaning='N/A = not explicit in checked assembly sources, not proven absent')


def write_csv(path, rows):
    with path.open('w', newline='', encoding='utf-8-sig') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    raw = json.loads((RAW / 'raw_methods.json').read_text())
    if len(raw) != 435 or any(r['status'] != 'ok' for r in raw):
        raise ValueError('Incomplete source audit')
    methods = {r['assembly_accession']: extract(r) for r in raw}
    base = json.loads(BASE.read_text())
    rows = [dict(r, **methods[r['assembly_accession']]) for r in base['rows']]
    full = [r for r in rows if r['included_in_primary_plots']]
    groups = sorted({r['sequencing_group'] for r in full}, key=lambda g: (-sum(r['sequencing_group'] == g for r in full), g))
    metrics = [('genome_size_mb', 'Genome FASTA size (Mb)'), ('gc_percent', 'GC content (%)'),
               ('number_of_proteins', 'Protein FASTA records'), ('contig_count', 'Contig count (log scale)')]
    summaries = []
    for group in groups:
        subset = [r for r in full if r['sequencing_group'] == group]
        record = dict(sequencing_group=group, n=len(subset))
        for key, _ in metrics:
            record[key + '_median'] = float(np.median([r[key] for r in subset]))
            record[key + '_q1'] = float(np.quantile([r[key] for r in subset], .25))
            record[key + '_q3'] = float(np.quantile([r[key] for r in subset], .75))
        summaries.append(record)
    counts = dict(sequencing_groups_all=dict(Counter(r['sequencing_group'] for r in rows)),
                  sequencing_missing=sum(r['sequencing_technology_raw'] == NA for r in rows),
                  hic_explicit=sum(r['hic_reported'] == 'Yes' for r in rows),
                  hifi_explicit=sum(r['pacbio_hifi_reported'] == 'Yes' for r in rows),
                  assembly_method_missing=sum(r['assembly_method_raw'] == NA for r in rows),
                  structural_tool_missing=sum(r['structural_annotation_tools'] == NA for r in rows),
                  structural_tools=dict(Counter(tool for r in rows if r['structural_annotation_tools'] != NA
                                                for tool in r['structural_annotation_tools'].split('; '))),
                  assembly_tools=dict(Counter(tool for r in rows if r['assembly_tools'] != NA for tool in r['assembly_tools'].split('; '))))
    notes = [
        '435 assembly records. Comparison plots use the same 434 Full-representation records as the original analysis.',
        'Raw NCBI sequencing/assembly labels are retained. N/A means not explicit in the checked sources, not confirmed absence.',
        'Hi-C is a library/contact method, recorded separately from platforms. Instrument names alone do not imply Hi-C.',
        'PacBio alone, RS II, Sequel, or Revio does not imply HiFi. HiFi/CCS requires explicit wording in sequencing metadata.',
        'Solexa and the explicit llumina spelling variant are normalized to Illumina; raw text is unchanged.',
        'Groups preserve all reported platforms. Platform-specific counts overlap for hybrid assemblies; group counts do not.',
        'Assembly method may list assemblers, scaffolders or polishers. Tool names do not imply pipeline order or role.',
        'Structural tools use explicit annotation fields or gene/CDS/transcript/exon source and ab-initio inference; provider names are not tools.',
        'No functional annotation tools are assigned. Gene products, GO terms, Pfam, BLAST/InterPro/eggNOG annotations are excluded.',
        'Sources: complete assembly reports and GFF structural metadata; GenBank inspection is limited to the first record header.',
        'SRA assays are not automatically attributed to an assembly: a shared BioSample alone does not establish reads used for assembly.',
        'No technologies/tools are inferred from publications, sequencing instruments, assembler requirements or annotation providers.',
        'Plots are descriptive. Sequencing groups are confounded with submitter, annotation pipeline and sample selection.',
    ]
    payload = dict(rows=rows, methods=methods, summary=summaries, counts=counts, notes=notes)
    (OUT / 'methods_analysis.json').write_text(json.dumps(payload, indent=2) + '\n')
    write_csv(OUT / 'genome_metadata_with_methods.csv', rows)
    write_csv(OUT / 'sequencing_group_summary.csv', summaries)
    write_csv(OUT / 'structural_annotation_tools.csv', [
        {key:r[key] for key in ('assembly_accession','strain','structural_annotation_tools',
                                'structural_annotation_tool_versions','structural_annotation_pipeline_raw',
                                'structural_annotation_method_raw','structural_annotation_evidence',
                                'annotation_gff_url','annotation_genbank_url')}
        for r in rows])
    write_csv(OUT / 'assembly_tool_counts.csv', [dict(reported_tool=name,assemblies=count)
              for name,count in sorted(counts['assembly_tools'].items(),key=lambda x:-x[1])])
    # Horizontal facets keep explicit hybrid-technology names readable.
    height = max(10.4, len(groups) * 1.0 + 3)
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'pdf.fonttype':42,
                         'axes.spines.top':False,'axes.spines.right':False})
    fig, axes = plt.subplots(2, 2, figsize=(16, height))
    rng = np.random.default_rng(4932)
    colors = plt.get_cmap('tab20').colors
    for ax, (metric, title) in zip(axes.flat, metrics):
        values = [[r[metric] for r in full if r['sequencing_group'] == group] for group in groups]
        bp = ax.boxplot(values, orientation='horizontal', showfliers=False, patch_artist=True,
                        widths=.55, whis=1.5, medianprops={'color':'#182634','linewidth':1.5})
        for i, (v, box) in enumerate(zip(values, bp['boxes'])):
            color = colors[i % len(colors)]
            box.set_facecolor(color)
            box.set_alpha(.35)
            ax.scatter(v, i + 1 + rng.uniform(-.15,.15,len(v)), color=color, s=10, alpha=.42, edgecolors='none')
        labels = [g.replace('PacBio (mode unspecified)', 'PacBio (unspecified)').replace(' + ', '\n+ ') + '\n(n=' + str(len(v)) + ')'
                  for g, v in zip(groups, values)]
        ax.set_yticks(range(1,len(groups)+1), labels, fontsize=9)
        ax.invert_yaxis()
        ax.set_title(title, loc='left', fontweight='bold', fontsize=13, pad=12)
        ax.set_xlabel(title)
        ax.set_ylabel('Reported sequencing group')
        ax.grid(axis='x', color='#E1E5E9', linewidth=.7)
        ax.set_axisbelow(True)
        if metric == 'contig_count':
            ax.set_xscale('log')
        if metric == 'number_of_proteins':
            ax.xaxis.set_major_formatter(StrMethodFormatter('{x:,.0f}'))
    fig.suptitle('S. cerevisiae assemblies by reported sequencing technology', x=.17, ha='left', fontsize=17, fontweight='bold', y=.985)
    fig.text(.17,.957,'434 Full-representation assemblies; explicit NCBI assembly metadata; unreported values are N/A',fontsize=11)
    fig.text(.17,.025,'Boxes: Q1–Q3; line: median; whiskers: 1.5 × IQR. Dots: all assemblies. No outlier removal.\n'
             'Hi-C is recorded separately; PacBio is not labeled HiFi unless explicitly stated. Associations are descriptive.',fontsize=9,linespacing=1.5)
    fig.subplots_adjust(left=.17,right=.985,top=.90,bottom=.10,hspace=.25,wspace=.65)
    fig.savefig(OUT / 'sequencing_technology_boxplots.png',dpi=180,facecolor='white')
    fig.savefig(OUT / 'sequencing_technology_boxplots.pdf',facecolor='white')
    plt.close(fig)
    print(json.dumps(counts,indent=2))


if __name__ == '__main__':
    main()
