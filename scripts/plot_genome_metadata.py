"""Generate reproducible assembly-level descriptive plots and research tables."""
import argparse
import csv
import json
import os
from pathlib import Path

os.environ.setdefault('MPLCONFIGDIR', '/tmp/genetrace-matplotlib')
os.environ.setdefault('XDG_CACHE_HOME', '/tmp/genetrace-cache')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import StrMethodFormatter
import numpy as np

LEVELS = ['Complete Genome', 'Chromosome', 'Scaffold', 'Contig']
COLORS = ['#355C7D', '#438A83', '#B6873D', '#8064A2']
METRICS = [('genome_size_mb', 'Genome FASTA size (Mb)'), ('gc_percent', 'GC content (%)'),
           ('number_of_proteins', 'Protein FASTA records'), ('contig_count', 'Number of contigs (log scale)')]


def summary(rows):
    result = []
    for level in ['All full genomes'] + LEVELS:
        group = rows if level == 'All full genomes' else [r for r in rows if r['assembly_level'] == level]
        record = dict(assembly_level=level, n=len(group))
        for key, _ in METRICS:
            values = np.array([r[key] for r in group if r[key] is not None], dtype=float)
            if not len(values):
                raise ValueError('Metric has no data: ' + key)
            for name, value in [('min', np.min(values)), ('q1', np.quantile(values, .25)),
                                ('median', np.median(values)), ('q3', np.quantile(values, .75)),
                                ('max', np.max(values)), ('mean', np.mean(values))]:
                record[key + '_' + name] = float(value)
            record[key + '_missing'] = len(group) - len(values)
        result.append(record)
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', required=True)
    parser.add_argument('--output', required=True)
    args = parser.parse_args()
    rows = json.loads(Path(args.input).read_text())
    if len(rows) != 435 or len({r['assembly_accession'] for r in rows}) != len(rows):
        raise ValueError('Unexpected cohort size or duplicate accessions')
    full = [r for r in rows if r['included_in_primary_plots']]
    if len(full) != 434:
        raise ValueError('Unexpected primary plot cohort')
    for r in rows:
        if not r.get('gc_method') or 'ncbi_gc_percent' not in r:
            raise ValueError('Recompute GC with collect_genome_metadata.py --recompute-gc before plotting')
        for key in ('genome_size_bp', 'gc_percent', 'number_of_proteins', 'contig_count'):
            if r[key] is None:
                raise ValueError('Missing requested metric: ' + r['assembly_accession'] + '/' + key)
        if not (0 <= r['gc_percent'] <= 100) or min(r['genome_size_bp'], r['number_of_proteins'], r['contig_count']) <= 0:
            raise ValueError('Metric out of bounds')
    output = Path(args.output)
    output.mkdir(parents=True, exist_ok=True)
    stats = summary(full)
    with (output / 'assembly_level_summary.csv').open('w', newline='', encoding='utf-8-sig') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(stats[0]))
        writer.writeheader()
        writer.writerows(stats)
    with (output / 'genome_metadata.csv').open('w', newline='', encoding='utf-8-sig') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10,
                         'axes.titlesize': 13, 'axes.labelsize': 11, 'pdf.fonttype': 42,
                         'axes.spines.top': False, 'axes.spines.right': False})
    fig, axes = plt.subplots(2, 2, figsize=(13.4, 9.8))
    rng = np.random.default_rng(4932)
    for panel, (ax, (key, label)) in enumerate(zip(axes.flat, METRICS)):
        groups = [[r[key] for r in full if r['assembly_level'] == level] for level in LEVELS]
        boxes = ax.boxplot(groups, positions=range(1, 5), widths=.52, patch_artist=True,
                           showfliers=False, whis=1.5,
                           medianprops={'color': '#182634', 'linewidth': 1.8},
                           whiskerprops={'color': '#52606D'}, capprops={'color': '#52606D'})
        for idx, (group, color, box) in enumerate(zip(groups, COLORS, boxes['boxes']), 1):
            box.set_facecolor(color)
            box.set_alpha(.32)
            box.set_edgecolor(color)
            ax.scatter(idx + rng.uniform(-.16, .16, len(group)), group, s=11,
                       color=color, alpha=.42, edgecolors='none', rasterized=True, zorder=2)
        ax.set_xticks(range(1, 5), [level.replace('Complete Genome', 'Complete\nGenome') + '\n(n=' + str(len(g)) + ')'
                                  for level, g in zip(LEVELS, groups)])
        ax.set_title(chr(65 + panel) + '   ' + label, loc='left', fontweight='bold', pad=12)
        ax.set_xlabel('Assembly level')
        ax.set_ylabel(label)
        ax.grid(axis='y', color='#E1E5E9', linewidth=.7)
        ax.set_axisbelow(True)
        ax.set_xlim(.5, 4.5)
        ax.margins(y=.08)
        if key == 'contig_count':
            ax.set_yscale('log')
        elif key == 'number_of_proteins':
            ax.yaxis.set_major_formatter(StrMethodFormatter('{x:,.0f}'))
    fig.suptitle('Saccharomyces cerevisiae: annotated NCBI assemblies', x=.075, ha='left',
                 fontsize=18, fontweight='bold', y=.98)
    fig.text(.075, .936, '434 full-representation assemblies | Latest-version snapshot: 30 Oct 2026',
             fontsize=11, color='#465361')
    fig.text(.075, .026, 'Boxes: Q1–Q3; center line: median; whiskers: 1.5 × IQR. Dots: all assemblies, including outliers.\n'
             'Protein counts and GC%: calculated from MD5-verified FASTA files. GC% excludes ambiguous bases.\n'
             'One partial assembly is retained in the table but excluded here. Full representation does not imply complete annotation.',
             fontsize=9, color='#465361', linespacing=1.6)
    fig.subplots_adjust(left=.075, right=.975, bottom=.20, top=.865, hspace=.53, wspace=.28)
    fig.savefig(output / 'assembly_metadata_boxplots.png', dpi=200, facecolor='white')
    fig.savefig(output / 'assembly_metadata_boxplots.pdf', facecolor='white')
    plt.close(fig)
    differences = [r for r in rows if r['number_of_proteins'] != r['protein_coding_gene_count']]
    notes = [
        'Snapshot: 2026-10-30. One row per assembly; GenBank/RefSeq paired duplicates removed.',
        'All 435 assemblies are in the table; the primary plots and summaries use 434 Full-representation assemblies.',
        'Full representation differs from Complete Genome assembly level. Partial GCA_000091065.2 is excluded from plots.',
        'Plotted genome size is total length of all genomic FASTA records. GC% is 100*(G+C)/(A+C+G+T) over the same FASTA.',
        'NCBI summary genome size covers the primary assembly only. It is preserved separately; FASTA may include other units.',
        'Rounded NCBI GC% is preserved separately. N and other ambiguity codes are excluded from the recomputed GC denominator.',
        'Number of proteins is the actual protein.faa.gz record count, verified against the NCBI-published compressed-file MD5.',
        'Protein records are not necessarily unique protein sequences or protein-coding gene loci. Both counts are kept.',
        'Genome size is in bp and Mb (1 Mb = 1,000,000 bp). GC is stored on the 0–100 percent scale.',
        'No outlier removal or inferential significance test. Assemblies are not necessarily independent biological strains.',
        'Genome representation, assembly type and annotation provider may explain some heterogeneity; causality was not tested.',
        'All source assembly_type fields read haploid; this label is not an independent ploidy measurement.',
        'Contig count is NCBI primary-assembly metadata. FASTA sequence count/N50 are separate and may describe scaffolds.',
        'GCA_903819185.2 has 16 protein FASTA records despite Full representation. File availability does not imply comprehensive annotation.',
        'Boxplots are static images; regenerate the analysis to reflect changes in the metadata table.',
    ]
    payload = dict(rows=rows, summary=stats, notes=notes,
                   protein_gene_count_differences=len(differences),
                   protein_download_bytes=sum(json.loads(p.read_text())['compressed_bytes'] for p in
                       Path(args.input).parent.joinpath('protein_counts').glob('*.json')))
    (output / 'analysis.json').write_text(json.dumps(payload, indent=2) + '\n')
    report = ['# S. cerevisiae assembly metadata', '',
              '435 assemblies tabulated; 434 Full-representation assemblies plotted.', '',
              '| Assembly level | n | Median size (Mb) | Median GC (%) | Median proteins |',
              '|---|---:|---:|---:|---:|']
    for s in stats:
        report.append('| %s | %d | %.3f | %.2f | %s |' %
                      (s['assembly_level'], s['n'], s['genome_size_mb_median'], s['gc_percent_median'],
                       format(s['number_of_proteins_median'], ',.1f')))
    report += ['', 'Protein FASTA records differ from protein-coding gene counts in %d assemblies.' % len(differences),
               '', '## Definitions and scope', ''] + ['- ' + n for n in notes]
    (output / 'analysis_notes.md').write_text('\n'.join(report) + '\n')
    print(json.dumps(dict(summary=stats, protein_gene_count_differences=len(differences)), indent=2))


if __name__ == '__main__':
    main()
