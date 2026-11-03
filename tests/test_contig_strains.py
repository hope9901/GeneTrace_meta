import importlib.util
import json
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location('contig_experiment', ROOT / 'scripts/test_contig_strains.py')
experiment = importlib.util.module_from_spec(spec)
spec.loader.exec_module(experiment)


class ContigExperimentGuards(unittest.TestCase):
    def test_n50_is_length_weighted(self):
        stats = experiment.assembly_metrics(dict(a='A' * 10, b='N' * 4, c='A' * 4))
        self.assertEqual(stats['n50_bp'], 10)
        self.assertEqual(stats['sequence_count'], 3)
        self.assertAlmostEqual(stats['n_fraction'], 4 / 18)

    def test_other_contig_or_strand_is_not_seeded_locus(self):
        locus = dict(seqid='ctg1', start=100, end=200, strand='+')
        hit = dict(sseqid='ctg2', start=100, end=200, strand='+')
        self.assertFalse(experiment.overlap(hit, locus))
        hit['sseqid'], hit['strand'] = 'ctg1', '-'
        self.assertFalse(experiment.overlap(hit, locus))

    def test_recorded_actual_blast_results(self):
        result = ROOT / 'runs/contig-strains-001/comparison.json'
        if not result.exists():
            self.skipTest('Run scripts/test_contig_strains.py first')
        summaries = json.loads(result.read_text())
        self.assertEqual({s['strain'] for s in summaries}, {'CLIB215', 'KACC30049'})
        for summary in summaries:
            with self.subTest(strain=summary['strain']):
                self.assertGreater(summary['metrics']['sequence_count'], 16)
                self.assertEqual(summary['best_hit']['pident'], 100.0)
                self.assertAlmostEqual(summary['best_hit']['query_coverage'], 0.992)
                cases = {s['scenario']: s for s in summary['scenarios']}
                self.assertEqual(cases['original']['verdict'], 'homolog_detected_annotation_unavailable')
                self.assertEqual(cases['locus_N_masked']['verdict'], 'assembly_limited')
                self.assertEqual(cases['locus_A_replaced']['verdict'], 'inconclusive')
                for scenario, case in cases.items():
                    self.assertEqual(case['controlled_perturbation'], scenario != 'original')
                    evidence = json.loads((result.parent / summary['strain'] / scenario / 'evidence.json').read_text())
                    self.assertTrue(all(e['status'] == 'ok' for e in evidence))
                    self.assertIn('tblastn', {e['tool'] for e in evidence})


if __name__ == '__main__':
    unittest.main()
