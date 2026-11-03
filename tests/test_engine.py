import json
import shutil
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from genetrace.engine import ROOT, classify, find_binary, investigate, validate_claims


class ScientificGuards(unittest.TestCase):
    def test_negative_search_is_not_gene_loss(self):
        verdict, _ = classify([], [], dict(seqid='chr', start=200, end=400, strand='+'),
                              dict(n_fraction=0, distance_to_contig_edge=200))
        self.assertEqual(verdict, 'inconclusive')

    def test_unknown_citation_rejected(self):
        with self.assertRaises(ValueError):
            validate_claims(dict(claims=[dict(evidence_ids=['E999'])]), [dict(id='E001')])


class PublicDataIntegration(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        if not (ROOT / 'data/demo/intact/manifest.json').exists():
            raise unittest.SkipTest('Run scripts/prepare_demo.py first')
        try:
            find_binary('blastp')
            find_binary('tblastn')
        except FileNotFoundError as error:
            raise unittest.SkipTest(str(error))

    def test_four_controls_with_real_blast(self):
        expected = dict(intact='annotated_homolog', annotation_removed='annotation_gap_candidate',
                        assembly_gap='assembly_limited', sequence_absent='inconclusive')
        for scenario, verdict in expected.items():
            with self.subTest(scenario=scenario), tempfile.TemporaryDirectory() as temp:
                output = Path(temp) / 'run'
                report = investigate(ROOT / 'data/demo' / scenario, output)
                self.assertEqual(report['verdict'], verdict)
                evidence = json.loads((output / 'evidence.json').read_text())
                self.assertTrue(all(e['status'] == 'ok' for e in evidence))
                validate_claims(report, evidence)
                if scenario == 'annotation_removed':
                    hit = next(e for e in evidence if e['tool'] == 'tblastn')['hits'][0]
                    self.assertEqual(hit['pident'], 100)
                    self.assertGreater(hit['query_coverage'], 0.98)
                    manifest = json.loads((output / 'input_manifest.json').read_text())
                    offset = manifest['reference_offset_0based']
                    # The gene span includes intron/stop bases absent from the protein HSP.
                    self.assertGreaterEqual(hit['start'] + offset, 53260)
                    self.assertLess(hit['start'] + offset, 53270)
                    self.assertLessEqual(hit['end'] + offset, 54696)
                    # Related actin proteins remain: family-wide absence must not be inferred.
                    protein_hits = next(e for e in evidence if e['tool'] == 'blastp')['hits']
                    self.assertTrue(any(h['passes_thresholds'] for h in protein_hits))

    def test_failed_tool_is_not_negative_evidence(self):
        with tempfile.TemporaryDirectory() as temp, patch('genetrace.engine.find_binary', side_effect=FileNotFoundError('missing')):
            report = investigate(ROOT / 'data/demo/intact', Path(temp) / 'run')
            self.assertEqual(report['verdict'], 'inconclusive')
            self.assertIn('failed', report['claims'][0]['text'])

    def test_modified_inputs_rejected(self):
        with tempfile.TemporaryDirectory() as temp:
            dataset = Path(temp) / 'dataset'
            shutil.copytree(ROOT / 'data/demo/intact', dataset)
            with (dataset / 'genome.fna').open('a') as handle:
                handle.write('A\n')
            with self.assertRaisesRegex(ValueError, 'checksum mismatch'):
                investigate(dataset, Path(temp) / 'run')

    def test_existing_output_is_preserved(self):
        with tempfile.TemporaryDirectory() as temp:
            with self.assertRaises(FileExistsError):
                investigate(ROOT / 'data/demo/intact', temp)


if __name__ == '__main__':
    unittest.main()
