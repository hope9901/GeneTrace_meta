import importlib.util
import json
import unittest
from pathlib import Path

from genetrace.orthology_evidence import scope_aware_presence

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'outputs/loss_artifact_tests_20261030'


class EvidenceGuards(unittest.TestCase):
    def test_nuclear_annotation_not_provided_is_not_a_zero_absence_call(self):
        self.assertEqual(scope_aware_presence(0,{'nuclear_or_unclassified'},{'nuclear_coding_loci':0}),'N/A')
        self.assertEqual(scope_aware_presence(1,{'nuclear_or_unclassified'},{'nuclear_coding_loci':0}),1)

    def test_missing_organelle_sequence_is_not_organelle_gene_loss(self):
        self.assertEqual(scope_aware_presence(0,{'organelle'},{'explicit_organelle_bp':0}),'N/A')
        self.assertEqual(scope_aware_presence(1,{'organelle'},{'explicit_organelle_bp':0}),1)

    def test_real_scope_audit_and_search_controls(self):
        if not (OUT/'adjudicated_cases.json').exists():self.skipTest('Run natural-case diagnostics first')
        audit={r['assembly_accession']:r for r in json.loads((OUT/'annotation_scope_audit.json').read_text())}
        self.assertEqual(len(audit),279)
        case=audit['GCA_903819185.2']
        self.assertEqual(case['nuclear_coding_loci'],0)
        self.assertEqual(case['organelle_coding_loci'],16)
        self.assertGreater(case['nuclear_or_unclassified_bp'],11_000_000)
        reports=json.loads((OUT/'adjudicated_cases.json').read_text())
        by_assembly={r['assembly_accession']:r for r in reports['cases']}
        self.assertEqual(by_assembly['GCA_000146045.2']['protein_positive'],64)
        self.assertEqual(by_assembly['GCA_000146045.2']['full_genome_positive'],64)
        self.assertEqual(by_assembly['GCA_903819185.2']['protein_positive'],0)
        self.assertGreaterEqual(by_assembly['GCA_903819185.2']['full_genome_positive'],60)
        self.assertTrue(all(r['biological_loss_calls']==0 for r in reports['cases']))
        self.assertTrue(all(not r['biological_loss_inferred'] for r in reports['loci']))
        swc3=next(r for r in reports['loci'] if r['target_assembly']=='GCA_000190255.1' and r['reference_gene_name']=='SWC3')
        self.assertEqual(swc3['classification'],'disrupted_or_pseudogene_model_unresolved')
        self.assertIn('frameshift',swc3['overlapping_gene_notes'])
        self.assertGreater(swc3['overlapping_sequence_alterations'],0)

    def test_flank_order_rejects_distant_paralog(self):
        spec=importlib.util.spec_from_file_location('syn',ROOT/'scripts/test_candidate_synteny.py')
        module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
        center=dict(seqid='chr1',start=10000,end=11000,strand='+')
        reference=dict(strand='+');flank=dict(strand='+')
        hit=dict(sseqid='chr1',start=7000,end=8000,strand='+',pident=99,query_coverage=.99,evalue=1e-100)
        self.assertTrue(module.compatible(hit,center,reference,flank,'left'))
        hit['sseqid']='chr2'
        self.assertFalse(module.compatible(hit,center,reference,flank,'left'))

    def test_same_sample_pairs_are_not_counted_twice(self):
        import csv
        file=OUT/'same_biosample_cross_level_pairs.csv'
        if not file.exists():self.skipTest('Run paired-assembly comparison first')
        with file.open(encoding='utf-8-sig') as handle:rows=list(csv.DictReader(handle))
        keys={(r['biosample'],r['level_a'],r['level_b']) for r in rows}
        self.assertEqual(len(keys),len(rows))
        self.assertEqual(len(rows),23)
        self.assertTrue(all(r['same_reported_sequencing']=='True' and r['same_annotation_provider']=='True' for r in rows))

    def test_synteny_does_not_override_recorded_disruption(self):
        file=OUT/'synteny_candidate_tests.json'
        if not file.exists():self.skipTest('Run neighborhood tests first')
        rows=json.loads(file.read_text())
        swc3=next(r for r in rows if r['reference_gene']=='SWC3')
        self.assertTrue(swc3['both_flanks_support'])
        self.assertEqual(swc3['source_classification'],'disrupted_or_pseudogene_model_unresolved')
        self.assertFalse(swc3['functional_gene_loss_inferred'])
        yju3=next(r for r in rows if r['reference_gene']=='YJU3')
        self.assertTrue(yju3['both_flanks_support'])


if __name__=='__main__':unittest.main()
