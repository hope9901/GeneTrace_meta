import importlib.util
import json
import unittest
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]


def module(name):
    spec=importlib.util.spec_from_file_location(name,ROOT/'scripts'/(name+'.py'))
    result=importlib.util.module_from_spec(spec);spec.loader.exec_module(result)
    return result


class CohortGuards(unittest.TestCase):
    def test_unknown_strains_never_have_a_merge_key(self):
        p=module('prepare_orthology_cohort')
        for value in (None,'','N/A','na','unknown','missing'):
            self.assertEqual(p.normalized(value),'')

    def test_derivative_strain_names_remain_distinct(self):
        p=module('prepare_orthology_cohort')
        self.assertNotEqual(p.normalized('BY4741'),p.normalized('BY4742'))
        self.assertNotEqual(p.normalized('Awamori-101'),p.normalized('Awamori 101-18'))
        self.assertEqual(p.normalized(' S288C '),p.normalized('s288c'))

    def test_reference_and_full_representation_win(self):
        p=module('prepare_orthology_cohort')
        base=dict(paired_refseq=None,genome_representation='Full',number_of_proteins=6000,
                  assembly_level='Contig',n_fraction=0,fasta_sequence_n50_bp=100000,
                  release_date='2026-01-01',assembly_accession='GCA_x')
        reference=dict(base,paired_refseq='GCF_000146045.2')
        self.assertGreater(p.rank(reference),p.rank(base))
        partial=dict(base,genome_representation='Partial',fasta_sequence_n50_bp=9999999)
        self.assertGreater(p.rank(base),p.rank(partial))

    def test_actual_cohort_accounting_and_sparse_genome_retention(self):
        out=ROOT/'outputs/orthology_20261030'
        if not (out/'all_gene_input_summary.json').exists():
            self.skipTest('Prepare the full cohort first')
        decisions=json.loads((out/'strain_deduplication.json').read_text())
        selected=json.loads((out/'selected_assemblies.json').read_text())
        inputs=json.loads((out/'all_gene_input_summary.json').read_text())
        ids={r['assembly_accession'] for r in selected}
        self.assertEqual(len(decisions),435)
        self.assertEqual(len(ids),279)
        self.assertEqual(sum(r['selected'] for r in decisions),len(ids))
        self.assertTrue(all(r['representative'] in ids for r in decisions))
        self.assertIn('GCA_903819185.2',ids)
        self.assertEqual({r['assembly_accession'] for r in inputs},ids)
        self.assertEqual(sum(r['orthology_coding_units'] for r in inputs),1541009)
        self.assertEqual(sum(r['unmapped_protein_units'] for r in inputs),0)


if __name__=='__main__':
    unittest.main()
