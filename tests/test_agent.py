import unittest
from unittest.mock import patch
from genetrace.agent import assessment,validate_model_step,resolve_gene,cohort,Investigation

class AgentEvidenceGuards(unittest.TestCase):
    def test_invalid_model_tool_is_rejected(self):
        with self.assertRaises(ValueError):validate_model_step({'next_tool':'execute_shell'},['E001'],['inspect_sequence'])
    def test_invented_model_citation_is_rejected(self):
        with self.assertRaises(ValueError):validate_model_step({'next_tool':None,'hypotheses':[{'text':'claim','evidence_ids':['E999']}]},['E001'],[])
    def test_uncited_model_claim_is_rejected(self):
        with self.assertRaises(ValueError):validate_model_step({'hypotheses':[{'text':'claim','evidence_ids':[]}]},[],[])
    def test_negative_search_is_not_gene_loss(self):
        r=assessment([dict(id='E001',tool='targeted_homology',status='ok',result=dict(best_genome_hit=None))])
        self.assertEqual(r['code'],'inconclusive');self.assertFalse(r['biological_loss_inferred'])
    def test_failed_annotation_cannot_be_an_annotation_gap(self):
        r=assessment([dict(id='E001',tool='targeted_homology',status='ok',result=dict(best_genome_hit={'start':1},protein_hit_count=0)),dict(id='E002',tool='inspect_annotation',status='error',result={}),dict(id='E003',tool='inspect_sequence',status='ok',result={})])
        self.assertEqual(r['code'],'inconclusive')
    def test_recorded_frameshift_is_not_an_annotation_omission(self):
        evidence=[dict(id='E001',tool='targeted_homology',status='ok',result=dict(best_genome_hit={'start':1},protein_hit_count=0)),dict(id='E002',tool='inspect_annotation',status='ok',result=dict(recorded_disruption=True,overlapping_cds=[])),dict(id='E003',tool='inspect_sequence',status='ok',result={})]
        self.assertEqual(assessment(evidence)['code'],'disrupted_model_unresolved')
    def test_scope_mismatch_is_not_individual_gene_loss(self):
        evidence=[dict(id='E001',tool='targeted_homology',status='ok',result=dict(best_genome_hit={'start':1},protein_hit_count=0)),dict(id='E002',tool='inspect_annotation',status='ok',result=dict(annotation_scope_mismatch=True)),dict(id='E003',tool='inspect_sequence',status='ok',result={})]
        self.assertEqual(assessment(evidence)['code'],'annotation_scope_mismatch')
    def test_all_cohort_and_reference_genes_are_available(self):
        self.assertEqual(sum(k.startswith('GCA_') for k in cohort()),279)
        self.assertEqual(sum(k.startswith('GCF_') for k in cohort()),100)
        self.assertEqual(resolve_gene('rpoB','GCF_000005845.2')['gene_name'],'rpoB')
        self.assertEqual(resolve_gene('YKL094W')['gene_name'],'YJU3')
        self.assertEqual(resolve_gene('ACT1')['gene_name'],'ACT1')
    def test_path_and_command_strings_cannot_be_gene_queries(self):
        for value in ['../../.env','ACT1; echo secret','$(env)']:
            with self.assertRaises(ValueError):resolve_gene(value)
