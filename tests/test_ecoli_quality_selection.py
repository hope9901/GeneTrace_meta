import importlib.util,pathlib,unittest
p=pathlib.Path(__file__).resolve().parents[1]/'scripts/prepare_ecoli_quality100.py'
spec=importlib.util.spec_from_file_location('ecoli_quality',p);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
class QualitySelectionGuards(unittest.TestCase):
 def record(self):return dict(accession='GCF_1.1',assembly_info={'assembly_level':'Complete Genome'},assembly_stats={'total_sequence_length':5000000,'total_ungapped_length':5000000,'contig_n50':4500000},checkm_info={'completeness':99.5,'contamination':.1},annotation_info={'stats':{'gene_counts':{'protein_coding':4300}}},average_nucleotide_identity={'taxonomy_check_status':'OK'})
 def test_missing_quality_is_not_good_quality(self):
  r=self.record();r['checkm_info']={};self.assertFalse(m.eligible(r))
 def test_contamination_excludes_candidate(self):
  r=self.record();r['checkm_info']['contamination']=2;self.assertFalse(m.eligible(r))
 def test_noncomplete_or_unannotated_excluded(self):
  r=self.record();r['assembly_info']['assembly_level']='Contig';self.assertFalse(m.eligible(r))
  r=self.record();r['annotation_info']={};self.assertFalse(m.eligible(r))
 def test_rank_prefers_less_contamination(self):
  a=self.record();b=self.record();b['checkm_info']['contamination']=.2
  self.assertLess(m.quality(a),m.quality(b));self.assertTrue(m.eligible(a))
