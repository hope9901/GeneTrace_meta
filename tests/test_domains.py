import unittest
from genetrace.domains import six_frame_orfs,domain_coordinates,genomic_assessment,parse_domtbl

class DomainEvidenceGuards(unittest.TestCase):
 def test_plus_and_minus_domain_coordinates(self):
  plus=dict(seqid='x',start=10,end=39,strand='+',length=10)
  minus=dict(seqid='x',start=10,end=39,strand='-',length=10)
  self.assertEqual(domain_coordinates(plus,2,4)['start'],13)
  self.assertEqual(domain_coordinates(minus,2,4)['start'],28)
  self.assertEqual(domain_coordinates(minus,2,4)['end'],36)
 def test_stop_separated_orfs_preserve_genomic_coordinates(self):
  peptides,mapping=six_frame_orfs({'x':'ATGAAATAAATGGGGTAA'},min_aa=2)
  plus=[(peptides[k],m['start'],m['end']) for k,m in mapping.items() if m['strand']=='+' and m['frame']==1]
  self.assertIn(('MK',1,6),plus);self.assertIn(('MG',10,15),plus)
 def test_pseudogene_is_not_an_annotation_gap(self):
  h=dict(profile='PF00005',domain_score=100,profile_coverage=1,coverage_class='near_complete');pos=dict(seqid='x',start=10,end=30,strand='+')
  fs=[dict(type='pseudogene',seqid='x',start=1,end=50,strand='+',attrs={'pseudo':'true','gene':'ycjV'})]
  r=genomic_assessment(h,pos,fs,[])
  self.assertEqual(r['classification'],'pseudogene_or_disrupted_model_unresolved');self.assertFalse(r['biological_loss_inferred'])
 def test_domain_hit_without_cds_remains_candidate(self):
  h=dict(profile='PF00005',domain_score=100,profile_coverage=1,coverage_class='near_complete');pos=dict(seqid='x',start=10,end=30,strand='+')
  self.assertEqual(genomic_assessment(h,pos,[],[])['classification'],'unannotated_genomic_domain_candidate')
 def test_domain_and_sequence_gathering_scores_both_required(self):
  m=[dict(accession='PF00005.1',ga='20 10;')]
  row='p - 100 ABC_tran PF00005.1 100 1e-20 25 0 1 1 1e-20 1e-20 9 0 1 100 1 100 1 100 1 desc'
  self.assertEqual(parse_domtbl(row,m),[])
 def test_partial_domain_is_retained_and_labeled(self):
  m=[dict(accession='PF00005.1',ga='20 10;')]
  row='p - 100 ABC_tran PF00005.1 100 1e-20 25 0 1 1 1e-20 1e-20 15 0 1 50 1 50 1 50 1 desc'
  self.assertEqual(parse_domtbl(row,m)[0]['coverage_class'],'partial')
