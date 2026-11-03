import unittest
from genetrace.regional_qc import analyze_sam,cigar_ops

def line(name,cigar,seq,start=1,flag=0,mapq=60,quality=None):
 return '\t'.join([name,str(flag),'x',str(start),str(mapq),cigar,'*','0','0',seq,quality or 'D'*len(seq)])

class RegionalAlignmentGuards(unittest.TestCase):
 def test_counts_only_local_mismatches_not_whole_read_nm(self):
  seq='C'+'A'*19;s=line('r','20M',seq)+'\tNM:i:1'
  r=analyze_sam(s,'A'*30,'x',10,15)
  self.assertEqual(r['metrics']['mismatches'],0)
 def test_reverse_sam_sequence_is_not_complemented_twice(self):
  r=analyze_sam(line('r','10M','ACGTACGTAA',flag=16),'ACGTACGTAA','x',1,10)
  self.assertEqual(r['metrics']['mismatches'],0)
 def test_reference_skip_is_not_a_deletion_or_spanning_read(self):
  r=analyze_sam(line('r','5M100N5M','A'*10),'A'*120,'x',1,110)
  self.assertEqual(r['metrics']['deletion_bases'],0);self.assertEqual(r['metrics']['spanning_plus100_records'],0)
 def test_indel_fraction_uses_reads_callable_at_both_anchors(self):
  text='\n'.join([line('a','5M2D5M','A'*10),line('b','2D5M','A'*5,start=6)])
  r=analyze_sam(text,'A'*20,'x',1,12)
  d=next(e for e in r['indels'] if e['kind']=='D')
  self.assertEqual(d['raw_supporting_records'],2);self.assertEqual(d['supporting_records'],1);self.assertEqual(d['support_fraction'],1)
 def test_zero_quality_is_unresolved_not_sequencing_failure(self):
  r=analyze_sam(line('a','10M','A'*10,quality='!'*10),'A'*20,'x',1,10)
  self.assertTrue(r['quality_information_unresolved']);self.assertEqual(r['causal_diagnosis'],'unresolved')
 def test_low_mapq_and_supplementary_not_counted_as_primary_support(self):
  text='\n'.join([line('a','10M','A'*10,mapq=0),line('b','10M','A'*10,flag=2048)])
  r=analyze_sam(text,'A'*20,'x',1,10)
  self.assertEqual(r['metrics']['usable_primary_records'],0);self.assertEqual(r['metrics']['supplementary_records'],1)
 def test_one_strand_only_cannot_create_coherent_substitution(self):
  text='\n'.join(line(str(i),'10M','AAAACAAAAA') for i in range(20))
  r=analyze_sam(text,'A'*20,'x',1,10)
  self.assertEqual(r['coherent_substitutions'],[])
 def test_malformed_cigar_is_rejected(self):
  with self.assertRaises(ValueError):cigar_ops('10Mgarbage')
 def test_partial_data_never_establish_biological_loss(self):
  r=analyze_sam('','A'*20,'x',1,10,partial_input=True)
  self.assertFalse(r['biological_loss_inferred']);self.assertEqual(r['causal_diagnosis'],'unresolved')
 def test_bidirectional_consensus_is_disagreement_not_causal_diagnosis(self):
  text='\n'.join(line(str(i),'10M','AAAACAAAAA',flag=16 if i%2 else 0) for i in range(20))
  r=analyze_sam(text,'A'*20,'x',1,10)
  self.assertEqual(len(r['coherent_substitutions']),1)
  self.assertEqual(r['causal_diagnosis'],'unresolved');self.assertFalse(r['variant_calls_made'])
 def test_mixed_alleles_are_not_coherent_reference_error(self):
  text='\n'.join(line(str(i),'10M','AAAACAAAAA' if i<20 else 'A'*10,flag=16 if i%2 else 0) for i in range(40))
  r=analyze_sam(text,'A'*20,'x',1,10)
  self.assertIn('mixed_allele_signal',r['signals']);self.assertEqual(r['coherent_substitutions'],[])
