import unittest
from genetrace.cohort_comparison import parse_cs,repeat_allele,repeat_bounds,summarize_marker,platform_group,distinct_loci,variable_reference_positions

class CohortComparisonGuards(unittest.TestCase):
 def test_uncallable_assemblies_do_not_count_as_reference_or_loss(self):
  r=summarize_marker([{'allele':'AAAA'},{'allele':None},{'allele':'AAAAA'}],'AAAAA')
  self.assertEqual(r['callable_assemblies'],2);self.assertEqual(r['uncallable_assemblies'],1)
  self.assertEqual(r['focal_allele_fraction'],.5);self.assertFalse(r['error_inferred'])
 def test_private_allele_is_not_automatically_an_error(self):
  r=summarize_marker([{'allele':'AAAA'}]*20,'AAAAA')
  self.assertEqual(r['focal_allele_count'],0);self.assertFalse(r['error_inferred'])
 def test_shifted_homopolymer_deletions_have_same_repeat_allele(self):
  ref='CAAAAAG';query='CAAAAG'
  a,ai=parse_cs('=CA-A=AAAG',0,0,query)
  b,bi=parse_cs('=CAAA-A=AG',0,0,query)
  self.assertEqual(repeat_allele(ref,3,a,ai)['allele'],repeat_allele(ref,3,b,bi)['allele'])
 def test_missing_outer_anchors_are_not_callable(self):
  self.assertIsNone(repeat_allele('CAAAAAG',3,{3:('A',2)},{}))
 def test_curated_publication_technology_preserves_multiplatform_context(self):
  self.assertEqual(platform_group({'publication_supported_platforms':'Illumina;PacBio','sequencing_group':'PacBio'}),'Illumina+PacBio')
 def test_short_cs_substitution_coordinates(self):
  m,ins=parse_cs('=AC*gt=TT+aa=CC',2,0,'ACTTTAACC')
  self.assertEqual(m[4][0],'T');self.assertEqual(ins[6],'AA')
 def test_unknown_cs_tokens_are_rejected(self):
  with self.assertRaises(ValueError):parse_cs('garbage',0,0,'AAAA')

 def test_variable_scan_collapses_repeat_deletion_and_respects_span(self):
  mapping,ins=parse_cs('=CA-A=AAAG',0,0,'CAAAAG')
  self.assertEqual(variable_reference_positions('CAAAAAG',[(mapping,ins)],1,6),[1])
  self.assertEqual(variable_reference_positions('CAAAAAG',[(mapping,ins)],6,7),[])
 def test_missing_comparison_does_not_create_a_variant(self):
  self.assertEqual(variable_reference_positions('CAAAAAG',[({}, {})],1,6),[])
