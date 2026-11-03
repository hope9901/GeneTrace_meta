import unittest
from genetrace.read_provenance import assess_read_provenance
class ReadProvenanceGuards(unittest.TestCase):
 def sam(self,regular=True):
  return '\n'.join(f'r{i}\t99\tchr\t{100+i*10 if regular else 100+i*i}\t60\t125M\t=\t{500+i*10}\t450\t'+('A'*125)+'\t'+('B'*125)+'\tNM:i:0' for i in range(150))
 def test_regular_tiling_requires_review_but_does_not_prove_simulation(self):
  r=assess_read_provenance(self.sam(),{'B':100000});self.assertTrue(r['provenance_review_required']);self.assertFalse(r['synthetic_data_confirmed']);self.assertFalse(r['independent_physical_read_validation_established'])
 def test_constant_quality_alone_does_not_trigger_tiling_review(self):
  self.assertFalse(assess_read_provenance(self.sam(False),{'B':100000})['provenance_review_required'])
 def test_small_sample_cannot_trigger_review(self):
  self.assertFalse(assess_read_provenance(self.sam().splitlines()[0],{'B':125})['provenance_review_required'])
