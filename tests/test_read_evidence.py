import unittest
from genetrace.read_evidence import quality_encoding_assessment,regional_evidence_assessment

class ReadEvidenceGuards(unittest.TestCase):
    def test_all_zero_encoded_qualities_are_not_a_sequencing_failure(self):
        r=quality_encoding_assessment({'!':83047354})
        self.assertTrue(r['bq_filter_can_remove_all_observed_bases'])
        self.assertFalse(r['sequencing_failure_inferred'])
        self.assertEqual(r['meaning'],'quality_information_unresolved')
    def test_prefix_zero_coverage_cannot_diagnose_absence_or_low_coverage(self):
        r=regional_evidence_assessment(0,0,True)
        self.assertEqual(r['classification'],'unresolved_read_support')
        self.assertFalse(r['low_coverage_diagnosis_allowed'])
        self.assertFalse(r['biological_loss_inferred'])
    def test_positive_prefix_support_is_candidate_presence_only(self):
        r=regional_evidence_assessment(1,17,True)
        self.assertEqual(r['classification'],'read_supported_sequence_presence_candidate')
        self.assertEqual(r['functional_status'],'not_assessed')
        self.assertEqual(r['independent_molecule_support'],'not_established')
    def test_unverified_full_input_cannot_enable_low_coverage_diagnosis(self):
        self.assertFalse(regional_evidence_assessment(.5,0,False)['low_coverage_diagnosis_allowed'])
    def test_empty_quality_data_are_not_normal_quality(self):
        self.assertEqual(quality_encoding_assessment({})['meaning'],'no_quality_observations')
