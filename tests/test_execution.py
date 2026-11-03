import copy,json,pathlib,tempfile,unittest
from unittest.mock import patch
from genetrace import execution as ex
class ExecutionGuards(unittest.TestCase):
 def candidate(self,signals=None):
  return dict(id='E003',tool='align_candidate_reads',result=dict(candidate=dict(signals=signals or ['coherent_indel_disagreement'])))
 def test_neighbor_check_is_exposed_after_disagreement_only(self):
  done={'inspect_annotation','align_candidate_reads'}
  self.assertIn('check_neighbor_regions',ex.available_tools(done,[self.candidate()]))
  self.assertNotIn('check_neighbor_regions',ex.available_tools(done,[self.candidate(['read_assembly_supported'])]))
 def test_cannot_skip_required_initial_checks(self):
  self.assertEqual(ex.available_tools(set(),[]),['inspect_annotation'])
  self.assertEqual(ex.available_tools({'inspect_annotation'},[]),['align_candidate_reads'])
 def test_candidate_alone_is_provisional_not_local_error_diagnosis(self):
  j=ex.verdict([self.candidate()]);self.assertEqual(j['status'],'provisional');self.assertFalse(j['biological_loss_inferred']);self.assertEqual(j['causal_diagnosis'],'unresolved')
 def test_neighbor_disagreement_changes_observation_but_not_causal_diagnosis(self):
  ev=[self.candidate(),dict(id='E004',tool='check_neighbor_regions',result=dict(neighbor_regions=[dict(signals=['coherent_substitution_disagreement'])])),dict(id='E005',tool='inspect_read_provenance',result={})]
  j=ex.verdict(ev);self.assertEqual(j['status'],'unresolved');self.assertFalse(j['local_only_explanation_supported']);self.assertIn('Reference-selected',j['limitations'][0])
 def test_input_tampering_is_rejected(self):
  with tempfile.TemporaryDirectory() as p:
   pack=pathlib.Path(p);(pack/'manifest.json').write_text(json.dumps(dict(files={'reads_1.fastq.gz':{'sha256':'0'*64}})));(pack/'reads_1.fastq.gz').write_bytes(b'tampered')
   with patch.object(ex,'PACK',pack),self.assertRaisesRegex(ValueError,'checksum'):ex.verify_pack()
 def test_unknown_input_path_rejected(self):
  with tempfile.TemporaryDirectory() as p:
   pack=pathlib.Path(p);(pack/'manifest.json').write_text(json.dumps(dict(files={'../../.env':{}})))
   with patch.object(ex,'PACK',pack),self.assertRaisesRegex(ValueError,'Unexpected'):ex.verify_pack()
 def test_unknown_tool_cannot_execute(self):
  obj=object.__new__(ex.Executor)
  with self.assertRaisesRegex(ValueError,'not allowed'):obj.execute('arbitrary_shell')
 def test_real_case_manifest_is_partial_and_hashed(self):
  m=ex.verify_pack();self.assertTrue(m['selection']['partial_input']);self.assertFalse(m['exact_assembly_run_use_verified']);self.assertEqual(m['source_fastq'][0]['selected_records'],756)
