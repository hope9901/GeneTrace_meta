import json,pathlib,tempfile,unittest
from unittest.mock import patch
from genetrace import demo
class SubmissionDemoGuards(unittest.TestCase):
 def test_unknown_case_is_rejected(self):
  with self.assertRaises(ValueError):demo.load_case('../../.env')
 def test_bundled_cases_have_hashed_real_sources(self):
  for key in ['nimt','provenance','contamination']:
   case=demo.load_case(key);self.assertFalse(case['model_called']);self.assertGreater(len(case['tools']),2)
   for t in case['tools']:
    for s in t['sources']:self.assertEqual(len(s['sha256']),64)
 def report(self,mode):return dict(id='a'*32,case_id='nimt',mode=mode,status='queued',model_calls=0,evidence=[],hypotheses=[],trace=[])
 def test_baseline_never_calls_model(self):
  with tempfile.TemporaryDirectory() as folder,patch.object(demo,'RUNS',pathlib.Path(folder)),patch.object(demo,'request_model') as model:
   r=self.report('baseline');demo.LOCK.acquire();demo.investigate(r);model.assert_not_called();self.assertEqual(r['status'],'complete');self.assertEqual(r['model_calls'],0);self.assertEqual(len(r['evidence']),4)
 def test_failed_meta_is_not_a_successful_baseline(self):
  with tempfile.TemporaryDirectory() as folder,patch.object(demo,'RUNS',pathlib.Path(folder)),patch.object(demo,'request_model',side_effect=RuntimeError('API failure')):
   r=self.report('meta');demo.LOCK.acquire();demo.investigate(r);self.assertEqual(r['status'],'failed');self.assertEqual(r['mode'],'meta');self.assertEqual(r['model_calls'],0)
 def test_mocked_model_invocations_are_recorded_and_bounded(self):
  def planner(context,available,evidence,config):return dict(next_tool=available[-1] if available else None,hypotheses=[]),{}
  with tempfile.TemporaryDirectory() as folder,patch.object(demo,'RUNS',pathlib.Path(folder)),patch.object(demo,'request_model',side_effect=planner),patch.object(demo,'model_config',return_value={'model':'test-only','key':'fake','base':'https://api.meta.ai/v1'}):
   r=self.report('meta');demo.LOCK.acquire();demo.investigate(r);self.assertEqual(r['status'],'complete');self.assertEqual(r['model_calls'],4);self.assertEqual(len(r['evidence']),4)
