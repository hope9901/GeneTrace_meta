import io,json,unittest
from unittest.mock import patch
from genetrace.agent import request_model
CONFIG={'key':'test-only-not-real','model':'muse-spark-1.3','base':'https://api.meta.ai/v1'}
class MetaResponseGuards(unittest.TestCase):
 def check(self,body):
  with patch('genetrace.agent.urllib.request.urlopen',return_value=io.StringIO(json.dumps(body))):return request_model({},[],[],CONFIG)
 def test_reasoning_budget_exhaustion_is_explicit_failure(self):
  with self.assertRaisesRegex(ValueError,'budget exhausted'):self.check({'choices':[{'finish_reason':'length','message':{'content':None}}]})
 def test_empty_content_is_not_a_valid_model_conclusion(self):
  with self.assertRaisesRegex(ValueError,'no text content'):self.check({'choices':[{'finish_reason':'stop','message':{'content':None}}]})
 def test_valid_json_records_provider_model_and_usage(self):
  step,usage=self.check({'model':'muse-spark-1.3','choices':[{'finish_reason':'stop','message':{'content':'{"next_tool":null,"hypotheses":[]}'}}],'usage':{'total_tokens':20}})
  self.assertIsNone(step['next_tool']);self.assertEqual(usage['response_model'],'muse-spark-1.3');self.assertEqual(usage['total_tokens'],20)
