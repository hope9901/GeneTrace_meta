import datetime,unittest
from unittest.mock import patch
from genetrace.background import job_snapshot,ROOT
class BackgroundStatusGuards(unittest.TestCase):
 def setUp(self):
  self.now=datetime.datetime.now(datetime.timezone.utc);self.config=('test','test','missing/status.json','missing/log','missing/core');self.status=dict(status='running',pid=42,updated_utc=self.now.isoformat(),stage='core_inference',core_size=2)
 def get(self,processes,error=None):
  with patch('genetrace.background.read_json',return_value=(self.status,None)):return job_snapshot(self.config,processes,error,self.now)
 def test_dead_process_does_not_appear_running_from_status_file(self):self.assertEqual(self.get({})['observed_status'],'not_running')
 def test_reused_pid_is_not_the_original_job(self):self.assertFalse(self.get({42:dict(ppid=1,stat='S',args='unrelated job',cpu=1)})['process_alive'])
 def test_process_check_failure_is_unknown_not_stopped(self):self.assertEqual(self.get({},'unavailable')['observed_status'],'unknown')
 def test_matching_live_process_has_no_invented_total_progress(self):
  r=self.get({42:dict(ppid=1,stat='S',args=str(ROOT/'tools/orthofinder')+' -b '+str(ROOT/'missing/core'),cpu=50)})
  self.assertTrue(r['process_alive']);self.assertEqual(r['cpu_percent'],50);self.assertIsNone(r['overall_percent'])
 def test_old_heartbeat_is_flagged(self):
  self.status['updated_utc']=(self.now-datetime.timedelta(seconds=120)).isoformat();self.assertTrue(self.get({})['heartbeat_stale'])
