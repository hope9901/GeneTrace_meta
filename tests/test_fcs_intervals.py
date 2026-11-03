import importlib.util,pathlib,unittest
p=pathlib.Path(__file__).resolve().parents[1]/'scripts/benchmark_fcs_ecoli_controls.py';spec=importlib.util.spec_from_file_location('fcs_benchmark',p);module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
class PopulationAlignmentCoverageGuards(unittest.TestCase):
 def test_overlapping_reference_matches_do_not_double_count_query(self):self.assertEqual(module.coverage([(0,100),(50,150),(20,30)]),150)
 def test_disjoint_intervals_retain_gaps(self):self.assertEqual(module.coverage([(0,10),(20,30)]),20)
 def test_no_alignment_is_zero_coverage(self):self.assertEqual(module.coverage([]),0)
