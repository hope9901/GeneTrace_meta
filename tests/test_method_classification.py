import importlib.util
from pathlib import Path
import unittest

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location('methods', ROOT / 'scripts/analyze_methods_metadata.py')
methods = importlib.util.module_from_spec(spec)
spec.loader.exec_module(methods)


class ClassificationGuards(unittest.TestCase):
    def test_hifiasm_does_not_imply_hifi(self):
        result = methods.technology('PacBio Sequel', 'hifiasm v. 0.19')
        self.assertEqual(result['pacbio_hifi_reported'], 'N/A')
        self.assertEqual(result['sequencing_group'], 'PacBio (mode unspecified)')

    def test_hiseq_does_not_imply_hic(self):
        self.assertEqual(methods.technology('Illumina HiSeq 2000', 'Canu')['hic_reported'], 'N/A')

    def test_explicit_hifi_and_hic(self):
        result = methods.technology('PacBio HiFi CCS; Illumina Hi-C', 'hifiasm')
        self.assertEqual(result['pacbio_hifi_reported'], 'Yes')
        self.assertEqual(result['hic_reported'], 'Yes')
        self.assertEqual(result['sequencing_group'], 'Illumina + PacBio HiFi/CCS + Hi-C')

    def test_missing_is_not_zero(self):
        result = methods.technology('N/A', 'N/A')
        self.assertEqual(result['sequencing_group'], 'N/A')
        self.assertEqual(result['hic_reported'], 'N/A')

    def test_functional_tools_and_provider_not_structural_tools(self):
        row = dict(report_fields={}, genbank_first_record_header='COMMENT     Functional annotations by InterProScan, BLAST and eggNOG.\n',
                   gff_structural_feature_sources={'Genbank': 12}, gff_header=['#!annotation-source SGD'],
                   report_url='report', gff_url='gff', gbff_url='gbff')
        result = methods.extract(row)
        self.assertEqual(result['structural_annotation_tools'], 'N/A')

    def test_explicit_structural_comment_and_version(self):
        row = dict(report_fields={}, genbank_first_record_header='COMMENT     Annotated by AUGUSTUS ver.3.3.3\n',
                   gff_structural_feature_sources={'Genbank': 12}, gff_header=[],
                   report_url='report', gff_url='gff', gbff_url='gbff')
        result = methods.extract(row)
        self.assertEqual(result['structural_annotation_tools'], 'AUGUSTUS')
        self.assertEqual(result['structural_annotation_tool_versions'], 'AUGUSTUS 3.3.3')


if __name__ == '__main__':
    unittest.main()
