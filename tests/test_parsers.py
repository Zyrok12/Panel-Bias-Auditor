import tempfile
import unittest
from pathlib import Path

from panel_bias_auditor.errors import InputFormatError
from panel_bias_auditor.parsers import parse_bed, parse_vcf


class ParserTests(unittest.TestCase):
    def test_parse_bed_reads_attrs(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "regions.bed"
            path.write_text("chr1\t10\t20\tregion_a\ttype=repeat,source=demo\n", encoding="utf-8")
            regions = parse_bed(path, source="test")
        self.assertEqual(len(regions), 1)
        self.assertEqual(regions[0].name, "region_a")
        self.assertEqual(regions[0].attrs["type"], "repeat")
        self.assertEqual(regions[0].attrs["source"], "demo")

    def test_parse_vcf_uses_end_for_sv(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "variants.vcf"
            path.write_text(
                "##fileformat=VCFv4.2\n"
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
                "chr1\t101\tsv1\tN\t<DEL>\t.\tPASS\tEND=250;SVTYPE=DEL\n",
                encoding="utf-8",
            )
            variants = parse_vcf(path)
        self.assertEqual(len(variants), 1)
        self.assertEqual(variants[0].start, 100)
        self.assertEqual(variants[0].end, 250)

    def test_parse_bed_rejects_end_before_start(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "bad.bed"
            path.write_text("chr1\t20\t10\tbad\n", encoding="utf-8")
            with self.assertRaisesRegex(InputFormatError, "end must be greater"):
                parse_bed(path, source="test")

    def test_parse_vcf_rejects_missing_alt(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "bad.vcf"
            path.write_text(
                "##fileformat=VCFv4.2\n"
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
                "chr1\t101\tv1\tA\t.\t.\tPASS\t.\n",
                encoding="utf-8",
            )
            with self.assertRaisesRegex(InputFormatError, "ALT must be present"):
                parse_vcf(path)


if __name__ == "__main__":
    unittest.main()
