import contextlib
import io
import tempfile
import unittest
from pathlib import Path

from panel_bias_auditor.parsers import parse_bed
from panel_bias_auditor.performance import (
    parse_performance_table,
    performance_rows_from_callability,
    write_performance_table,
)
from research.benchmark_regions_to_performance import main as benchmark_main


class PublicBenchmarkConversionTests(unittest.TestCase):
    def test_callability_beds_become_performance_rows(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            base = Path(tmpdir)
            evaluated = base / "evaluated.bed"
            callable_bed = base / "callable.bed"
            performance = base / "performance.tsv"
            evaluated.write_text(
                "chr1\t0\t100\tmostly_callable\n"
                "chr1\t100\t200\tpartly_callable\n"
                "chr1\t200\t300\tmostly_uncallable\n",
                encoding="utf-8",
            )
            callable_bed.write_text(
                "chr1\t0\t95\ttruth_1\n"
                "chr1\t120\t200\ttruth_2\n"
                "chr1\t250\t260\ttruth_3\n",
                encoding="utf-8",
            )
            rows = performance_rows_from_callability(
                parse_bed(evaluated, source="evaluated"),
                parse_bed(callable_bed, source="callable"),
                min_callable_fraction=0.9,
            )
            write_performance_table(rows, performance)
            parsed = parse_performance_table(performance)

        self.assertEqual([row.status for row in parsed], ["pass", "low_callability", "low_callability"])
        self.assertAlmostEqual(parsed[0].callable_fraction or 0.0, 0.95)
        self.assertAlmostEqual(parsed[1].callable_fraction or 0.0, 0.80)
        self.assertAlmostEqual(parsed[2].callable_fraction or 0.0, 0.10)

    def test_research_script_writes_conversion_and_optional_enrichment(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            base = Path(tmpdir)
            evaluated = base / "evaluated.bed"
            callable_bed = base / "callable.bed"
            risk_bed = base / "risk.bed"
            out = base / "out"
            evaluated.write_text(
                "chr1\t0\t100\tpass_region\n"
                "chr1\t100\t200\tfail_region\n",
                encoding="utf-8",
            )
            callable_bed.write_text("chr1\t0\t100\tcallable_pass\nchr1\t100\t150\tcallable_partial\n", encoding="utf-8")
            risk_bed.write_text("chr1\t100\t200\trisk\ttype=gc_extreme,gc_class=high_gc\n", encoding="utf-8")

            with contextlib.redirect_stdout(io.StringIO()):
                exit_code = benchmark_main(
                    [
                        "--evaluated-bed",
                        str(evaluated),
                        "--callable-bed",
                        str(callable_bed),
                        "--risk-bed",
                        str(risk_bed),
                        "--output-dir",
                        str(out),
                        "--prefix",
                        "demo",
                        "--min-callable-fraction",
                        "0.9",
                    ]
                )

            self.assertEqual(exit_code, 0)
            self.assertTrue((out / "demo_public_benchmark_performance.tsv").exists())
            self.assertTrue((out / "demo_callability_conversion.md").exists())
            self.assertTrue((out / "demo_public_benchmark_enrichment.json").exists())


if __name__ == "__main__":
    unittest.main()
