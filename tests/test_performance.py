import contextlib
import io
import tempfile
import unittest
from pathlib import Path

from panel_bias_auditor.models import Interval
from panel_bias_auditor.performance import analyze_performance_enrichment, parse_performance_table
from research.assay_performance_enrichment import main as performance_main


class PerformanceAnalysisTests(unittest.TestCase):
    def test_parse_performance_table_and_classify_failures(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            table = Path(tmpdir) / "performance.tsv"
            table.write_text(
                "chrom\tstart\tend\tname\tmean_depth\tcallable_fraction\tstatus\tfalse_negative_count\tfalse_positive_count\ttotal_variant_count\n"
                "chr1\t0\t100\tok\t500\t0.99\tpass\t0\t0\t10\n"
                "chr1\t100\t200\tlow\t40\t0.80\tpass\t0\t0\t10\n"
                "chr1\t200\t300\tfn\t500\t0.99\tpass\t1\t0\t10\n",
                encoding="utf-8",
            )
            rows = parse_performance_table(table)
        self.assertEqual(len(rows), 3)
        self.assertFalse(rows[0].is_failure(min_depth=100, min_callable_fraction=0.9))
        self.assertTrue(rows[1].is_failure(min_depth=100, min_callable_fraction=0.9))
        self.assertTrue(rows[2].is_failure(min_depth=100, min_callable_fraction=0.9))

    def test_enrichment_detects_failures_in_risk_regions(self):
        performance = [
            "chrom\tstart\tend\tname\tmean_depth\tcallable_fraction\tstatus\n",
            "chr1\t0\t100\tfail_risk_1\t20\t0.50\tpass\n",
            "chr1\t100\t200\tfail_risk_2\t30\t0.60\tpass\n",
            "chr1\t200\t300\tpass_nonrisk_1\t500\t0.99\tpass\n",
            "chr1\t300\t400\tpass_nonrisk_2\t500\t0.99\tpass\n",
        ]
        with tempfile.TemporaryDirectory() as tmpdir:
            table = Path(tmpdir) / "performance.tsv"
            table.write_text("".join(performance), encoding="utf-8")
            rows = parse_performance_table(table)
        report = analyze_performance_enrichment(
            rows,
            [Interval("chr1", 0, 200, "risk", attrs={"type": "gc_extreme"})],
            min_depth=100,
            min_callable_fraction=0.9,
        )
        self.assertEqual(report["summary"]["failure_interval_count"], 2)
        self.assertEqual(report["overall_enrichment"]["risk_fail_bases"], 200)
        self.assertEqual(report["overall_enrichment"]["nonrisk_fail_bases"], 0)
        self.assertGreater(report["overall_enrichment"]["odds_ratio_haldane"], 1)

    def test_research_script_writes_enrichment_reports(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            base = Path(tmpdir)
            performance = base / "performance.tsv"
            risk = base / "risk.bed"
            out = base / "out"
            performance.write_text(
                "chrom\tstart\tend\tname\tmean_depth\tcallable_fraction\tstatus\n"
                "chr1\t0\t100\tfail\t20\t0.5\tpass\n"
                "chr1\t100\t200\tok\t500\t0.99\tpass\n",
                encoding="utf-8",
            )
            risk.write_text(
                "chr1\t0\t100\trisk\ttype=gc_extreme,gc_class=high_gc,source=fasta_gcderive\n",
                encoding="utf-8",
            )
            with contextlib.redirect_stdout(io.StringIO()):
                exit_code = performance_main(
                    [
                        "--performance",
                        str(performance),
                        "--risk-bed",
                        str(risk),
                        "--output-dir",
                        str(out),
                        "--prefix",
                        "demo",
                        "--min-depth",
                        "100",
                        "--min-callable-fraction",
                        "0.9",
                    ]
                )
            self.assertEqual(exit_code, 0)
            self.assertTrue((out / "demo_assay_performance_enrichment.md").exists())
            self.assertTrue((out / "demo_assay_performance_enrichment.json").exists())


if __name__ == "__main__":
    unittest.main()
