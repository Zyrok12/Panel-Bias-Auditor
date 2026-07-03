import tempfile
import unittest
from pathlib import Path

from panel_bias_auditor.cli import main
from panel_bias_auditor.models import Interval
from panel_bias_auditor.technical import (
    merge_gc_windows,
    merge_sequence_risk_windows,
    parse_fasta,
    scan_gc_windows,
    scan_sequence_risk_windows,
)


class TechnicalTrackTests(unittest.TestCase):
    def test_gc_scan_detects_high_and_low_gc_blocks(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            fasta = Path(tmpdir) / "toy.fa"
            fasta.write_text(">chr1\n" + "G" * 100 + "A" * 100 + "C" * 100 + "T" * 100 + "\n", encoding="utf-8")
            sequences = parse_fasta(fasta)
            windows = scan_gc_windows(
                sequences,
                regions=[Interval("1", 0, 400, "toy")],
                window_size=50,
                step_size=50,
                low_gc=0.2,
                high_gc=0.8,
            )
            intervals = merge_gc_windows(windows)
        classes = [interval.attrs["gc_class"] for interval in intervals]
        loci = {(interval.start, interval.end, interval.attrs["gc_class"]) for interval in intervals}
        self.assertEqual(classes.count("high_gc"), 2)
        self.assertEqual(classes.count("low_gc"), 2)
        self.assertIn((0, 100, "high_gc"), loci)
        self.assertIn((100, 200, "low_gc"), loci)

    def test_gcderive_cli_writes_bed_and_reports(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            base = Path(tmpdir)
            fasta = base / "toy.fa"
            regions = base / "regions.bed"
            out_bed = base / "gc.bed"
            report = base / "gc.md"
            json_report = base / "gc.json"
            fasta.write_text(">chr1\n" + "G" * 100 + "A" * 100 + "\n", encoding="utf-8")
            regions.write_text("chr1\t0\t200\ttoy\n", encoding="utf-8")
            exit_code = main(
                [
                    "gcderive",
                    "--fasta",
                    str(fasta),
                    "--regions",
                    str(regions),
                    "--window-size",
                    "50",
                    "--step-size",
                    "50",
                    "--low-gc",
                    "0.2",
                    "--high-gc",
                    "0.8",
                    "--out-bed",
                    str(out_bed),
                    "--report",
                    str(report),
                    "--json",
                    str(json_report),
                ]
            )
            bed_text = out_bed.read_text(encoding="utf-8")
            self.assertEqual(exit_code, 0)
            self.assertIn("high_gc", bed_text)
            self.assertIn("low_gc", bed_text)
            self.assertTrue(report.exists())
            self.assertTrue(json_report.exists())

    def test_sequence_risk_scan_detects_homopolymer_and_low_complexity(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            fasta = Path(tmpdir) / "toy.fa"
            fasta.write_text(">chr1\n" + "A" * 80 + "ACGT" * 20 + "C" * 70 + "GATTACA" * 10 + "\n", encoding="utf-8")
            sequences = parse_fasta(fasta)
            windows = scan_sequence_risk_windows(
                sequences,
                regions=[Interval("chr1", 0, 300, "toy")],
                window_size=40,
                step_size=40,
                min_homopolymer=10,
                low_complexity_fraction=0.8,
            )
            intervals = merge_sequence_risk_windows(windows)
        classes = {interval.attrs["sequence_risk_class"] for interval in intervals}
        self.assertIn("homopolymer", classes)
        self.assertTrue(any(int(interval.attrs["max_homopolymer"]) >= 10 for interval in intervals))

    def test_seqderive_cli_writes_bed_and_reports(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            base = Path(tmpdir)
            fasta = base / "toy.fa"
            regions = base / "regions.bed"
            out_bed = base / "seq.bed"
            report = base / "seq.md"
            json_report = base / "seq.json"
            fasta.write_text(">chr1\n" + "A" * 80 + "CGTA" * 20 + "\n", encoding="utf-8")
            regions.write_text("chr1\t0\t160\ttoy\n", encoding="utf-8")
            exit_code = main(
                [
                    "seqderive",
                    "--fasta",
                    str(fasta),
                    "--regions",
                    str(regions),
                    "--window-size",
                    "40",
                    "--step-size",
                    "40",
                    "--min-homopolymer",
                    "10",
                    "--out-bed",
                    str(out_bed),
                    "--report",
                    str(report),
                    "--json",
                    str(json_report),
                ]
            )
            bed_text = out_bed.read_text(encoding="utf-8")
            self.assertEqual(exit_code, 0)
            self.assertIn("homopolymer", bed_text)
            self.assertTrue(report.exists())
            self.assertTrue(json_report.exists())


if __name__ == "__main__":
    unittest.main()
