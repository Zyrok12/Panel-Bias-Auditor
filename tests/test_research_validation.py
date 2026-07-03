import json
import tempfile
import unittest
from pathlib import Path

from research.validate_empirical_tracks import validate_outputs


class ResearchValidationTests(unittest.TestCase):
    def write_outputs(self, base: Path, combined_count: int = 2) -> None:
        (base / "demo_gc_extremes.bed").write_text(
            "# gc\n"
            "chr1\t0\t50\tgc1\ttype=gc_extreme,gc_class=high_gc,curation_status=source_verified,evidence_level=assay_computable_gc,source=fasta_gcderive,min_gc=0.900,max_gc=0.950\n",
            encoding="utf-8",
        )
        (base / "demo_sequence_risk.bed").write_text(
            "# seq\n"
            "chr1\t50\t100\tseq1\ttype=homopolymer,sequence_risk_class=homopolymer,curation_status=source_verified,evidence_level=assay_computable_sequence_complexity,source=fasta_seqderive,max_homopolymer=12,max_dominant_base_fraction=0.900\n",
            encoding="utf-8",
        )
        (base / "demo_combined_technical_risk.bed").write_text(
            "# combined\n"
            "chr1\t0\t50\tgc1\ttype=gc_extreme,gc_class=high_gc,curation_status=source_verified,evidence_level=assay_computable_gc,source=fasta_gcderive,min_gc=0.900,max_gc=0.950\n"
            "chr1\t50\t100\tseq1\ttype=homopolymer,sequence_risk_class=homopolymer,curation_status=source_verified,evidence_level=assay_computable_sequence_complexity,source=fasta_seqderive,max_homopolymer=12,max_dominant_base_fraction=0.900\n",
            encoding="utf-8",
        )
        (base / "demo_gc_report.json").write_text(json.dumps({"merged_interval_count": 1}), encoding="utf-8")
        (base / "demo_sequence_risk_report.json").write_text(
            json.dumps({"merged_interval_count": 1}), encoding="utf-8"
        )
        (base / "demo_research_summary.json").write_text(
            json.dumps(
                {
                    "gc_interval_count": 1,
                    "sequence_risk_interval_count": 1,
                    "combined_interval_count": combined_count,
                    "combined_merged_bases": 100,
                }
            ),
            encoding="utf-8",
        )

    def test_validate_outputs_passes_on_consistent_artifacts(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            base = Path(tmpdir)
            self.write_outputs(base)
            report = validate_outputs(base, "demo")
        self.assertEqual(report["status"], "pass")
        self.assertEqual(report["fail_count"], 0)

    def test_validate_outputs_fails_on_summary_mismatch(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            base = Path(tmpdir)
            self.write_outputs(base, combined_count=99)
            report = validate_outputs(base, "demo")
        self.assertEqual(report["status"], "fail")
        self.assertTrue(any(item["check"] == "summary_combined_count_mismatch" for item in report["findings"]))


if __name__ == "__main__":
    unittest.main()
