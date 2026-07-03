import unittest

from panel_bias_auditor.core import audit_panel, coverage_by_region, merged_interval_bases
from panel_bias_auditor.models import Interval, Variant


class CoreTests(unittest.TestCase):
    def test_merged_interval_bases_merges_overlaps_by_chrom(self):
        intervals = [
            Interval("chr1", 10, 20),
            Interval("1", 15, 30),
            Interval("chr2", 10, 20),
        ]
        self.assertEqual(merged_interval_bases(intervals), 30)

    def test_coverage_by_region_counts_unique_overlap(self):
        panel = [
            Interval("chr1", 10, 20),
            Interval("chr1", 15, 25),
        ]
        critical = [Interval("1", 0, 30, "region")]
        coverage = coverage_by_region(panel, critical)
        self.assertEqual(coverage[0].covered_bases, 15)
        self.assertAlmostEqual(coverage[0].coverage_fraction, 0.5)

    def test_audit_panel_flags_outside_variants(self):
        panel = [Interval("chr7", 100, 200, "target")]
        critical = [Interval("chr7", 100, 250, "critical")]
        difficult = [Interval("chr7", 150, 180, "repeat", attrs={"type": "repeat"})]
        variants = [
            Variant("7", 160, "A", "T", "inside"),
            Variant("7", 300, "A", "T", "outside"),
        ]
        report = audit_panel(panel, critical, difficult, variants, genome_build="test")
        self.assertEqual(report["variants"]["outside_panel_count"], 1)
        self.assertEqual(report["critical_regions"]["covered_bases"], 100)
        self.assertEqual(report["difficult_regions"]["unique_overlap_bases"], 30)
        self.assertEqual(report["metadata"]["genome_build"], "test")
        self.assertTrue(report["recommendations"])

    def test_audit_panel_recommends_track_manifest_when_raw_tracks_used(self):
        report = audit_panel(
            panel=[Interval("chr1", 10, 20, "panel")],
            critical_regions=[Interval("chr1", 10, 20, "critical")],
            difficult_regions=[],
            variants=[],
            genome_build="test",
        )
        categories = {item["category"] for item in report["recommendations"]}
        self.assertIn("provenance", categories)


if __name__ == "__main__":
    unittest.main()
