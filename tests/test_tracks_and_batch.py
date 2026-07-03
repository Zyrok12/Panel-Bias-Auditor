import json
import tempfile
import unittest
from pathlib import Path

from panel_bias_auditor.batch import compare_panels, load_panel_manifest
from panel_bias_auditor.errors import InputFormatError
from panel_bias_auditor.models import Interval
from panel_bias_auditor.tracks import load_track_manifest, load_tracks_from_bundle


class TrackAndBatchTests(unittest.TestCase):
    def test_track_manifest_loads_relative_paths(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            base = Path(tmpdir)
            (base / "critical.bed").write_text("chr1\t10\t20\tcritical_a\n", encoding="utf-8")
            manifest = base / "tracks.json"
            manifest.write_text(
                json.dumps(
                    {
                        "name": "demo",
                        "version": "1",
                        "genome_build": "GRCh38-test",
                        "tracks": [
                            {
                                "name": "critical",
                                "role": "critical",
                                "path": "critical.bed",
                            }
                        ],
                    }
                ),
                encoding="utf-8",
            )
            bundle = load_track_manifest(manifest)
            critical, difficult = load_tracks_from_bundle(bundle)
        self.assertEqual(bundle.genome_build, "GRCh38-test")
        self.assertEqual(len(critical), 1)
        self.assertEqual(len(difficult), 0)
        self.assertEqual(critical[0].name, "critical_a")

    def test_compare_panels_ranks_lower_risk_first(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            base = Path(tmpdir)
            (base / "weak.bed").write_text("chr1\t10\t15\tweak\n", encoding="utf-8")
            (base / "strong.bed").write_text("chr1\t10\t20\tstrong\n", encoding="utf-8")
            manifest = base / "panels.tsv"
            manifest.write_text(
                "name\tpath\tassay_type\tnotes\n"
                "Weak\tweak.bed\tdemo\tpartial\n"
                "Strong\tstrong.bed\tdemo\tcomplete\n",
                encoding="utf-8",
            )
            panels = load_panel_manifest(manifest)
            report = compare_panels(
                panels,
                critical_regions=[Interval("chr1", 10, 20, "critical")],
                difficult_regions=[],
                genome_build="test",
            )
        self.assertEqual(report["results"][0]["name"], "Strong")
        self.assertLess(report["results"][0]["risk_score"], report["results"][1]["risk_score"])

    def test_panel_manifest_reports_missing_panel_bed(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            base = Path(tmpdir)
            manifest = base / "panels.tsv"
            manifest.write_text("name\tpath\nMissing\tmissing.bed\n", encoding="utf-8")
            with self.assertRaisesRegex(InputFormatError, "panel BED does not exist"):
                load_panel_manifest(manifest)


if __name__ == "__main__":
    unittest.main()
