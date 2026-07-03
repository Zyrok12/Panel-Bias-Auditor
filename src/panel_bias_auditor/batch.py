from __future__ import annotations

import csv
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

from .core import audit_panel, pct
from .errors import InputFormatError
from .models import Interval
from .parsers import parse_bed


@dataclass(frozen=True)
class PanelDefinition:
    name: str
    path: Path
    assay_type: str = ""
    notes: str = ""


def _detect_delimiter(path: Path) -> str:
    if not path.exists():
        raise InputFormatError(f"Panel manifest does not exist: {path}")
    if not path.is_file():
        raise InputFormatError(f"Panel manifest path is not a file: {path}")
    lines = path.read_text(encoding="utf-8").splitlines()
    if not lines:
        raise InputFormatError(f"Panel manifest is empty: {path}")
    first_line = lines[0]
    return "\t" if "\t" in first_line else ","


def load_panel_manifest(path: str | Path) -> list[PanelDefinition]:
    manifest_path = Path(path)
    delimiter = _detect_delimiter(manifest_path)
    panels: list[PanelDefinition] = []
    with manifest_path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter=delimiter)
        required = {"name", "path"}
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise InputFormatError(f"Panel manifest is missing required columns: {', '.join(sorted(missing))}")
        for row_number, row in enumerate(reader, start=2):
            name = (row.get("name") or "").strip()
            path_raw = (row.get("path") or "").strip()
            if not name or not path_raw:
                raise InputFormatError(f"{manifest_path}:{row_number}: panel rows require name and path")
            panel_path = Path(path_raw)
            if not panel_path.is_absolute():
                panel_path = manifest_path.parent / panel_path
            if not panel_path.exists():
                raise InputFormatError(f"{manifest_path}:{row_number}: panel BED does not exist: {panel_path}")
            panels.append(
                PanelDefinition(
                    name=name,
                    path=panel_path,
                    assay_type=(row.get("assay_type") or "").strip(),
                    notes=(row.get("notes") or "").strip(),
                )
            )
    if not panels:
        raise InputFormatError("Panel manifest did not contain any panels")
    return panels


def compare_panels(
    panels: Iterable[PanelDefinition],
    critical_regions: list[Interval],
    difficult_regions: list[Interval],
    genome_build: str,
    track_metadata: dict[str, object] | None = None,
) -> dict[str, object]:
    results: list[dict[str, object]] = []
    for panel_def in panels:
        panel_intervals = parse_bed(panel_def.path, source=f"panel:{panel_def.name}")
        report = audit_panel(
            panel_intervals,
            critical_regions=critical_regions,
            difficult_regions=difficult_regions,
            variants=[],
            genome_build=genome_build,
            track_metadata=track_metadata,
            panel_name=panel_def.name,
        )
        critical = report["critical_regions"]
        difficult = report["difficult_regions"]
        panel = report["panel"]
        risk = report["risk"]
        results.append(
            {
                "name": panel_def.name,
                "path": str(panel_def.path),
                "assay_type": panel_def.assay_type,
                "notes": panel_def.notes,
                "risk_score": risk["score"],
                "risk_label": risk["label"],
                "merged_target_bases": panel["merged_target_bases"],
                "interval_count": panel["interval_count"],
                "critical_coverage_fraction": critical["coverage_fraction"],
                "missing_critical_bases": critical["missing_bases"],
                "difficult_overlap_bases": difficult["unique_overlap_bases"],
                "difficult_region_fraction": risk["components"]["difficult_region_fraction"],
                "full_report": report,
            }
        )
    results.sort(key=lambda item: (float(item["risk_score"]), -float(item["critical_coverage_fraction"] or 0), str(item["name"])))
    for rank, item in enumerate(results, start=1):
        item["rank"] = rank
    return {
        "genome_build": genome_build,
        "track_bundle": track_metadata,
        "panel_count": len(results),
        "results": results,
    }


def compact_batch_report(report: dict[str, object]) -> dict[str, object]:
    compact = dict(report)
    compact["results"] = [
        {key: value for key, value in item.items() if key != "full_report"}
        for item in report["results"]
    ]
    return compact


def render_compare_markdown(report: dict[str, object]) -> str:
    lines: list[str] = []
    lines.append("# Panel Comparison Report")
    lines.append("")
    lines.append(f"- Genome build: `{report['genome_build']}`")
    if report.get("track_bundle"):
        bundle = report["track_bundle"]
        lines.append(f"- Track bundle: `{bundle.get('name', 'unknown')}` version `{bundle.get('version', '')}`")
    lines.append(f"- Panels compared: **{report['panel_count']}**")
    lines.append("")
    lines.append("| Rank | Panel | Risk | Critical coverage | Missing critical bp | Difficult overlap bp | Target bp |")
    lines.append("|---:|---|---:|---:|---:|---:|---:|")
    for item in report["results"]:
        lines.append(
            "| {rank} | {name} | {score} ({label}) | {coverage} | {missing:,} | {difficult:,} | {target:,} |".format(
                rank=item["rank"],
                name=item["name"],
                score=item["risk_score"],
                label=item["risk_label"],
                coverage=pct(item["critical_coverage_fraction"]),
                missing=int(item["missing_critical_bases"]),
                difficult=int(item["difficult_overlap_bases"]),
                target=int(item["merged_target_bases"]),
            )
        )
    lines.append("")
    lines.append("## Interpretation")
    lines.append("")
    lines.append("Lower risk scores indicate better coverage of supplied critical regions with less difficult-region overlap.")
    lines.append("This comparison is only as good as the supplied track bundle; production use requires curated, versioned tracks.")
    lines.append("")
    return "\n".join(lines)


def write_compare_markdown(report: dict[str, object], path: str | Path) -> None:
    Path(path).write_text(render_compare_markdown(report), encoding="utf-8")


def write_compare_json(report: dict[str, object], path: str | Path, include_full_reports: bool = False) -> None:
    payload = report if include_full_reports else compact_batch_report(report)
    Path(path).write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")


def write_compare_tsv(report: dict[str, object], path: str | Path) -> None:
    fields = [
        "rank",
        "name",
        "assay_type",
        "risk_score",
        "risk_label",
        "critical_coverage_fraction",
        "missing_critical_bases",
        "difficult_overlap_bases",
        "difficult_region_fraction",
        "merged_target_bases",
        "interval_count",
        "path",
        "notes",
    ]
    with Path(path).open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields)
        writer.writeheader()
        for item in report["results"]:
            writer.writerow({field: item.get(field, "") for field in fields})
