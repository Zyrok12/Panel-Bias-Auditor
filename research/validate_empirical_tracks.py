from __future__ import annotations

import argparse
import json
from dataclasses import dataclass
from pathlib import Path

from panel_bias_auditor.core import merged_interval_bases
from panel_bias_auditor.parsers import parse_bed


EXPECTED_FILES = {
    "gc_bed": "{prefix}_gc_extremes.bed",
    "gc_json": "{prefix}_gc_report.json",
    "sequence_bed": "{prefix}_sequence_risk.bed",
    "sequence_json": "{prefix}_sequence_risk_report.json",
    "combined_bed": "{prefix}_combined_technical_risk.bed",
    "summary_json": "{prefix}_research_summary.json",
}

GC_REQUIRED_ATTRS = {"type", "gc_class", "curation_status", "evidence_level", "source", "min_gc", "max_gc"}
SEQ_REQUIRED_ATTRS = {
    "type",
    "sequence_risk_class",
    "curation_status",
    "evidence_level",
    "source",
    "max_homopolymer",
    "max_dominant_base_fraction",
}


@dataclass(frozen=True)
class ValidationFinding:
    severity: str
    check: str
    message: str


def load_json(path: Path) -> dict[str, object]:
    return json.loads(path.read_text(encoding="utf-8"))


def expected_paths(output_dir: str | Path, prefix: str) -> dict[str, Path]:
    base = Path(output_dir)
    return {label: base / template.format(prefix=prefix) for label, template in EXPECTED_FILES.items()}


def count_noncomment_bed_lines(path: Path) -> int:
    return sum(1 for line in path.read_text(encoding="utf-8").splitlines() if line and not line.startswith("#"))


def validate_attrs(intervals, required_attrs: set[str], label: str) -> list[ValidationFinding]:
    findings: list[ValidationFinding] = []
    for interval in intervals:
        missing = sorted(attr for attr in required_attrs if attr not in interval.attrs)
        if missing:
            findings.append(
                ValidationFinding(
                    "fail",
                    f"{label}_attrs",
                    f"{interval.name} is missing required attrs: {', '.join(missing)}",
                )
            )
    return findings


def validate_outputs(output_dir: str | Path, prefix: str = "demo") -> dict[str, object]:
    paths = expected_paths(output_dir, prefix)
    findings: list[ValidationFinding] = []
    for label, path in paths.items():
        if not path.exists():
            findings.append(ValidationFinding("fail", "missing_file", f"Missing {label}: {path}"))

    if findings:
        return build_report(paths, findings, {}, {}, [], [], [])

    gc_intervals = parse_bed(paths["gc_bed"], source="gc_validation")
    seq_intervals = parse_bed(paths["sequence_bed"], source="sequence_validation")
    combined_intervals = parse_bed(paths["combined_bed"], source="combined_validation")
    gc_json = load_json(paths["gc_json"])
    seq_json = load_json(paths["sequence_json"])
    summary_json = load_json(paths["summary_json"])

    findings.extend(validate_attrs(gc_intervals, GC_REQUIRED_ATTRS, "gc"))
    findings.extend(validate_attrs(seq_intervals, SEQ_REQUIRED_ATTRS, "sequence"))

    if int(gc_json.get("merged_interval_count", -1)) != len(gc_intervals):
        findings.append(ValidationFinding("fail", "gc_count_mismatch", "GC JSON merged_interval_count does not match GC BED."))
    if int(seq_json.get("merged_interval_count", -1)) != len(seq_intervals):
        findings.append(
            ValidationFinding(
                "fail",
                "sequence_count_mismatch",
                "Sequence JSON merged_interval_count does not match sequence BED.",
            )
        )
    if int(summary_json.get("gc_interval_count", -1)) != len(gc_intervals):
        findings.append(ValidationFinding("fail", "summary_gc_count_mismatch", "Summary GC count does not match GC BED."))
    if int(summary_json.get("sequence_risk_interval_count", -1)) != len(seq_intervals):
        findings.append(
            ValidationFinding(
                "fail",
                "summary_sequence_count_mismatch",
                "Summary sequence-risk count does not match sequence BED.",
            )
        )
    if count_noncomment_bed_lines(paths["combined_bed"]) != len(gc_intervals) + len(seq_intervals):
        findings.append(
            ValidationFinding(
                "fail",
                "combined_line_count_mismatch",
                "Combined BED line count does not equal GC plus sequence-risk BED line count.",
            )
        )
    if int(summary_json.get("combined_interval_count", -1)) != len(combined_intervals):
        findings.append(
            ValidationFinding(
                "fail",
                "summary_combined_count_mismatch",
                "Summary combined interval count does not match combined BED.",
            )
        )
    combined_bases = merged_interval_bases(combined_intervals)
    if int(summary_json.get("combined_merged_bases", -1)) != combined_bases:
        findings.append(
            ValidationFinding(
                "fail",
                "summary_combined_bases_mismatch",
                "Summary combined merged bases does not match parsed combined BED.",
            )
        )
    if not gc_intervals:
        findings.append(ValidationFinding("warn", "no_gc_intervals", "GC BED contains no empirical GC-risk intervals."))
    if not seq_intervals:
        findings.append(
            ValidationFinding("warn", "no_sequence_intervals", "Sequence-risk BED contains no empirical sequence-risk intervals.")
        )

    return build_report(paths, findings, gc_json, seq_json, gc_intervals, seq_intervals, combined_intervals)


def build_report(
    paths: dict[str, Path],
    findings: list[ValidationFinding],
    gc_json: dict[str, object],
    seq_json: dict[str, object],
    gc_intervals,
    seq_intervals,
    combined_intervals,
) -> dict[str, object]:
    fail_count = sum(1 for finding in findings if finding.severity == "fail")
    warn_count = sum(1 for finding in findings if finding.severity == "warn")
    status = "fail" if fail_count else "pass_with_warnings" if warn_count else "pass"
    return {
        "status": status,
        "fail_count": fail_count,
        "warning_count": warn_count,
        "paths": {label: str(path) for label, path in paths.items()},
        "metrics": {
            "gc_interval_count": len(gc_intervals),
            "sequence_risk_interval_count": len(seq_intervals),
            "combined_interval_count": len(combined_intervals),
            "combined_merged_bases": merged_interval_bases(combined_intervals) if combined_intervals else 0,
            "gc_json_interval_count": gc_json.get("merged_interval_count") if gc_json else None,
            "sequence_json_interval_count": seq_json.get("merged_interval_count") if seq_json else None,
        },
        "findings": [finding.__dict__ for finding in findings],
    }


def render_validation_markdown(report: dict[str, object]) -> str:
    lines = [
        "# Empirical Track Validation Report",
        "",
        f"- Status: **{report['status']}**",
        f"- Failures: **{report['fail_count']}**",
        f"- Warnings: **{report['warning_count']}**",
        "",
        "## Metrics",
        "",
        "| Metric | Value |",
        "|---|---:|",
    ]
    for key, value in report["metrics"].items():
        rendered = "n/a" if value is None else str(value)
        lines.append(f"| {key} | {rendered} |")
    lines.append("")
    lines.append("## Findings")
    lines.append("")
    if report["findings"]:
        lines.append("| Severity | Check | Message |")
        lines.append("|---|---|---|")
        for finding in report["findings"]:
            lines.append(f"| {finding['severity']} | {finding['check']} | {finding['message']} |")
    else:
        lines.append("- No validation findings.")
    lines.append("")
    return "\n".join(lines)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Validate empirical technical-track research outputs.")
    parser.add_argument("--output-dir", required=True, help="Research output directory")
    parser.add_argument("--prefix", default="demo", help="Output filename prefix")
    parser.add_argument("--out", help="Markdown validation report path")
    parser.add_argument("--json", help="JSON validation report path")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    report = validate_outputs(args.output_dir, args.prefix)
    if args.out:
        output_path = Path(args.out)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        output_path.write_text(render_validation_markdown(report), encoding="utf-8")
    if args.json:
        json_path = Path(args.json)
        json_path.parent.mkdir(parents=True, exist_ok=True)
        json_path.write_text(json.dumps(report, indent=2, sort_keys=True), encoding="utf-8")
    if not args.out and not args.json:
        print(render_validation_markdown(report))
    return 0 if report["status"] != "fail" else 1


if __name__ == "__main__":
    raise SystemExit(main())
