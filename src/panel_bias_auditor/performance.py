from __future__ import annotations

import csv
import json
from dataclasses import dataclass
from pathlib import Path

from .core import merged_interval_bases, overlap_bases
from .errors import InputFormatError
from .models import Interval, normalize_chrom
from .parsers import parse_bed


REQUIRED_PERFORMANCE_COLUMNS = {"chrom", "start", "end", "name"}
FAIL_STATUSES = {"fail", "failed", "uncallable", "low_coverage", "low_callability", "no_call"}
PERFORMANCE_COLUMNS = [
    "chrom",
    "start",
    "end",
    "name",
    "mean_depth",
    "callable_fraction",
    "status",
    "false_negative_count",
    "false_positive_count",
    "total_variant_count",
]


@dataclass(frozen=True)
class PerformanceInterval:
    interval: Interval
    mean_depth: float | None = None
    callable_fraction: float | None = None
    status: str = ""
    false_negative_count: int = 0
    false_positive_count: int = 0
    total_variant_count: int = 0

    def is_failure(self, min_depth: float | None = None, min_callable_fraction: float | None = None) -> bool:
        status_value = self.status.strip().lower()
        if status_value in FAIL_STATUSES:
            return True
        if min_depth is not None and self.mean_depth is not None and self.mean_depth < min_depth:
            return True
        if (
            min_callable_fraction is not None
            and self.callable_fraction is not None
            and self.callable_fraction < min_callable_fraction
        ):
            return True
        return self.false_negative_count > 0 or self.false_positive_count > 0


def detect_delimiter(path: Path) -> str:
    first_line = path.read_text(encoding="utf-8").splitlines()[0]
    return "\t" if "\t" in first_line else ","


def parse_optional_float(raw: str | None, path: Path, line_number: int, column: str) -> float | None:
    if raw is None or raw.strip() == "":
        return None
    try:
        return float(raw)
    except ValueError as exc:
        raise InputFormatError(f"{path}:{line_number}: {column} must be numeric") from exc


def parse_optional_int(raw: str | None, path: Path, line_number: int, column: str) -> int:
    if raw is None or raw.strip() == "":
        return 0
    try:
        return int(raw)
    except ValueError as exc:
        raise InputFormatError(f"{path}:{line_number}: {column} must be an integer") from exc


def parse_performance_table(path: str | Path) -> list[PerformanceInterval]:
    source_path = Path(path)
    if not source_path.exists():
        raise InputFormatError(f"Assay performance table does not exist: {source_path}")
    if not source_path.is_file():
        raise InputFormatError(f"Assay performance path is not a file: {source_path}")
    delimiter = detect_delimiter(source_path)
    rows: list[PerformanceInterval] = []
    with source_path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter=delimiter)
        missing = REQUIRED_PERFORMANCE_COLUMNS - set(reader.fieldnames or [])
        if missing:
            raise InputFormatError(f"Assay performance table is missing required columns: {', '.join(sorted(missing))}")
        for line_number, row in enumerate(reader, start=2):
            chrom = (row.get("chrom") or "").strip()
            name = (row.get("name") or "").strip()
            if not chrom or not name:
                raise InputFormatError(f"{source_path}:{line_number}: chrom and name are required")
            try:
                start = int((row.get("start") or "").strip())
                end = int((row.get("end") or "").strip())
            except ValueError as exc:
                raise InputFormatError(f"{source_path}:{line_number}: start/end must be integers") from exc
            if start < 0 or end <= start:
                raise InputFormatError(f"{source_path}:{line_number}: invalid interval {chrom}:{start}-{end}")
            rows.append(
                PerformanceInterval(
                    interval=Interval(chrom, start, end, name=name, source="assay_performance"),
                    mean_depth=parse_optional_float(row.get("mean_depth"), source_path, line_number, "mean_depth"),
                    callable_fraction=parse_optional_float(
                        row.get("callable_fraction"), source_path, line_number, "callable_fraction"
                    ),
                    status=(row.get("status") or "").strip(),
                    false_negative_count=parse_optional_int(
                        row.get("false_negative_count"), source_path, line_number, "false_negative_count"
                    ),
                    false_positive_count=parse_optional_int(
                        row.get("false_positive_count"), source_path, line_number, "false_positive_count"
                    ),
                    total_variant_count=parse_optional_int(
                        row.get("total_variant_count"), source_path, line_number, "total_variant_count"
                    ),
                )
            )
    if not rows:
        raise InputFormatError(f"Assay performance table contains no records: {source_path}")
    return rows


def performance_rows_from_callability(
    evaluated_regions: list[Interval],
    callable_regions: list[Interval],
    min_callable_fraction: float = 0.95,
) -> list[PerformanceInterval]:
    if not 0 <= min_callable_fraction <= 1:
        raise InputFormatError("min_callable_fraction must be between 0 and 1")
    if not evaluated_regions:
        raise InputFormatError("At least one evaluated region is required")
    if not callable_regions:
        raise InputFormatError("At least one callable/high-confidence region is required")

    rows: list[PerformanceInterval] = []
    for region in evaluated_regions:
        callable_bases = overlap_bases(region, callable_regions)
        callable_fraction = callable_bases / region.length if region.length else 0.0
        status = "pass" if callable_fraction >= min_callable_fraction else "low_callability"
        rows.append(
            PerformanceInterval(
                interval=Interval(region.chrom, region.start, region.end, name=region.name, source="public_benchmark"),
                callable_fraction=callable_fraction,
                status=status,
            )
        )
    return rows


def write_performance_table(rows: list[PerformanceInterval], path: str | Path) -> None:
    output_path = Path(path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=PERFORMANCE_COLUMNS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow(
                {
                    "chrom": row.interval.chrom,
                    "start": row.interval.start,
                    "end": row.interval.end,
                    "name": row.interval.name,
                    "mean_depth": "" if row.mean_depth is None else f"{row.mean_depth:.6g}",
                    "callable_fraction": "" if row.callable_fraction is None else f"{row.callable_fraction:.6f}",
                    "status": row.status,
                    "false_negative_count": row.false_negative_count,
                    "false_positive_count": row.false_positive_count,
                    "total_variant_count": row.total_variant_count,
                }
            )


def summarize_callability_rows(rows: list[PerformanceInterval]) -> dict[str, object]:
    if not rows:
        raise InputFormatError("No callability-derived performance rows were generated")
    callable_fractions = [row.callable_fraction or 0.0 for row in rows]
    failure_rows = [row for row in rows if row.is_failure()]
    total_bases = sum(row.interval.length for row in rows)
    failed_bases = sum(row.interval.length for row in failure_rows)
    return {
        "interval_count": len(rows),
        "failure_interval_count": len(failure_rows),
        "pass_interval_count": len(rows) - len(failure_rows),
        "total_bases": total_bases,
        "failed_bases": failed_bases,
        "mean_callable_fraction": sum(callable_fractions) / len(callable_fractions),
        "min_callable_fraction": min(callable_fractions),
        "max_callable_fraction": max(callable_fractions),
    }


def render_callability_conversion_markdown(
    rows: list[PerformanceInterval],
    source_label: str,
    min_callable_fraction: float,
) -> str:
    summary = summarize_callability_rows(rows)
    lines = [
        "# Public Benchmark Callability Conversion",
        "",
        f"- Source label: **{source_label}**",
        f"- Evaluated intervals: **{summary['interval_count']}**",
        f"- Low-callability intervals: **{summary['failure_interval_count']}**",
        f"- Total evaluated bases: **{int(summary['total_bases']):,}**",
        f"- Low-callability bases: **{int(summary['failed_bases']):,}**",
        f"- Failure threshold: callable_fraction < **{min_callable_fraction:.3f}**",
        f"- Mean callable fraction: **{float(summary['mean_callable_fraction']):.3f}**",
        "",
        "## Interval Preview",
        "",
        "| Region | Locus | Callable fraction | Status |",
        "|---|---|---:|---|",
    ]
    sorted_rows = sorted(rows, key=lambda item: ((item.callable_fraction or 0.0), item.interval.chrom, item.interval.start))
    for row in sorted_rows[:20]:
        lines.append(
            f"| {row.interval.name} | {row.interval.chrom}:{row.interval.start}-{row.interval.end} | "
            f"{(row.callable_fraction or 0.0):.3f} | {row.status} |"
        )
    lines.extend(
        [
            "",
            "## Interpretation",
            "",
            "This table converts public benchmark high-confidence or callable regions into the same interval-level schema used for assay-performance enrichment. A low-callability row means the evaluated region is poorly represented by the benchmark/callable set, not that a specific clinical assay failed there.",
            "",
        ]
    )
    return "\n".join(lines)


def clip_risk_to_universe(risk_intervals: list[Interval], universe: list[Interval]) -> list[Interval]:
    clipped: list[Interval] = []
    by_chrom: dict[str, list[Interval]] = {}
    for region in universe:
        by_chrom.setdefault(region.norm_chrom, []).append(region)
    for risk in risk_intervals:
        for region in by_chrom.get(risk.norm_chrom, []):
            span = risk.overlap_span(region)
            if span is None:
                continue
            clipped.append(
                Interval(
                    risk.chrom,
                    span[0],
                    span[1],
                    name=risk.name,
                    source=risk.source,
                    attrs=risk.attrs,
                )
            )
    return clipped


def category_for(interval: Interval) -> str:
    return (
        interval.attrs.get("type")
        or interval.attrs.get("gc_class")
        or interval.attrs.get("sequence_risk_class")
        or "technical_risk"
    )


def enrichment_from_counts(risk_fail: int, risk_pass: int, nonrisk_fail: int, nonrisk_pass: int) -> dict[str, object]:
    risk_total = risk_fail + risk_pass
    nonrisk_total = nonrisk_fail + nonrisk_pass
    risk_failure_rate = risk_fail / risk_total if risk_total else None
    nonrisk_failure_rate = nonrisk_fail / nonrisk_total if nonrisk_total else None
    odds_ratio = ((risk_fail + 0.5) * (nonrisk_pass + 0.5)) / ((risk_pass + 0.5) * (nonrisk_fail + 0.5))
    relative_risk = None
    if risk_failure_rate is not None and nonrisk_failure_rate not in (None, 0):
        relative_risk = risk_failure_rate / nonrisk_failure_rate
    return {
        "risk_fail_bases": risk_fail,
        "risk_pass_bases": risk_pass,
        "nonrisk_fail_bases": nonrisk_fail,
        "nonrisk_pass_bases": nonrisk_pass,
        "risk_failure_rate": risk_failure_rate,
        "nonrisk_failure_rate": nonrisk_failure_rate,
        "odds_ratio_haldane": odds_ratio,
        "relative_risk": relative_risk,
    }


def analyze_performance_enrichment(
    performance_rows: list[PerformanceInterval],
    risk_intervals: list[Interval],
    min_depth: float | None = None,
    min_callable_fraction: float | None = None,
) -> dict[str, object]:
    universe = [row.interval for row in performance_rows]
    failures = [row.interval for row in performance_rows if row.is_failure(min_depth, min_callable_fraction)]
    passes = [row.interval for row in performance_rows if not row.is_failure(min_depth, min_callable_fraction)]
    universe_bases = merged_interval_bases(universe)
    failure_bases = merged_interval_bases(failures)
    risk_in_universe = clip_risk_to_universe(risk_intervals, universe)
    risk_bases = merged_interval_bases(risk_in_universe)
    risk_fail_bases = sum(overlap_bases(failure, risk_in_universe) for failure in failures)
    risk_pass_bases = max(0, risk_bases - risk_fail_bases)
    nonrisk_fail_bases = max(0, failure_bases - risk_fail_bases)
    nonrisk_pass_bases = max(0, universe_bases - risk_fail_bases - risk_pass_bases - nonrisk_fail_bases)

    category_results: dict[str, dict[str, object]] = {}
    for category in sorted({category_for(interval) for interval in risk_in_universe}):
        category_intervals = [interval for interval in risk_in_universe if category_for(interval) == category]
        category_risk_bases = merged_interval_bases(category_intervals)
        category_fail = sum(overlap_bases(failure, category_intervals) for failure in failures)
        category_pass = max(0, category_risk_bases - category_fail)
        category_results[category] = {
            "risk_bases": category_risk_bases,
            **enrichment_from_counts(
                category_fail,
                category_pass,
                max(0, failure_bases - category_fail),
                max(0, universe_bases - failure_bases - category_pass),
            ),
        }

    row_summaries = []
    for row in performance_rows:
        overlap = overlap_bases(row.interval, risk_in_universe)
        row_summaries.append(
            {
                "name": row.interval.name,
                "chrom": row.interval.chrom,
                "start": row.interval.start,
                "end": row.interval.end,
                "length": row.interval.length,
                "mean_depth": row.mean_depth,
                "callable_fraction": row.callable_fraction,
                "status": row.status,
                "false_negative_count": row.false_negative_count,
                "false_positive_count": row.false_positive_count,
                "is_failure": row.is_failure(min_depth, min_callable_fraction),
                "risk_overlap_bases": overlap,
                "risk_overlap_fraction": overlap / row.interval.length if row.interval.length else 0.0,
            }
        )
    row_summaries.sort(key=lambda item: (not bool(item["is_failure"]), -int(item["risk_overlap_bases"]), str(item["name"])))

    return {
        "parameters": {
            "min_depth": min_depth,
            "min_callable_fraction": min_callable_fraction,
        },
        "summary": {
            "performance_interval_count": len(performance_rows),
            "failure_interval_count": len(failures),
            "pass_interval_count": len(passes),
            "universe_bases": universe_bases,
            "failure_bases": failure_bases,
            "risk_bases_in_universe": risk_bases,
        },
        "overall_enrichment": enrichment_from_counts(
            risk_fail_bases,
            risk_pass_bases,
            nonrisk_fail_bases,
            nonrisk_pass_bases,
        ),
        "by_risk_category": category_results,
        "performance_rows": row_summaries,
    }


def render_performance_markdown(report: dict[str, object]) -> str:
    summary = report["summary"]
    overall = report["overall_enrichment"]
    lines = [
        "# Assay Performance Enrichment Report",
        "",
        f"- Performance intervals: **{summary['performance_interval_count']}**",
        f"- Failure intervals: **{summary['failure_interval_count']}**",
        f"- Universe bases: **{int(summary['universe_bases']):,}**",
        f"- Failure bases: **{int(summary['failure_bases']):,}**",
        f"- Technical-risk bases in universe: **{int(summary['risk_bases_in_universe']):,}**",
        "",
        "## Overall Enrichment",
        "",
        "| Metric | Value |",
        "|---|---:|",
    ]
    for key, value in overall.items():
        if isinstance(value, float):
            rendered = f"{value:.4f}"
        elif value is None:
            rendered = "n/a"
        else:
            rendered = f"{int(value):,}"
        lines.append(f"| {key} | {rendered} |")
    lines.extend(["", "## By Risk Category", "", "| Category | Risk bases | Risk fail bases | Odds ratio | Relative risk |", "|---|---:|---:|---:|---:|"])
    for category, item in report["by_risk_category"].items():
        odds = item["odds_ratio_haldane"]
        rr = item["relative_risk"]
        lines.append(
            f"| {category} | {int(item['risk_bases']):,} | {int(item['risk_fail_bases']):,} | {odds:.4f} | {'n/a' if rr is None else f'{rr:.4f}'} |"
        )
    if not report["by_risk_category"]:
        lines.append("| none | 0 | 0 | n/a | n/a |")
    lines.extend(["", "## Top Performance Rows", "", "| Row | Locus | Failure | Mean depth | Callable fraction | Risk overlap |", "|---|---|---|---:|---:|---:|"])
    for row in report["performance_rows"][:20]:
        locus = f"{row['chrom']}:{row['start']}-{row['end']}"
        depth = "n/a" if row["mean_depth"] is None else f"{float(row['mean_depth']):.2f}"
        callable_fraction = "n/a" if row["callable_fraction"] is None else f"{float(row['callable_fraction']):.3f}"
        lines.append(
            f"| {row['name']} | {locus} | {row['is_failure']} | {depth} | {callable_fraction} | {int(row['risk_overlap_bases']):,} |"
        )
    lines.extend(
        [
            "",
            "## Interpretation",
            "",
            "This analysis tests whether empirical technical-risk intervals overlap observed assay failures more than expected over evaluated assay regions. It is an exploratory research statistic, not proof of causal failure.",
            "",
        ]
    )
    return "\n".join(lines)


def write_performance_report(report: dict[str, object], markdown_path: str | Path, json_path: str | Path) -> None:
    markdown_output = Path(markdown_path)
    json_output = Path(json_path)
    markdown_output.parent.mkdir(parents=True, exist_ok=True)
    json_output.parent.mkdir(parents=True, exist_ok=True)
    markdown_output.write_text(render_performance_markdown(report), encoding="utf-8")
    json_output.write_text(json.dumps(report, indent=2, sort_keys=True), encoding="utf-8")


def load_risk_bed(path: str | Path) -> list[Interval]:
    return parse_bed(path, source="technical_risk")
