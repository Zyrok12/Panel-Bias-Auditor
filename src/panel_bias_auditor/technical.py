from __future__ import annotations

import json
from dataclasses import asdict, dataclass
from pathlib import Path

from .core import merged_interval_bases
from .errors import InputFormatError
from .models import Interval, normalize_chrom


VALID_BASES = {"A", "C", "G", "T"}


@dataclass(frozen=True)
class GCWindow:
    chrom: str
    start: int
    end: int
    gc_fraction: float
    callable_bases: int
    category: str


@dataclass(frozen=True)
class SequenceRiskWindow:
    chrom: str
    start: int
    end: int
    category: str
    longest_homopolymer: int
    dominant_base_fraction: float
    callable_bases: int


def parse_fasta(path: str | Path) -> dict[str, str]:
    fasta_path = Path(path)
    if not fasta_path.exists():
        raise InputFormatError(f"FASTA file does not exist: {fasta_path}")
    if not fasta_path.is_file():
        raise InputFormatError(f"FASTA path is not a file: {fasta_path}")

    sequences: dict[str, list[str]] = {}
    current_name: str | None = None
    with fasta_path.open("r", encoding="utf-8") as handle:
        for line_number, raw in enumerate(handle, start=1):
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                current_name = line[1:].split()[0]
                if not current_name:
                    raise InputFormatError(f"{fasta_path}:{line_number}: FASTA header is empty")
                sequences.setdefault(current_name, [])
                continue
            if current_name is None:
                raise InputFormatError(f"{fasta_path}:{line_number}: FASTA sequence appears before the first header")
            sequences[current_name].append(line.upper())

    parsed = {name: "".join(parts) for name, parts in sequences.items()}
    if not parsed:
        raise InputFormatError(f"FASTA file contains no sequences: {fasta_path}")
    return parsed


def sequence_for_chrom(sequences: dict[str, str], chrom: str) -> str | None:
    wanted = normalize_chrom(chrom)
    for name, sequence in sequences.items():
        if normalize_chrom(name) == wanted:
            return sequence
    return None


def gc_fraction(sequence: str) -> tuple[float, int]:
    counts = {base: 0 for base in VALID_BASES}
    for base in sequence.upper():
        if base in counts:
            counts[base] += 1
    callable_bases = sum(counts.values())
    if callable_bases == 0:
        return 0.0, 0
    return (counts["G"] + counts["C"]) / callable_bases, callable_bases


def base_counts(sequence: str) -> dict[str, int]:
    counts = {base: 0 for base in VALID_BASES}
    for base in sequence.upper():
        if base in counts:
            counts[base] += 1
    return counts


def longest_homopolymer_run(sequence: str) -> int:
    longest = 0
    current_base = ""
    current_length = 0
    for base in sequence.upper():
        if base not in VALID_BASES:
            current_base = ""
            current_length = 0
            continue
        if base == current_base:
            current_length += 1
        else:
            current_base = base
            current_length = 1
        longest = max(longest, current_length)
    return longest


def default_scan_regions(sequences: dict[str, str]) -> list[Interval]:
    return [Interval(chrom, 0, len(sequence), name=f"{chrom}:0-{len(sequence)}") for chrom, sequence in sequences.items()]


def validate_scan_region(region: Interval, sequence_length: int) -> Interval:
    start = max(0, region.start)
    end = min(region.end, sequence_length)
    if end <= start:
        raise InputFormatError(f"Region {region.chrom}:{region.start}-{region.end} is outside FASTA sequence bounds")
    return Interval(region.chrom, start, end, name=region.name, source=region.source, attrs=region.attrs)


def scan_gc_windows(
    sequences: dict[str, str],
    regions: list[Interval] | None = None,
    window_size: int = 120,
    step_size: int = 60,
    low_gc: float = 0.25,
    high_gc: float = 0.75,
    min_callable_fraction: float = 0.8,
) -> list[GCWindow]:
    if window_size <= 0:
        raise InputFormatError("window_size must be > 0")
    if step_size <= 0:
        raise InputFormatError("step_size must be > 0")
    if not 0 <= low_gc <= 1 or not 0 <= high_gc <= 1 or low_gc >= high_gc:
        raise InputFormatError("GC thresholds must satisfy 0 <= low_gc < high_gc <= 1")
    if not 0 < min_callable_fraction <= 1:
        raise InputFormatError("min_callable_fraction must satisfy 0 < value <= 1")

    scan_regions = regions or default_scan_regions(sequences)
    windows: list[GCWindow] = []
    for region in scan_regions:
        sequence = sequence_for_chrom(sequences, region.chrom)
        if sequence is None:
            raise InputFormatError(f"Region chromosome is not present in FASTA: {region.chrom}")
        bounded = validate_scan_region(region, len(sequence))
        if bounded.length < window_size:
            continue
        last_start = bounded.end - window_size
        for start in range(bounded.start, last_start + 1, step_size):
            end = start + window_size
            fraction, callable_bases = gc_fraction(sequence[start:end])
            if callable_bases / window_size < min_callable_fraction:
                continue
            if fraction <= low_gc:
                windows.append(GCWindow(bounded.chrom, start, end, fraction, callable_bases, "low_gc"))
            elif fraction >= high_gc:
                windows.append(GCWindow(bounded.chrom, start, end, fraction, callable_bases, "high_gc"))
    return windows


def scan_sequence_risk_windows(
    sequences: dict[str, str],
    regions: list[Interval] | None = None,
    window_size: int = 80,
    step_size: int = 40,
    min_homopolymer: int = 8,
    low_complexity_fraction: float = 0.85,
    min_callable_fraction: float = 0.8,
) -> list[SequenceRiskWindow]:
    if window_size <= 0:
        raise InputFormatError("window_size must be > 0")
    if step_size <= 0:
        raise InputFormatError("step_size must be > 0")
    if min_homopolymer <= 1:
        raise InputFormatError("min_homopolymer must be > 1")
    if not 0 < low_complexity_fraction <= 1:
        raise InputFormatError("low_complexity_fraction must satisfy 0 < value <= 1")
    if not 0 < min_callable_fraction <= 1:
        raise InputFormatError("min_callable_fraction must satisfy 0 < value <= 1")

    scan_regions = regions or default_scan_regions(sequences)
    windows: list[SequenceRiskWindow] = []
    for region in scan_regions:
        sequence = sequence_for_chrom(sequences, region.chrom)
        if sequence is None:
            raise InputFormatError(f"Region chromosome is not present in FASTA: {region.chrom}")
        bounded = validate_scan_region(region, len(sequence))
        if bounded.length < window_size:
            continue
        last_start = bounded.end - window_size
        for start in range(bounded.start, last_start + 1, step_size):
            end = start + window_size
            window_sequence = sequence[start:end]
            counts = base_counts(window_sequence)
            callable_bases = sum(counts.values())
            if callable_bases / window_size < min_callable_fraction:
                continue
            longest_run = longest_homopolymer_run(window_sequence)
            dominant_fraction = max(counts.values()) / callable_bases if callable_bases else 0.0
            if longest_run >= min_homopolymer:
                windows.append(
                    SequenceRiskWindow(
                        bounded.chrom,
                        start,
                        end,
                        "homopolymer",
                        longest_run,
                        dominant_fraction,
                        callable_bases,
                    )
                )
            elif dominant_fraction >= low_complexity_fraction:
                windows.append(
                    SequenceRiskWindow(
                        bounded.chrom,
                        start,
                        end,
                        "low_complexity",
                        longest_run,
                        dominant_fraction,
                        callable_bases,
                    )
                )
    return windows


def merge_gc_windows(windows: list[GCWindow]) -> list[Interval]:
    if not windows:
        return []
    sorted_windows = sorted(windows, key=lambda item: (normalize_chrom(item.chrom), item.category, item.start, item.end))
    merged: list[dict[str, object]] = []
    for window in sorted_windows:
        if (
            merged
            and normalize_chrom(str(merged[-1]["chrom"])) == normalize_chrom(window.chrom)
            and merged[-1]["category"] == window.category
            and int(merged[-1]["end"]) >= window.start
        ):
            merged[-1]["end"] = max(int(merged[-1]["end"]), window.end)
            merged[-1]["min_gc"] = min(float(merged[-1]["min_gc"]), window.gc_fraction)
            merged[-1]["max_gc"] = max(float(merged[-1]["max_gc"]), window.gc_fraction)
            merged[-1]["window_count"] = int(merged[-1]["window_count"]) + 1
            continue
        merged.append(
            {
                "chrom": window.chrom,
                "start": window.start,
                "end": window.end,
                "category": window.category,
                "min_gc": window.gc_fraction,
                "max_gc": window.gc_fraction,
                "window_count": 1,
            }
        )

    intervals: list[Interval] = []
    for item in merged:
        category = str(item["category"])
        name = f"{category}_{item['chrom']}_{item['start']}_{item['end']}"
        attrs = {
            "type": "gc_extreme",
            "gc_class": category,
            "scope": "assay_validation",
            "curation_status": "source_verified",
            "evidence_level": "assay_computable_gc",
            "source": "fasta_gcderive",
            "min_gc": f"{float(item['min_gc']):.3f}",
            "max_gc": f"{float(item['max_gc']):.3f}",
            "window_count": str(item["window_count"]),
        }
        intervals.append(
            Interval(
                str(item["chrom"]),
                int(item["start"]),
                int(item["end"]),
                name=name,
                source="gcderive",
                attrs=attrs,
            )
        )
    return intervals


def merge_sequence_risk_windows(windows: list[SequenceRiskWindow]) -> list[Interval]:
    if not windows:
        return []
    sorted_windows = sorted(windows, key=lambda item: (normalize_chrom(item.chrom), item.category, item.start, item.end))
    merged: list[dict[str, object]] = []
    for window in sorted_windows:
        if (
            merged
            and normalize_chrom(str(merged[-1]["chrom"])) == normalize_chrom(window.chrom)
            and merged[-1]["category"] == window.category
            and int(merged[-1]["end"]) >= window.start
        ):
            merged[-1]["end"] = max(int(merged[-1]["end"]), window.end)
            merged[-1]["max_homopolymer"] = max(int(merged[-1]["max_homopolymer"]), window.longest_homopolymer)
            merged[-1]["max_dominant_base_fraction"] = max(
                float(merged[-1]["max_dominant_base_fraction"]), window.dominant_base_fraction
            )
            merged[-1]["window_count"] = int(merged[-1]["window_count"]) + 1
            continue
        merged.append(
            {
                "chrom": window.chrom,
                "start": window.start,
                "end": window.end,
                "category": window.category,
                "max_homopolymer": window.longest_homopolymer,
                "max_dominant_base_fraction": window.dominant_base_fraction,
                "window_count": 1,
            }
        )

    intervals: list[Interval] = []
    for item in merged:
        category = str(item["category"])
        name = f"{category}_{item['chrom']}_{item['start']}_{item['end']}"
        attrs = {
            "type": category,
            "sequence_risk_class": category,
            "scope": "assay_validation",
            "curation_status": "source_verified",
            "evidence_level": "assay_computable_sequence_complexity",
            "source": "fasta_seqderive",
            "max_homopolymer": str(item["max_homopolymer"]),
            "max_dominant_base_fraction": f"{float(item['max_dominant_base_fraction']):.3f}",
            "window_count": str(item["window_count"]),
        }
        intervals.append(
            Interval(
                str(item["chrom"]),
                int(item["start"]),
                int(item["end"]),
                name=name,
                source="seqderive",
                attrs=attrs,
            )
        )
    return intervals


def attrs_to_bed_field(attrs: dict[str, str]) -> str:
    return ",".join(f"{key}={value}" for key, value in attrs.items())


def write_gc_bed(intervals: list[Interval], path: str | Path) -> None:
    output_path = Path(path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", encoding="utf-8", newline="\n") as handle:
        handle.write("# Generated by panel-bias-auditor gcderive from FASTA sequence.\n")
        for interval in sorted(intervals, key=lambda item: (normalize_chrom(item.chrom), item.start, item.end, item.name)):
            handle.write(
                f"{interval.chrom}\t{interval.start}\t{interval.end}\t{interval.name}\t{attrs_to_bed_field(interval.attrs)}\n"
            )


def write_sequence_risk_bed(intervals: list[Interval], path: str | Path) -> None:
    output_path = Path(path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", encoding="utf-8", newline="\n") as handle:
        handle.write("# Generated by panel-bias-auditor seqderive from FASTA sequence.\n")
        for interval in sorted(intervals, key=lambda item: (normalize_chrom(item.chrom), item.start, item.end, item.name)):
            handle.write(
                f"{interval.chrom}\t{interval.start}\t{interval.end}\t{interval.name}\t{attrs_to_bed_field(interval.attrs)}\n"
            )


def build_gc_report(
    fasta_path: str | Path,
    windows: list[GCWindow],
    intervals: list[Interval],
    window_size: int,
    step_size: int,
    low_gc: float,
    high_gc: float,
    regions: list[Interval] | None = None,
) -> dict[str, object]:
    by_class: dict[str, dict[str, object]] = {}
    for interval in intervals:
        gc_class = interval.attrs.get("gc_class", "unknown")
        bucket = by_class.setdefault(gc_class, {"interval_count": 0, "merged_bases": 0})
        bucket["interval_count"] = int(bucket["interval_count"]) + 1
    for gc_class, bucket in by_class.items():
        bucket_intervals = [interval for interval in intervals if interval.attrs.get("gc_class") == gc_class]
        bucket["merged_bases"] = merged_interval_bases(bucket_intervals)
    return {
        "fasta": str(fasta_path),
        "window_size": window_size,
        "step_size": step_size,
        "low_gc": low_gc,
        "high_gc": high_gc,
        "regions_scanned": len(regions) if regions is not None else "all_contigs",
        "extreme_window_count": len(windows),
        "merged_interval_count": len(intervals),
        "merged_bases": merged_interval_bases(intervals),
        "by_class": dict(sorted(by_class.items())),
        "intervals": [
            {
                "chrom": interval.chrom,
                "start": interval.start,
                "end": interval.end,
                "name": interval.name,
                "attrs": interval.attrs,
            }
            for interval in intervals
        ],
        "windows": [asdict(window) for window in windows],
    }


def build_sequence_risk_report(
    fasta_path: str | Path,
    windows: list[SequenceRiskWindow],
    intervals: list[Interval],
    window_size: int,
    step_size: int,
    min_homopolymer: int,
    low_complexity_fraction: float,
    regions: list[Interval] | None = None,
) -> dict[str, object]:
    by_class: dict[str, dict[str, object]] = {}
    for interval in intervals:
        risk_class = interval.attrs.get("sequence_risk_class", "unknown")
        bucket = by_class.setdefault(risk_class, {"interval_count": 0, "merged_bases": 0})
        bucket["interval_count"] = int(bucket["interval_count"]) + 1
    for risk_class, bucket in by_class.items():
        bucket_intervals = [interval for interval in intervals if interval.attrs.get("sequence_risk_class") == risk_class]
        bucket["merged_bases"] = merged_interval_bases(bucket_intervals)
    return {
        "fasta": str(fasta_path),
        "window_size": window_size,
        "step_size": step_size,
        "min_homopolymer": min_homopolymer,
        "low_complexity_fraction": low_complexity_fraction,
        "regions_scanned": len(regions) if regions is not None else "all_contigs",
        "risk_window_count": len(windows),
        "merged_interval_count": len(intervals),
        "merged_bases": merged_interval_bases(intervals),
        "by_class": dict(sorted(by_class.items())),
        "intervals": [
            {
                "chrom": interval.chrom,
                "start": interval.start,
                "end": interval.end,
                "name": interval.name,
                "attrs": interval.attrs,
            }
            for interval in intervals
        ],
        "windows": [asdict(window) for window in windows],
    }


def render_gc_report_markdown(report: dict[str, object]) -> str:
    lines: list[str] = []
    lines.append("# Empirical GC Track Report")
    lines.append("")
    lines.append(f"- FASTA: `{report['fasta']}`")
    lines.append(f"- Window size: **{report['window_size']} bp**")
    lines.append(f"- Step size: **{report['step_size']} bp**")
    lines.append(f"- Low-GC threshold: **{float(report['low_gc']):.2f}**")
    lines.append(f"- High-GC threshold: **{float(report['high_gc']):.2f}**")
    lines.append(f"- Regions scanned: **{report['regions_scanned']}**")
    lines.append(f"- Extreme windows: **{report['extreme_window_count']}**")
    lines.append(f"- Merged intervals: **{report['merged_interval_count']}**")
    lines.append(f"- Merged bases: **{int(report['merged_bases']):,}**")
    lines.append("")
    lines.append("## By GC Class")
    lines.append("")
    lines.append("| Class | Intervals | Merged bases |")
    lines.append("|---|---:|---:|")
    for gc_class, bucket in report["by_class"].items():
        lines.append(f"| {gc_class} | {bucket['interval_count']} | {int(bucket['merged_bases']):,} |")
    if not report["by_class"]:
        lines.append("| none | 0 | 0 |")
    lines.append("")
    lines.append("## Largest Intervals")
    lines.append("")
    lines.append("| Region | Class | Locus | Length | GC range |")
    lines.append("|---|---|---|---:|---|")
    intervals = sorted(report["intervals"], key=lambda item: (int(item["end"]) - int(item["start"])), reverse=True)
    for item in intervals[:15]:
        attrs = item["attrs"]
        locus = f"{item['chrom']}:{item['start']}-{item['end']}"
        length = int(item["end"]) - int(item["start"])
        gc_range = f"{attrs.get('min_gc', '')}-{attrs.get('max_gc', '')}"
        lines.append(f"| {item['name']} | {attrs.get('gc_class', '')} | {locus} | {length:,} | {gc_range} |")
    if not intervals:
        lines.append("| none | none | none | 0 | none |")
    lines.append("")
    return "\n".join(lines)


def render_sequence_risk_report_markdown(report: dict[str, object]) -> str:
    lines: list[str] = []
    lines.append("# Empirical Sequence Risk Track Report")
    lines.append("")
    lines.append(f"- FASTA: `{report['fasta']}`")
    lines.append(f"- Window size: **{report['window_size']} bp**")
    lines.append(f"- Step size: **{report['step_size']} bp**")
    lines.append(f"- Homopolymer threshold: **{report['min_homopolymer']} bp**")
    lines.append(f"- Low-complexity dominant-base threshold: **{float(report['low_complexity_fraction']):.2f}**")
    lines.append(f"- Regions scanned: **{report['regions_scanned']}**")
    lines.append(f"- Risk windows: **{report['risk_window_count']}**")
    lines.append(f"- Merged intervals: **{report['merged_interval_count']}**")
    lines.append(f"- Merged bases: **{int(report['merged_bases']):,}**")
    lines.append("")
    lines.append("## By Sequence Risk Class")
    lines.append("")
    lines.append("| Class | Intervals | Merged bases |")
    lines.append("|---|---:|---:|")
    for risk_class, bucket in report["by_class"].items():
        lines.append(f"| {risk_class} | {bucket['interval_count']} | {int(bucket['merged_bases']):,} |")
    if not report["by_class"]:
        lines.append("| none | 0 | 0 |")
    lines.append("")
    lines.append("## Largest Intervals")
    lines.append("")
    lines.append("| Region | Class | Locus | Length | Max homopolymer | Max dominant base |")
    lines.append("|---|---|---|---:|---:|---:|")
    intervals = sorted(report["intervals"], key=lambda item: (int(item["end"]) - int(item["start"])), reverse=True)
    for item in intervals[:15]:
        attrs = item["attrs"]
        locus = f"{item['chrom']}:{item['start']}-{item['end']}"
        length = int(item["end"]) - int(item["start"])
        lines.append(
            f"| {item['name']} | {attrs.get('sequence_risk_class', '')} | {locus} | {length:,} | {attrs.get('max_homopolymer', '')} | {attrs.get('max_dominant_base_fraction', '')} |"
        )
    if not intervals:
        lines.append("| none | none | none | 0 | 0 | 0 |")
    lines.append("")
    return "\n".join(lines)


def write_gc_report_markdown(report: dict[str, object], path: str | Path) -> None:
    output_path = Path(path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(render_gc_report_markdown(report), encoding="utf-8")


def write_gc_report_json(report: dict[str, object], path: str | Path) -> None:
    output_path = Path(path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(json.dumps(report, indent=2, sort_keys=True), encoding="utf-8")


def write_sequence_risk_report_markdown(report: dict[str, object], path: str | Path) -> None:
    output_path = Path(path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(render_sequence_risk_report_markdown(report), encoding="utf-8")


def write_sequence_risk_report_json(report: dict[str, object], path: str | Path) -> None:
    output_path = Path(path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(json.dumps(report, indent=2, sort_keys=True), encoding="utf-8")
