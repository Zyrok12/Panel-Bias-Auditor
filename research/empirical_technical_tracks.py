from __future__ import annotations

import argparse
import json
from pathlib import Path

from panel_bias_auditor.core import merged_interval_bases
from panel_bias_auditor.parsers import parse_bed
from panel_bias_auditor.technical import (
    build_gc_report,
    build_sequence_risk_report,
    merge_gc_windows,
    merge_sequence_risk_windows,
    parse_fasta,
    scan_gc_windows,
    scan_sequence_risk_windows,
    write_gc_bed,
    write_gc_report_json,
    write_gc_report_markdown,
    write_sequence_risk_bed,
    write_sequence_risk_report_json,
    write_sequence_risk_report_markdown,
)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Research runner for empirical assay technical-risk tracks from FASTA sequence."
    )
    parser.add_argument("--fasta", required=True, help="Reference or targeted FASTA")
    parser.add_argument("--regions", help="Optional BED regions to scan")
    parser.add_argument("--output-dir", required=True, help="Output directory")
    parser.add_argument("--prefix", default="empirical", help="Output filename prefix")
    parser.add_argument("--gc-window", type=int, default=120)
    parser.add_argument("--gc-step", type=int, default=60)
    parser.add_argument("--low-gc", type=float, default=0.25)
    parser.add_argument("--high-gc", type=float, default=0.75)
    parser.add_argument("--seq-window", type=int, default=80)
    parser.add_argument("--seq-step", type=int, default=40)
    parser.add_argument("--min-homopolymer", type=int, default=8)
    parser.add_argument("--low-complexity-fraction", type=float, default=0.85)
    return parser


def attrs_to_bed_field(attrs: dict[str, str]) -> str:
    return ",".join(f"{key}={value}" for key, value in attrs.items())


def write_merged_bed(gc_bed: Path, seq_bed: Path, output_path: Path) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", encoding="utf-8", newline="\n") as handle:
        handle.write("# Combined empirical technical-risk track.\n")
        for source_path in (gc_bed, seq_bed):
            for line in source_path.read_text(encoding="utf-8").splitlines():
                if not line or line.startswith("#"):
                    continue
                handle.write(line + "\n")


def render_summary(summary: dict[str, object]) -> str:
    lines = [
        "# Empirical Technical Track Research Summary",
        "",
        f"- FASTA: `{summary['fasta']}`",
        f"- Regions: `{summary['regions']}`",
        f"- Combined BED: `{summary['combined_bed']}`",
        f"- GC intervals: **{summary['gc_interval_count']}**",
        f"- Sequence-risk intervals: **{summary['sequence_risk_interval_count']}**",
        f"- Combined raw intervals: **{summary['combined_interval_count']}**",
        f"- Combined merged bases: **{int(summary['combined_merged_bases']):,}**",
        "",
        "## Outputs",
        "",
    ]
    for label, path in summary["outputs"].items():
        lines.append(f"- {label}: `{path}`")
    lines.append("")
    lines.append("## Research Caveat")
    lines.append("")
    lines.append(
        "These tracks are empirical sequence-risk annotations. They identify regions that deserve assay-validation review; they do not prove a region is uncallable."
    )
    lines.append("")
    return "\n".join(lines)


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    sequences = parse_fasta(args.fasta)
    regions = parse_bed(args.regions, source="research_regions") if args.regions else None

    gc_windows = scan_gc_windows(
        sequences,
        regions=regions,
        window_size=args.gc_window,
        step_size=args.gc_step,
        low_gc=args.low_gc,
        high_gc=args.high_gc,
    )
    gc_intervals = merge_gc_windows(gc_windows)
    gc_bed = output_dir / f"{args.prefix}_gc_extremes.bed"
    gc_md = output_dir / f"{args.prefix}_gc_report.md"
    gc_json = output_dir / f"{args.prefix}_gc_report.json"
    write_gc_bed(gc_intervals, gc_bed)
    gc_report = build_gc_report(
        fasta_path=args.fasta,
        windows=gc_windows,
        intervals=gc_intervals,
        window_size=args.gc_window,
        step_size=args.gc_step,
        low_gc=args.low_gc,
        high_gc=args.high_gc,
        regions=regions,
    )
    write_gc_report_markdown(gc_report, gc_md)
    write_gc_report_json(gc_report, gc_json)

    seq_windows = scan_sequence_risk_windows(
        sequences,
        regions=regions,
        window_size=args.seq_window,
        step_size=args.seq_step,
        min_homopolymer=args.min_homopolymer,
        low_complexity_fraction=args.low_complexity_fraction,
    )
    seq_intervals = merge_sequence_risk_windows(seq_windows)
    seq_bed = output_dir / f"{args.prefix}_sequence_risk.bed"
    seq_md = output_dir / f"{args.prefix}_sequence_risk_report.md"
    seq_json = output_dir / f"{args.prefix}_sequence_risk_report.json"
    write_sequence_risk_bed(seq_intervals, seq_bed)
    seq_report = build_sequence_risk_report(
        fasta_path=args.fasta,
        windows=seq_windows,
        intervals=seq_intervals,
        window_size=args.seq_window,
        step_size=args.seq_step,
        min_homopolymer=args.min_homopolymer,
        low_complexity_fraction=args.low_complexity_fraction,
        regions=regions,
    )
    write_sequence_risk_report_markdown(seq_report, seq_md)
    write_sequence_risk_report_json(seq_report, seq_json)

    combined_bed = output_dir / f"{args.prefix}_combined_technical_risk.bed"
    write_merged_bed(gc_bed, seq_bed, combined_bed)
    combined_intervals = gc_intervals + seq_intervals
    summary = {
        "fasta": args.fasta,
        "regions": args.regions or "all_contigs",
        "combined_bed": str(combined_bed),
        "gc_interval_count": len(gc_intervals),
        "sequence_risk_interval_count": len(seq_intervals),
        "combined_interval_count": len(combined_intervals),
        "combined_merged_bases": merged_interval_bases(combined_intervals),
        "outputs": {
            "gc_bed": str(gc_bed),
            "gc_markdown": str(gc_md),
            "gc_json": str(gc_json),
            "sequence_risk_bed": str(seq_bed),
            "sequence_risk_markdown": str(seq_md),
            "sequence_risk_json": str(seq_json),
            "combined_bed": str(combined_bed),
        },
    }
    summary_md = output_dir / f"{args.prefix}_research_summary.md"
    summary_json = output_dir / f"{args.prefix}_research_summary.json"
    summary_md.write_text(render_summary(summary), encoding="utf-8")
    summary_json.write_text(json.dumps(summary, indent=2, sort_keys=True), encoding="utf-8")
    print(render_summary(summary))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
