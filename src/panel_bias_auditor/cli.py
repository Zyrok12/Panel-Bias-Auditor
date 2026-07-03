from __future__ import annotations

import argparse
import sys
from pathlib import Path

from .batch import (
    compare_panels,
    load_panel_manifest,
    render_compare_markdown,
    write_compare_json,
    write_compare_markdown,
    write_compare_tsv,
)
from .core import audit_panel, render_markdown, write_html_report, write_json_report, write_markdown_report
from .errors import AuditorError
from .parsers import parse_bed, parse_vcf
from .technical import (
    build_gc_report,
    build_sequence_risk_report,
    merge_gc_windows,
    merge_sequence_risk_windows,
    parse_fasta,
    render_gc_report_markdown,
    render_sequence_risk_report_markdown,
    scan_gc_windows,
    scan_sequence_risk_windows,
    write_gc_bed,
    write_gc_report_json,
    write_gc_report_markdown,
    write_sequence_risk_bed,
    write_sequence_risk_report_json,
    write_sequence_risk_report_markdown,
)
from .tracks import load_track_manifest, load_tracks_from_bundle


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="panel-bias-auditor",
        description="Audit genomic panel blind spots against critical and difficult region tracks.",
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    audit = subparsers.add_parser("audit", help="Run a panel audit")
    audit.add_argument("--panel", required=True, help="Panel BED file")
    audit.add_argument("--critical-regions", help="Critical regions BED file")
    audit.add_argument("--difficult-regions", help="Difficult regions BED file")
    audit.add_argument("--track-manifest", help="JSON track bundle manifest")
    audit.add_argument("--variants", help="Optional VCF file to annotate")
    audit.add_argument("--genome-build", default="unknown", help="Genome build label for report metadata")
    audit.add_argument("--panel-name", help="Optional panel name for report metadata")
    audit.add_argument("--out", help="Markdown report path")
    audit.add_argument("--json", help="JSON report path")
    audit.add_argument("--html", help="Standalone HTML report path")
    audit.add_argument("--stdout", action="store_true", help="Print Markdown report to stdout")

    compare = subparsers.add_parser("compare", help="Compare multiple panels under the same tracks")
    compare.add_argument("--panels", required=True, help="CSV/TSV manifest with name,path columns")
    compare.add_argument("--critical-regions", help="Critical regions BED file")
    compare.add_argument("--difficult-regions", help="Difficult regions BED file")
    compare.add_argument("--track-manifest", help="JSON track bundle manifest")
    compare.add_argument("--genome-build", default="unknown", help="Genome build label for report metadata")
    compare.add_argument("--out", help="Markdown comparison report path")
    compare.add_argument("--json", help="JSON comparison report path")
    compare.add_argument("--tsv", help="TSV comparison table path")
    compare.add_argument("--stdout", action="store_true", help="Print Markdown comparison to stdout")

    gcderive = subparsers.add_parser("gcderive", help="Derive empirical GC-extreme difficult regions from FASTA")
    gcderive.add_argument("--fasta", required=True, help="Reference or targeted FASTA file")
    gcderive.add_argument("--regions", help="Optional BED regions to scan instead of all FASTA contigs")
    gcderive.add_argument("--window-size", type=int, default=120, help="Sliding window size in bp")
    gcderive.add_argument("--step-size", type=int, default=60, help="Sliding window step size in bp")
    gcderive.add_argument("--low-gc", type=float, default=0.25, help="Low-GC threshold, inclusive")
    gcderive.add_argument("--high-gc", type=float, default=0.75, help="High-GC threshold, inclusive")
    gcderive.add_argument(
        "--min-callable-fraction",
        type=float,
        default=0.8,
        help="Minimum fraction of A/C/G/T bases required inside a window",
    )
    gcderive.add_argument("--out-bed", required=True, help="Output BED file for empirical GC-extreme regions")
    gcderive.add_argument("--report", help="Markdown GC derivation report path")
    gcderive.add_argument("--json", help="JSON GC derivation report path")
    gcderive.add_argument("--stdout", action="store_true", help="Print Markdown GC derivation report to stdout")

    seqderive = subparsers.add_parser(
        "seqderive",
        help="Derive empirical homopolymer and low-complexity difficult regions from FASTA",
    )
    seqderive.add_argument("--fasta", required=True, help="Reference or targeted FASTA file")
    seqderive.add_argument("--regions", help="Optional BED regions to scan instead of all FASTA contigs")
    seqderive.add_argument("--window-size", type=int, default=80, help="Sliding window size in bp")
    seqderive.add_argument("--step-size", type=int, default=40, help="Sliding window step size in bp")
    seqderive.add_argument("--min-homopolymer", type=int, default=8, help="Minimum homopolymer run to flag")
    seqderive.add_argument(
        "--low-complexity-fraction",
        type=float,
        default=0.85,
        help="Dominant-base fraction threshold for low-complexity windows",
    )
    seqderive.add_argument(
        "--min-callable-fraction",
        type=float,
        default=0.8,
        help="Minimum fraction of A/C/G/T bases required inside a window",
    )
    seqderive.add_argument("--out-bed", required=True, help="Output BED file for empirical sequence-risk regions")
    seqderive.add_argument("--report", help="Markdown sequence-risk derivation report path")
    seqderive.add_argument("--json", help="JSON sequence-risk derivation report path")
    seqderive.add_argument("--stdout", action="store_true", help="Print Markdown sequence-risk report to stdout")

    return parser


def load_tracks(args: argparse.Namespace) -> tuple[list, list, dict[str, object] | None, str]:
    if args.track_manifest:
        bundle = load_track_manifest(args.track_manifest)
        critical, difficult = load_tracks_from_bundle(bundle)
        genome_build = args.genome_build if args.genome_build != "unknown" else bundle.genome_build
        return critical, difficult, bundle.metadata, genome_build
    critical = parse_bed(args.critical_regions, source="critical") if args.critical_regions else []
    difficult = parse_bed(args.difficult_regions, source="difficult") if args.difficult_regions else []
    return critical, difficult, None, args.genome_build


def run_audit(args: argparse.Namespace) -> int:
    panel = parse_bed(args.panel, source="panel")
    critical, difficult, track_metadata, genome_build = load_tracks(args)
    variants = parse_vcf(args.variants) if args.variants else []

    report = audit_panel(
        panel=panel,
        critical_regions=critical,
        difficult_regions=difficult,
        variants=variants,
        genome_build=genome_build,
        track_metadata=track_metadata,
        panel_name=args.panel_name,
    )

    if args.out:
        write_markdown_report(report, Path(args.out))
    if args.json:
        write_json_report(report, Path(args.json))
    if args.html:
        write_html_report(report, Path(args.html))
    if args.stdout or not (args.out or args.json or args.html):
        print(render_markdown(report))
    return 0


def run_compare(args: argparse.Namespace) -> int:
    panels = load_panel_manifest(args.panels)
    critical, difficult, track_metadata, genome_build = load_tracks(args)
    report = compare_panels(
        panels=panels,
        critical_regions=critical,
        difficult_regions=difficult,
        genome_build=genome_build,
        track_metadata=track_metadata,
    )
    if args.out:
        write_compare_markdown(report, Path(args.out))
    if args.json:
        write_compare_json(report, Path(args.json))
    if args.tsv:
        write_compare_tsv(report, Path(args.tsv))
    if args.stdout or not (args.out or args.json or args.tsv):
        print(render_compare_markdown(report))
    return 0


def run_gcderive(args: argparse.Namespace) -> int:
    sequences = parse_fasta(args.fasta)
    regions = parse_bed(args.regions, source="gcderive_regions") if args.regions else None
    windows = scan_gc_windows(
        sequences=sequences,
        regions=regions,
        window_size=args.window_size,
        step_size=args.step_size,
        low_gc=args.low_gc,
        high_gc=args.high_gc,
        min_callable_fraction=args.min_callable_fraction,
    )
    intervals = merge_gc_windows(windows)
    write_gc_bed(intervals, Path(args.out_bed))
    report = build_gc_report(
        fasta_path=args.fasta,
        windows=windows,
        intervals=intervals,
        window_size=args.window_size,
        step_size=args.step_size,
        low_gc=args.low_gc,
        high_gc=args.high_gc,
        regions=regions,
    )
    if args.report:
        write_gc_report_markdown(report, Path(args.report))
    if args.json:
        write_gc_report_json(report, Path(args.json))
    if args.stdout or not (args.report or args.json):
        print(render_gc_report_markdown(report))
    return 0


def run_seqderive(args: argparse.Namespace) -> int:
    sequences = parse_fasta(args.fasta)
    regions = parse_bed(args.regions, source="seqderive_regions") if args.regions else None
    windows = scan_sequence_risk_windows(
        sequences=sequences,
        regions=regions,
        window_size=args.window_size,
        step_size=args.step_size,
        min_homopolymer=args.min_homopolymer,
        low_complexity_fraction=args.low_complexity_fraction,
        min_callable_fraction=args.min_callable_fraction,
    )
    intervals = merge_sequence_risk_windows(windows)
    write_sequence_risk_bed(intervals, Path(args.out_bed))
    report = build_sequence_risk_report(
        fasta_path=args.fasta,
        windows=windows,
        intervals=intervals,
        window_size=args.window_size,
        step_size=args.step_size,
        min_homopolymer=args.min_homopolymer,
        low_complexity_fraction=args.low_complexity_fraction,
        regions=regions,
    )
    if args.report:
        write_sequence_risk_report_markdown(report, Path(args.report))
    if args.json:
        write_sequence_risk_report_json(report, Path(args.json))
    if args.stdout or not (args.report or args.json):
        print(render_sequence_risk_report_markdown(report))
    return 0


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        if args.command == "audit":
            return run_audit(args)
        if args.command == "compare":
            return run_compare(args)
        if args.command == "gcderive":
            return run_gcderive(args)
        if args.command == "seqderive":
            return run_seqderive(args)
    except AuditorError as exc:
        print(f"panel-bias-auditor: error: {exc}", file=sys.stderr)
        return 1
    parser.error(f"Unknown command: {args.command}")
    return 2
