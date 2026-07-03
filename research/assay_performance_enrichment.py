from __future__ import annotations

import argparse
from pathlib import Path

from panel_bias_auditor.performance import (
    analyze_performance_enrichment,
    load_risk_bed,
    parse_performance_table,
    write_performance_report,
)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Analyze whether empirical technical-risk tracks are enriched for observed assay-performance failures."
    )
    parser.add_argument("--performance", required=True, help="Assay performance CSV/TSV")
    parser.add_argument("--risk-bed", required=True, help="Empirical technical-risk BED")
    parser.add_argument("--output-dir", required=True, help="Output directory")
    parser.add_argument("--prefix", default="assay_performance", help="Output prefix")
    parser.add_argument("--min-depth", type=float, help="Optional failure threshold for mean_depth")
    parser.add_argument(
        "--min-callable-fraction",
        type=float,
        help="Optional failure threshold for callable_fraction",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    performance_rows = parse_performance_table(args.performance)
    risk_intervals = load_risk_bed(args.risk_bed)
    report = analyze_performance_enrichment(
        performance_rows,
        risk_intervals,
        min_depth=args.min_depth,
        min_callable_fraction=args.min_callable_fraction,
    )
    markdown = output_dir / f"{args.prefix}_assay_performance_enrichment.md"
    json_path = output_dir / f"{args.prefix}_assay_performance_enrichment.json"
    write_performance_report(report, markdown, json_path)
    print(markdown.read_text(encoding="utf-8"))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
