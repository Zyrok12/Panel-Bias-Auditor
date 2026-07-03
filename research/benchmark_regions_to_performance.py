from __future__ import annotations

import argparse
import json
from pathlib import Path

from panel_bias_auditor.parsers import parse_bed
from panel_bias_auditor.performance import (
    analyze_performance_enrichment,
    load_risk_bed,
    performance_rows_from_callability,
    render_callability_conversion_markdown,
    summarize_callability_rows,
    write_performance_report,
    write_performance_table,
)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Convert public benchmark/callability BED files into the assay-performance "
            "schema used by the technical-risk enrichment workflow."
        )
    )
    parser.add_argument("--evaluated-bed", required=True, help="BED intervals to evaluate, such as panel targets")
    parser.add_argument("--callable-bed", required=True, help="High-confidence/callable benchmark BED")
    parser.add_argument("--output-dir", required=True, help="Output directory")
    parser.add_argument("--prefix", default="public_benchmark", help="Output prefix")
    parser.add_argument("--source-label", default="public_benchmark", help="Label shown in the conversion report")
    parser.add_argument(
        "--min-callable-fraction",
        type=float,
        default=0.95,
        help="Intervals below this callable fraction are marked low_callability",
    )
    parser.add_argument(
        "--risk-bed",
        help="Optional empirical technical-risk BED for immediate enrichment analysis",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    evaluated_regions = parse_bed(args.evaluated_bed, source="evaluated_regions")
    callable_regions = parse_bed(args.callable_bed, source="callable_benchmark")
    rows = performance_rows_from_callability(
        evaluated_regions,
        callable_regions,
        min_callable_fraction=args.min_callable_fraction,
    )

    performance_tsv = output_dir / f"{args.prefix}_public_benchmark_performance.tsv"
    conversion_md = output_dir / f"{args.prefix}_callability_conversion.md"
    conversion_json = output_dir / f"{args.prefix}_callability_conversion.json"

    write_performance_table(rows, performance_tsv)
    conversion_md.write_text(
        render_callability_conversion_markdown(rows, args.source_label, args.min_callable_fraction),
        encoding="utf-8",
    )
    conversion_json.write_text(
        json.dumps(
            {
                "source_label": args.source_label,
                "min_callable_fraction": args.min_callable_fraction,
                "summary": summarize_callability_rows(rows),
                "performance_table": str(performance_tsv),
            },
            indent=2,
            sort_keys=True,
        ),
        encoding="utf-8",
    )

    if args.risk_bed:
        risk_intervals = load_risk_bed(args.risk_bed)
        enrichment = analyze_performance_enrichment(rows, risk_intervals)
        write_performance_report(
            enrichment,
            output_dir / f"{args.prefix}_public_benchmark_enrichment.md",
            output_dir / f"{args.prefix}_public_benchmark_enrichment.json",
        )

    print(conversion_md.read_text(encoding="utf-8"))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
