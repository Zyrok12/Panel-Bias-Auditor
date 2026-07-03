from __future__ import annotations

from collections import defaultdict
from dataclasses import asdict
from datetime import datetime, timezone
from html import escape
from pathlib import Path
from typing import Iterable

from .models import Interval, RegionCoverage, Variant


def index_intervals(intervals: Iterable[Interval]) -> dict[str, list[Interval]]:
    indexed: dict[str, list[Interval]] = defaultdict(list)
    for interval in intervals:
        indexed[interval.norm_chrom].append(interval)
    for chrom in indexed:
        indexed[chrom].sort(key=lambda item: (item.start, item.end, item.name))
    return dict(indexed)


def merge_spans(spans: Iterable[tuple[int, int]]) -> list[tuple[int, int]]:
    sorted_spans = sorted(spans)
    if not sorted_spans:
        return []
    merged = [sorted_spans[0]]
    for start, end in sorted_spans[1:]:
        last_start, last_end = merged[-1]
        if start <= last_end:
            merged[-1] = (last_start, max(last_end, end))
        else:
            merged.append((start, end))
    return merged


def total_span_bases(spans: Iterable[tuple[int, int]]) -> int:
    return sum(end - start for start, end in merge_spans(spans))


def merged_interval_bases(intervals: Iterable[Interval]) -> int:
    by_chrom: dict[str, list[tuple[int, int]]] = defaultdict(list)
    for interval in intervals:
        by_chrom[interval.norm_chrom].append((interval.start, interval.end))
    return sum(total_span_bases(spans) for spans in by_chrom.values())


def overlap_bases(target: Interval, query_intervals: Iterable[Interval]) -> int:
    spans: list[tuple[int, int]] = []
    for query in query_intervals:
        span = target.overlap_span(query)
        if span is not None:
            spans.append(span)
    return total_span_bases(spans)


def coverage_by_region(panel: list[Interval], critical_regions: list[Interval]) -> list[RegionCoverage]:
    panel_by_chrom = index_intervals(panel)
    result: list[RegionCoverage] = []
    for region in critical_regions:
        covered = overlap_bases(region, panel_by_chrom.get(region.norm_chrom, []))
        result.append(RegionCoverage(region=region, covered_bases=covered))
    return result


def difficult_overlap(panel: list[Interval], difficult_regions: list[Interval]) -> dict[str, object]:
    difficult_by_chrom = index_intervals(difficult_regions)
    total_spans_by_chrom: dict[str, list[tuple[int, int]]] = defaultdict(list)
    by_category: dict[str, int] = defaultdict(int)
    by_region: list[dict[str, object]] = []

    for panel_interval in panel:
        overlapping_spans: list[tuple[int, int]] = []
        for difficult in difficult_by_chrom.get(panel_interval.norm_chrom, []):
            span = panel_interval.overlap_span(difficult)
            if span is None:
                continue
            category = difficult.attrs.get("type") or difficult.attrs.get("category") or difficult.name
            span_bases = span[1] - span[0]
            by_category[category] += span_bases
            overlapping_spans.append(span)
            total_spans_by_chrom[panel_interval.norm_chrom].append(span)
            by_region.append(
                {
                    "panel_region": panel_interval.name,
                    "difficult_region": difficult.name,
                    "category": category,
                    "overlap_bases": span_bases,
                    "chrom": panel_interval.chrom,
                    "start": span[0],
                    "end": span[1],
                }
            )
        # Keep the loop explicit for future per-panel-region risk features.
        _ = overlapping_spans

    total_unique = sum(total_span_bases(spans) for spans in total_spans_by_chrom.values())
    return {
        "unique_overlap_bases": total_unique,
        "by_category": dict(sorted(by_category.items(), key=lambda item: (-item[1], item[0]))),
        "by_region": sorted(by_region, key=lambda item: (-int(item["overlap_bases"]), str(item["panel_region"])))[:50],
    }


def find_overlaps(query: Interval, targets: list[Interval]) -> list[Interval]:
    return [target for target in targets if query.overlaps(target)]


def annotate_variants(
    variants: list[Variant],
    panel: list[Interval],
    critical_regions: list[Interval],
    difficult_regions: list[Interval],
) -> list[dict[str, object]]:
    panel_by_chrom = index_intervals(panel)
    critical_by_chrom = index_intervals(critical_regions)
    difficult_by_chrom = index_intervals(difficult_regions)
    annotated: list[dict[str, object]] = []

    for variant in variants:
        interval = variant.interval
        panel_hits = find_overlaps(interval, panel_by_chrom.get(interval.norm_chrom, []))
        critical_hits = find_overlaps(interval, critical_by_chrom.get(interval.norm_chrom, []))
        difficult_hits = find_overlaps(interval, difficult_by_chrom.get(interval.norm_chrom, []))
        annotated.append(
            {
                "variant": f"{variant.chrom}:{variant.pos}:{variant.ref}>{variant.alt}",
                "id": variant.identifier,
                "inside_panel": bool(panel_hits),
                "panel_regions": [hit.name for hit in panel_hits],
                "critical_regions": [hit.name for hit in critical_hits],
                "difficult_regions": [hit.name for hit in difficult_hits],
                "filter": variant.filt,
                "qual": variant.qual,
            }
        )
    return annotated


def summarize_coverage(region_coverages: list[RegionCoverage]) -> dict[str, object]:
    total_bases = sum(item.region.length for item in region_coverages)
    covered_bases = sum(item.covered_bases for item in region_coverages)
    missing_bases = total_bases - covered_bases
    status_counts: dict[str, int] = defaultdict(int)
    for item in region_coverages:
        status_counts[item.status] += 1
    return {
        "region_count": len(region_coverages),
        "total_bases": total_bases,
        "covered_bases": covered_bases,
        "missing_bases": missing_bases,
        "coverage_fraction": covered_bases / total_bases if total_bases else None,
        "status_counts": dict(status_counts),
    }


def risk_label(score: float) -> str:
    if score <= 20:
        return "low"
    if score <= 50:
        return "moderate"
    if score <= 75:
        return "high"
    return "severe"


def compute_risk_score(
    panel_bases: int,
    critical_summary: dict[str, object],
    hard_unique_overlap: int,
    variant_annotations: list[dict[str, object]],
) -> dict[str, object]:
    coverage_fraction = critical_summary.get("coverage_fraction")
    critical_missing_fraction = 0.0 if coverage_fraction is None else 1.0 - float(coverage_fraction)
    hard_fraction = hard_unique_overlap / panel_bases if panel_bases else 0.0
    if variant_annotations:
        outside = sum(1 for item in variant_annotations if not item["inside_panel"])
        outside_fraction = outside / len(variant_annotations)
    else:
        outside_fraction = 0.0

    hard_component = min(hard_fraction * 5.0, 1.0)
    score = 100.0 * (
        0.65 * critical_missing_fraction
        + 0.25 * hard_component
        + 0.10 * outside_fraction
    )
    score = max(0.0, min(100.0, score))
    return {
        "score": round(score, 2),
        "label": risk_label(score),
        "components": {
            "critical_missing_fraction": round(critical_missing_fraction, 4),
            "difficult_region_fraction": round(hard_fraction, 4),
            "variant_outside_panel_fraction": round(outside_fraction, 4),
        },
    }


def generate_recommendations(
    critical_summary: dict[str, object],
    difficult_summary: dict[str, object],
    variant_annotations: list[dict[str, object]],
    missing_regions: list[dict[str, object]],
    track_metadata: dict[str, object] | None,
) -> list[dict[str, str]]:
    recommendations: list[dict[str, str]] = []
    coverage_fraction = critical_summary.get("coverage_fraction")
    missing_bases = int(critical_summary.get("missing_bases") or 0)

    if coverage_fraction is None:
        recommendations.append(
            {
                "priority": "high",
                "category": "critical-track",
                "finding": "No critical-region track was supplied.",
                "action": "Define a use-case-specific critical-region BED before interpreting assay quality.",
            }
        )
    elif coverage_fraction < 0.95:
        worst = ", ".join(str(item["name"]) for item in missing_regions[:3]) or "none"
        recommendations.append(
            {
                "priority": "high",
                "category": "critical-coverage",
                "finding": f"{missing_bases:,} critical bases are missing or partially covered. Top affected regions: {worst}.",
                "action": "Review whether these regions are biologically required for the panel claim; add probes/amplicons or narrow the assay claim.",
            }
        )

    hard_by_category = difficult_summary.get("by_category") or {}
    if hard_by_category:
        top_category, top_bases = next(iter(hard_by_category.items()))
        recommendations.append(
            {
                "priority": "medium",
                "category": "callability",
                "finding": f"The largest difficult-region burden is {top_category} ({int(top_bases):,} bp overlap).",
                "action": "Treat these regions as validation priorities; confirm sensitivity, specificity, and false-negative behavior with orthogonal or truth-set data.",
            }
        )

    outside_variants = [item for item in variant_annotations if not item["inside_panel"]]
    if outside_variants:
        examples = ", ".join(str(item["variant"]) for item in outside_variants[:3])
        recommendations.append(
            {
                "priority": "medium",
                "category": "variant-check",
                "finding": f"{len(outside_variants)} supplied variants fall outside the panel. Examples: {examples}.",
                "action": "If these variants represent expected positives or customer requirements, revise targets or document them as out-of-scope.",
            }
        )

    if track_metadata is None:
        recommendations.append(
            {
                "priority": "medium",
                "category": "provenance",
                "finding": "Tracks were supplied as raw BED files without a versioned manifest.",
                "action": "Use a track manifest so reports preserve track names, versions, source descriptions, and genome build.",
            }
        )

    if not recommendations:
        recommendations.append(
            {
                "priority": "low",
                "category": "review",
                "finding": "No major design-risk flags were detected against the supplied tracks.",
                "action": "Before any clinical or commercial claim, validate against real truth sets and confirm the track bundle matches the assay indication.",
            }
        )
    return recommendations


def audit_panel(
    panel: list[Interval],
    critical_regions: list[Interval] | None = None,
    difficult_regions: list[Interval] | None = None,
    variants: list[Variant] | None = None,
    genome_build: str = "unknown",
    track_metadata: dict[str, object] | None = None,
    panel_name: str | None = None,
) -> dict[str, object]:
    critical_regions = critical_regions or []
    difficult_regions = difficult_regions or []
    variants = variants or []

    panel_bases = merged_interval_bases(panel)
    region_coverages = coverage_by_region(panel, critical_regions)
    critical_summary = summarize_coverage(region_coverages)
    hard = difficult_overlap(panel, difficult_regions)
    variant_annotations = annotate_variants(variants, panel, critical_regions, difficult_regions)
    risk = compute_risk_score(
        panel_bases=panel_bases,
        critical_summary=critical_summary,
        hard_unique_overlap=int(hard["unique_overlap_bases"]),
        variant_annotations=variant_annotations,
    )

    missing_regions = [
        {
            "name": item.region.name,
            "chrom": item.region.chrom,
            "start": item.region.start,
            "end": item.region.end,
            "length": item.region.length,
            "covered_bases": item.covered_bases,
            "missing_bases": item.missing_bases,
            "coverage_fraction": round(item.coverage_fraction, 4),
            "status": item.status,
        }
        for item in region_coverages
        if item.status != "covered"
    ]
    missing_regions.sort(key=lambda item: (-int(item["missing_bases"]), str(item["name"])))
    recommendations = generate_recommendations(
        critical_summary=critical_summary,
        difficult_summary=hard,
        variant_annotations=variant_annotations,
        missing_regions=missing_regions,
        track_metadata=track_metadata,
    )

    return {
        "metadata": {
            "generated_at": datetime.now(timezone.utc).isoformat(),
            "genome_build": genome_build,
            "tool": "panel-bias-auditor-lite",
            "track_bundle": track_metadata,
        },
        "panel": {
            "name": panel_name,
            "interval_count": len(panel),
            "merged_target_bases": panel_bases,
        },
        "critical_regions": critical_summary,
        "difficult_regions": hard,
        "variants": {
            "count": len(variant_annotations),
            "outside_panel_count": sum(1 for item in variant_annotations if not item["inside_panel"]),
            "annotations": variant_annotations,
        },
        "risk": risk,
        "recommendations": recommendations,
        "missing_or_partial_critical_regions": missing_regions[:50],
    }


def write_json_report(report: dict[str, object], path: str | Path) -> None:
    import json

    Path(path).write_text(json.dumps(report, indent=2, sort_keys=True), encoding="utf-8")


def pct(value: float | None) -> str:
    if value is None:
        return "n/a"
    return f"{value * 100:.1f}%"


def render_markdown(report: dict[str, object]) -> str:
    metadata = report["metadata"]
    panel = report["panel"]
    critical = report["critical_regions"]
    difficult = report["difficult_regions"]
    variants = report["variants"]
    risk = report["risk"]

    lines: list[str] = []
    lines.append("# Panel Bias Audit Report")
    lines.append("")
    lines.append(f"- Generated: `{metadata['generated_at']}`")
    lines.append(f"- Genome build: `{metadata['genome_build']}`")
    lines.append(f"- Tool: `{metadata['tool']}`")
    if metadata.get("track_bundle"):
        bundle = metadata["track_bundle"]
        lines.append(f"- Track bundle: `{bundle.get('name', 'unknown')}` version `{bundle.get('version', '')}`")
    lines.append("")
    lines.append("## Executive Summary")
    lines.append("")
    lines.append(f"- Risk score: **{risk['score']} / 100** ({risk['label']})")
    lines.append(f"- Panel intervals: **{panel['interval_count']}**")
    if panel.get("name"):
        lines.append(f"- Panel name: **{panel['name']}**")
    lines.append(f"- Merged target footprint: **{panel['merged_target_bases']:,} bp**")
    lines.append(f"- Critical-region coverage: **{pct(critical['coverage_fraction'])}**")
    lines.append(f"- Missing critical bases: **{critical['missing_bases']:,} bp**")
    lines.append(f"- Difficult-region overlap: **{difficult['unique_overlap_bases']:,} bp**")
    lines.append(f"- Variants outside panel: **{variants['outside_panel_count']} / {variants['count']}**")
    lines.append("")
    lines.append("## Risk Components")
    lines.append("")
    lines.append("| Component | Value |")
    lines.append("|---|---:|")
    for key, value in risk["components"].items():
        lines.append(f"| {key} | {value} |")
    lines.append("")
    lines.append("## Recommendations")
    lines.append("")
    lines.append("| Priority | Category | Finding | Suggested action |")
    lines.append("|---|---|---|---|")
    for item in report["recommendations"]:
        lines.append(f"| {item['priority']} | {item['category']} | {item['finding']} | {item['action']} |")
    lines.append("")
    lines.append("## Critical Regions")
    lines.append("")
    lines.append(f"- Critical regions assessed: **{critical['region_count']}**")
    lines.append(f"- Covered bases: **{critical['covered_bases']:,} / {critical['total_bases']:,} bp**")
    lines.append(f"- Status counts: `{critical['status_counts']}`")
    lines.append("")

    missing = report["missing_or_partial_critical_regions"]
    if missing:
        lines.append("### Missing Or Partial Critical Regions")
        lines.append("")
        lines.append("| Region | Locus | Status | Covered | Missing |")
        lines.append("|---|---|---|---:|---:|")
        for item in missing[:25]:
            locus = f"{item['chrom']}:{item['start']}-{item['end']}"
            lines.append(
                f"| {item['name']} | {locus} | {item['status']} | {item['covered_bases']:,} | {item['missing_bases']:,} |"
            )
        lines.append("")
    else:
        lines.append("All supplied critical regions were fully covered.")
        lines.append("")

    lines.append("## Difficult Region Overlap")
    lines.append("")
    if difficult["by_category"]:
        lines.append("| Category | Overlap bases |")
        lines.append("|---|---:|")
        for category, bases in difficult["by_category"].items():
            lines.append(f"| {category} | {bases:,} |")
        lines.append("")
    else:
        lines.append("No difficult-region overlaps were found in the supplied track.")
        lines.append("")

    if variants["count"]:
        lines.append("## Variant Checks")
        lines.append("")
        lines.append("| Variant | In panel | Panel regions | Critical regions | Difficult regions |")
        lines.append("|---|---|---|---|---|")
        for item in variants["annotations"][:50]:
            lines.append(
                "| {variant} | {inside} | {panel_regions} | {critical_regions} | {difficult_regions} |".format(
                    variant=item["variant"],
                    inside="yes" if item["inside_panel"] else "no",
                    panel_regions=", ".join(item["panel_regions"]) or "-",
                    critical_regions=", ".join(item["critical_regions"]) or "-",
                    difficult_regions=", ".join(item["difficult_regions"]) or "-",
                )
            )
        lines.append("")

    lines.append("## Interpretation Notes")
    lines.append("")
    lines.append("- This report evaluates assay-design risk, not clinical validity.")
    lines.append("- A region can be technically covered but still hard to call if it overlaps repeats, homologous sequence, GC extremes, or reference-biased loci.")
    lines.append("- The quality of the audit depends on the quality of the critical-region and difficult-region tracks supplied.")
    lines.append("- For commercial or clinical use, tracks must be curated, versioned, and validated against appropriate truth sets.")
    lines.append("")
    return "\n".join(lines)


def write_markdown_report(report: dict[str, object], path: str | Path) -> None:
    Path(path).write_text(render_markdown(report), encoding="utf-8")


def render_html(report: dict[str, object]) -> str:
    """Render a compact standalone HTML report without external dependencies."""
    metadata = report["metadata"]
    panel = report["panel"]
    critical = report["critical_regions"]
    difficult = report["difficult_regions"]
    variants = report["variants"]
    risk = report["risk"]
    risk_class = str(risk["label"])

    missing_rows = []
    for item in report["missing_or_partial_critical_regions"][:30]:
        locus = f"{item['chrom']}:{item['start']}-{item['end']}"
        missing_rows.append(
            "<tr>"
            f"<td>{escape(str(item['name']))}</td>"
            f"<td>{escape(locus)}</td>"
            f"<td>{escape(str(item['status']))}</td>"
            f"<td>{int(item['covered_bases']):,}</td>"
            f"<td>{int(item['missing_bases']):,}</td>"
            "</tr>"
        )

    difficult_rows = []
    for category, bases in difficult["by_category"].items():
        difficult_rows.append(f"<tr><td>{escape(str(category))}</td><td>{int(bases):,}</td></tr>")

    variant_rows = []
    for item in variants["annotations"][:50]:
        variant_rows.append(
            "<tr>"
            f"<td>{escape(str(item['variant']))}</td>"
            f"<td>{'yes' if item['inside_panel'] else 'no'}</td>"
            f"<td>{escape(', '.join(item['panel_regions']) or '-')}</td>"
            f"<td>{escape(', '.join(item['critical_regions']) or '-')}</td>"
            f"<td>{escape(', '.join(item['difficult_regions']) or '-')}</td>"
            "</tr>"
        )

    recommendation_rows = []
    for item in report["recommendations"]:
        recommendation_rows.append(
            "<tr>"
            f"<td>{escape(str(item['priority']))}</td>"
            f"<td>{escape(str(item['category']))}</td>"
            f"<td>{escape(str(item['finding']))}</td>"
            f"<td>{escape(str(item['action']))}</td>"
            "</tr>"
        )

    return f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <title>Panel Bias Audit Report</title>
  <style>
    :root {{
      color-scheme: light;
      --ink: #17202a;
      --muted: #5d6d7e;
      --line: #d5dde5;
      --soft: #f5f7fa;
      --low: #247a3d;
      --moderate: #946200;
      --high: #a94700;
      --severe: #9f1d35;
    }}
    body {{
      margin: 0;
      font-family: Arial, Helvetica, sans-serif;
      color: var(--ink);
      background: #ffffff;
      line-height: 1.45;
    }}
    main {{
      max-width: 1080px;
      margin: 0 auto;
      padding: 32px 24px 48px;
    }}
    h1, h2, h3 {{ margin: 0 0 12px; }}
    h2 {{ margin-top: 28px; border-bottom: 1px solid var(--line); padding-bottom: 6px; }}
    .meta, .note {{ color: var(--muted); }}
    .summary {{
      display: grid;
      grid-template-columns: repeat(auto-fit, minmax(180px, 1fr));
      gap: 12px;
      margin: 20px 0;
    }}
    .metric {{
      border: 1px solid var(--line);
      border-radius: 8px;
      padding: 14px;
      background: var(--soft);
    }}
    .metric strong {{ display: block; font-size: 1.35rem; margin-top: 4px; }}
    .risk {{
      color: white;
      background: var(--{risk_class});
      border-radius: 8px;
      padding: 18px;
      margin: 16px 0;
    }}
    table {{
      width: 100%;
      border-collapse: collapse;
      margin: 12px 0 20px;
      font-size: 0.95rem;
    }}
    th, td {{
      border: 1px solid var(--line);
      padding: 8px 10px;
      text-align: left;
      vertical-align: top;
    }}
    th {{ background: var(--soft); }}
    code {{ background: var(--soft); padding: 2px 4px; border-radius: 4px; }}
  </style>
</head>
<body>
<main>
  <h1>Panel Bias Audit Report</h1>
  <p class="meta">Generated: <code>{escape(str(metadata['generated_at']))}</code> | Genome build: <code>{escape(str(metadata['genome_build']))}</code> | Tool: <code>{escape(str(metadata['tool']))}</code></p>

  <section class="risk">
    <h2 style="border:0;margin-top:0;padding:0;">Risk Score: {risk['score']} / 100 ({escape(risk_class)})</h2>
    <p>This score summarizes missing critical-region coverage, difficult-region overlap, and optional variant checks.</p>
  </section>

  <section class="summary">
    <div class="metric">Panel<strong>{escape(str(panel.get('name') or 'single panel'))}</strong></div>
    <div class="metric">Panel intervals<strong>{panel['interval_count']}</strong></div>
    <div class="metric">Merged target footprint<strong>{int(panel['merged_target_bases']):,} bp</strong></div>
    <div class="metric">Critical coverage<strong>{pct(critical['coverage_fraction'])}</strong></div>
    <div class="metric">Missing critical bases<strong>{int(critical['missing_bases']):,} bp</strong></div>
    <div class="metric">Difficult overlap<strong>{int(difficult['unique_overlap_bases']):,} bp</strong></div>
    <div class="metric">Variants outside panel<strong>{variants['outside_panel_count']} / {variants['count']}</strong></div>
  </section>

  <h2>Risk Components</h2>
  <table>
    <tr><th>Component</th><th>Value</th></tr>
    {''.join(f"<tr><td>{escape(str(k))}</td><td>{v}</td></tr>" for k, v in risk['components'].items())}
  </table>

  <h2>Recommendations</h2>
  <table>
    <tr><th>Priority</th><th>Category</th><th>Finding</th><th>Suggested action</th></tr>
    {''.join(recommendation_rows)}
  </table>

  <h2>Missing Or Partial Critical Regions</h2>
  <table>
    <tr><th>Region</th><th>Locus</th><th>Status</th><th>Covered</th><th>Missing</th></tr>
    {''.join(missing_rows) if missing_rows else '<tr><td colspan="5">All supplied critical regions were fully covered.</td></tr>'}
  </table>

  <h2>Difficult Region Overlap</h2>
  <table>
    <tr><th>Category</th><th>Overlap bases</th></tr>
    {''.join(difficult_rows) if difficult_rows else '<tr><td colspan="2">No difficult-region overlaps were found.</td></tr>'}
  </table>

  <h2>Variant Checks</h2>
  <table>
    <tr><th>Variant</th><th>In panel</th><th>Panel regions</th><th>Critical regions</th><th>Difficult regions</th></tr>
    {''.join(variant_rows) if variant_rows else '<tr><td colspan="5">No VCF was supplied.</td></tr>'}
  </table>

  <h2>Interpretation Notes</h2>
  <p class="note">This report evaluates assay-design risk, not clinical validity. A region can be technically covered but still hard to call if it overlaps repeats, homologous sequence, GC extremes, or reference-biased loci. For commercial or clinical use, tracks must be curated, versioned, and validated against appropriate truth sets.</p>
</main>
</body>
</html>
"""


def write_html_report(report: dict[str, object], path: str | Path) -> None:
    Path(path).write_text(render_html(report), encoding="utf-8")


def report_to_plain_dict(report: dict[str, object]) -> dict[str, object]:
    """Hook for future dataclass-heavy reports; currently returns the report unchanged."""
    return asdict(report) if hasattr(report, "__dataclass_fields__") else report
