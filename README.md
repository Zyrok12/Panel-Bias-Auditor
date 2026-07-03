# Panel Bias Auditor Lite

`panel-bias-auditor-lite` is a dependency-light Python CLI for auditing genomic panel designs against user-provided biological risk tracks.

## What It Does

Given a panel BED file, the auditor can:

- Calculate target footprint and merged covered bases.
- Measure how well the panel covers a critical-region BED track.
- Flag overlap with difficult-region tracks such as repeats, low mappability, high-GC proxies, homologous genes, pseudogene-like regions, or pangenome/reference-bias tracks.
- Optionally annotate VCF variants as inside/outside the panel and inside/outside critical or difficult regions.
- Generate practical recommendations about missing critical coverage, difficult-region burden, provenance, and supplied variant checks.
- Produce a Markdown report and machine-readable JSON.

The tool is intentionally track-driven. The repo includes small demo tracks only; real analyses should use curated, versioned tracks selected for the assay and biological question.

## Quick Start

From this directory:

```bash
PYTHONPATH=src python3 -m panel_bias_auditor audit \
  --panel examples/demo_oncology_panel.bed \
  --track-manifest examples/demo_track_manifest.json \
  --variants examples/demo_variants.vcf \
  --out examples/demo_report.md \
  --json examples/demo_report.json \
  --html examples/demo_report.html \
  --genome-build GRCh38-demo
```

On Windows through WSL from this repository:

```powershell
wsl.exe -e bash -lc "cd /path/to/panel-bias-auditor-lite && PYTHONPATH=src python3 -m panel_bias_auditor audit --panel examples/demo_oncology_panel.bed --track-manifest examples/demo_track_manifest.json --variants examples/demo_variants.vcf --out examples/demo_report.md --json examples/demo_report.json --html examples/demo_report.html"
```

## Compare Multiple Panels

```bash
PYTHONPATH=src python3 -m panel_bias_auditor compare \
  --panels examples/demo_panels.tsv \
  --track-manifest examples/demo_track_manifest.json \
  --out examples/demo_comparison.md \
  --json examples/demo_comparison.json \
  --tsv examples/demo_comparison.tsv
```

## Inputs

### Panel BED

Required. Standard BED-like file:

```text
chrom  start  end  name
```

Coordinates are interpreted as zero-based half-open intervals.

### Critical Regions BED

Optional but strongly recommended. Regions that should be confidently covered for the panel use case.

Examples:

- Cancer hotspots.
- Full coding regions of high-value genes.
- Fusion-relevant introns.
- ctDNA or tumor-informed marker candidate regions.
- Clinically important founder or recurrent variants.

### Difficult Regions BED

Optional. Regions where "covered" does not necessarily mean "confidently callable."

Examples:

- Low mappability.
- Segmental duplications.
- Repeats.
- GC extremes.
- Homologous/pseudogene-prone regions.
- Pangenome/reference-bias regions.

### VCF

Optional. Used to check whether observed or test variants fall inside the panel and/or overlap known risk tracks.

## Outputs

### Markdown Report

Human-readable report for technical review.

### JSON Report

Machine-readable output for downstream analysis, dashboards, APIs, or batch workflows.

### HTML Report

Standalone browser-viewable report for demos and technical review.

## Recommendation Layer

The report includes a small rule-based recommendation layer. This is deliberately transparent:

- Missing or partial critical regions produce assay-design recommendations.
- Difficult-region burden produces validation-priority recommendations.
- Supplied variants outside the panel are flagged as potential scope or design mismatches.
- Raw BED use without a track manifest produces a provenance recommendation.

This is not clinical decision support. It is a research-use assay-design explanation layer.

## Important Caveat

This public version is not a medical device and is not intended for clinical diagnosis or treatment selection. It is a research-use and assay-design aid.

The demo tracks in `examples/` are deliberately small and illustrative. Do not use them as production truth sets.

## Possible Extensions

- Additional curated oncology, germline, or rare disease demonstration tracks.
- Pangenome-aware callability tracks.
- PDF report export.
- Larger parser and CLI test coverage.
- Example notebooks or visualizations.
