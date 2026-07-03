# Panel Bias Auditor Lite

`panel-bias-auditor-lite` is a dependency-light Python CLI and research toolkit
for auditing genomic panel designs against user-provided biological risk tracks.

It is built for research-use assay design review: where does a panel miss
important biology, and where might apparently covered regions still be hard to
call because of genomic context?

## What It Does

Given a panel BED file, the auditor can:

- Calculate target footprint and merged covered bases.
- Measure how well the panel covers a critical-region BED track.
- Flag overlap with difficult-region tracks such as repeats, low mappability, high-GC proxies, homologous genes, pseudogene-like regions, or pangenome/reference-bias tracks.
- Optionally annotate VCF variants as inside/outside the panel and inside/outside critical or difficult regions.
- Generate practical recommendations about missing critical coverage, difficult-region burden, provenance, and supplied variant checks.
- Produce a Markdown report and machine-readable JSON.
- Derive empirical GC-extreme, homopolymer and low-complexity technical-risk
  tracks from FASTA sequence.
- Test whether technical-risk tracks are enriched in assay-performance failures
  or low-callability public benchmark regions.

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

## Derive Empirical Technical-Risk Tracks

The repo includes public-safe research commands for generating technical-risk
tracks directly from sequence.

GC extremes:

```bash
PYTHONPATH=src python3 -m panel_bias_auditor gcderive \
  --fasta examples/demo_reference.fa \
  --regions examples/demo_gc_scan_regions.bed \
  --window-size 50 \
  --step-size 50 \
  --low-gc 0.2 \
  --high-gc 0.8 \
  --out-bed examples/demo_empirical_gc_extremes.bed \
  --report examples/demo_empirical_gc_report.md \
  --json examples/demo_empirical_gc_report.json
```

Combined GC plus sequence-complexity research workflow:

```bash
PYTHONPATH=src python research/empirical_technical_tracks.py \
  --fasta examples/demo_reference.fa \
  --regions examples/demo_gc_scan_regions.bed \
  --output-dir research_outputs/demo_empirical_tracks \
  --prefix demo \
  --gc-window 50 \
  --gc-step 50 \
  --low-gc 0.2 \
  --high-gc 0.8 \
  --seq-window 40 \
  --seq-step 40 \
  --min-homopolymer 10 \
  --low-complexity-fraction 0.8
```

## Validate Against Performance Or Benchmark Data

Assay-performance enrichment:

```bash
PYTHONPATH=src python research/assay_performance_enrichment.py \
  --performance research_inputs/demo_assay_performance.tsv \
  --risk-bed research_outputs/demo_empirical_tracks/demo_combined_technical_risk.bed \
  --output-dir research_outputs/demo_assay_performance \
  --prefix demo \
  --min-depth 100 \
  --min-callable-fraction 0.9
```

Public benchmark/callability conversion:

```bash
PYTHONPATH=src python research/benchmark_regions_to_performance.py \
  --evaluated-bed research_inputs/demo_public_benchmark_regions.bed \
  --callable-bed research_inputs/demo_public_callable_regions.bed \
  --risk-bed research_outputs/demo_empirical_tracks/demo_combined_technical_risk.bed \
  --output-dir research_outputs/demo_public_benchmark \
  --prefix demo \
  --source-label demo_public_benchmark \
  --min-callable-fraction 0.9
```

The public dataset manifest is in
`research_inputs/public_validation_sources.tsv`. It points to resources such as
NIST Genome in a Bottle, GIAB stratification BEDs and Cancer Genome in a Bottle.

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

### Research Outputs

Small demo outputs are included under `research_outputs/` so expected scientific
artifacts can be inspected without downloading large public datasets.

Generated research artifacts include:

- empirical GC and sequence-risk BED files;
- validation summaries for generated tracks;
- assay-performance enrichment reports;
- public benchmark callability conversion reports.

## Recommendation Layer

The report includes a small rule-based recommendation layer. This is deliberately transparent:

- Missing or partial critical regions produce assay-design recommendations.
- Difficult-region burden produces validation-priority recommendations.
- Supplied variants outside the panel are flagged as potential scope or design mismatches.
- Raw BED use without a track manifest produces a provenance recommendation.

This is not clinical decision support. It is a research-use assay-design explanation layer.

## Important Caveat

This public version is not a medical device and is not intended for clinical diagnosis or treatment selection. It is a research-use and assay-design aid.

The demo tracks in `examples/` and `research_outputs/` are deliberately small
and illustrative. Do not use them as production truth sets.

Public benchmark callability is useful validation evidence, but it is not the
same thing as a hospital or commercial lab's assay-performance table. Clinical
or regulated use requires assay-specific validation against appropriate truth
sets, coverage metrics and orthogonal confirmation.

## Tests

```bash
PYTHONPATH=src python -m unittest discover -s tests
```

Current tests cover parsers, panel audits, track bundles, empirical technical
track derivation, public benchmark conversion and assay-performance enrichment.

## Possible Extensions

- Additional curated oncology, germline, or rare disease demonstration tracks.
- Pangenome-aware callability tracks.
- Raw BAM/CRAM coverage extraction into the public performance schema.
- PDF report export.
- Larger parser and CLI test coverage.
- Example notebooks or visualizations.
