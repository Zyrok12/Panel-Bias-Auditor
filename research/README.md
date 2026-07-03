# Research Scripts

This folder is for scientific exploration scripts that derive, validate and
stress-test technical-risk tracks before they are used in broader analyses.

## `empirical_technical_tracks.py`

Builds empirical technical-risk tracks from FASTA sequence:

- GC extremes.
- Homopolymer risk.
- Low-complexity risk.
- Combined difficult-region BED.

Example:

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

These outputs are research annotations. They identify regions that deserve
assay-validation review; they do not prove clinical uncallability.

## `validate_empirical_tracks.py`

Validates the output directory from `empirical_technical_tracks.py`.

It checks:

- expected files exist;
- BED files parse;
- required BED attributes are present;
- JSON interval counts match BED outputs;
- combined BED counts match component BEDs;
- combined merged bases match the research summary.

Example:

```bash
PYTHONPATH=src python research/validate_empirical_tracks.py \
  --output-dir research_outputs/demo_empirical_tracks \
  --prefix demo \
  --out research_outputs/demo_empirical_tracks/demo_validation_report.md \
  --json research_outputs/demo_empirical_tracks/demo_validation_report.json
```

## `assay_performance_enrichment.py`

Tests whether empirical technical-risk tracks overlap observed assay failures.

Example:

```bash
PYTHONPATH=src python research/assay_performance_enrichment.py \
  --performance research_inputs/demo_assay_performance.tsv \
  --risk-bed research_outputs/demo_empirical_tracks/demo_combined_technical_risk.bed \
  --output-dir research_outputs/demo_assay_performance \
  --prefix demo \
  --min-depth 100 \
  --min-callable-fraction 0.9
```

The input table can contain depth, callable fraction, failure status and
false-positive/false-negative counts.

## `benchmark_regions_to_performance.py`

Converts a public benchmark/high-confidence BED into the same performance table
schema used by `assay_performance_enrichment.py`.

Example:

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

This is the current bridge from public GIAB-style datasets into our scientific
validation workflow. It treats poor benchmark callability as exploratory
evidence, not as clinical assay failure.
