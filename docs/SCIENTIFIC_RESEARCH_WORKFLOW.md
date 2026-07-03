# Scientific Research Workflow

This project now has two layers:

- audit/report commands;
- research scripts for deriving and testing scientific track ideas before they
  are used in broader panel-design analyses.

The current scientific focus is empirical technical-risk annotation from FASTA
sequence.

## Empirical Technical-Risk Tracks

Use:

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

The script produces:

- GC-extreme BED and reports.
- Homopolymer/low-complexity BED and reports.
- Combined empirical technical-risk BED.
- Research summary Markdown/JSON.

Validate the outputs:

```bash
PYTHONPATH=src python research/validate_empirical_tracks.py \
  --output-dir research_outputs/demo_empirical_tracks \
  --prefix demo \
  --out research_outputs/demo_empirical_tracks/demo_validation_report.md \
  --json research_outputs/demo_empirical_tracks/demo_validation_report.json
```

Analyze assay-performance enrichment:

```bash
PYTHONPATH=src python research/assay_performance_enrichment.py \
  --performance research_inputs/demo_assay_performance.tsv \
  --risk-bed research_outputs/demo_empirical_tracks/demo_combined_technical_risk.bed \
  --output-dir research_outputs/demo_assay_performance \
  --prefix demo \
  --min-depth 100 \
  --min-callable-fraction 0.9
```

Convert public benchmark callability into assay-performance evidence:

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

Demo outputs are in:

```text
research_outputs/demo_empirical_tracks/
```

## Biological Meaning

GC extremes can affect capture efficiency, amplification balance, coverage
uniformity and variant-calling confidence.

Homopolymers and low-complexity sequence can create alignment, indel-calling
and sequencing-error problems, especially in short-read assays and low-VAF
settings.

These tracks should be interpreted as assay-validation prompts. They do not
prove that a locus is uncallable.

## Next Research Direction

For real oncology use, run the research script on:

- GRCh38 FASTA with oncology regions supplied as BED; or
- target-extracted FASTA from an assay design.

Then compare empirical technical-risk density against:

- existing assay targets;
- failed/low-depth regions;
- variant-level false positives/false negatives;
- orthogonal validation calls.

The next useful research scripts are:

- empirical mappability import and summarization;
- segmental duplication/homology import;
- fusion breakpoint interval curation;
- variant-level assay failure enrichment.

## Current Validation Status

The local scientific workflow is validated for toy/demo data:

- unit tests cover FASTA parsing, GC scanning, sequence-risk scanning, CLIs and
  research-output validation;
- demo empirical outputs pass `validate_empirical_tracks.py`;
- demo assay-performance enrichment runs against empirical technical-risk tracks;
- demo public benchmark callability conversion runs against empirical
  technical-risk tracks;
- generated BED, Markdown and JSON outputs are internally consistent.

It is not yet scientifically validated on real oncology assay data. That needs
real GRCh38 or target-extracted FASTA plus observed assay performance data.
