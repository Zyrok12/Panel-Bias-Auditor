# Assay Performance Validation

The scientific loop is:

1. Derive empirical technical-risk tracks from sequence.
2. Ingest observed assay-performance windows.
3. Test whether technical-risk regions are enriched for assay failures.

This is the step that turns a plausible technical-risk annotation into a
testable scientific hypothesis.

## Performance Table Schema

Use CSV or TSV. Required columns:

| Column | Meaning |
|---|---|
| `chrom` | Chromosome or contig name. |
| `start` | Zero-based interval start. |
| `end` | Half-open interval end. |
| `name` | Stable assay-performance window name. |

Optional columns:

| Column | Meaning |
|---|---|
| `mean_depth` | Mean observed depth in the interval. |
| `callable_fraction` | Fraction of bases passing the assay callability rule. |
| `status` | `pass`, `low_coverage`, `low_callability`, `uncallable`, etc. |
| `false_negative_count` | Number of known positives missed in the interval. |
| `false_positive_count` | Number of false positives observed in the interval. |
| `total_variant_count` | Number of truth/test variants assessed. |

Rows are classified as failures if:

- `status` is a failure-like label;
- `mean_depth` is below `--min-depth`;
- `callable_fraction` is below `--min-callable-fraction`;
- `false_negative_count` or `false_positive_count` is greater than zero.

## Demo Assay Performance Data

Demo input:

```text
research_inputs/demo_assay_performance.tsv
```

Run:

```bash
PYTHONPATH=src python research/assay_performance_enrichment.py \
  --performance research_inputs/demo_assay_performance.tsv \
  --risk-bed research_outputs/demo_empirical_tracks/demo_combined_technical_risk.bed \
  --output-dir research_outputs/demo_assay_performance \
  --prefix demo \
  --min-depth 100 \
  --min-callable-fraction 0.9
```

Outputs:

- `demo_assay_performance_enrichment.md`
- `demo_assay_performance_enrichment.json`

## Interpreting The Enrichment

The analysis computes a base-level 2x2 table:

| | Failure bases | Passing bases |
|---|---:|---:|
| Technical-risk bases | risk/fail | risk/pass |
| Non-risk bases | non-risk/fail | non-risk/pass |

It reports failure rates and a Haldane-corrected odds ratio.

This is exploratory. Bases within genomic intervals are not independent, so the
odds ratio should be interpreted as an effect-size signal, not a clinical proof.

## Real-Data Sources

For real scientific validation, use assay-specific data whenever possible:

- per-target or per-window depth;
- callable-base fraction;
- no-call intervals;
- false negatives and false positives against truth variants;
- orthogonal validation results.

Public benchmark resources can help bootstrap this:

- NIST Genome in a Bottle provides benchmark calls, high-confidence regions,
  sequencing data and difficult-region stratifications for benchmarking variant
  calling pipelines: https://www.nist.gov/programs-projects/genome-bottle
- GIAB release data are hosted at:
  https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/
- Cancer Genome in a Bottle provides public tumor-normal resources and draft
  somatic benchmarks for cancer sequencing validation:
  https://www.nist.gov/programs-projects/cancer-genome-bottle

GIAB is not a substitute for a lab's assay performance table. It is a reference
benchmark source you can use to construct truth/callability labels or evaluate
variant calling around difficult sequence contexts.

## Public Benchmark Conversion

Use `research/benchmark_regions_to_performance.py` when the validation source is
a public high-confidence or callable BED rather than an assay-performance table.

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

This emits an assay-performance-compatible TSV and, if `--risk-bed` is supplied,
an immediate technical-risk enrichment report.

See `docs/PUBLIC_VALIDATION_DATASETS.md` for dataset sources and interpretation
limits.

## What Is Finished

- Performance table parser.
- Failure classification by status, depth, callability and error counts.
- Base-level enrichment against empirical technical-risk tracks.
- Public benchmark/callability BED conversion into the performance schema.
- Markdown/JSON reporting.
- Demo assay-performance dataset.
- Unit tests for parser, enrichment and research script.

## What Still Requires Real Biology

- Real GRCh38 or target-extracted FASTA.
- Real assay windows and coverage metrics.
- Truth variants or orthogonal validation calls.
- Analysis across multiple samples/runs to estimate reproducibility.
