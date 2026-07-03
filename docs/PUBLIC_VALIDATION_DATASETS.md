# Public Validation Datasets

Yes: we can use public experimental and benchmark datasets for validation.

The safest first use is not to download huge raw BAM/FASTQ files. It is to use
public benchmark/high-confidence BED files as a callability proxy, then test
whether our empirical technical-risk tracks overlap regions with poor benchmark
callability.

## Best Immediate Sources

Machine-readable source manifest:

```text
research_inputs/public_validation_sources.tsv
```

| Source | Use | URL |
|---|---|---|
| NIST Genome in a Bottle | Germline benchmark VCF/BED, high-confidence regions, sequencing data indexes. | https://www.nist.gov/programs-projects/genome-bottle |
| GIAB release FTP | Direct benchmark releases for GRCh37/GRCh38 and reference files. | https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/ |
| GIAB data indexes | Public raw reads and alignments across platforms. | https://github.com/genome-in-a-bottle/giab_data_indexes |
| GIAB/GA4GH stratifications | Difficult-region BEDs for GC, low-complexity, mappability and other contexts. | https://github.com/genome-in-a-bottle/genome-stratifications |
| Cancer Genome in a Bottle | Public tumor-normal HG008 resources and draft somatic benchmarks. | https://www.nist.gov/programs-projects/cancer-genome-bottle |

## What These Data Give Us

Public benchmark resources can provide:

- high-confidence/callable regions;
- truth variants;
- difficult-region stratifications;
- raw reads and alignments for deeper assay simulation;
- tumor-normal somatic benchmark regions for cancer workflows.

They usually do not directly provide a commercial lab-style assay performance
table with every target's depth, no-call rate, false negatives and false
positives. We derive that table from benchmark/callability files, or we compute
coverage and error metrics from raw alignments when needed.

## Current Implemented Path

Convert public benchmark/callable regions into a performance table:

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

Outputs:

- `demo_public_benchmark_performance.tsv`
- `demo_callability_conversion.md`
- `demo_callability_conversion.json`
- `demo_public_benchmark_enrichment.md`
- `demo_public_benchmark_enrichment.json`

## Scientific Interpretation

This gives us public benchmark validation, not final clinical validation.

If a region has low overlap with a high-confidence/callable BED, that means the
benchmark resource does not confidently assess that region. It is a strong clue
that the region may be technically difficult, but it is not proof that a
specific hospital assay fails there.

The evidence ladder is:

1. Sequence-derived risk: GC, homopolymer and low-complexity tracks.
2. Public benchmark support: risk regions are less callable or more error-prone
   in GIAB/GIAB Cancer resources.
3. Platform/run support: the same regions show low coverage or no-calls in real
   sequencing runs.
4. Clinical assay support: validated false negatives, false positives or
   orthogonal discordance in a real panel.

## Next Real Dataset Run

For a real public run:

1. Choose the evaluated BED: a panel design, oncology loci BED or benchmark
   target universe.
2. Download the matching GIAB or Cancer GIAB high-confidence/callable BED for
   the same reference build.
3. Run `benchmark_regions_to_performance.py`.
4. Run or reuse `empirical_technical_tracks.py` on matching FASTA/targets.
5. Compare technical-risk enrichment against low-callability regions.

For raw read validation, add a separate coverage extractor over BAM/CRAM files
and emit the same performance TSV schema. Large raw datasets should stay out of
the repository; only manifests, checksums and derived summaries should be
committed.
