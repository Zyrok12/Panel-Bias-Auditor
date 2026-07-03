# Empirical Technical Tracks

Handwritten proxy tracks are useful for small demos, but serious assay-design
review needs technical-risk tracks that can be reproduced from data.

The first empirical technical-track command is `gcderive`.

## GC Extremes

`gcderive` scans a FASTA file with sliding windows and emits merged high-GC and
low-GC regions as BED intervals. These intervals can be used as difficult-region
tracks in `audit`, `compare`, or downstream research workflows.

Example:

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

## Inputs

- FASTA sequence with chromosome/contig names.
- Optional BED file limiting the scan to panel, gene, or track-pack regions.
- Window size, step size, and GC thresholds.

If no BED file is supplied, all FASTA contigs are scanned.

## Outputs

- BED difficult-region track with `type=gc_extreme`.
- Markdown report.
- JSON report with windows and merged intervals.

The generated BED rows include:

- `gc_class`: `high_gc` or `low_gc`.
- `curation_status=source_verified`.
- `evidence_level=assay_computable_gc`.
- `min_gc` and `max_gc` across merged windows.
- `window_count`.

## Biological Meaning

GC extremes affect capture efficiency, amplification balance, local coverage
uniformity and variant-calling confidence. They do not automatically make a
region uncallable, but they identify places where assay validation and coverage
review deserve extra attention.

## Next Empirical Tracks

- `seqderive` homopolymer and low-complexity scanner.
- Mappability import from bigWig/BED sources.
- Segmental duplication/homology import.
- Fusion-breakpoint evidence tracks.

## Sequence Complexity

`seqderive` scans FASTA sequence for homopolymer-heavy and low-complexity
windows.

Example:

```bash
PYTHONPATH=src python3 -m panel_bias_auditor seqderive \
  --fasta examples/demo_reference.fa \
  --regions examples/demo_gc_scan_regions.bed \
  --window-size 40 \
  --step-size 40 \
  --min-homopolymer 10 \
  --low-complexity-fraction 0.8 \
  --out-bed examples/demo_sequence_risk.bed \
  --report examples/demo_sequence_risk_report.md \
  --json examples/demo_sequence_risk_report.json
```

For research use, `research/empirical_technical_tracks.py` runs GC and sequence
complexity together and writes a combined empirical technical-risk BED.

Use `research/validate_empirical_tracks.py` to validate that the generated BED
and JSON outputs are internally consistent.
