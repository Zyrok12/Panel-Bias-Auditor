# Track Metadata And Batch Comparison

## Why Track Metadata Matters

An assay audit is only credible if the tracks are credible.

For a real analysis, it is not enough to say:

> This panel misses 20,000 bases.

The reviewer will ask:

- Which biological regions were considered important?
- Who curated those regions?
- Which genome build was used?
- Which version of the track was used?
- What license or source did the track come from?
- Can the same audit be reproduced next quarter?

The track manifest answers those questions.

## Track Bundle Format

`demo_track_manifest.json` is an example:

```json
{
  "name": "demo-oncology-assay-risk-bundle",
  "version": "0.1.0",
  "genome_build": "GRCh38-demo",
  "description": "Illustrative oncology assay-risk tracks for demos only.",
  "tracks": [
    {
      "name": "demo critical oncology regions",
      "role": "critical",
      "path": "demo_critical_regions.bed",
      "description": "Toy cancer hotspot, hereditary cancer, and fusion proxy regions.",
      "source": "synthetic demo coordinates",
      "version": "0.1.0",
      "license": "MIT demo data"
    }
  ]
}
```

Supported roles:

- `critical`: regions the assay should cover confidently.
- `difficult`: regions where coverage may not equal callability.

## Batch Comparison

The `compare` command ranks panels under the same track bundle.

This is useful because users often want to compare:

- Their panel against a competitor panel.
- Version 1 versus version 2 of their panel.
- A small hotspot panel versus a broader comprehensive genomic profiling panel.
- A proposed ctDNA marker set versus alternative designs.

Example:

```bash
PYTHONPATH=src python3 -m panel_bias_auditor compare \
  --panels examples/demo_panels.tsv \
  --track-manifest examples/demo_track_manifest.json \
  --out examples/demo_comparison.md \
  --json examples/demo_comparison.json \
  --tsv examples/demo_comparison.tsv
```

The output table ranks panels by lower assay-design risk.
