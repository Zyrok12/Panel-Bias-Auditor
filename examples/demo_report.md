# Panel Bias Audit Report

- Generated: `2026-07-03T10:27:59.642295+00:00`
- Genome build: `GRCh38-demo`
- Tool: `panel-bias-auditor-lite`
- Track bundle: `demo-oncology-assay-risk-bundle` version `0.1.0`

## Executive Summary

- Risk score: **81.75 / 100** (severe)
- Panel intervals: **7**
- Merged target footprint: **42,600 bp**
- Critical-region coverage: **20.4%**
- Missing critical bases: **164,000 bp**
- Difficult-region overlap: **9,300 bp**
- Variants outside panel: **2 / 4**

## Risk Components

| Component | Value |
|---|---:|
| critical_missing_fraction | 0.7961 |
| difficult_region_fraction | 0.2183 |
| variant_outside_panel_fraction | 0.5 |

## Recommendations

| Priority | Category | Finding | Suggested action |
|---|---|---|---|
| high | critical-coverage | 164,000 critical bases are missing or partially covered. Top affected regions: BRCA1_core_demo, ALK_fusion_intron_proxy_demo, MET_exon14_and_amp_proxy_demo. | Review whether these regions are biologically required for the panel claim; add probes/amplicons or narrow the assay claim. |
| medium | callability | The largest difficult-region burden is homology (5,000 bp overlap). | Treat these regions as validation priorities; confirm sensitivity, specificity, and false-negative behavior with orthogonal or truth-set data. |
| medium | variant-check | 2 supplied variants fall outside the panel. Examples: chr2:29480000:A>G, chr13:32350000:C>T. | If these variants represent expected positives or customer requirements, revise targets or document them as out-of-scope. |

## Critical Regions

- Critical regions assessed: **8**
- Covered bases: **42,000 / 206,000 bp**
- Status counts: `{'covered': 4, 'partial': 2, 'missing': 2}`

### Missing Or Partial Critical Regions

| Region | Locus | Status | Covered | Missing |
|---|---|---|---:|---:|
| BRCA1_core_demo | chr17:43044000-43125000 | missing | 0 | 81,000 |
| ALK_fusion_intron_proxy_demo | chr2:29440000-29500000 | missing | 0 | 60,000 |
| MET_exon14_and_amp_proxy_demo | chr7:116770000-116800000 | partial | 10,000 | 20,000 |
| TP53_core_extended_demo | chr17:7668000-7675500 | partial | 4,500 | 3,000 |

## Difficult Region Overlap

| Category | Overlap bases |
|---|---:|
| homology | 5,000 |
| low_mappability | 2,500 |
| gc_extreme | 1,800 |

## Variant Checks

| Variant | In panel | Panel regions | Critical regions | Difficult regions |
|---|---|---|---|---|
| chr7:55141478:C>T | yes | EGFR_exon20_demo | EGFR_exon20_demo | EGFR_exon20_gc_proxy |
| chr17:7670700:G>A | yes | TP53_core_demo | TP53_core_extended_demo | TP53_mapping_complexity_proxy |
| chr2:29480000:A>G | no | - | ALK_fusion_intron_proxy_demo | ALK_large_intron_capture_risk |
| chr13:32350000:C>T | no | - | - | - |

## Interpretation Notes

- This report evaluates assay-design risk, not clinical validity.
- A region can be technically covered but still hard to call if it overlaps repeats, homologous sequence, GC extremes, or reference-biased loci.
- The quality of the audit depends on the quality of the critical-region and difficult-region tracks supplied.
- For commercial or clinical use, tracks must be curated, versioned, and validated against appropriate truth sets.
