# Assay Performance Enrichment Report

- Performance intervals: **5**
- Failure intervals: **2**
- Universe bases: **500**
- Failure bases: **200**
- Technical-risk bases in universe: **400**

## Overall Enrichment

| Metric | Value |
|---|---:|
| risk_fail_bases | 200 |
| risk_pass_bases | 200 |
| nonrisk_fail_bases | 0 |
| nonrisk_pass_bases | 100 |
| risk_failure_rate | 0.5000 |
| nonrisk_failure_rate | 0.0000 |
| odds_ratio_haldane | 201.0000 |
| relative_risk | n/a |

## By Risk Category

| Category | Risk bases | Risk fail bases | Odds ratio | Relative risk |
|---|---:|---:|---:|---:|
| gc_extreme | 400 | 200 | 201.0000 | n/a |
| homopolymer | 400 | 200 | 201.0000 | n/a |

## Top Performance Rows

| Row | Locus | Failure | Mean depth | Callable fraction | Risk overlap |
|---|---|---|---:|---:|---:|
| demo_high_gc_region_2 | chrDemo:200-300 | True | n/a | 0.800 | 100 |
| demo_low_gc_region | chrDemo:100-200 | True | n/a | 0.400 | 100 |
| demo_high_gc_region | chrDemo:0-100 | False | n/a | 1.000 | 100 |
| demo_low_gc_region_2 | chrDemo:300-400 | False | n/a | 1.000 | 100 |
| demo_neutral_region | chrDemo:400-500 | False | n/a | 1.000 | 0 |

## Interpretation

This analysis tests whether empirical technical-risk intervals overlap observed assay failures more than expected over evaluated assay regions. It is an exploratory research statistic, not proof of causal failure.
