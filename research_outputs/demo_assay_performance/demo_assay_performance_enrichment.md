# Assay Performance Enrichment Report

- Performance intervals: **6**
- Failure intervals: **4**
- Universe bases: **600**
- Failure bases: **400**
- Technical-risk bases in universe: **400**

## Overall Enrichment

| Metric | Value |
|---|---:|
| risk_fail_bases | 400 |
| risk_pass_bases | 0 |
| nonrisk_fail_bases | 0 |
| nonrisk_pass_bases | 200 |
| risk_failure_rate | 1.0000 |
| nonrisk_failure_rate | 0.0000 |
| odds_ratio_haldane | 321201.0000 |
| relative_risk | n/a |

## By Risk Category

| Category | Risk bases | Risk fail bases | Odds ratio | Relative risk |
|---|---:|---:|---:|---:|
| gc_extreme | 400 | 400 | 321201.0000 | n/a |
| homopolymer | 400 | 400 | 321201.0000 | n/a |

## Top Performance Rows

| Row | Locus | Failure | Mean depth | Callable fraction | Risk overlap |
|---|---|---|---:|---:|---:|
| demo_high_gc_dropout | chrDemo:0-100 | True | 35.00 | 0.520 | 100 |
| demo_high_gc_fp_cluster | chrDemo:200-300 | True | 85.00 | 0.880 | 100 |
| demo_low_gc_dropout | chrDemo:100-200 | True | 42.00 | 0.610 | 100 |
| demo_low_gc_low_callability | chrDemo:300-400 | True | 60.00 | 0.700 | 100 |
| demo_clean_control_1 | chrDemo:400-500 | False | 520.00 | 0.990 | 0 |
| demo_clean_control_2 | chrDemo:500-600 | False | 510.00 | 0.980 | 0 |

## Interpretation

This analysis tests whether empirical technical-risk intervals overlap observed assay failures more than expected over evaluated assay regions. It is an exploratory research statistic, not proof of causal failure.
