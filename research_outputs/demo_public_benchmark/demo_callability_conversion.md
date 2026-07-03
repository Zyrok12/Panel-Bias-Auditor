# Public Benchmark Callability Conversion

- Source label: **demo_public_benchmark**
- Evaluated intervals: **5**
- Low-callability intervals: **2**
- Total evaluated bases: **500**
- Low-callability bases: **200**
- Failure threshold: callable_fraction < **0.900**
- Mean callable fraction: **0.840**

## Interval Preview

| Region | Locus | Callable fraction | Status |
|---|---|---:|---|
| demo_low_gc_region | chrDemo:100-200 | 0.400 | low_callability |
| demo_high_gc_region_2 | chrDemo:200-300 | 0.800 | low_callability |
| demo_high_gc_region | chrDemo:0-100 | 1.000 | pass |
| demo_low_gc_region_2 | chrDemo:300-400 | 1.000 | pass |
| demo_neutral_region | chrDemo:400-500 | 1.000 | pass |

## Interpretation

This table converts public benchmark high-confidence or callable regions into the same interval-level schema used for assay-performance enrichment. A low-callability row means the evaluated region is poorly represented by the benchmark/callable set, not that a specific clinical assay failed there.
