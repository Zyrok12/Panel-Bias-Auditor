# Panel Comparison Report

- Genome build: `GRCh38-demo`
- Track bundle: `demo-oncology-assay-risk-bundle` version `0.1.0`
- Panels compared: **2**

| Rank | Panel | Risk | Critical coverage | Missing critical bp | Difficult overlap bp | Target bp |
|---:|---|---:|---:|---:|---:|---:|
| 1 | Demo improved oncology panel | 25.0 (moderate) | 100.0% | 0 | 82,300 | 206,600 |
| 2 | Demo baseline oncology panel | 76.75 (severe) | 20.4% | 164,000 | 9,300 | 42,600 |

## Interpretation

Lower risk scores indicate better coverage of supplied critical regions with less difficult-region overlap.
This comparison is only as good as the supplied track bundle; production use requires curated, versioned tracks.
