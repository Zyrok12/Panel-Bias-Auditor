# Biological Context

## Why Panel Bias Matters

Targeted sequencing panels are often described by what they "cover": genes, exons, hotspots, introns, or promoter regions. But biological and clinical usefulness depends on more than nominal coverage.

A panel may include a region in its BED file while still producing weak or unreliable calls because of:

- Repetitive sequence.
- Homologous genes or pseudogenes.
- Segmental duplications.
- GC extremes.
- Large introns needed for fusion detection.
- Structural variation not represented well by a linear reference.
- Population-specific sequence absent from the reference genome.
- Assay chemistry limitations.

This distinction is the central biological idea behind the tool:

> Covered does not always mean confidently callable.

## Oncology Examples

In cancer genomics, missing or weakly calling a region can have direct downstream consequences.

Examples:

- **EGFR**: NSCLC treatment selection and resistance monitoring depend on small variants in specific exons. Some resistance variants are low-frequency and require strong assay sensitivity.
- **MET**: Exon 14 skipping and amplification biology can be difficult to capture with a simple SNV-focused panel.
- **ALK/ROS1/RET/NTRK fusions**: Fusion detection may require intronic or RNA-aware assay design. A panel that covers coding exons may still miss relevant rearrangements.
- **TP53**: Important across many cancers, but variant interpretation and detection can be affected by technical and biological complexity.
- **TERT promoter**: Clinically relevant promoter hotspots are outside coding regions, so exon-only designs can miss them.
- **BRCA1/BRCA2 and mismatch-repair genes**: Hereditary cancer risk and tumor treatment decisions can depend on confident coverage of clinically important regions.

## Pangenome And Reference Bias

Most clinical sequencing still uses GRCh37 or GRCh38-style linear references. These references are useful, but they do not represent all human diversity equally. If reads come from a locus where the patient's haplotype differs substantially from the reference, alignment and variant calling can become biased.

Pangenome resources try to reduce this by representing multiple haplotypes and structural forms of the genome. For an assay audit, the practical question is:

> Are there regions where a panel appears adequate against one reference but may perform unevenly across ancestries or haplotypes?

The public lite tool does not solve graph-genome analysis yet. Instead, it is built so pangenome-derived "risk tracks" can be supplied as difficult-region BED files.

## Why A Track-Driven Design Is Powerful

The auditor accepts user-supplied tracks:

- Critical regions: places the assay should cover well.
- Difficult regions: places where technical confidence may be reduced.
- Variants: known or test variants to check against the design.

That means the same engine can be used for:

- Oncology panels.
- Germline cancer panels.
- Rare disease panels.
- Pharmacogenomics panels.
- ctDNA or tumor-informed marker panels.
- Internal validation tracks from a lab or company.

The same engine can support different use cases when the track bundle is curated and versioned carefully.

## Connection To ctDNA Marker Design

Tumor-informed ctDNA assays select a set of patient-specific tumor markers and then track them in plasma. This depends on marker quality.

Good tumor-informed markers should ideally be:

- Tumor-specific.
- Clonal or truncal.
- Technically callable at very low allele fractions.
- Not likely to be clonal hematopoiesis.
- Not in regions with high sequencing or alignment noise.
- Distributed enough to reduce single-marker failure risk.

A marker-design workflow can reuse this foundation:

1. Audit candidate markers against difficult-region tracks.
2. Penalize CHIP-prone contexts.
3. Simulate detection at low tumor fraction.
4. Generate a design-risk report.

## What This Tool Is Not

This tool is not:

- A clinical diagnostic.
- A replacement for validation.
- A variant interpretation engine.
- A treatment recommendation system.
- A substitute for molecular pathology review.

It is an assay-design and research-use quality layer.
