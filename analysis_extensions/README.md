# JAHA v5 analysis extensions

This directory contains the final-lineage extension code and aggregate evidence package for the measurement-context analysis. It includes Phase 1.6 lineage checks, Phase 2A measurement-context analyses, Phase 2B analytic-context landscape analyses, and Phase 2C evidence-lock assembly.

## Public-release boundary

The public package contains code, aggregate tables, evidence registries, and candidate figures only. It does not contain MIMIC-IV, eICU-CRD, or INSPIRE data; patient/stay-level rows; real-data bootstrap replicates; local manifests; logs; manuscript files; or the private reviewer data package. See `provenance/release_scope.csv`.

## Reproduction

1. Copy `config.example.yaml` to `config.yaml`.
2. Replace placeholders with authorized local input and code paths.
3. Run the Phase 1.6 checks before any extension analysis.
4. Run the Phase 2A, Phase 2B, and Phase 2C R scripts in their numbered order.

The locked aggregate outputs in `results/` are provenance-bearing records, not substitutes for the controlled source data. No result should be manually retyped into a manuscript.

## Key public outputs

- `results/machine_readable/phase2b_analytic_context_landscape.csv`: conceptual, SD-GV-only context landscape.
- `results/machine_readable/phase2b_log_scale_agreement.csv`: multiplicative-scale sensitivity.
- `results/machine_readable/phase2b_source_sampling_imbalance.csv`: descriptive count/span closure, including the declared NOT_ESTIMABLE span boundary.
- `results/evidence_lock/master_evidence_registry.csv`: claim-to-result provenance registry.
- `figures/`: candidate source-dependence and analytic-context figures.

## Interpretation boundary

The evidence supports a bounded measurement-context framing. It does not establish a reference-standard hierarchy between POCT and laboratory measurements, does not identify a causal sampling mechanism, does not treat count/span as confounders, and does not pool the non-equivalent database estimands.
