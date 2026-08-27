# JAHA v5 analysis extensions

This directory contains the final-lineage extension code and aggregate evidence package for the measurement-context analysis. It includes Phase 1.6 lineage checks, Phase 2A measurement-context analyses, Phase 2B analytic-context landscape analyses, Phase 2C evidence-lock assembly, the Phase 3B public aggregate production layer, and the v1.1.1 event-limited sensitivity bundle. The v1.1.1 release preserves the v1.1.0 materials and adds the event-limited sensitivity bundle.

## Public-release boundary

The public package contains code, aggregate tables, evidence registries, and candidate figures only. It does not contain MIMIC-IV, eICU-CRD, or INSPIRE data; patient/stay-level rows; real-data bootstrap replicates; local manifests; logs; manuscript files; or the private reviewer data package. See `provenance/release_scope.csv`.

## Reproduction

1. Copy `config.example.yaml` to `config.yaml`.
2. Replace placeholders with authorized local input and code paths.
3. Run the Phase 1.6 checks before any extension analysis.
4. Run the Phase 2A, Phase 2B, and Phase 2C R scripts in their numbered order.
5. To rebuild the public aggregate Figure 3, run `Rscript scripts/phase3b/01_public_aggregate_landscape.R`.

The locked aggregate outputs in `results/` are provenance-bearing records, not substitutes for the controlled source data. No result should be manually retyped into a manuscript.

## Key public outputs

- `results/machine_readable/phase2b_analytic_context_landscape.csv`: conceptual, SD-GV-only context landscape.
- `results/machine_readable/phase2b_log_scale_agreement.csv`: multiplicative-scale sensitivity.
- `results/machine_readable/phase2b_source_sampling_imbalance.csv`: descriptive count/span closure, including the declared NOT_ESTIMABLE span boundary.
- `results/evidence_lock/master_evidence_registry.csv`: claim-to-result provenance registry.
- `results/phase3b/`: final aggregate source-dependence closure, landscape rows,
  and public figure source data.
- `figures/phase3b/`: aggregate-safe Figure 1 and Figure 3 production files.
- `figures/phase3b_public/`: output location for the public aggregate builder.
- `provenance/phase3b/FINAL_FIGURE_SOURCE_MAP.csv`: panel-level source map.

## Event-limited sensitivity bundle

- `results/phase3b/anaesthesia_new_sensitivities/`: aggregate sensitivity
  results, bootstrap summaries, the at-least-three-measurement MICE result and
  the QC record.
- `scripts/phase3b/run_anaesthesia_new_sensitivities.R`: locked runner. It
  requires authorised local inputs specified through a private `config.yaml`;
  those inputs are not part of this release.
- `provenance/phase3b/ANAESTHESIA_NEW_SENSITIVITY_LOCK.md`: analysis lock and
  interpretation boundary.

## Interpretation boundary

The evidence supports a bounded measurement-context framing. It does not establish a reference-standard hierarchy between POCT and laboratory measurements, does not identify a causal sampling mechanism, does not treat count/span as confounders, and does not pool the non-equivalent database estimands.

The controlled local production package contains the patient-level Figure 2
scatter/Bland–Altman source values. Those rows and the final manuscript files
are intentionally excluded from this public repository. Public Figure 3 is
reproducible from aggregate rows only.
