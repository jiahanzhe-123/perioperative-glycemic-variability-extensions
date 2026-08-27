# v1.1.1 — Event-limited sensitivity bundle

This release preserves the public v1.1.0 Phase 3B extension package and adds
the locked post hoc sensitivity bundle used for the Anaesthesia submission
candidate.

## Included

- The pre-run analysis lock documenting the three sensitivity analyses.
- The reproducible R runner for the parsimonious same-patient Cox, ridge-
  penalised same-patient Cox, and MIMIC priority-series analyses requiring at
  least three retained glucose measurements.
- Aggregate sensitivity estimates, paired-bootstrap summaries, the MIMIC
  multiple-imputation output, a run manifest, and a QC record.
- A release file manifest with SHA-256 hashes and sanitized provenance.

## Boundary

The repository contains no patient-level or stay-level records, controlled
inputs, real-data bootstrap replicates, credentials, local absolute paths, or
manuscript submission files. The released runner requires separately authorised
local inputs and does not bypass dataset access controls.

No versioned Zenodo DOI is claimed for v1.1.1 without an independently verified
archive record. The v1.1.0 Zenodo record remains the archival record for the
earlier public code snapshot.
