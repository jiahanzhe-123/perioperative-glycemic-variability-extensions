# v1.1.0 — JAHA v5 measurement-context extension

This release adds the final-lineage Phase 2A–2C and Phase 3B public-safe
extension package for the JAHA v5 revision.

Included:

- aggregate same-patient source-agreement summaries, log-scale sensitivity,
  source-count closure, and source-defined coefficient summaries;
- the locked 18-row SD-GV analytic-context landscape with provenance;
- aggregate-safe Figure 1 and Figure 3 files plus a public Figure 3 builder;
- release-scope, privacy, figure-source, and file-hash manifests;
- tests and synthetic checks inherited from the original public code.

Not included:

- MIMIC-IV, INSPIRE, or eICU-CRD raw/controlled data;
- patient/stay-level rows, real bootstrap replicates, or the Figure 2
  patient-level source values/derived point-level figures;
- manuscript/submission files, credentials, logs, or local absolute paths.

The aggregate Figure 3 builder can be run with:

```bash
Rscript analysis_extensions/scripts/phase3b/01_public_aggregate_landscape.R
```

The public package preserves the measurement-context interpretation: no
source is designated as a reference standard, source-count imbalance is not
treated as a causal explanation, and MIMIC/INSPIRE HRs are not pooled with
eICU RRs.
