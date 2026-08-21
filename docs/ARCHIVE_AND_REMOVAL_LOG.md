# Archive and Removal Log

Nothing was silently deleted. Superseded scripts were moved to
`archive/not_used_in_final_analysis/`; scratch artifacts were excluded from
the repository entirely (not archived). Full old→new mapping of migrated
files: `docs/FILE_MIGRATION_MAP.csv`.

## Archived (kept for audit, not used in any final result)

| File | Origin | Reason archived | Superseded by |
|---|---|---|---|
| `archive/not_used_in_final_analysis/mimic_03_primary_models_superseded_by_03b.R` | internal rebuild workspace | first-pass primary-model script; replaced by the MICE rerun used for frozen results | `analyses/03_primary_mimic/03_run_primary_models_mice.R` |
| `archive/not_used_in_final_analysis/inspire_11_superseded_by_15_16_17.R` | INSPIRE workspace | pre-repair INSPIRE model script (outcome/time structure pre-dated the reconciliation) | `analyses/06_inspire/02_outcome_time_repairs.R`, `03_sensitivity_joint_repairs.R`, `04_final_models_v3.R` |
| `archive/not_used_in_final_analysis/inspire_12_superseded_by_16.R` | INSPIRE workspace | intermediate sensitivity script | `analyses/06_inspire/03_sensitivity_joint_repairs.R` |
| `archive/not_used_in_final_analysis/inspire_13_superseded_by_16.R` | INSPIRE workspace | intermediate sensitivity script | `analyses/06_inspire/03_sensitivity_joint_repairs.R` |

## Excluded from the repository (not archived)

| Category | Examples | Reason |
|---|---|---|
| Legacy public draft | `github_release/` inside the old supplementary-analysis workspace | superseded in full by this repository; kept locally |
| Internal Word/Supplement tooling | replication-workspace Python helpers for DOCX tables | manuscript-production utilities, not analysis code |
| Patient-level frames and outputs | analysis frames, MICE objects, bootstrap outputs | DUA-restricted; must never be published |
| Local full archive packages | desktop `*_final_package_*` bundles | contain patient-level data; local-only |
| OS/session noise | `.DS_Store`, `__pycache__`, `.Rhistory`, `.pytest_cache` | noise |

## Policy going forward

- Archive entries are read-only. If a fix is needed, create a new versioned
  file in the module and log the change in `CHANGELOG.md`.
- Any future removal must be appended to this log with a reason.
