# Archive: not used in the final analysis

Scripts in this directory were superseded during development and are kept
**only for audit**. They are not referenced by any final result, are not part
of any run target, and must not be executed to regenerate results.

| File | Why it exists here | Replacement |
|---|---|---|
| `mimic_03_primary_models_superseded_by_03b.R` | First-pass MIMIC primary-model script; the frozen results come from the later MICE rerun | `analyses/03_primary_mimic/03_run_primary_models_mice.R` |
| `inspire_11_superseded_by_15_16_17.R` | INSPIRE models built before the outcome/time-structure reconciliation | `analyses/06_inspire/02_outcome_time_repairs.R`, `03_sensitivity_joint_repairs.R`, `04_final_models_v3.R` |
| `inspire_12_superseded_by_16.R` | Intermediate INSPIRE sensitivity version | `analyses/06_inspire/03_sensitivity_joint_repairs.R` |
| `inspire_13_superseded_by_16.R` | Intermediate INSPIRE sensitivity version | `analyses/06_inspire/03_sensitivity_joint_repairs.R` |

Archived files may still contain outdated assumptions (e.g., the pre-repair
INSPIRE outcome definition). Read them only alongside
`docs/ARCHIVE_AND_REMOVAL_LOG.md` and the final reports.
