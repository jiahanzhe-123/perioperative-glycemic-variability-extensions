# Anaesthesia candidate: new sensitivity-analysis lock

**Lock status:** `LOCKED_BEFORE_RUN`
**Lock date:** 2026-08-27
**Workspace:** `JAHA_v5_analysis_extensions`
**Scope:** Anaesthesia candidate only. The locked manuscript, authoritative public code, controlled inputs, historical branches and other journal candidates are not modified by these analyses.

## Purpose

These analyses are registered to address the event-limited, high-dimensional
same-patient mortality comparison and the instability of a two-measurement
standard deviation. They are sensitivity analyses, not new primary tests,
not a specification search, and not a basis for changing the evidence-lock
chronology. No manuscript number is replaced before the run outputs pass QC.

## Locked inputs

The analysis uses the final same-lineage MIMIC-derived inputs already named in
the extension configuration. The following SHA-256 values were recorded before
the run:

| Input | SHA-256 |
|---|---|
| `analysis_base_bmi_repaired.csv` | `2163168d1c2cf07828dd10d3b79392ca94545fbf36e7fd68dcd9d2675eeefba4` |
| `features_priority.csv` | `649dcf92d9f2fa35bab2521b123d54b1390ec508c1724ab8911ba1d0e599d906` |
| `samepatient_poct_lab.csv` | `9ff80386a5145f4cdc664a1ffe67db0e43a2bda08c93351bd81d6a9623b8774c` |
| `standardization_constants.json` | `f67cda9ad42193e07472d979f44984284545b70fa8932b02e2fc99bae4d9e27e` |
| locked primary `mice_pooled_models.csv` | `41c6a64f707f40518b5c809e67b26b23e9b3ade65ddea82708d5c5c1ce1c4d0a` |
| locked source `measurement_source_results.csv` | `ebb489c3f31bc879aa24e9ce98b6b5da66232a0312700c9ba237153c87180909` |

Data are merged by `stay_id`. Boolean landmark eligibility is interpreted
using the final-lineage truth values. The primary final target is restricted
to landmark-eligible patients with non-missing priority-series GV and at least
two retained measurements. The source-paired analysis starts from the
validated same-patient source file and uses the identical complete-case
intersection for both source models.

## Locked analyses

### 1. Parsimonious same-patient Cox sensitivity

- Outcome: 30-day all-cause mortality from the final day-1 landmark frame.
- Cohort: the same complete-case paired source cohort used by the locked
  source models; expected `N=409`, `events=49`.
- Exposures: source-specific GV and source-specific mean glucose, fitted
  separately for POCT-derived and laboratory-derived values.
- Model: Efron-tie Cox model with linear terms for GV (reported per
  `10 mg.dl⁻¹`), mean glucose, age, sex, diabetes and six-level procedure
  category. No BMI, Charlson, lactate, creatinine, SOFA, spline expansion,
  measurement-count terms or source-process terms are added.
- Inference: model-based 95% Wald confidence intervals for each source
  coefficient; a paired coefficient contrast is the POCT minus laboratory
  log-coefficient. A patient-level bootstrap with 2000 requested successful
  replicates is used for the paired contrast and its percentile interval
  (`seed=20260928`).

### 2. Ridge-penalised same-patient Cox sensitivity

- Outcome and cohort: identical to Analysis 1 (`N=409`, `events=49`).
- Design: the locked Model B clinical covariate frame, with source-specific GV
  and source-specific mean glucose, fixed restricted-cubic-spline mean-glucose
  knots from `standardization_constants.json`, and the locked Model B clinical
  covariates. No Model C measurement-process terms are added.
- Penalisation: ridge Cox (`alpha=0`) using the same fixed patient-level
  10-fold cross-validation folds for both sources (`fold seed=20260929`);
  `lambda.1se` is selected separately for each source from the full paired
  cohort. The coefficient is
  reported per `10 mg.dl⁻¹` higher source-defined GV.
- Inference: the full-cohort ridge coefficients and source-specific HRs are
  reported. For paired comparison, the full-cohort lambdas are held fixed and
  a patient-level bootstrap with 2000 requested successful replicates produces
  percentile intervals for each source coefficient and for the paired
  POCT-minus-laboratory log-coefficient difference (`seed=20260930`). Ridge
  coefficients and
  their contrast are descriptive stability summaries, not causal interactions.

### 3. MIMIC priority-series sensitivity requiring at least three measurements

- Cohort: final priority-series target restricted to `glucose_count >= 3`.
- Expected pre-run frame: `N=10398`, `30-day events=272`; the run must verify
  these values rather than assume them.
- Model: the locked Model B formula and fixed mean-glucose spline knots,
  including the same clinical covariate set and no post hoc significance-based
  sample selection.
- Missing data: the same MICE architecture as the locked primary analysis,
  `m=50`, `maxit=20`, predictive mean matching, Rubin pooling, with the new
  deterministic seed `20260827` recorded in the run manifest.
- Effect scale: GV per `0.555 mmol.l⁻¹ (10 mg.dl⁻¹)` higher. This is a
  measurement-count eligibility sensitivity, not a claim that three values
  establish a stable individual-level GV construct.

## Fixed interpretation rules

1. These analyses do not alter the prespecified primary analysis or promote
   post hoc same-patient findings to primary evidence.
2. A stable direction or magnitude across the parsimonious and ridge models
   supports robustness to covariate complexity only; it does not remove
   overfit risk or establish a causal source interaction.
3. The `>=3` measurement analysis addresses a minimum-measurement boundary;
   it does not make routine-care sampling equivalent to protocolised or CGM
   measurement.
4. Failed fits, non-convergence, insufficient bootstrap replicates or input
   hash mismatches are reported as QC failures and are not silently repaired.

## Locked outputs

The run writes only inside this extension workspace:

- `new_sensitivity_results.csv`
- `new_sensitivity_paired_bootstrap.csv`
- `new_sensitivity_mimic_ge3_mice.csv`
- `new_sensitivity_table.csv`
- `new_sensitivity_run_manifest.json`
- `new_sensitivity_qc.md`

The Anaesthesia manuscript and its submitted supporting files are not changed
by the analysis runner.
