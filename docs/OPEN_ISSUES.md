# Open Issues

Status: **PUBLIC CODE RELEASE v1.0.1 (2026-08-04)**. The canonical
repository is `/perioperative-glycemic-variability` on `main` and has been
pushed to the public remote
`jiahanzhe-123/perioperative-glycemic-variability`. All three GitHub Actions
workflows pass on the frozen review commit.

## Resolved in the current public release

- The canonical repository is the single editable copy. The desktop copy is a
  read-only backup and is not part of this repository.
- All analysis-of-record code is present, including the BMI plausibility rule,
  365-day pooled-interval and covariate time-transform sensitivities, INSPIRE
  v5 administrative-censoring analysis, eICU random-intercept sensitivity, and
  the final 13-exhibit figure pipeline.
- Frozen validation is **35/35 PASS** in the authorized environment. The
  repository also records 23/23 Python tests, passing R tests, a passing
  synthetic workflow, and syntax checks for the migrated R and Python files.
- `docs/FINAL_OUTPUT_PROVENANCE.csv` contains all 43 manuscript and
  supplementary exhibits. Figure source entries use repository-relative paths;
  no pending migration marker remains.
- The staged repository security review found no clinical data, credentials,
  symlinks, local configuration, or oversized analysis artifacts. Synthetic
  identifiers are simulated and are explicitly labelled as such.
- The GitHub repository was created privately, then prepared for public
  release. On the DOI-metadata commit `6d9aad939d518faae7e2f690bf5a3097ff3fe8f8`,
  the `documentation` (run `30915721274`), `lint` (run `30915718898`), and
  `tests` plus synthetic workflow (run `30915721243`) checks all passed.
- The author-confirmed manuscript work copy records that the INSPIRE secondary
  analysis was exempt from ethical review and that ChatGPT was used only for
  language polishing and code review without access to patient-level data.
- The local PostgreSQL credential was rotated on 2026-08-04 for role `postgres`
  on `inspire-pg` (`inspire_v142`); the password is stored only in ignored local
  configuration and is not tracked here.
- The Zenodo archive is complete. The concept DOI is
  [`10.5281/zenodo.21791846`](https://doi.org/10.5281/zenodo.21791846); the
  version DOI for `v1.0.1` is
  [`10.5281/zenodo.21791847`](https://doi.org/10.5281/zenodo.21791847).

## Remaining author-controlled items after public release

1. Add ORCID identifiers to the citation metadata if the authors wish to use
   them; none were supplied for this release, and this does not block code
   availability.
2. Confirm any remaining INSPIRE data-governance, coverage, and access-contact
   facts for the eventual manuscript submission.

The local PostgreSQL credential rotation was completed on 2026-08-04; no new
credential is tracked here.

## Non-blocking technical notes

- The eICU random-intercept model retains its documented convergence warning
  and is a sensitivity analysis rather than the primary harmonized estimate.
- Royston–Parmar diagnostics may depend on the local `flexsurv` availability;
  this does not change the frozen manuscript estimates.
- Superseded scripts remain under `archive/` for auditability and are not part
  of the analysis-of-record execution path.

See [RELEASE_CHECKLIST.md](RELEASE_CHECKLIST.md) for the gate-by-gate release
checklist and [PRIVACY_REVIEW.md](PRIVACY_REVIEW.md) for the data and secret
handling review.
