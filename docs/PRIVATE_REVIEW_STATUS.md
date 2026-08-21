# Private GitHub Review Status

**Historical status: READY FOR PRIVATE GITHUB REVIEW — superseded by the
public `v1.0.0` code release on 2026-08-04.**

This document records the private release-candidate state. It is a hand-off
summary for private repository review; it is not a public-release declaration.

## Verified locally

- One canonical `main` repository with a fresh history: 12 commits total (10
  project analysis/QC commits plus two private-review documentation syncs), pushed to
  the private remote [`jiahanzhe-123/perioperative-glycemic-variability`](https://github.com/jiahanzhe-123/perioperative-glycemic-variability).
- Analysis-of-record scripts are present, including the BMI repair path, the
  365-day pooled-interval and covariate time-transform sensitivities, INSPIRE
  v5 administrative censoring, eICU random-intercept sensitivity, and the
  13-exhibit figure pipeline.
- Frozen validation: **35/35 PASS** in the authorized environment.
- Python tests: **23/23 PASS**; R tests, synthetic workflow, and migrated-file
  syntax checks also pass.
- Exhibit provenance: **43/43** manuscript and supplementary exhibits mapped
  to repository-relative scripts or documented controlled-environment inputs.
- Staged security review: no patient-level clinical data, credentials, local
  configuration, symlinks, or oversized analysis artifacts.
- GitHub Actions on the current commit `7f5c81124b580588401067a2f602e03ecb356fb2`:
  **documentation PASS** (run `30887264620`), **lint PASS** (run
  `30887264657`), and **tests + synthetic workflow PASS** (run
  `30887264746`). The preceding analysis commit was also green.

## Deliberately not performed

- No public release, version tag, or Zenodo deposit has been made.
- No DOI or ORCID facts were invented. The MIT license, copyright holder,
  author names, contact e-mails, the author-confirmed INSPIRE statement that
  this secondary analysis was exempt from ethical review, and the factual
  ChatGPT-use declaration are recorded in the manuscript work copy. The local
  PostgreSQL password was rotated on 2026-08-04; the new value remains only in
  ignored local configuration and is not recorded here.

## Before public release

Use [RELEASE_CHECKLIST.md](RELEASE_CHECKLIST.md) and
[OPEN_ISSUES.md](OPEN_ISSUES.md) to close the author-controlled items. In
particular, confirm any remaining INSPIRE governance details, any ORCID
identifiers, and the eventual version/Zenodo decision.
