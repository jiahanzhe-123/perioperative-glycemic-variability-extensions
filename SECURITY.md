# Security Policy

## Scope

This repository contains analysis code only. It must never contain:

- patient-level data from MIMIC-IV, eICU-CRD, or INSPIRE;
- database credentials, API keys, tokens, or SSH keys;
- local absolute paths, server names, or internal network addresses.

## Reporting a vulnerability or data exposure

If you believe patient-level data, credentials, or any other sensitive material
has been committed to this repository (including in git history), please report
it privately to the maintainers at **jhz0223@outlook.com**
instead of opening a public issue.

We will acknowledge receipt as soon as possible, remove the material, rewrite
history if necessary, and revoke any exposed credentials.

## For contributors

Run `make lint` before every commit. It scans the working tree for credential
patterns, absolute local paths, private IP addresses, and oversized files.
Never commit files under `data/` (other than `data/synthetic/`) or
`config/paths.yml` — both are excluded by `.gitignore`.
