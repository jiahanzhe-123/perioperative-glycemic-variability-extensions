#!/usr/bin/env python3
"""Run the Phase 1.6 final-lineage hard gate and write the audit report."""

from __future__ import annotations

import csv
import json
import re
import subprocess
import sys
from collections import Counter
from pathlib import Path
from typing import Any

WORKSPACE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(WORKSPACE))

from src.config_loader import load_config, output_root, resolve_config_path, workspace_root  # noqa: E402
from src.final_lineage import CONTEXT_FIELDS, reconstruct_all  # noqa: E402
from src.path_manifest import INPUT_KEYS, sha256_file  # noqa: E402
from src.schema_checks import read_csv, read_header, utc_now_iso, write_json  # noqa: E402


def _number(value: object) -> float | None:
    if value is None or str(value).strip().lower() in {"", "na", "nan", "null", "none"}:
        return None
    return float(str(value))


def _equal(value: object, expected: float, tolerance: float = 1e-6) -> bool:
    parsed = _number(value)
    return parsed is not None and abs(parsed - expected) <= tolerance


def _find(rows: list[dict[str, str]], key: str, expected: str) -> dict[str, str]:
    for row in rows:
        if row.get(key) == expected:
            return row
    raise KeyError(f"No row with {key}={expected}")


def _find_prefix(rows: list[dict[str, str]], key: str, prefix: str) -> dict[str, str]:
    for row in rows:
        if row.get(key, "").startswith(prefix):
            return row
    raise KeyError(f"No row with {key} starting with {prefix}")


def _git_head(path: Path) -> str | None:
    try:
        result = subprocess.run(["git", "-C", str(path), "rev-parse", "HEAD"], capture_output=True, text=True, check=True)
    except (OSError, subprocess.CalledProcessError):
        return None
    return result.stdout.strip() or None


def _modeling_tokens_present() -> list[str]:
    banned = (
        r"\bcoxph\s*\(",
        r"\bglm\s*\(",
        r"\bglmer\s*\(",
        r"\blifelines\b",
        r"\bstatsmodels\b",
        r"\bbootstrap\s*\(",
    )
    hits: list[str] = []
    for path in list((WORKSPACE / "src").glob("*.py")) + list((WORKSPACE / "scripts").glob("*.py")):
        text = path.read_text(encoding="utf-8")
        for pattern in banned:
            if re.search(pattern, text, flags=re.IGNORECASE):
                hits.append(f"{path.name}:{pattern}")
    return hits


def _historical_reference(config: dict[str, Any]) -> dict[str, Any]:
    expected = {
        "outputs/machine_readable/measurement_context_patient_level.csv": "6967b5e3e220f9c861f1a5e38717f21d5a7fa8fc2e80f4c1ab646bb201fdd091",
        "outputs/qc/measurement_context_qc_summary.json": "5ee2843525b300a13bd293fea326ec71c2cbb6a12e7832cce71be20303a75c41",
    }
    root = resolve_config_path(config, "historical_extension_root")
    current: dict[str, str | None] = {}
    for relative in expected:
        path = root / relative
        current[relative] = sha256_file(path) if path.exists() and path.is_file() else None
    return {
        "root": str(root.resolve()),
        "expected_sha256": expected,
        "current_sha256": current,
        "status": "PASS" if current == expected else "FAIL",
    }


def _write_report(
    config: dict[str, Any],
    decision: str,
    checks: list[dict[str, Any]],
    reconstruction: dict[str, Any],
    inventory_rows: list[dict[str, str]],
    locked: dict[str, Any],
    historical: dict[str, Any],
) -> Path:
    m = reconstruction["mimic"]
    s = reconstruction["samepatient"]
    e = reconstruction["eicu"]
    i = reconstruction["inspire"]
    status_by = Counter(row.get("fit_status", "NA") for row in inventory_rows)
    check_table = "\n".join(
        f"| {check['name']} | {check['status']} | {check.get('detail', '')} |" for check in checks
    )
    inventory_table = "\n".join(f"| {key} | {value} |" for key, value in sorted(status_by.items()))
    lines = [
        "# Phase 1.6 final-lineage extension QC",
        "",
        f"- Generated: `{utc_now_iso()}`",
        f"- Workspace: `{workspace_root(config).resolve()}`",
        f"- Random seed recorded in config: `{config['random_seed']}`",
        "- Modeling performed: **NO**; all operations are deterministic reconstruction, cohort-mask validation, provenance, or locked-result fingerprint checks.",
        "",
        "## Final gate",
        "",
        f"`{decision}`",
        "",
        "## Check matrix",
        "",
        "| Check | Status | Detail |",
        "|---|---|---|",
        check_table,
        "",
        "## MIMIC-IV final-lineage reconstruction",
        "",
        f"- Final target: **N={m['target_n']}**, 30-day events **{m['events_30d']}**, 365-day events **{m['events_365d']}**.",
        f"- Measurement count Q1/median/Q3: **{m['count_q1']:.0f} / {m['count_median']:.0f} / {m['count_q3']:.0f}**.",
        f"- Reconstructed measurement span Q1/median/Q3: **{m['span_q1']:.10f} / {m['span_median']:.10f} / {m['span_q3']:.10f} h**.",
        f"- Reconstructed GV mean: **{m['gv_mean']:.10f}**; reconstructed mean-glucose mean: **{m['mean_glucose_mean']:.10f}**.",
        f"- Stored-vs-reconstructed maximum absolute differences: mean glucose **{m['max_mean_absolute_difference']:.3g}**, GV **{m['max_gv_absolute_difference']:.3g}**, span **{m['max_span_absolute_difference']:.3g} h**.",
        f"- Stored/recomputed span ratio deviation maximum: **{m['max_stored_recomputed_span_ratio_difference']:.3g}**; no x1000 correction was applied.",
        "",
        "### Locked MIMIC fingerprints",
        "",
        f"- Model B 30d HR: **{locked['mimic_model_b_hr']:.10f}**.",
        f"- Model C 30d HR: **{locked['mimic_model_c_hr']:.10f}**.",
        f"- POCT source 30d HR: **{locked['poct_hr']:.10f}**.",
        f"- Laboratory source 30d HR: **{locked['lab_hr']:.10f}**.",
        "",
        "## Same-patient source context",
        "",
        f"- Candidate pairs in final source file: **{s['candidate_n']}**.",
        f"- Exact agreement-eligible pairs with non-missing POCT and laboratory GV: **{s['agreement_eligible_n']}**; within final MIMIC target: **{s['agreement_eligible_final_mimic_n']}**.",
        f"- Outcome-model-aligned complete-case cohort: **N={s['outcome_model_n']}**, 30-day events **{s['outcome_model_events_30d']}**, 365-day events **{s['outcome_model_events_365d']}**.",
        "- Agreement eligibility and outcome-model eligibility are separate flags in the long-format context table; no correlation or Bland–Altman analysis was run.",
        "",
        "## eICU final-input limitation",
        "",
        f"- Aggregate input rows: **{e['aggregate_input_n']}**; deterministic main/landmark mask: **{e['main_landmark_n']}** rows and **{e['main_landmark_events']}** events.",
        f"- M3: **N={e['m3_n']}**, events **{e['m3_events']}**, locked RR **{locked['eicu_m3_rr']:.10f}**.",
        f"- M4: **N={e['m4_n']}**, events **{e['m4_events']}**, locked RR **{locked['eicu_m4_rr']:.10f}**.",
        "- No separate event-level eICU glucose file is preserved. This workspace reconstructs aggregate eligibility masks only and leaves event-level reconstruction explicitly NOT FEASIBLE.",
        "- eICU day-30/day-365 event fields in `final_measurement_context.csv` are `NA`, not fabricated.",
        "",
        "## INSPIRE final cohort and coverage",
        "",
        f"- 24h landmark: **N={i['24h_n']}**, events **{i['24h_events']}**, locked I2 HR **{locked['inspire_24h_hr']:.10f}**.",
        f"- 48h landmark: **N={i['48h_n']}**, events **{i['48h_events']}**, locked I2 HR **{locked['inspire_48h_hr']:.10f}**.",
        f"- Coverage metadata preserved as `{i['coverage_status']}` with `patient_level_coverage_end_field={str(i['patient_level_coverage_end_field']).upper()}`.",
        "- Historical/withdrawn 365-day endpoint is not reintroduced.",
        "",
        "## Phase 2 estimand inventory",
        "",
        "| Fit status | Rows |",
        "|---|---:|",
        inventory_table,
        "",
        "The inventory keeps MIMIC effects as HR, eICU effects as RR, and INSPIRE effects as HR. It explicitly prohibits cross-database pooling/interchangeability and records the eICU event-level request as NOT FEASIBLE until a valid event-level input exists.",
        "",
        "## Historical workspace and provenance",
        "",
        f"- Historical extension root checked read-only: `{historical['root']}`; reference hash status: **{historical['status']}**.",
        f"- Input manifest: `{(workspace_root(config) / 'manifests' / 'final_input_manifest.json').resolve()}`.",
        f"- Rebase manifest: `{(workspace_root(config) / 'manifests' / 'rebase_manifest.json').resolve()}`.",
        "- Authoritative roots were not moved, copied, or rewritten.",
        "",
    ]
    path = output_root(config) / "qc" / "final_lineage_extension_qc.md"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines), encoding="utf-8")
    return path


def main() -> int:
    config = load_config()
    checks: list[dict[str, Any]] = []
    failures: list[str] = []
    try:
        reconstruction = reconstruct_all(config)
    except Exception as exc:
        failure = {"name": "deterministic_reconstruction", "status": "FAIL", "detail": f"{type(exc).__name__}: {exc}"}
        checks.append(failure)
        failures.append(failure["detail"])
        reconstruction = {
            "mimic": {"target_n": 0, "events_30d": 0, "events_365d": 0, "count_q1": 0, "count_median": 0, "count_q3": 0, "span_q1": 0, "span_median": 0, "span_q3": 0, "gv_mean": 0, "mean_glucose_mean": 0, "max_mean_absolute_difference": 0, "max_gv_absolute_difference": 0, "max_span_absolute_difference": 0, "max_stored_recomputed_span_ratio_difference": 0},
            "samepatient": {"candidate_n": 0, "agreement_eligible_n": 0, "agreement_eligible_final_mimic_n": 0, "outcome_model_n": 0, "outcome_model_events_30d": 0, "outcome_model_events_365d": 0},
            "eicu": {"aggregate_input_n": 0, "main_landmark_n": 0, "main_landmark_events": 0, "m3_n": 0, "m3_events": 0, "m4_n": 0, "m4_events": 0, "event_level_input_available": False},
            "inspire": {"24h_n": 0, "24h_events": 0, "48h_n": 0, "48h_events": 0, "coverage_status": "FAIL", "patient_level_coverage_end_field": False},
        }
    else:
        checks.append({"name": "deterministic_reconstruction", "status": "PASS", "detail": "MIMIC, same-patient, eICU, and INSPIRE masks reconstructed from declared final inputs"})

    summary = reconstruction
    m = summary["mimic"]
    s = summary["samepatient"]
    e = summary["eicu"]
    i = summary["inspire"]
    target_checks = [
        ("MIMIC target/events", m.get("target_n") == 10561 and m.get("events_30d") == 296 and m.get("events_365d") == 745, f"N={m.get('target_n')}, events30={m.get('events_30d')}, events365={m.get('events_365d')}"),
        ("MIMIC span/reconstruction", _equal(m.get("span_median"), 14.833333333333334, 1e-10) and _equal(m.get("max_span_absolute_difference"), 0.0, 1e-12) and (m.get("max_gv_absolute_difference") or 1) <= 1e-8, f"median={m.get('span_median')}, max span diff={m.get('max_span_absolute_difference')}, max GV diff={m.get('max_gv_absolute_difference')}"),
        ("Same-patient candidate/agreement/model", s.get("candidate_n") == 453 and s.get("agreement_eligible_n") == 453 and s.get("outcome_model_n") == 409 and s.get("outcome_model_events_30d") == 49, f"candidate={s.get('candidate_n')}, agreement={s.get('agreement_eligible_n')}, model={s.get('outcome_model_n')}/{s.get('outcome_model_events_30d')}"),
        ("eICU aggregate masks", e.get("m3_n") == 7115 and e.get("m3_events") == 130 and e.get("m4_n") == 7115 and e.get("m4_events") == 130, f"M3={e.get('m3_n')}/{e.get('m3_events')}, M4={e.get('m4_n')}/{e.get('m4_events')}"),
        ("eICU event-level limitation", e.get("event_level_input_available") is False, "event-level input unavailable; aggregate-only validation"),
        ("INSPIRE 24h/48h masks", i.get("24h_n") == 1353 and i.get("24h_events") == 27 and i.get("48h_n") == 1511 and i.get("48h_events") == 31, f"24h={i.get('24h_n')}/{i.get('24h_events')}, 48h={i.get('48h_n')}/{i.get('48h_events')}"),
        ("INSPIRE coverage metadata", i.get("coverage_status") == "SOURCE_PROVENANCE_ATTESTED_TO_DAY30" and i.get("patient_level_coverage_end_field") is False, f"coverage={i.get('coverage_status')}, patient_level_coverage_end_field={i.get('patient_level_coverage_end_field')}"),
    ]
    for name, passed, detail in target_checks:
        check = {"name": name, "status": "PASS" if passed else "FAIL", "detail": detail}
        checks.append(check)
        if not passed:
            failures.append(f"{name}: {detail}")

    locked: dict[str, Any] = {}
    try:
        primary = read_csv(resolve_config_path(config, "mimic_primary_results"))
        mice = read_csv(resolve_config_path(config, "mimic_mice_results"))
        source = read_csv(resolve_config_path(config, "mimic_source_results"))
        eicu = read_csv(resolve_config_path(config, "eicu_locked_results"))
        inspire24 = read_csv(resolve_config_path(config, "inspire_primary_results"))
        inspire48 = read_csv(resolve_config_path(config, "inspire_48h_results"))
        model_b = _find_prefix(primary, "item", "PRIMARY TEST: GV per 10")
        model_c = _find(mice, "model_id", "MICE_C_30d")
        poct = _find(source, "model_id", "SRC_samepatient_poct_30d")
        lab = _find(source, "model_id", "SRC_samepatient_lab_30d")
        eicu_m3 = _find(eicu, "model_id", "HARM_eICU-CRD_M3")
        eicu_m4 = _find(eicu, "model_id", "HARM_eICU-CRD_M4")
        insp24 = _find(inspire24, "model_id", "ADMINV5_I2_30d")
        insp48 = _find(inspire48, "model_id", "ADMINV5_I2_48h_landmark_30d")
        locked = {
            "mimic_model_b_hr": float(model_b["HR"]),
            "mimic_model_c_hr": float(model_c["HR_per10"]),
            "poct_hr": float(poct["HR_per10"]),
            "lab_hr": float(lab["HR_per10"]),
            "eicu_m3_rr": float(eicu_m3["RR_per10"]),
            "eicu_m4_rr": float(eicu_m4["RR_per10"]),
            "inspire_24h_hr": float(insp24["HR_per10"]),
            "inspire_48h_hr": float(insp48["HR_per10"]),
        }
        fingerprint_checks = [
            ("MIMIC Model B fingerprint", locked["mimic_model_b_hr"], 0.979380602603274),
            ("MIMIC Model C fingerprint", locked["mimic_model_c_hr"], 1.02857335963722),
            ("POCT source fingerprint", locked["poct_hr"], 0.758857245048576),
            ("laboratory source fingerprint", locked["lab_hr"], 0.970897006154936),
            ("eICU M3 fingerprint", locked["eicu_m3_rr"], 1.15774790083148),
            ("eICU M4 fingerprint", locked["eicu_m4_rr"], 1.13279595096984),
            ("INSPIRE 24h I2 fingerprint", locked["inspire_24h_hr"], 0.904554212678046),
            ("INSPIRE 48h I2 fingerprint", locked["inspire_48h_hr"], 1.10438443373211),
        ]
        for name, actual, expected in fingerprint_checks:
            passed = abs(actual - expected) <= 1e-8
            check = {"name": name, "status": "PASS" if passed else "FAIL", "detail": f"actual={actual:.15g}, expected={expected:.15g}"}
            checks.append(check)
            if not passed:
                failures.append(f"{name}: actual={actual}, expected={expected}")
    except Exception as exc:
        check = {"name": "locked-result fingerprints", "status": "FAIL", "detail": f"{type(exc).__name__}: {exc}"}
        checks.append(check)
        failures.append(check["detail"])

    context_path = output_root(config) / "machine_readable" / "final_measurement_context.csv"
    inventory_path = output_root(config) / "machine_readable" / "phase2_estimand_inventory.csv"
    inventory_rows: list[dict[str, str]] = []
    if context_path.exists() and inventory_path.exists():
        context_rows = read_csv(context_path)
        inventory_rows = read_csv(inventory_path)
        context_header = read_header(context_path)
        context_pass = all(field in context_header for field in CONTEXT_FIELDS)
        eicu_context = [row for row in context_rows if row.get("database") == "eICU-CRD"]
        no_fake_eicu_events = all(row.get("event_30d") == "NA" and row.get("event_365d") == "NA" for row in eicu_context)
        no_modeling_tokens = _modeling_tokens_present()
        checks_more = [
            ("measurement context schema", context_pass, f"rows={len(context_rows)}"),
            ("eICU event fields remain unavailable", no_fake_eicu_events, f"eICU context rows={len(eicu_context)}"),
            ("no inferential fitting tokens", not no_modeling_tokens, ", ".join(no_modeling_tokens) or "none"),
            ("estimand inventory present", len(inventory_rows) >= 10, f"rows={len(inventory_rows)}"),
        ]
        for name, passed, detail in checks_more:
            check = {"name": name, "status": "PASS" if passed else "FAIL", "detail": detail}
            checks.append(check)
            if not passed:
                failures.append(f"{name}: {detail}")
    else:
        failures.append("Required generated context or estimand inventory is missing")
        checks.append({"name": "generated outputs present", "status": "FAIL", "detail": f"context={context_path.exists()}, inventory={inventory_path.exists()}"})

    historical = _historical_reference(config)
    checks.append({"name": "historical extension unchanged", "status": historical["status"], "detail": historical["root"]})
    if historical["status"] != "PASS":
        failures.append("Historical extension hash differs from the recorded Phase 1.5 reference")

    decision = "GO_PHASE2_MEASUREMENT_CONTEXT_ANALYSES" if not failures else "NO_GO_PHASE2_FINAL_LINEAGE_QC_FAIL"
    report_path = _write_report(config, decision, checks, reconstruction, inventory_rows, locked or {key: 0.0 for key in ("mimic_model_b_hr", "mimic_model_c_hr", "poct_hr", "lab_hr", "eicu_m3_rr", "eicu_m4_rr", "inspire_24h_hr", "inspire_48h_hr")}, historical)
    output_files = [
        report_path,
        output_root(config) / "machine_readable" / "final_measurement_context.csv",
        output_root(config) / "machine_readable" / "phase2_estimand_inventory.csv",
        output_root(config) / "qc" / "final_reconstruction_summary.json",
        output_root(config) / "qc" / "mimic_reconstruction_comparison.csv",
    ]
    rebase_manifest = {
        "manifest_version": "phase1.6-final-lineage-rebase-manifest-v1",
        "generated_at_utc": utc_now_iso(),
        "workspace_root": str(workspace_root(config).resolve()),
        "decision": decision,
        "modeling_performed": False,
        "public_code_git_head": _git_head(resolve_config_path(config, "public_code_root")),
        "final_input_manifest": str((workspace_root(config) / "manifests" / "final_input_manifest.json").resolve()),
        "authoritative_roots": {key: str(resolve_config_path(config, key).resolve()) for key in ("public_code_root", "analysis_record_root", "controlled_package_root", "manuscript_provenance_root")},
        "locked_fingerprints": locked,
        "historical_extension_reference": historical,
        "output_sha256": {str(path.relative_to(workspace_root(config))): sha256_file(path) for path in output_files if path.exists()},
        "checks_failed": failures,
    }
    rebase_path = workspace_root(config) / "manifests" / "rebase_manifest.json"
    write_json(rebase_path, rebase_manifest)
    final_gate = {
        "generated_at_utc": utc_now_iso(),
        "decision": decision,
        "modeling_performed": False,
        "failure_count": len(failures),
        "failures": failures,
        "qc_report": str(report_path.resolve()),
        "rebase_manifest": str(rebase_path.resolve()),
    }
    write_json(output_root(config) / "qc" / "final_gate.json", final_gate)
    log_path = workspace_root(config) / "logs" / "03_final_lineage_qc.log"
    log_path.parent.mkdir(parents=True, exist_ok=True)
    log_path.write_text(json.dumps(final_gate, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps(final_gate, indent=2, ensure_ascii=False))
    return 0 if not failures else 2


if __name__ == "__main__":
    raise SystemExit(main())
