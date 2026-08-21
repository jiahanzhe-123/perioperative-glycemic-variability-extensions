"""Input and output provenance records for the final-lineage extension."""

from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from .config_loader import output_root, resolve_config_path, workspace_root


INPUT_KEYS = (
    "mimic_analysis_base",
    "mimic_features_priority",
    "mimic_series_priority",
    "mimic_samepatient_source",
    "mimic_primary_results",
    "mimic_mice_results",
    "mimic_source_results",
    "eicu_aggregate_input",
    "eicu_locked_results",
    "inspire_base",
    "inspire_outcome_30d",
    "inspire_glucose_features_anchor",
    "inspire_cohort_48h",
    "inspire_comorbidity",
    "inspire_primary_results",
    "inspire_48h_results",
    "inspire_time_rule_qc",
    "inspire_analysis_manifest",
)

CODE_KEYS = (
    "public_primary_model_script",
    "public_source_builder_script",
    "public_source_model_script",
    "public_eicu_model_script",
    "controlled_inspire_runner",
)

ROOT_KEYS = (
    "public_code_root",
    "analysis_record_root",
    "controlled_package_root",
    "manuscript_provenance_root",
    "historical_extension_root",
)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _mtime(path: Path) -> str | None:
    if not path.exists() and not path.is_symlink():
        return None
    return datetime.fromtimestamp(path.stat().st_mtime, tz=timezone.utc).isoformat()


def path_record(path: Path, role: str, inferential_input: bool) -> dict[str, Any]:
    resolved = path.expanduser().resolve(strict=False)
    exists = resolved.exists()
    is_file = exists and resolved.is_file()
    record: dict[str, Any] = {
        "role": role,
        "absolute_path": str(resolved),
        "exists": exists,
        "kind": "file" if is_file else ("directory" if exists and resolved.is_dir() else "missing"),
        "size_bytes": resolved.stat().st_size if is_file else None,
        "modified_time_utc": _mtime(resolved),
        "sha256": sha256_file(resolved) if is_file else None,
        "sha256_status": "computed" if is_file else "not_applicable",
        "inferential_input": inferential_input,
        "status": "AVAILABLE" if is_file else ("AVAILABLE_DIRECTORY" if exists else "MISSING"),
    }
    return record


def build_input_manifest(config: dict[str, Any]) -> dict[str, Any]:
    inputs = [
        path_record(resolve_config_path(config, key), key, True)
        for key in INPUT_KEYS
    ]
    code = [
        path_record(resolve_config_path(config, key), key, False)
        for key in CODE_KEYS
    ]
    roots = [
        path_record(resolve_config_path(config, key), key, False)
        for key in ROOT_KEYS
    ]
    return {
        "manifest_version": "phase1.6-final-lineage-input-manifest-v1",
        "generated_at_utc": datetime.now(timezone.utc).replace(microsecond=0).isoformat(),
        "random_seed": config["random_seed"],
        "workspace_root": str(workspace_root(config).resolve()),
        "output_root": str(output_root(config).resolve()),
        "raw_input_policy": "authoritative roots remain in place; this workspace reads only declared final CSVs and locked result files",
        "historical_extension_policy": "historical extension is provenance-only and is not copied, imported, or modified",
        "roots": roots,
        "inputs": inputs,
        "code_contracts": code,
        "status": "READY" if all(record["exists"] for record in inputs + code + roots if record["role"] != "manuscript_provenance_root") else "BLOCKED_MISSING_INPUT",
    }


def write_input_manifest(config: dict[str, Any]) -> Path:
    path = workspace_root(config) / "manifests" / "final_input_manifest.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(build_input_manifest(config), indent=2, ensure_ascii=False, sort_keys=True) + "\n", encoding="utf-8")
    return path
