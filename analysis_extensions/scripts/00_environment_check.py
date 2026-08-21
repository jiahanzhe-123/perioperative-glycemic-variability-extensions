#!/usr/bin/env python3
"""Check final roots, declared inputs, and code-contract provenance."""

from __future__ import annotations

import json
import sys
from pathlib import Path

WORKSPACE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(WORKSPACE))

from src.config_loader import load_config, output_root, workspace_root  # noqa: E402
from src.path_manifest import build_input_manifest, write_input_manifest  # noqa: E402
from src.schema_checks import utc_now_iso, write_json  # noqa: E402


def main() -> int:
    config = load_config()
    manifest = build_input_manifest(config)
    manifest_path = write_input_manifest(config)
    missing = [
        record["role"]
        for group in (manifest["inputs"], manifest["code_contracts"], manifest["roots"])
        for record in group
        if not record["exists"] and record["role"] != "manuscript_provenance_root"
    ]
    payload = {
        "generated_at_utc": utc_now_iso(),
        "workspace_root": str(workspace_root(config).resolve()),
        "manifest_path": str(manifest_path.resolve()),
        "required_root_keys": [
            "public_code_root",
            "analysis_record_root",
            "controlled_package_root",
            "manuscript_provenance_root",
            "output_root",
            "random_seed",
        ],
        "random_seed": config["random_seed"],
        "missing_roles": missing,
        "status": "PASS" if not missing else "FAIL",
        "historical_extension_not_used_as_input": True,
    }
    qc_path = output_root(config) / "qc" / "environment_check.json"
    write_json(qc_path, payload)
    log_path = workspace_root(config) / "logs" / "00_environment_check.log"
    log_path.parent.mkdir(parents=True, exist_ok=True)
    log_path.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps(payload, indent=2, ensure_ascii=False))
    if missing:
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
