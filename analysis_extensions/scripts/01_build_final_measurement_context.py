#!/usr/bin/env python3
"""Build final-lineage measurement context without inferential fitting."""

from __future__ import annotations

import json
import sys
from pathlib import Path

WORKSPACE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(WORKSPACE))

from src.config_loader import load_config, workspace_root  # noqa: E402
from src.final_lineage import reconstruct_all, write_context_outputs  # noqa: E402


def main() -> int:
    config = load_config()
    try:
        reconstruction = reconstruct_all(config)
        summary = write_context_outputs(config, reconstruction)
        status = {"status": "PASS", **summary}
    except Exception as exc:
        status = {"status": "FAIL", "error_type": type(exc).__name__, "error": str(exc), "modeling_performed": False}
        log_path = workspace_root(config) / "logs" / "01_build_final_measurement_context.log"
        log_path.parent.mkdir(parents=True, exist_ok=True)
        log_path.write_text(json.dumps(status, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
        print(json.dumps(status, indent=2, ensure_ascii=False))
        return 2
    log_path = workspace_root(config) / "logs" / "01_build_final_measurement_context.log"
    log_path.parent.mkdir(parents=True, exist_ok=True)
    log_path.write_text(json.dumps(status, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps(status, indent=2, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
