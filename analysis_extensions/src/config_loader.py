"""Dependency-free loader for the flat Phase 1.6 configuration."""

from __future__ import annotations

import re
from pathlib import Path
from typing import Any


REQUIRED_KEYS = (
    "public_code_root",
    "analysis_record_root",
    "controlled_package_root",
    "manuscript_provenance_root",
    "output_root",
    "random_seed",
)


def _parse_scalar(value: str) -> Any:
    value = value.strip()
    if (value.startswith("'") and value.endswith("'")) or (
        value.startswith('"') and value.endswith('"')
    ):
        return value[1:-1]
    lowered = value.lower()
    if lowered in {"true", "yes"}:
        return True
    if lowered in {"false", "no"}:
        return False
    if lowered in {"null", "none", "~"}:
        return None
    if re.fullmatch(r"[-+]?\d+", value):
        return int(value)
    if re.fullmatch(r"[-+]?(?:\d+\.\d*|\d*\.\d+)(?:[eE][-+]?\d+)?", value):
        return float(value)
    return value


def _parse_flat_yaml(path: Path) -> dict[str, Any]:
    values: dict[str, Any] = {}
    for line_number, raw_line in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
        line = raw_line.strip()
        if not line or line.startswith("#"):
            continue
        if raw_line[:1].isspace():
            raise ValueError(f"Nested YAML is not supported at {path}:{line_number}")
        if ":" not in line:
            raise ValueError(f"Expected key: value at {path}:{line_number}")
        key, raw_value = line.split(":", 1)
        key = key.strip()
        if not key:
            raise ValueError(f"Empty configuration key at {path}:{line_number}")
        values[key] = _parse_scalar(raw_value)
    return values


def load_config(config_path: str | Path | None = None) -> dict[str, Any]:
    config_file = (
        Path(config_path).expanduser().resolve()
        if config_path is not None
        else Path(__file__).resolve().parents[1] / "config.yaml"
    )
    if not config_file.exists():
        raise FileNotFoundError(f"Configuration file not found: {config_file}")
    values = _parse_flat_yaml(config_file)
    missing = [key for key in REQUIRED_KEYS if key not in values]
    if missing:
        raise ValueError(f"Configuration is missing required keys: {', '.join(missing)}")
    for key in REQUIRED_KEYS[:-2]:
        if not Path(str(values[key])).expanduser().is_absolute():
            raise ValueError(f"{key} must be an absolute path")
    if not isinstance(values["random_seed"], int):
        raise ValueError("random_seed must be an integer")
    values["_config_path"] = config_file
    values["_workspace_root"] = config_file.parent
    return values


def resolve_config_path(config: dict[str, Any], key: str) -> Path:
    if key not in config:
        raise KeyError(f"Unknown configuration key: {key}")
    value = Path(str(config[key])).expanduser()
    if value.is_absolute():
        return value
    return Path(config["_workspace_root"]) / value


def workspace_root(config: dict[str, Any]) -> Path:
    return Path(config["_workspace_root"])


def output_root(config: dict[str, Any]) -> Path:
    return resolve_config_path(config, "output_root")
