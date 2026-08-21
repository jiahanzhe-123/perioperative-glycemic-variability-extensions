"""Small standard-library schema, parsing, and output helpers."""

from __future__ import annotations

import csv
import json
import math
from datetime import datetime, timezone
from pathlib import Path
from statistics import mean, median, stdev
from typing import Any, Iterable, Iterator, Mapping, Sequence


class SchemaError(RuntimeError):
    """Raised when a final-lineage input violates its declared contract."""


def iter_csv(path: Path) -> Iterator[dict[str, str]]:
    with path.open(newline="", encoding="utf-8", errors="replace") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None:
            raise SchemaError(f"CSV has no header: {path}")
        for row in reader:
            yield row


def read_csv(path: Path) -> list[dict[str, str]]:
    return list(iter_csv(path))


def read_header(path: Path) -> list[str]:
    with path.open(newline="", encoding="utf-8", errors="replace") as handle:
        reader = csv.reader(handle)
        try:
            header = next(reader)
        except StopIteration as exc:
            raise SchemaError(f"CSV is empty: {path}") from exc
    if not header or any(not str(column).strip() for column in header):
        raise SchemaError(f"CSV has an invalid header: {path}")
    return [str(column) for column in header]


def require_columns(path: Path, required: Sequence[str], label: str | None = None) -> None:
    available = set(read_header(path))
    missing = [column for column in required if column not in available]
    if missing:
        raise SchemaError(f"{label or path} is missing required columns: {', '.join(missing)}")


def is_missing(value: object) -> bool:
    return value is None or str(value).strip().lower() in {"", "na", "nan", "null", "none"}


def parse_float(value: object) -> float | None:
    if is_missing(value):
        return None
    try:
        number = float(str(value).strip())
    except ValueError as exc:
        raise SchemaError(f"Expected numeric value, received {value!r}") from exc
    return number if math.isfinite(number) else None


def parse_bool(value: object) -> bool | None:
    if is_missing(value):
        return None
    text = str(value).strip().lower()
    if text in {"true", "1", "yes", "y", "t"}:
        return True
    if text in {"false", "0", "no", "n", "f"}:
        return False
    raise SchemaError(f"Unrecognized boolean value: {value!r}")


def parse_datetime(value: object) -> datetime | None:
    if is_missing(value):
        return None
    text = str(value).strip().replace("Z", "+00:00")
    try:
        return datetime.fromisoformat(text)
    except ValueError as exc:
        raise SchemaError(f"Invalid ISO-like timestamp: {value!r}") from exc


def nonmissing(value: object) -> bool:
    return not is_missing(value)


def all_nonmissing(row: Mapping[str, object], columns: Sequence[str]) -> bool:
    return all(nonmissing(row.get(column)) for column in columns)


def sample_sd(values: Sequence[float]) -> float | None:
    return float(stdev(values)) if len(values) >= 2 else None


def finite_mean(values: Sequence[float]) -> float | None:
    return float(mean(values)) if values else None


def finite_median(values: Sequence[float]) -> float | None:
    return float(median(values)) if values else None


def quantile(values: Sequence[float], probability: float) -> float:
    if not values:
        raise ValueError("Cannot calculate a quantile for an empty sequence")
    ordered = sorted(float(value) for value in values)
    position = (len(ordered) - 1) * probability
    lower = math.floor(position)
    upper = math.ceil(position)
    if lower == upper:
        return ordered[lower]
    return ordered[lower] + (ordered[upper] - ordered[lower]) * (position - lower)


def utc_now_iso() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


def format_value(value: Any) -> str:
    if value is None:
        return "NA"
    if isinstance(value, bool):
        return "TRUE" if value else "FALSE"
    if isinstance(value, float):
        return "NA" if not math.isfinite(value) else format(value, ".15g")
    return str(value)


def write_csv(path: Path, rows: Iterable[Mapping[str, Any]], fieldnames: Sequence[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(fieldnames), extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({key: format_value(row.get(key)) for key in fieldnames})


def write_json(path: Path, payload: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, ensure_ascii=False, sort_keys=True) + "\n", encoding="utf-8")
