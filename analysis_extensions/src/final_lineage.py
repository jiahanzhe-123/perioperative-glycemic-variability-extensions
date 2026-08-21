"""Deterministic, non-inferential reconstruction for the final JAHA v5 lineage.

This module deliberately stops at cohort masks, measurement summaries, and
provenance. It never calls a statistical fitting routine and never writes to
an authoritative root.
"""

from __future__ import annotations

from collections import defaultdict
from pathlib import Path
from typing import Any, Iterable, Mapping

from .config_loader import output_root, resolve_config_path
from .schema_checks import (
    SchemaError,
    all_nonmissing,
    finite_mean,
    is_missing,
    iter_csv,
    nonmissing,
    parse_bool,
    parse_datetime,
    parse_float,
    quantile,
    read_csv,
    require_columns,
    sample_sd,
    utc_now_iso,
    write_csv,
    write_json,
)


MIMIC_COVARIATES = (
    "age_at_admission",
    "gender",
    "bmi",
    "diabetes",
    "charlson_without_diabetes",
    "procedure_cat6",
    "lactate_postop_first",
    "creat_postop_first",
    "sofa_24h",
)

CONTEXT_FIELDS = (
    "database",
    "patient_or_episode_id",
    "analysis_context",
    "source_definition",
    "gv_sd",
    "mean_glucose",
    "measurement_count",
    "measurement_span_hours",
    "outcome_type",
    "event_30d",
    "event_365d",
    "hospital_id",
    "agreement_eligible",
    "outcome_model_eligible",
    "lineage_version",
    "input_source",
    "metric_role",
    "event_hospital_post_landmark",
    "coverage_status",
    "patient_level_coverage_end_field",
)

MIMIC_COMPARISON_FIELDS = (
    "stay_id",
    "stored_count",
    "recomputed_count",
    "count_match",
    "stored_mean_glucose",
    "recomputed_mean_glucose",
    "mean_absolute_difference",
    "stored_gv_sd",
    "recomputed_gv_sd",
    "gv_absolute_difference",
    "stored_span_hours",
    "recomputed_span_hours",
    "span_absolute_difference",
    "span_match",
    "event_30d",
    "event_365d",
)


def _path(config: Mapping[str, Any], key: str) -> Path:
    return resolve_config_path(dict(config), key)


def _input_label(config: Mapping[str, Any], *keys: str) -> str:
    return " | ".join(f"{key}:{_path(config, key)}" for key in keys)


def _event_value(row: Mapping[str, object], key: str) -> int | None:
    value = parse_bool(row.get(key))
    return None if value is None else int(value)


def _assert_unique(rows: Iterable[Mapping[str, object]], key: str, label: str) -> dict[str, Mapping[str, object]]:
    indexed: dict[str, Mapping[str, object]] = {}
    for row in rows:
        value = str(row.get(key, "")).strip()
        if not value:
            raise SchemaError(f"{label} contains an empty {key}")
        if value in indexed:
            raise SchemaError(f"{label} contains duplicate {key}: {value}")
        indexed[value] = row
    return indexed


def _target_mimic_rows(config: Mapping[str, Any]) -> tuple[list[dict[str, str]], dict[str, dict[str, str]], dict[str, dict[str, str]]]:
    base_path = _path(config, "mimic_analysis_base")
    feature_path = _path(config, "mimic_features_priority")
    series_path = _path(config, "mimic_series_priority")
    require_columns(
        base_path,
        ("stay_id", "landmark_eligible", "t_lm_30", "t_lm_365", "event_lm_30", "event_lm_365", *MIMIC_COVARIATES),
        "MIMIC analysis base",
    )
    require_columns(feature_path, ("stay_id", "glucose_count", "mean_glucose", "gv_sd", "span_hours"), "MIMIC priority features")
    require_columns(series_path, ("stay_id", "minute", "value", "final_source_class"), "MIMIC priority series")
    base = [dict(row) for row in iter_csv(base_path)]
    features = _assert_unique(iter_csv(feature_path), "stay_id", "MIMIC priority features")
    base_index = _assert_unique(base, "stay_id", "MIMIC analysis base")
    target: list[dict[str, str]] = []
    for row in base:
        stay_id = str(row["stay_id"])
        feature = features.get(stay_id)
        if (
            parse_bool(row.get("landmark_eligible")) is True
            and feature is not None
            and parse_float(feature.get("gv_sd")) is not None
            and (parse_float(feature.get("glucose_count")) or 0) >= 2
        ):
            target.append(row)
    if len(target) != 10561:
        raise SchemaError(f"Final MIMIC target N={len(target)}; expected 10561")
    return target, {str(row["stay_id"]): row for row in target}, features


def reconstruct_mimic(config: Mapping[str, Any]) -> dict[str, Any]:
    target_rows, target_index, features = _target_mimic_rows(config)
    series_path = _path(config, "mimic_series_priority")
    grouped: dict[str, list[dict[str, Any]]] = defaultdict(list)
    seen_minute: set[tuple[str, str]] = set()
    target_ids = set(target_index)
    for row in iter_csv(series_path):
        stay_id = str(row.get("stay_id", "")).strip()
        if stay_id not in target_ids:
            continue
        minute = str(row.get("minute", "")).strip()
        key = (stay_id, minute)
        if key in seen_minute:
            raise SchemaError(f"MIMIC priority series has duplicate stay_id/minute: {key}")
        seen_minute.add(key)
        value = parse_float(row.get("value"))
        timestamp = parse_datetime(minute)
        if value is None or timestamp is None:
            raise SchemaError(f"MIMIC priority series has a non-numeric value/timestamp for {key}")
        grouped[stay_id].append({"value": value, "timestamp": timestamp})

    context_rows: list[dict[str, Any]] = []
    comparisons: list[dict[str, Any]] = []
    counts: list[float] = []
    means: list[float] = []
    gvs: list[float] = []
    spans: list[float] = []
    event30 = 0
    event365 = 0
    for row in target_rows:
        stay_id = str(row["stay_id"])
        points = grouped.get(stay_id, [])
        if len(points) < 2:
            raise SchemaError(f"Final MIMIC target stay {stay_id} has fewer than two priority-series points")
        values = [float(point["value"]) for point in points]
        timestamps = [point["timestamp"] for point in points]
        rebuilt_count = len(values)
        rebuilt_mean = finite_mean(values)
        rebuilt_gv = sample_sd(values)
        rebuilt_span = (max(timestamps) - min(timestamps)).total_seconds() / 3600.0
        feature = features[stay_id]
        stored_count = parse_float(feature.get("glucose_count"))
        stored_mean = parse_float(feature.get("mean_glucose"))
        stored_gv = parse_float(feature.get("gv_sd"))
        stored_span = parse_float(feature.get("span_hours"))
        count_match = stored_count == rebuilt_count
        mean_diff = abs(float(stored_mean) - float(rebuilt_mean)) if stored_mean is not None and rebuilt_mean is not None else None
        gv_diff = abs(float(stored_gv) - float(rebuilt_gv)) if stored_gv is not None and rebuilt_gv is not None else None
        span_diff = abs(float(stored_span) - rebuilt_span) if stored_span is not None else None
        span_match = span_diff is not None and span_diff <= 1e-8
        if not count_match or mean_diff is None or mean_diff > 1e-8 or gv_diff is None or gv_diff > 1e-8 or not span_match:
            raise SchemaError(f"MIMIC stored/recomputed measurement mismatch for stay_id={stay_id}")
        e30 = _event_value(row, "event_lm_30")
        e365 = _event_value(row, "event_lm_365")
        event30 += int(e30 or 0)
        event365 += int(e365 or 0)
        counts.append(float(rebuilt_count))
        means.append(float(rebuilt_mean))
        gvs.append(float(rebuilt_gv))
        spans.append(float(rebuilt_span))
        comparisons.append(
            {
                "stay_id": stay_id,
                "stored_count": stored_count,
                "recomputed_count": rebuilt_count,
                "count_match": count_match,
                "stored_mean_glucose": stored_mean,
                "recomputed_mean_glucose": rebuilt_mean,
                "mean_absolute_difference": mean_diff,
                "stored_gv_sd": stored_gv,
                "recomputed_gv_sd": rebuilt_gv,
                "gv_absolute_difference": gv_diff,
                "stored_span_hours": stored_span,
                "recomputed_span_hours": rebuilt_span,
                "span_absolute_difference": span_diff,
                "span_match": span_match,
                "event_30d": e30,
                "event_365d": e365,
            }
        )
        context_rows.append(
            {
                "database": "MIMIC-IV",
                "patient_or_episode_id": stay_id,
                "analysis_context": "MIMIC_primary_24h_day1_landmark",
                "source_definition": "final source-priority series: central_lab > blood_gas > POCT > ICU_charted",
                "gv_sd": rebuilt_gv,
                "mean_glucose": rebuilt_mean,
                "measurement_count": rebuilt_count,
                "measurement_span_hours": rebuilt_span,
                "outcome_type": "all-cause mortality after day-1 landmark",
                "event_30d": e30,
                "event_365d": e365,
                "hospital_id": None,
                "agreement_eligible": None,
                "outcome_model_eligible": True,
                "lineage_version": "JAHA_v5_final_lineage",
                "input_source": _input_label(config, "mimic_analysis_base", "mimic_features_priority", "mimic_series_priority"),
                "metric_role": "primary_measurement_context",
                "event_hospital_post_landmark": None,
                "coverage_status": None,
                "patient_level_coverage_end_field": None,
            }
        )
    if event30 != 296 or event365 != 745:
        raise SchemaError(f"Final MIMIC events are {event30}/{event365}; expected 296/745")
    stored_gv = [float(row["stored_gv_sd"]) for row in comparisons]
    stored_mean = [float(row["stored_mean_glucose"]) for row in comparisons]
    stored_span = [float(row["stored_span_hours"]) for row in comparisons]
    return {
        "context_rows": context_rows,
        "comparisons": comparisons,
        "target_n": len(target_rows),
        "events_30d": event30,
        "events_365d": event365,
        "count_q1": quantile(counts, 0.25),
        "count_median": quantile(counts, 0.5),
        "count_q3": quantile(counts, 0.75),
        "span_q1": quantile(spans, 0.25),
        "span_median": quantile(spans, 0.5),
        "span_q3": quantile(spans, 0.75),
        "gv_mean": finite_mean(gvs),
        "mean_glucose_mean": finite_mean(means),
        "max_mean_absolute_difference": max(float(row["mean_absolute_difference"]) for row in comparisons),
        "max_gv_absolute_difference": max(float(row["gv_absolute_difference"]) for row in comparisons),
        "max_span_absolute_difference": max(float(row["span_absolute_difference"]) for row in comparisons),
        "max_stored_recomputed_span_ratio_difference": max(
            abs(float(row["recomputed_span_hours"]) / float(row["stored_span_hours"]) - 1.0)
            for row in comparisons
        ),
        "stored_gv_mean": finite_mean(stored_gv),
        "stored_mean_glucose_mean": finite_mean(stored_mean),
        "stored_span_median": quantile(stored_span, 0.5),
    }


def reconstruct_samepatient(config: Mapping[str, Any], mimic: Mapping[str, Any]) -> dict[str, Any]:
    path = _path(config, "mimic_samepatient_source")
    require_columns(path, ("stay_id", "gv_poct", "mean_glu_poct", "n_src_poct", "gv_lab", "mean_glu_lab", "n_src_lab"), "MIMIC same-patient source")
    source_rows = [dict(row) for row in iter_csv(path)]
    base = _assert_unique(iter_csv(_path(config, "mimic_analysis_base")), "stay_id", "MIMIC analysis base")
    target_ids = {str(row["patient_or_episode_id"]) for row in mimic["context_rows"]}
    agreement_n = 0
    agreement_target_n = 0
    model_rows: list[dict[str, str]] = []
    context_rows: list[dict[str, Any]] = []
    for row in source_rows:
        stay_id = str(row["stay_id"])
        pair_rule = (
            (parse_float(row.get("n_src_poct")) or 0) >= 2
            and (parse_float(row.get("n_src_lab")) or 0) >= 2
            and parse_float(row.get("gv_poct")) is not None
            and parse_float(row.get("gv_lab")) is not None
        )
        if pair_rule:
            agreement_n += 1
            if stay_id in target_ids:
                agreement_target_n += 1
        base_row = base.get(stay_id)
        outcome_model_eligible = bool(
            stay_id in target_ids
            and base_row is not None
            and all_nonmissing(base_row, ("t_lm_30", "t_lm_365", "event_lm_30", "event_lm_365", *MIMIC_COVARIATES))
        )
        if outcome_model_eligible:
            model_rows.append(row)
        e30 = _event_value(base_row or {}, "event_lm_30")
        e365 = _event_value(base_row or {}, "event_lm_365")
        for source_name, gv_key, mean_key, n_key in (
            ("POCT", "gv_poct", "mean_glu_poct", "n_src_poct"),
            ("central_laboratory", "gv_lab", "mean_glu_lab", "n_src_lab"),
        ):
            context_rows.append(
                {
                    "database": "MIMIC-IV",
                    "patient_or_episode_id": stay_id,
                    "analysis_context": "MIMIC_samepatient_paired_source_agreement",
                    "source_definition": source_name,
                    "gv_sd": parse_float(row.get(gv_key)),
                    "mean_glucose": parse_float(row.get(mean_key)),
                    "measurement_count": parse_float(row.get(n_key)),
                    "measurement_span_hours": None,
                    "outcome_type": "all-cause mortality after day-1 landmark",
                    "event_30d": e30,
                    "event_365d": e365,
                    "hospital_id": None,
                    "agreement_eligible": pair_rule,
                    "outcome_model_eligible": outcome_model_eligible,
                    "lineage_version": "JAHA_v5_final_lineage",
                    "input_source": _input_label(config, "mimic_samepatient_source", "mimic_analysis_base"),
                    "metric_role": f"samepatient_{source_name.lower()}_measurement_context",
                    "event_hospital_post_landmark": None,
                    "coverage_status": None,
                    "patient_level_coverage_end_field": None,
                }
            )
    model_events = sum(int(_event_value(base.get(str(row["stay_id"]), {}), "event_lm_30") or 0) for row in model_rows)
    model_events_365 = sum(int(_event_value(base.get(str(row["stay_id"]), {}), "event_lm_365") or 0) for row in model_rows)
    if len(source_rows) != 453 or agreement_n != 453 or len(model_rows) != 409 or model_events != 49:
        raise SchemaError(
            f"Same-patient counts candidate/agreement/model/events={len(source_rows)}/{agreement_n}/{len(model_rows)}/{model_events}; expected 453/453/409/49"
        )
    return {
        "context_rows": context_rows,
        "candidate_n": len(source_rows),
        "agreement_eligible_n": agreement_n,
        "agreement_eligible_final_mimic_n": agreement_target_n,
        "outcome_model_n": len(model_rows),
        "outcome_model_events_30d": model_events,
        "outcome_model_events_365d": model_events_365,
    }


def _eicu_outcome(row: Mapping[str, object]) -> int:
    return int(parse_bool(row.get("hosp_mortality")) is True and parse_bool(row.get("died_within_24h")) is not True)


def _eicu_complete(row: Mapping[str, object], columns: tuple[str, ...]) -> bool:
    categorical = {"sex", "procedure_category", "diabetes"}
    if not all_nonmissing(row, ("hosp_mortality", "died_within_24h")):
        return False
    if not all_nonmissing(row, columns):
        return False
    numeric = [column for column in columns if column not in categorical]
    return all(parse_float(row.get(column)) is not None for column in numeric)


def reconstruct_eicu(config: Mapping[str, Any]) -> dict[str, Any]:
    path = _path(config, "eicu_aggregate_input")
    require_columns(
        path,
        (
            "patientunitstayid", "hospitalid", "procedure_category", "age", "sex", "diabetes", "creatinine", "apachescore",
            "glucose_n", "glucose_mean_24h", "glucose_sd_24h", "measurement_span_minutes", "hosp_mortality", "died_within_24h", "in_main", "in_landmark",
        ),
        "eICU aggregate input",
    )
    rows = [dict(row) for row in iter_csv(path)]
    main = [row for row in rows if parse_bool(row.get("in_main")) is True and parse_bool(row.get("in_landmark")) is True]
    m3_columns = ("glucose_sd_24h", "age", "sex", "procedure_category", "diabetes", "creatinine", "apachescore", "glucose_mean_24h")
    m4_columns = m3_columns + ("glucose_n", "measurement_span_minutes")
    m3 = [row for row in main if _eicu_complete(row, m3_columns)]
    m4 = [row for row in main if _eicu_complete(row, m4_columns)]
    if len(m3) != 7115 or len(m4) != 7115 or sum(_eicu_outcome(row) for row in m3) != 130 or sum(_eicu_outcome(row) for row in m4) != 130:
        raise SchemaError(f"eICU masks are M3={len(m3)}/{sum(_eicu_outcome(row) for row in m3)}, M4={len(m4)}/{sum(_eicu_outcome(row) for row in m4)}; expected 7115/130")
    context_rows: list[dict[str, Any]] = []
    for context_name, selected in (("eICU_harmonized_M3", m3), ("eICU_harmonized_M4", m4)):
        for row in selected:
            context_rows.append(
                {
                    "database": "eICU-CRD",
                    "patient_or_episode_id": str(row["patientunitstayid"]),
                    "analysis_context": context_name,
                    "source_definition": "final aggregate all-source glucose features; no event-level eICU glucose file preserved",
                    "gv_sd": parse_float(row.get("glucose_sd_24h")),
                    "mean_glucose": parse_float(row.get("glucose_mean_24h")),
                    "measurement_count": parse_float(row.get("glucose_n")),
                    "measurement_span_hours": (parse_float(row.get("measurement_span_minutes")) or 0.0) / 60.0,
                    "outcome_type": "post-landmark hospital mortality",
                    "event_30d": None,
                    "event_365d": None,
                    "hospital_id": row.get("hospitalid"),
                    "agreement_eligible": None,
                    "outcome_model_eligible": True,
                    "lineage_version": "JAHA_v5_final_lineage",
                    "input_source": _input_label(config, "eicu_aggregate_input"),
                    "metric_role": context_name,
                    "event_hospital_post_landmark": _eicu_outcome(row),
                    "coverage_status": None,
                    "patient_level_coverage_end_field": None,
                }
            )
    return {
        "context_rows": context_rows,
        "aggregate_input_n": len(rows),
        "main_landmark_n": len(main),
        "main_landmark_events": sum(_eicu_outcome(row) for row in main),
        "m3_n": len(m3),
        "m3_events": sum(_eicu_outcome(row) for row in m3),
        "m4_n": len(m4),
        "m4_events": sum(_eicu_outcome(row) for row in m4),
        "event_level_input_available": False,
    }


def _inspire_event(row: Mapping[str, object], landmark_key: str) -> int:
    death = parse_float(row.get("death_time_composite"))
    landmark = parse_float(row.get(landmark_key))
    cutoff = parse_float(row.get("day30_cutoff"))
    return int(death is not None and landmark is not None and cutoff is not None and death > landmark and death <= cutoff)


def _inspire_post_landmark(row: Mapping[str, object], landmark_key: str) -> bool:
    death = parse_float(row.get("death_time_composite"))
    landmark = parse_float(row.get(landmark_key))
    return landmark is not None and (death is None or death > landmark)


def reconstruct_inspire(config: Mapping[str, Any]) -> dict[str, Any]:
    base_path = _path(config, "inspire_base")
    outcome_path = _path(config, "inspire_outcome_30d")
    anchor_path = _path(config, "inspire_glucose_features_anchor")
    cohort48_path = _path(config, "inspire_cohort_48h")
    comorbidity_path = _path(config, "inspire_comorbidity")
    require_columns(base_path, ("subject_id", "op_id", "age", "mean_glucose", "gv_sd", "n_glucose_0_24h", "span_hours"), "INSPIRE base")
    require_columns(outcome_path, ("subject_id", "op_id", "landmark_24h", "day30_cutoff", "death_time_composite", "event_30d_reconciled"), "INSPIRE reconciled outcome")
    require_columns(anchor_path, ("subject_id", "anchor", "n", "mean_glucose", "gv_sd", "span_hours"), "INSPIRE anchor features")
    require_columns(cohort48_path, ("subject_id", "landmark_48h", "day30_cutoff", "death_time_composite", "age", "mean_glucose", "gv_sd", "n_glucose_0_48h", "span_hours"), "INSPIRE 48h cohort")
    require_columns(comorbidity_path, ("subject_id", "charlson_without_diabetes", "diabetes"), "INSPIRE comorbidity")
    base = _assert_unique(iter_csv(base_path), "subject_id", "INSPIRE base")
    base_by_key = {(str(row["subject_id"]), str(row["op_id"])): row for row in base.values()}
    anchor_rows = [dict(row) for row in iter_csv(anchor_path)]
    if not anchor_rows:
        raise SchemaError("INSPIRE anchor feature input is empty")
    outcome_rows = [dict(row) for row in iter_csv(outcome_path)]
    d24: list[dict[str, str]] = []
    for row in outcome_rows:
        key = (str(row["subject_id"]), str(row["op_id"]))
        base_row = base_by_key.get(key)
        if base_row is None:
            continue
        merged = {
            "subject_id": key[0],
            "op_id": key[1],
            "landmark_24h": row.get("landmark_24h", ""),
            "day30_cutoff": row.get("day30_cutoff", ""),
            "death_time_composite": row.get("death_time_composite", ""),
            "event_30d_reconciled": row.get("event_30d_reconciled", ""),
            **base_row,
        }
        if _inspire_post_landmark(merged, "landmark_24h"):
            d24.append(merged)
    required24 = ("landmark_24h", "day30_cutoff", "age", "gv_sd", "mean_glucose", "n_glucose_0_24h", "span_hours")
    if len(d24) != 1353 or not all(all_nonmissing(row, required24) for row in d24) or sum(_inspire_event(row, "landmark_24h") for row in d24) != 27:
        raise SchemaError(f"INSPIRE 24h frame is N={len(d24)}, events={sum(_inspire_event(row, 'landmark_24h') for row in d24)}; expected 1353/27")
    d48 = [dict(row) for row in iter_csv(cohort48_path) if _inspire_post_landmark(row, "landmark_48h")]
    required48 = ("landmark_48h", "day30_cutoff", "age", "gv_sd", "mean_glucose", "n_glucose_0_48h", "span_hours")
    if len(d48) != 1511 or not all(all_nonmissing(row, required48) for row in d48) or sum(_inspire_event(row, "landmark_48h") for row in d48) != 31:
        raise SchemaError(f"INSPIRE 48h frame is N={len(d48)}, events={sum(_inspire_event(row, 'landmark_48h') for row in d48)}; expected 1511/31")
    coverage_rows = read_csv(_path(config, "inspire_time_rule_qc"))
    coverage = next((row for row in coverage_rows if row.get("check") == "coverage_assumption_recorded"), None)
    if coverage is None or coverage.get("status") != "ATTESTED":
        raise SchemaError("INSPIRE day-30 coverage attestation is missing or not ATTESTED")
    coverage_status = "SOURCE_PROVENANCE_ATTESTED_TO_DAY30"
    context_rows: list[dict[str, Any]] = []
    for context_name, selected, landmark_key, count_key in (
        ("INSPIRE_24h_landmark", d24, "landmark_24h", "n_glucose_0_24h"),
        ("INSPIRE_48h_landmark", d48, "landmark_48h", "n_glucose_0_48h"),
    ):
        for row in selected:
            context_rows.append(
                {
                    "database": "INSPIRE",
                    "patient_or_episode_id": f"{row['subject_id']}:{row.get('op_id', 'NA')}",
                    "analysis_context": context_name,
                    "source_definition": "INSPIRE final anchor glucose features with uniform administrative day-30 censoring",
                    "gv_sd": parse_float(row.get("gv_sd")),
                    "mean_glucose": parse_float(row.get("mean_glucose")),
                    "measurement_count": parse_float(row.get(count_key)),
                    "measurement_span_hours": parse_float(row.get("span_hours")),
                    "outcome_type": "all-cause mortality by postoperative day 30",
                    "event_30d": _inspire_event(row, landmark_key),
                    "event_365d": None,
                    "hospital_id": None,
                    "agreement_eligible": None,
                    "outcome_model_eligible": True,
                    "lineage_version": "JAHA_v5_final_lineage",
                    "input_source": _input_label(config, "inspire_base", "inspire_outcome_30d", "inspire_cohort_48h", "inspire_time_rule_qc"),
                    "metric_role": context_name,
                    "event_hospital_post_landmark": None,
                    "coverage_status": coverage_status,
                    "patient_level_coverage_end_field": False,
                }
            )
    return {
        "context_rows": context_rows,
        "24h_n": len(d24),
        "24h_events": sum(_inspire_event(row, "landmark_24h") for row in d24),
        "48h_n": len(d48),
        "48h_events": sum(_inspire_event(row, "landmark_48h") for row in d48),
        "coverage_status": coverage_status,
        "patient_level_coverage_end_field": False,
        "anchor_input_n": len(anchor_rows),
    }


def reconstruct_all(config: Mapping[str, Any]) -> dict[str, Any]:
    mimic = reconstruct_mimic(config)
    samepatient = reconstruct_samepatient(config, mimic)
    eicu = reconstruct_eicu(config)
    inspire = reconstruct_inspire(config)
    return {"mimic": mimic, "samepatient": samepatient, "eicu": eicu, "inspire": inspire}


def write_context_outputs(config: Mapping[str, Any], reconstruction: Mapping[str, Any]) -> dict[str, Any]:
    output = output_root(dict(config))
    machine = output / "machine_readable"
    qc = output / "qc"
    context_rows: list[dict[str, Any]] = []
    for key in ("mimic", "samepatient", "eicu", "inspire"):
        context_rows.extend(reconstruction[key]["context_rows"])
    write_csv(machine / "final_measurement_context.csv", context_rows, CONTEXT_FIELDS)
    write_csv(qc / "mimic_reconstruction_comparison.csv", reconstruction["mimic"]["comparisons"], MIMIC_COMPARISON_FIELDS)
    summary = {
        "generated_at_utc": utc_now_iso(),
        "phase": "Phase 1.6",
        "modeling_performed": False,
        "random_seed": config["random_seed"],
        "mimic": {key: value for key, value in reconstruction["mimic"].items() if key not in {"context_rows", "comparisons"}},
        "samepatient": {key: value for key, value in reconstruction["samepatient"].items() if key != "context_rows"},
        "eicu": {key: value for key, value in reconstruction["eicu"].items() if key != "context_rows"},
        "inspire": {key: value for key, value in reconstruction["inspire"].items() if key != "context_rows"},
        "context_row_count": len(context_rows),
        "outputs": {
            "final_measurement_context": str((machine / "final_measurement_context.csv").resolve()),
            "mimic_reconstruction_comparison": str((qc / "mimic_reconstruction_comparison.csv").resolve()),
            "summary": str((qc / "final_reconstruction_summary.json").resolve()),
        },
    }
    write_json(qc / "final_reconstruction_summary.json", summary)
    return summary
