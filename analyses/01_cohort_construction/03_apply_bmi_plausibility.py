#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""03_apply_bmi_plausibility.py — BMI 合理性规则 BMI-DQ-1(analysis of record 数据准备步骤)。

规则(2026-07-30 在任何结局模型重跑前冻结;见 docs/RESULT_PROVENANCE.md):
    bmi 为有限值且 10 <= bmi <= 80 kg/m2  -> 保留原值
    bmi < 10 或 > 80 或非有限值           -> 设为缺失(NA),交冻结 MICE 插补
    原本缺失                              -> 保持缺失
仅作用于 analysis_base.csv 的 bmi 列(MIMIC-IV 侧);height_cm/weight_kg 不进模型,不修复。

输入 : PGV("mimic_derived_data")/data/analysis_base.csv        (只读授权数据)
输出 : PGV("mimic_record_work")/data/analysis_base_bmi_repaired.csv
QC   : 聚合 manifest + checksums 写入 PGV("mimic_record_work")/qc/;
        患者级 change log 写入 PGV("mimic_record_work")/qc_private/(本地受限,绝不入库)。

本脚本是 phase3_5 BMI 修复工程 _work/_task_d.py 的配置化迁移(逻辑一致,
路径改由 config/paths.yml 驱动),并额外校验输出与冻结 checksum 的一致性。
"""
import os, sys, hashlib, json, datetime
import pandas as pd
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "src", "common"))
from paths import PGV

BMI_LO, BMI_HI = 10.0, 80.0           # BMI-DQ-1 主规则边界(冻结)
SENS_LO, SENS_HI = 12.0, 60.0         # 边界敏感性(只计数,不用于任何结局模型)
EXPECTED_ROWS = 12992                 # 冻结分析框行数
# 冻结的 analysis-of-record 派生框 checksum(phase3_5 BMI_REPAIR_CHECKSUMS.csv)
FROZEN_SHA256 = "2163168d1c2cf07828dd10d3b79392ca94545fbf36e7fd68dcd9d2675eeefba4"
FROZEN_CHANGED = 35                   # 冻结的规则命中行数(28 例 >80 + 7 例 <10)

def sha256(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()

def main():
    src = os.path.join(PGV("mimic_derived_data"), "data", "analysis_base.csv")
    work = PGV("mimic_record_work")
    d_data = os.path.join(work, "data")
    d_qc = os.path.join(work, "qc")
    d_priv = os.path.join(work, "qc_private")
    for d in (d_data, d_qc, d_priv):
        os.makedirs(d, exist_ok=True)
    dst = os.path.join(d_data, "analysis_base_bmi_repaired.csv")

    df = pd.read_csv(src, low_memory=False)
    assert df["stay_id"].is_unique, "stay_id not unique in source frame"
    assert len(df) == EXPECTED_ROWS, f"source frame rows {len(df)} != frozen {EXPECTED_ROWS}"

    bmi = df["bmi"]
    bad = bmi.notna() & ((bmi < BMI_LO) | (bmi > BMI_HI) | ~np.isfinite(bmi))
    bad_sens = bmi.notna() & ((bmi < SENS_LO) | (bmi > SENS_HI))
    extra_sens = int((bad_sens & ~bad).sum())
    n_bad = int(bad.sum())
    print(f"[BMI-DQ-1] rule [{BMI_LO:g},{BMI_HI:g}]: set-missing n = {n_bad}")
    print(f"[BMI-DQ-1] boundary sensitivity [{SENS_LO:g},{SENS_HI:g}]: extra flags = {extra_sens} (count only)")
    if n_bad != FROZEN_CHANGED:
        print(f"FATAL: rule hit {n_bad} rows != frozen {FROZEN_CHANGED}; aborting (no files written).")
        sys.exit(1)

    rep = df.copy()
    rep.loc[bad, "bmi"] = np.nan
    rep.to_csv(dst, index=False)

    # QC 1: 除 bmi 外所有列逐值一致
    chk = df.drop(columns=["bmi"]).equals(rep.drop(columns=["bmi"]))
    assert chk, "non-bmi columns differ after rule application"
    # QC 2: 行数/唯一性/缺失增量
    reread = pd.read_csv(dst, low_memory=False)
    assert len(reread) == len(df) and reread["stay_id"].is_unique
    assert int(reread["bmi"].isna().sum()) == int(df["bmi"].isna().sum()) + n_bad
    print("[QC] non-bmi columns identical:", chk, "; rows:", len(reread),
          "; newly missing:", n_bad)

    # QC 3: 与冻结 checksum 比对(证明复现 analysis-of-record 输入框)
    dst_sha = sha256(dst)
    sha_match = (dst_sha == FROZEN_SHA256)
    print("[QC] derived frame sha256:", dst_sha)
    print("[QC] frozen analysis-of-record sha256 match:", sha_match)
    if not sha_match:
        ref = os.path.join(PGV("bmi_repair_work"), "analysis_base_bmi_repaired.csv")
        if os.path.exists(ref):
            a = pd.read_csv(dst, low_memory=False)
            b = pd.read_csv(ref, low_memory=False)
            same = a.equals(b)
            cell_diff = int((a.fillna(-999).astype(str) != b.fillna(-999).astype(str)).sum().sum()) if not same else 0
            print(f"[QC] byte differs but semantic-equal vs frozen frame: {same} (cell diffs: {cell_diff})")
            if not same:
                print("FATAL: derived frame semantically differs from frozen analysis-of-record frame.")
                sys.exit(1)
        else:
            print("WARN: frozen reference frame unavailable; checksum match could not be confirmed.")

    # 患者级 change log(本地受限目录;绝不入库)
    log = df.loc[bad, ["subject_id", "hadm_id", "stay_id", "bmi"]].rename(columns={"bmi": "old_bmi"})
    log["new_bmi"] = "NA (set missing)"
    log["rule_id"] = "BMI-DQ-1"
    log["reason"] = np.where(log["old_bmi"] < BMI_LO, "bmi < 10 kg/m2 (implausible)",
                     np.where(log["old_bmi"] > BMI_HI, "bmi > 80 kg/m2 (implausible)", "nonfinite"))
    log.to_csv(os.path.join(d_priv, "BMI_REPAIR_CHANGE_LOG_PRIVATE.csv"), index=False)

    # 聚合 manifest + checksums(可入库格式的本地副本)
    man = pd.DataFrame([
        {"item": "rule_id", "value": "BMI-DQ-1"},
        {"item": "rule", "value": f"bmi<{BMI_LO:g} or >{BMI_HI:g} or nonfinite -> NA -> frozen MICE"},
        {"item": "source_sha256", "value": sha256(src)},
        {"item": "derived_sha256", "value": dst_sha},
        {"item": "frozen_sha256_match", "value": sha_match},
        {"item": "rows_total", "value": len(df)},
        {"item": "records_changed_set_missing", "value": n_bad},
        {"item": "records_originally_missing", "value": int(df["bmi"].isna().sum())},
        {"item": "non_bmi_columns_identical", "value": chk},
        {"item": "stay_id_unique", "value": bool(reread["stay_id"].is_unique)},
        {"item": "boundary_sensitivity_extra_flags_count_only", "value": extra_sens},
        {"item": "created", "value": datetime.datetime.now().isoformat(timespec="seconds")},
    ])
    man.to_csv(os.path.join(d_qc, "bmi_plausibility_manifest.csv"), index=False)
    print("BMI_PLAUSIBILITY_DONE", dst)

if __name__ == "__main__":
    main()
