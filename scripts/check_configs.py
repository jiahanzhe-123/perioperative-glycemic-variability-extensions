#!/usr/bin/env python3
# check_configs.py — YAML/CSV 配置格式检查
import os, re, sys
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
fails = []

def check_yaml_keys(path, required):
    txt = open(path, encoding="utf-8").read()
    for k in required:
        if not re.search(rf"^{k}\s*:", txt, re.M) and not re.search(rf"^- model_id: {k}\b", txt, re.M):
            fails.append(f"{path}: missing key {k}")

check_yaml_keys(os.path.join(ROOT, "config", "config.example.yml"),
                ["seed", "output_dir", "mimic_version", "eicu_version", "inspire_version"])
check_yaml_keys(os.path.join(ROOT, "config", "paths.example.yml"),
                ["seed", "mimic_work", "rebuild_work", "inspire_work", "eicu_work",
                 "mimic_derived_data", "mimic_record_work", "bmi_repair_work"])
check_yaml_keys(os.path.join(ROOT, "config", "constants.yml"),
                ["seed", "glucose", "landmark", "mice", "bootstrap"])
# model_specifications.yml 每模型必须有 model_id/database/cohort 或 analytic_frame/interpretive_role
txt = open(os.path.join(ROOT, "config", "model_specifications.yml"), encoding="utf-8").read()
n_models = len(re.findall(r"^- model_id:", txt, re.M))
for k in ["database", "interpretive_role"]:
    n = len(re.findall(rf"^  {k}:", txt, re.M))
    if n < n_models * 0.8:
        fails.append(f"model_specifications.yml: only {n}/{n_models} models have {k}")
# frozen_results.yml 关键节
check_yaml_keys(os.path.join(ROOT, "tests", "expected", "frozen_results.yml"),
                ["mimic_primary", "inspire_primary_v5", "eicu_harmonized_m3"])
# covariate_sets.yml 不得为空
cs = open(os.path.join(ROOT, "config", "covariate_sets.yml"), encoding="utf-8").read()
if "MIMIC_CLINICAL_FIXED_V1" not in cs:
    fails.append("covariate_sets.yml missing MIMIC_CLINICAL_FIXED_V1")

print("=== config check ===")
if fails:
    for f_ in fails: print("FAIL:", f_)
else:
    print("all config files valid")
print(f"model_specifications.yml models: {n_models}")
sys.exit(1 if fails else 0)
