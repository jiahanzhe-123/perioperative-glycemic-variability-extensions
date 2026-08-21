# test_core_rules.py — 核心规则单元测试(不依赖真实患者数据)
import numpy as np
import pandas as pd

PRIO = {"central_lab": 1, "blood_gas": 2, "poct": 3, "icu_charted": 4}

def mgdl_to_mmoll(x):
    return x / 18.018

def minute_median(recs):
    return recs.groupby(["stay", "minute"])["value"].median()

def select_priority(df):
    df = df.copy()
    df["pr"] = df["source"].map(PRIO)
    best = df.groupby(["stay", "minute"])["pr"].transform("min")
    return df[df["pr"] == best]

def gv_sd(vals):
    vals = np.asarray(vals, dtype=float)
    return vals.std(ddof=1) if len(vals) >= 2 else np.nan

def eAG(a1c):
    return 28.7 * a1c - 46.7

def in_window(t, t0, span=1440):
    return (t >= t0) & (t < t0 + span)

def bmi_plausible(bmi):
    return (bmi >= 10) & (bmi <= 80)

def interval_label(days):
    return np.where(days <= 7, "1-7", np.where(days <= 30, "8-30", "31-365"))

def test_mmoll_conversion():
    assert abs(mgdl_to_mmoll(10) - 0.555) < 0.001

def test_right_open_window():
    t0 = 100
    assert in_window(np.array([99]), t0).tolist() == [False]
    assert in_window(np.array([100]), t0).tolist() == [True]
    assert in_window(np.array([100 + 1439]), t0).tolist() == [True]
    assert in_window(np.array([100 + 1440]), t0).tolist() == [False]

def test_exact_duplicate_removal():
    df = pd.DataFrame({"stay": [1, 1, 1, 2], "minute": [5, 5, 6, 7],
                       "value": [100.0, 100.0, 110.0, 120.0],
                       "source": ["poct", "poct", "poct", "central_lab"]})
    dd = df.drop_duplicates(["stay", "minute", "value", "source"])
    assert len(dd) == 3

def test_source_priority_per_minute():
    df = pd.DataFrame({"stay": [1]*4, "minute": [5, 5, 5, 6],
                       "value": [100, 102, 104, 120],
                       "source": ["icu_charted", "poct", "blood_gas", "poct"]})
    sel = select_priority(df)
    assert (sel["source"].tolist() == ["blood_gas", "poct"])

def test_same_minute_median():
    df = pd.DataFrame({"stay": [1]*3, "minute": [5, 5, 5], "value": [100, 110, 400]})
    m = minute_median(df)
    assert m.iloc[0] == 110

def test_gv_requires_two_timepoints():
    assert np.isnan(gv_sd([120]))
    assert gv_sd([100, 120]) > 0

def test_gv_and_mean_same_set():
    vals = [90, 110, 130]
    assert abs(np.mean(vals) - 110) < 1e-9
    assert gv_sd(vals) == np.std(vals, ddof=1)

def test_shr_and_eag_positive():
    eag = eAG(6.0)
    assert abs(eag - 125.5) < 1e-6
    assert eag > 0
    assert abs((140 / eag) - 140/125.5) < 1e-9

def test_hba1c_strict_window():
    days = np.array([-91, -90, -45, -1, 0, 1])
    ok = (days < 0) & (days >= -90)
    assert ok.tolist() == [False, True, True, True, False, False]

def test_bmi_plausibility_rule():
    x = np.array([9.9, 10, 45, 80, 80.1, np.nan])
    out = bmi_plausible(x)
    assert out[~np.isnan(x)].tolist() == [False, True, True, True, False]

def test_post_landmark_interval_labels():
    assert interval_label(np.array([1, 7, 8, 30, 31, 365])).tolist() == \
        ["1-7", "1-7", "8-30", "8-30", "31-365", "31-365"]

def test_per10_equals_half_mmoll():
    assert abs(10/18.018 - 0.555) < 0.001
