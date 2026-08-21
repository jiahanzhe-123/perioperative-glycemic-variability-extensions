"""PGV 配置访问(Python)。

规则:
- 外部数据键(患者级/授权目录)必须来自 config/paths.yml(本机,gitignored);
  模板见 config/paths.example.yml。未配置时明确报错,绝不回退到仓库内目录
  (正式代码不依赖仓库内 symlink 或默认数据目录)。
- 仓库内部键(synthetic/home)有仓库相对默认值,无需配置。
"""
import os
import re

_REPO_INTERNAL = {
    "synthetic": os.path.join("data", "synthetic"),
    "home": ".",
}

def _repo_root(start=None):
    d = os.path.abspath(start or os.getcwd())
    while True:
        if os.path.exists(os.path.join(d, "config", "paths.example.yml")):
            return d
        nd = os.path.dirname(d)
        if nd == d:
            raise RuntimeError(f"repository root not found above {start}")
        d = nd

def _read_yml(path):
    out = {}
    if not os.path.exists(path):
        return out
    for ln in open(path, encoding="utf-8"):
        ln = re.sub(r"#.*$", "", ln).rstrip()
        m = re.match(r"^([A-Za-z0-9_]+)\s*:\s*(.*)$", ln)
        if m:
            out[m.group(1)] = m.group(2).strip().strip('"').strip("'")
    return out

_ROOT = _repo_root()
def _is_placeholder(v):
    return (not v) or v.startswith("/path/to/")

_CFG_USER = {k: os.path.expanduser(v) for k, v in
             _read_yml(os.path.join(_ROOT, "config", "paths.yml")).items()
             if not _is_placeholder(v)}

def PGV(key):
    if key in _CFG_USER:
        return _CFG_USER[key]
    if key in _REPO_INTERNAL:
        return os.path.join(_ROOT, _REPO_INTERNAL[key])
    raise KeyError(
        f"PGV: key '{key}' is not configured. External data paths must be set in "
        f"config/paths.yml (copy config/paths.example.yml and edit). "
        f"Formal code never falls back to in-repo data directories.")

def PGV_OUT(sub="results"):
    return os.path.join(_ROOT, sub)

PGV_SEED = 20260726
