#!/usr/bin/env python3
# sensitive_scan.py — 仓库敏感信息扫描(凭证/患者 ID 字面量/绝对路径/大文件)
import os, re, sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SKIP_DIRS = {".git", "data", "results", "__pycache__", ".ipynb_checkpoints", ".pytest_cache", "node_modules", "archive/.git"}
SKIP_EXTS = (".pyc", ".pyo", ".png", ".pdf", ".rds", ".pkl", ".zip", ".gz", ".RData", ".rda")
PATTERNS = [
    ("credential_password", re.compile(r"(?i)(password|passwd|pwd|token|api[_-]?key|secret)\s*[:=]\s*['\"]?[^\s'\"]{4,}")),
    ("pg_password_literal", re.compile(r"(?i)(postgres|pg)[^\n]{0,40}(password|pwd)\s*[:=]\s*['\"]?[^\s'\"]{4,}")),
    ("abs_user_path", re.compile(r"/Users/[A-Za-z0-9_]+/")),
    ("home_path", re.compile(r"/home/[A-Za-z0-9_]+/")),
    ("volume_path", re.compile(r"/Volumes/[^/]+/")),
    ("server_path", re.compile(r"(?i)(srv|server|nas|internal)[A-Za-z0-9\.\-]*[:/]{2}")),
    ("ipv4", re.compile(r"\b(?:10|172\.(?:1[6-9]|2\d|3[01])|192\.168)\.\d{1,3}\.\d{1,3}\b")),
    ("ssh_key", re.compile(r"-----BEGIN (RSA|OPENSSH|EC|DSA) PRIVATE KEY-----")),
]
PLACEHOLDER_OK = re.compile(r"(PGV\(|/path/to/|your-path|example|placeholder|authorized|re\.compile)")
findings = []
big = []
for dp, dns, fns in os.walk(ROOT):
    dns[:] = [d for d in dns if d not in SKIP_DIRS]
    for fn in fns:
        if fn.endswith(SKIP_EXTS) or fn == "sensitive_scan.py":
            continue
        # config/paths.yml 是本机授权路径的指定存放处(已 gitignore),不扫描
        if fn == "paths.yml" and os.path.basename(dp) == "config":
            continue
        p = os.path.join(dp, fn)
        rp = os.path.relpath(p, ROOT)
        try:
            sz = os.path.getsize(p)
            if sz > 2_000_000:
                big.append((rp, sz))
            if sz > 500_000:
                continue
            with open(p, encoding="utf-8", errors="ignore") as f:
                for i, line in enumerate(f, 1):
                    for name, pat in PATTERNS:
                        if pat.search(line) and not PLACEHOLDER_OK.search(line):
                            findings.append((rp, i, name, line.strip()[:120]))
        except Exception:
            continue
print("=== sensitive scan ===")
if findings:
    for f_ in findings[:40]:
        print("FIND:", f_)
else:
    print("no credential/path/ID literals found in tracked source files")
if big:
    print("large files (>2MB):")
    for b in big: print("  ", b)
ok = len(findings) == 0
print("scan result:", "PASS" if ok else "FAIL")
sys.exit(0 if ok else 1)
