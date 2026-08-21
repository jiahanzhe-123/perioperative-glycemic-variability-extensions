#!/usr/bin/env bash
# 全量入口(Level 1 可运行;Level 2/3 需授权数据与 config/paths.yml)
set -euo pipefail
cd "$(dirname "$0")/.."
echo "[pgv] lint"; make lint
echo "[pgv] tests"; make test
echo "[pgv] synthetic"; make synthetic
echo "[pgv] validate (frozen-config format)"; make validate
