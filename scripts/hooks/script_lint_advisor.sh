#!/usr/bin/env bash
# script_lint_advisor.sh — PostToolUse(Edit|Write) advisory lint for .py / .sh under scripts/
#
# WHY: C++ has a compile Hard Gate (pre_commit_compile_check.sh) but .py/.sh had ZERO syntax
# check. The harness's own governance helpers (task_graph.py, number_provenance.py,
# harness_health.py, _build_full_tree.py) ship unlinted — a parse error in an "auditor" silently
# produces wrong output (insights friction #4 root cause + builder-vs-helper hazard).
#
# DESIGN: ADVISORY only (exit 0 always) + fail-OPEN. Lint != test; never blocks legitimate WIP.
# Scoped to scripts/ to avoid noise on research analysis scratch files. py_compile for .py,
# bash -n for .sh (+ shellcheck if available). Mirrors the advisory-exit-0 pattern of
# provenance_stamp_advisor.sh / health_drift_advisor.sh.
set -u

# Read hook stdin JSON; extract edited file path (tool_input.file_path). Fail-open on any parse issue.
payload="$(cat 2>/dev/null || true)"
fp=""
if command -v python3 >/dev/null 2>&1; then
    fp="$(printf '%s' "$payload" | python3 -c 'import sys,json
try:
    d=json.load(sys.stdin); ti=d.get("tool_input",{})
    print(ti.get("file_path") or ti.get("path") or "")
except Exception:
    print("")' 2>/dev/null || true)"
fi
[ -z "$fp" ] && exit 0

# Scope: only scripts/ (the governance-helper + pipeline home). Skip elsewhere.
case "$fp" in
    */scripts/*|scripts/*) : ;;
    *) exit 0 ;;
esac
[ -f "$fp" ] || exit 0

case "$fp" in
    *.py)
        if command -v python3 >/dev/null 2>&1; then
            err="$(python3 -m py_compile "$fp" 2>&1)" || {
                echo "[script-lint] ⚠ py_compile FAILED — $fp" 1>&2
                echo "$err" | head -5 1>&2
                echo "[script-lint] advisory only（未阻擋）；治理 helper 解析錯會默默產出錯誤結果，請修正後再跑。" 1>&2
            }
        fi
        ;;
    *.sh)
        if command -v bash >/dev/null 2>&1; then
            err="$(bash -n "$fp" 2>&1)" || {
                echo "[script-lint] ⚠ bash -n syntax error — $fp" 1>&2
                echo "$err" | head -5 1>&2
            }
        fi
        if command -v shellcheck >/dev/null 2>&1; then
            shellcheck -S error "$fp" 1>&2 2>/dev/null || true
        fi
        ;;
esac
exit 0
