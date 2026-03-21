#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "${ROOT_DIR}"

PASS_COUNT=0
WARN_COUNT=0
FAIL_COUNT=0

ok() {
    echo "[OK] $*"
    PASS_COUNT=$((PASS_COUNT + 1))
}

warn() {
    echo "[WARN] $*"
    WARN_COUNT=$((WARN_COUNT + 1))
}

fail() {
    echo "[FAIL] $*"
    FAIL_COUNT=$((FAIL_COUNT + 1))
}

echo "=== AI Agent Readiness Check ==="
echo "Repo: ${ROOT_DIR}"
echo "Time: $(date '+%Y-%m-%d %H:%M:%S')"
echo

# 1) Core entry files
for f in README.md QUICKSTART.md docs/README.md docs/CURRENT_FOCUS.md docs/standards/README.md AGENTS.md .mcp.json; do
    if [[ -f "${f}" ]]; then
        ok "Found file: ${f}"
    else
        fail "Missing file: ${f}"
    fi
done

echo

# 2) MCP knowledge server config
if python3 - <<'PY'
import json
from pathlib import Path
cfg = json.loads(Path(".mcp.json").read_text())
srv = cfg.get("mcpServers", {}).get("knowledge", {})
cmd = srv.get("command")
args = srv.get("args", [])
arg0 = args[0].lower() if len(args) >= 1 else ""
ok = cmd == "python3" and len(args) >= 1 and arg0.endswith("/knowledge/scripts/mcp/knowledge_server.py")
raise SystemExit(0 if ok else 1)
PY
then
    ok ".mcp.json knowledge server configuration looks correct"
else
    fail ".mcp.json knowledge server is missing or path mismatch"
fi

echo

# 3) output symlink policy
if [[ -L output ]]; then
    target="$(readlink output || true)"
    if [[ "${target}" == "/big8_disk/liaoyoyo2001/InterSubMod_runs/output" ]]; then
        ok "output symlink target is canonical: ${target}"
    else
        warn "output symlink target is non-canonical: ${target}"
    fi
else
    fail "output is not a symlink"
fi

echo

# 4) docs structural validation
if [[ -f scripts/analysis/validate_docs_structure.py ]]; then
    tmp="$(mktemp)"
    python3 scripts/analysis/validate_docs_structure.py > "${tmp}" || true

    total_md="$(awk -F': ' '/^Total markdown files:/ {print $2}' "${tmp}")"
    name_issues="$(awk -F': ' '/^Name issues:/ {print $2}' "${tmp}")"
    metadata_issues="$(awk -F': ' '/^Metadata issues:/ {print $2}' "${tmp}")"
    link_issues="$(awk -F': ' '/^Link issues:/ {print $2}' "${tmp}")"

    if [[ "${name_issues}" == "0" ]]; then
        ok "Name issues = 0"
    else
        fail "Name issues = ${name_issues}"
    fi

    if [[ "${metadata_issues}" == "0" ]]; then
        ok "Metadata issues = 0"
    else
        fail "Metadata issues = ${metadata_issues}"
    fi

    if [[ "${link_issues}" == "0" ]]; then
        ok "Link issues = 0"
    else
        warn "Link issues = ${link_issues} (expected if only archive/deep immutable)"
    fi

    ok "Total markdown files = ${total_md}"
    rm -f "${tmp}"
else
    fail "Missing scripts/analysis/validate_docs_structure.py"
fi

echo

# 5) Legacy-path signals that confuse agents
legacy_patterns=(
    "experiments/outputs"
    "output/bip8_disk_output"
    "docs/reports/2026"
    "docs/analysis/"
)

scan_paths=(
    "README.md"
    "QUICKSTART.md"
    "README_PROJECT_SUMMARY.md"
    "docs/README.md"
    "docs/CURRENT_FOCUS.md"
    "docs/standards"
    "scripts"
)

for p in "${legacy_patterns[@]}"; do
    count="$(
        (
            rg -n "${p}" "${scan_paths[@]}" -S \
                --glob '!scripts/analysis/check_ai_agent_readiness.sh' \
                --glob '!docs/references/manual/20260301_AI_Agent_快速操作手冊_01.md' \
                --glob '!docs/reports/validated/2026/03/20260301_AI_Agent可操作性盤點與改善方案_01.md' \
                || true
        ) | wc -l | tr -d ' '
    )"
    if [[ "${count}" == "0" ]]; then
        ok "Legacy pattern '${p}' = 0"
    else
        warn "Legacy pattern '${p}' occurrences = ${count}"
    fi
done

echo

# 6) Empty directories (potentially confusing)
empty_count="$(find docs -maxdepth 4 -type d -empty | wc -l | tr -d ' ')"
if [[ "${empty_count}" == "0" ]]; then
    ok "No empty directories under docs (maxdepth=4)"
else
    warn "Empty directories under docs (maxdepth=4): ${empty_count}"
fi

echo
echo "=== Summary ==="
echo "PASS: ${PASS_COUNT}"
echo "WARN: ${WARN_COUNT}"
echo "FAIL: ${FAIL_COUNT}"

if [[ "${FAIL_COUNT}" -gt 0 ]]; then
    exit 1
fi

exit 0
