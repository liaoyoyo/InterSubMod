#!/usr/bin/env bash
# Trigger Routing Hook — PostToolUse on Edit|Write
# Reads tool_input from stdin (JSON), checks file_path, outputs context-aware guidance.
# Exit 0 always (advisory only, never blocks).

set -euo pipefail

# Read tool input from stdin
INPUT=$(cat)
FILE_PATH=$(echo "$INPUT" | jq -r '.tool_input.file_path // empty' 2>/dev/null)

# Skip if no file path
[ -z "$FILE_PATH" ] && exit 0

case "$FILE_PATH" in
    */src/core/*.cpp|*/src/core/*.hpp|*/include/core/*.hpp|*/include/core/*.h)
        echo "[Routing] C++ core modified — consider: /methodology-audit before changes, /cpp-change for implementation protocol"
        ;;
    */CMakeLists.txt)
        echo "[Routing] CMakeLists.txt modified — run: mkdir -p build && cd build && cmake .. -DCMAKE_BUILD_TYPE=Release && make -j\$(nproc)"
        ;;
    */scripts/analysis/*.py|*/research/*.py)
        echo "[Routing] Analysis script modified — confirm: evidence_ledger update needed?"
        ;;
    */.claude/CLAUDE.md)
        echo "[Routing] CLAUDE.md modified — confirm consistency with AGENTS.md"
        ;;
    */AGENTS.md)
        echo "[Routing] AGENTS.md modified — confirm consistency with .claude/CLAUDE.md"
        ;;
    */.claude/skills/*/SKILL.md)
        echo "[Routing] Skill modified — verify: description has USE WHEN trigger, < 50 lines"
        ;;
    */.claude/settings*.json)
        echo "[Routing] Settings modified — verify JSON syntax: jq . < file"
        ;;
    */docs/experiments/*.md)
        echo "[Routing] Experiment doc modified — verify: docs/experiments/INDEX.md updated?"
        ;;
    */research/autoresearch/hypothesis_queue.json)
        echo "[Routing] Hypothesis queue modified — verify: tombstone check passed, priority scores correct"
        ;;
    */research/autoresearch/evidence_ledger.jsonl)
        echo "[Routing] Evidence ledger modified — verify: all required fields present (conditions, conclusion_stability)"
        ;;
esac

exit 0
