#!/usr/bin/env bash
# Standalone HTML Trigger Hook — PostToolUse on Edit|Write
# Detects when LLM writes/edits a markdown report in paths that warrant a
# standalone HTML companion, then prints a reminder for Claude to invoke
# the `html-report-build` skill in `standalone` mode.
#
# Tier-A path triggers (validated reports, concluded experiments, etc.):
#   docs/reports/validated/{YYYY}/{MM}/*.md
#   docs/reports/research_landscape/*.md
#   docs/reports/pi_reports/*.md
#   docs/experiments/concluded/{YYYY}/{MM}/*.md
#   docs/concepts/*.md
#   docs/data_specs/*.md
#   docs/references/manual/*.md
#
# SKIP paths (LLM-consumption / source-of-truth / state):
#   CLAUDE.md / AGENTS.md / MEMORY.md / README.md
#   state/*.json / *.yaml / *.csv / *.tsv
#   docs/experiments/in_progress/*.md  (drafts; use report mode if needed)
#
# Exit 0 always (advisory only, never blocks).

set -euo pipefail

INPUT=$(cat)
FILE_PATH=$(echo "$INPUT" | jq -r '.tool_input.file_path // empty' 2>/dev/null)

# Skip if no file path
[ -z "$FILE_PATH" ] && exit 0

# Only act on .md files
case "$FILE_PATH" in
    *.md) ;;
    *) exit 0 ;;
esac

# Skip LLM-consumption / source-of-truth docs
case "$FILE_PATH" in
    */CLAUDE.md|*/AGENTS.md|*/MEMORY.md|*/README.md)
        exit 0
        ;;
    */docs/experiments/in_progress/*)
        # Drafts — companion not warranted yet
        exit 0
        ;;
    */hypothesis_queue.json|*/evidence_ledger.jsonl|*/state.json|*/active.json)
        exit 0
        ;;
esac

# Check companion HTML existence — if .standalone.html already exists & is newer, skip
BASENAME="${FILE_PATH%.md}"
STANDALONE_HTML="${BASENAME}.standalone.html"
if [ -f "$STANDALONE_HTML" ] && [ "$STANDALONE_HTML" -nt "$FILE_PATH" ]; then
    # Standalone newer than .md → already up-to-date
    exit 0
fi

# Path-based mode recommendation
case "$FILE_PATH" in
    */docs/reports/validated/*)
        echo "[html-report-build] Validated report modified: $(basename "$FILE_PATH")"
        echo "[html-report-build] → 建議跑 'html-report-build' skill in standalone mode"
        echo "[html-report-build] → Output: ${STANDALONE_HTML#*/InterSubMod/}"
        ;;
    */docs/reports/research_landscape/*)
        echo "[html-report-build] Research landscape doc modified: $(basename "$FILE_PATH")"
        echo "[html-report-build] → 建議跑 'html-report-build' skill in standalone mode (主軸大圖景，PI 反覆讀)"
        ;;
    */docs/reports/pi_reports/*)
        echo "[html-report-build] PI report modified: $(basename "$FILE_PATH")"
        echo "[html-report-build] → 建議跑 'html-report-build' skill in standalone mode (給 PI 終版)"
        ;;
    */docs/experiments/concluded/*)
        echo "[html-report-build] Concluded experiment modified: $(basename "$FILE_PATH")"
        echo "[html-report-build] → 建議跑 'html-report-build' skill in standalone mode (tier ⭐4-5 收尾)"
        ;;
    */docs/concepts/*)
        echo "[html-report-build] Concept doc modified: $(basename "$FILE_PATH")"
        echo "[html-report-build] → 建議跑 'html-report-build' skill in standalone mode (技術說明文章，sticky TOC 受益大)"
        ;;
    */docs/data_specs/*)
        echo "[html-report-build] Data spec modified: $(basename "$FILE_PATH")"
        echo "[html-report-build] → 建議跑 'html-report-build' skill in standalone mode (規格類技術文章)"
        ;;
    */docs/references/manual/*)
        echo "[html-report-build] Manual / handbook modified: $(basename "$FILE_PATH")"
        echo "[html-report-build] → 建議跑 'html-report-build' skill in standalone mode (技術 handbook)"
        ;;
    */docs/presentations/*/02_slide_outline.md|*/docs/presentations/*/03_slide_layout_script.md)
        echo "[html-report-build] Slide outline modified: $(basename "$FILE_PATH")"
        echo "[html-report-build] → 建議跑 'html-report-build' skill in **slide** mode"
        ;;
    *)
        # Other .md paths — no specific recommendation, exit silently
        exit 0
        ;;
esac

exit 0
