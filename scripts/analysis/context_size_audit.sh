#!/usr/bin/env bash
# context_size_audit.sh
# Measure InterSubMod's "always loaded" Claude Code context size.
# Reports CLAUDE.md + 6 mandatory research files + N landscape files,
# with line/word/estimated-token counts for each + totals.
#
# Usage:
#   bash context_size_audit.sh                 # stdout report
#   bash context_size_audit.sh --snapshot LBL  # also writes JSON to .claude/state/
#
# Env overrides:
#   REPO_ROOT             explicit repo root (else `git rev-parse --show-toplevel`)
#   CLAUDE_MD_OVERRIDE    explicit path to CLAUDE.md (for testing missing-file behavior)
#
# Token estimate: words * 1.3 (integer math: words * 13 / 10).
# Real BPE tokens may differ by ~20%.

set -euo pipefail

# --- Argument parsing ---------------------------------------------------------
SNAPSHOT_LABEL=""
while [ $# -gt 0 ]; do
    case "$1" in
        --snapshot)
            SNAPSHOT_LABEL="${2:-}"
            if [ -z "$SNAPSHOT_LABEL" ]; then
                echo "ERROR: --snapshot requires a LABEL argument" >&2
                exit 1
            fi
            shift 2
            ;;
        -h|--help)
            sed -n '2,16p' "$0"
            exit 0
            ;;
        *)
            echo "ERROR: unknown argument: $1" >&2
            exit 1
            ;;
    esac
done

# --- Repo root resolution -----------------------------------------------------
if [ -z "${REPO_ROOT:-}" ]; then
    REPO_ROOT="$(git rev-parse --show-toplevel 2>/dev/null || pwd)"
fi

CLAUDE_MD_PATH="${CLAUDE_MD_OVERRIDE:-${REPO_ROOT}/.claude/CLAUDE.md}"

# --- Mandatory file list (hardcoded — resilient to CLAUDE.md format drift) ---
MANDATORY_FILES=(
    "docs/CURRENT_FOCUS.md"
    "docs/experiments/INDEX.md"
    "docs/README.md"
    "docs/references/manual/20260424_AI啟動壓縮上下文與研究索引_01.md"
    "docs/concepts/2026/04/20260409_研究構想總索引_01.md"
    "docs/reports/research_landscape/00_INDEX.md"
)

# Landscape files: discover dynamically (resilient to additions/removals)
LANDSCAPE_DIR="${REPO_ROOT}/docs/reports/research_landscape"
mapfile -t LANDSCAPE_FILES < <(
    if [ -d "$LANDSCAPE_DIR" ]; then
        find "$LANDSCAPE_DIR" -maxdepth 1 -type f -name "*.md" \
            ! -name "00_INDEX.md" | sort
    fi
)

# --- Helpers ------------------------------------------------------------------
# count_file <abs_path>  → echoes "lines words tokens"
# missing → "MISSING MISSING MISSING"
count_file() {
    local f="$1"
    if [ -f "$f" ]; then
        local lines words
        lines=$(wc -l < "$f" | tr -d ' ')
        words=$(wc -w < "$f" | tr -d ' ')
        local tokens=$(( words * 13 / 10 ))
        echo "$lines $words $tokens"
    else
        echo "MISSING MISSING MISSING"
    fi
}

# truncate string to N chars with trailing "..." if over
truncate() {
    local s="$1" n="$2"
    if [ "${#s}" -gt "$n" ]; then
        echo "${s:0:$((n - 3))}..."
    else
        echo "$s"
    fi
}

# --- Header -------------------------------------------------------------------
TIMESTAMP="$(date -u +%Y-%m-%dT%H:%M:%SZ)"
echo "=== InterSubMod Always-Loaded Context Audit ==="
echo "Repo root: $REPO_ROOT"
echo "Audit time: $TIMESTAMP"
echo

# --- CLAUDE.md ----------------------------------------------------------------
echo "--- CLAUDE.md ---"
read -r CMD_LINES CMD_WORDS CMD_TOKENS <<< "$(count_file "$CLAUDE_MD_PATH")"
if [ "$CMD_LINES" = "MISSING" ]; then
    echo "  MISSING: $CLAUDE_MD_PATH"
    CMD_TOKENS_NUM=0
    CMD_LINES_NUM=0
    CMD_WORDS_NUM=0
else
    echo "Lines: $CMD_LINES  Words: $CMD_WORDS  Est. tokens: $CMD_TOKENS"
    CMD_TOKENS_NUM=$CMD_TOKENS
    CMD_LINES_NUM=$CMD_LINES
    CMD_WORDS_NUM=$CMD_WORDS
fi
echo

# --- Mandatory files ----------------------------------------------------------
echo "--- ${#MANDATORY_FILES[@]} mandatory files ---"
printf "%-55s %8s %8s %8s\n" "File" "Lines" "Words" "Tokens"

MAND_LINES_TOTAL=0
MAND_WORDS_TOTAL=0
MAND_TOKENS_TOTAL=0
MAND_JSON_ITEMS=()

for rel in "${MANDATORY_FILES[@]}"; do
    abs="${REPO_ROOT}/${rel}"
    read -r L W T <<< "$(count_file "$abs")"
    disp=$(truncate "$rel" 55)
    if [ "$L" = "MISSING" ]; then
        printf "%-55s %8s %8s %8s\n" "$disp" "MISSING" "-" "-"
        MAND_JSON_ITEMS+=("{\"path\":\"$rel\",\"missing\":true}")
    else
        printf "%-55s %8d %8d %8d\n" "$disp" "$L" "$W" "$T"
        MAND_LINES_TOTAL=$(( MAND_LINES_TOTAL + L ))
        MAND_WORDS_TOTAL=$(( MAND_WORDS_TOTAL + W ))
        MAND_TOKENS_TOTAL=$(( MAND_TOKENS_TOTAL + T ))
        MAND_JSON_ITEMS+=("{\"path\":\"$rel\",\"lines\":$L,\"words\":$W,\"tokens\":$T}")
    fi
done
printf "%-55s %8d %8d %8d\n" "SUBTOTAL" "$MAND_LINES_TOTAL" "$MAND_WORDS_TOTAL" "$MAND_TOKENS_TOTAL"
echo

# --- Landscape files ----------------------------------------------------------
echo "--- ${#LANDSCAPE_FILES[@]} research landscape files (excl 00_INDEX) ---"
printf "%-55s %8s %8s %8s\n" "File" "Lines" "Words" "Tokens"

LAND_LINES_TOTAL=0
LAND_WORDS_TOTAL=0
LAND_TOKENS_TOTAL=0
LAND_JSON_ITEMS=()

for abs in "${LANDSCAPE_FILES[@]}"; do
    rel="${abs#${REPO_ROOT}/}"
    read -r L W T <<< "$(count_file "$abs")"
    disp=$(truncate "$rel" 55)
    if [ "$L" = "MISSING" ]; then
        printf "%-55s %8s %8s %8s\n" "$disp" "MISSING" "-" "-"
        LAND_JSON_ITEMS+=("{\"path\":\"$rel\",\"missing\":true}")
    else
        printf "%-55s %8d %8d %8d\n" "$disp" "$L" "$W" "$T"
        LAND_LINES_TOTAL=$(( LAND_LINES_TOTAL + L ))
        LAND_WORDS_TOTAL=$(( LAND_WORDS_TOTAL + W ))
        LAND_TOKENS_TOTAL=$(( LAND_TOKENS_TOTAL + T ))
        LAND_JSON_ITEMS+=("{\"path\":\"$rel\",\"lines\":$L,\"words\":$W,\"tokens\":$T}")
    fi
done
printf "%-55s %8d %8d %8d\n" "SUBTOTAL" "$LAND_LINES_TOTAL" "$LAND_WORDS_TOTAL" "$LAND_TOKENS_TOTAL"
echo

# --- Totals -------------------------------------------------------------------
TOTAL_TOKENS=$(( CMD_TOKENS_NUM + MAND_TOKENS_TOTAL + LAND_TOKENS_TOTAL ))
TOTAL_LINES=$(( CMD_LINES_NUM + MAND_LINES_TOTAL + LAND_LINES_TOTAL ))
TOTAL_WORDS=$(( CMD_WORDS_NUM + MAND_WORDS_TOTAL + LAND_WORDS_TOTAL ))
BUDGET=300000
# Headroom (1 decimal): budget * 10 / total → e.g. 64 → "6.4x"
if [ "$TOTAL_TOKENS" -gt 0 ]; then
    HEADROOM_X10=$(( BUDGET * 10 / TOTAL_TOKENS ))
    HEADROOM="$(( HEADROOM_X10 / 10 )).$(( HEADROOM_X10 % 10 ))x"
else
    HEADROOM="n/a"
fi

echo "=== TOTAL ==="
echo "Always-loaded: ~${TOTAL_TOKENS} tokens"
echo "Recommended budget for context rot: <${BUDGET} tokens"
echo "Headroom: ${HEADROOM}"
echo
echo "[NOTE: Token estimates use words * 1.3 ratio. Real BPE tokens may differ +/- 20%.]"

# --- Snapshot JSON ------------------------------------------------------------
if [ -n "$SNAPSHOT_LABEL" ]; then
    STATE_DIR="${REPO_ROOT}/.claude/state"
    mkdir -p "$STATE_DIR"
    TS_FILE="$(date -u +%Y%m%dT%H%M%SZ)"
    JSON_PATH="${STATE_DIR}/context_audit_${SNAPSHOT_LABEL}_${TS_FILE}.json"

    # Build JSON with jq for robust escaping
    MAND_JSON="[$(IFS=,; echo "${MAND_JSON_ITEMS[*]}")]"
    LAND_JSON="[$(IFS=,; echo "${LAND_JSON_ITEMS[*]}")]"

    if [ "$CMD_LINES" = "MISSING" ]; then
        CLAUDE_MD_JSON="{\"path\":\"$CLAUDE_MD_PATH\",\"missing\":true}"
    else
        CLAUDE_MD_JSON="{\"path\":\"$CLAUDE_MD_PATH\",\"lines\":$CMD_LINES_NUM,\"words\":$CMD_WORDS_NUM,\"tokens\":$CMD_TOKENS_NUM}"
    fi

    jq -n \
        --arg ts "$TIMESTAMP" \
        --arg snap_label "$SNAPSHOT_LABEL" \
        --arg repo "$REPO_ROOT" \
        --argjson cmd "$CLAUDE_MD_JSON" \
        --argjson mand "$MAND_JSON" \
        --argjson land "$LAND_JSON" \
        --argjson tot_lines "$TOTAL_LINES" \
        --argjson tot_words "$TOTAL_WORDS" \
        --argjson tot_tokens "$TOTAL_TOKENS" \
        --argjson mand_tokens "$MAND_TOKENS_TOTAL" \
        --argjson land_tokens "$LAND_TOKENS_TOTAL" \
        --argjson budget "$BUDGET" \
        '{
            timestamp: $ts,
            label: $snap_label,
            repo_root: $repo,
            claude_md: $cmd,
            mandatory_files: $mand,
            landscape_files: $land,
            totals: {
                lines: $tot_lines,
                words: $tot_words,
                tokens: $tot_tokens,
                mandatory_tokens: $mand_tokens,
                landscape_tokens: $land_tokens,
                budget: $budget
            }
        }' > "$JSON_PATH"

    echo
    echo "Snapshot written: $JSON_PATH"
fi
