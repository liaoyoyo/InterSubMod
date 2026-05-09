#!/usr/bin/env bash
# tmpdir_guard.sh — Forced TMPDIR export + pre-flight free-space check
#
# Source: prior art D-3 (Snakemake `resources.tmpdir`) + issue #1059 (wrapper may not honor TMPDIR)
# Plan ref: ~/.claude/plans/agent-harness-langgraph-resilient-waterfall.md §12.2 T1-2
# Lessons: feedback_tmp_disk_full_pipeline_pitfall.md (longphase-to wrote BAM to /tmp → 800GB disaster)
#
# Usage (source from any pipeline launcher; do NOT execute):
#   source /big7_disk/liaoyoyo2001/InterSubMod/scripts/lib/tmpdir_guard.sh
#
# Effect:
#   1. export TMPDIR to per-process directory under $TMPDIR_BASE (default /big7_disk/liaoyoyo2001/tmp)
#   2. mkdir -p the path
#   3. trap EXIT cleanup
#   4. Pre-flight: abort if $TMPDIR_BASE has <10GB free
#   5. echo TMPDIR to stderr (caller can verify wrapper actually consumes it)
#
# Tools known to silently ignore TMPDIR (must verify per Snakemake #1059):
#   - longphase-to (use -o explicit output dir)
#   - some samtools subcommands (use -T <prefix>)
#   - clairs / clairs-to (use --output_prefix)

# Do NOT use 'set -e' here — sourcing scripts must not abort caller on minor issues.

# --- Config ---
: "${TMPDIR_BASE:=/big7_disk/liaoyoyo2001/tmp}"
: "${TMPDIR_MIN_FREE_KB:=10000000}"   # 10 GB
: "${TMPDIR_GUARD_KEEP_ON_EXIT:=0}"   # set 1 to skip cleanup (debugging)

# --- Build per-process tmpdir ---
TMPDIR_BASE="$TMPDIR_BASE"  # in case of typo
USER_NAME="$(whoami)"
PROC_TMPDIR="${TMPDIR_BASE}/${USER_NAME}/$$"

if ! mkdir -p "$PROC_TMPDIR" 2>/dev/null; then
    echo "[tmpdir_guard] ERROR: cannot create $PROC_TMPDIR" >&2
    return 1 2>/dev/null || exit 1
fi

# --- Pre-flight free-space check ---
AVAIL_KB="$(df --output=avail "$TMPDIR_BASE" 2>/dev/null | tail -1 | tr -d ' ')"
if [[ "$AVAIL_KB" =~ ^[0-9]+$ ]] && [ "$AVAIL_KB" -lt "$TMPDIR_MIN_FREE_KB" ]; then
    echo "[tmpdir_guard] ABORT: $TMPDIR_BASE has only ${AVAIL_KB}KB free (<${TMPDIR_MIN_FREE_KB}KB threshold)" >&2
    return 1 2>/dev/null || exit 1
fi

# --- Export to env ---
export TMPDIR="$PROC_TMPDIR"
export TMP="$PROC_TMPDIR"
export TEMP="$PROC_TMPDIR"

# --- Cleanup trap ---
if [ "$TMPDIR_GUARD_KEEP_ON_EXIT" != "1" ]; then
    trap 'rm -rf "$TMPDIR" 2>/dev/null' EXIT
fi

# --- Audit line to stderr (lets caller verify) ---
echo "[tmpdir_guard] TMPDIR=$TMPDIR free=$((AVAIL_KB/1024))MB" >&2
