#!/usr/bin/env bash
# pipeline_block_check.sh — PreToolUse hook: hard-block long pipelines when disk is critical
#
# Source: prior art D-1+D-2 integration
# Plan ref: ~/.claude/plans/agent-harness-langgraph-resilient-waterfall.md §12.2 T1-3
#
# Behavior:
#   - Read .PIPELINE_BLOCKED flag (created by disk_guard.sh when worst mount ≥95%)
#   - If present AND tool is invoking a known long-running pipeline script → exit 2 (block)
#   - Otherwise pass through
#
# Triggers (matchers configured in .claude/settings.local.json):
#   - Bash invoking scripts/run_*.sh
#   - Bash invoking longphase / clairs / samtools sort with high-volume input
#
# Override (when user accepts risk):
#   rm $HOME/.intersubmod_pipeline_blocked    # explicit removal
#   OR
#   PIPELINE_BLOCK_OVERRIDE=1 <command>       # one-shot bypass

set -u

BLOCK_FLAG="${DISK_GUARD_BLOCK_FLAG:-$HOME/.intersubmod_pipeline_blocked}"

# Read tool input from stdin (Claude Code hook protocol)
TOOL_INPUT="$(cat 2>/dev/null || true)"

# Override env: user explicitly accepted risk
if [ "${PIPELINE_BLOCK_OVERRIDE:-0}" = "1" ]; then
    exit 0
fi

# No flag → no action
if [ ! -f "$BLOCK_FLAG" ]; then
    exit 0
fi

# Flag present: check if command looks like a long pipeline
COMMAND="$(echo "$TOOL_INPUT" | grep -oE '"command"\s*:\s*"[^"]*"' | head -1 || true)"
if echo "$COMMAND" | grep -qE 'run_(vcf|batch|random)|longphase|clairs|samtools.*sort|build_master|run_pipeline'; then
    REASON="$(cat "${BLOCK_FLAG}.reason" 2>/dev/null || echo 'disk usage critical')"
    cat <<EOF >&2
[pipeline_block_check] BLOCKED: long pipeline execution rejected.
  Flag: $BLOCK_FLAG
  Reason: $REASON
  Resolve options:
    1. Free disk space (recommended) — flag auto-clears when disk drops below 95%
    2. Move data off /big7_disk if /big7_disk is the offender
    3. Manual override (use sparingly):
         rm $BLOCK_FLAG
       or one-shot:
         PIPELINE_BLOCK_OVERRIDE=1 <your-command>
EOF
    exit 2
fi

exit 0
