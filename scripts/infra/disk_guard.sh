#!/usr/bin/env bash
# disk_guard.sh — Periodic disk usage monitor for /tmp safety
#
# Source: prior art D-1 (Prometheus textfile collector) + D-2 (bash df cron)
# Plan ref: ~/.claude/plans/agent-harness-langgraph-resilient-waterfall.md §12.2 T1-1
# Lessons: feedback_tmp_disk_full_pipeline_pitfall.md (2026-05-08 800GB disaster)
#
# Behavior:
#   1. Sample df for / and /tmp each run
#   2. Append timestamped JSONL to $LOG_FILE (default /var/log/disk_guard.log; falls back to ~/.disk_guard.log if not writable)
#   3. WARN  threshold (default 90%) → touch $WARN_FLAG + logger
#   4. BLOCK threshold (default 95%) → touch $BLOCK_FLAG (read by pipeline_block_check.sh hook)
#   5. When usage drops below thresholds, remove flags
#
# Install (cron, no sudo):
#   crontab -e
#   */10 * * * * /big7_disk/liaoyoyo2001/InterSubMod/scripts/infra/disk_guard.sh
#
# Install (systemd timer, requires sudo — see disk_guard.timer + disk_guard.service templates):
#   sudo cp scripts/infra/disk_guard.{service,timer} /etc/systemd/system/
#   sudo systemctl enable --now disk_guard.timer

set -euo pipefail

# --- Config (override via env) ---
WARN_PCT="${DISK_GUARD_WARN_PCT:-90}"
BLOCK_PCT="${DISK_GUARD_BLOCK_PCT:-95}"
LOG_FILE="${DISK_GUARD_LOG:-/var/log/disk_guard.log}"
# Flags written to user-writable location (default: $HOME). Override via env if shared between users.
WARN_FLAG="${DISK_GUARD_WARN_FLAG:-$HOME/.intersubmod_disk_warn}"
BLOCK_FLAG="${DISK_GUARD_BLOCK_FLAG:-$HOME/.intersubmod_pipeline_blocked}"
MOUNTS_TO_WATCH="${DISK_GUARD_MOUNTS:-/ /tmp /big7_disk}"

# --- Fallback log path if /var/log not writable ---
if ! touch "$LOG_FILE" 2>/dev/null; then
    LOG_FILE="$HOME/.disk_guard.log"
fi

# --- Sample disk usage; record JSONL ---
TS="$(date -Iseconds)"
declare -A USAGE
WORST_PCT=0
WORST_MOUNT=""

for MOUNT in $MOUNTS_TO_WATCH; do
    if [ -d "$MOUNT" ]; then
        # df --output=pcent gives percentage; strip '%' and whitespace
        PCT="$(df --output=pcent "$MOUNT" 2>/dev/null | tail -1 | tr -d ' %')"
        if [[ "$PCT" =~ ^[0-9]+$ ]]; then
            USAGE[$MOUNT]="$PCT"
            if [ "$PCT" -gt "$WORST_PCT" ]; then
                WORST_PCT="$PCT"
                WORST_MOUNT="$MOUNT"
            fi
        fi
    fi
done

# Build JSON line
JSON_LINE="{\"ts\":\"$TS\""
for MOUNT in "${!USAGE[@]}"; do
    KEY="$(echo "$MOUNT" | tr '/' '_' | sed 's/^_//;s/^$/root/')"
    JSON_LINE+=",\"$KEY\":${USAGE[$MOUNT]}"
done
JSON_LINE+=",\"worst_pct\":$WORST_PCT,\"worst_mount\":\"$WORST_MOUNT\"}"
echo "$JSON_LINE" >> "$LOG_FILE"

# --- Threshold actions ---
if [ "$WORST_PCT" -ge "$BLOCK_PCT" ]; then
    touch "$BLOCK_FLAG"
    echo "$TS BLOCK reason=$WORST_MOUNT@${WORST_PCT}%" > "${BLOCK_FLAG}.reason"
    command -v logger >/dev/null && logger -t disk_guard "CRITICAL: $WORST_MOUNT at ${WORST_PCT}% — pipeline BLOCKED"
    exit 0
fi

# Below block threshold: clear block flag if present (was an earlier alert that recovered)
[ -f "$BLOCK_FLAG" ] && rm -f "$BLOCK_FLAG" "${BLOCK_FLAG}.reason"

if [ "$WORST_PCT" -ge "$WARN_PCT" ]; then
    touch "$WARN_FLAG"
    echo "$TS WARN reason=$WORST_MOUNT@${WORST_PCT}%" > "${WARN_FLAG}.reason"
    command -v logger >/dev/null && logger -t disk_guard "WARN: $WORST_MOUNT at ${WORST_PCT}%"
    exit 0
fi

# All clear: remove warn flag
[ -f "$WARN_FLAG" ] && rm -f "$WARN_FLAG" "${WARN_FLAG}.reason"
