#!/bin/bash
# watch_session.sh — Real-time observability dashboard for long cycles.
#
# Purpose: 對齊 Walking Labs L11 「Observability belongs inside the harness」+
#          Anthropic Harness Design 結構化交付.
# Source:  P4 Top 3 high-ROI fix (T1)
# Usage:
#   bash InterSubMod/scripts/watch_session.sh               # default 5s refresh
#   bash InterSubMod/scripts/watch_session.sh 10            # 10s refresh
#   bash InterSubMod/scripts/watch_session.sh --no-clear    # append (no clear screen)
# Stop: Ctrl+C
#
# Streams to terminal:
#   - CURRENT_FOCUS.md head (top 30 lines)
#   - Recent git commits (last 5)
#   - Cache telemetry summary (latest session)
#   - Hook audit log (skill_audit + hook_failures + subagent_completion latest)
#   - Active /tmp pipeline files (size + mtime)
#   - Recent figures in research/ (last 5)

set -o pipefail

REFRESH="${1:-5}"
NO_CLEAR=""
case "${2:-}" in
    --no-clear) NO_CLEAR="1" ;;
esac
# Handle position swap
[ "$REFRESH" = "--no-clear" ] && { NO_CLEAR="1"; REFRESH="5"; }

ROOT="/big7_disk/liaoyoyo2001/InterSubMod"
YYYYMM=$(date +%Y%m)

trap 'echo ""; echo "[watch_session] stopped."; exit 0' INT TERM

while true; do
    [ -z "$NO_CLEAR" ] && clear

    cat <<HEADER
╔════════════════════════════════════════════════════════════════════════╗
║ InterSubMod Session Watch · $(date +'%Y-%m-%d %H:%M:%S')  · refresh=${REFRESH}s ║
╚════════════════════════════════════════════════════════════════════════╝
HEADER

    # 1. CURRENT_FOCUS.md
    echo ""
    echo "─── §1 CURRENT_FOCUS.md (head 12 lines) ───"
    [ -f "$ROOT/docs/CURRENT_FOCUS.md" ] && head -12 "$ROOT/docs/CURRENT_FOCUS.md" || echo "  (missing)"

    # 2. Recent commits
    echo ""
    echo "─── §2 Recent commits (last 5) ───"
    cd "$ROOT" && git log -5 --oneline 2>/dev/null || echo "  (git unavailable)"

    # 3. Cache telemetry (latest session)
    echo ""
    echo "─── §3 Cache effectiveness (latest session) ───"
    bash "$ROOT/scripts/hooks/cache_telemetry.sh" 2>/dev/null | head -16 | tail -12

    # 4. Hook logs (skill_audit + hook_failures + subagent_completion)
    echo ""
    echo "─── §4 Hook logs (last 3 entries each) ───"
    for log in "skill_audit_${YYYYMM}.log" "hook_failures.log" "subagent_completion_${YYYYMM}.log"; do
        path="$ROOT/docs/postmortems/$log"
        if [ -f "$path" ]; then
            echo "  [$log]"
            tail -3 "$path" | sed 's/^/    /'
        fi
    done

    # 5. Active /tmp pipeline outputs (recent)
    echo ""
    echo "─── §5 Active /tmp pipeline files (recent 5, by mtime) ───"
    ls -t /tmp/ism_*/ /tmp/skill_*/ 2>/dev/null | head -10 | sed 's/^/  /' || echo "  (none)"

    # 6. Recent figures
    echo ""
    echo "─── §6 Recent research figures (last 5 *.png) ───"
    find "$ROOT/research" "$ROOT/docs/experiments/in_progress" -name "*.png" -newer "$ROOT/.git/HEAD" 2>/dev/null | head -5 | sed 's/^/  /' || true
    echo "  (newer than HEAD)"

    # 7. Active disk usage warning
    echo ""
    echo "─── §7 Disk status ───"
    df -h /big7_disk /big8_disk /tmp 2>/dev/null | awk 'NR==1 || /\/big|\/tmp/' | head -4

    echo ""
    echo "─── (Ctrl+C to stop) ───"

    [ -z "$NO_CLEAR" ] && sleep "$REFRESH" || { sleep "$REFRESH"; echo ""; }
done
