#!/usr/bin/env bash
# install_disk_guard.sh — Install disk_guard via cron (preferred, no sudo) or systemd timer
#
# Usage:
#   ./install_disk_guard.sh cron      # add to current user's crontab (no sudo)
#   ./install_disk_guard.sh systemd   # copy unit files to /etc/systemd/system (requires sudo)
#   ./install_disk_guard.sh status    # show current install state

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
GUARD="$REPO_ROOT/scripts/infra/disk_guard.sh"
CRON_LINE="*/10 * * * * $GUARD"

case "${1:-status}" in
    cron)
        # Check if entry already present
        if crontab -l 2>/dev/null | grep -F "$GUARD" >/dev/null; then
            echo "[skip] disk_guard already in crontab"
        else
            (crontab -l 2>/dev/null; echo "$CRON_LINE") | crontab -
            echo "[ok] added to crontab: $CRON_LINE"
        fi
        ;;
    systemd)
        if [ "$EUID" -ne 0 ]; then
            echo "ERROR: systemd install requires sudo. Re-run with: sudo ./install_disk_guard.sh systemd"
            exit 1
        fi
        cp "$REPO_ROOT/scripts/infra/disk_guard.service" /etc/systemd/system/
        cp "$REPO_ROOT/scripts/infra/disk_guard.timer"   /etc/systemd/system/
        systemctl daemon-reload
        systemctl enable --now disk_guard.timer
        echo "[ok] systemd timer enabled. Status:"
        systemctl status disk_guard.timer --no-pager
        ;;
    status)
        echo "== crontab =="
        if crontab -l 2>/dev/null | grep -F "$GUARD" >/dev/null; then
            echo "  [installed] $(crontab -l | grep -F "$GUARD")"
        else
            echo "  [not installed]"
        fi
        echo "== systemd =="
        if systemctl is-enabled disk_guard.timer 2>/dev/null | grep -q enabled; then
            echo "  [enabled]"
            systemctl status disk_guard.timer --no-pager 2>&1 | head -5
        else
            echo "  [not installed]"
        fi
        echo "== last log =="
        for LOG in /var/log/disk_guard.log "$HOME/.disk_guard.log"; do
            if [ -f "$LOG" ]; then
                echo "  $LOG (last 3 lines):"
                tail -3 "$LOG" | sed 's/^/    /'
                break
            fi
        done
        ;;
    *)
        echo "Usage: $0 {cron|systemd|status}"
        exit 1
        ;;
esac
