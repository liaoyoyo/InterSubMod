#!/usr/bin/env bash
# Serve /big7_disk/liaoyoyo2001 over http://localhost:8765/
# (the canonical layered_workstation lives under this root — see WS_URL below).
#
# Bound to 127.0.0.1 ONLY (this machine) — never 0.0.0.0, so the confidential
# research files are not exposed to the local network. To view from another
# computer, SSH port-forward instead:  ssh -L 8765:localhost:8765 <host>
#
# Usage:
#   ./serve_layered_workstation.sh [start|stop|restart|status]
#   PORT=9000 ./serve_layered_workstation.sh start          # override port
#   SERVE_DIR=/some/dir ./serve_layered_workstation.sh start # override root
#
# Default action is "start".
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SERVE_DIR="${SERVE_DIR:-/big7_disk/liaoyoyo2001}"
WS_REL="InterSubMod/docs/methodology/_assets/layered_workstation/index.html"
PORT="${PORT:-8765}"
HOST="127.0.0.1"
PID_FILE="/tmp/layered_workstation_http_${PORT}.pid"
LOG_FILE="/tmp/layered_workstation_http_${PORT}.log"
PYTHON="$(command -v python3 || command -v python)"

is_running() {  # echoes PID if a server we started is alive
  if [[ -f "$PID_FILE" ]]; then
    local pid; pid="$(cat "$PID_FILE" 2>/dev/null || true)"
    if [[ -n "$pid" ]] && kill -0 "$pid" 2>/dev/null; then echo "$pid"; return 0; fi
  fi
  return 1
}

port_owner() {  # PID currently listening on $PORT (any owner), if any
  (ss -ltnp 2>/dev/null || netstat -ltnp 2>/dev/null) \
    | grep -oE "127.0.0.1:${PORT}\b.*pid=[0-9]+" | grep -oE 'pid=[0-9]+' | head -1 | cut -d= -f2
}

start() {
  if pid="$(is_running)"; then
    echo "already running (PID $pid) → http://localhost:${PORT}/"; return 0
  fi
  local owner; owner="$(port_owner || true)"
  if [[ -n "${owner:-}" ]]; then
    echo "ERROR: port ${PORT} already in use by PID ${owner} (not started by this script)." >&2
    echo "       Use a different PORT=... or stop that process." >&2
    return 1
  fi
  if [[ ! -d "$SERVE_DIR" ]]; then
    echo "ERROR: serve root $SERVE_DIR not found." >&2
    return 1
  fi
  cd "$SERVE_DIR"
  nohup "$PYTHON" -m http.server "$PORT" --bind "$HOST" > "$LOG_FILE" 2>&1 &
  echo $! > "$PID_FILE"
  sleep 1
  if pid="$(is_running)"; then
    echo "started (PID $pid) → http://localhost:${PORT}/"
    echo "  serving:     $SERVE_DIR"
    echo "  workstation: http://localhost:${PORT}/${WS_REL}"
    echo "  log:         $LOG_FILE"
  else
    echo "ERROR: server failed to start; see $LOG_FILE" >&2
    tail -5 "$LOG_FILE" >&2 || true
    return 1
  fi
}

stop() {
  if pid="$(is_running)"; then
    kill "$pid" 2>/dev/null || true
    sleep 1
    kill -9 "$pid" 2>/dev/null || true
    rm -f "$PID_FILE"
    echo "stopped (PID $pid)"
  else
    rm -f "$PID_FILE"
    echo "not running"
  fi
}

status() {
  if pid="$(is_running)"; then
    echo "RUNNING (PID $pid) → http://localhost:${PORT}/"
    local code; code="$(curl -sS -o /dev/null -w '%{http_code}' "http://localhost:${PORT}/${WS_REL}" 2>/dev/null || echo '---')"
    echo "  workstation HTTP $code · serving $SERVE_DIR"
  else
    echo "NOT running (port ${PORT})"
    local owner; owner="$(port_owner || true)"
    [[ -n "${owner:-}" ]] && echo "  note: port ${PORT} held by another PID ${owner}"
  fi
}

case "${1:-start}" in
  start)   start ;;
  stop)    stop ;;
  restart) stop; start ;;
  status)  status ;;
  *) echo "usage: $0 [start|stop|restart|status]   (PORT=nnnn to override)" >&2; exit 2 ;;
esac
