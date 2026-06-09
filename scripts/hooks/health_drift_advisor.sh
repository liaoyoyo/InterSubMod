#!/bin/bash
# SessionStart hook — C7: change-gated read-only health-drift advisor.
#
# Purpose:  Fix the cross-cutting "invocation-dependence" gap (harness_health /
#           cycle-state only run when someone remembers). Detect whether harness
#           config changed since last check, and nudge the user to run
#           /harness-health — WITHOUT running it eagerly (no latency, no snapshot
#           spam, no hang risk).
#
# Design (per 20260609 Loop-Engineering ADR §5, ADOPT_WITH_GUARDRAILS):
#   - CHANGE-GATED: only speaks when .claude/{skills,agents,rules,hooks}/ or
#     settings/CLAUDE.md/AGENTS.md changed since the marker. Silent otherwise
#     (most pure code-edit/single-doc sessions get nothing — aligns CLAUDE.md §6).
#   - ADVISORY-ONLY (advisory only): always returns 0, never blocks session
#     start (this is NOT a hard-gate; it has no blocking exit code).
#   - ZERO research numbers: harness-maintenance telemetry only; no report/ledger/
#     cycle writes, no findings. Does NOT invoke harness_health.py.
#   - NOT /loop|/goal|cron: a plain SessionStart event hook (precedent:
#     session_start_inject_focus.sh). Marker debounces to ~once per change-batch.
#
# Side effect: touches state/health_snapshots/.drift_advisor_marker after advising.

set -uo pipefail

python3 << 'PYEOF'
import json, os, time

REPO = "/big7_disk/liaoyoyo2001/InterSubMod"
MARKER = os.path.join(REPO, "state", "health_snapshots", ".drift_advisor_marker")
WATCH_DIRS = [
    ".claude/skills", ".claude/agents", ".claude/rules", "scripts/hooks",
    ".claude/workflows",
]
WATCH_FILES = [".claude/settings.local.json", ".claude/CLAUDE.md", "AGENTS.md"]
WATCH_EXT = (".md", ".sh", ".json", ".js")
STALE_DAYS = 7          # cycle no-advance threshold
DEBOUNCE_H = 20         # don't re-evaluate stale-only nudges more than ~once/20h


def emit(ctx):
    print(json.dumps({
        "hookSpecificOutput": {
            "hookEventName": "SessionStart",
            "additionalContext": ctx,
        }
    }, ensure_ascii=False))


def safe_exit():
    raise SystemExit(0)


# --- marker mtime (0 if absent → treat as first run = changed) ---
marker_mtime = os.path.getmtime(MARKER) if os.path.exists(MARKER) else 0.0
first_run = marker_mtime == 0.0

# --- count harness files changed since marker (scoped, light dirs only) ---
changed = []
try:
    for d in WATCH_DIRS:
        full = os.path.join(REPO, d)
        if not os.path.isdir(full):
            continue
        for root, _dirs, files in os.walk(full):
            for fn in files:
                if not fn.endswith(WATCH_EXT):
                    continue
                p = os.path.join(root, fn)
                try:
                    if os.path.getmtime(p) > marker_mtime:
                        changed.append(os.path.relpath(p, REPO))
                        if len(changed) >= 60:   # cap scan cost
                            raise StopIteration
                except OSError:
                    continue
    for f in WATCH_FILES:
        p = os.path.join(REPO, f)
        if os.path.exists(p):
            try:
                if os.path.getmtime(p) > marker_mtime:
                    changed.append(f)
            except OSError:
                pass
except StopIteration:
    pass
except Exception:
    # fail-OPEN: any scan error → stay silent, never block startup
    safe_exit()

# --- stale active cycles (no advance > STALE_DAYS) ---
stale_cycles = []
try:
    aj = os.path.join(REPO, "state", "active.json")
    if os.path.exists(aj):
        active = json.load(open(aj, encoding="utf-8"))
        now = time.time()
        for c in active.get("cycles", []):
            la = c.get("last_advanced_at", "")
            # robust parse: take the YYYY-MM-DD prefix (ignore time/tz)
            try:
                tt = time.mktime(time.strptime(la[:10], "%Y-%m-%d"))
                age_d = (now - tt) / 86400.0
                if age_d > STALE_DAYS:
                    stale_cycles.append((c.get("cycle_id", "?"), int(age_d),
                                         c.get("phase", "?")))
            except Exception:
                continue
except Exception:
    stale_cycles = []

# --- decide whether to speak ---
harness_changed = first_run or bool(changed)
marker_fresh = (not first_run) and (time.time() - marker_mtime) < DEBOUNCE_H * 3600

# If harness unchanged AND we advised recently → stay silent (debounce stale nudges)
if not harness_changed and marker_fresh:
    safe_exit()
# If nothing to say at all → silent
if not harness_changed and not stale_cycles:
    safe_exit()

# --- build advisory (1-3 lines) ---
lines = ["[harness-drift advisor]"]
if harness_changed:
    n = len(changed)
    if first_run:
        lines.append("首次啟用 drift advisor。建議跑 `/harness-health` 建立基準 8 燈儀表板。")
    else:
        sample = ", ".join(os.path.basename(c) for c in changed[:4])
        more = f" 等 {n} 檔" if n > 4 else ""
        lines.append(f"harness 設定自上次稽核後有變動（{sample}{more}）→ 跑 `/harness-health` 確認無 count-drift / Hard-Gate neutering / state↔doc 漂移。")
if stale_cycles:
    cid, age, ph = stale_cycles[0]
    extra = f"（另 {len(stale_cycles)-1} 個）" if len(stale_cycles) > 1 else ""
    lines.append(f"active cycle `{cid}` 已 {age} 天未推進（{ph}）{extra} → `/cycle-state` 檢視是否該推進或 archive。")
lines.append("（唯讀提醒；如不需要可在 settings 移除 health_drift_advisor）")

emit("\n".join(lines))

# --- touch marker so next session is quiet until new change ---
try:
    os.makedirs(os.path.dirname(MARKER), exist_ok=True)
    with open(MARKER, "w", encoding="utf-8") as f:
        f.write(time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()) + "\n")
except OSError:
    pass
PYEOF

# Always succeed — advisory hook must never block session start.
exit 0
