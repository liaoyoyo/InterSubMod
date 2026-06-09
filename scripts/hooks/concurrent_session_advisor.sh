#!/bin/bash
# SessionStart hook — concurrent-session collision advisor (git governance §G).
#
# Purpose:  Detect when >1 Claude session is active in the MAIN repo working dir
#           at once (they share ONE git HEAD → one session's `git checkout -b`
#           silently re-homes another session's commits onto the wrong branch).
#           Real incident 2026-06-09: 2 stray commits landed on a harness branch.
#
# Signal:   count of recently-modified *.jsonl in THIS repo's transcript project
#           dir. Worktree sessions get a DIFFERENT transcript dir (cwd-derived),
#           so a fresh .jsonl here == a session sharing the main dir+HEAD.
#           >=2 fresh (mine + >=1 other) → concurrent → advise git worktree.
#
# Design:   advisory only (always returns 0, no blocking code), fail-OPEN, fast
#           (one dir stat). Not /loop|/goal|cron. Pairs with C7 health advisor.

set -uo pipefail

# Read the hook's SessionStart JSON (stdin) BEFORE the heredoc consumes stdin as
# the python script; pass it via env so python can exclude THIS session's own
# transcript precisely (robust even if our own .jsonl is not yet flushed).
HOOK_INPUT="$(cat 2>/dev/null || true)"
export HOOK_INPUT

python3 << 'PYEOF'
import json, os, time, glob

# Transcript project dir for the MAIN repo (cwd /big7_disk/.../InterSubMod).
# Worktree sessions resolve to a different (cwd-derived) project dir, so counting
# fresh .jsonl HERE isolates "sessions sharing the main working dir + HEAD".
TDIR = "/bip7_disk/liaoyoyo2001/.claude/projects/-big7-disk-liaoyoyo2001-InterSubMod"
WINDOW_S = 180   # "active" = touched within this window

def emit(ctx):
    print(json.dumps({
        "hookSpecificOutput": {"hookEventName": "SessionStart", "additionalContext": ctx}
    }, ensure_ascii=False))

# Identify our own transcript so we count only OTHER sessions.
my_tx = ""
try:
    hi = json.loads(os.environ.get("HOOK_INPUT", "") or "{}")
    my_tx = os.path.realpath(hi.get("transcript_path", "")) if hi.get("transcript_path") else ""
except Exception:
    my_tx = ""

try:
    if not os.path.isdir(TDIR):
        raise SystemExit(0)
    now = time.time()
    others = 0
    for f in glob.glob(os.path.join(TDIR, "*.jsonl")):
        try:
            if os.path.realpath(f) == my_tx:
                continue   # skip self
            if now - os.path.getmtime(f) < WINDOW_S:
                others += 1
        except OSError:
            continue
except Exception:
    raise SystemExit(0)   # fail-OPEN: never block session start

# Fallback when transcript_path is unavailable (my_tx==""): one of the fresh
# files is ours, so require >=2 fresh to claim a concurrent OTHER session.
if not my_tx and others >= 1:
    others -= 1

if others >= 1:
    emit(
        f"[concurrent-session advisor] 偵測到主 repo working dir 另有 ~{others} "
        f"個近期活躍 Claude session（共用同一 git HEAD）。\n"
        f"⚠ 並行 session 共用主 dir 會互撞：一個 session 的 `git checkout -b` / commit "
        f"會把別 session 的 commit 落到錯 branch（2026-06-09 已發生）。\n"
        f"→ 並行工作請各開 git worktree（`git worktree add ../wt-<topic> <branch>`）後在該 dir 起 session；"
        f"見 git governance §G。本 session 若要動 branch/commit 前先確認無並行 session。"
    )
raise SystemExit(0)
PYEOF

# Advisory only — never block session start.
exit 0
