#!/usr/bin/env bash
# gh-write.sh — GitHub CLI wrapper for write operations (PR create, merge, release, etc.)
#
# Problem: Claude Code runs in a non-interactive shell that does not source ~/.bashrc,
# so the read/write token-switching wrapper defined there is never activated.
# Direct `gh` calls fall back to GH_TOKEN (read-only), causing 403 on write operations.
#
# Solution: This script reads GH_TOKEN_WRITE from ~/.config/gh/tokens.env and injects
# it as GH_TOKEN before calling gh, ensuring write operations always succeed.
#
# Usage (anywhere gh write is needed):
#   ./scripts/utils/gh-write.sh pr create --base main --head develop --title "..."
#   ./scripts/utils/gh-write.sh pr merge 4 --squash
#   ./scripts/utils/gh-write.sh release create v1.0.0

TOKEN_FILE="$HOME/.config/gh/tokens.env"

if [[ ! -f "$TOKEN_FILE" ]]; then
    echo "[gh-write] ERROR: Token file not found: $TOKEN_FILE" >&2
    exit 1
fi

WRITE_TOKEN=$(grep '^GH_TOKEN_WRITE=' "$TOKEN_FILE" | cut -d= -f2-)

if [[ -z "$WRITE_TOKEN" ]]; then
    echo "[gh-write] ERROR: GH_TOKEN_WRITE is not set in $TOKEN_FILE" >&2
    echo "[gh-write] Add: GH_TOKEN_WRITE=ghp_yourtoken" >&2
    exit 1
fi

GH_TOKEN="$WRITE_TOKEN" gh "$@"
