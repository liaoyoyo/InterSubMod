#!/usr/bin/env bash
set -uo pipefail

# disk_audit.sh - audit InterSubMod project disk usage
# Read-only. Reports sizes, never modifies.

REPO_ROOT="${REPO_ROOT:-/big7_disk/liaoyoyo2001/InterSubMod}"

echo "=== InterSubMod disk audit ==="
echo "Repo root: $REPO_ROOT"
echo ""

echo "--- Top-level directories ---"
du -sh "$REPO_ROOT"/* 2>/dev/null | sort -h | tail -20

echo ""
echo "--- output/ subdirectories ---"
if [ -d "$REPO_ROOT/output" ]; then
    du -sh "$REPO_ROOT/output"/* 2>/dev/null | sort -h | tail -10
else
    echo "(no output/ dir)"
fi

echo ""
echo "--- build/ ---"
if [ -d "$REPO_ROOT/build" ]; then
    du -sh "$REPO_ROOT/build" 2>/dev/null
else
    echo "(no build/ dir)"
fi

echo ""
echo "--- research/ subdirectories ---"
if [ -d "$REPO_ROOT/research" ]; then
    du -sh "$REPO_ROOT/research"/* 2>/dev/null | sort -h | tail -10
else
    echo "(no research/ dir)"
fi

echo ""
echo "--- Disk summary ---"
df -h "$REPO_ROOT"
