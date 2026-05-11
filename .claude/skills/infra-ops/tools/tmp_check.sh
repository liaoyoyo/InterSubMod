#!/usr/bin/env bash
set -uo pipefail

# tmp_check.sh - report /tmp state, mount type, TMPDIR safety
# Read-only. Never modifies anything.

echo "=== Root partition ==="
df -h /

echo ""
echo "=== /tmp partition ==="
df -h /tmp

echo ""
echo "=== /tmp mount type ==="
if mount | grep -E '^[^ ]+ on /tmp ' >/dev/null; then
    mount | grep -E '^[^ ]+ on /tmp '
else
    echo "/tmp is NOT a separate mount - uses root volume (DANGEROUS for GB-scale intermediates)"
fi

echo ""
echo "=== Current TMPDIR ==="
echo "TMPDIR=${TMPDIR:-<unset, defaults to /tmp>}"
if [ -z "${TMPDIR:-}" ]; then
    echo "WARNING: TMPDIR unset - pipeline tools may write to /tmp by default"
    echo "Recommended: export TMPDIR=/big7_disk/liaoyoyo2001/tmp && mkdir -p \"\$TMPDIR\""
fi

echo ""
echo "=== Largest /tmp occupants (top 10) ==="
du -sh /tmp/* 2>/dev/null | sort -h | tail -10 || echo "(no /tmp contents accessible)"
