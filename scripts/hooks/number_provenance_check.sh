#!/bin/bash
# number_provenance_check.sh — PreToolUse Edit|Write anti-fabrication gate (Layer B).
#
# Tiered (per 2026-06-01 fabricated-metric postmortem, action items A4/A5/A7):
#   - validated/ or pi_reports/ paths -> STRICT: unsourced metric => exit 2 (block write)
#   - other paths (in_progress, etc.) -> ADVISORY: remind, exit 0
#
# FAIL-OPEN: if python3 is missing or errors, exit 0 (a broken hook must NEVER block
# all doc writes). This is NOT the `|| exit 0` neutering bug — genuine detections in
# strict paths still propagate exit 2 because we do NOT mask the python exit code here.
#
# Wiring (settings.local.json) MUST be a plain call (no `|| true` / `|| exit 0`):
#   bash scripts/hooks/number_provenance_check.sh
set -o pipefail

ROOT="/big7_disk/liaoyoyo2001/InterSubMod"
PY="$ROOT/scripts/number_provenance.py"

# fail-open if interpreter or script absent
command -v python3 >/dev/null 2>&1 || exit 0
[ -f "$PY" ] || exit 0

# pass the hook JSON (stdin) straight through; python decides advisory vs exit 2.
python3 "$PY" gate
rc=$?
# only 0 (pass/advisory) and 2 (strict block) are meaningful; map anything else to 0.
if [ "$rc" -eq 2 ]; then
    exit 2
fi
exit 0
