#!/usr/bin/env bash
# archive_candidate.sh — auto-build the archive manifest ROW for candidate run dir(s).
#
# WHY: the archive workflow (grep-check → safe mtime/BAM stat → manifest row → printed mv)
#      was done by hand ≥2× (batch1+batch2) with ≥3 more batches pending (2026-06-15 audit
#      D6-1). This encodes the repeated recipe so each batch isn't re-derived manually.
#
# SAFETY:
#   - NEVER moves/deletes anything — only PRINTS the mv command for you to run after review.
#   - §12-safe: ls -la for BAM (no du), grep scoped to scripts/+docs/ (no repo-root recursion).
#   - Deletion stays a Hard Gate (keep >=1 week + user confirm), per the manifests.
#
# Usage: bash scripts/infra/archive_candidate.sh <run_dir> [<run_dir> ...]
set -uo pipefail
REPO="/big7_disk/liaoyoyo2001/InterSubMod"
OUT_ROOT="/big7_disk/liaoyoyo2001/big7_disk_output"
ARCHIVE_ROOT="${OUT_ROOT}/_ARCHIVE_pending_cleanup_202606"
BIG_BAM_BYTES=104857600   # 100 MB — "big BAM" worth keeping out of archive
[ "$#" -eq 0 ] && { echo "usage: archive_candidate.sh <run_dir> [<run_dir> ...]"; exit 1; }

printf '%-58s | %-10s | %-7s | %-4s | %s\n' "run" "mtime" "big_BAM" "refs" "verdict"
printf '%s\n' "----------------------------------------------------------------------------------------"
declare -a HOLD_NOTES=()
for d in "$@"; do
    d="${d%/}"
    [ -d "$d" ] || { echo "❌ not a dir: $d"; continue; }
    base=$(basename "$d")
    mtime=$(ls -dl --time-style=+%Y-%m-%d "$d" 2>/dev/null | awk '{print $6}')
    # big-BAM count via ls -la (NEVER du): top-level + longphase_s/
    nbam=$(ls -la "$d"/*.bam "$d"/longphase_s/*.bam 2>/dev/null | awk -v t="$BIG_BAM_BYTES" '$5>t' | wc -l)
    # §12-safe scoped grep-check: is the run basename a FUNCTIONAL dep (referenced by script/doc)?
    refs=$(grep -rIl --include='*.py' --include='*.sh' --include='*.md' -- "$base" \
             "$REPO/scripts" "$REPO/docs" 2>/dev/null | wc -l)
    if [ "$refs" -gt 0 ]; then
        verdict="⏸ HOLD (referenced)"
        HOLD_NOTES+=("$(grep -rIln --include='*.py' --include='*.sh' --include='*.md' -- "$base" "$REPO/scripts" "$REPO/docs" 2>/dev/null | sed "s|^|      $base ref: |")")
    else
        verdict="✅ archive"
    fi
    printf '%-58s | %-10s | %-7s | %-4s | %s\n' "$base" "${mtime:-?}" "$nbam" "$refs" "$verdict"
done

if [ "${#HOLD_NOTES[@]}" -gt 0 ]; then
    echo; echo "# HOLD references (review before archiving — may be a live functional dep):"
    printf '%s\n' "${HOLD_NOTES[@]}"
fi

echo
echo "# ---- After review, mv the ✅ rows (same-fs = instant, recoverable, NOT delete) ----"
for d in "$@"; do
    d="${d%/}"; [ -d "$d" ] || continue
    rel="${d#"${OUT_ROOT}"/}"
    echo "#   mkdir -p \"${ARCHIVE_ROOT}/$(dirname "$rel")\" && mv \"$d\" \"${ARCHIVE_ROOT}/${rel}\""
done
echo "# Deletion = Hard Gate: keep >=1 week in _ARCHIVE_pending_cleanup_202606, user confirm before rm."
