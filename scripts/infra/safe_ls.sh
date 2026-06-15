#!/usr/bin/env bash
# safe_ls.sh — §12-safe stat of canonical ISM run dir(s) WITHOUT du/find recursion.
#
# WHY: ISM region dirs (intersubmod_{tp,fp}, filtered_snv_{tp,fp}) hold >10,000 subdirs
#      each → `du`/`find`/`ls -R`/`wc -l <glob>/*/...` HANG (2026-06-15 audit D3-1, the #1
#      "卡住" cause). This encodes the safe triple — run mtime / BAM sizes / one-level region
#      count — so the inventory commands can't be mistyped into a hang. Generalizes the
#      `ls -d ... | wc -l` idiom already used in scripts/verify_output.sh (D6-4).
#
# Usage: bash scripts/infra/safe_ls.sh <run_dir> [<run_dir> ...]
#   e.g. bash scripts/infra/safe_ls.sh /big7_disk/.../canonical/HCC1395/paired_full/*_complete_matrix/
set -uo pipefail
[ "$#" -eq 0 ] && { echo "usage: safe_ls.sh <run_dir> [<run_dir> ...]"; exit 1; }

REGION_PARENTS=(intersubmod_tp intersubmod_fp filtered_snv_tp filtered_snv_fp)

for d in "$@"; do
    d="${d%/}"
    if [ ! -d "$d" ]; then echo "❌ not a dir: $d"; continue; fi
    echo "=== $d ==="
    # 1) run mtime — single dir stat, NO recursion
    echo "  mtime: $(ls -dl --time-style=+%Y-%m-%d_%H:%M "$d" 2>/dev/null | awk '{print $6}')"
    # 2) BAM sizes — ls -la ONLY (never du); haplotagged BAM lives in longphase_s/
    bams=$(ls -la "$d"/*.bam "$d"/longphase_s/*.bam 2>/dev/null | awk '{printf "    %12s  %s\n",$5,$NF}')
    if [ -n "$bams" ]; then echo "  BAMs (bytes):"; echo "$bams"; else echo "  BAMs: (none at top-level / longphase_s)"; fi
    # 3) region subdir count — ONE level only (ls -d .../*/), never find / ls -R / recursive wc
    for p in "${REGION_PARENTS[@]}"; do
        if [ -d "$d/$p" ]; then
            n=$(ls -d "$d/$p"/*/ 2>/dev/null | wc -l)
            echo "  $p: $n region dirs (one-level)"
        fi
    done
done
