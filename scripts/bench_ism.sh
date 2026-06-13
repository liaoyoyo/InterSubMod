#!/bin/bash
# ISM benchmark wrapper — 每次跑 ISM 記錄 wall-time / max-RSS / disk / commit
# 用途：改動前後比較效能（time/memory/space），對齊參數實驗框架的 cost 欄。
# 用法: bench_ism.sh <label> <out_dir> <benchmark.jsonl> -- <ism_binary 與所有參數，不含 -o>
#   e.g. bench_ism.sh "len1000_chr1" /tmp/out /path/bench.jsonl -- \
#          build/bin/inter_sub_mod -t ... -n ... -r ... -v ... -w 5000 -j 16
set -uo pipefail

LABEL="$1"; OUTDIR="$2"; BENCH="$3"; shift 3
[ "${1:-}" = "--" ] && shift
mkdir -p "$OUTDIR" "$(dirname "$BENCH")"

TL=$(mktemp)
# /usr/bin/time -v 給 Elapsed(wall) + Maximum resident set size(max RSS)
/usr/bin/time -v "$@" -o "$OUTDIR" > "${OUTDIR%/}.runlog" 2>"$TL"
RC=$?

WALL=$(grep -i "Elapsed (wall" "$TL" | grep -oE "[0-9:.]+$" | tail -1)
RSS_KB=$(grep -i "Maximum resident" "$TL" | grep -oE "[0-9]+$" | tail -1)
DISK_KB=$(du -sk "$OUTDIR" 2>/dev/null | cut -f1)
COMMIT=$(git -C "$(cd "$(dirname "$0")/.." && pwd)" rev-parse --short HEAD 2>/dev/null || echo NA)
TS=$(date +%Y-%m-%dT%H:%M:%S)

# append 一行 JSON（每 run 一筆，可累積比較）
printf '{"ts":"%s","label":"%s","exit":%d,"wall":"%s","max_rss_kb":%s,"disk_kb":%s,"commit":"%s"}\n' \
  "$TS" "$LABEL" "$RC" "${WALL:-NA}" "${RSS_KB:-0}" "${DISK_KB:-0}" "$COMMIT" >> "$BENCH"
echo "[bench] $LABEL: exit=$RC wall=${WALL} maxRSS=$((${RSS_KB:-0}/1024))MB disk=$((${DISK_KB:-0}/1024))MB commit=$COMMIT"
rm -f "$TL"
exit $RC
