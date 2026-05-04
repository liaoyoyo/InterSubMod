#!/bin/bash
# scan_weekly.sh — weekly-report v2 W1 raw data 自動掃描
# 整合 8 個 grep/find 命令為單一入口
#
# Usage:
#   scripts/analysis/scan_weekly.sh                       # 預設 7 天 + 全主題
#   scripts/analysis/scan_weekly.sh --days 14             # 14 天範圍
#   scripts/analysis/scan_weekly.sh --keyword "self-phasing|methylation"  # 主題過濾
#   scripts/analysis/scan_weekly.sh --raw                 # 顯示原始 grep output（C0 expand mode，修正點 5）

set -e
cd "$(git rev-parse --show-toplevel)"

DAYS=7
KEYWORD=""
RAW_OUTPUT=0

while [[ $# -gt 0 ]]; do
  case "$1" in
    --days) DAYS="$2"; shift 2 ;;
    --keyword) KEYWORD="$2"; shift 2 ;;
    --raw) RAW_OUTPUT=1; shift ;;
    *) echo "Unknown option: $1"; exit 1 ;;
  esac
done

echo "=========================================="
echo "weekly-report W1 raw data 掃描"
echo "Range: 過去 $DAYS 天"
echo "Keyword: ${KEYWORD:-全主題}"
echo "Raw output: $([ $RAW_OUTPUT -eq 1 ] && echo 'YES' || echo 'NO (摘要)')"
echo "=========================================="

echo -e "\n=== 1. Git log (last $DAYS days) ==="
# v2.2 fix: 加 --extended-regexp (-E) 支援 OR pattern (e.g. "self.phasing|V5")
# git log --grep 預設用 BRE，| 不被當 OR；需 --extended-regexp 才能正確處理
if [ -n "$KEYWORD" ]; then
  git log --since="$DAYS days ago" --grep="$KEYWORD" -i --extended-regexp --pretty=format:"%h %ad %s" --date=short
else
  git log --since="$DAYS days ago" --pretty=format:"%h %ad %s" --date=short
fi

if [ -n "$KEYWORD" ]; then
  echo -e "\n\n=== 2. Git log (keyword '$KEYWORD' all-time, last 30) ==="
  git log --all --grep="$KEYWORD" -i --extended-regexp --pretty=format:"%h %ad %s" --date=short | head -30
fi

echo -e "\n\n=== 3. experiments/in_progress 最近修改 (top 15) ==="
ls -t docs/experiments/in_progress/2026/*/2026*.md 2>/dev/null | head -15

if [ -n "$KEYWORD" ]; then
  echo -e "\n\n=== 4. keyword '$KEYWORD' in experiments (last $DAYS days) ==="
  for f in $(find docs/experiments/in_progress/2026/ -name "2026*.md" -mtime -$DAYS 2>/dev/null); do
    grep -l -i "$KEYWORD" "$f" 2>/dev/null
  done | head -15
fi

echo -e "\n\n=== 5. evidence_ledger 最近條目 ==="
if [ -n "$KEYWORD" ]; then
  grep -i "$KEYWORD" research/autoresearch/evidence_ledger.jsonl 2>/dev/null | tail -10
else
  tail -10 research/autoresearch/evidence_ledger.jsonl 2>/dev/null
fi

echo -e "\n\n=== 6. CURRENT_FOCUS 主軸切換偵測 (last $DAYS days commits) ==="
git log --follow --since="$DAYS days ago" --pretty=format:"%h %ad %s" --date=short docs/CURRENT_FOCUS.md 2>/dev/null | head -5

echo -e "\n\n=== 7. 上週 weekly report (top 3) ==="
ls -t docs/reports/validated/2026/*/2026*.md 2>/dev/null | head -3

echo -e "\n\n=== 8. Memory 中相關條目 ==="
MEMORY_DIR="/bip7_disk/liaoyoyo2001/.claude/projects/-big7-disk-liaoyoyo2001-InterSubMod/memory"
if [ -d "$MEMORY_DIR" ]; then
  if [ -n "$KEYWORD" ]; then
    grep -l -i "$KEYWORD" "$MEMORY_DIR"/*.md 2>/dev/null | head -10
  else
    ls -t "$MEMORY_DIR"/*.md 2>/dev/null | head -5
  fi
fi

if [ $RAW_OUTPUT -eq 1 ]; then
  echo -e "\n\n=========================================="
  echo "=== RAW OUTPUT MODE (C0 expand) ==="
  echo "=========================================="
  if [ -n "$KEYWORD" ]; then
    echo -e "\n=== 9. 全文 grep '$KEYWORD' in experiments (raw lines) ==="
    grep -rn -i "$KEYWORD" docs/experiments/in_progress/2026/ 2>/dev/null | head -30
    echo -e "\n=== 10. 全文 grep '$KEYWORD' in reports/validated 2026 ==="
    grep -rn -i "$KEYWORD" docs/reports/validated/2026/ 2>/dev/null | head -30
  fi
fi

echo -e "\n\n=========================================="
echo "掃描完成。AI 應將上述 1-8 (+9-10 if --raw) 整理為 W1 結構化清單，提交 C0 確認。"
echo "=========================================="
