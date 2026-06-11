#!/usr/bin/env bash
# stale_redirect_banner.sh — 在已知 stale 索引/狀態快照檔插 redirect banner，指向現役 CURRENT_FOCUS。
# 改進 ⑥（D3 / pitfalls P-17）：pivot 後不重寫 stale 內容，只標「過時 + 指路」，把「誤導」降級為「明示過時」。
# frontmatter-aware：YAML(---) / HTML-comment(<!--) 起頭的檔，banner 插在 frontmatter 之後（不破壞 KB MCP 解析 last_verified）；純文字插最上面。
# Idempotent（marker 檢查，可重跑）。
# 用法：bash scripts/stale_redirect_banner.sh            # 跑預設已知 stale 清單
#       bash scripts/stale_redirect_banner.sh <rel1> ..  # 指定相對 repo root 的檔
set -euo pipefail
REPO="$(git rev-parse --show-toplevel 2>/dev/null || echo /big7_disk/liaoyoyo2001/InterSubMod)"
MARKER="<!-- STALE-REDIRECT-BANNER (scripts/stale_redirect_banner.sh) -->"
POINTER="InterSubMod/docs/CURRENT_FOCUS.md"
NOTE="主軸已於 2026-06-11 pivot 至 Subclonal reconstruction（取代 G6）"

DEFAULTS=(
  "knowledge/10_research_status/01_current_focus_snapshot.md"
  "knowledge/10_research_status/02_active_hypotheses.md"
  "knowledge/10_research_status/04_blockers_and_risks.md"
  "knowledge/10_research_status/05_next_milestones.md"
  "knowledge/09_conclusions/05_hypothesis_queue_snapshot.md"
  "docs/reports/INDEX.md"
)
files=("$@"); [ "${#files[@]}" -eq 0 ] && files=("${DEFAULTS[@]}")

for rel in "${files[@]}"; do
  f="$REPO/$rel"
  [ -f "$f" ] || { echo "SKIP (missing): $rel"; continue; }
  if grep -qF "$MARKER" "$f"; then echo "OK (already bannered): $rel"; continue; fi
  mtime="$(date -r "$f" +%Y-%m-%d 2>/dev/null || echo unknown)"
  text="> ⚠ **此檔為 ${mtime} 前後快照，可能已過時** — 現役主軸/狀態以 \`${POINTER}\` 為準（${NOTE}）。本檔僅供歷史對照，勿據此判斷現況。"
  tmp="$(mktemp)"
  awk -v marker="$MARKER" -v text="$text" '
    function emit(){ print ""; print marker; print text; print "" }
    BEGIN{ inserted=0; mode=""; n=0 }
    { n++
      if(n==1){
        if($0=="---"){ mode="yaml"; print; next }
        else if($0 ~ /^<!--/){ mode="html"; print; if($0 ~ /-->/){ emit(); inserted=1 } next }
        else { emit(); inserted=1; print; next }
      }
      print
      if(!inserted){
        if(mode=="yaml" && $0=="---"){ emit(); inserted=1 }
        else if(mode=="html" && $0 ~ /-->/){ emit(); inserted=1 }
      }
    }
    END{ if(!inserted){ emit() } }   # 安全網：未偵測到 frontmatter 結尾也補上
  ' "$f" > "$tmp"
  mv "$tmp" "$f"
  echo "BANNERED: $rel (snapshot ~${mtime})"
done
