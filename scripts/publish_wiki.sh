#!/usr/bin/env bash
# 把 docs/wiki/ 的內容發布到 GitHub Wiki。
#
# 前提（只需做一次）：GitHub 的 wiki 是一個獨立 git repo，而它要等到「網頁上存過
# 第一個頁面」之後才會被建立 —— 這一步沒有 API 可以代勞。若尚未做過：
#   1. 開 https://github.com/liaoyoyo/InterSubMod/wiki
#   2. 點 "Create the first page"，隨便輸入內容按 Save（本腳本之後會覆蓋它）
#
# 用法:
#   bash scripts/publish_wiki.sh            # 預覽會做什麼，不推送
#   bash scripts/publish_wiki.sh --push     # 實際推送
set -euo pipefail

REPO_SLUG="liaoyoyo/InterSubMod"
WIKI_URL="https://github.com/${REPO_SLUG}.wiki.git"
SRC="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)/docs/wiki"
TMP="$(mktemp -d -t ism-wiki-XXXXXX)"
DO_PUSH=0
[ "${1:-}" = "--push" ] && DO_PUSH=1

cleanup() { rm -rf "$TMP"; }
trap cleanup EXIT

echo "來源目錄 : $SRC"
echo "Wiki repo: $WIKI_URL"
echo

# --- 1. 圖片必須已經在 develop 上，否則 wiki 會是一堆破圖 ---
echo "[1/4] 檢查圖片是否已在 develop 上（wiki 用 raw URL 引用）"
missing=0
# 只掃真正會進 wiki 的頁面；README.md 是給開發者看的說明，內含 <name> 佔位符範例
for f in "$SRC"/*.md; do
  [ "$(basename "$f")" = "README.md" ] && continue
  while IFS= read -r name; do
    if ! git cat-file -e "origin/develop:docs/images/${name}.png" 2>/dev/null; then
      echo "  ✘ docs/images/${name}.png 不在 origin/develop"
      missing=$((missing + 1))
    fi
  done < <(grep -oE 'raw\.githubusercontent\.com/[^)]*/docs/images/([^)/]+)\.png' "$f" \
           | sed -E 's#.*/docs/images/([^)]+)\.png#\1#' | sort -u)
done
if [ "$missing" -gt 0 ]; then
  echo
  echo "🔴 有 $missing 個圖片尚未推到 develop —— 先推圖再發布，否則 wiki 上會是破圖。"
  exit 1
fi
echo "  ✔ 所有引用的圖片都在 develop 上"

# --- 2. clone wiki ---
echo "[2/4] clone wiki repo"
if ! GIT_TERMINAL_PROMPT=0 git clone -q "$WIKI_URL" "$TMP/wiki" 2>/dev/null; then
  echo
  echo "🔴 clone 失敗 —— wiki repo 還不存在。"
  echo "   請先到 https://github.com/${REPO_SLUG}/wiki 點 'Create the first page' 存一次，再重跑。"
  exit 1
fi
echo "  ✔ 已取得（現有 $(find "$TMP/wiki" -maxdepth 1 -name '*.md' | wc -l) 個頁面）"

# --- 3. 同步檔案（README.md 是給開發者看的說明，不進 wiki）---
echo "[3/4] 同步頁面"
cp "$SRC"/*.md "$TMP/wiki/"
rm -f "$TMP/wiki/README.md"
cd "$TMP/wiki"
git add -A
if git diff --cached --quiet; then
  echo "  沒有變更，不需發布。"
  exit 0
fi
git --no-pager diff --cached --stat | sed 's/^/  /'

# --- 4. 推送 ---
if [ "$DO_PUSH" -eq 0 ]; then
  echo
  echo "[4/4] 這是預覽模式。確認以上變更無誤後，加 --push 實際發布："
  echo "      bash scripts/publish_wiki.sh --push"
  exit 0
fi

echo "[4/4] 推送"
git -c user.name="InterSubMod Research" -c user.email="noreply@anthropic.com" \
    commit -q -m "docs(wiki): sync from docs/wiki/ (source of truth: docs/explain/)"
git push -q origin HEAD
echo "  ✔ 已發布 → https://github.com/${REPO_SLUG}/wiki"
echo
echo "請開啟上面的網址確認：每頁的圖都有顯示、側邊欄有出現、頁間連結可點。"
