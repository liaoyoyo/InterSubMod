#!/usr/bin/env bash
# hook_path_rewrite.sh
# 把拷貝來的 hook scripts 內 hardcoded ISM 絕對路徑改為 ${CLAUDE_PROJECT_DIR}
#
# Usage: bash hook_path_rewrite.sh <hooks_dir>
# Example: bash hook_path_rewrite.sh /new-project/scripts/hooks

set -euo pipefail

HOOKS_DIR="${1:?Usage: $0 <hooks_dir>}"
[ -d "$HOOKS_DIR" ] || { echo "ERROR: $HOOKS_DIR not a directory"; exit 1; }

# ISM 絕對路徑映射表
declare -A REPLACEMENTS=(
  ["/big7_disk/liaoyoyo2001/InterSubMod"]="\${CLAUDE_PROJECT_DIR}"
  ["/bip7_disk/liaoyoyo2001"]="\${HOME}"
  ["/big8_disk/liaoyoyo2001"]="\${HOME}/big8_disk"  # 若新機沒有 big8，手動改
  ["/big7_disk/liaoyoyo2001"]="\${HOME}/big7_disk"  # 若新機沒有 big7，手動改
)

# Backup
BACKUP_DIR="${HOOKS_DIR}/.backup_$(date +%s)"
mkdir -p "$BACKUP_DIR"
cp -r "${HOOKS_DIR}"/*.sh "$BACKUP_DIR/" 2>/dev/null || true

echo "[hook-rewrite] backup → $BACKUP_DIR"

# Rewrite
for f in "${HOOKS_DIR}"/*.sh; do
  [ -f "$f" ] || continue
  for old_path in "${!REPLACEMENTS[@]}"; do
    new_path="${REPLACEMENTS[$old_path]}"
    # 用 | 作分隔避免 / 衝突
    sed -i "s|${old_path}|${new_path}|g" "$f"
  done
  echo "[hook-rewrite] $(basename "$f") rewritten"
done

# 找出殘留的可疑路徑（人工檢視）
echo ""
echo "=== 殘留可疑絕對路徑（請人工檢視）==="
grep -nE "/(big[0-9]+|bip[0-9]+)_disk|/home/liaoyoyo|HCC1395|COLO829|HCC1937|HCC1954|H1437|H2009|longphase|ClairS|InterSubMod|evidence_ledger|cycle_id" "${HOOKS_DIR}"/*.sh 2>/dev/null || echo "  (無殘留)"

echo ""
echo "[hook-rewrite] 完成。下一步："
echo "  1. 人工檢視上述殘留並酌情改寫"
echo "  2. 測試：bash <hook>.sh （某些 hook 需 stdin JSON，可跑 echo '{}' | bash <hook>.sh）"
echo "  3. 確認 hooks 加入 .claude/settings.local.json 的 hooks 區段"
echo "  4. 開新 Claude session 觸發 hooks 驗證"
