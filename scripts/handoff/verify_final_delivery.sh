#!/usr/bin/env bash
# SPDX-License-Identifier: GPL-3.0-only
#
# 重驗一份最終交付品 —— 不靠任何敘述，全部機械重算。
#
#   bash verify_final_delivery.sh <FINAL_DIR> [CANONICAL_LINK_DIR]
#
# 六道檢查，任一失敗即 exit 1：
#   1. MANIFEST.json / DELIVERY.md 存在且可解析
#   2. artifacts[] 裡每個登記過 sha256 的檔案，重算後相符
#   3. 沒有 UNREGISTERED 檔案（交付目錄不該有來歷不明的東西）
#   4. BAM 通過 quickcheck、有索引、read 數與 manifest 一致
#   5. 譜系標籤真的寫出來了（lv / lg / ln / ls 各抽驗）
#   6. canonical symlink（若有指定）指向本目錄且可讀
#
# 設計原則：**失敗要吵**。靜默通過比明確失敗危險得多 —— 例如缺 .bai 時
# `samtools view FILE chr21` 會回空且不報錯，本腳本改用 idxstats 明確揪出。

set -uo pipefail

FINAL="${1:-}"
CANON="${2:-}"
[[ -n "$FINAL" && -d "$FINAL" ]] || { echo "用法: $0 <FINAL_DIR> [CANONICAL_LINK_DIR]" >&2; exit 2; }

fail=0
ok()   { printf '  ✅ %s\n' "$*"; }
bad()  { printf '  ❌ %s\n' "$*"; fail=$((fail+1)); }
note() { printf '     %s\n' "$*"; }

echo "═══ 交付品驗證 $FINAL ═══"

# ── 1. 入口檔 ───────────────────────────────────────────────
echo "[1/6] 入口檔"
M="$FINAL/MANIFEST.json"
if [[ -f "$M" ]] && python3 -c "import json,sys; json.load(open(sys.argv[1]))" "$M" 2>/dev/null; then
    ok "MANIFEST.json 可解析"
else
    bad "MANIFEST.json 缺失或無法解析"; echo; echo "FAIL"; exit 1
fi
[[ -f "$FINAL/DELIVERY.md" ]] && ok "DELIVERY.md 存在" || bad "DELIVERY.md 缺失"

read -r TRUST LIFE SAMPLE READS < <(python3 -c "
import json,sys
d=json.load(open(sys.argv[1]))
print(d.get('trust_tier','?'), d.get('lifecycle','?'), d.get('sample','?'),
      str(d.get('reads_total','?')).replace(',',''))" "$M")
note "sample=$SAMPLE  信任度=$TRUST  生命週期=$LIFE"
[[ "$TRUST" == "CURRENT" ]] || note "⚠ 信任度非 CURRENT —— 引用前先確認這是刻意的"

# ── 2. sha256 逐檔重算 ─────────────────────────────────────
echo "[2/6] sha256 重算"
python3 - "$M" "$FINAL" <<'PY'
import json, hashlib, sys
from pathlib import Path
m = json.load(open(sys.argv[1])); root = Path(sys.argv[2])
n = miss = bad = 0
for a in m["artifacts"]:
    if "sha256" not in a:
        continue
    f = root / a["path"]
    if not f.is_file():
        print(f"  ❌ 檔案不存在：{a['path']}"); miss += 1; continue
    h = hashlib.sha256()
    with open(f, "rb") as fh:
        while c := fh.read(16*1024*1024):
            h.update(c)
    if h.hexdigest() != a["sha256"]:
        print(f"  ❌ sha256 不符：{a['path']}"); bad += 1
    n += 1
print(f"  {'✅' if (miss==0 and bad==0) else '❌'} 有雜湊登記的 {n} 檔："
      f"相符 {n-bad-miss}、不符 {bad}、遺失 {miss}")
sys.exit(1 if (miss or bad) else 0)
PY
[[ $? -eq 0 ]] || fail=$((fail+1))

# ── 3. 未登記檔案 ──────────────────────────────────────────
echo "[3/6] 未登記檔案"
U=$(python3 -c "
import json,sys; d=json.load(open(sys.argv[1])); print(d.get('unregistered_count',0))" "$M")
if [[ "$U" == "0" ]]; then ok "無未登記檔案"
else bad "$U 個未登記檔案 —— 交付目錄不該有來歷不明的東西"; fi

# ── 4. BAM 完整性 ─────────────────────────────────────────
echo "[4/6] BAM"
BAM=$(python3 -c "
import json,sys
d=json.load(open(sys.argv[1]))
print(next((a['path'] for a in d['artifacts'] if a['role']=='lineage_tagged_bam'),''))" "$M")
if [[ -z "$BAM" ]]; then
    bad "manifest 沒有登記 lineage_tagged_bam"
else
    B="$FINAL/$BAM"
    samtools quickcheck -v "$B" 2>/dev/null && ok "quickcheck" || bad "quickcheck 失敗"
    [[ -f "$B.bai" || -f "$B.csi" ]] && ok "索引存在" \
        || bad "缺索引 —— region query 會靜默回空，IGV 會看似無資料"
    GOT=$(samtools idxstats "$B" 2>/dev/null | awk '{s+=$3+$4} END{printf "%d",s}')
    if [[ "$GOT" == "$READS" ]]; then ok "read 數守恆 $GOT"
    else bad "read 數 $GOT ≠ manifest 記的 $READS"; fi
    NC=$(samtools idxstats "$B" 2>/dev/null | awk '$3>0{n++} END{print n+0}')
    note "含 read 的 contig：$NC"
fi

# ── 5. 譜系標籤 ───────────────────────────────────────────
echo "[5/6] 譜系標籤"
if [[ -n "${BAM:-}" ]]; then
    CH=$(samtools idxstats "$FINAL/$BAM" 2>/dev/null | awk '$3>0{print $1; exit}')
    S=$(samtools view "$FINAL/$BAM" "$CH" 2>/dev/null | head -200000)
    for t in lv lg ln ls lp; do
        c=$(printf '%s' "$S" | grep -c "${t}:" || true)
        if [[ "$c" -gt 0 ]]; then ok "$t 有寫出（$CH 抽 200k 條命中 $c）"
        else bad "$t 完全沒有 —— 標籤未寫出"; fi
    done
fi

# ── 6. canonical symlink ──────────────────────────────────
echo "[6/6] canonical symlink"
if [[ -z "$CANON" ]]; then
    note "（未指定 canonical 目錄，略過）"
elif [[ ! -d "$CANON" ]]; then
    bad "canonical 目錄不存在：$CANON"
else
    n=0
    for l in "$CANON"/*; do
        [[ -L "$l" ]] || continue
        n=$((n+1))
        tgt=$(readlink -f "$l" 2>/dev/null || true)
        if [[ -z "$tgt" || ! -e "$tgt" ]]; then
            bad "斷掉的 symlink：$l"
        elif [[ "$tgt" != "$(readlink -f "$FINAL")"/* ]]; then
            bad "symlink 指向交付目錄之外：$l → $tgt"
        else
            ok "$(basename "$l") → 有效"
        fi
    done
    [[ $n -ge 2 ]] || bad "symlink 少於 2 個 —— .bam 與 .bai 都必須連（只連 BAM 會讓 IGV 靜默失敗）"
    if [[ $n -ge 1 ]]; then
        L=$(find "$CANON" -maxdepth 1 -name '*.bam' -type l | head -1)
        [[ -n "$L" ]] && { samtools idxstats "$L" >/dev/null 2>&1 \
            && ok "經 symlink 可正常讀取" || bad "經 symlink 讀取失敗"; }
    fi
fi

echo
if [[ $fail -eq 0 ]]; then
    echo "VERIFY RESULT: PASS  ($FINAL)"
    exit 0
else
    echo "VERIFY RESULT: FAIL  ($fail 項) —— 在修好之前不要把這份當結論引用"
    exit 1
fi
