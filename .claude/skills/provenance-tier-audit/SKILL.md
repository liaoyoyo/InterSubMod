---
name: provenance-tier-audit
description: 研究證據鏈審計與 tier 分佈報告。跨 hypothesis_queue / evidence_ledger / cycles / MEMORY.md 四處檢查一致性。產出週報素材、偵測 orphan cycles / over-claimed tiers / stale entries。觸發：「週報 evidence」「tier 分佈」「audit ledger」「證據鏈檢查」「結論盤點」。
allowed-tools: Read, Write, Bash, Grep
user-invocable: true
---

# Provenance Tier Audit（研究證據鏈審計）

跨 artifact 掃描研究證據鏈完整性，產出週報/landscape 可用素材，偵測結構性問題。

## 何時使用

- 撰寫研究週報前 → 取得「本週 tier 分佈」段落
- 撰寫 landscape 報告前 → 驗證結論穩定性分級
- 懷疑某結論被 over-claim 時
- 每週一次例行盤點（建議）

## 資料來源

| 來源 | 用途 |
|------|------|
| `research/autoresearch/evidence_ledger.jsonl` | 結論 truth source |
| `research/autoresearch/hypothesis_queue.json` | 待驗證假說 |
| `research/autoresearch/cycles/*/` | 執行證據（腳本、raw result） |
| `MEMORY.md` | 人類可讀結論索引 |

## 核心命令

```bash
# 完整 Markdown 報告
python3 scripts/analysis/audit_provenance_chain.py

# 週報素材（精簡版）
python3 scripts/analysis/audit_provenance_chain.py --weekly

# 機器可讀（供進一步處理）
python3 scripts/analysis/audit_provenance_chain.py --json

# 調整 staleness 閾值（預設 365 天）
python3 scripts/analysis/audit_provenance_chain.py --stale-days 180
```

## 五項檢查

### Check 1: Provenance Coverage
每筆 evidence_ledger 條目的 `artifacts_path` 是否可解析？
- **失敗情境**：artifacts 被刪除、搬移但 ledger 未更新
- **修正**：補寫 artifacts_path 或在 ledger 註記「artifacts migrated to `<new_path>`」

### Check 2: Orphan Cycles
`cycles/` 下的資料夾是否都有對應 ledger 條目？
- **失敗情境**：跑了 cycle 但沒寫入 ledger（結論被遺失）
- **修正**：補寫 ledger entry；或將 cycle 搬到 `cycles/archive/` 並註記放棄原因

### Check 3: Tier Distribution
統計 tier 1-5 的條目分佈；偵測 over-claimed（tier > confidence_cap）
- **失敗情境**：宣告 tier 3 但 tier_flags 只支撐 tier 2
- **修正**：補實驗達到宣告的 tier（不得直接降級覆寫歷史條目）

### Check 4: Staleness
tier ≤ 2 的 positive 條目是否超過 365 天未升級？
- **失敗情境**：早期 pilot 結論未驗證就滲透到結論區
- **修正**：排入升級實驗佇列，或改標 `annotation_only`

### Check 5: MEMORY.md ↔ Ledger
MEMORY.md Concluded 區段的 POSITIVE/NEGATIVE/NO-GO 數量與 ledger verdict 數量是否一致？
- **失敗情境**：ledger 新增但 MEMORY 未同步（對應 L1 hook 應已警告）
- **修正**：依 L1 hook 提示更新 MEMORY.md

## 週報段落範本

執行 `--weekly` 後可直接貼入週報：

```markdown
## 研究證據盤點（本週）

- **總條目**: 15
- Tier 1: 1 (7%)
- Tier missing: 14 (93%)

**需處理 issues: 2**
  - 1 個 artifacts 失效
  - 14 個條目缺 tier（需套用 tier rubric）
```

## Landscape 報告用法

全報告格式（預設）可直接引用至 `docs/reports/research_landscape/`：

```bash
python3 scripts/analysis/audit_provenance_chain.py \
  > docs/reports/research_landscape/07_Tier_Audit_$(date +%Y%m%d).md
```

## 建議執行頻率

| 場景 | 頻率 |
|------|------|
| 日常研究 | 無需主動跑 |
| 撰寫週報 | 每週一次 |
| 撰寫 landscape | 每次更新 |
| 用戶要求「盤點結論」 | 即時 |

## 與其他 skill 的關係

- `/auc-confound-guard` → 產出 tier_flags 後，本 skill 可驗證 tier 宣告
- `/multi-sample-consistency` → `multi_sample_consistent` flag 計入 tier 計算
- `/conclude-research` → 結論歸檔前應先跑本審計
- L1 `evidence_ledger_sync.sh` → 同步提醒；本 skill 做深度盤點

## 退出碼

- `0`: 審計完成（即使有 warnings）
- `2`: ledger 不存在或無條目

## 輸出範例（Markdown 格式）

```
# Provenance Chain Audit Report

**Generated**: 2026-04-18T...
**Ledger**: research/autoresearch/evidence_ledger.jsonl

## 1. Provenance Coverage
- Entries with valid artifacts_path: 14
- Entries with missing artifacts: 1

## 2. Orphan Cycles (no ledger entry)
None

## 3. Tier Distribution
  - Tier 1: 1
  - Tier missing: 14

## 4. Staleness (tier ≤ 2 positive, > 365 days)
None

## 5. MEMORY.md vs Ledger Cross-Check
  MEMORY.md markers: {'POSITIVE': 3, 'NEGATIVE': 7, ...}
  Ledger verdicts: {'positive_pilot': 3, 'negative': 8, ...}
```
