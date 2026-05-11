---
name: provenance-tier-audit
description: **跨 cycle 證據鏈一致性審計（system-level，非單 cycle 升級判定）**。掃描 hypothesis_queue / evidence_ledger / state/cycles / docs/experiments/INDEX / MEMORY.md 五處 artifact，偵測：(1) orphan cycle (無 ledger 對應)、(2) over-claimed tier (ledger stability < tier_used)、(3) stale entry (>365 天 active)、(4) 跨 artifact 不一致（如 INDEX 標 ⭐4 但 ledger verdict=NEGATIVE）、(5) tier 分佈異常（⭐4-5 過多）。**與 /run-evaluator 分工**：本 skill = 全域審計（週級別跑），/run-evaluator = 單 cycle P5 升 tier 前 retraction risk（cycle 級別跑）；兩者互補不重疊。觸發：「週報 evidence」「tier 分佈」「audit ledger」「證據鏈檢查」「結論盤點」「provenance 審計」「跨 cycle 一致性」。SKIP WHEN 單 cycle P5 tier 升級判定（用 run-evaluator）、單實驗結果分析（用 results-analysis）、查單一假說歷史（用 review-evidence）、純 build / commit / docs。
allowed-tools: Read, Write, Bash, Grep
user-invocable: true
version: 0.2.0
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

---

## Phase & Chain Position

- **Phase**: **Governance / Cross-phase**（非單一 phase；定期或事件驅動跨 cycle 審計）
- **Chain**: 服務 chain #9 (週報路徑) — `weekly-report` 觸發前先 `provenance-tier-audit` 提供素材
  ```
  weekly-report (cron weekly OR 「週報」觸發)
      ↓ (Layer 1 raw data 收集前)
  provenance-tier-audit ← (本 skill: 跨 5 處 artifact 一致性審計)
      ↓
  weekly-report 整合素材 → 17 段母稿
      ↓
  [選 A] pptx-build / [選 B] conclude-research
  ```
- **與 /run-evaluator 互補（非重疊）**：
  - `/run-evaluator`: **per-cycle** P5 升 tier 前 retraction risk（cycle scope, forward-looking）
  - `provenance-tier-audit`: **全域** 跨 cycle 一致性（system scope, backward-looking）
  - 兩者結合：cycle 內 evaluator 把關升級；全域 audit 把關歷史紀錄一致性
- **上游觸發**: weekly-report 流程開始 / 用戶手動「週報 evidence」「跨 cycle 審計」 / cron 週級別自動 / pivot-direction 後（檢查歷史是否誤判）
- **下游 skill**: `weekly-report`（消費素材）/ `memory-consolidation`（清理 stale entry）/ `pivot-direction`（基於 stale 標記決定是否轉向）

## Dependencies

| 類別 | 項目 |
|---|---|
| **Uses** (本 skill 內部呼叫) | Bash（python3 inline scripts: jsonl 解析 / glob 統計 / 跨 artifact diff）、Read（5 處 artifact 全讀）、Grep（搜 ⭐4-5 markers in INDEX.md） |
| **Used by** (誰會觸發本 skill) | `weekly-report` Layer 1 / 用戶手動「結論盤點」 / cron weekly hook / `pivot-direction` 換方向前的歷史檢查 |
| **Reads** | `research/autoresearch/hypothesis_queue.json`、`research/autoresearch/evidence_ledger.jsonl`（全 jsonl）、`state/cycles/*/state.json`（每 cycle）、`state/cycles_archived/*/state.json`、`docs/experiments/INDEX.md`（grep ⭐ markers）、`~/.claude/projects/*/memory/MEMORY.md` 索引 |
| **Writes** | **不直接寫永久檔案**；輸出文字報告（5 段：orphan / over-claim / stale / 跨 artifact 不一致 / tier 分佈異常）給呼叫者（weekly-report 或 stdout）。**例外**：偵測到嚴重 over-claim 時可 append 警告到 `state/invalidation/audit_warnings.jsonl`（新建，與 stale_marks.jsonl / tier_overrides.jsonl 並列） |

## Failure Mode & Diagnostics

| # | 失敗症狀 | 先看哪 | 排查步驟 |
|---|---|---|---|
| 1 | Orphan cycle (state/cycles/X/ 存在但 ledger 無 hypothesis_id 對應) | state/cycles/{id}/state.json `hypothesis_id` 欄位 vs ledger `hypothesis_id` jq 比對 | 補 ledger 條目（用 cycle artifacts 回推填寫）OR 移到 cycles_archived 並標 incomplete |
| 2 | Over-claimed tier（INDEX.md ⭐4 但 ledger stability="2"） | `evidence_ledger.jsonl` `stability` 字串解析 vs INDEX.md ⭐ markers | 個人風格 anchor #1「L4 多層必驗」要求；強制降級 INDEX 標記 OR 補跑 multi-sample-consistency 升 stability |
| 3 | Stale entry（last_relevant > 365 天，狀態仍 active） | MEMORY.md / state/active.json / hypothesis_queue.json `status: in_progress` + timestamp | 移到 cycles_archived；MEMORY 改成 Concluded；個人風格 anchor #7「pivot 容忍」允許但需明確 archive |
| 4 | 跨 artifact 不一致（INDEX 標 ⭐4 但 MEMORY.md 標 NEGATIVE） | grep ⭐ in INDEX.md vs MEMORY.md `## Concluded` 段 | 用戶決定 ground truth；通常 ledger > MEMORY > INDEX 優先順序 |
| 5 | Tier 分佈異常（⭐4-5 佔 >40% 全部 cycle） | jq count by tier in cycles + cycles_archived | 警告可能 over-claim 系統性；建議跑 evaluator 復驗 top-tier cycles |
| 6 | Cycle 缺 evaluation.json 但已標 ⭐4-5 | state/cycles/{id}/ ls + state.json `tier` | pre_tier_upgrade_check.sh hook 應已擋；歷史 cycle 補跑 /run-evaluator OR 加 tier-upgrade-override 注解 |
| 7 | Schema mismatch（state.json 用舊 schema_version） | state.json `schema_version` vs `state/schemas/state.schema.json` const | 寫遷移腳本；不直接改既有 cycle 的 state.json（plan §10 #1 規則） |

**何時升級到別的 skill / agent / 人工審查**：
- Over-claim ≥3 cycles 同時觸發 → 升級 `/run-evaluator` 對每個 cycle 重新評估
- Stale entry >10 個 → 升級 `memory-consolidation` skill 清理
- 跨 artifact 不一致涉及論文層級結論 → Hard Gate 暫停問用戶
- Tier 分佈異常（>40% ⭐4-5）→ 系統性問題；建議降低 evaluator 閾值或加 negative control（Drill 1 §6 改進建議）

**個人風格適配**（依 `feedback_*` memory）：
- **Anchor #1 「L4 多層驗證必建」** → over-claim 偵測必須交叉檢查 ledger.stability + multi-sample passed_count + 跨 artifact 一致性，不單看 INDEX 標籤
- **Anchor #7 「Pivot 容忍 + ⭐4-5 lock**」 → ⭐1-3 stale 寬容（reopen 容易），⭐4-5 stale 嚴查（lock 後撤回成本高）
- **Anchor #3 「報告 5 段骨架」** → 輸出格式對應 weekly-report Layer 1，便於直接餵入週報 §1 證據盤點段

## DO NOT USE WHEN（v1.7 batch A）

- **單 cycle 升 tier 判定** — 用 `/run-evaluator`（per-cycle retraction risk）；本 skill 是**全域**而非單 cycle
- **想跑 pilot 分析** — 用 `/feature-layered-observation`
- **寫週報** — 用 `/weekly-report`（本 skill 可作 Layer 1 input 但不取代）
- **修復某 cycle artifact** — 用該 cycle 對應 phase skill；本 skill 純 read-only audit
- **mass tier downgrade** — 必須逐 cycle 過 run-evaluator + reviewer override，不可批次改 INDEX

## Quality Checklist — 交付 audit 報告前自我檢查（v1.7 batch B）

- [ ] **5 處 artifact 全掃**：hypothesis_queue / evidence_ledger / state/cycles / experiments INDEX / MEMORY.md
- [ ] **5 種異常分類偵測**：(1) orphan cycle (2) over-claimed tier (3) stale entry (4) cross-artifact 不一致 (5) tier 分佈異常
- [ ] over-claim 交叉檢查 stability + passed_count + 跨 artifact 一致（anchor #1）
- [ ] **⭐4-5 stale 嚴查**（>30 天即警示），⭐1-3 寬容（>365 天）
- [ ] 輸出格式對應 weekly-report Layer 1（§1 證據盤點直接套）
- [ ] 每異常含：cycle_id + 證據 file path + 建議 action + severity
- [ ] **不修狀態，只報告**（read-only audit；任何 write 都是 bug）
- [ ] tier 分佈異常閾值：⭐4-5 占比 > 40% 提 systematic alarm
