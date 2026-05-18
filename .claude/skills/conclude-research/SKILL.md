---
name: conclude-research
description: 研究專案收尾與知識傳播。驗證 artifacts 完整性 → 生成報告骨架 → 自動更新 INDEX/MEMORY/evidence_ledger/CURRENT_FOCUS。觸發詞：conclude-research、結束研究、撰寫結論、P5→P6 transition、cycle 達到 ⭐4-5 verdict 收尾。SKIP WHEN cycle 仍在 P1-P4 探索階段（未到 verdict）、單實驗結果報告（用 results-report）、單一工程修改紀錄（用 structured-tech-report）、週進度多主題彙整（用 weekly-report）、AI 對話 session 紀錄（用 report）。
allowed-tools: Read, Write, Edit, Bash, Glob, Grep
user-invocable: true
---

# 研究收尾與知識傳播 (Conclude Research)

完成研究後的標準化收尾流程：驗證產出物 → 撰寫報告 → 傳播結論到所有追蹤系統。

## 觸發時機

`/conclude-research {project_name} {conclusion}`

- `project_name`：research/ 下的目錄名
- `conclusion`：`POSITIVE` | `NEGATIVE` | `CONDITIONAL`

## 執行流程

### Phase A：CONCLUDE（收尾）

#### Step 1：讀取 manifest.yaml

```bash
MANIFEST="research/${project_name}/manifest.yaml"
if [ ! -f "$MANIFEST" ]; then
    echo "ERROR: manifest.yaml not found. Run /init-research first or create manually."
    exit 1
fi
```

讀取 manifest 取得 artifacts 清單、假說列表、數據來源。

#### Step 2：Artifact 完整性驗證

對 manifest 中列出的所有 artifacts 執行檢查：

```python
checks = []

# 1. Scripts exist
for script in manifest['artifacts']['scripts']:
    path = f"research/{project_name}/{script}"
    checks.append(("script_exists", path, os.path.exists(path)))

# 2. Figures exist and not empty (>1KB)
for fig in manifest['artifacts']['figures']:
    path = f"research/{project_name}/{fig}"
    exists = os.path.exists(path)
    size_ok = os.path.getsize(path) > 1024 if exists else False
    checks.append(("figure_valid", path, exists and size_ok))

# 3. Data files exist
for data in manifest['artifacts']['data']:
    path = f"research/{project_name}/{data}"
    checks.append(("data_exists", path, os.path.exists(path)))

# 4. Plan exists
plan_path = f"research/{project_name}/00_PLAN.md"
checks.append(("plan_exists", plan_path, os.path.exists(plan_path)))
```

**失敗處理**：
- 任何 check 失敗 → 列出缺失項目，詢問用戶是否繼續
- 全部通過 → 自動繼續

#### Step 3：結論判定確認

向用戶展示結論摘要：

```
研究結論：{conclusion}

假說驗證結果：
  H1: {verdict} — {one-line summary}
  H2: {verdict} — {one-line summary}

是否確認此結論？ [Y/n]
```

#### Step 4：撰寫正式報告

在 `docs/experiments/validated/YYYY/MM/` 或 `in_progress/` 建立報告。

報告必須包含（參考 `results-report` skill 格式）：
- 背景與動機
- 方法（含篩選條件、公式、變數定義 — 足以完全複現）
- 結果（含統計數值、圖表引用）
- 結論判定（POSITIVE/NEGATIVE/CONDITIONAL）
- Limitations
- Next steps

#### Step 5：更新 manifest.yaml

```yaml
project:
  status: "completed"
  conclusion: "{conclusion}"

report:
  path: "docs/experiments/validated/YYYY/MM/{report_filename}"
  index_entry: true
```

---

### Phase B：PROPAGATE（知識傳播）

以下步驟全部自動執行，無需用戶確認。

#### Step 6：更新 INDEX.md

在 `docs/experiments/INDEX.md` 的適當位置加入新條目：

```markdown
### {研究標題} — {結論標記} {conclusion} {YYYY-MM-DD}

- **報告**: [{report_filename}](validated/YYYY/MM/{report_filename})
- **結論**: {one-line conclusion}
- **核心數據**: {key metrics}
- **建議後續**: {next steps or "已完成"}
```

結論標記：✅ = POSITIVE、❌ = NEGATIVE、⚠ = CONDITIONAL

#### Step 7：更新 MEMORY

在 memory 目錄建立或更新對應的 memory 檔案：

```markdown
---
name: {project title} {conclusion}
description: {one-line description for relevance matching}
type: project
status: active
last_relevant: {today}
---

{conclusion summary}

**Why:** {motivation}

**How to apply:**
- {key findings}
- {implications for future work}

**報告**: {report path}
**數據**: {data directory path}
```

更新 `MEMORY.md` index。

#### Step 8：更新 evidence_ledger.jsonl

追加結論 record（支援兩種類型）：

**Type 1: hypothesis_test（F1 驗證型）**
```json
{
  "cycle_id": "{YYYYMMDD}_{project_name}",
  "type": "hypothesis_test",
  "hypothesis_id": "{H_ID}",
  "delta_f1": {...},
  "verdict": "{verdict}"
}
```

**Type 2: characterization（特徵化研究型）** ← 新增
```json
{
  "cycle_id": "{YYYYMMDD}_{project_name}_{step}",
  "type": "characterization",
  "hypothesis": "{hypothesis statement}",
  "method": "{statistical method}",
  "result": {
    "effect_size": "{value}",
    "p_value": "{value}",
    "cross_sample": "{N/7}",
    "confound_control": "{method}"
  },
  "verdict": "{verdict}",
  "project": "{project_name}",
  "report_path": "{report path}"
}
```

#### Step 9：更新 CURRENT_FOCUS.md（條件性）

若研究結論影響當前工作方向：
- POSITIVE → 加入「新發現」區段
- NEGATIVE → 加入「已排除方向」提醒
- 否則不動

#### Step 10：完成確認

```
研究收尾完成：

✓ Artifacts 驗證通過 ({N}/{N} checks)
✓ 報告已寫入: {report_path}
✓ INDEX.md 已更新
✓ Memory 已建立: {memory_file}
✓ evidence_ledger 已追加 {M} 條 records
{✓|—} CURRENT_FOCUS.md {已更新|未需更新}
```

## DECLARE 7-letter checklist — 收尾前最後一道 cognitive forcing（v1.8 T2-5）

> 設計來源：DECLARE / Croskerry cognitive forcing strategies（PMC10149772 / PMC3786644）；改寫為 InterSubMod 研究結論語境。**Phase A Step 4「撰寫正式報告」之後、Phase B Step 6「更新 INDEX.md」之前必跑**。

| 字母 | 提問 | 對 InterSubMod 怎麼答 |
|---|---|---|
| **D** ifferential | 「除了我提的這個結論，至少還有哪 2 種競爭假說可以解釋同樣數據？」 | 列 ≥2 個 alternative explanations，且明確說明本 cycle 為何排除（非「沒想過」）|
| **E** vidence | 「L1-L4 4 軌證據齊全嗎？」 | Statistical / Cross-sample / Mechanism / Orthogonal — ⭐4-5 cycle 必須四軌齊（anchor #1 L4 mandatory）|
| **C** onfounders | 「n_reads / AF / LOH / spatial / OLS bias 5 大 confound 都已 sweep 過？」 | plan.confound_sweep_plan + pilot.confound_guard.* 對照；任一缺 → 補做或降 tier |
| **L** ikelihood | 「tier_used 與 evidence stability_grade 對齊？」 | 不可 over-claim（ledger stability=B 但 tier_used=⭐4 → over-claim） |
| **A** lternatives | 「失敗 / NEGATIVE 結論已對照嗎？」 | 至少引用 1 個歷史 NEGATIVE case（從 review-evidence 拉）做為對比；不可只列正面 |
| **R** eassess | 「跨 cycle 結論與本 cycle 矛盾嗎？」 | 與 INDEX.md / MEMORY.md 既有 ⭐4-5 結論 cross-check；衝突 → 暫停升 tier，啟動 provenance-tier-audit |
| **E** ngage | 「結論已主動 propagate 進 INDEX/MEMORY/CURRENT_FOCUS？」 | Phase B Step 6-9 完成；不留結論孤兒（orphan cycle）|

**全自動模式**：DECLARE 7 字逐項 self-check；任一不過 → 暫停回報用戶；7 字全過 → 進 Phase B propagate。
**互動模式**：列 7 字答覆給用戶，用戶 ack 後進 propagate。

## 注意事項

1. **不重寫報告**：若報告已存在，只更新 manifest 和 propagate；不覆蓋已有報告
2. **結論不可逆**：一旦標記 completed，不應改回 executing（如需重開，建新專案）
3. **NEGATIVE 同樣重要**：失敗結論也必須完整記錄，防止未來重複
4. **手動已完成的也可用**：即使研究是手動完成的（如 LOH-AF），也可用此 skill 補做 PROPAGATE
5. **DECLARE 不取代 evidence**：checklist 是 cognitive forcing；不能替代 L4 多軌證據鏈本身（PMID 24423999 有警示 CFS 對新手 randomized trial 不顯效）

---

## Phase Chain Position & Dependencies

- **Phase**: P6 CONCLUDE（cycle 終點）
- **Upstream**: `/run-evaluator` P5（tier 已升至 ⭐4-5）
- **Downstream**: `/memory-consolidation`（寫 MEMORY.md）+ `/weekly-report`（report PI）
- **Reads**: cycle artifacts + INDEX.md + evidence_ledger.jsonl
- **Writes**: MEMORY.md Concluded 新增條目 + INDEX.md status flip + CURRENT_FOCUS 更新

## 與 /scientific-rigor 元方法論的關係

- **§9.2 SRE Postmortem**: NEGATIVE 結論強制走 `templates/postmortem.md` 8 段範本
- **§8 任務壓縮**: 寫 MEMORY.md Concluded = semantic memory consolidation
- **§9 PDCA Act 階段**: 本 skill 為 cycle 的 Act 收尾（升級 / 修正 / 廢棄判定）

## Failure Mode & Diagnostics

| 症狀 | 可能原因 | 修法 |
|------|---------|------|
| Conclude 但缺多樣本驗證 | 未跑 multi-sample-consistency | 回 P4 補 7-sample cross-check |
| ⭐4-5 但 ledger stability 低 | over-claim — `/provenance-tier-audit` 會抓 | 降級或補實驗 |
| Concluded 條目超 200 行 cap | MEMORY.md 過載 | archive 至 `memory/archive/`（P4 fix） |

