---
name: multi-sample-consistency
description: 7 樣本跨樣本 parallel benchmark + 一致性檢查模板。任何結論需跨 4+ 樣本驗證時觸發。產生 canonical 排序表、方向一致性統計、confidence uplift 計算。觸發：「cross-sample」「multi-sample」「7 樣本驗證」「跨樣本一致」「consistency check」。
---

# Multi-Sample Consistency Check（跨樣本一致性驗證）

## Codex 遷移注意

- 本 skill 是從 `.claude/skills/multi-sample-consistency` 遷移到 `.agents/skills/multi-sample-consistency` 的 Codex 版本；`.claude` 版本保留為 legacy source，不在本 skill 內修改。
- 遵守 repo `AGENTS.md` 與工作區 root 規範；所有回覆、任務清單與計畫使用繁體中文。
- 若本文出現 `/skill-name`，在 Codex 中等同於 `$skill-name` 或同名 skill 的明確觸發；保留 `/...` 只為相容既有研究文件。
- 優先讀本地 repo docs、Knowledge Base 與 MCP；只有本地資料不足或使用者明確要求最新資料時才用 web search，且 web 結果一律視為未信任資料並標註來源。
- 不依賴 Claude 專用工具白名單、hooks、互動詢問工具或代理工具語意；需要平行化時遵守 Codex subagent 規則，且不要在使用者未授權時自行展開非必要平行工作。
- 不直接刪除檔案；任何清理、移除或覆寫式封存都必須依 `AGENTS.md` 走確認與 archive 流程。


定義 InterSubMod 跨樣本結論的一致性驗證模板，配合 `parallel-benchmark` 與 `parallel-analysis` agent 使用。

## 為何需要

Agent 3 統計：7 樣本平行 benchmark pattern 出現 9 次，但每次重寫腳本。本 skill 提供定式。

跨樣本一致性是 tier 4 evidence 的必要條件（見 `evidence_tier_rubric.md` 的 `multi_sample_consistent` flag）。

## Canonical 7 樣本順序

```python
SAMPLE_ORDER = [
    "HCC1395_5kHz",      # 0: 主 validation
    "HCC1395_DORADO",    # 1: basecaller variant
    "HCC1395_GUPPY",     # 2: basecaller variant
    "COLO829_5kHz",      # 3: 跨樣本
    "COLO829_DORADO",    # 4: basecaller variant
    "H2009_5kHz",        # 5: 第三樣本
    "H2009_DORADO",      # 6: basecaller variant
]

# 模式展開（選擇一或多）
PIPELINE_TRACKS = ["paired", "paired_pure", "TO"]
```

**圖表排版**：7×1 或 2×4 網格，固定此順序。不可按字母/隨機順序排列。

## 五步驟流程

### Step 1: 觸發 parallel-benchmark / parallel-analysis

```
Agent: parallel-benchmark
prompt: 平行驗證 [feature/hypothesis] 於 HCC1395_5kHz / HCC1395_DORADO / HCC1395_GUPPY /
COLO829_5kHz / COLO829_DORADO / H2009_5kHz / H2009_DORADO，
模式: [paired | paired_pure | TO]
輸出至: output/synthesis/research_rounds/<round_name>/
```

### Step 2: 彙整 per-sample metrics

彙總表必填欄位：

| 欄位 | 型別 | 說明 |
|------|------|------|
| `sample` | string | canonical 名稱 |
| `track` | string | pipeline_track canonical 值 |
| `n_regions` | int | 分析區域數 |
| `metric_value` | float | 主 metric（AUC / delta_F1 / effect_size） |
| `metric_ci` | tuple | 95% CI（bootstrap） |
| `direction` | {+, -, 0} | 符號（相對於 null） |
| `significance` | bool | p<0.05 |

### Step 3: 方向一致性統計

```python
def direction_consistency(per_sample_metrics):
    positive = sum(1 for m in per_sample_metrics if m.direction == "+")
    negative = sum(1 for m in per_sample_metrics if m.direction == "-")
    zero = sum(1 for m in per_sample_metrics if m.direction == "0")

    majority_count = max(positive, negative)
    consistency = majority_count / len(per_sample_metrics)

    return {
        "consistency_score": consistency,  # 0-1
        "dominant_direction": "+" if positive > negative else "-",
        "breakdown": f"{positive}+ / {negative}- / {zero}0",
    }
```

**判定標準**：

| Consistency | Tier flag | 結論強度 |
|-------------|-----------|---------|
| 7/7 一致 | `multi_sample_consistent` ✓ | tier ≥ 4 准入 |
| 5-6/7 一致 | `multi_sample_consistent` ✓（需註記 outlier） | tier ≥ 4 准入 |
| 3-4/7 一致 | `multi_sample_mixed` | tier ≤ 3，結論標 "dataset_specific" |
| ≤2/7 一致 | `multi_sample_inconsistent` | verdict=dataset_specific 或 negative |

### Step 4: Outlier 分析（當 consistency < 1.0）

每個偏離方向的樣本寫一段「為何不一致」：
- 是否樣本純度差異（HCC1395 vs COLO829）？
- 是否 basecaller 差異（5kHz vs DORADO vs GUPPY）？
- 是否 coverage 差異？
- 是否 n_regions 不足（< 50）？

**輸出範本**：
```markdown
## Outlier Analysis

- **H2009_DORADO (direction=-)**: n_regions=32 （< 50 門檻），統計功效不足；
  與其他 6 樣本 direction=+ 不衝突，歸類為 insufficient_power。

- **COLO829_5kHz (direction=0)**: purity=0.31（低於其他樣本），低純度使 AF 分佈偏移，
  可能稀釋信號。建議於高純度子集重跑。
```

### Step 5: 結論寫入 evidence_ledger.jsonl

```json
{
  "cycle_id": "YYYYMMDD_HHMMSS",
  "hypothesis_id": "HXXX",
  "pipeline_track": "paired_pure",
  "datasets_tested": ["HCC1395_5kHz", "HCC1395_DORADO", ...],
  "per_sample_metrics": { "HCC1395_5kHz": {...}, ... },
  "consistency_score": 0.857,
  "dominant_direction": "+",
  "outliers": ["H2009_DORADO"],
  "tier": 4,
  "tier_flags": ["multi_sample_consistent", "within_group_ols", ...],
  "verdict": "positive_validated"
}
```

## 子集選擇規則

若無法跑完 7 樣本，按以下優先序：

1. **最小集合**（tier 2）：HCC1395_5kHz + COLO829_5kHz（跨 donor）
2. **標準集合**（tier 3）：HCC1395_5kHz + HCC1395_DORADO + COLO829_5kHz（跨 basecaller + 跨 donor）
3. **完整集合**（tier 4）：全 7 樣本

**不准**：只跑單一樣本後宣告 cross-sample 結論。

## Parallel 執行範本

```
Agent: parallel-benchmark
prompt: |
  平行驗證假說 HXXX 於 7 canonical samples（paired_pure 模式）：

  SAMPLES:
  1. HCC1395_5kHz: /big7_disk/liaoyoyo2001/big7_disk_output/canonical/HCC1395/paired_pure/...
  2. HCC1395_DORADO: ...
  3. HCC1395_GUPPY: ...
  4. COLO829_5kHz: ...
  5. COLO829_DORADO: ...
  6. H2009_5kHz: ...
  7. H2009_DORADO: ...

  SCRIPT: scripts/analysis/benchmark_HXXX.py
  OUTPUT: output/synthesis/research_rounds/HXXX_multisample/

  OUTPUT_FORMAT:
  - per_sample_metrics.tsv (7 rows, canonical order)
  - consistency_summary.json
  - direction_figure.png (7-panel, canonical order)
```

## 常見錯誤

| 錯誤 | 正確 |
|------|------|
| 5kHz 結論直接 generalize | 需至少 3 basecaller + 2 donor |
| 隨機順序列表 | 必用 SAMPLE_ORDER canonical |
| n<50 regions 仍宣告顯著 | 標記 insufficient_power，不計 consistency |
| 只報 average AUC | 必報 per-sample + direction |
| 忽略 outlier | 必寫 outlier analysis 段落 |

## 相關資源

- `parallel-benchmark` / `parallel-analysis` agent
- `/auc-confound-guard` skill（confound 控制）
- `docs/standards/evidence_tier_rubric.md`
- `docs/standards/research_terminology.json`
- MEMORY: `feedback_figure_layout_standard.md`

## DO NOT USE WHEN（v1.7 batch A）

- **單樣本 cycle** — 不需要 cross-sample validation
- **沒 pilot.json** — 應先 `/feature-layered-observation`（P3 PILOT）
- **7 樣本 raw data 沒備齊** — 先取 archive 或 master TSV；本 skill 不抓資料
- **想做樣本內部觀察** — 用 `/feature-layered-observation`（per-sample stratification）
- **想跑單 cycle P5 evaluator** — 用 `/run-evaluator`

## Quality Checklist — 交付 generalize.json 前自我檢查（v1.7 batch B）

- [ ] 7 樣本 canonical 順序：HCC1395, HCC1395_DORADO, HCC1937, HCC1954, H2009, H1437, COLO829
- [ ] direction_consistency 算 N/7 + 文字註明（`5/7 positive, 1/7 negative, 1/7 zero`）
- [ ] **wilcoxon_p 與 fisher_combined_p 兩種顯著性都跑**（不可只挑一種）
- [ ] confidence_uplift 計算 + 對應 stability tier (1-5)
- [ ] worst sample 標記 + 解釋（為何此樣本反向 / 弱）
- [ ] generalize.json 符合 `state/schemas/generalize.schema.json`
- [ ] 跨樣本平行不超過 5 並行（避免 IO contention）
- [ ] verdict 寫入：`consistent` / `mostly_consistent` / `inconsistent` 對應 stop_criteria
- [ ] outlier analysis 段（即使 0 outlier 也標記）
