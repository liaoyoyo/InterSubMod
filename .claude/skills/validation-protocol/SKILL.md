---
name: validation-protocol
description: 假說驗證協議。根據假說類型自動建議 L1→L2→L3→L4 漸進驗證層級，每層含通過/失敗閾值與決策規則。觸發：「驗證假說」「validation protocol」「怎麼驗證」
allowed-tools: Read, Write, Glob, Grep
user-invocable: true
---

# Validation Protocol Skill

根據假說類型，自動建議漸進式驗證層級（L1→L2→L3→L4），每層有明確的通過/失敗閾值。

---

## 術語消歧義

本專案有四套分類系統，回答不同問題：

| 系統 | 代號 | 回答什麼 | 定義位置 |
|------|------|----------|---------|
| 修改層級 | **ML1/ML2/ML3** | 用什麼語言/方式修改程式碼 | `.claude/skills/research-loop/references/SCALE_LADDER.md` |
| 規模梯度 | **S1/S2/S3** | 用多少資料測試 | `.claude/skills/research-loop/references/SCALE_LADDER.md` |
| 驗證層級 | **L1/L2/L3/L4** | 信號是否可信（本文件） | `.claude/skills/validation-protocol/SKILL.md` |
| 證據卡層級 | **Tier 1/2/3** | 週報中報告多詳細 | `.claude/skills/weekly-report/references/LAYER_STRUCTURE.md` |

**容易混淆的組合**：
- ML3 ≠ L3：ML3 = C++ 修改，L3 = 跨樣本驗證
- S3 ≠ L3：S3 = 全量資料，L3 = 跨樣本驗證
- Tier 1（證據卡）≠ ML1（修改層級）：前者是報告格式，後者是程式語言選擇

---

## 觸發時機

- research-loop Step 3 DESIGN TEST 時自動觸發
- 用戶問「這個假說怎麼驗證」
- 新假說注入時建議驗證計劃

---

## 四層驗證梯度

### L1: AUC Screening（必做）

**目標**：快速確認信號是否存在。

| 指標 | 通過 | 失敗 | 灰色地帶 |
|------|------|------|----------|
| AUC | > 0.58 | < 0.52 | 0.52-0.58 |
| Delta F1 | > +0.001 | < -0.001 | -0.001 ~ +0.001 |

**灰色地帶處理**：記錄觀察，不做判斷，進入 L2 嘗試。

**數據集**：Pilot (HCC1395 only)，通過後擴展。

---

### L2: Confound Check（L1 通過後必做）

**目標**：確認信號不是已知 confounders 的假象。

**已知 confounders 清單**（依研究歷史）：
1. **AF (Allele Frequency)** — 最常見，幾乎所有特徵與 AF 相關
2. **n_reads / NumReads** — read count 影響統計量
3. **LOH status** — LOH region 行為完全不同
4. **HP tags** — haplotype assignment 可能帶入 self-phasing bias

**方法**：
- **Within-group OLS residualization**（不用 pooled OLS！參見 memory `feedback_pooled_ols_residualization_trap.md`）
- **AF-bin stratification**（L3 交叉驗證，防 L2 collider bias）
- 特徵分 Paired/TO 分開分析

| 指標 | 通過 | 失敗 |
|------|------|------|
| Residualized AUC | > 0.58 | <= 0.58 |
| AF-bin stratification | >=3/5 bins 同方向 | <3/5 bins |

---

### L3: Cross-sample Validation（L2 通過後）

**目標**：確認不是樣本特異性。

**固定順序**：HCC1395, HCC1395_DORADO, HCC1937, HCC1954, H2009, H1437, COLO829

| 指標 | 通過 | 失敗 |
|------|------|------|
| 方向一致 | >= 5/7 同方向 | < 5/7 |
| 最差樣本 AUC | > 0.55 | < 0.52 |

**注意**：H2009 和 H1437 行為可能與 HCC 系列不同（非 HCC 樣本）。

---

### L4: Beyond-AUC Validation（L3 通過後，目標為結論性判定）

**目標**：確認有生物學合理性，排除所有替代解釋。

檢查項目：
1. **Mechanism plausibility**：是否有生物學解釋？
2. **Literature support**：是否有論文支持此現象？
3. **Alternative explanations**：是否已排除所有替代解釋？
4. **Practical utility**：能否轉化為實際的過濾/分類規則？

| 結果 | 判定 |
|------|------|
| 所有 4 項通過 | **POSITIVE** — 結論性確認 |
| 1-3 通過但實用性差 | **CONDITIONAL** — 科學發現但非工具 |
| mechanism 不合理 | **NEGATIVE** — 即使 AUC 高也拒絕 |

---

### L4 升 ⭐4/⭐5 的硬性條件（2026-05-07 個人風格 anchor #1 硬化）

**規則**：cycle 升至 ⭐4 (validated) 或 ⭐5 (canonical) 之**前**，必須提交「4 軌獨立證據鏈」覆蓋表，且每軌至少 1 件 artifact 引用：

| # | 軌 | 定義 | artifact 範例 |
|---|---|---|---|
| (i) | **Statistical 統計** | AUC + p-value + CI（含 within-group residualization 後的數值） | pilot.json `metric_results` + auc-confound-guard report |
| (ii) | **Cross-sample 跨樣本** | ≥5/7 樣本方向一致 + 最差樣本 AUC > 0.55 | generalize.json `consistency` |
| (iii) | **Mechanism 機制/文獻** | 生物學解釋 + 至少 1 篇文獻支持（or 推理鏈到既有結論） | structured-tech-report §Background+Mechanism |
| (iv) | **Orthogonal 正交資料** | 來自不同 pipeline / archive / caller 的相同方向證據（避免 self-phasing 與 single-track artifact） | archive 對比 / replicate 跑 / 替代 caller |

**升 ⭐4/⭐5 PR 必附覆蓋表**（範本見下）；缺軌 → 降至 ⭐3 (described, single-track)，**不寫 evaluation.json 即觸 `pre_tier_upgrade_check.sh` hook 阻擋**。

**4-Track Coverage Table 範本**（直接複製到報告 §Evidence Layers 結尾）：

| Track | Status | Artifact Reference | Notes |
|---|---|---|---|
| (i) Statistical | ✅/❌ | `<path/to/auc-table.csv>` | AUC=...; residualized AUC=...; p=... |
| (ii) Cross-sample | ✅/❌ | `<path/to/generalize.json>` | n_passed/total = .../7; worst-sample AUC=... |
| (iii) Mechanism | ✅/❌ | `<path/to/report.md>` §X | 1 篇文獻 + 推理鏈 |
| (iv) Orthogonal | ✅/❌ | `<path/to/replicate>` | 替代 pipeline / archive / caller 結果 |

**為什麼是 4 軌而非 3 軌或 5 軌**：
- 3 軌缺 orthogonal → 易被 single-pipeline self-phasing 誤導（如 V3-Fixed haplotag 若沒 archive 對比就升 ⭐4 會犯 04-26 thread B 撤回）
- 5 軌會把「文獻」與「推理鏈」拆開，但本專案實務上常合一，5 軌反而拖延升級
- 4 軌是「最少獨立證據維度」+「經 04 月 6 個撤回事件回溯驗證」（plan v1.6 §4.5.4-A anchor #1）

**Connection to Harness P5 EVALUATE**：
- `/run-evaluator` 計算 6 components 中 `multi_sample_consistency` 是 (ii) cross-sample 的 surrogate
- `pitfall_coverage_score` 命中 P-07「single-track validated」時直接降為 0（即 (iv) 缺軌）
- 其餘 (i) (iii) 仍由人工填 4-track table 確認；Path B 將加 component 7 `multi_track_corroboration` 自動算覆蓋率

---

## 假說類型 → 驗證路徑建議

| 假說類型 | 建議路徑 | 說明 |
|----------|---------|------|
| `rule_change` | L1 → L2 → L3 | 不需 L4（規則修改非生物發現） |
| `feature_hypothesis` | L1 → L2 → L3 → L4 | 完整路徑 |
| `literature_feature` | L1 → L2 → L4 → L3 | 先驗機制再跨樣本 |
| `cpp_improvement` | L1 → L3（全量） | 效能改進直接全量測試 |
| `param_combo` | L1 → L1（多參數） → L3 | 參數搜索後直接跨樣本 |

---

## 決策規則

```
L1 通過? ──No──→ REJECT（記錄，下一個假說）
    │
   Yes
    ↓
L2 通過? ──No──→ CONFOUND（記錄 confound 來源，探索 residualization）
    │
   Yes
    ↓
L3 通過? ──No──→ SAMPLE-SPECIFIC（記錄，標記例外條件）
    │
   Yes
    ↓
L4 通過? ──No──→ CONDITIONAL（科學發現，非工具）
    │
   Yes
    ↓
  POSITIVE（確認，整合進 pipeline）
```

**每層失敗都需填寫 evidence_ledger + 更新 hypothesis_queue.json status。**

---

## Verdict 術語對照

本專案有兩層 verdict，各有其用途：

| 層次 | 欄位 | 可選值 | 寫入位置 | 定義文件 |
|------|------|--------|---------|---------|
| **單輪測試** | `verdict` | `positive_pilot` / `negative` / `dataset_specific` / `annotation_only` | evidence_ledger.jsonl | `.claude/skills/research-loop/references/TEST_DISCIPLINE.md` |
| **假說最終判定** | `Verdict` | `positive` / `negative` / `NO-GO` / `conditional` | Hypothesis Tracker | 本文件 + `.claude/skills/research-loop/references/HYPOTHESIS_TRACKER_TEMPLATE.md` |

**驗證決策 → 最終判定 對照**：

| 驗證決策 (本文件) | 對應最終判定 | 何時使用 |
|-------------------|-------------|---------|
| REJECT | negative 或 NO-GO | L1 失敗 → negative；多輪全 REJECT → NO-GO |
| CONFOUND | negative | L2 失敗，信號來自 confounders |
| SAMPLE-SPECIFIC | conditional | L3 失敗，僅部分樣本有效 |
| CONDITIONAL | conditional | L4 部分通過 |
| POSITIVE | positive | L4 全通過 |

---

## 結論生命週期 Checklist

當假說到達最終判定時，依判定類型更新以下檔案：

### POSITIVE / CONDITIONAL

1. `hypothesis_queue.json` — status → `concluded`
2. `evidence_ledger.jsonl` — 最終記錄（含 verdict）
3. Hypothesis Tracker 報告（使用 `HYPOTHESIS_TRACKER_TEMPLATE.md`）
4. Memory — 建立/更新 `project_*.md`（status: active）
5. `MEMORY.md` — 加入 Active Research 區段
6. `docs/CURRENT_FOCUS.md` — 如影響當前方向則更新

### NEGATIVE

1-3 同上
4. Memory — 建立 `project_*.md`（status: concluded）
5. `MEMORY.md` — 加入 Concluded 區段

### NO-GO（額外步驟）

1-5 同 NEGATIVE，再加：
6. NO-GO 報告（使用 `NOGO_REPORT_TEMPLATE.md`）
7. `scripts/hooks/research_direction_guard.sh` — **新增守門 regex**（見該檔案頂部擴展指南）
8. `docs/experiments/INDEX.md` — 更新實驗結論

---

## 擴展指南

### 新增 Confounder（L2 清單）

當發現新的系統性 confound 時：
1. 在本檔案 L2 的「已知 confounders 清單」中新增條目
2. 說明：名稱 + 一句話影響描述 + 發現來源
3. 如有對應的 memory，在條目後標註 memory 檔名

### 新增假說類型

在 `.claude/skills/research-loop/SKILL.md`「假說類型與優先分數」表中：
1. 新增類型行，含基礎分和加分條件
2. 在本檔案「假說類型 → 驗證路徑建議」表中新增對應路徑
3. 如需跳過某層（如 `rule_change` 不需 L4），明確說明原因

---

## 輸出

對用戶展示：
1. 建議的驗證路徑（基於假說類型）
2. 每層的具體通過/失敗閾值
3. 當前應執行的層級和預期結果

模板路徑：`.claude/skills/research-loop/references/HYPOTHESIS_TRACKER_TEMPLATE.md`
NO-GO 報告模板：`.claude/skills/research-loop/references/NOGO_REPORT_TEMPLATE.md`
