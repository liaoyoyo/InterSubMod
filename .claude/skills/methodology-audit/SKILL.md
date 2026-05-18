---
name: methodology-audit
description: 修改任何 C++ 方法前的方法學審查協議。審查現有方法的合理性、量化問題影響、產出決策文件供用戶選擇。觸發條件：「審查 XX 方法」、「評估是否要改 XX」、「XX 是否合理」、修改 .cpp/.hpp 前、發現統計失效率過高/分類邊界模糊/特徵設計可疑時。SKIP WHEN 純 comment / format / typo edits、純 docs 寫作、無 .cpp/.hpp/.h 觸及、已通過審查的方案進入實作階段（用 cpp-change）、純 build 驗證（用 verification-loop）。
user-invocable: true
---

# Methodology Audit Skill — PDD C++ 審查協議

## 觸發條件

- 用戶說「審查某方法」、「評估是否要改 XX」、「XX 是否合理」
- 即將修改任何 `.cpp` / `.hpp` 前
- 發現統計失效率過高、分類邊界模糊、或特徵設計可疑時

## 執行流程

### Step 1：IDENTIFY — 定位問題

```
目標：確認問題位置（檔案:行號）、影響範圍、相關統計指標
```

1. 使用 Grep/Glob 找到相關 `.cpp`/`.hpp` 檔案與對應行號
2. 確認問題影響的 pipeline 環節：
   - `passed_gating`（global chi-square）
   - `cluster_structure`（PERMANOVA cluster-first）
   - `label_structure`（PERMANOVA label-first）
   - `VerificationClass`（分類決策）
   - `heuristic_score` / `quality_score`（評分）

### Step 2：QUANTIFY — 量化現況

```
目標：從 label_first_metrics.tsv 或 significance.json 量化當前行為
```

使用 Python 腳本（或呼叫 `methodology-reviewer` subagent）：

```python
import pandas as pd

# 載入 7 個樣本的 label_first_metrics.tsv（路徑於各 round 的 step05_intersubmod）
# 計算：失效率、分佈、AUROC、class 比例
```

**必須量化的指標**（依問題類型選擇）：
- 若問題是「失效率」：計算 `*PermanovaValid=False` 的比例（分 TP/FP）
- 若問題是「分類邊界」：計算各 VerificationClass 的特徵分佈（中位數、四分位）
- 若問題是「特徵貢獻」：計算 AUROC（sklearn.metrics.roc_auc_score）
- 若問題是「閾值」：畫 precision-recall 曲線或計算 F1 delta per threshold

### Step 3：OPTIONS — 列出選項 + SWOT 2x2 matrix

列出至少 3 個選項（含「不改」選項）：

| 方案 | 名稱 | 摘要 | 修改位置 | 預估成本 | 預估 F1 影響 |
|------|------|------|---------|--------|------------|
| A | 不修改 | 接受現狀，補文件 | 無 | 低 | 0 |
| B | [具體名稱] | [一句話] | `src/core/XX.cpp:LL` | 中 | ±? |
| C | [備選] | [一句話] | `src/core/YY.cpp:LL` | 高 | ±? |

**SWOT 4 象限分析**（2026-05-18 P2 audit M4 強制 — 每個非 A 方案必填）：

對方案 B（同樣對 C, D...）：

|  | Helpful（有助達成目標）| Harmful（阻礙達成目標）|
|---|---|---|
| **Internal**（可控因子） | **S** Strengths：本方案內部優勢（如：實作簡單、可重用既有 module、test coverage 高）| **W** Weaknesses：本方案內部劣勢（如：增加 coupling、降低 testability、引入新依賴）|
| **External**（外部因子） | **O** Opportunities：外部機會（如：對齊業界 best practice、開啟下游研究方向、可發 paper）| **T** Threats：外部風險（如：未來 API breaking、SEQC2 標準變動、ONT 升版相容）|

**強制規則**：
- **每個非 A 方案至少填 1 個 S + 1 個 W + 1 個 O + 1 個 T**（4 象限不可空）
- 若 W + T > S + O 數量 → 自動降級為「條件性方案」，需 Step 2 補強量化才可推進
- A 方案（不改）的 SWOT 是「對照組」：S=穩定 / W=不解問題 / O=資源投他處 / T=問題擴大

### Step 4：WRITE — 輸出審查文件

輸出到 `docs/methodology/YYYYMMDD_{主題}_{流水號}.md`（使用下方模板）：

```markdown
<!--
建立時間: YYYY-MM-DD
問題類型: C++ 方法 | 統計方法 | 特徵設計
影響 track: TO | paired | 兩者
狀態: pending_decision
-->

# [方法名稱] 審查文件

## 問題描述
[清晰說明問題，引用具體程式碼位置]

## 量化影響
- 失效率/異常比例：XX%（來源：label_first_metrics.tsv）
- 受影響的樣本：[列表]
- 當前 F1 影響評估：[說明]

## 修改選項
### 方案 A：不修改（接受現狀）
### 方案 B：[具體修改名稱]
### 方案 C：[備選方案]

## 驗收標準
- [ ] test-quick 通過
- [ ] test-full F1 delta ≥ 0 或 [具體條件]
- [ ] 跨 3 個樣本確認無回退

## 用戶決策
**選擇**：[ ] A / [ ] B / [ ] C
**日期**：
**理由**：
```

### Step 5：DECISION — 等待用戶決策

展示審查文件後，等待用戶選擇 A/B/C/skip。

若用戶選定非 A 方案 → 記錄決策至文件，呼叫 `cpp-change` skill 開始實作。

---

## 關鍵數據路徑

```
# label_first_metrics.tsv（各樣本）
/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/{round}/step05_intersubmod/{tp|fp}/label_first_metrics.tsv

# significance.json（單位點）
{step05_intersubmod}/{tp|fp}/filtered_snv_{tp|fp}/{chr}/{anchor}/{window}/clustering/significance.json

# 跨樣本基線
python3 scripts/analysis/collect_baseline_metrics.py
```

## 相關工具

- `methodology-reviewer` subagent：深度分析 C++ 實作細節
- `cpp-change` skill：取得決策後的 6 步驟實作協議
- `scripts/analysis/classify_boundary_analysis.py`：VerificationClass 邊界分析

---

## 與 /scientific-rigor 元方法論的關係

本 skill 為 `/scientific-rigor §6 消融實驗設計` 的**具體實作**:
- Step 1 IDENTIFY / Step 2 QUANTIFY / Step 3 OPTIONS 對應 `/scientific-rigor §6 消融原則 4 步紀錄`（改了什麼 / 預期 / 實際 / 差異解讀）
- Step 3 OPTIONS 表（3+ 方案 + 成本估算）為 `/scientific-rigor §6 最小單元改動` 的實質執行格式
- Step 5 DECISION 必經 `/scientific-rigor §1 Hard Gate` 與 §11.6 雙環學習（質疑根本假設）對齊

**級聯觸發**: `/scientific-rigor §0.5 最小可用子集` 對「中影響」任務必跑 §6 → 觸發本 skill → Step 5 DECISION 後進 `/cpp-change` 6 步 PDD

---

## Phase Chain Position & Dependencies

- **Phase**: C++ change 前置（在 /cpp-change Step 0 之前）
- **Upstream**: 用戶提出 C++ 修改需求
- **Downstream**: `/cpp-change` 6-step（用戶選方案後）
- **Reads**: src/ / include/ / tests/ + benchmark baseline
- **Writes**: `InterSubMod/docs/methodology/YYYYMMDD_*.md`（含 Step 1-3 OPTIONS）

## Failure Mode & Diagnostics

| 症狀 | 可能原因 | 修法 |
|------|---------|------|
| 方案列出但未量化 effect size | Step 2 不完整 | 補 baseline F1 + expected delta + cohen's d |
| OPTIONS 只列 1-2 方案 | 未盡 §6 「列 ≥3 方案」 | 補替代方案 + SWOT 比較（P2 M4 fix） |
| 方案 A「不做」未顯式列 | 潛在偏好 bias | 始終列 A=「維持現狀」為對照 |

