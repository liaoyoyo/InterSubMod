---
name: observation-analysis
description: 標準化觀察分析腳本生成 — 產生符合 O1-O15 系列模式的 Python 分析腳本。USE WHEN：建立新觀察分析（O-系列）、產生 .py 分析腳本含數據載入/圖表/統計。涉及 scripts/analysis/*.py、docs/experiments/ 輸出。SKIP WHEN 既有 O-script 微調（直接 Edit）、純圖表參數調整、純 .md 寫作、純 build / commit、程式碼級驗證（用 verification-loop）。
allowed-tools: ["Read", "Write", "Bash", "Glob", "Grep"]
user-invocable: false
tags: ["analysis", "observation", "research"]
---

# 觀察分析腳本生成（Observation Analysis Builder）

依據 O1-O15 系列分析腳本的標準模式，生成新的觀察分析 Python 腳本。

## ⭐ 2026-05-26 新增 — Task Type + Scope Gate（A3 落地）

**腳本生成前必先做** task type classification（A-F per `AGENTS.md §15.3` / CLAUDE.md §0）：

1. **AskUserQuestion 必含 scope 維度**：
   > Q: 分析腳本 scope？
   > A1: 全基因組 + 全樣本（Recommended for B/C/D）— 預估時間長但證據完整
   > A2: Subset（chr1+chr8+chr19 / 1-2 樣本）— 必標 partial flag
   > A3: 自訂
2. **task type ∈ {B, C, D}** → 腳本 default 跑全 scope；用戶明示 subset 才縮小
3. **task type ∈ {A, E, F}** → subset OK 但腳本輸出 .md / .tsv 必加 `# PARTIAL_SCOPE: ...` header
4. **腳本 docstring 必含** `task_type` 與 `scope` 兩個元信息：
   ```python
   """[主題描述]
   task_type: B (Comprehensive validation)
   scope: full genome × 7 samples
   """
   ```

## 標準模式

所有分析腳本遵循以下結構：

### 1. 匯入與基礎設定

```python
#!/usr/bin/env python3
"""[主題描述]"""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from observation_common import (
    load_master_dataset, setup_plot_style, save_figure,
    compute_auc, compute_effect_size, compute_enrichment,
    bootstrap_ci, write_round_context, write_feature_statistics,
    ensure_dir, safe_div,
    SAMPLE_ORDER, MODE_ORDER, OUTPUT_ROOT,
    TRUTH_PALETTE, MODE_PALETTE,
    COLOR_TP, COLOR_FP, COLOR_PAIRED, COLOR_TO,
)

OUT_DIR = ensure_dir(OUTPUT_ROOT / "{YYYYMMDD}_{topic_slug}")
```

### 2. 數據載入

```python
df = load_master_dataset()  # 748K rows × 116 cols
# 依需求篩選
df_to = df[df['mode'] == 'to'].copy()
df_paired = df[df['mode'] == 'paired'].copy()
```

### 3. 圖表生成

每張圖遵循：
- `setup_plot_style()` 統一風格
- `save_figure(fig, OUT_DIR / "NN_description.png")` 統一保存
- 7 samples 使用 `SAMPLE_ORDER` 固定順序
- TP/FP 使用 `TRUTH_PALETTE`，Paired/TO 使用 `MODE_PALETTE`
- 圖表編號用 `{NN}_` 前綴（01, 02, ...）

### 4. 統計計算

可用工具：
- `compute_auc(y_true, y_score)` — ROC AUC
- `compute_effect_size(g1, g2)` — Cohen's d + bootstrap CI
- `compute_enrichment(tp, fp, tp_total, fp_total)` — enrichment ratio
- `bootstrap_ci(fn, data)` — bootstrap 信賴區間
- `scipy.stats.mannwhitneyu`, `scipy.stats.fisher_exact`, `scipy.stats.chi2_contingency`

### 5. 輸出

```python
# 統計摘要
write_feature_statistics(df, columns, OUT_DIR / "feature_statistics.tsv")

# 上下文記錄
write_round_context(
    OUT_DIR,
    observation_id="loh_concordance",
    title="LOH Definition Concordance",
    description="...",
    n_figures=N,
    n_tables=M,
)
```

## 輸入參數

生成腳本時需要：

| 參數 | 說明 | 範例 |
|------|------|------|
| topic_slug | 主題英文代稱 | `loh_concordance` |
| title | 中文標題 | `LOH 定義一致性分析` |
| description | 分析目標描述 | `比較 HP Imbalance 與 LOH.bed 的差異` |
| figure_specs | 圖表清單 | 見下方 |
| filter_mode | 數據篩選 | `to` / `paired` / `all` |
| additional_imports | 額外匯入 | `from sklearn.metrics import cohen_kappa_score` |

## 圖表規格格式

```python
FIGURE_SPECS = [
    {
        "id": "01",
        "name": "confusion_matrix",
        "title": "HP Imbalance vs LOH.bed Confusion Matrix",
        "type": "heatmap",  # heatmap / violin / bar / scatter / roc / distribution
        "stratify_by": "sample",  # sample / mode / quadrant / none
        "columns": ["Potential_LOH", "to_loh_bed_hit"],
    },
    # ...
]
```

## LOH.bed 路徑常數

如需使用 LOH.bed 路徑，直接從 O15b 複製：

```python
LOH_BED_PATHS = {
    "HCC1395": "/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260315_hcc1395_to_pilot/step03_longphase_to/tumor_phased_LOH.bed",
    "HCC1395_DORADO": "...20260315_hcc1395_dorado_to_pilot/step03_longphase_to/tumor_phased_LOH.bed",
    "COLO829": "...20260317_colo829_to_pilot_1/step03_longphase_to/tumor_phased_LOH.bed",
    "H1437": "...20260318_h1437_to_pilot_fastresume/step03_longphase_to/tumor_phased_LOH.bed",
    "H2009": "...20260318_h2009_to_pilot_fastresume/step03_longphase_to/tumor_phased_LOH.bed",
    "HCC1937": "...20260318_hcc1937_to_pilot_fastresume/step03_longphase_to/tumor_phased_LOH.bed",
    "HCC1954": "...20260318_hcc1954_to_pilot_fastresume/step03_longphase_to/tumor_phased_LOH.bed",
}
```

## 版本基準

- **LongPhase 版本**：所有分析以 baseline 為基準，PON-only 僅觀察
- **ISM Potential_LOH**：報告中改稱「HP Imbalance（HP 不平衡）」
- **Master dataset**：2026-03-27 HP fix 後重建版本

## 驗證標準

每個分析必須包含：
1. 數值一致性檢查（與已知基準比對）
2. 跨 7 samples 方向一致性
3. 微弱信號（AUC 0.55-0.65）必須排除 confound（O11-O13 教訓）

---

## Phase Chain Position & Dependencies

- **Phase**: P3 PILOT 子任務（產生 O-系列 Python 腳本，feed feature-layered-observation）
- **Upstream**: `/inject-hypothesis`（hypothesis_id 註冊後）
- **Downstream**: `/feature-layered-observation`（執行 Step 1-6）+ `/results-report`（產報告）
- **Writes**: `InterSubMod/scripts/analysis/O*.py` + `InterSubMod/docs/experiments/in_progress/YYYY/MM/<topic>.md`

## 與 /scientific-rigor 元方法論的關係

- **§4 因果 DAG**: O-script 產生 AUC/相關性前必畫 DAG（標 confounder / collider）→ 殘差過 collider 自動降級 characterization
- **§5 對照組 + 多方驗證**: 微弱信號 AUC 0.55-0.65 強制走 `/auc-confound-guard` 3-gate（within-group OLS / AF-bin / permutation）
- **§3 Effect Size**: O-script 報 AUC delta 時必標 Cohen's d / NNT / 臨床意義 ribbon

## Failure Mode & Diagnostics

| 症狀 | 可能原因 | 修法 |
|------|---------|------|
| AUC 跨 7 樣本方向不一致 | confound 或 sample-specific artifact | 跑 `/multi-sample-consistency` + 分層分析 |
| 殘差 OLS 後 AUC 升高 | collider bias (P-01) | 自動降級 verdict 為 "characterization only" |
| 微弱信號 chr+pos 聚合特徵 | spatial autocorrelation (P-09) | 必跑 mid-TP-rate window 驗證 |
