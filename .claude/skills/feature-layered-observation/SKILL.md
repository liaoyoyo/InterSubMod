---
name: feature-layered-observation
description: **P3 PILOT 主要分析 skill** — 標準化 TP/FP 特徵分層觀察協議。對任一 ISM/VCF/BAM-QC 特徵或特徵群組執行 Step 1-6（global AUC → LOH×AF×CN 32-cell → 7-sample 一致性 → 分層 AUC → confound guard → spatial autocorr），產出 10 章節 .md + 6 類圖 + verdict + pilot.json（schema 對齊 `state/schemas/pilot.schema.json`）供 P4 multi-sample-consistency 接手。封裝 Phase A-E G1-G10 方法學。觸發：「觀察特徵」「feature layered」「feature × LOH×AF×CN」「TP FP 分層」「跨樣本一致性」「單一特徵觀察」「特徵群組觀察」「run Step 1-6」「P3 PILOT」。USE WHEN：新增功能群組、對現有特徵重跑、延伸到新 pipeline 的特徵驗證。
allowed-tools: ["Read", "Write", "Edit", "Bash", "Grep", "Glob"]
user-invocable: true
version: 0.2.0
tags: ["feature", "observation", "tp-fp", "auc", "loh", "stratified", "confound", "P3-PILOT"]
---

# Feature Layered Observation（特徵分層觀察協議）

封裝 Phase A-E (2026-04-23 完成) 對 137 個 ISM+VCF 特徵執行過的標準觀察流程，讓未來可一鍵對任意新特徵或群組重跑完整 Step 1-6 + verdict + 10 章節報告。

## 為何需要

G1-G10 十個功能群組驗證了同一條流程穩定可行（見 `research/feature_layered_observation/00_main_observation.md`）；未來：
- 新增 G11+ 群組（例：normal BAM 特徵、cross-region 特徵）時，避免重發明輪子
- 既有特徵在新 pipeline/新樣本下的重驗
- 單特徵 deep-dive（如 `HPFineN_HP1S` 後續 flag=on 重跑）

每次重建流程 = 約 1.5-3h (single feature) / 6-10h (full group 10+ features)；skill 化後 = 一鍵 spec + auto-layout。

## 輸入模式

| 模式 | 說明 | 參數 | 輸出 |
|------|------|------|------|
| **single** | 對單一特徵跑完整 Step 1-6 | `feature=<name>`，可選 `sample=<subset>`、`mode=<paired\|to\|all>` | 1 個 .md (10 章節) + 6 圖 + AUC/cell_delta/confound tsv |
| **group** | 對一組功能相關特徵（2-15 個）跑完整流程 | `group_id=G11`、`features=[f1,f2,...]`、`source=<c++ file pattern>` | 1 個 group .md + features/ 每特徵子章節 + 群組彙總圖 |
| **all** | 對 registry 中標記 `pending` 的所有特徵跑 | `filter=tier:raw/smoke` | 每特徵單獨 md + global index 更新 |
| **orthogonal** | 同時觀察 2-4 個互補特徵（跨群組）找 interaction | `features=[NumReads, HPMergedDelta, SampleASM_Delta]` | 1 個 orthogonal .md + pairwise scatter + residual AUC |

**預設模式**：未指定時 = `single`。

## 資料源契約

- **主資料**：`research/feature_layered_observation/data/merged_with_vcf.tsv.gz`（748,676 rows × 60 cols）
- **Registry**：`research/feature_layered_observation/data/feature_registry.tsv`（137 features；欄位定義見 `references/feature_registry_template.tsv`）
- **Extended master（per-group re-join）**：`data/G{N}/master_g{N}.tsv.gz`（當特徵需要 extra joins 如 HP_Ratio_norm、SampleASM_Delta）
- **⚠ 陷阱**：master 欄位 `AF` = `|AlleleDelta|`（非 caller VAF）— 下游切分必用 `vcf_AF` 作 caller VAF。若使用 `AF_class` 欄位，已預先用 `vcf_AF` 切好（Extreme/Near-half/Intermediate）。

## 執行流程（Step 1-6 + skill 分派）

### Step 0 · 前置檢查（必做）

1. 確認 feature 存在於 `merged_with_vcf.tsv.gz` 或指定 master；不存在 → 暫停並詢問來源 pipeline
2. 查 `feature_registry.tsv` 是否已有 `prior_conclusion` — 若已有 NEGATIVE/CONFOUND_COLLAPSED → **一行告知後繼續**或用戶明示 skip
3. 讀 `docs/reports/research_landscape/03_ISM分析價值界定.md` 確認該特徵未被標為 `annotation_only`（若已標：降為 characterization-only observation）
4. 呼叫 `/known-pitfalls`：涉及 OLS/residualization/AF collider 場景必讀

### Step 1 · 全域分佈

- **圖**：`fig01_global_distribution.png`（violin + KDE overlay）
- **統計**：mean±std、median、Cohen's d、Mann-Whitney U p、AUC + Wilson 95% CI
- **輸出表**：`{group}_global_stats.tsv`（欄位：feature, auc, lo, hi, cohen_d, mwu_p, mean_tp, mean_fp, n_tp, n_fp）
- **判定分流**：
  - AUC ≥ 0.58 → 進 Step 5 confound guard
  - AUC < 0.53 → 標 NEGATIVE，Step 2-6 簡化（僅 Step 2 分層 heatmap）
  - 0.53 ≤ AUC < 0.58 → 邊緣，跑完 Step 2-6 但 verdict 預設 NEGATIVE

### Step 2 · LOH × AF × CN 32-cell heatmap

- **切分**：`LOH_Subtype(5) × AF_class(3 from vcf_AF) × cn_tier_F(5 from CovM)` = 75 cells；有效（n≥20）約 32
- **四張圖**：`fig02_{tp_rate, feat_tp, feat_fp, delta}.png`
  - A: per-cell TP rate（色階 0-1）
  - B: feature mean for TP only
  - C: feature mean for FP only
  - D: Δ(TP−FP) + Cohen's d（signed divergent 色階）
- **每 cell 標註 n**（≥20 顯示數字，<20 顯示 `·`）

### Step 3 · 跨樣本一致性

- 7 samples × 核心 stratum grid heatmap（固定順序見下）
- Spearman concordance matrix（7×7）on per-cell TP rate 或 feature mean
- **圖**：`fig03_per_sample_consistency.png`
- **關鍵門檻**：median ρ ≥ 0.5 → 一致；同向 ≥5/7 樣本 → cross-sample pass
- **呼叫**：`/multi-sample-consistency`（當結論依賴 cross-sample 時強制觸發）

### Step 4 · 分層 AUC

- Bar chart：Global / per-LOH(5) / per-AF(3) / per-CN(5)
- 兩條虛線：0.50 (random)、0.58 (Beyond-AUC ceiling)
- **圖**：`fig04_stratified_auc.png`
- **表**：`{group}_auc_table.tsv`

### Step 5 · Confound guard（auc-confound-guard 整合）

- **前置**：AUC ≥ 0.58 才需要此步；否則輸出 `{group}_confound.tsv` 標 `skipped_reason=low_raw_auc`
- **必做**（整合 `/auc-confound-guard` 三關）：
  - Gate 1: within-group OLS 殘差化 on `NumReads + vcf_AF + Coverage_Multiple`（G10 加 AlleleDelta）
  - Gate 2: AF-bin 交叉驗證（Extreme/Near-half/Intermediate 各自計 AUC）
  - Gate 3: permutation test (1000 reps, p<0.05)
- **圖**：`fig05_confound_residualized.png`（raw vs resid AUC bar + AF-bin AUC trace）
- **輸出**：`{group}_confound.tsv`（feature, raw_auc, resid_auc, af_bin_aucs dict, perm_p）
- **呼叫**：`/auc-confound-guard`（必呼叫，無例外）

### Step 6 · Spatial autocorrelation

- chr + pos 5Mb bin 聚合，per-bin AUC
- **圖**：`fig06_spatial_autocorrelation.png`（genome-wide track + high-TP-rate-bin overlay）
- **Artifact warning**：若 AUC 只在 baseline TP rate > 80% 的 bin 出現 → spatial artifact，verdict 降級
- **參考**：memory `feedback_spatial_autocorrelation_confound`

## Agent 分派拓撲（Opus 4.7 subagent 觸發規則）

預設單回合完成；下列情境**才明確觸發** subagent：

| 情境 | Agent 模式 | 分派規則 |
|------|-----------|---------|
| `mode=all` 跑 10+ features | parallel per-feature | spawn 1 agent / 3-5 features（避免 overhead） |
| `mode=group` + features≥8 | parallel per-feature inside group | spawn 1 agent / 4 features |
| Step 3 per-sample consistency | sequential（single agent）| 不分派，因 Spearman matrix 需整體 |
| Step 5 AF-bin + permutation | sequential | 不分派 |
| Cross-group orthogonal analysis | parallel per-pair | 若 ≥6 pairs 才分派 |

**明示觸發語範例**：「spawn parallel feature-layered-observation for [f1, f2, f3, f4] under G11」。

## 整合 Skill 呼叫點

| Skill | 何時呼叫 | 為何 |
|-------|---------|------|
| `/auc-confound-guard` | Step 5 | 三關驗證；raw AUC≥0.58 強制呼叫 |
| `/multi-sample-consistency` | Step 3 | 跨樣本 canonical 排序、一致性統計 |
| `/known-pitfalls` | Step 0 + Step 5 | OLS/AF collider/n_reads 警告必讀 |
| `/citation-verification` | 第 9 章（論文背景） | bibliography 寫入前必驗 |
| `/observation-analysis` | 生成 Python 腳本時 | 已有標準模板 |
| `/grill-me` | 第 10 章（結論與質疑） | 3 個質疑題 |

## 標準輸出格式（feature-level 10 章節 .md）

骨架見 `references/observation_report_template.md`。章節固定：

1. **特徵定義與來源**（C++ file:line + 計算公式 + upstream deps）
2. **觀察目標**（1-2 段 hypothesis + prior expectation）
3. **全域分佈**（Step 1 stats table + fig01）
4. **LOH×AF×CN 分層**（Step 2 32-cell highlights + fig02×4）
5. **跨樣本一致性**（Step 3 concordance + fig03）
6. **分層 AUC**（Step 4 bar + fig04）
7. **Confound 檢查**（Step 5 residualized AUC + AF-bin + fig05）
8. **Spatial autocorrelation**（Step 6 warning or pass + fig06）
9. **論文與知識庫背景**（knowledge MCP search + bibliography，至少 3 refs）
10. **結論與質疑**（verdict + 邏輯鏈 + 3 質疑 + 後續建議）

## 圖表規範

- **標題**：一句話目標 + sample scope + mode（例：`NumReads TP/FP discrimination across 7 samples (paired + to)`）
- **座標軸**：帶單位（counts / rate / log10 / Δ）
- **圖例**：TP = 藍（`#1f77b4`），FP = 橘（`#ff7f0e`）；Paired/TO 見 `MODE_PALETTE`
- **Heatmap**：colorbar + threshold line（0.58、0.50）；divergent 用 TwoSlopeNorm(0 center)
- **n 標註**：每 cell / bar；`n<20` 顯示 `·`，否則顯示整數
- **字型**：Latin + CJK fallback（Noto Sans CJK）避免方塊
- **固定樣本順序**：`HCC1395, HCC1395_DORADO, HCC1937, HCC1954, H2009, H1437, COLO829`
- **LOH 順序**：`LOH_None, LOH_Weak, LOH_Noise, LOH_Strong, LOH_Subclone`
- **AF 順序**：`Extreme, Near-half, Intermediate`
- **CN tier 順序**：`CN_Loss, CN_Near1, CN_Diploid, CN_Gain, CN_HighGain`
- **DPI**：150；`bbox_inches='tight'`；等比縮放不擠壓

## Verdict 規則（決策樹見 `references/verdict_rubric.md`）

| Raw AUC | Resid AUC | Cross-sample (≥5/7) | Spatial | Verdict |
|---------|-----------|---------------------|---------|---------|
| ≥0.65 | ≥0.55 | Yes | non-artifact | **POSITIVE** |
| 0.58–0.65 | ≥0.55 | Yes | non-artifact | **CONDITIONAL_POSITIVE** |
| ≥0.58 | <0.55 | any | any | **CONFOUND_COLLAPSED** |
| ≥0.58 | ≥0.55 | <5/7 | any | **SAMPLE_SPECIFIC** |
| 0.50–0.58 | any | any | any | **NEGATIVE** |
| <0.50 | any | any | any | **ANTI_SIGNAL**（標註 direction flipped） |
| any | any | any | artifact only | **SPATIAL_ARTIFACT** |

**Tier 寫入 evidence_ledger.jsonl**：POSITIVE=tier 3-4, CONDITIONAL_POSITIVE=tier 3, 其餘 ≤ tier 2（characterization_only / annotation_only）。

## 知名陷阱清單（觀察此領域前必看）

1. **L2 collider bias**（AlleleDelta ↔ AF）— memory `feedback_L2_collider_bias`；殘差化 AF 會產生虛假反向信號
2. **Pooled OLS residualization trap** — memory `feedback_pooled_ols_residualization_trap`；必用 within-group OLS
3. **n_reads confound**（O11 epipolymorphism）— AUC<0.60 微弱信號近全為此所致
4. **Sample-specific fingerprint**（HP tag driven）— HCC1395 特異，其他樣本散亂 → SAMPLE_SPECIFIC 非 generic
5. **Pipeline-dependent**（paired vs TO）— paired 是 HP-rich，TO 是 HP-self-phased；同特徵跨 mode 不可比
6. **Spatial autocorrelation artifact** — chr + pos 聚合特徵必 mid-TP-rate 驗證
7. **Diploid_Coverage_Used 是 assay-level**（G1 教訓）— 全樣本同值，AUC 0.79 是 "which sample?" shortcut
8. **Master `AF` 欄位 ≠ caller VAF** — master `AF` = `|AlleleDelta|`；要用 `vcf_AF`

## 質疑模板（第 10 章 grill-me 3 問）

對每個 verdict POSITIVE / CONDITIONAL_POSITIVE 的結論，必問：

1. **「這個信號在 NumReads + vcf_AF + Coverage_Multiple 三個殘差後還有多少？若 <0.05 就只是這三者的 proxy。」** → 強制 Step 5 Gate 1 數字
2. **「Cross-sample 一致性在 high-TP-rate bin 還成立嗎？去掉 HCC1395 (TP 92%) 後呢？」** → mid-TP-rate window 驗證，排除 spatial artifact
3. **「這個特徵在 paired 與 TO 模式下是否同向？若只在 paired 有效，是否為 HP-tag 人工產物？」** → pipeline-dependent 檢查（memory `project_readparser_germline_hp_only_phase1_negative`）

## 啟動範例

### 範例 A（single）

> 「觀察 `HPFineN_HP1S` 特徵，在 7 樣本 paired 模式下跑 Step 1-6。」

預期動作：讀 registry → 確認 HPFineN_HP1S 存在於 G6 extended master → 跑 `scripts/quick_observe.py --feature HPFineN_HP1S --mode paired` → 輸出 `features/HPFineN_HP1S.md`。

### 範例 B（group）

> 「建立 G11 normal_bam_features 群組，包含 `NNormalReads_Unmethylated_Ratio` 與 `Normal_ASM_Delta` 兩個特徵。」

預期動作：Step 0 確認兩特徵 data source → 生成 group-level observation script → Step 1-6 per-feature → synthesize G11.md。

### 範例 C（orthogonal）

> 「同時觀察 `NumReads`, `HPMergedDelta`, `SampleASM_Delta` 的 3-way interaction。」

預期動作：spawn 3 parallel agents (Step 1-4 per feature) + 1 agent for pairwise scatter / partial AUC / 3-way residual。

## 快速檢核（AI 自我審查用）

輸出 feature.md 前先問自己：

```
☐ 10 章節全寫了嗎（不是 3 章節）
☐ 6 張圖都生成了嗎（Step 2 是 4 張）
☐ fig01 violin 有 KDE + n 標註嗎
☐ Step 5 residualize 用 within-group OLS（不是 pooled）嗎
☐ Step 5 AF-bin 極差 <0.10 嗎
☐ Step 6 spatial artifact 檢查過嗎
☐ Cross-sample 固定順序對嗎
☐ 第 9 章至少 3 個 citation 且 /citation-verification 驗過嗎
☐ 第 10 章 3 個 grill-me 質疑寫了嗎
☐ Verdict 有根據 rubric（非主觀）嗎
☐ evidence_ledger.jsonl 登記了嗎
```

## 回報範本

> 本觀察（feature=X）通過 feature-layered-observation Step 1-6：raw AUC={A}；resid AUC={B}；cross-sample {n}/7 同向；verdict={V}；tier={T}；圖存於 `figures/G{N}/`；報告 `features/{feature_id}.md`。

## 相關資源

- **方法學 SOP**：`research/feature_layered_observation/02_methodology.md`
- **Feature registry**：`research/feature_layered_observation/data/feature_registry.tsv`
- **Phase E 彙整**：`research/feature_layered_observation/00_main_observation.md`
- **範例（G1 完整）**：`research/feature_layered_observation/features/G1_coverage.md` + `scripts/observe_G1_coverage.py`
- **References 子目錄**：
  - `references/feature_registry_template.tsv` — 欄位定義
  - `references/observation_report_template.md` — 10 章節骨架
  - `references/step_1_6_sop.md` — Step 1-6 詳細 SOP
  - `references/verdict_rubric.md` — verdict 決策樹
- **Quick script**：`scripts/quick_observe.py` — 通用單特徵觀察器
- **相關 skills**：`/auc-confound-guard`, `/multi-sample-consistency`, `/known-pitfalls`, `/citation-verification`, `/observation-analysis`, `/grill-me`
- **Memory**：
  - `feedback_L2_collider_bias`
  - `feedback_pooled_ols_residualization_trap`
  - `feedback_spatial_autocorrelation_confound`
  - `feedback_feature_name_vs_definition_rule`
  - `project_beyond_auc_exhaustion_confirmed`

---

## Phase & Chain Position

- **Phase**: **P3 PILOT**（Resilient Waterfall harness 第 4 phase）— 主分析 skill
- **Chain**: forward-link chain #3 第 3 環
  ```
  /check-staleness（P2 PRECHECK）→ /test-quick (P3 啟動測試)
      ↓
  feature-layered-observation ← (本 skill: Step 1-6 P3 PILOT 主分析)
      ↓
  auc-confound-guard（若 AUC>0.58 則 confound 三關）
      ↓
  results-analysis（統計與圖表收尾）
      ↓
  multi-sample-consistency（P4 GENERALIZE 跨 7 樣本擴展）
      ↓
  /run-evaluator（P5 EVALUATE retraction risk）
  ```
- **上游觸發**: `/check-staleness` PASS 後 + `research-loop` 完成 P1 plan.json + 用戶手動「觀察特徵」
- **下游 skill**: `auc-confound-guard`（必走，當 AUC > 0.58）→ `multi-sample-consistency`（必走，cross-sample 擴展）→ `/run-evaluator`

## Dependencies

| 類別 | 項目 |
|---|---|
| **Uses** (本 skill 內部呼叫) | Bash（python3 inline scripts: AUC / 32-cell / Wilcoxon / OLS residualization / chr-shuffle null）、Read（讀 merged_with_vcf.tsv.gz / feature_registry.tsv / 各 G{N} master）、Write（寫 10 章節 .md + 6 類圖 + pilot.json）、`/known-pitfalls` skill（Step 0 觸發）、`/auc-confound-guard` skill（Step 5 內呼叫） |
| **Used by** (誰會觸發本 skill) | `research-loop` Step 4 forward-link / `/test-quick` 後接續 / 用戶手動「觀察特徵」「P3 PILOT」 |
| **Reads** | `research/feature_layered_observation/data/merged_with_vcf.tsv.gz`、`research/feature_layered_observation/data/feature_registry.tsv`、`research/feature_layered_observation/data/G{N}/master_g{N}.tsv.gz`、`docs/reports/research_landscape/03_ISM分析價值界定.md`、`state/cycles/{id}/plan.json`（讀 preconditions） |
| **Writes** | **`state/cycles/{id}/pilot.json`** (核心輸出，schema 對齊 `state/schemas/pilot.schema.json`)、`research/feature_layered_observation/G{N}/00_main_observation.md` (10 章節報告)、`research/feature_layered_observation/G{N}/figures/*.png` (6 類圖) |

## Failure Mode & Diagnostics

| # | 失敗症狀 | 先看哪 | 排查步驟 |
|---|---|---|---|
| 1 | 特徵不存在於 merged 或 master | `feature_registry.tsv` 是否有該 entry；`merged_with_vcf.tsv.gz` columns 用 `head -1` 確認 | 檢查 source pipeline；若是新特徵，先把 C++ 計算結果重 join 到 master_g{N}.tsv.gz |
| 2 | `prior_conclusion` 已標 NEGATIVE / CONFOUND_COLLAPSED | `feature_registry.tsv` `prior_conclusion` 欄位 | 一行告知後繼續（用戶明示重驗）或跳過；若重驗，必須在新 cycle plan.json 標 `revival_reason` |
| 3 | 特徵已標 `annotation_only` | `docs/reports/research_landscape/03_ISM分析價值界定.md` | 降為 characterization-only observation；不要試圖找 F1 改善（個人風格 anchor #5 「One-turn freeze」避免徒勞 iteration） |
| 4 | Step 5 confound guard CONFOUND_COLLAPSED（L2/L3 AUC 差距 >0.10） | `pilot.confound_guard.af_bin_check` / `within_group_ols` 兩欄；`/known-pitfalls` P-01/P-02 | 標記 P-01 (L2 collider) 或 P-02 (pooled OLS)；evaluator 會在 P5 自動 catch；若用戶仍要推進，須 override + 完整證據鏈 |
| 5 | Step 6 spatial autocorr fail（chr+pos aggregation AUC 假象） | `pilot.confound_guard.spatial_autocorr` + `feedback_spatial_autocorrelation_confound` memory | 必跑 mid-TP-rate window 驗證；通不過則 verdict=ARTIFACT |
| 6 | Step 3 跨 7 樣本一致性 <5/7 | `pilot.metric_results` × 7 sample 表 | 個人風格 anchor #1「L4 多層必驗」；<5/7 時不該升級 tier，仍可 characterization-only |
| 7 | data 含 caller AF 與 master AF 混用 (merged AF trap) | `vcf_AF` vs `AF` 欄位差異；`feedback_merged_dataset_af_and_loh_pitfalls` | 強制用 `vcf_AF` 作 caller VAF 切分；`AF_class` 已預先用 `vcf_AF` 切好可直用 |

**何時升級到別的 skill / agent / 人工審查**：
- AUC > 0.58 → 強制呼叫 `/auc-confound-guard`（Step 5 已內建，但外顯為獨立 gate）
- Step 3 一致性 5-6/7（邊緣）→ 升級 `multi-sample-consistency` 跑完整 7-sample fan-out
- 結論將升 ⭐4-5 → P5 必呼叫 `/run-evaluator`（pre_tier_upgrade_check.sh hook 強制）
- 涉及 BAM 計算改動（如 HP tag 重新生成） → `methodology-audit` + `cpp-change`

**個人風格適配**（依 `feedback_*` memory）：
- **Anchor #1 「L4 多層驗證必建」** → Step 5 confound guard 是 L4 必經；單跑 Step 1 global AUC 不夠
- **Anchor #5 「One-turn mechanism freeze」** → Step 0 prior_conclusion 已標 NEGATIVE 的特徵，**不在迴圈中試「換 filter 重做」**；要重驗必須在新 plan.json 寫 `revival_reason` 並過 inject-hypothesis
- **Anchor #6 「5 層圖表細節」** → 6 類圖 caption 須含 4 元素：發現 / 樣本 n + 統計量 / 對應表格 / 下一步預期；avoid 紅綠單色（色盲無關）
- **Anchor #7 「Pivot 容忍」** → Step 3 一致性 0/7 時，建議 `pivot-direction` 而非 `keep retrying`
