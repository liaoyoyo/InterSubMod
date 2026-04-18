<!--
建立時間: 2026-04-08 02:00
目標: InterSubMod 階段性研究報告 — 系統性否定 + 正面成果 + 可執行後續
處理範圍: 2025-11 ~ 2026-04-08 全部研究方向
關聯檔案:
  - docs/reports/research_landscape/00_INDEX.md
  - docs/experiments/INDEX.md
  - docs/CURRENT_FOCUS.md
-->

# InterSubMod 階段性研究報告 (v1)

**日期**: 2026-04-08
**版本**: v1.1 (review-revised)
**涵蓋期間**: 2025-11 ~ 2026-04-08

---

## Executive Summary

InterSubMod 是一個以長讀取測序甲基化模式分析腫瘤亞克隆結構的 C++17 工具。本報告整合了約 5 個月的系統性研究，涵蓋：

- **122 個分析觀察**（O1-O13, G1-G7, R1-R5, Wave 1-3, O9, TO-pure 等）
- **748,391 TP/FP regions + 122,790 FN regions**（7 cancer cell lines × 2 modes）
- **60+ ISM 特徵**系統評估

### 核心結論

1. **ISM 甲基化特徵無法作為 variant filter**（全方向 NEGATIVE，AUC 天花板 0.58）
2. **ISM 有明確的生物學價值**（ASM 定量 32-66%、somatic heterogeneity 標記）
3. **Self-phasing 循環依賴是 TO 模式 HP 特徵失效的根因**（已因果驗證）
4. **PON-only phasing 是唯一已知的修正路徑**（N50 +99.7%，待全量驗證）
5. **Caller AF 是 TO 模式唯一有效判別器**（AUC=0.654，超越全部 ISM）

### 研究定位轉向

```
原始目標: ISM = somatic variant filter → 提升 F1
    ↓ 系統否定後
現在定位: ISM = read-level epigenetic characterizer
          → ASM 定量、subclone 結構解析、二次打擊推論
```

---

## 1. 方法論概覽

### 1.1 系統架構

```
ONT Long-Read BAM (tumor ± normal)
    ↓
LongPhase (phasing + haplotagging)
    ↓
InterSubMod C++17
    ├── ReadParser: MM/ML tag 解碼 → CpG 甲基化狀態
    ├── DistanceCalculator: 6 種距離度量（Bernoulli 為主）
    ├── TreeCutter: 階層式 clustering
    ├── SignificanceAnalyzer: PERMANOVA + χ² + LabelTest
    └── QualityScorer: 綜合評分
    ↓
60+ ISM 特徵 per region
```

### 1.2 數據規模

| 資料集 | Rows | 描述 |
|--------|------|------|
| Master TP/FP | 748,391 | 7 samples × 2 modes (paired/TO) × TP+FP |
| FN manifest | 122,790 | 7 samples × 2 modes × FN (2026-04-08 新增) |
| 特徵數 | 116 cols | ISM 60+ + Caller + metadata |
| Samples | 7 | HCC1395, HCC1395_DORADO, COLO829, H1437, H2009, HCC1937, HCC1954 |

### 1.3 評估框架

- **LOSO (Leave-One-Sample-Out)**: 7-fold cross-validation，防止 sample-specific overfitting
- **Per-sample AUC**: 每個特徵 × 每個 sample 單獨計算，確認跨樣本一致性
- **三層分析**: L1 pooled → L2 per-sample → L3 AF-bin stratified

---

## 2. 系統性否定結果 — 完整清單

### 2.1 Variant Filter 方向（全部 NEGATIVE）

| # | 方向 | 核心數字 | 穩定度 | 日期 |
|---|------|---------|--------|------|
| 1 | O1-O10 系統觀察 | TO max AUC=0.566 (caller_gq)；QS=0.497 | ⭐4 | 03-31 |
| 2 | O11 Within-group heterogeneity | raw AUC=0.845 → corrected 0.530 (n_reads confound) | ⭐5 | 04-01 |
| 3 | O12 LOH 場景區分 | 3 場景不可分；AlleleDelta=AF confound；L2 collider bias | ⭐4 | 04-02 |
| 4 | O13 跨區域 Correlation | shared read count confound；OLS p=0.464, d=-0.071 | ⭐5 | 04-03 |
| 5 | G1-G7 Germline FP 鑑別 | 60+ 特徵全 AUC<0.64；TP loss≤2% → FP removal=0% | ⭐4 | 04-04 |
| 6 | Read-level Germline FP | LOSO AUC=0.721 但安全約束 FP removal=0% | ⭐3 | 04-04 |
| 7 | Wave 3 LOH 分層 J11-J16 | combo AUC=0.577；cnLOH 0.587=Simpson's Paradox | ⭐4 | 04-06 |
| 8 | R1-R5 特徵設計 | identifiability problem（非設計問題） | ⭐4 | 04-07 |
| 9 | Option C HP-free 雙路 | ClusterPermanovaF=0.512；HP-free combo=0.564 | ⭐4 | 04-07 |
| 10 | **O9 FN Rescue** | **HP-free AUC<0.53；TO QS=0.338 反轉** | ⭐4 | 04-08 |
| 11 | **TO-pure 獨立建模** | **ISM 增量 +0.003~+0.030 over caller-only** | ⭐4 | 04-08 |

### 2.2 否定結論的邏輯結構

```
根本原因: Germline FP ≈ Somatic TP in methylation space
    │
    ├── HP-free features → AUC 0.50-0.53 (random)
    │   ├── NumReads, NumCpGs: coverage, not biology
    │   ├── CramersV: 93% zero (2×2 framework limitation)
    │   ├── PairwiseMeanDist/MedianDist: methylation distance, no TP/FP signal
    │   └── ClusterPermanovaF: 0.512 (random)
    │
    ├── HP-dependent features → AUC 0.55-0.72 但有 confounds
    │   ├── HP2FamilyN: 0.72 → 0.54 after circular artifact removal
    │   ├── HPFineNGroups: 0.617 residualized (subclone marker, not filter)
    │   ├── AlleleDelta: AF proxy, not methylation
    │   └── All TO HP features: self-phasing contaminated
    │
    └── Caller features → best available
        ├── caller_af: 0.654 (single best)
        ├── caller_ad_alt: 0.655
        └── ISM adds only +0.003-0.030 on top
```

### 2.3 為什麼每個 confound 被排除

| Confound | 發現方式 | 排除方法 | 結果 | 來源 |
|----------|---------|---------|------|------|
| n_reads | Spearman r=0.79 with epipolymorphism | Residualize on n_reads | AUC 0.845→0.530 | `.../20260401_O11_heterogeneity/` |
| Shared read count | FP median=36 vs TP=21 reads shared | 分層+OLS+matching | p=0.464, d=-0.071 | `.../20260403_O13_cross_region/` |
| AF (allele frequency) | AlleleDelta=AF proxy | AF-bin stratification | All within-bin AUC<0.55 | `.../20260402_O12_loh_scenarios/` |
| Self-phasing circularity | 94.6% somatic→HP1, bias 17.3:1 | PON-only phasing | N50 +99.7%, bias eliminated | `.../20260402_self_phasing_causal_analysis/` |
| Simpson's Paradox | cnLOH pooled 0.587 | Per-sample decomposition | Mean=0.50 | `.../20260406_wave3_loh_stratification/` |
| L2 Collider bias | Near-constant feature residualized on AF | L3 AF-bin cross-validation | Spurious signal disappears | `.../20260402_O12_loh_scenarios/` |

---

## 3. 正面成果

### 3.1 ASM (Allele-Specific Methylation) 定量 — ⭐4 POSITIVE

**發現**: 32-66% 的體細胞突變位點顯示等位基因特異性甲基化。

| 方法 | ASM 比例 | 驗證方式 |
|------|---------|---------|
| AlleleDelta > 0.1 | 66% | 直接閾值 |
| AlleleDelta > 0.2 | 45% | 嚴格閾值 |
| PERMANOVA p<0.05 | 32% | 統計檢定 |
| Haplotype methylation KDE | 一致 | 視覺化確認 |
| Within-haplotype entropy | 一致 | 獨立度量 |

**意義**: ISM 能定量測量每個 somatic SNV 位點的等位基因甲基化差異，這是其他工具難以提供的 read-level 表觀遺傳學資訊。

**限制**: FP 也顯示 ASM（germline variants 也有 ASM），因此不能用於 TP/FP 區分。

**主要度量**: 報告中引用的 ASM 比例以 **AlleleDelta > 0.1**（66%）為寬鬆估計、**PERMANOVA p<0.05**（32%）為嚴格估計。來源：`.../20260405_snv_methylation_association/`。

### 3.2 HPFineNGroups Subclone Marker — ⭐4 POSITIVE (Paired mode)

**發現**: HPFineNGroups ≥ 4（多甲基化亞群）是 somatic heterogeneity 的可靠標記。

| 指標 | 數值 | 來源 |
|------|------|------|
| N≥4 + NR≥80 → TP rate | **89.1%** | `.../20260407_hpfinengroups_subclone_marker/` |
| 低 AF (0.1-0.2) 信號增量 | **+50 percentage points** | AF-bin stratified analysis |
| 跨樣本一致性 | **7/7 方向一致** | per-sample AUC table |
| Residualized AUC | **0.617** | AF-residualized, paired mode |
| Cohen's d | 0.35 | pooled paired data |

**意義**: 當 ISM 偵測到 4+ 個甲基化亞群且有足夠 reads 時，該位點幾乎確定是真實體細胞突變，而非 germline leak。這不是 filter（不能用它排除 FP），但是有力的 somatic heterogeneity 解釋工具。

**重要限制**: HPFineNGroups 是 HP-dependent 特徵（依賴 haplotype assignment 進行分群）。上述數字來自 **paired mode**（使用 normal BAM 的可靠 germline phasing）。在 TO mode 下，self-phasing 汙染可能膨脹分群數量。因此 TO mode 的 HPFineNGroups 結果**待 PON-only phasing 重跑後驗證**——若重跑後趨勢維持，則確認其跨模式穩健性；若消失，則限縮為 paired-only 工具。

### 3.3 Self-Phasing 因果鏈 — ⭐4 CONFIRMED

**完整五步因果鏈** (2026-04-02/03 驗證):

```
Step 1: LongPhase-TO 用 somatic VCF 作 phasing input
Step 2: 94.6% somatic variants 被 assign 到 HP1 (bias 17.3:1)
Step 3: HP_Ratio 在 TO 模式失真 (paired-TO r=0.001)
Step 4: 62% TO LOH 是 artifact (d=-1.20)
Step 5: 31.2% self-phasing LOH detected
```

**量化影響**: 23 TSV + 7 PNG 完整證據鏈。來源：`.../20260402_self_phasing_causal_analysis/`（94.6% 偏差數字來自 `step2_hp_distribution.tsv`）。

### 3.4 PON-Only Phasing — ⭐4 PARTIAL SUCCESS

| 指標 | Baseline TO | PON-only | 改善 |
|------|------------|----------|------|
| Somatic bias | 17.3:1 | **消除** | 完全修正 |
| Phase block N50 | — | **+99.7%** | 大幅提升 |
| Phased rate | — | **+23.6pp** | 顯著 |
| LOH.bed | Jaccard=1.0 | 不變 | LOH 定義穩健 |

**待辦**: 需 haplotag + ISM 全量重跑以驗證 HP-dependent 特徵是否突破。

---

## 4. Quality Score 問題

### 4.1 TO QS 完全失效

| 指標 | Paired | TO |
|------|--------|-----|
| QS AUC (TP vs FP) | 0.754 | **0.497 (random)** |
| QS AUC (FN vs TP) | 0.551 | **0.338 (inverted!)** |
| LOH penalty 效果 | 正確 | **反向** |

### 4.2 根因分析

TO QS 失效的根因鏈：
1. LOH penalty 基於 HP_Ratio 判斷 LOH
2. TO mode 的 HP_Ratio 被 self-phasing 汙染
3. 大量 TP 被錯誤標記為 LOH → QS 被壓低
4. FP 反而因為不在 LOH 區域 → QS 較高

### 4.3 已實施修正

`QualityScorer.cpp` 已實作 mode-aware 邏輯：TO mode 停用 LOH penalty 和 verify bonus。但根本修正需要 PON-only phasing 後的正確 HP tags。

---

## 5. TO 模式的完整診斷

### 5.1 問題全貌

```
ClairS-TO 輸出
├── TP: ~310K (7 samples)
├── FP: ~132K (FP rate 8.7% ~ 74.6%, 8.6× spread)
│   └── 98.6% = raw_absent (paired 從未 call)
│       └── 多數是 germline leak
└── FN: ~78K (TO recall 0.59 ~ 0.92)

PON 過濾: 移除 99.48% raw FP → 殘餘 FP = ISM 的目標 [來源: TO FP provenance analysis]
ISM 對殘餘 FP: AUC < 0.58 → 無法過濾
```

### 5.2 TO-pure LOSO 建模結果 (2026-04-08)

| Model | Features | LR AUC | RF AUC |
|-------|----------|--------|--------|
| A: HP-free only | 8 methylation | 0.529 | 0.535 |
| B: All ISM | 60+ features | 0.601 | 0.635 |
| C: Caller only | AF, DP, GQ, AD | **0.632** | — |
| D: ISM+Caller | All combined | 0.636 | **0.662** |

**ISM 增量**: D - C = +0.004 (LR) ~ +0.030 (RF) — 可忽略。

> **注**: Model B (All ISM, RF) AUC=0.635 略高於 Model C (Caller only, LR) AUC=0.632，看似 ISM 超越 caller。但兩者使用不同學習器（RF vs LR），不可直接比較。正確比較應控制方法：D vs C (LR) = +0.004，D vs C (RF, D=0.662 vs C 未測 RF) 也需同方法。Model B 的 RF 0.635 反映的是 RF 對 ISM 特徵的非線性擬合能力（含 HP-dependent 殘餘信號），非 ISM 本身的判別力。來源：`.../20260408_to_pure_independent_modeling/loso_results.tsv`。

**Per-fold LOSO 變異**: 各 fold AUC 標準差 0.03-0.08（H2009 為 outlier fold，AUC 偏高因 FP rate=8.7% → 極不平衡）。完整 7-fold 結果見 `.../20260408_to_pure_independent_modeling/loso_results.tsv`。

### 5.3 Per-Sample FP Rate Heterogeneity

| Sample | FP Rate | 特徵 |
|--------|---------|------|
| HCC1954 | **74.6%** | 極高 FP |
| HCC1937 | 48.8% | 高 FP |
| COLO829 | 34.6% | 中等 |
| HCC1395 | 28.9% | 中等 |
| HCC1395_DORADO | 28.6% | 中等 |
| H1437 | 22.8% | 較低 |
| H2009 | **8.7%** | 極低 FP |

8.6× spread — 這意味著任何 LOSO 模型需在極不同的 class balance 下泛化。

---

## 6. FN (False Negative) 特徵觀察 (O9, 2026-04-08)

### 6.1 數據

- 122,790 FN regions（44,415 paired + 78,375 TO）
- 7 samples 完整 ISM 執行

### 6.2 FN vs TP Feature AUC

| Feature | Paired AUC | TO AUC | Category |
|---------|-----------|--------|----------|
| LabelAllelePermanovaF | 0.664 | 0.636 | Allele (AF proxy) |
| AlleleDelta | 0.642 | 0.607 | Allele (AF proxy) |
| HPFineF | 0.578 | 0.558 | HP-dependent |
| Quality_Score | 0.551 | **0.338** | HP-dependent |
| NumReads | 0.507 | 0.503 | **HP-free (random)** |
| CramersV | 0.507 | 0.503 | **HP-free (random)** |
| PairwiseMeanDist | 0.492 | 0.489 | **HP-free (random)** |

### 6.3 結論

**NO-GO**: 甲基化空間中 FN ≡ TP（HP-free 特徵 AUC 全 <0.53，Cohen's d < 0.05）。

**FN 被遺漏的推論原因**: ISM 資料中唯一能區分 FN vs TP 的信號來自 AF-proxy 特徵（LabelAllelePermanovaF=0.664），而非甲基化特徵。這暗示 FN 被 caller 遺漏的原因是讀取證據不足（低 AF / 低 coverage），而非表觀遺傳學差異。此推論基於間接證據（ISM 特徵的 AUC 分布），非直接 AF 比較——直接驗證需比對 truth set AF 與 caller AF threshold，屬 caller 層面分析，超出 ISM 範疇。

來源：`.../20260408_O09_fn_characterization/O9_FN_characterization_report.md`，`fn_vs_tp_summary_auc.tsv`。

---

## 7. 結論穩定性總覽

### 7.1 穩定度分布

| 等級 | 數量 | 結論 |
|------|------|------|
| ⭐5 堅固 | 3 | Paired/TO 分離、O11 n_reads、O13 shared reads |
| ⭐4 穩固 | 8+ | TO AUC<0.58、O12、G1-G7、Self-phasing、ASM、O9、TO-pure |
| ⭐3 需注意 | 5 | PON 99.48%、Read-level、Phase 1A、LOH Jaccard、QS |

### 7.2 不會被 PON-only 重跑翻轉的結論

以下結論**不依賴 HP tags**，即使 self-phasing 修正後也不會改變：

1. **O11 (n_reads confound)** — 與 phasing 無關
2. **O13 (shared read confound)** — 與 phasing 無關
3. **G1-G7 (VCF 特徵 NO-GO)** — VCF 特徵不受甲基化影響
4. **FP Provenance (98.6% raw_absent)** — Caller 行為，不受 ISM 影響
5. **O9 FN Rescue (HP-free random)** — HP-free 特徵不依賴 phasing
6. **TO-pure ISM 增量微弱** — caller_af 優勢是結構性的

### 7.3 可能被 PON-only 重跑改善的方向

1. **HP-dependent 特徵 AUC** — paired HP1FamilyN=0.834 暗示有真實信號
2. **within_dom_alt_frac** — read-level feature，HP 修正後可能改善
3. **Phase 1A TO 負增益** — 目前 TO ISM 是負增益，修正後可能轉正
4. **QS 公式** — 正確 HP 後 LOH 判斷改善 → QS 可能恢復區分力

---

## 8. 現在可執行的工作（不依賴 PON-only 重跑）

### 8.1 C++ 程式碼改進

| 優先 | 項目 | 說明 | 狀態 |
|------|------|------|------|
| P0 | ReadParser HP tag 過濾 | 忽略 somatic HP tags (11/21/33)，只用 germline 1/2 | 待做 |
| P0 | QS TO mode 公式重設計 | 基於 caller_af 而非 HP-based LOH | mode-aware 已做，根本修正待做 |
| P1 | CramersV 2×2 框架修正 | 93% 為零，改為多群框架 | 設計完成，待實作 |
| P1 | Output 格式擴充 | 增加 caller feature passthrough 欄位 | 待做 |
| P2 | Gating 邏輯 TO 適配 | TO mode 22% passed_gate，考慮放寬條件 | 待評估 |

### 8.2 Python 分析腳本

| 優先 | 項目 | 說明 |
|------|------|------|
| P1 | PON-only haplotag pipeline | 自動化 LongPhase haplotag + ISM 重跑腳本 |
| P1 | 重跑後 before/after 比較框架 | 自動比較 HP-dependent 特徵 AUC 變化 |
| P2 | H2009 異質性診斷腳本 | Phase 1A 唯一負向 sample 的深度分析 |
| P2 | Platform normalization | 5kHz / DORADO / PAO 跨平台特徵對齊 |

### 8.3 研究方向

| 優先 | 方向 | 依賴 | 預期收益 |
|------|------|------|---------|
| P0 | PON-only haplotag + ISM 全量重跑 | PON-only BAM（已有） | 解鎖 HP-dependent 特徵 |
| P1 | 重跑後 HP-dependent 重評估 | 重跑完成 | 可能突破 AUC 0.58 |
| P2 | 低純度臨床樣本驗證 | 臨床數據 | 驗證 read-level FP 潛力 |
| P2 | 論文撰寫 | 以上全部 | 發表 |

---

## 9. 關鍵數字速查表

### ISM 特徵 AUC（TO mode, 最佳值）

| 特徵 | TP vs FP | FN vs TP | 用途 |
|------|---------|---------|------|
| caller_af | **0.654** | — | 最佳判別器 |
| caller_ad_alt | **0.655** | — | ≈ caller_af |
| LabelAllelePermanovaF | 0.558 | 0.636 | Allele proxy |
| HPFineNGroups | 0.571 | 0.521 | Subclone marker |
| HP2FamilyN | 0.523 (circular removed) | — | Self-phasing artifact |
| PairwiseMeanDist | 0.539 | 0.489 | 弱/random |
| CramersV | 0.509 | 0.503 | Random |
| ClusterPermanovaF | 0.522 | 0.514 | Random |
| Quality_Score | 0.518 | **0.338** (inverted) | 失效 |

### 模型 AUC 比較

| 設定 | AUC |
|------|-----|
| TO ISM+Caller (RF, LOSO) | 0.662 |
| TO Caller-only (LR, LOSO) | 0.632 |
| TO All ISM (RF, LOSO) | 0.635 |
| TO HP-free only (RF, LOSO) | 0.535 |
| Paired All ISM (RF, LOSO) | 0.574 |
| Paired HP-free (LR, LOSO) | 0.553 |

### 關鍵常數

| 指標 | 值 |
|------|-----|
| ASM 比例 | 32-66% |
| HPFineNGroups N≥4 TP rate | 89.1% |
| Self-phasing HP bias | 17.3:1 (HP1) |
| TO LOH artifact 比例 | 62% |
| PON FP removal rate | 99.48% |
| Paired QS AUC | 0.754 |
| TO QS AUC | 0.497 |

---

## 10. 報告、數據與圖表位置

### 觀察工作區

| 工作區 | 路徑 |
|--------|------|
| O1-O10 系統觀察 | `.../20260318_per_sample_tp_fp_analysis/` |
| O9 FN 特徵 | `.../20260408_O09_fn_characterization/` |
| TO-pure 建模 | `.../20260408_to_pure_independent_modeling/` |
| LOH Round 1-4 | `.../20260327_loh_round1_cross_sample_audit/` ~ `20260328_loh_round4_final_validation/` |
| Self-phasing 因果鏈 | `.../20260402_self_phasing_causal_analysis/` |
| Germline FP G1-G7 | documented in research_landscape reports |

### 研究全景文件

| 文件 | 說明 |
|------|------|
| `docs/reports/research_landscape/00_INDEX.md` | 全局入口 |
| `docs/reports/research_landscape/01_TO_FP問題全貌.md` | FP 來源分析 |
| `docs/reports/research_landscape/02_Self_Phasing根因.md` | Self-phasing 因果鏈 |
| `docs/reports/research_landscape/03_ISM分析價值界定.md` | 特徵 HP 分類 + 價值重定位 |
| `docs/reports/research_landscape/04_暫停判定與重評估.md` | 哪些需重測 |
| `docs/reports/research_landscape/05_證據鏈總覽.md` | 8 條證據鏈 |
| `docs/reports/research_landscape/06_結論穩定性審查.md` | 14 結論穩定度 |

---

## 附錄 A: 研究時間線

| 日期 | 里程碑 |
|------|--------|
| 2025-11 | InterSubMod 基礎設施建設（MM/ML 解碼、距離計算、PERMANOVA） |
| 2025-12 | 初步 TP/FP 特徵富集分析、F1 最佳化 |
| 2026-01~02 | 雙向驗證、Subsample 分析、Phase 1 baseline |
| 2026-03-11~18 | TO pilot 7 samples、Master dataset 建構 |
| 2026-03-20~27 | Phase 1A ML、LOH 深度調查 Round 1-4 |
| 2026-03-30 | HP integer tag bug fix + 全量重跑 |
| 2026-03-31 | O1-O10 系統觀察（82 圖表） |
| 2026-04-01 | O11 Heterogeneity (NEGATIVE) |
| 2026-04-02 | O12 LOH 場景 (NEGATIVE)、Self-phasing 因果鏈 (CONFIRMED) |
| 2026-04-03 | O13 跨區域 (NEGATIVE)、PON-only phasing (PARTIAL SUCCESS) |
| 2026-04-04 | G1-G7 Germline FP (NO-GO)、Research Landscape 6 文件 |
| 2026-04-05 | ASM 定量 (POSITIVE)、Read-level (CONDITIONAL NO-GO) |
| 2026-04-06 | Wave 3 LOH 分層 (CLOSURE)、LOH concordance |
| 2026-04-07 | R1-R5 特徵設計、Option C HP-free (NEGATIVE)、HPFineNGroups (POSITIVE) |
| **2026-04-08** | **O9 FN (NO-GO)、TO-pure 建模 (NEGATIVE)、本報告** |
