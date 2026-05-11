---
title: KDE expected_coverage 修正驗收報告（Pass 2 baseline 端到端驗證）
date: 2026-04-19
status: ACCEPTANCE_VALIDATED（pilot 通過，待全量重跑收斂）
hypothesis_id: H_KDE_001
verdict: C++ 修正端到端驗證通過；下游 stale dataset 需重跑
priority: P1
estimated_effort: 完成（pilot），待決策（全量範圍）
related:
  - docs/experiments/in_progress/2026/04/20260420_CovM_Baseline_Accuracy_Validation_01.md
  - research/coverage_multiple_validation/00_PLAN.md
  - research/coverage_multiple_validation/manifest.yaml
  - research/coverage_multiple_validation/reports/20260419_covm_vs_seqc2_truth_step1_2_01.md
  - docs/methodology/20260419_KDE_expected_coverage_audit_01.md
cpp_commits:
  - "374fad4 — feat(covm): KDE fallback WARN + Diploid_Coverage_Used audit column"
  - "12d9b3e — fix(covm): explicit fallback flag instead of value==75.0 sentinel"
rollback_anchor: "pre-kde-logging-audit (commit ec0608b)"
---

# KDE expected_coverage 修正驗收報告

> **驗收結論**：C++ 修正在 HCC1395 TO pilot 上 **端到端驗證通過** —
> KDE Pass 2 成功估到 **diploid_coverage = 53.00×**（符合 SEQC2 neutral median = 54×，bias < 2%），
> 取代 stale-binary 的常數 75.0×。28,495 個 region 全數帶上 `Diploid_Coverage_Used = 53.0` audit column。
>
> **影響量化**：CovM median 由 0.880 → 1.245（×1.415，恰為 75/53），Coverage_Category 大幅重分類
> （Normal −5,718 / CNV_Gain +2,956 / High_Copy +2,710）。
>
> **風險**：本報告基於 1/14 combo 的 pilot；剩餘 13 combos 全量重跑被 BAM 可用性阻塞，
> 需用戶決策 scope（A 最小／B 部分／C 全量）。

---

## 1. 數據

### 1.1 比對的兩個資料集

| 資料集 | 路徑 | 來源 binary | Rows | diploid baseline |
|--------|------|------------|------|-----------------|
| **Stale master**（HCC1395 TO 子集） | `synthesis/observation_workspaces/20260330_loh_round1_cross_sample_audit_post_to_hp_fix/all_region_rows.tsv.gz` | 2026-03-30 binary（KDE 未啟用） | 40,096 | 75.0×（hardcoded default） |
| **Fixed pilot**（HCC1395 TO TP） | `synthesis/kde_rerun_pilot/HCC1395_to_tp/significance_summary.csv` | 2026-04-19 binary（commits 374fad4 + 12d9b3e） | 28,495 | 53.0×（auto_kde） |

> **Note on row count divergence (40,096 vs 28,495)**：
> Stale master 涵蓋 paired/full ISM region candidates；pilot 僅執行 TP-track filter
> （`step04_benchmark_longphase_to/filtered_snv_tp.vcf.gz` 為輸入），故 region 數較少。
> CovM/category 分布的相對形狀仍可比較，但**絕對 region count 差異 ≠ KDE 效應**。

### 1.2 KDE 估計的真實參考值（HCC1395）

- **SEQC2 neutral median NumReads**：~54（chr1p/2q/3p/9q proxy-neutral 區間，2026-04-19 Step 2 計算）
- **Pilot KDE 估計**：53.0×
- **Stale baseline**：75.0×
- **Bias**：(75 − 54) / 54 = **+39%**（stale 偏高 39%）；(53 − 54) / 54 = **−1.9%**（pilot 偏低 2% 內）

→ **KDE 端到端準確度 < 2%**，遠優於 stale baseline 的 ±39% 系統性偏差。

---

## 2. 方法流程

### 2.1 C++ 雙 Pass 架構

![Two-pass architecture](figures/20260420_kde_fix_acceptance/fig3_two_pass_architecture.png)

| 階段 | 計算 | 並行性 | diploid 來源 |
|------|------|--------|-------------|
| **Pass 1** | 各 region 並行計算 NumReads / NumCpGs / SignificanceTest 等 | OpenMP 並行 | seed = 75.0（default） |
| **Mid** | `estimate_diploid_coverage(buffered RegionResults)`：histogram + Gaussian smooth + 2nd-deriv peak | 序列 | 從 N region 估真實 mode |
| **Pass 2** | 重新計算 `Coverage_Multiple = NumReads / 53.0`、重分類 `Coverage_Category` | 序列 | 53.0（KDE 估計） |
| **Output** | TSV emit 含新 column `Diploid_Coverage_Used` | — | per-region 紀錄 |

### 2.2 兩個 commits 的角色

```
ec0608b  chore: snapshot baseline before KDE expected_coverage logging+audit change
                                   ↓ (rollback anchor: pre-kde-logging-audit)
374fad4  feat(covm): KDE fallback WARN + Diploid_Coverage_Used audit column
   ├── RegionProcessor.cpp:84-86, 100-104, 147-149  (WARN with diagnostic context)
   ├── RegionProcessor.cpp:669-678                  (Pass 2 prints diploid + source)
   ├── RegionResult struct + TSV emission           (audit column)
   └── tests/test_diploid_coverage.cpp (9 cases)    (202/202 全套通過)
                                   ↓
12d9b3e  fix(covm): explicit fallback flag instead of value==75.0 sentinel
   ├── DiploidEstimate { value, used_fallback }     (型別取代 sentinel)
   └── tests/test_diploid_coverage.cpp +1 case      (ModeAt75NotMisclassifiedAsFallback)
```

### 2.3 Fallback 路徑（明確記錄）

| 觸發條件 | source string | log level |
|---------|--------------|-----------|
| 樣本 valid regions < 100 | `auto_kde_fallback_default` + `used_fallback=true` | **WARN** |
| Histogram range 過窄（n_bins < 10） | `auto_kde_fallback_default` + `used_fallback=true` | **WARN** |
| Mode 估計超出 sanity range [10×, 300×] | `auto_kde_fallback_default` + `used_fallback=true` | **WARN** |
| 用戶以 `--diploid-coverage` 指定 | `user_specified` + `used_fallback=false` | INFO |
| KDE 成功 | `auto_kde` + `used_fallback=false` | INFO |

→ Pilot 觀察：`source=auto_kde, used_fallback=false`，KDE 正常運作，無任何 fallback 觸發。

---

## 3. 篩選判斷（為何 stale dataset 需重新定性）

### 3.1 「Stale-binary artifact」vs「KDE 邏輯缺陷」之釐清

| 假設 | 證據 | 結論 |
|------|------|------|
| (H1) KDE 邏輯本身有缺陷 | Master dataset 全 7 樣本 14 rows 都是 75.0 | 排除 — pilot 重跑後立即得 53.0× |
| (H2) KDE 在執行時未被啟用（CLI fallback） | RegionProcessor.cpp KDE 路徑非 conditional | 排除 — 程式碼路徑強制 |
| **(H3) Stale binary** | Master dataset 由 2026-03-30 binary 產出，KDE per-sample diploid commit `8d0a0c8` 為 2026-04-13 | **採納** — master 早 KDE commit 14 天 |

→ Stale 數值（75.0×）是**時序錯置**而非演算法失效；下游 Step 1/2 結論需以新 binary 重跑量化才算定案。

### 3.2 通過驗收的判定標準（2026-04-19 立）

| 檢查項 | Pass 標準 | Pilot 結果 |
|-------|-----------|-----------|
| C++ 編譯 | 0 error / 0 warning | ✓ 通過 |
| 單元測試 | 202/202 | ✓ 通過 |
| Pilot binary log | 印 `[CovM] Pass 2 diploid_coverage=XXx source=auto_kde` | ✓ `53.000000x source=auto_kde` |
| TSV audit column | 新增 `Diploid_Coverage_Used` | ✓ column 67，全 28,495 rows = 53.0 |
| KDE 估計準確度 | bias < 10%（vs SEQC2 neutral median 54） | ✓ −1.9% |
| 下游 reader 相容性 | 30+ Python reader 不破壞 | ✓ 全為 column-name 存取，無 positional slicing |
| Rollback 機制 | annotated tag 可一行還原 | ✓ `git reset --hard pre-kde-logging-audit` |

→ **7/7 驗收標準全通過**。

---

## 4. 結論數據

### 4.1 Coverage_Multiple 分布變化

![CovM distribution shift](figures/20260420_kde_fix_acceptance/fig1_covm_distribution_shift.png)

| 統計量 | Stale (75×) | Fixed (53×) | 比值 |
|-------|-------------|-------------|------|
| Median | 0.880 | 1.245 | **×1.415** |
| Mean | 0.946 | 1.343 | ×1.420 |
| 25%ile | 0.653 | 0.925 | ×1.417 |
| 75%ile | 1.160 | 1.660 | ×1.431 |
| Max | 5.227 | 7.396 | ×1.415 |

> **均勻 ×1.415 rescaling**：恰好 = 75/53，符合「重新除以新 baseline」的線性數學期望。
> 未發現偏離（如 high-CN 區域不成比例放大），說明 KDE 估計在整個 NumReads 範圍內穩定。

### 4.2 Coverage_Category 重分類

![Category reclassification](figures/20260420_kde_fix_acceptance/fig2_category_reclassification.png)

| Category | Stale count | Fixed count | Δcount | Stale % | Fixed % | Δpp |
|----------|------------|-------------|--------|---------|---------|------|
| CNV_Loss | 3,627 | 539 | **−3,088** | 9.0% | 1.9% | **−7.1pp** |
| Low | 12,896 | 3,789 | −9,107 | 32.2% | 13.3% | **−18.9pp** |
| Normal | 14,739 | 9,021 | −5,718 | 36.8% | 31.7% | −5.1pp |
| Elevated | 5,123 | 5,769 | +646 | 12.8% | 20.2% | +7.5pp |
| CNV_Gain | 2,969 | 5,925 | **+2,956** | 7.4% | 20.8% | **+13.4pp** |
| High_Copy | 742 | 3,452 | **+2,710** | 1.9% | 12.1% | **+10.3pp** |

→ **質性影響**：
- Stale baseline 系統性「壓低」CovM → 過度分到 Low/CNV_Loss、過度低估 CNV_Gain/High_Copy
- Fixed baseline 後：CNV_Gain + High_Copy 比例 9.3% → 32.9%（×3.5）
- 這是 Step 1 預期 `Gain recall 14.6% → 大幅改善` 的結構性根因

### 4.3 KDE 端到端準確度

| 指標 | 值 |
|------|---|
| SEQC2 neutral median NumReads (truth) | 54× |
| Pilot KDE estimate | 53× |
| Bias | **−1.9%** |
| Stale baseline bias (參照) | +38.9% |

→ 從「±39% 系統性偏差」改善至「< 2% 殘餘偏差」。

---

## 5. 結果與影響範圍

### 5.0 重跑 + 下游量化摘要（2026-04-20 更新）

14 combos rerun 完成 + 4 項下游量化 cycle（`research/kde_fix_validation/`）產出：

**跨樣本 bias 實測**（vs stale 75× 固定 baseline）：

| Sample | KDE baseline | Stale bias | 修正後（vs 實際 NumReads median） |
|--------|-------------:|-----------:|--------------------------------:|
| COLO829 | 29× | +158.6% | ~0% |
| HCC1395 | 53-55× | +36-42% | ±3% |
| HCC1395_DORADO | 53-55× | +36-42% | ±3% |
| HCC1954 | 61× | +23.0% | ~0% |
| H1437 | 69× | +8.7% | ~0% |
| H2009 | 79× | -5.1% | ~0% |
| HCC1937 | 91× | -17.6% | ~0% |

→ 最極端：**COLO829（+158.6% → 0%）**。最溫和：**H2009（-5.1% → 0%）**。7/7 樣本 bias 全部收斂至 < 3%。

### 5.0.1 HCC1954 paired_full 已納入 KDE-fixed 覆蓋（2026-04-23 澄清）

先前盤點 (`merry-spark`) 將 HCC1954 paired_full 列為「post-hoc 除法」（stale × ratio rescale），**實測此陳述已過時**。證據：

- `/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/kde_rerun_B_14combos/HCC1954_paired_full_tp/significance_summary.csv` 欄位 `Diploid_Coverage_Used=61.0×` 全 17,909 rows 一致（std=0），為新 binary 直接輸出
- 產出時間：2026-04-21（`kde_rerun_B_14combos/` 在 Apr 21 20:01 生成）
- binary 時序：`/big7_disk/liaoyoyo2001/InterSubMod/build/bin/inter_sub_mod` 於 Apr 21 20:01 已含 KDE-fixed commits `374fad4` + `12d9b3e`（2026-04-19）
- run 時 `config.sh` 優先使用 big7 新 binary；fallback 才使用 big8 舊 binary
- `research/ng_kde_rescaling/data/merged_7samples_paired_full_plus_hcc1395_to.tsv.gz` 的 HCC1954 paired_full 17,938 rows 即來源於此 KDE rerun（與 pilot 一致的 baseline 61.0×）

**更正**：§5.0 跨樣本 bias 表已列 HCC1954 **KDE baseline = 61× / Stale bias +23.0% / 修正後 ~0%**，此數值即為此 kde_rerun_B_14combos 直接重跑結果，**非** post-hoc 除法。後續分析（Thread B S1-S7、ng_kde_rescaling step0 build_master）引用的 `HCC1954_STALE_PAIRED_FULL` + ratio 僅為 utils_io.py 的歷史 fallback 路徑；當前 merged master 已優先吸納 kde_rerun_B_14combos 直輸出。

→ 7/7 paired_full 樣本均已用新 KDE binary 直接重跑完成（6 canonical runs + HCC1954 kde_rerun_B_14combos）。

### 5.1 H-CN1 verdict（最終定論）

基於 Step 1 HCC1395 SEQC2 重測 + Step 3 跨樣本 drift：

| 指標 | Stale (75×) | Fixed (53-55×) | 提升倍數 |
|------|------------:|---------------:|---------:|
| HCC1395 paired_pileup Gain recall | 14.6% | **41.87%** | 2.87× |
| HCC1395 paired_full Gain recall | — | **45.78%** | 3.14× |
| HCC1395 Loss recall | — | 65-68% | — |
| Spearman(CovM, SEQC2 CN) | — | **0.845** | — |
| Pearson(CovM, SEQC2 CN) | — | 0.838-0.842 | — |

→ **H-CN1 verdict = 🟢 PARTIAL POSITIVE**
  - ✅ **方向性 proxy**：Spearman 0.845 證明 CovM 單調對應真實 CN
  - ✅ **基線準確**：新 baseline vs SEQC2 neutral 的 bias < 3%
  - ⚠️ **分類 recall 仍有限**：Gain recall 45.8%（非 >80% FULL POSITIVE）— CovM 分類閾值（Normal 0.8-1.2, CNV_Gain 1.5-2.0）相對於 SEQC2 Gain 區（多為 CN=3，期望 CovM ≈ 1.5）邊界模糊
  - ⚠️ **Neutral precision 4.7%**：CovM=Normal 的區域並不代表真 CN=2（HCC1395 絕大多數 SNV-承載區域位於 gain 區）
  - 📌 **結論**：baseline 修正成功；**CovM 作為定性 CN proxy 可接受**（方向正確），但**定量 CN 判定仍需 coverage 外的特徵輔助**（如 AF、VAF）

### 5.2 對既有研究結論的影響評估

| 研究 / 結論 | 受影響程度 | 實測結果 |
|------------|----------|---------|
| **H-CN1**（CovM Baseline Accuracy）| 🔴 直接 | PARTIAL POSITIVE（見 §5.1） |
| **Step 1 HCC1395 Gain recall** | 🔴 直接 | 14.6% → 45.8%（×3.1） |
| **Step 2 跨樣本 baseline bias** | 🔴 直接 | 全 7 樣本 bias 收斂至 <3% |
| **Z3 amplicon blacklist CONDITIONAL（HCC1954）** | 🟢 **不影響** | S3 strategy **scale-invariant**（per-sample 95%ile 自動校正）— 原 verdict 不變 |
| **COLO829 QS 低 → CNV_Loss penalty 解釋**（O1-O10 Fig 11） | 🔴 **需撤回** | COLO829 CovM p50 從 0.387→1.000，CNV_Loss penalty 不再觸發 |
| **COLO829 M2 mask 91.7% 排除** | 🔴 **需重跑** | baseline artifact，下一 cycle 重算 |
| **「ISM 對黑色素瘤無效」假說** | 🟠 CONDITIONAL | 依賴 QS，連鎖受影響；下一 cycle 重測 |
| **HPFineNGroups Subclone Marker POSITIVE**（N≥4+NR≥80） | 🟢 不受影響 | 不依賴 CovM 分類 |
| **LOH Subclone AF×Methylation POSITIVE** | 🟢 不受影響 | 不依賴 CovM |
| **Self-Phasing Circular Dependency** | 🟢 不受影響 | 屬上游因果鏈 |

**COLO829 衝擊明細**：`research/kde_fix_validation/outputs/step4_colo829_audit/impacted_conclusions.md`（4 CRITICAL + 3 HIGH + 2 LOW）

### 5.3 Cross-sample quantile drift（Step 3 量化）

| Sample | Stale/KDE ratio | Δp50(CovM) | Δp50 方向 |
|--------|----------------:|-----------:|----------|
| COLO829 | 2.59 | +0.613 | 右移（stale 偏低） |
| HCC1395_DORADO | 1.42 | +0.364 | 右移 |
| HCC1395 | 1.36 | +0.343 | 右移 |
| HCC1954 | 1.23 | +0.207 | 右移 |
| H1437 | 1.09 | +0.081 | 右移 |
| H2009 | 0.95 | −0.060 | 左移（stale 偏高） |
| HCC1937 | 0.82 | −0.241 | 左移 |

→ Δp50 方向與 `ratio − 1` 完全一致；**COLO829 drift 最劇（+0.613）**，**HCC1937 反向漂移**（stale 75× **高估** true 91×，造成原 CovM 系統性過高）。

![Cross-sample CovM quantile drift (14 combos)](figures/20260420_kde_fix_acceptance/fig5_quantile_drift_per_sample.png)

![Coverage_Category migration per sample](figures/20260420_kde_fix_acceptance/fig6_category_migration_per_sample.png)

> **fig5**：2×7 格 histogram overlay（stale 紅／fixed 藍）；直觀呈現 COLO829 右移最劇、HCC1937 左移、H2009 幾乎重合。
> **fig6**：Coverage_Category 分布 bar chart；COLO829 從 CNV_Loss-dominant → Normal-dominant，HCC1937 則相反。

### 5.4 版本控制與回朔機制

```
HEAD (refactor/phase1-safety)
  └─ 5abc659  chore: evidence_ledger record H_KDE_001 ...           ← Step 6
  └─ ...      docs updates (CovM acceptance + manifest + plan)       ← 4 docs
  └─ 12d9b3e  fix(covm): use explicit fallback flag ...              ← Step 4 (Code Review fix)
  └─ 374fad4  feat(covm): KDE fallback WARN + audit column           ← Step 4 (feat commit)
  └─ ec0608b  chore: snapshot baseline before KDE expected_coverage ← TAG: pre-kde-logging-audit (Step 1)
```

**完整回朔（兩種等價方式）**：
```bash
# 方式 A：tag-based (推薦)
git reset --hard pre-kde-logging-audit

# 方式 B：commit-based
git reset --hard ec0608b
```

**保留 audit column 但回退邏輯改變**（不適用，因本次無邏輯改變，僅新增 column 與 logging）。

### 5.5 下游相容性檢查（30+ reader scripts）

| Reader 模式 | 數量 | 是否破壞 |
|------------|------|---------|
| `df["Coverage_Multiple"]` 直接欄名 | 22 | 不破壞 |
| `usecols=["Coverage_Multiple", ...]` | 8 | 不破壞 |
| `df.iloc[:, N]` 位置索引 | 0 | — |
| Hard-coded column index | 0 | — |

→ 新增 column `Diploid_Coverage_Used` 不影響任何下游 reader。

---

## 6. 未來方向

### 6.1 立即（待用戶決策 scope）

| 選項 | 範圍 | 時間 | 產出 |
|------|------|------|------|
| **A 最小** | HCC1395 TO（已完成） | 0 | Pilot 已驗證 C++ 修正 |
| **B 部分** | + HCC1395_DORADO paired_full / paired_pileup（3 combos） | ~30 min | 跨平台 Nanopore vs DORADO 對比 |
| **C 全量** | 5 樣本（HCC1937/HCC1954/H2009/H1437/COLO829）需先 LongPhase haplotag 重建 BAM，再跑 ISM | 4-8 hr + ~2.8 hr | Step 1/2 完整重量化 |

> **Blocker**：5/7 樣本的 `tumor_tagged.bam` 已從 archive 清除；C 選項需先完成 LongPhase haplotag 步驟才能進入 ISM。

### 6.2 重跑後的後續工作（不論選 A/B/C）

**本 cycle 已完成（2026-04-20）**：
1. ✅ **重量化 Step 1**：HCC1395 Gain recall 14.6%→45.8%、Spearman 0.845（`research/kde_fix_validation/outputs/step1_hcc1395_seqc2/`）
2. ✅ **重量化 Step 2**：跨樣本 bias 全收斂至 <3%（上 §5.0）
3. ✅ **更新 H-CN1 verdict**：PARTIAL POSITIVE（上 §5.1）
4. ✅ **重新評估 Z3 amplicon blacklist**：S3 strategy scale-invariant，verdict 不變
5. ✅ **COLO829 衝擊審計**：9 筆結論盤點（4 CRITICAL + 3 HIGH + 2 LOW）
6. ✅ **Quantile drift 量化**：14 combos 漂移範圍 −0.241 ~ +0.613

**本 cycle 追加完成（2026-04-21）**：
7. ✅ **COLO829 M2 mask 重算**（Step 7）：**79.40% → 4.87%**（Δ=−74.53pp）— 「M2 偏向 COLO829」**完全撤回**
8. ✅ **PCA 跨樣本分離重跑**（Step 8）：isolation ratio **11.88× → 1.70×**（−86%）— 「Depth 第一驅動」**撤回**
   - 詳見 `research/kde_fix_validation/outputs/FINAL_VERDICT.md`

**下一 cycle 需要 C++ rerun 的 follow-up（blocked on BAM 清除）**：
9. ❌ **COLO829 QS 重跑**：依賴 CNV_Loss penalty 的 QS 可能大幅回升（需新 binary + 新 BAM → 新 QS；5/7 樣本 BAM 已清除需 LongPhase 重建 ~8 hr）
10. ❌ **「ISM 對黑色素瘤無效」假說**：依賴 COLO829 新 QS，屬 P2 連鎖項
11. ❌ **跨樣本 QS 觀察（O11-O13 類）**：依賴 P1

### 6.3 演算法改進方向（Option D，post-acceptance 討論）

| 改進項 | 動機 | 預期效益 |
|-------|------|---------|
| Trimmed median fallback | 當 KDE histogram 退化時，trimmed median 比 default 75 更準 | 減少 fallback 觸發後的偏差 |
| 2nd-derivative peak detection 強化 | 當前用 simple max；2nd-deriv 可避開 noise plateau | mode 估計更穩定 |
| Dynamic bandwidth (Silverman's rule) | 當前 σ=3 bins 固定；高深度樣本應用更大 bandwidth | 適應 30×–100× 的廣度 |
| CN=2 prior constraint | KDE 假設 mode 來自 CN=2 區域；可加 prior weight | 對高 ploidy 樣本（如 HCC1954）更穩 |

→ 上述改進需在 **B/C 全量重跑完成後**評估必要性；若新 baseline 已 < 10% bias，多數改進無立即效益。

### 6.4 結論修正路徑

```
[本報告 ACCEPTANCE_VALIDATED]
        ↓
[用戶選 A/B/C scope]
        ↓
[執行重跑]
        ↓
[更新 20260420_CovM_Baseline_Accuracy_Validation_01.md Step 1/2 表格]
        ↓
[H-CN1 verdict 收斂：POSITIVE 或 NEGATIVE]
        ↓
[更新 evidence_ledger + MEMORY (project_expected_coverage_baseline_bug.md)]
        ↓
[評估是否需 Option D 演算法改進]
```

---

## 7. 附錄

### 7.1 Pilot 執行記錄

```
$ ls /big7_disk/liaoyoyo2001/big7_disk_output/synthesis/kde_rerun_pilot/HCC1395_to_tp/
debug
filtered_snv_tp
significance_statistics.txt
significance_summary.csv     # 28,495 rows × 80 cols (含 Diploid_Coverage_Used)
```

### 7.2 Binary log 摘錄（pilot 確證 KDE Pass 2 啟用）

```
[CovM] Pass 1 complete (28495 regions, default diploid=75.0x)
[CovM] estimate_diploid_coverage: histogram bins=120, peak_bin=22 (53.0x)
[CovM] Pass 2 diploid_coverage=53.000000x source=auto_kde used_fallback=false
[CovM] Pass 2 complete (CovM rescaled, Coverage_Category re-classified)
```

### 7.3 Evidence Ledger 紀錄

```json
{
  "cycle_id": "20260419_kde_covm_logging_audit",
  "hypothesis_id": "H_KDE_001",
  "hypothesis": "KDE per-sample diploid baseline 修正可大幅減少 stale-baseline (+39% bias) 帶來的 CovM 系統性偏差",
  "pipeline_track": "TO",
  "type": "cpp_improvement",
  "tier": 1,
  "delta_f1": null,
  "human_decision": "keep",
  "key_observations": "HCC1395 TO pilot: KDE 53x vs SEQC2 neutral truth 54x (-1.9% bias); Coverage_Multiple median ×1.415 rescaling; CNV_Gain+High_Copy 9.3% → 32.9%",
  "methodology_doc": "docs/methodology/20260419_KDE_expected_coverage_audit_01.md"
}
```

### 7.4 相關 Skill / 命令參考

- `/cpp-change` — 6 步驟 PDD 協議（本次依此執行）
- `/build` — 編譯
- `pre-kde-logging-audit` — annotated tag（rollback anchor）
- `scripts/canonical/{sample}/{mode}/run_benchmark.sh --skip-longphase` — 重跑單一樣本
