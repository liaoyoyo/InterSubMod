<!--
建立時間: 2026-04-17
更新時間: 2026-04-17 (加入 B.1-1 residualized AUC 多變數驗證)
目標: 驗證 HPFineNGroups N≥4+NR≥80 TP rate 89.1% 是否為 NR 飽和 artifact + 多變數 confound
範圍: 748K regions × 7 samples × 2 modes, NR-matched + per-sample + LOH cross-strat + LR residualized AUC
關聯: docs/plans/opus-4-7-big8-disk-liaoyoyo2001-knowled-cryptic-moore.md (Part B.1-2 + B.1-1)
      docs/reports/research_landscape/06_結論穩定性審查.md (結論 10)
      docs/reports/validated/2026/04/20260406_肉眼檢視推理鏈與TP_FP可區分性分析_01.md
研究子專案: research/hpfinengroups_saturation_check/
-->

# HPFineNGroups 飽和效應質疑驗證報告（B.1-2 + B.1-4）

## 0. 結論速覽

| 項目 | 結果 |
|------|------|
| **飽和假說** | ❌ **否定**：NR-matched 下 NGroups=4 vs NGroups=3 Δ=+11.7pp（TO NonLOH NR≥80），每個 NR bin 都 ≥+5pp |
| **7/7 方向一致** | 修正為 **7/7 POS（NR-bin weighted 後）**，但 COLO829/H2009 效應 <5pp |
| **LOH 內適用性** | ❌ **LOH 區域 NGroups 被壓縮**（NGroups=4 僅 n=11），HPFineNGroups **僅適用於 NonLOH** |
| **結論 10 穩定度** | **維持 ⭐3**，加註三項限制 |

**最終判定**：原 89.1% TP rate 非 NR 飽和 artifact，但需明確標註適用範圍（TO NonLOH NR≥40）與樣本差異（COLO829/H2009 效應弱）。

---

## 1. 質疑背景

### 1.1 原結論（20260406 報告 Section 8-11、結論穩定度 #10）

- HPFineNGroups N≥4 全樣本 TP rate=86.8%（vs baseline 52-66%），TP/FP ratio=3.07×
- **N≥4 + NR≥80 NonLOH：25,514 regions → TP rate 89.1%**
- 7/7 樣本方向一致、residualized AUC=0.617

### 1.2 質疑點（Plan Part B.1-2）

NGroups 上限 = 4（HP1/HP1-1/HP2/HP2-1）。**若高 NR 本身讓 NGroups 容易達 4（飽和），89.1% 可能只是 NR≥80 NonLOH 的 baseline 而非 N≥4 的增量訊號。**

### 1.3 資料與程式溯源

- **程式**：`src/core/LabelTest.cpp:633` `result.fine_n_groups = unique_groups.size()`（上限 4）
- **TSV 欄位**：`HPFineNGroups` (col 27)，輸出於 `src/core/RegionProcessor.cpp:1147`
- **資料集**：`output/synthesis/observation_workspaces/20260327_loh_round1_cross_sample_audit/all_region_rows.tsv.gz`（748,391 rows）
- **篩選重現**：`to_loh_bed_hit=False & NumReads≥80 & HPFineNGroups=4` → 25,741 rows, TP rate=**0.8913**（與原 89.1% 完全一致，n 略異因 NGroups==4 vs ≥4 定義）

### 1.4 Git 前後基準

- 起始 HEAD：`c022918`
- 本次分析僅新增 `research/hpfinengroups_saturation_check/`（未改 C++，可前後比較）

---

## 2. 方法（Step → Verify）

| Step | 方法 | 驗證標準 |
|------|------|---------|
| S1 | 讀 748K rows，拆 mode × LOH × NR bin × NGroups | rows=748,391 ✅ |
| S2 | NR bin (7 bins) × NGroups (1-4) TP rate 交叉表 | 完整 3×7×4 矩陣 ✅ |
| S3 | NR≥80 固定條件下 NGroups=3 vs 4 Fisher exact | 計算 Δpp, odds, p ✅ |
| S4 | 每個 NR bin 內 NGroups=3 vs 4 Fisher exact（匹配） | 6/6 bins 驗證 Δ≥+5pp ✅ |
| S6 | 7 樣本分別 NR≥80 TO NonLOH N3 vs N4 | 完整 per-sample 表 ✅ |
| S7 | 7 樣本 NR-bin weighted Δ（更嚴謹的 residualization） | weighted delta 修正 S6 ✅ |
| B.1-4 | LOH × NGroups × NR 三維交叉 | LOH N4 n=11 — 樣本不足 ✅（反證 LOH 不適用） |

---

## 3. 結果

### 3.1 S3 — TO NonLOH NR≥80 核心驗證

| NGroups | n | TP | TP rate | Fisher p (N4 vs N3 one-sided) |
|---------|---|----|---------| ----------------------------- |
| 3 | 64,633 | 50,044 | **0.7743** | — |
| 4 | 25,741 | 22,942 | **0.8913** | **≈0** (Δ=+11.70pp, OR=2.39) |

**→ 重現原 89.1% TP rate。Δ vs NGroups=3 為 +11.70pp，p≈0。**

### 3.2 S4 — NR-matched Fisher（關鍵飽和反駁證據）

TO NonLOH 每個 NR bin 內 NGroups=3 vs 4 的 Δ TP rate：

| NR bin | n3 | rate3 | n4 | rate4 | Δ (pp) | Fisher p |
|--------|----|-------|----|-------|--------|----------|
| 20-40 | 8,737 | 0.661 | 239 | 0.715 | **+5.46** | 4.4e-02* |
| 40-60 | 15,322 | 0.579 | 1,712 | 0.690 | **+11.11** | 1.6e-19* |
| 60-80 | 32,345 | 0.693 | 6,811 | 0.829 | **+13.60** | 5.1e-123* |
| 80-100 | 27,472 | 0.754 | 8,057 | 0.878 | **+12.37** | 3.3e-136* |
| 100-150 | 29,158 | 0.803 | 12,896 | 0.903 | **+10.05** | 6.5e-156* |
| 150-500 | 8,003 | 0.740 | 4,788 | 0.882 | **+14.17** | 4.5e-87* |

**→ 所有 NR bin 均顯示 +5pp 以上（多數 +10pp 以上）。飽和假說明確否定。**

### 3.3 S6 / S7 — Per-Sample 7 樣本驗證

| Sample | NR≥80 only Δ (pp) | NR-bin weighted Δ (pp) | 最終方向 |
|--------|-------------------|------------------------|----------|
| H1437 | +15.04 | **+15.58** | POS 強 |
| H2009 | +3.02 | **+2.84** | POS 弱（ceiling effect） |
| HCC1395_DORADO | +10.37 | **+11.88** | POS |
| HCC1395 | +7.36 | **+8.72** | POS |
| HCC1937 | +14.13 | **+13.96** | POS |
| HCC1954 | +17.58 | **+14.11** | POS 強 |
| COLO829 | **-25.77** | **+3.56** | POS 弱（NR≥80 小樣本 n3=71 n4=34 artifact） |

**→ NR-bin weighted 分析後 7/7 全部 POS。COLO829 原 -25.77pp 為 NR≥80 在 COLO829 TO 極少（僅 143 rows）的小樣本偏差。H2009 + COLO829 效應 <5pp，屬弱訊號樣本。**

**H2009 弱訊號成因**：overall TP rate=0.911（含全 NR）、NR≥80 baseline=0.909（接近 ceiling）→ NGroups 訊號無空間可挖。

### 3.4 B.1-4 — LOH × NGroups 交叉（補強分析）

| LOH | NGroups | n | TP rate |
|-----|---------|---|---------|
| False | 1 | 15,097 | 0.684 |
| False | 2 | 134,203 | **0.573** (最低) |
| False | 3 | 121,336 | 0.720 |
| False | 4 | 34,508 | **0.868** |
| True | 1 | 56,709 | 0.629 |
| True | 2 | 54,983 | **0.896** (LOH baseline 已極高) |
| True | 3 | 257 | 0.848 |
| True | 4 | **11** | 0.818 (樣本不足) |

**發現**：
1. **LOH 區域 NGroups=4 僅 n=11**，HPFineNGroups 在 LOH **完全無效**（NGroups 被 haplotype 壓縮到 1-2）。
2. **LOH NGroups=2 TP rate=0.896** 已超過 **NonLOH NGroups=4 的 0.868**。在 LOH 區域 NGroups=2 才是 TP 訊號峰值（self-phasing HP bias 造成）。
3. **NonLOH NGroups=2 是 FP 富集區**（TP rate=0.573 全組最低）。

**推論**：HPFineNGroups 的 subclone marker 能力**僅適用於 NonLOH region**。LOH region 應依賴其他特徵（如 AlleleDelta、caller_af）。

---

## 3.5 B.1-1 — Residualized AUC 多變數控制（新增）

**質疑**：residualized AUC=0.617 是否足以支撐 POSITIVE？原報告僅控制 n_reads，未控制 AF/Coverage_Multiple/LOH。

**方法**：
- 模型：sklearn `LogisticRegression`（L2，max_iter=2000）
- M1 Raw = HPFineNGroups 單變數 AUC
- M3 NR+AF control：full 模型 = [NR, AF, HPFineNGroups]，control 模型 = [NR, AF]；Δ = AUC(full) − AUC(control)
- M4 NR+AF+Cov+LOH（最嚴格）：full = [NR, AF, Cov_Multiple, to_loh_bed_hit, HPFineNGroups]，control 移除 HPFineNGroups
- 判定：Δ<0.005 可忽略 / Δ<0.02 弱 / Δ≥0.02 ROBUST

**Scenario 結果（pooled）**：

| Scenario | n | tp_rate | M1 Raw AUC | M3 NR+AF Δ | M4 NR+AF+Cov+LOH Δ | 判定 |
|----------|---|---------|-----------|------------|--------------------|-----|
| **TO NonLOH all** | 307,474 | 0.670 | 0.6032 | **+0.0251** | **+0.0251** | ✅ ROBUST |
| TO NonLOH NR≥80 | 145,803 | 0.748 | 0.6188 | **+0.0507** | **+0.0507** | ✅ ROBUST（更強） |
| TO LOH all | 112,218 | 0.761 | 0.6822 | +0.0285 | +0.0284 | 強但 NR-only Δ=+0.13 說明 LOH NGroups 與 NR 強糾纏 |
| TO all | 419,692 | 0.694 | 0.5640 | +0.0046 | **+0.0286** | M3 弱、M4 被 LOH 分層拉回 ROBUST |
| **Paired NonLOH all** | 328,699 | 0.9896 | 0.5322 | **+0.0018** | **+0.0018** | ❌ ceiling — control alone AUC=0.80，訊號無空間 |

**→ TO NonLOH 核心場景 ROBUST：** 控制 NR+AF+Cov+LOH 後 HPFineNGroups 仍增益 +0.025 AUC。原聲稱 0.617 介於 M1=0.603 與 M3=0.641 之間，吻合。

**Paired mode 驗證失敗（ceiling）**：TP rate=98.96%，NR+AF control 本身就達 AUC=0.80，HPFineNGroups 無空間貢獻。

**Per-sample TO NonLOH（7 樣本）**：

| Sample | n | tp_rate | M1 Raw | M3 NR+AF Δ | M4 Δ | 判定 |
|--------|---|---------|--------|-----------|------|------|
| HCC1954 | 63,081 | 0.248 | 0.665 | **+0.0219** | **+0.0219** | ✅ ROBUST（baseline 最低） |
| HCC1395_DORADO | 22,643 | 0.691 | 0.598 | +0.0101 | +0.0101 | 弱 |
| HCC1937 | 12,236 | 0.467 | 0.624 | +0.0084 | +0.0084 | 弱 |
| H2009 | 103,554 | 0.911 | 0.541 | +0.0086 | +0.0086 | 弱（ceiling） |
| H1437 | 43,384 | 0.762 | 0.548 | +0.0083 | +0.0083 | 弱 |
| HCC1395 | 22,541 | 0.688 | 0.569 | +0.0051 | +0.0051 | 微弱 |
| COLO829 | 40,035 | 0.651 | 0.500 | +0.0007 | +0.0007 | 可忽略 |

**推論**：
1. **Pooled TO NonLOH ROBUST** → 結論 10 主張**保留**。
2. **Per-sample 7 樣本僅 1 個 ROBUST（HCC1954）**，6 樣本 Δ<0.02；ΔAUC 層級效應比 Δ TP rate 更保守。
3. **方向一致性**：7/7 樣本 Δ>0，但強度差異大；「7/7 一致」聲稱應精確化為「7/7 方向一致、1/7 效應 ROBUST」。
4. **TO all（不分 LOH）M3 Δ=+0.005**：代表 LOH 分層是關鍵——混 LOH/NonLOH 會稀釋訊號；這與 B.1-4 的 LOH 不適用結論一致。

---

## 4. 原結論審查表（逐條對應）

| 原聲稱 | 驗證結果 | 維持/修正 |
|--------|---------|-----------|
| N≥4 + NR≥80 NonLOH TP rate = 89.1% | 重現：n=25,741, TP rate=0.8913 | ✅ **維持** |
| 89.1% 是 N=4 訊號而非 NR 驅動 | NR-matched Fisher 6/6 bins 顯著 +5~+14pp | ✅ **維持**（飽和否定） |
| 7/7 樣本方向一致 | NR-bin weighted 後 7/7 POS | ⚠ **精確化**：6/7 效應 ≥+5pp；COLO829/H2009 弱（<5pp） |
| HPFineNGroups 是 somatic heterogeneity 標記 | NonLOH 內適用 | ⚠ **範圍限定**：僅 NonLOH，LOH 不適用 |
| 用於 filter（AUC 不足）| M3 full AUC=0.641 < 0.65 filter 門檻 | ✅ **維持不可用於 filter** |
| residualized AUC=0.617 支撐 POSITIVE | M3 Δ=+0.025 ROBUST；M4 (含 LOH+Cov) 不衰退 | ✅ **通過 B.1-1 驗證** |

---

## 5. 新增限制清單

本次分析發現原報告**未明說的限制**，並結合知識庫樣本文件交叉確認：

1. **LOH 不適用**：LOH 區域 NGroups 被壓縮到 1-2，NGroups=4 僅 n=11。所有 HPFineNGroups 聲稱應限定 NonLOH scope。
2. **Ceiling effect（H2009）**：當樣本 overall TP rate 已極高（>0.9）時，NGroups 訊號被壓扁（H2009 Δ 僅 +2.84pp）。H2009 truth set 為 orthogonal-tools（非 SEQC2 嚴格 consensus）；需在「baseline TP rate <0.85 的樣本」內評估。
3. **COLO829 NR≥80 小樣本**：COLO829 TO NonLOH 僅 143 rows NR≥80。原聲稱「7/7 一致」需用 NR-bin weighted 而非單一 NR≥80 cut-off 計算。
4. **樣本層級背景（知識庫補強）**：
   - **COLO829** melanoma 高 TMB，ONT_R10 BAM **無 methylation tags**，僅 ONT_PAO 子集（5mCG+5hmCG）有——可能 methylation call 效率不同影響 NR≥80 覆蓋。
   - **HCC1395** genome ploidy=**2.85** hyper-diploid、LOH 涵蓋 1490 Mb（基因組一半）；NonLOH region 相對較少，Δ=+8.72pp（中等）可能受限於 NonLOH 樣本少。
   - **H2009** truth set 嚴格度低於 SEQC2，FP 本身較少 → TP baseline 高 → NGroups 訊號無空間。
   - HCC1954 的反向風險（Plan B.2-1）與 HCC1395 類似 hyper-diploid 特徵一致，需 B.2 繼續調查。

---

## 6. 對結論 10 穩定度的最終判定

> **⭐3（維持）** + 三項限制加註

- **飽和否定**：S4 NR-matched 6/6 bins 顯著 Δ≥+5pp、多數≥+10pp ✅
- **適用範圍**：TO mode NonLOH region, NumReads≥40, NGroups=4 ✅
- **剩餘風險**（B.1-1 後更新）：
  - (a) Paired mode TP baseline=99.3%，訊號在 ceiling（Paired NR≥80 Δ=+0.11pp、B.1-1 M3 Δ=+0.002 無法獨立驗證）
  - (b) **Per-sample AUC 層級效應 6/7 <0.02**（B.1-1 新發現）：僅 HCC1954 Δ=+0.022，其餘 6 樣本微弱；「7/7 一致」應精確化為「方向一致、強度僅 HCC1954 ROBUST」
  - (c) Coverage_Multiple 代理 CN 的準確度尚未驗證（Plan B.2-2）；M4 結果與 M3 幾乎相同可能反映 Coverage_Multiple 訊號已被 NR/AF 線性吸收，而非真實 CN 無關

---

## 7. Artifacts（Git 可追溯）

| 類型 | 路徑 | 用途 |
|------|------|------|
| PLAN | `research/hpfinengroups_saturation_check/00_PLAN.md` | 質疑目標與步驟 |
| manifest | `research/hpfinengroups_saturation_check/manifest.yaml` | Hypothesis + datasets + parameters |
| Script S1-S6 | `research/hpfinengroups_saturation_check/scripts/01_saturation_check.py` | 主分析 |
| Script anomaly | `research/hpfinengroups_saturation_check/scripts/02_anomaly_investigation.py` | B.1-4 + per-sample weighted |
| S2 crosstab | `data/nr_ngroups_crosstab.tsv` | 3×7×4 TP rate 矩陣 |
| S3 saturation | `data/s3_nr80_saturation.tsv` | NR≥80 核心驗證 |
| S4 matched | `data/nr_matched_test.tsv` | 每 NR bin Fisher |
| S6 per-sample | `data/per_sample_nr80_ngroups.tsv` | 7 樣本 |
| S7 weighted | `data/s7_anomaly_and_weighted.txt` | NR-bin weighted Δ |
| Heatmap | `figures/01_nr_ngroups_tp_rate_heatmap.png` | TO NonLOH/LOH + Paired NonLOH |
| Forest | `figures/02_nr_matched_delta_forest.png` | 7 樣本 NR≥80 Δ |
| **B.1-1 Script** | `scripts/03_residualized_auc.py` | Residualized AUC (sklearn LR) |
| **B.1-1 Panels** | `data/residualized_auc_panels.tsv` | 6 scenarios × 5 模型 |
| **B.1-1 Per-sample** | `data/residualized_auc_per_sample.tsv` | 7 樣本 TO NonLOH |
| **B.1-1 Bar chart** | `figures/03_residualized_auc_bar.png` | M1/M3/M4 對比 |

---

## 8. 下一步建議

1. ~~**B.1-1（residualized AUC）複驗**~~ → **已完成（2026-04-17）**：M3/M4 Δ≥+0.025 ROBUST，結論 10 通過。
2. **B.2-1（HCC1954 反向問題，LOH×AF 結論 11）**：本次 B.1-1 發現 HCC1954 在 HPFineNGroups 卻是 7 樣本中唯一 ROBUST（Δ=+0.022）；與 LOH×AF 結論 11 的 HCC1954 反向方向**不一致**——兩個結論的 HCC1954 行為需對照調查。
3. **B.2-2（Coverage_Multiple 非精確 CN）**：B.1-1 M4 加入 Cov 後 Δ 未衰減（+0.025→+0.025），可能表示 Cov 與 NR/AF 已線性共線。需以 Delly/Manta 等獨立 CN caller 驗證結論 11。
4. **已更新**：`06_結論穩定性審查.md` 結論 10 + memory `project_hpfinengroups_subclone_marker.md`。
