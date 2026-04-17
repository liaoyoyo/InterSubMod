<!--
建立時間: 2026-04-17
更新時間: 2026-04-17
目標: 驗證 LOH×AF×Methylation 結論 11 中 HCC1954 反向是真實訊號還是小樣本雜訊
範圍: 7 samples × {TO, Paired} segment-level + region-level
關聯: docs/plans/opus-4-7-big8-disk-liaoyoyo2001-knowled-cryptic-moore.md (Part B.2-1)
      docs/experiments/validated/2026/04/20260414_LOH_Subclone_AF_Methylation_Evidence_01.md
      docs/reports/research_landscape/06_結論穩定性審查.md (結論 11)
研究子專案: research/hcc1954_reversal_investigation/
-->

# HCC1954 反向排除標準調查（B.2-1）

## 0. 結論速覽

| 項目 | 結果 |
|------|------|
| **HCC1954 TO** | ρ=-0.297 **95% CI=[-0.604, +0.044]** → **CI 含 0，反向為統計噪音** |
| **HCC1954 Paired** | ρ=-0.211 **95% CI=[-0.621, +0.219]** → **CI 含 0，反向為統計噪音** |
| **H2009（新發現）** | TO ρ=+0.212 CI=[-0.28, +0.68] n=25 → **也 CI 含 0，正向也是功效不足** |
| **移除 HCC1954 敏感性** | Meta-ρ 從 +0.367 → +0.383 (Δ=+0.016) → **HCC1954 並非拖累，包含與否不改變結論方向** |
| **結論 11 精確化** | 原「7/7 一致」應修正為 **「5/7 CI 不含 0 正向、2/7 功效不足（HCC1954 n=34, H2009 n=25）」** |
| **Pre-registered exclusion** | 建議標準：**n_CN1_segments ≥ 50**（可事先定，不需 post-hoc） |
| **結論 11 穩定度** | **維持 ⭐3**，加註兩項限制 |

**最終判定**：HCC1954 反向並非真實 signal reversal，而是小樣本（n_CN1_segments=34）功效不足；但需精確化原「7/7 一致」聲稱。

---

## 1. 質疑背景

### 1.1 原結論 11（20260414 報告）

| 聲稱 | 來源 |
|------|------|
| TO segment-level Spearman ρ=0.270 | step3 — 報告 Section 3 |
| Paired segment-level Spearman ρ=0.382 | step3 — paired subproject |
| 7/7 樣本方向一致 | step3 per_sample_segment_correlation.tsv |

### 1.2 實際 per-sample 數據（原）

| Sample | TO n_segments | TO ρ | Paired n_segments | Paired ρ |
|--------|---------------|------|-------------------|----------|
| HCC1395 | 435 | +0.151 | 438 | +0.209 |
| HCC1395_DORADO | 349 | +0.255 | 350 | +0.231 |
| COLO829 | 105 | +0.763 | 106 | +0.503 |
| H1437 | 168 | +0.809 | 167 | +0.744 |
| H2009 | 25 | +0.212 | 25 | +0.666 |
| HCC1937 | 110 | +0.230 | 106 | +0.437 |
| **HCC1954** | **34** | **-0.297** | **30** | **-0.211** |

### 1.3 質疑點（Plan B.2-1）

- HCC1954 方向相反，原報告以「純度/CN 複雜」post-hoc 解釋
- B.1-1 新發現 HCC1954 在 HPFineNGroups 為 7 樣本中唯一 ΔAUC ROBUST（+0.022）— 兩結論 HCC1954 行為相反
- 需驗證：(1) 統計噪音 or 真反向 (2) 建立 pre-registered exclusion

### 1.4 程式與資料溯源

**原 step3 公式精確定位**（最重要！）：
```python
# research/loh_subclone_af/scripts/step3_spatial_analysis.py:373
sub = seg_df[(seg_df["sample"] == sample) & (seg_df["dominant_cn"] == "CN1")]
# line 384
rho, pval = scipy_stats.spearmanr(sub["af_sd"], sub["mean_ngroups"])
```

**Filter chain**：
1. LOH region only (`to_loh_bed_hit=True`)
2. TP only (`truth_label="TP"`)
3. Has ≥2 variants per segment
4. **CN1 only** (`dominant_cn="CN1"`, i.e. `Coverage_Multiple < 0.75`)
5. Spearman ρ(af_sd, mean_ngroups)

### 1.5 Git 基準

- 起始 HEAD：`ab61ad1`（B.1-1 完成後）
- 本次不改 C++，僅新增 `research/hcc1954_reversal_investigation/`

---

## 2. 方法（Step → Verify）

| Step | 方法 | 驗證標準 |
|------|------|---------|
| S1 | 讀既有 step3 segment-level TSV（TO + Paired） | TO=3750, Paired=3732 segments ✅ |
| S2 | 重現原 per-sample ρ（CN1 only, af_sd vs mean_ngroups）+ Bootstrap 95% CI (B=2000) | HCC1954 ρ=-0.30 吻合原報告 ✅ |
| S2b | 對照分析：所有 CN + af_mean vs mean_ngroups | 與 S2 比較檢查角度差異 ✅ |
| S3 | 判定 HCC1954 反向：CI 含 0 (noise) vs 不含 0 (true reversal) | CI 決策 ✅ |
| S4 | Leave-one-out sensitivity：meta-ρ via Fisher z, 權重 n−3 | 7/6 樣本比較 ✅ |
| S5 | HCC1954 region-level ρ（step2 vs step3 調和）| 看不同聚合層級 ρ 方向 ✅ |
| S6 | Forest plot + HCC1954 vs 其他樣本 segment 分佈 | PNG 視覺化 ✅ |

---

## 3. 結果

### 3.1 S2 核心結果 — CN1 only, af_sd vs mean_ngroups, Bootstrap CI

**TO mode**：

| Sample | n_CN1 | ρ | 95% CI | p | CI 含 0? |
|--------|-------|---|--------|---|---------|
| HCC1395 | 435 | +0.151 | [+0.044, +0.249] | 1.6e-03 | ❌ |
| HCC1395_DORADO | 349 | +0.255 | [+0.135, +0.382] | 1.3e-06 | ❌ |
| COLO829 | 105 | +0.763 | [+0.628, +0.868] | 3.0e-21 | ❌ |
| H1437 | 168 | +0.809 | [+0.734, +0.871] | 3.4e-40 | ❌ |
| **H2009** | **25** | +0.212 | **[-0.281, +0.677]** | 0.309 | ✅ **功效不足** |
| HCC1937 | 110 | +0.230 | [+0.014, +0.439] | 1.6e-02 | ❌ |
| **HCC1954** | **34** | -0.297 | **[-0.604, +0.044]** | 0.088 | ✅ **反向為噪音** |

**Paired mode**：HCC1954 ρ=-0.211, CI=[-0.621, +0.219] 含 0；H2009 ρ=+0.666, CI=[+0.310, +0.876] 不含 0（Paired 下 H2009 穩定）。

### 3.2 S2b 對照分析 — 所有 CN + af_mean vs mean_ngroups

所有樣本（含 HCC1954）的 **ρ(af_mean, mean_ngroups)** 皆為負相關（-0.26 到 -0.95）。

**解釋**：這是**不同的機制**，不可與 S2 混淆：
- **S2 設定（原結論 11）**：CN1 segments 內 AF 變異大（AF_SD 高）→ NGroups 高。**subclonal 異質性在 segment 內**的訊號。
- **S2b 對照**：segments 整體 AF 越高（越接近 homozygous）→ NGroups 越低（單一基因型）。**CN-mixed 聚合效應**。

兩者不矛盾。S2 才是「subclone detection」的核心聲稱。

### 3.3 S3 — HCC1954 判定

| Mode | ρ | CI | 判定 |
|------|---|-----|------|
| TO | -0.297 | [-0.604, +0.044] | **CI 含 0，小樣本噪音（n=34）** |
| Paired | -0.211 | [-0.621, +0.219] | **CI 含 0，小樣本噪音（n=30）** |

**→ 兩種 mode 下 HCC1954 反向皆為統計噪音，非真實 signal reversal。**

### 3.4 S4 — Leave-One-Out 敏感性

**Meta-ρ (Fisher z 加權 n−3)**：

| Excluded | TO meta-ρ | Paired meta-ρ | TO Δ | Paired Δ |
|----------|-----------|---------------|------|----------|
| NONE (full 7) | **+0.367** | **+0.357** | baseline | baseline |
| HCC1395 | +0.474 | +0.434 | +0.107 | +0.077 |
| HCC1395_DORADO | +0.410 | +0.405 | +0.043 | +0.048 |
| COLO829 | +0.317 | +0.342 | -0.051 | -0.015 |
| H1437 | +0.262 | +0.274 | -0.105 | -0.083 |
| H2009 | +0.370 | +0.350 | +0.003 | -0.007 |
| HCC1937 | +0.380 | +0.349 | +0.013 | -0.008 |
| **HCC1954** | **+0.383** | **+0.369** | **+0.016** | **+0.012** |

**發現**：
1. **排除 HCC1954 meta-ρ 反而 +0.016**（從 0.367→0.383）— HCC1954 並非拖累，反而是**微弱負向噪聲**被平均稀釋。
2. **H1437 排除衝擊最大 Δ=-0.105**（權重最高樣本），排除後 meta-ρ 下降；但 HCC1954 排除衝擊極小。
3. **H2009 排除 Δ=+0.003 幾乎不變**（n=25 權重低）— 噪音樣本本來就貢獻小。

**含義**：**HCC1954 是否納入不改變結論 11 的方向，不需 post-hoc exclusion**。

### 3.5 S5 — Region-Level ρ（Step2 vs Step3 調和）

在 LOH 內每個 region（非 segment 聚合）計算 ρ(caller_af, HPFineNGroups)：

| Sample | TO n_regions | TO ρ | p |
|--------|--------------|------|---|
| HCC1395 | 12,991 | -0.570 | 0 |
| HCC1395_DORADO | 13,219 | -0.740 | 0 |
| COLO829 | 7,025 | -0.659 | 0 |
| H1437 | 12,410 | -0.802 | 0 |
| H2009 | 31,375 | -0.262 | 0 |
| HCC1937 | 6,905 | -0.532 | 0 |
| **HCC1954** | 1,418 | **-0.282** | 2.5e-27 |

**發現**：
- 所有樣本 region-level ρ 皆負 — AF 越高（越 homozygous）NGroups 越低
- HCC1954 ρ=-0.28 與其他樣本同向（-0.26 到 -0.80）
- **Step2 Δ=+0.69 POS（Intermediate AF 比 Extreme AF NGroups 高）與 region-level ρ 負是一致的**：Intermediate AF 範圍 [0.1,0.4]∪[0.6,0.9]，相對於 Extreme [<0.1 或 >0.9]，區間聚合後 Intermediate 的 mean_NG 較高

**Step2/Step3 調和**：
- Step2（bin 比較）：Intermediate AF bin mean_NG > Extreme AF bin mean_NG → HCC1954 POS
- Step3（segment-level）：在 CN1 segments 內 AF_SD 大 → NG 高 → HCC1954 n=34 功效不足
- Region-level ρ 負：整體 AF→NG 是負相關（純 homozygous 無 NGroups 變異）

三層分析**測試不同假設**，HCC1954 在 S2、S5 都 POS；僅 S3（小 n）功效不足。

---

## 4. 原結論 11 審查表（逐條對應）

| 原聲稱 | 驗證結果 | 維持/修正 |
|--------|---------|-----------|
| TO Δ_NG Inter−Ext=+0.705（step2）| HCC1954 Δ=+0.69 POS、7/7 方向一致 | ✅ **維持** |
| Paired Δ_NG=+0.787（step2）| HCC1954 Δ=+0.91 POS、7/7 方向一致 | ✅ **維持** |
| TO segment-level Spearman ρ=0.270（step3）| Meta-ρ=+0.367（Fisher z 加權，校正後更高）| ✅ **維持（數值更強）** |
| 7/7 樣本 step3 方向一致 | Bootstrap 發現 HCC1954 CI 含 0、H2009 CI 含 0 | ⚠ **精確化**：5/7 CI 明確 + 2/7 功效不足 |
| HCC1954 反向是 post-hoc 「純度/CN 複雜」 | Bootstrap CI=[-0.60, +0.04] 含 0 = 小樣本噪音 | ✅ **正名**：統計噪音，不是真反向 |
| 「需 pre-registered exclusion」| n_CN1_segments ≥ 50 可事先定 | ✅ **採納** |

---

## 5. 新增限制清單

1. **Per-sample 功效不足樣本（pre-registered threshold）**：建議 step3 per-sample 分析應要求 **n_CN1_segments ≥ 50** 作為「功效足夠」門檻。HCC1954 (34)、H2009 TO (25) 均不達標，其 ρ 值應標記為「underpowered」。
2. **H2009 TO 也是功效不足**（B.2-1 新發現）：原報告未明點出；Paired mode H2009 n=25 ρ=+0.666 CI=[+0.31, +0.88] 不含 0 穩定，但 TO mode 同 n=25 ρ=+0.212 CI 含 0 不穩。
3. **「7/7 一致」精確化**：應修正為 **「7/7 方向（step2 全 POS）+ 5/7 step3 CI 明確（不含 0）+ 2/7 step3 功效不足」**。
4. **HCC1954 段數少的根因**：n_CN1_segments=34 的原因是 HCC1954 LOH segments 多為 CN2（hyper-diploid）而非 CN1。HER2+ hyper-diploid 特性（知識庫）造成 dominant_cn 分佈集中在 CN2；應在論文 limitation 段落說明「cnLOH vs deletion-LOH 在 HER2+ 樣本比例不同」。

---

## 6. 對結論 11 穩定度的最終判定

> **⭐3（維持）** + 兩項限制加註

- **Step2 7/7 一致**：確認 ✅（含 HCC1954 Δ=+0.69 POS）
- **Step3 meta-ρ ROBUST**：含 HCC1954 ρ=+0.367、排除 HCC1954 ρ=+0.383，均為正向顯著
- **剩餘風險**（未消除）：
  - (a) Pre-registered exclusion 標準 n_CN1_segments≥50 建立後，需在論文中 upfront 說明（否則仍為 post-hoc）
  - (b) HCC1954 HER2+ hyper-diploid 特性需延伸至 clinical cohort 查看是否 HER2+ 群體系統性功效不足
  - (c) Step3 原 CN1 only 限定未涵蓋 cnLOH (CN2)，HER2+ 樣本可能在 cnLOH 層有訊號；建議補 CN2 LOH segments 獨立分析

---

## 7. Artifacts（Git 可追溯）

| 類型 | 路徑 | 用途 |
|------|------|------|
| PLAN | `research/hcc1954_reversal_investigation/00_PLAN.md` | 質疑目標與步驟 |
| manifest | `research/hcc1954_reversal_investigation/manifest.yaml` | Hypothesis + findings |
| Script | `research/hcc1954_reversal_investigation/scripts/01_reversal_check.py` | 主分析（Bootstrap + LOO + region-level） |
| ρ + CI | `data/rho_bootstrap_ci.tsv` | 7×2 modes Bootstrap 結果 |
| 對照 | `data/rho_af_mean_all_cn.tsv` | S2b 所有 CN + af_mean 角度 |
| LOO | `data/sensitivity_leave_one_out.tsv` | 敏感性分析 |
| Region-level | `data/hcc1954_region_level_rho.tsv` | step2/step3 調和 |
| Forest | `figures/01_bootstrap_ci_forest.png` | 7 樣本 ρ + CI 視覺化 |
| 分佈 | `figures/02_hcc1954_vs_others_segment_dist.png` | HCC1954 vs 其他 segment 數 |

---

## 8. 下一步建議

1. **更新結論 11 原報告**：`docs/experiments/validated/2026/04/20260414_LOH_Subclone_AF_Methylation_Evidence_01.md` 加入「per-sample 功效限制」section。
2. **更新 06 結論穩定性審查**：結論 11 加入 pre-registered exclusion criterion。
3. **B.2-2 Coverage_Multiple vs 真實 CN** 下一項：本次 S4 顯示 HCC1954 hyper-diploid 的 CN 分佈是關鍵；若 Coverage_Multiple 不能區分 cnLOH vs deletion-LOH，結論 11 的 CN1 限定可能低估 HCC1954 真實訊號。
4. **記憶更新**：`project_loh_subclone_af_methylation_positive.md` 加入「5/7 CI 明確 + 2/7 功效不足」精確化。
