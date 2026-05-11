---
title: 研究問題與假說
status: in_progress
last_updated: 2026-04-22
---

# 01. 研究問題與假說

## 1.1 背景脈絡（為什麼做這件事）

### 起點：v1 偏題

v1 研究（`20260421_NG_KDE_Rescaled_Multi_CN_Analysis_01.md`，已 superseded）在新 KDE baseline 下重跑了 20260414 的 ΔNG cross-section 觀察，發現 6/7 ΔNG 方向反轉。用戶指出此分析偏題，因為不直接回答「LOH 與 AF 是否能分出 TP 與 FP」的核心問題。

### v2 轉向：TP/FP 判別

v2 研究（`20260422_LOH_AF_KDE_TPFP_Discrimination_02.md`）重新框架為：
- 測試 LOH_Subtype × AF_class × CN-tier 是否能切分出 TP 富集區與 FP 富集區
- 在 HCC1395 TO mode（baseline TP=71.1%）為主測試場
- 建立 S1-S8 biology-defined scheme，測得 S3 Diploid Het=95.5% TP、S5 Combo=91.8% TP

### E1-E5 深度驗證（§11-14）

- E1：相對 caller baseline 的 TP:FP fold-improvement（S3=8.69×）
- E2：跨樣本穩健性（paired_full 全樣本 S3/S5 ≥97%）
- E3：統計力（36 STRONG + 10 NO_USE cells）
- E4：TP vs FP 特徵差（AlleleDelta Cliff's δ=-0.49 最強）
- E5：Sub-scheme 嘗試（S4d 最佳 fold=1.17×，其他 <1.17×）

### E6 綜觀突破（§17）

放棄「在 scheme 內加單特徵」框架，採用**聯合 5D cube**：
- LOH_Subtype (5) × AF_class (3) × cn_tier (5) × NG (1-4) × NR_band (3)
- 發現 17-cell envelope=96.1% purity + 3.7% recall（fold 10×），**Pareto-dominates** S3
- 28-cell envelope=90.4% purity + 7.4% recall（fold 4.73×），**Pareto-dominates** S5

---

## 1.2 核心問題（五層）

### Q1：baseline 比值驗證
每 sample × mode 的 caller baseline TP:FP ratio 是多少？哪些樣本有足夠統計力做 FP 改進檢定？
- **回答**：v2 §11.1 - §11.2（見 `tpfp_baseline_ratio.tsv`）
- **結論**：HCC1395 TO (2.47:1, 11606 FP) + COLO829 pf (15.5:1, 2273 FP) 是唯二有足夠 FP 的測試場

### Q2：scheme 切分有效性
LOH × AF × CN 定義的 S1-S7 scheme 是否切分出 TP:FP 兩極？
- **回答**：v2 §3 + §11（見 `tpfp_stratified_filter_schemes_TO.tsv`）
- **結論**：S3 (95.5%)、S5 (91.8%)、S1 (90.1%) 為高 purity 白名單；S4 (71.1%) 為無辨別力桶

### Q3：跨樣本泛化性
這些 scheme 在所有樣本都成立嗎？FP 是否集中在某樣本？
- **回答**：v2 §12（見 `tpfp_per_sample_scheme_full.tsv`）
- **結論**：paired_full S3/S5 TP rate ≥97% 全樣本一致；但 S1 在 COLO829 fold=0.59× 異常，LOH_Strong 跨樣本定義不同

### Q4：特徵是否可再細分
同一 scheme 內 TP vs FP 的特徵分佈是否不同？若有差，可否加特徵 threshold 提高 purity？
- **回答**：v2 §14（見 `tpfp_feature_diffs.tsv`、`tpfp_subschemes.tsv`）
- **結論**：AlleleDelta δ=-0.49 為最強特徵差，但單特徵 refinement 最大 fold=1.17×（S4d）→ **在 scheme 內飽和**

### Q5：綜觀聯合空間是否飽和
若放棄「在 scheme 內加單特徵」框架，綜觀 5D 聯合 bin，是否存在超越 S3/S5 的 operating point？
- **回答**：v2 §17（見 `tpfp_5d_pareto_HCC1395_TO.tsv`、`tpfp_5d_cumulative_envelope_HCC1395_TO.tsv`）
- **結論**：**未飽和**。17-cell envelope (96.1%, 3.7%) Pareto-dominates S3；28-cell envelope (90.4%, 7.4%) Pareto-dominates S5

---

## 1.3 假說與可否證條件

### H1：LOH × AF × CN 能切分 TP/FP

- **內容**：biology-informed 切分比隨機切分更能區分 TP 與 FP
- **證據層級**：Tier 1（P1）
- **支持**：v2 §3 scheme S3/S5 TP rate 遠高於 baseline；§12 跨樣本穩健
- **否證條件**：若 HCC1395 TO 的 S3/S5 TP rate ≤ 80%（僅比 baseline 71% 高 9%）→ H1 被弱化
- **目前狀態**：✅ 強支持

### H2：在預定 scheme 內加單特徵 threshold 會繼續提高 purity

- **內容**：若 TP vs FP 在 S4 內有特徵差（如 AlleleDelta），加 threshold 可切出 sub-high-TP
- **證據層級**：Tier 2（P4）
- **否證條件**：若所有 sub-scheme fold≤1.2× → H2 被駁斥
- **目前狀態**：❌ **駁斥**。最佳 S4d fold=1.17× → 單特徵 refinement 在此框架飽和

### H3：聯合 5D cube 存在超越單 scheme 的 Pareto operating point

- **內容**：若 cube 的聯合 bin 能捕捉 biology 同時表達的多個維度，top-k cells 可 Pareto-dominate S3/S5
- **證據層級**：Tier 1（P1）
- **否證條件**：若 cumulative envelope 緊貼 S3/S5（以 ±1% purity 帶）→ H3 被駁斥
- **目前狀態**：✅ **強支持**。17-cell envelope 在 HCC1395 TO 下 purity 高 + recall 高（96.1%, 3.7%）同時勝過 S3 (95.5%, 1.3%)

### H4：top-k cells 跨樣本泛化（**尚未驗證**）

- **內容**：在 HCC1395 TO 找到的 top-17 cells 在 COLO829 pf 也呈現高 purity
- **證據層級**：Tier 1（P1）
- **否證條件**：若 HCC1395 TO 的 top cells 在 COLO829 TP rate ≤ 85%（劣於 COLO829 baseline 94%）→ H4 駁斥
- **目前狀態**：🔴 **未驗證**。初步觀察（v2 §17.3）顯示 COLO829 top cells 座標異於 HCC1395 TO（NG=1 vs NG=2-3；Intermediate vs Near-half）→ 泛化性疑慮
- **本資料夾要做**：`obs08_overlap_venn.py` + 後續 LOSO 檢驗

---

## 1.4 測試框架矩陣

| 測試場 | 目的 | 統計力 | 限制 |
|--------|------|--------|------|
| **HCC1395 TO mode** | 主 FP 辨別力檢定 | 極高（11,606 FP, 71.1% baseline）| 單樣本，無跨樣本驗證力 |
| **COLO829 paired_full** | 次要 FP 檢定（high baseline）| 高（2,273 FP, 94.1% baseline）| baseline 已高，fold 上限 1.3× |
| HCC1395 paired_full | TP purity sanity check | 中（627 FP, 97.9% baseline）| FP 稀少，scheme 改進無意義 |
| HCC1937/H2009/H1437 pf | Cross-sample sanity | 低-極低（FP<200）| 無辨別力量測力 |

---

## 1.5 非本研究範疇

- ❌ 跑其他 6 樣本 TO mode（需 C++ pipeline rerun）
- ❌ 多變量 ML 模型（logistic / RF / XGBoost）— 單特徵 threshold 先用
- ❌ 真 F1 計算（需 SEQC2 FN 資料，不在 master TSV）
- ❌ Bootstrap CI 全面重算（Wilson CI 已足夠當前 n）
- ❌ flag=on 重跑驗證 NG 貢獻（屬獨立 HP-only 研究）

---

## 1.6 成功標準

✅ 本研究已達成：
- 回答 Q1-Q5 五層問題
- 建立 5D Pareto envelope 框架
- 提供 17-cell/28-cell empirical white-list candidates

🟡 未完待續：
- H4（top-k 跨樣本泛化）需 LOSO 驗證
- 論文定位（read-level epigenetic characterization vs variant filter）需進一步對焦
