<!--
建立時間: 2026-04-26
狀態: validated
類型: retraction declaration
撤回對象: Thread B (LOH × AF × CN cross-sample whitelist filter)
觸發證據: docs/experiments/in_progress/2026/04/20260424_X6_Caller_AF_S3S5_CrossSample_01.md
受眾: PI、外部讀者、未來 AI agent（避免重訪已撤回方向）
-->

# Thread B 撤回宣告 — LOH × AF × CN 跨樣本 Whitelist Filter

> **狀態**：⚠ **RETRACTED 2026-04-26**（filter 用途）／✅ **保留**（HCC1395 single-sample case study + per-sample characterization）
> **撤回理由**：X6 caller_af 重 merge + KDE-corrected master 跨樣本驗證證實 S3/S5 跨樣本不穩定，原 v2 HCC1395 TO S3 95.5% 為 stale-binary CN_tier artifact。
> **論文主軸切換**：自本日起，TO 層論文主軸鎖定 **Thread D LOH-constrained phasing signatures**（佔位連結：[InterSubMod/docs/reports/validated/2026/04/20260426_Thread_D_LOH_constrained_phasing_main_axis_01.md](20260426_Thread_D_LOH_constrained_phasing_main_axis_01.md)）。

---

## 1. 撤回摘要

**原宣稱**：Thread B v2 報告（`InterSubMod/docs/experiments/in_progress/2026/04/20260422_LOH_AF_KDE_TPFP_Discrimination_02.md`）依 HCC1395 TO 單樣本 40,115 sites 推論 **S3 Diploid Het (LOH=None ∩ Near-half AF ∩ CN∈T1/T2) TP rate 95.5%**、**S5 combo TP 91.8% / FP reduction 99.37%**，並計畫推廣為 TO 全樣本 **biology-informed stratified whitelist filter**（白名單型 cross-sample filter）。

**為何撤回**：X6（2026-04-24）以 caller AF 重 merge + KDE-corrected master 在 6 TO 樣本（HCC1395 / HCC1395_DORADO / H1437 / H2009 / HCC1937 / HCC1954）跨樣本重算 S3/S5 scheme，得到決定性反證——**S3 TP≥0.85 達成比例 1/6**（僅 H2009，且為 baseline 飽和 artifact）、**Wilcoxon S3 > baseline one-sided W=0、p=1**（S3 系統性低於 baseline）、原 v2 HCC1395 TO S3 95.5% 在 KDE-corrected master 完全無法複現（n 從 380 → 2,200，TP 從 95.5% → 58.3%）。

**撤回範圍**：跨樣本 LOH × AF × CN whitelist filter 用途的所有宣稱、CL-S3-001 ⭐4 tier 升級、論文 whitelist 落地工具定位。**保留**：HCC1395 single-sample TO case study 描述、per-sample biology characterization、Thread B 的 LOH_Subtype × AF_class × CN-tier 切分框架本身（作 characterization annotation）。

---

## 2. X6 決定性證據

### 2.1 跨樣本 S3 Diploid Het scheme TP rate（X6 caller_af verified）

| Sample | S3 n | S3 TP rate | vs baseline |
|--------|-----:|-----------:|-------------|
| HCC1395 | 2,200 | **0.583** | ↓（baseline ~0.71） |
| HCC1395_DORADO | 1,805 | **0.597** | ↓ |
| H1437 | 6,588 | **0.738** | flat |
| H2009 | 31,436 | **0.903** | ≈baseline 0.93（飽和） |
| HCC1937 | 1,920 | **0.358** | ↓↓ |
| HCC1954 | 10,729 | **0.129** | ↓↓↓ |

- **S3 TP≥0.85 & n≥20 達成比例**：**1/6**（僅 H2009，且 H2009 baseline 已 0.93 → S3 屬飽和 artifact，非 S3 切分力證據）
- **Wilcoxon one-sided S3 > baseline**：**W=0.0, p=1**（S3 系統性**低於** baseline）
- **S5 combo TP≥0.85 & n≥50 達成比例**：1/6（同 S3，僅 H2009 飽和）

### 2.2 原 v2 HCC1395 TO S3 95.5% 的 stale-binary CN_tier 歸因

| 指標 | 原 v2（HCC1395 TO） | X6 post-KDE（同樣本） |
|------|-------------------:|--------------------:|
| S3 n | 380 | **2,200**（5.8×） |
| S3 TP rate | **95.5%** | **58.3%**（−37 pp） |
| S3 FP reduction | 99.85% | 92.1% |
| S5 n | 886 | 10,652 |
| S5 TP rate | 91.8% | 72.6% |
| S5 FP reduction | 99.37% | 74.9% |

**根因**：原 v2 master 由 KDE commit `8d0a0c8`（2026-04-13）**之前**的 stale binary 產出（`expected_coverage` 共用 75.0 hardcoded default），導致 Coverage_Multiple 系統偏低 → CN_tier 分類落入錯誤 bucket → S3 cell 範圍人工縮小到僅含高純度 canonical het 子集 → TP rate 人工偏高至 95.5%。KDE-corrected 後 CovM median 0.880 → 1.245（×1.415），S3 cell 重新涵蓋更大範圍的 None × Near-half × CN∈T1/T2 sites，TP rate 收斂到接近 baseline。

### 2.3 數據來源

- **撤回觸發實驗**：`InterSubMod/docs/experiments/in_progress/2026/04/20260424_X6_Caller_AF_S3S5_CrossSample_01.md`
- **原 v2 報告（被撤回）**：`InterSubMod/docs/experiments/in_progress/2026/04/20260422_LOH_AF_KDE_TPFP_Discrimination_02.md`
- **caller AF merge script**：`InterSubMod/scripts/analysis/20260424_X6_merge_caller_af_S3S5.py`
- **JSON artifact**：`InterSubMod/research/tpfp_loh_af_kde_discrimination/data/X6_caller_AF_S3S5.json`

---

## 3. 降級範圍對照表

| 宣稱／用途 | 撤回前 | 撤回後（2026-04-26） | 證據 |
|-----------|-------|------------------|------|
| **跨樣本 whitelist filter（論文落地工具）** | ⭐4 POSITIVE，計畫推廣 6 TO 樣本 | ❌ **NEGATIVE 撤回**（5/6 樣本 < baseline，Wilcoxon p=1） | X6 §2.2 |
| **HCC1395 TO S3 95.5% TP rate** | 主結果 ⭐ | ❌ **stale-binary artifact**（KDE-corrected 為 58.3%） | X6 §2.4 |
| **S5 combo (S1∨S2∨S3 \ S4) FP reduction 99.37%** | TP recall 2.85% 高 precision 白名單 | ❌ KDE-corrected FP reduction 74.9%（−24.5pp） | X6 §2.3 |
| **CL-S3-001 ⭐4 tier**（週報 0423） | ⭐4 升級 | ⬇ **降級為 ⭐2 characterization-only** | 本文件 |
| **HCC1395 single-sample TO case study** | — | ✅ **保留**（已知 caller_af + LOH bed + KDE 條件下的單樣本切分樣態） | X6 §3.1 |
| **per-sample biology characterization** | — | ✅ **保留**（LOH_Subtype × AF_class × cn_tier_F 切分框架本身） | 本文件 §5 |
| **S4 ambiguous bucket（LOH=None ∩ Extreme AF）** | 75% TP + 76% FP 無辨別力 | ✅ 保留結論（X6 多樣本同樣表現 baseline） | X6 §2.4 |

---

## 4. 受影響舊報告清單與修正建議

| 檔案 | 影響段落 | 修正動作 |
|------|---------|---------|
| `InterSubMod/docs/reports/validated/2026/03/20260329_LOH_round1_cross_sample_audit_v2_01.md` | 全文（v2 提出 LOH × AF cross-sample 審計） | 已加 retraction banner（單樣本 case study 保留） |
| `InterSubMod/docs/reports/validated/2026/04/20260423_研究週報_20260416_20260423_NG2_LOH_constrained_phasing與TO_pivot_01.md` | Thread B 段、S3/S5 證據卡（lines 39, 72-73, 84, 365-488） | 已於 S3/S5 段加 caveat（caveat 連結本撤回宣告） |
| `InterSubMod/docs/experiments/in_progress/2026/04/20260422_LOH_AF_KDE_TPFP_Discrimination_02.md` | 全文（v2 報告本體） | 已於 INDEX.md 入口加 `[RETRACTED 2026-04-26]` 標記 |
| `InterSubMod/docs/experiments/INDEX.md` | line 236 (LOH × AF × CN) / line 242 (週報 0423) | 已加 `[RETRACTED 2026-04-26]` 前綴 |
| `InterSubMod/docs/CURRENT_FOCUS.md` | 主軸描述、P2 行動 | 已切換為 Thread D 主軸 |
| `InterSubMod/research/feature_layered_observation/00_main_observation.md` | §4「與 Thread B S1-S7 scheme 的 feature coverage 對照」 | 加 retraction caveat（S1-S7 scheme 仍可作 characterization 框架，不再宣稱 cross-sample filter） |

---

## 5. 未撤回項目（明確區分）

### 5.1 ❌ 撤回（filter 用途 / cross-sample whitelist）

- S3 Diploid Het 為 TO 全樣本 white-list filter 的宣稱
- S5 combo 為高 precision white-list 的宣稱
- 「biology-informed stratified filter framework POSITIVE」（cross-sample filter 含義）
- CL-S3-001 ⭐4 tier 升級
- 推廣 Thread B 為論文落地工具的計畫
- HCC1395 TO S3 = 95.5% TP rate 的數值宣稱（單樣本亦撤回此具體數值，因屬 stale-binary）

### 5.2 ✅ 保留（characterization / single-sample case study）

- LOH_Subtype × AF_class × cn_tier_F **切分框架本身**（作 per-variant annotation，不作 filter）
- **HCC1395 TO single-sample case study**：在 KDE-corrected master 下 S3 cell n=2,200 / TP=58.3% / FP reduction 92.1% 的單樣本切分描述
- **S4 ambiguous bucket** 結論：LOH=None ∩ Extreme AF 在多樣本上一致表現為「baseline 級無辨別力 bucket」（X6 跨樣本 reconfirmed）
- **LOH 任何註解 × Extreme AF TP 88-96%** 觀察（per-sample characterization，未推廣為 filter）
- **NG≥3 邊際貢獻 <1pp** 觀察（與 Thread D 主軸一致，NG 由 phasing 機制吸收）
- Thread A（CN KDE 校準方法學）獨立保留，**不受本撤回影響**
- Thread C（`--germline-hp-only` flag）獨立保留，**不受本撤回影響**
- Thread D（LOH-constrained phasing）獨立保留並升級為論文主軸

---

## 6. 連結到 Thread D 主軸報告（佔位）

自 2026-04-26 起，TO 層論文主軸正式鎖定為：

> **"LOH-constrained phasing signatures distinguish somatic from germline-like variants in tumor-only sequencing"**

主軸報告（佔位）：[InterSubMod/docs/reports/validated/2026/04/20260426_Thread_D_LOH_constrained_phasing_main_axis_01.md](20260426_Thread_D_LOH_constrained_phasing_main_axis_01.md)

關鍵證據：obs18 跨 6 TO 樣本 Inner × NG=2 ≥93% same-haplotype（6/6 一致，median 97%）、Inner same-hap vs Outer cross-het TP gap median +0.37（6/6 正向）。詳見 `InterSubMod/research/tpfp_loh_af_kde_discrimination/09_TO_sample_af_lohside_ng.md` 與 obs18 stacked composition data。

---

## 7. 後續行動

- [x] 撤回宣告本檔案（2026-04-26 完成）
- [ ] Thread D 主軸 validated 報告完稿（P0-1）
- [ ] HCC1954 standalone case panel — 於 Thread D 內處理 outlier（P0-3）
- [ ] AutoResearch direction 重寫，移除 Thread B whitelist 推廣（P2-3）
- [ ] Wakhan / SAVANA external CN pilot 設計（P2-2，作 Thread B 替代方案，非繼承 whitelist 宣稱）
