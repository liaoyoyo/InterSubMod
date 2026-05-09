<!--
build_date: 2026-05-10
agent: V6-C Phase A theory derivation
status: validated
report_class: theory-analysis
audience: PI / lab member / 自己未來
parent_plan: InterSubMod/research/paired_priority_bug_audit/03_V6C_HPFineNGroups_remand_plan.md
inputs:
  - InterSubMod/research/F_hpfinengroups_deepening/data/step1_per_sample.tsv
  - InterSubMod/research/F_hpfinengroups_deepening/data/step4a_cohens_d_af04.tsv
  - InterSubMod/research/tpfp_loh_af_kde_discrimination/data/obs18_NG2_composition_by_sample.tsv
  - InterSubMod/research/v5_provenance_followup/T1_2_read_level_audit/T1_2_F1_genome_wide_audit.md
  - InterSubMod/src/core/ReadParser.cpp:121-153 + LabelTest.cpp:230-282
outputs:
  - 本檔（Phase A 結論）
verdict: Phase A 揭露關鍵機制 — flag=on 下 NGroups ≥3=0 是 bucket schema collapse 數學必然（demote hp=11/21/33 為 unphased → somatic-tag buckets HP1-1/HP2-1/HP3 全消失），**非 marker artifact 證據**。Phase A 無法獨立判斷 marker 真實性；必須升級至 Phase B region-level TP rate cross-flag 驗證
last_verified: 2026-05-10
decision: Phase A 結束（無法定論）→ 直接進 Phase B chr19 子集 ISM 重跑 flag=on/off 對比，看 region-level TP rate 是否保留
-->

# V6-C Phase A 理論推導結論 — 升級至 Phase B 必要

## 0. TL;DR

> Phase A 從既有數據推導發現：**flag=on 下 NGroups ≥3 完全消失（22,732 → 0 TP regions）**不是 marker artifact 證據，而是 **bucket schema collapse 的數學必然**（demote hp=11/21/33 為 unphased → somatic-tag buckets 全消失，NG count 機制上必降）。原計劃 Phase A 想用 inflation rate 估計直接判斷 marker 真實性，但發現此判斷需要 **region-level cross-flag TP rate tracking**（即同一 region 在 flag=off 與 flag=on 下 TP rate 比較），無法從現有 aggregated data 推導。**升級至 Phase B 必要。**

---

## 1. Phase A 動機與目標

依 V6-C plan §2.A：
- 用既有 5/9 Step D + obs18 + F pilot 數據做機制推導
- 估計 priority bug 對 HPFineNGroups marker 影響度上下界
- Decision gate：< 10% 影響 → 結束 / 10-50% → Phase B / > 50% → Phase C

## 2. 機制推導

### 2.1 HPFineNGroups 計算規則

從 `InterSubMod/src/core/LabelTest.cpp:265-302` `hp_to_fine_labels()`：

每個 region 統計 4-bucket occupancy（read tag 各類別的存在與否）：

| Bucket | 對應 BAM HP tag | 對應 ISM `hp_tag_raw` | flag=on 後 |
|---|---|---|---|
| HP1（germline）| HP:i:1 | "1" | **保留** |
| HP2（germline）| HP:i:2 | "2" | **保留** |
| HP1-1（somatic on HP1）| HP:i:11 | "1-1" | **demote 為 "0"** → 消失 |
| HP2-1（somatic on HP2）| HP:i:21 | "2-1" | **demote 為 "0"** → 消失 |
| HP3（somatic ambiguous）| HP:i:33 | "3" | **demote 為 "0"** → 消失 |

NG = number of distinct buckets with ≥1 read。

### 2.2 flag=on 下 NG count 的數學上界

| flag=off NG | bucket pattern (例) | flag=on bucket（demote hp=11/21/33）| flag=on NG |
|---|---|---|---|
| 4 | {HP1, HP1-1, HP2, HP2-1} | {HP1, HP2} + 一些 "0" | **2** |
| 3 | {HP1, HP1-1, HP2} | {HP1, HP2} | **2** |
| 3 | {HP1-1, HP2, HP2-1} | {HP2} + 多 "0" | **1** |
| 2 same_HP1 | {HP1, HP1-1} | {HP1} + 多 "0" | **1** |
| 2 same_HP2 | {HP2, HP2-1} | {HP2} + 多 "0" | **1** |
| 2 cross_het | {HP1, HP2-1} | {HP1} + 多 "0" | **1** |
| 1 germline-only | {HP1} | {HP1} | **1** |
| 1 somatic-only | {HP1-1} | 全 "0" | **0** |

→ **flag=on 下 NG max = 2**（只有 HP1, HP2 兩 buckets）。**NGroups ≥3 完全消失是數學必然**。

### 2.3 Phase 1 audit 的「NGroups ≥3=0」應如何解讀

| 解讀 | 是否成立 |
|---|---|
| flag=on 下 NG≥3=0 是 priority bug artifact 證據 | **不成立** — 此為 schema collapse 必然結果，與 priority bug 無關 |
| flag=on 下 NG≥3=0 等於 marker 信號消失 | **不成立** — bucket schema 改變但 region 物理特徵（reads / 位置 / methylation）不變 |
| flag=on 下 NGroups distribution 改變反映 marker 對 somatic-tag 依賴 | **部分成立** — 但無法區分「真實生物學依賴」 vs 「priority bug 副產物依賴」 |

→ Phase 1 audit 結論「flag=on NGroups ≥3=0」**不能直接證明 marker 是 artifact**。

## 3. 真正需要驗證的問題（Phase B/C 必要）

### 3.1 重定義驗證問題

V6-C 真正要回答：

> 對 flag=off 下被 marker filter 抓到的 regions（如 NG=4+AF<0.4+NR≥80 NonLOH 共 14,197 regions），這些 regions 在 **flag=on 重定義後**（NG count 機制上會變）的 **per-region TP rate 是否保留 ~93%**？

而**不是**：

> flag=on 下 NG distribution 是否改變（→ 必然改變，schema collapse）

### 3.2 Cross-flag region-level tracking 需求

| 量化任務 | 需要的資料 |
|---|---|
| Region 在 flag=off 認定為「NG=4+AF<0.4+NR≥80 NonLOH」的 set A | F pilot canonical filter 已有（14,197 regions）|
| Region 在 flag=on 重新計算後的 NG / TP / FP | **需新跑 ISM with flag=on**（per-region significance_summary） |
| set A 在 flag=on 下的 per-region TP rate | 上面兩者 join |

### 3.3 為什麼 Phase A 既有資料無法答

| 既有資料 | 限制 |
|---|---|
| F pilot per-sample TP rate（flag=off）| 沒 flag=on 對比 |
| obs18 NG=2 composition（flag=off）| 沒 flag=on 對比 |
| Phase 1 audit summary（aggregated AUC）| 沒 region-level cross-flag tracking |
| 5/9 Step D paired germline-absent | 是 read-level 不是 region-level |
| T1.2-F1 全基因組 vote dump | 是 read-level，無法直接對應到 region 的 TP rate |

→ **Phase A 無法獨立定論**，必須升級。

## 4. inflation rate 概略估計（lower bound 評估）

雖然不能定論，Phase A 可給 priority bug 對 marker 信號污染的 **lower bound** 估計：

### 4.1 read-level inflation

- 全基因組 V5 Layer 1.5 額外 tag = 560,881 reads（5/9 Step D）
- 全基因組 V3F/V5 tagged total = 18,895,432 reads
- V5 vs V3F germline-absent inflation rate = 560,881 / 18,895,432 = **2.97% reads**

### 4.2 region-level inflation 估計

若 F pilot canonical filter regions（NG=4+AF<0.4+NR≥80 NonLOH，14,197 個）內 reads 分布類似全基因組：
- 每 region 平均約 80+ NR
- 每 region 中 germline-absent reads 占比約 3%（同全基因組比例）
- 完全靠 germline-absent reads 才達到 NG=4 的 region 占 < 5% 估計

→ **Lower bound：~5% F pilot regions 可能受 priority bug inflation**（即 NG count 純靠 germline-absent reads）

### 4.3 Upper bound 限制

但 **bucket schema collapse 機制** 使得 flag=on 下這些 regions 統一變成 NG ≤ 2，即使原 region 的「真實 sub-clone signal」存在。所以「flag=on 下這些 regions 失去 marker 標籤」並非「marker 是 artifact」的證據。

→ Upper bound：~5% inflation 是合理估計。

## 5. Phase A 結論

| 項目 | 結論 |
|---|---|
| Phase 1 audit「flag=on NGroups ≥3=0」是 marker artifact 證據？ | **否** — schema collapse 數學必然 |
| Priority bug 對 marker 信號的 lower bound | ~5%（read-level inflation）|
| Phase A 既有資料能定論嗎？ | **不能** — 需 region-level cross-flag tracking |
| 是否進 Phase B？ | **是** — 必要，否則 V6-C 無法定論 |

## 6. Phase B 修訂 design

原 Phase B（chr19 子集 ISM 重跑 + NG distribution 對比）**需要修訂**為 region-level 對比：

### 6.1 修訂後 Phase B 步驟

B1. **chr19 子集 ISM 跑兩 mode**：
- flag=off：產 `output/v6c_phaseB_off/significance_summary.csv`
- flag=on：產 `output/v6c_phaseB_on/significance_summary.csv`

B2. **Region-level join**：
- 用 region 識別碼（chr_pos_window）對齊兩 mode 結果
- 計算 set A = flag=off 內 NG=4+AF<0.4+NR≥80 NonLOH regions
- 對 set A 在 flag=on 下：(a) NG 重新計算結果 (b) TP / FP 標籤是否一致 (c) per-region TP rate

B3. **Decision gate**:

| flag=on 下 set A 的 TP rate | 結論 |
|---|---|
| ≥ 0.85（保留 marker level）| 真實生物學 marker；schema collapse 不 invalidate |
| 0.70-0.85 | 部分 schema 依賴；中間 |
| < 0.70 | marker 主要靠 schema artifact；需 Phase C 7 樣本確認 |

B4. **預期時間**：~5 hr（chr19 ISM 5 min × 2 + 分析 30 min + write-up 2 hr）

### 6.2 Phase C 仍按 plan 執行

若 Phase B 顯示中間 / 矛盾結果，Phase C 7 樣本全量重驗按 plan §2.C 不變。

## 7. 推薦 next step

**Decision**：直接進 Phase B 修訂版（chr19 子集 ISM × 2 flag region-level cross-tab）。

**啟動條件**：
- KDE-fixed binary 仍可用 ✅
- HCC1395 5kHz BAM + ClairS-TO VCF + REF 仍可用 ✅
- germline_hp_only flag CLI 介面 ✅（既有）

**估時**：~5 hr 完整 Phase B 含 write-up。

User 確認後啟動 Phase B 修訂版。

## 8. 對 V6 evaluation 的補充

V6 evaluation （`02_V6_proposal_evaluation.md`）的核心結論不變：
- V6 binary 不必要（ISM flag 已等效）✅
- 但其中「Phase 1 audit 已驗證 marker 對 filter 無 ≥+0.02 改善」的解讀需精確化：
  - **AUC 衰退**（HPFineNGroups -0.026, etc.）是 schema collapse 後 NG count 自然改變的結果
  - **不能直接等同於「marker 是 artifact」**
  - 真實 artifact rate 需 Phase B region-level cross-flag tracking 才能確認

→ 5/8 主報告 §8.6 的 V6-C 描述可加 caveat（可選），但本檔已記錄機制釐清。
