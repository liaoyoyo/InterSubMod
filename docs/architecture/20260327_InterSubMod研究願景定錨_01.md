<!--
建立時間: 2026-03-27 00:00
目標: InterSubMod 研究願景正式定錨文件，作為後續所有計劃的北極星
處理範圍: 五個研究目標、主軸優先順序、成功標準、LOH 定位
關聯檔案:
  - docs/CURRENT_FOCUS.md
  - docs/references/20260324_InterSubMod研究方法與相關文獻突破方向全域分析_01.md
  - docs/reports/validated/2026/03/20260327_LOH_round1_cross_sample_audit_01.md
-->

# InterSubMod 研究願景定錨文件

> **性質**：正式 north star 文件。後續所有計劃、實驗、報告都以此為參照基準。
> 除非有強證據或後續研究發現有無法解決的困難，造成需要修正與重新排序重要性，否則以此方向為目標。
> 若需修正，修改此文件並記錄原因。

---

## 一、研究定位

InterSubMod 是一個「**可解釋 epigenetic evidence 整合框架**」，不是另一個 somatic caller。

核心貢獻定位：

1. 提供 per-CpG 甲基化的多標籤關聯性分析，作為 caller 無法直接觀察的維度
2. 利用 read-level epigenetic context 補充 variant 的生物背景資訊
3. 特別針對 tumor-only 與低 purity 情境的困難問題，提供補充 evidence

**不是**：替換現有的 somatic variant caller（ClairS、DeepSomatic 等）。

---

## 二、五個研究目標

### 目標 1 — per-CpG 甲基位點多標籤關聯性評分

**核心問題**：某個 CpG 位點的甲基狀態，是否與 ALT/REF、HP1/HP2 等標籤具有統計關聯？

**輸入**：
- BAM MM/ML tags（甲基化狀態）
- phased VCF（allele 資訊）
- HP tags（haplotype 標籤）

**輸出**：
- 每個 CpG 位點的「關聯標籤集合」
- 每個標籤的「關聯強度分數」（例如：HP imbalance score、allele-specific methylation score）

**適用條件**：
- 需要足夠 read coverage 的位點（有效 read threshold 待定義）
- coverage 不足的位點，關聯強度分數要標記為「不可信」

**生物動機**：
- 甲基狀態受多種因素影響，不穩定
- 癌症甲基化存在 cis 效應，常導致抑癌基因表達下調
- 甲基演化可能與 SNV 演化分離，部分甲基異常是獨立發展的，需要細分與分析

---

### 目標 2 — 單 region 與跨 region 的 clone 結構分析

**核心問題**：

1. 單位點：某個 sSNV region 是否存在「clone 後再 clone」的第二層突變結構（sub-clone within a clone）？
2. 跨位點：多個 sSNV 是否同屬一個 clone 事件（共演化的 clone 同步性）？

**分析方式**：
- 單位點：分析 ALT read 中的甲基與 HP 分佈，找出 within-ALT 的異質性
- 跨位點：分析多個 sSNV region 的 read overlap 與甲基模式，推論共演化關係

**輸出**：
- subclone 結構圖
- 突變 co-occurrence 分析

**注意**：相關二代定序研究已有先例，但長讀取在 read-level 的解析力更高。推論需謹慎，需清楚說明可信度邊界。

---

### 目標 3 — 二次打擊與事件發生順序推論

**核心問題**：某個基因/區域是否存在兩次打擊（two-hit hypothesis）？發生順序是什麼？

**輸入**：
- 甲基演化差異（目標 1 的輸出）
- LOH、CNV、HP、SNV、INDEL 等標籤
- Haplotype 資訊（在 PS block 範圍內）

**分析方式**：
- 組合分析：哪些標籤的組合暗示先後順序
- 條件機率或模擬方法

**輸出**：
- 事件發生順序的推論結果（含可信度說明）
- 二次打擊候選位點清單

**前提條件**：此目標依賴目標 1 和目標 2 的基礎分析，不能直接跳入。

---

### 目標 4 — Tumor-Only 的 normal 資訊補強

**核心問題**：在無配對 normal 的情況下，如何利用甲基與 HP 資訊推估 normal 背景？

**挑戰**：
- TO 的 normal 資料（若有）來自血液，含有不同身體組織的甲基差異
- Purity 問題：tumor 與 normal 混合的比例影響甲基與 allele 訊號的解讀

**輸入補充**：
- 是否提供 `normal.bam`
- longphase-s 或 longphase-to 的 purity 估計輸出

**分析方向**：
- 利用甲基與 germline SNP 的關聯性，推估 normal tissue methylation 背景
- 配合 purity 資訊做 deconvolution

**目標**：提升 TO 在低 purity 情境下的 F1，特別是 TO 的 FP 過濾

**說明**：TO 是獨立主線課題（不是 paired 的簡化版），但可以用 paired 作對照與驗證。

---

### 目標 5 — 整合 evidence panel 提升最終 somatic 判斷的 F1

**核心問題**：整合目標 1–4 的資訊，能否提升 caller 輸出的最終 F1？

**策略**：

1. 分析各種特徵在 TP/FP 之間的分佈差異，且確認是否跨樣本一致
2. 對 TP/FP 分佈無清楚差異的特徵，嘗試加入分類標籤（例如 LOH、HP 分層）後再分析
3. 利用 germline SNP/PON 資訊做連鎖判斷（經過的 read 是否更像 germline）

**輸出**：
- 可保留更多 TP 的 rescue 規則
- 可過濾 FP 的 filter 規則
- 最終 F1 benchmark 與 baseline 比較

**注意**：不替換 caller 輸出，而是在 caller 輸出的基礎上做 evidence panel annotation 與過濾。

---

## 三、研究主軸優先順序

```
Phase 1（當前活躍）：方向 E — ML read classification
                      within-tumor alt-support read 的分類模型
                      目標 1 + 目標 5 的第一步落地

Phase 2：方向 A + D  — normal methylation reference + CN/Purity-aware
                      目標 4 的核心方向

Phase 3：方向 B      — gene-level / mechanism-level evidence integration
                      目標 3 的基礎

Phase 4：方向 C      — CpG 功能分層與智慧選點
                      目標 1 的精細化
```

**LOH evidence panel** 屬於跨 Phase 的基礎觀察層：
- 在 Phase 1 的 feature 工程中就需要
- 先做好 region-level LOH 診斷，再往 per-CpG 展開

**當前 Phase 1A 進度**（截至 2026-03-25）：
- 已建立 per-read training/export layer
- exporter 已納入 Phase 1A/1B label schema
- 141,014 個 regions 已納入 manifest（4 個 baseline datasets）
- Phase 1B（is_tumor 分類）目前受限於 `normal_reader` 尚未接入 `process_single_region`

---

## 四、成功標準（拆分定義）

成功標準拆為兩類，**兩者都是合法的成功**，不需要都達成才算有貢獻：

| 成功類型 | 定義 | 驗證方式 |
|---------|------|---------|
| **可解釋性成功** | 能清楚說明某個 CpG/region 的甲基模式與 HP/SNV 的關聯理由，提供可追溯的 evidence | case study 驗證、evidence panel 可追溯 |
| **F1 提升成功** | paired 或 TO 的 precision/recall 有可量化的改善 | benchmark 與 baseline 對比，需達統計顯著 |

**不能混用的陷阱**：
- 「evidence panel 解釋了很多位點」≠「F1 提升」
- 「某個 feature 在 TP/FP 有差異」≠「用它來過濾一定能提升 F1」

---

## 五、LOH 在研究願景中的定位

### 層次定義

LOH/HP imbalance 屬於 **evidence panel 第一層（基礎觀察），非主判斷器**。

### 在 paired 端的定位

- `LOH-like` 具有輕度 FP enrichment（整體 1.194×）
- 但有顯著 sample heterogeneity：HCC1954 (3.185×) vs HCC1395_HKU_5kHz (1.016×)
- 最合理用法：**sample-aware risk feature**，搭配 support quality 分層後才能作為 filter candidate

### 在 TO 端的定位

- `LOH-like` 呈 **TP 富集**（overall ~~0.912×~~ → **0.805×**）（[修正 2026-03-30] HP integer tag bug 修正後，TO LOH-like 在 TP 更常見，與 paired 方向相反。詳見 `20260330_TO_LOH_enrichment_post_hp_fix_01.md`）
- TO LOH-like 不是 FP 的 discriminative signal，而是 TP enrichment marker
- 必須配合 `effective_hp_reads`、`hp0_ratio`、`PS availability` 才能解讀

### 關鍵限制：PS block 邊界

- LOH 只在同一個 PS block 範圍內有效
- 不同 PS block 的 HP1 / HP2 **無法**互相對應到父本或母本
- 不能跨 block 串聯成長程 haplotype 敘事

### 里程碑

| Round | 目標 | 狀態 |
|-------|------|------|
| Round 1 | 確立第一層診斷底圖（LOH-like 分佈、TP/FP 差異、same-locus compare） | ✅ 完成 |
| Round 2 | 確立「有效 LOH evidence」的品質條件（support 分層 + HP0 來源分析） | ⏳ 待執行 |
| Round 3 | LOH × methylation 聯合分析（在 support quality 分層確立後） | 🔜 待規劃 |

---

## 六、paired 與 TO 的互補策略

**TO 是獨立主線**：
- TO 的 FP 多、課題性質不同（沒有配對 normal）
- 不是 paired 的簡化版，需要獨立的分析框架

**paired 是對照與驗證來源**：
- 提供較乾淨的 reference，幫助定義「可解釋 evidence panel」
- 提供 cross-mode 一致性驗證

**互補方向（不對稱）**：
- `paired → TO`：paired 先定義哪些 evidence 在高品質條件下是可信的，再移植到 TO
- `TO → paired`：TO 的大量 FP 提供壓力測試，驗證 evidence 是否真能處理最難的場景

**同位點比較（Round 1 發現）**：
- `both_tp`：62.44%
- `TO-only FP`：27.59%（主要不一致來源）
- `paired_only FP`：0.42%（少）
- `both_fp`：0.33%（極少）

---

## 七、分析進展路徑

```
region-first → per-CpG
     ↓
「可解釋 evidence panel」→「事件發生順序 inference」
     ↓
tumor-only（獨立主線）   paired（對照驗證）
         ↕ 互相輔助（不對稱）
```

每一步的前置條件：
1. region-level LOH 診斷品質確認 → 才能往 per-CpG 展開
2. 可解釋 evidence panel 建立 → 才能進入事件順序推論
3. paired 的 evidence panel → 才能移植到 TO 主線
