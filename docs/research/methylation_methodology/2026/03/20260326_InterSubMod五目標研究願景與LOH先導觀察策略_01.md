<!--
建立時間: 2026-03-26 18:40
目標: 定義 InterSubMod 的五個最終研究願景，對齊現行 phase roadmap，並評估 LOH / HP / PS 是否適合作為第一步先導觀察
處理範圍:
  - 五個最終研究願景的正式定義
  - 願景與 Phase 1~4 的對應關係
  - paired / tumor-only 的分線與互補
  - LOH / HP family / HP0 / HP3 / PS block 的分析邊界
  - 可快速驗證與小步試錯的研究路線
關聯檔案:
  - /big7_disk/liaoyoyo2001/InterSubMod/docs/CURRENT_FOCUS.md
  - /big7_disk/liaoyoyo2001/InterSubMod/docs/references/20260324_InterSubMod研究方法與相關文獻突破方向全域分析_01.md
  - /big7_disk/liaoyoyo2001/InterSubMod/docs/methodology/20260324_方法學審查全域結論報告_01.md
  - /big7_disk/liaoyoyo2001/InterSubMod/docs/research/methylation_methodology/2026/03/20260307_5kHz主實驗與方法學驗證藍圖_01.md
  - /big7_disk/liaoyoyo2001/InterSubMod/docs/plans/2026/03/20260307_純樣本甲基研究執行計畫_01.md
  - /big7_disk/liaoyoyo2001/InterSubMod/docs/references/manual/20260316_HP_LOH術語校正手冊_01.md
  - /big7_disk/liaoyoyo2001/InterSubMod/docs/methodology/20260323_LOH_HP偏移TP_FP區分能力分析_01.md
  - /big7_disk/liaoyoyo2001/InterSubMod/docs/reports/validated/2026/03/20260322_cross_sample_TO_ISM_gradient_analysis_01.md
-->

# InterSubMod 五目標研究願景與 LOH 先導觀察策略

## 1. 一句總結

InterSubMod 後續應明確分成兩層來思考：

1. **五個目標**是最終研究願景，不等於近期執行順序。
2. **Phase 1~4** 是為了逐步實現這五個願景的技術路線。

在這個前提下，近期研究應採：

1. `region-first`，之後再往 `per-CpG` 展開。
2. 先建立**可解釋的 evidence panel**，再追求**事件發生順序 inference**。
3. `paired` 與 `tumor-only` 視為**不同課題**，但應保留方法與證據互相借鏡的設計。
4. `F1 提升` 與 `可解釋性提升` 視為兩條並列成功標準，不互相取代。

---

## 2. 為什麼需要這份文件

`2026-03-24` 之後的主線文件，已經把近期技術優先序收斂為：

1. `Phase 1`: `ML read classification`
2. `Phase 2`: `normal methylation reference + CN/Purity-aware`
3. `Phase 3`: `gene-level evidence integration`
4. `Phase 4`: `CpG functional stratification`

這套 phase 規劃解決的是「**先做什麼技術補強**」，但尚未完整回答「**InterSubMod 最終到底想成為什麼**」。

因此本文件的角色不是取代現有 phase roadmap，而是把 InterSubMod 的**最終願景**獨立定義清楚，讓後續的實驗、文件、簡報與證據鏈都有同一個高層主軸。

若要把這份高層願景落到可執行的 LOH round，請直接搭配：

1. [LOH 盤點執行規格](/big7_disk/liaoyoyo2001/InterSubMod/docs/plans/2026/03/20260326_LOH盤點執行規格_01.md)

該文件負責固定：

1. summary table
2. decision ledger
3. figure naming
4. case panel
5. open questions 與追溯欄位

---

## 3. 五個最終研究願景

### 3.1 願景 1：建立 region-first 的甲基關聯標籤與強度輸出

**定義**

對每個以 somatic SNV 為錨點的 region，分析該區域內 read-level 甲基模式是否與：

1. `ALT / REF`
2. `HP family`
3. `Tumor / Normal`
4. 後續新增的其他結構性標籤

存在穩定且可量化的關聯，並輸出：

1. 關聯標籤
2. 關聯強度
3. 適用邊界
4. 可信度或 QC 註記

**近期定義**

1. 先以 `region` 作為主輸出單位。
2. `per-CpG` 展開屬於後續精化層，不作為第一輪必要條件。

**這個願景的價值**

1. 提供後續 second-hit、clone、tumor-only、F1 改善的共同底層證據。
2. 讓「哪些甲基位點值得信、哪些只是噪音或混雜」可以先被形式化。

### 3.2 願景 2：建立 clone / subclone 的 region 與多位點證據面板

**定義**

不只回答某個位點有沒有甲基差異，而是回答：

1. 此差異是否呈現 clone / subclone 結構
2. 是單一 region 的局部現象，還是多位點可重複的結構
3. 是否存在「clone 後又 clone」或不同 locus 上的共現模式

**近期定義**

1. 先建立可解釋的 `region-level evidence panel`
2. 暫不直接宣稱完成真正的跨位點 clone tree 或演化順序

### 3.3 願景 3：建立 second-hit / biallelic event 的解釋框架

**定義**

整合：

1. SNV / INDEL
2. LOH / CNV
3. promoter methylation
4. haplotype evidence
5. 多個 region 的甲基關聯結果

去回答：

1. 是否存在 second hit
2. 哪一種 second-hit archetype 最合理
3. 哪些 case 已足夠形成機制層的 evidence panel

**近期定義**

先建立**可解釋的 gene-level / locus-level evidence panel**，再考慮事件先後順序 inference。

### 3.4 願景 4：建立 tumor-only 與 purity-aware 的獨立研究主線

**定義**

`tumor-only` 不只是 paired 的簡化版，而是獨立課題：

1. 無 matched normal
2. PON 與 caller verdict 的角色不同
3. LOH / purity / normal contamination 的判讀方式不同
4. FP 型態更複雜，也更適合當方法壓力測試場景

**定位**

1. `tumor-only` 應被視為獨立主線。
2. 但 paired 的方法學概念，例如 normal baseline、purity-aware correction、LOH / CN 解釋框架，仍可作為 TO 的設計靈感。
3. 反過來，TO 中較多的 FP 與 hard cases，也可拿來反驗 paired 方法是否真的有效。

### 3.5 願景 5：整合多層資訊，提升 somatic 判讀能力

**定義**

把前四個願景產出的資訊整合後，用來回答：

1. 某個位點更像真實 somatic 還是 artifact / germline-like
2. 哪些 TP 可以被保留
3. 哪些 FP 可以被移除
4. 哪些 case 雖然不提升 F1，但可被更清楚地解釋與分層

**成功標準**

1. `F1 提升` 是一種成功
2. `可解釋性提升` 也是一種成功
3. 兩者應分開回報，不應互相覆蓋

---

## 4. 五個願景之間如何互相輔助

| 願景 | 直接幫助的後續願景 | 互補關係 |
|---|---|---|
| 願景 1 | 2, 3, 4, 5 | 若連 region-level 關聯都不穩，後續 clone、second hit、TO 與 F1 整合都無從建立 |
| 願景 2 | 3, 5 | clone / subclone 結構是 second-hit 與整合判讀的重要背景層 |
| 願景 3 | 5 | second-hit evidence panel 能把 region-level 訊號升級成生物學解釋 |
| 願景 4 | 1, 5 | TO 提供更高 FP 壓力，可用來驗證哪些方法真的有效 |
| 願景 5 | 1, 2, 3, 4 | 整合評分不是取代前面四個願景，而是把它們的輸出變成最終決策層 |

**核心原則**

1. 五個願景不是五條互斥路線。
2. 它們應被設計成共享資料層、共享 diagnostics、共享 hard cases。
3. 真正需要分開的是**問題定義**與**驗收標準**，不是完全拆掉彼此的資料與方法連結。

---

## 5. 五個願景與現行 Phase 路線的對應

| 現行 Phase | 技術主題 | 主要服務的願景 |
|---|---|---|
| `Phase 1` | `ML read classification` | 願景 1、4、5 |
| `Phase 2` | `normal baseline + CN/Purity-aware` | 願景 1、3、4、5 |
| `Phase 3` | `gene-level evidence integration` | 願景 2、3、5 |
| `Phase 4` | `CpG functional stratification` | 願景 1、2、3、5 |

**解讀**

1. 當前 roadmap 並沒有偏離你的五個願景。
2. 差別只在於：phase 是**技術工程順序**，五個願景是**最終研究定位**。

---

## 6. LOH 是否適合作為第一步觀察

## 結論

**適合，但只能當「第一步觀察 / 校正 / failure-mode mapping」，不適合直接當第一個主證明。**

### 6.1 為什麼適合先觀察

1. 現有輸出已經有 `HP_Ratio`、`Potential_LOH`、`LOH_Subtype`，不需要先大改 C++ 才能開始分析。
2. LOH 與 HP family imbalance 本身就是後續：
   - germline ASM
   - second hit
   - purity / CNV
   - TO/paired 差異
   的共同交界點。
3. 先看 LOH 能快速知道目前哪些訊號其實是在反映：
   - 真實結構性事件
   - 胚系 ASM
   - phase / HP 不完整
   - 或根本只是噪音

### 6.2 為什麼不適合當第一個主證明

根據既有分析：

1. `TO mode` 中 `Potential_LOH` 對 TP/FP 幾乎無區分力，`AUROC ≈ 0.50`。
2. `AF>=0.9` 的 LOH 區域在跨樣本 TO 分析中幾乎全部落在 `Noise`，ISM 對這類區域普遍失效。
3. `paired mode` 的 LOH / HP imbalance 雖有訊號，但很大一部分其實在解釋 germline ASM FP，不是直接證明 methylation 對所有問題都有效。
4. 目前 InterSubMod 尚未把 `PS / phase block` 正式輸出到 read-level 主表，因此 HP-driven 解釋仍有跨 block 誤判風險。

### 6.3 因此最合理的定位

LOH 應先作為：

1. **觀察軸**
2. **失效模式拆解軸**
3. **paired / TO 對照軸**
4. **Phase 2 前的前置 diagnostics 軸**

而不是第一個直接宣稱「可提升 F1」或「已完成 second-hit inference」的主證據。

---

## 7. LOH / HP / PS 的正式分析邊界

### 7.1 HP family 的建議定義

若要先做 LOH / HP imbalance 觀察，建議採：

1. `HP1 family = HP1 + HP1-1`
2. `HP2 family = HP2 + HP2-1`

這與目前核心程式一致，因為 `LabelTest` 的 merged HP test 就是這樣分組。

### 7.2 HP ratio 的建議定義

建議主指標：

```text
hp_ratio_core = HP1_family_reads / (HP1_family_reads + HP2_family_reads)
```

並同時保留：

1. `effective_hp_reads = HP1_family_reads + HP2_family_reads`
2. `hp_assign_rate = effective_hp_reads / total_reads`
3. `hp0_ratio = HP0_or_unphased_reads / total_reads`
4. `hp3_ratio = HP3_reads / total_reads`

### 7.3 HP0 與 HP3 的角色

`HP0` 與 `HP3` 不應直接併入 HP1 / HP2 分母。

原因：

1. `HP0` 更接近 phase incompleteness / unphased 背景
2. `HP3` 在 LongPhase-S 是 ALT-supporting 但 germline HP 無法唯一解析
3. 兩者若硬併入 HP1 / HP2，會把 phase/QC 問題誤寫成生物訊號

**建議**

1. `HP0` 與 `HP3` 要獨立統計
2. 作為 `phase/QC annotation`
3. 需要時再做 `unassigned affinity` 類型分析，不直接當 LOH 主量尺

### 7.4 read 數量對分析的重要性

LOH / HP imbalance 觀察不能只看比例，也要看讀數：

1. `hp_ratio = 1.0` 若只有 `1 vs 0` 讀，不應被解讀成強 LOH 線索
2. `hp_ratio = 1.0` 若是 `48 vs 0`，才有較強的結構性意義
3. 因此所有 LOH 觀察都應至少同時回報：
   - `HP1_family_reads`
   - `HP2_family_reads`
   - `effective_hp_reads`
   - `total_reads`

### 7.5 PS / phase block 的邊界

這是後續所有文件都必須遵守的硬限制：

1. 同一個 `PS` / phase block 內，`HP1` 與 `HP2` 才有可直接比較的親緣意義。
2. **不同 phase block 之間，HP1 不保證對應同一親代來源。**
3. 因此跨 block 直接把 `HP1` 串成同一條生物學線索，風險很高。

**目前系統狀態**

1. `ReadParser.hpp` 註解提到會處理 `HP, PS`
2. 但目前 `ReadParser.cpp` 實作只實際抽出 `HP`
3. `reads.tsv` 也沒有 `PS` 欄位

因此目前 InterSubMod 的主輸出**尚未 phase-block-aware**。

這代表：

1. 目前所有 HP-driven 結論都應被視為 `local haplotype label` 層級
2. 若要升級成更強的 allelic / parental interpretation，必須先補 `PS` 匯出與 phase-block-aware diagnostics

---

## 8. paired 與 tumor-only 的研究關係

### 8.1 為什麼要拆成兩條主線

1. paired 有 matched normal，可用來研究 germline ASM、normal baseline、purity-aware correction
2. TO 沒有 matched normal，需要把 PON、LOH、purity、artifact 都重新思考
3. 因此兩者不能共用同一個「有效 / 無效」結論

### 8.2 為什麼又不能完全切開

1. paired 可用來建立 mechanistic reference
2. TO 可用來壓力測試 hard cases 與 FP
3. paired 中有效的概念，未必能直接移植到 TO
4. 但 TO 中失敗得最嚴重的 case，常能告訴 paired 哪些特徵其實只是脆弱訊號

### 8.3 實務原則

1. paired 與 TO 要分開建模、分開結論
2. 但應共用：
   - region schema
   - evidence panel schema
   - hard-case inventory
   - diagnostics template

---

## 9. 哪些問題能先快速驗證，哪些需要完整路線

### 9.1 可快速驗證的問題

| 問題 | 為什麼快 | 建議輸出 |
|---|---|---|
| LOH / HP ratio 在 TP/FP 的分佈 | 既有 CSV 與 summary 已有欄位 | `cross_sample_loh_distribution.tsv` |
| `HP1+HP1-1` vs `HP2+HP2-1` 的 imbalance 是否有樣本差異 | 無需重跑主流程，只需後處理 | `hp_family_balance_summary.tsv` |
| `HP0 / HP3` 是否在特定失敗 case 富集 | 既有 reads / case diagnostics 可補觀察 | `phase_qc_case_notes.md` |
| 同樣本同位點在 paired / TO 是否呈現相同 LOH-like 行為 | 可用現有 canonical / synthesis 產物先做 inventory | `paired_to_same_locus_loh_compare.tsv` |
| `LOH_Subtype` 在 Subclone / Strong / Noise 的分佈 | 已有欄位 | `loh_subtype_by_class.tsv` |

### 9.2 需要中等成本的驗證

| 問題 | 主要缺口 | 建議輸出 |
|---|---|---|
| PS / phase block 是否造成 HP-driven case 誤判 | 主輸出目前沒有 `PS` | `phase_block_qc.tsv` |
| paired / TO 同 locus 是否能形成共用 evidence panel | 需跨 run 對齊 region key 與上游 phasing 輸出 | `same_locus_evidence_panel.tsv` |
| LOH 與 methylation pattern 是否在特定樣本 truly coupled | 需 case-level heatmap + BAM/VCF 驗證 | `loh_case_panel/` |

### 9.3 需要完整路線的問題

| 問題 | 原因 |
|---|---|
| second-hit 順序 inference | 需 normal baseline、CN/purity-aware、gene-level evidence integration 同時到位 |
| TO 的 robust purity-aware correction | 目前缺少足夠的 TO 專用 baseline 與 normal-like reference |
| per-CpG 最終輸出層 | 需先確立 region-level schema 與 evidence panel 邏輯 |
| 單一整合模型穩定提升所有樣本 F1 | 目前跨樣本異質性仍高，paired 與 TO 也不應混成同一結論 |

---

## 10. 建議的最小試錯順序

### Round A：LOH / HP 分佈盤點

目標：

1. 先把 paired 與 TO 的 LOH / HP ratio / HP0 / HP3 分佈整理清楚
2. 確認哪些現象是跨樣本成立，哪些只是單一樣本

### Round B：LOH 先導 evidence panel

目標：

1. 針對代表性的 TP / FP / paired / TO / LongPhase-S / LongPhase-TO case
2. 建出一份可以直接閱讀的 case evidence panel

### Round C：phase-block-aware 缺口盤點

目標：

1. 確認哪些高風險 HP-driven case 需要 `PS`
2. 定義後續要不要把 `PS` 輸出到 `reads.tsv` 或獨立 QC 表

### Round D：paired / TO 分線但互驗

目標：

1. paired 線看 mechanistic clarity
2. TO 線看 hard FP / stress-test
3. 回頭判斷哪些特徵可共用、哪些必須分線

---

## 11. 目前可直接採納的研究原則

1. 五個目標是最終研究願景，不等於近期 phase 順序。
2. 先 `region-first`，後 `per-CpG`。
3. 先 evidence panel，後事件順序 inference。
4. `tumor-only` 是獨立主線，但應設計成可與 paired 互相借鏡。
5. `F1` 與 `可解釋性` 分開回報。
6. LOH 適合作為第一步觀察與 failure-mode mapping，不適合作為第一個主證明。
7. `HP1-1` 歸入 `HP1 family`、`HP2-1` 歸入 `HP2 family` 進行 merged HP 觀察是合理的。
8. `HP0` 與 `HP3` 不應硬併入 LOH 主量尺，應獨立作 phase/QC annotation。
9. 不同 `PS` / phase block 的 `HP1` 與 `HP2` 不能直接串成同一親緣敘事。
10. 若之後要強化 HP-driven 解釋，`PS` 匯出是必要缺口，而不是可忽略細節。

---

## 12. 對接下來文件與簡報的影響

後續若要整理正式簡報或長報告，建議章節主軸依序寫成：

1. InterSubMod 的五個最終研究願景
2. 為什麼當前技術 phase 需要先走 `1 → 2 → 3 → 4`
3. paired 與 TO 為何要分線，但又不能完全分家
4. LOH 為何適合當第一步 observation axis
5. 目前哪些證據已成立，哪些仍是待驗證假說
6. 哪些可以快速試錯，哪些需要完整路線

這樣能把「目標定義」、「方法路線」與「證據層級」分清楚，不再互相打架。
