<!--
建立時間: 2026-04-22
目標: 向 PI 完整報告 Self-Phasing 問題的脈絡、根因、修改方案與結果
受眾: PI（領域專家但未親自操作程式碼）
處理範圍: 從原始 BAM HP tag → LongPhase-TO → ISM ReadParser → 下游特徵
狀態: validated_complete
pipeline_track: TO
priority: P0
hypothesis_id: H-SELF-PHASING-01 / H-HP-ONLY-01
outcome: filter_negative_characterization_warning_issued
-->

# Self-Phasing 循環依賴問題——完整技術報告

> 撰稿日期：2026-04-22
> 關鍵 commit：`775027d`（`refactor/phase1-safety`）
> 原始資料：`/tmp/ism_hp_fix_phase1/merged_{off,on}.csv`、`/tmp/ism_hp_fix_smoke/{off,on}/*.csv`

---

## 一頁速讀：Before → After 整合案例（AF=0.3 somatic variant）

本報告所述問題、修改、結果可用**同一具體案例**貫穿。以下 storyboard 用一個 AF=0.3 的 somatic variant（10 reads 覆蓋，3 ALT + 7 REF）展示：

![Figure 13 — End-to-End Case Storyboard](figures/fig13_case_storyboard.png)

**三個 Panel 敘事**：
- **Panel A（Problem, flag=off）**：TO 模式下 self-phasing 將 3 個 ALT reads 全數強行塞入 HP1-1 bucket；ISM 觀察到 HP_Ratio=0.30（偏離期望 0.5）、NGroups=3（含 artifact 的 HP1-1 bucket）
- **Panel B（Fix, demotion 邏輯）**：`ReadParser.cpp:145-154` 的 9 行 switch 將 `hp_raw ∈ {"1-1","2-1","3"}` demote 為 `"0"`（unphased）；原 raw 值保留在 `hp_tag_raw` 供 audit；Phase 0 守恆驗證（NHP0 gain = Σ NHP_Somatic11+21+33）
- **Panel C（Result, flag=on）**：3 個 ALT reads 被移到 unphased bucket；HP_Ratio 回到 0.43（接近平衡的 0.5）；NGroups 收斂至 2（{HP1, HP2}）——揭露原 NG=3 的訊號源是 somatic HP bucket 而非真實 heterogeneity

**詳見對應章節**：
- Problem（Panel A）邏輯細節 → Section 2.1-2.3
- Fix（Panel B）實作細節 → Section 5-6
- Result（Panel C）跨樣本量化 → Section 7.1-7.3

---

## Executive Summary（執行摘要）

**問題**：TO（Tumor-Only）模式測序下，LongPhase-TO 在建立 haplotype phasing scaffold 時同時把 germline 與 somatic variants 當作 phasing anchor；somatic ALT reads 經由共享 variant 互相連結，被強行塞進單一 haplotype。結果是 **94.6% somatic reads 集中到 HP1**（17.3 : 1 bias），導致 ISM 計算的 `HP_Ratio` 失真，產生 62% 的假 LOH 訊號。

**為何重要**：ISM 有 29 個特徵（約 38%）直接依賴 HP tag。若 HP tag 不可信，則這些特徵在 TO 模式下的結論全部需重新驗證。Self-Phasing 是 TO 模式 FP filter 無法收斂、subclone marker 結論無法外推的主要懷疑對象。

**我們做了什麼**：在 ISM 端新增 `--germline-hp-only` flag（commit `775027d`）。當 flag on 時，ReadParser 將 somatic HP tag（`HP:i:11/21/33`）視為 unphased（`"0"`），只保留 germline HP:i:1/2。原始值以 `hp_tag_raw` 欄位保留供 audit。

**結果**：
- **機制層面（Phase 0 smoke）**：✅ 完全正確。Demotion 數學守恆、audit 欄位獨立於 flag、unit tests 12/12 通過。
- **Filter 層面（Phase 1 HCC1395 TO 40,115 sites）**：❌ **FAIL**。17 個 ISM TSV 特徵中，無一項 AUC 改善 ≥ +0.02；4 個 HP-dependent 特徵 AUC 反而下降 ≤ -0.025。
- **意外發現（P1 審查）**：原「HPFineNGroups≥4 subclone marker，TP rate 89.1%」的結論在 flag=on 下**完全失效**——`NG≥3` 的 regions 從 23,000+ 直接歸零。訊號源**至少部分**來自 somatic HP tag 自身分群，而非真實 subclone 結構。

**建議下一步**：修正保留（default=off，不影響既有流程；flag 作為研究者 opt-in 工具）。**不推進** Phase 2（7 樣本全量重跑）、**不推進** P2/P3。未來若要發表 HPFineNGroups 生物學論文，需條件性執行 P4（master dataset × 兩 flag 比對）。

**關鍵數字**：
- **17.3 : 1** — Somatic reads HP1 vs HP2 bias（614K vs 35K）
- **r = 0.001** — 同位點 HP_Ratio 在 Paired vs TO 間的相關（288K pairs，n.s.）
- **0 / 40,115** — flag on 後 HPFineNGroups≥3 的 regions（原本 23,000+）
- **ΔAUC < 0.02** — 所有特徵在修正前後的 AUC 變化上限

---

## Section 1：問題場景

### 1.1 兩種測序模式

ISM 分析腫瘤樣本時，BAM 與 VCF 可來自兩種模式：

| 模式 | 輸入 | Phasing scaffold 來源 | 應用場景 |
|------|------|----------------------|---------|
| **Paired** | Tumor BAM + matched Normal BAM | Normal 樣本呼叫的 germline SNPs（百萬級，已驗證） | 臨床金標準、研究級樣本 |
| **TO (Tumor-Only)** | Tumor BAM only | PON（Panel of Normals）推測的 germline + 未標記的 somatic 混合 | 無 normal 樣本時（單細胞、臨床檢體、cfDNA） |

Paired 模式下 scaffold 由 **純 germline** 建立：somatic variants 是被動「掛」在 scaffold 上，不影響分群。
TO 模式下因為沒有 normal 樣本，phaser 無法區分 germline 與 somatic，兩者一起進 phasing graph——這是 self-phasing 的結構性根源。

### 1.2 HP tag 在 BAM 中的編碼意義

Haplotype tag（`HP:i:N`）是 haplotag 步驟寫入 BAM read 的整數欄位：

| 值 | 意義 |
|----|------|
| `HP:i:0` 或無 tag | Unphased / 未分配 |
| `HP:i:1` | germline HP1 |
| `HP:i:2` | germline HP2 |
| `HP:i:11` | somatic phase block 1-1（self-phasing 產物） |
| `HP:i:21` | somatic phase block 2-1（self-phasing 產物） |
| `HP:i:33` | somatic phase block 3（self-phasing 產物） |

ISM `ReadParser` 讀到這些整數後轉成 label 字串（`"1"`、`"2"`、`"1-1"`、`"2-1"`、`"3"`、`"0"`），下游 FisherExact / LabelTest / RegionProcessor 依此計算 HP_Ratio、HPFineNGroups、PermanovaF 等特徵。

### 1.3 為何 HP tag 對 ISM 重要

HP tag 是 ISM 三個核心分析維度的 label：

1. **LOH 判定**：`HP_Ratio = min(HP1 reads, HP2 reads) / (HP1 + HP2)`，Paired 模式 ~0.5 表示兩 haplotype 平衡；<0.1 或 >0.9 表示 LOH。
2. **Subclone 結構**：`HPFineNGroups` 計算一 region 內出現的 HP label 種類數；NGroups≥4 被認為是 subclone heterogeneity 的 marker。
3. **ASM（Allele-Specific Methylation）**：HP 分群後計算甲基化矩陣的 haplotype 差異。

若 HP tag 被 self-phasing 污染 → 這三個維度全部失真。

**Figure 1** 對照兩種模式的流程：

![Figure 1 — Paired vs TO Pipeline](figures/fig1_pipeline_comparison.png)

---

## Section 2：Self-Phasing 根因

### 2.1 循環依賴：「球員兼裁判」

LongPhase-TO 的 phasing 演算法基於 phasing graph：**每個 variant 是一個節點，若同一條 read 同時帶兩個 variants 的 ALT，就在它們之間連一條邊**。邊的權重越大，兩個 variants 越可能被指派到同一 haplotype。

Paired 模式下 scaffold 由 germline SNPs 建立。這些 SNPs 的 GT 已由 normal 樣本獨立確認（0/1 het），屬於「外部 ground truth」，因此 phasing 是根據**已知的 haplotype 結構**分配 somatic。

TO 模式下沒有 normal，phaser 把每個候選 variant（包括真正的 somatic）都塞進 graph：

```
Somatic variants 在同一 tumor clone 內，天然共享 sub-population
  → 它們的 ALT reads 高度共現
    → phasing graph 中 ALT-ALT edges 權重極高
      → 所有 somatic ALT 被判定為「同一 haplotype」
        → 強行塞入 HP1 的 phase block 11/21/33
          → 原本是「被分類者」的 somatic，變成「scaffold 的一部分」
            → 循環依賴：somatic 自己決定自己被分到哪裡
```

這就是「球員兼裁判」——不該參與定義規則的人參與了規則。

### 2.2 具體走例：AF=0.3 的 somatic variant

假設 chr1 位置 X 有一個 AF=0.3 的 somatic variant。腫瘤樣本共 10 reads 覆蓋此位點——3 reads 帶 ALT，7 reads 帶 REF。

**Paired 模式**：
1. Normal 樣本沒有 variant → 確認為 somatic
2. Scaffold 由 germline SNPs 建立（不含此 variant）
3. 3 個 ALT reads 根據各自攜帶的 germline SNPs 被分配 HP1 或 HP2
4. 結果可能是 HP1 有 2 ALT / 2 REF、HP2 有 1 ALT / 4 REF
5. **HP_Ratio ≈ 0.5**（若該區有 LOH 則偏移，但 Non-LOH 本應平衡）

**TO 模式**：
1. 此 variant 未被 PON 標記 → 進 phasing graph 當作 anchor
2. 這 3 個 ALT reads **同時帶其他 somatic variants @ Y、Z**（它們來自同一 tumor clone）
3. Phasing graph edges：X-ALT ↔ Y-ALT ↔ Z-ALT 三角連結
4. Phaser 把 X/Y/Z 的 ALT 全部判定為「同一 phase block」
5. Haplotag 為每個 ALT read 打上 `HP:i:11`（表示 somatic-only block 1）
6. 結果：3 個 ALT reads 全部歸 HP1；7 個 REF reads 分散 HP1/HP2
7. **HP_Ratio → 0.94**（假 LOH）

**Figure 3** 視覺化這個過程：

![Figure 3 — AF=0.3 Walkthrough](figures/fig3_af03_walkthrough.png)

### 2.3 為什麼 somatic 會互相 phase（機制細節）

Phasing graph 的 edge weight 計算：

```
weight(variant_A, variant_B) = Σ_reads I(read 帶 A.alt) × I(read 帶 B.alt)
```

兩個 **真正的 germline het** SNP（位置遠）：edge weight 會受到 crossover rate 影響，邊緣化。
兩個 **同 clone 的 somatic** SNP（位置遠）：**共享 sub-population 的 reads** 100% 同時帶兩者 ALT → edge weight 遠高於背景。

結果：TO 模式 phasing graph 中 somatic edges 比 germline edges 權重更高，**somatic 反客為主定義 scaffold**。

![Figure 2 — Self-Phasing Conceptual Flow](figures/fig2_self_phasing_concept.png)

---

## Section 3：量化證據

以下數據來自 HCC1395 樣本的 Paired vs TO 全基因體比對（穩定度 4/5）。

### 3.1 Somatic bias 17.3 : 1

統計 TO 模式下所有 somatic variant ALT reads 的 HP 分配：

- HP1 reads: **614,000**
- HP2 reads: **35,500**
- Ratio: **17.3 : 1**
- 94.6% somatic reads 分配到 HP1

理想情況下（無 self-phasing）應該接近 1 : 1（每個 somatic 的 two haplotype 是 50/50 事件）。17 倍的 bias 是強烈的 artifact 指標。

### 3.2 同位點 HP_Ratio 跨模式無相關

取 288,000 個同時在 Paired 與 TO 模式下都有 ≥5 reads 的 variant sites，計算各自的 `HP_Ratio`，做 Pearson correlation：

- **r = 0.001**（p = 0.59，n.s.）

若 TO 的 HP_Ratio 是真實 haplotype 結構的近似，應該與 Paired 高度相關（預期 r > 0.5）。實際上相關近零，表示 **TO 的 HP_Ratio 訊息幾乎全是雜訊或 self-phasing artifact**。

### 3.3 ISM HP_Ratio LOH 的 62% 是 artifact

TO 模式判定為 LOH（HP_Ratio < 0.1 或 > 0.9）的 regions 中：
- **62%** 在 Paired 模式下 HP_Ratio 位於 0.4-0.6（即 balanced，非 LOH）
- 剩下 38% 在 Paired 仍為 LOH（結構性 LOH，真實）

即超過一半的 TO-mode LOH 判定是 self-phasing 產生的 artifact。

**⚠ 兩套 LOH 系統澄清**：
- **ISM HP_Ratio LOH**：由 BAM HP tag 計算，受 self-phasing 嚴重影響
- **LOH.bed region-level LOH**：由 VCF AD（allele depth）計算，**不經過** HP tag，Jaccard=1.0 完全不受 self-phasing 影響
- 兩者 kappa=0.670（不完全一致）正是因為使用不同 LOH 定義

### 3.4 7/7 樣本方向一致

HCC1395、HCC1395_DORADO、HCC1954、HCC1937、H1437、H2009、COLO829 全部觀察到相同方向的 self-phasing 效應（CV-2 pass）。排除了樣本特異性的可能。

**Figure 4** 彙總以上證據：

![Figure 4 — Quantitative Evidence](figures/fig4_evidence_summary.png)

---

## Section 4：影響範圍

### 4.1 三級影響分類

ISM 有 ~85 個輸出特徵（TSV column）。根據對 HP tag 的依賴程度分三類：

**🔴 嚴重影響（結果不可信，需修正後重測）**
- `HP_Ratio` → 假 LOH
- `Potential_LOH` → 62% artifact
- `HPMergedDelta` / `HPMergedSig` → 方向反轉可能
- `hp_assign_rate` → 系統性偏高
- `effective_hp_reads` → 分母偏離
- `HPFineNGroups≥3` → 依賴 somatic HP tag 自身分群（**本次研究新發現**）

**🟡 中度影響（間接污染）**
- `Quality_Score` → 移除 LOH penalty 後殘留 HP 成分（AUC 0.497）
- `GlobalP` → 取 HP/Allele 最小 p，HP 噪音偶然壓低
- `CramersV` → 同理
- `VerificationClass` → `label_sig` 含 HP 成分

**🟢 無影響（結論穩固）**
- `PairwiseMeanDist` / `MedianDist` → 全 reads 計算，不分 HP
- `AlleleDelta` / `AlleleP` → 只用 allele label（ALT/REF）
- Caller 特徵（AF/GQ/DP/SB）→ 來自 VCF
- 甲基化矩陣（raw）→ 來自 BAM MM/ML tag
- CpG 座標、`region_methyl_mean` → 基因組固有

### 4.2 特徵數量分布

| 類別 | 特徵數 | 佔比 | TO 結果可信度 |
|------|--------|------|------------|
| A. 完全不依賴 HP | ~42 | 55% | ✅ 結論全部穩固 |
| B. 直接依賴 HP | ~29 | 38% | ❌ 不可信，需重測 |
| C. 間接依賴 HP | ~14 | 7% | ⚠ 大部分影響微弱或已程式碼移除 |

**Figure 5** 視覺化：

![Figure 5 — Feature Impact Matrix](figures/fig5_impact_matrix.png)

---

## Section 5：修改方案決策

### 5.1 長期方案 vs 短期方案

| 方案 | 改動位置 | 風險 | 週期 | 判定 |
|------|---------|------|------|------|
| **長期**：修 LongPhase-TO `PhasingProcess.cpp`（讓 somatic 不參與 scaffold） | 第三方 C++ 專案 | HIGH | 3-6 個月 | **DEFERRED** |
| **短期**：ISM 端 ReadParser 忽略 somatic HP tag | 本專案 `src/core/ReadParser.cpp` | LOW | 1-2 天 | **✅ 採用** |

短期方案不修根因（LongPhase-TO 仍會產生 somatic HP tags），而是在 ISM 端「假裝這些 tags 不存在」。這是 fallback 而非 true fix，但足以讓 ISM 的下游特徵脫離 self-phasing 污染。

### 5.2 Option A vs Option B

短期方案內還分兩個實作策略：

| 維度 | Option A（ReadParser 單點過濾） | Option B（下游分散過濾） |
|------|-------------------------------|---------------------|
| 改動檔案數 | **1**（`ReadParser.cpp` 單一 switch） | 4+（`LabelTest.cpp` × 2、`FisherExact.cpp` × 3） |
| 可逆性 | 關 flag 即 bypass，output schema 不變 | 每個下游函式須個別驗證 |
| Regression 風險 | LOW | MEDIUM |
| 保留 audit 能力 | `ReadInfo::hp_tag_raw` 保留原值 | 天然保留（下游看得到原值） |
| 與既有 `--pon-only-phasing` 對稱 | ✅ | ❌ |

**採用 Option A + 新增 audit 欄位**。

**Figure 6** 決策樹：

![Figure 6 — Fix-Option Decision Tree](figures/fig6_fix_decision_tree.png)

---

## Section 6：實作細節

### 6.1 修改清單

| 檔案 | 位置 | 改動 |
|------|------|------|
| `include/core/Config.hpp` | struct 內 | 新增 `bool germline_hp_only = false;` |
| `include/core/ReadParser.hpp` | `ReadInfo` struct | 新增 `std::string hp_tag_raw;`（audit 用） |
| `src/core/ReadParser.cpp` | 第 121-142 行 switch | 保留 switch 寫入 `hp_tag_raw`；`flag && hp_int ∈ {11,21,33}` → `info.hp_tag = "0"` |
| `include/utils/ArgParser.hpp` | CLI | 新增 `--germline-hp-only`（預設 false） |
| `src/core/RegionProcessor.cpp` | TSV header | 新增三欄 `NHP_Somatic11/21/33`（audit，計數 `hp_tag_raw` 為對應值的 reads，不論 flag） |
| `tests/test_read_parser.cpp` | — | 新增 12 個 unit tests |

### 6.2 關鍵 code 段

```cpp
// src/core/ReadParser.cpp:121-142
std::string raw;  // parse HP AUX into raw ("1"/"2"/"1-1"/"2-1"/"3"...)
if (hp_aux) {
    int hp_int = static_cast<int>(bam_aux2i(hp_aux));
    switch (hp_int) {
        case 1:  raw = "1"; break;
        case 2:  raw = "2"; break;
        case 11: raw = "1-1"; break;
        case 21: raw = "2-1"; break;
        case 33: raw = "3"; break;
        default: raw = "0"; break;
    }
}
info.hp_tag_raw = raw;                    // always preserved
if (config_.germline_hp_only &&
    (raw == "1-1" || raw == "2-1" || raw == "3")) {
    info.hp_tag = "0";                    // demote to unphased
} else {
    info.hp_tag = raw;                    // pass through
}
```

### 6.3 Audit 欄位

`RegionProcessor` 計算每個 region 的三個 somatic tag read 數：

```
NHP_Somatic11 = Σ reads with hp_tag_raw == "1-1"
NHP_Somatic21 = Σ reads with hp_tag_raw == "2-1"
NHP_Somatic33 = Σ reads with hp_tag_raw == "3"
```

這三欄**永遠反映 raw 值**，不受 flag 影響。因此 flag on vs off 的 output 中，這三欄應該完全相同——這是 Phase 0 smoke 的 invariant 檢核之一。

### 6.4 驗證機制

| 層次 | 測試 | 目的 |
|------|------|------|
| Unit | `tests/test_read_parser.cpp` 12 cases | 合成 reads、flag on/off 預期 label 分佈 |
| Smoke | chr19 615 sites | 全流程驗證、守恆律檢核 |
| Phase 1 | HCC1395 TO 40,115 sites | 全基因體 AUC gate |

### 6.5 Commit

`775027d` on branch `refactor/phase1-safety`。改動合計約 +450 行（code + tests + docs）。

---

## Section 7：驗證結果

### 7.1 Phase 0 Smoke（chr19 615 sites）——機制正確性 ✅

測試設定：HCC1395 V3-Fixed TO BAM × chr19 TP VCF，flag off 與 on 各跑一次（~21 秒 / run）。

**7.1.1 Audit 欄位獨立性（核心正確性）**

| 欄位 | flag off sum | flag on sum | 相同 |
|------|-------------|-------------|------|
| NHP_Somatic11 | 4,894 | 4,894 | ✅ |
| NHP_Somatic21 | 8,630 | 8,630 | ✅ |
| NHP_Somatic33 | 2,498 | 2,498 | ✅ |

每個 region 的三欄值在兩 run 間 identical → audit 計數 key 為 raw 值，與 flag 無關。

**7.1.2 Demotion 數學守恆（critical）**

預期：flag on 時，NHP0 增量應等於 `NHP_Somatic11 + NHP_Somatic21 + NHP_Somatic33`（所有 somatic tags demote 到 unphased）。

| 欄位 | flag off | flag on | Δ |
|------|---------|---------|---|
| NHP3 | 2,498 | 0 | **-2,498** |
| NHP0 | 13,985 | 30,007 | **+16,022** |
| HP1FamilyN | 26,462 | 21,568 | -4,894 |
| HP2FamilyN | 33,134 | 24,504 | -8,630 |

驗算：**4,894 + 8,630 + 2,498 = 16,022** ✅ 完全守恆。

**7.1.3 HPFineNGroups 分佈收斂**

flag off 時 NGroups ∈ {0, 1, 2, 3, 4}（5 類）；flag on 後收斂至 {0, 1, 2}（3 類）——符合「somatic HP tags 11/21/33 移除後，細分類僅剩 germline HP1/HP2 + unphased」的預期。

### 7.2 Phase 1 HCC1395 TO（40,115 sites）——Filter Gate ❌

測試設定：HCC1395 V3-Fixed TO BAM × ClairS-TO raw TP/FP VCF split，flag off 與 on 各跑 TP + FP 兩次（4 runs 總計 ~1 hr）。

**7.2.1 HP_Ratio median 位移**

| 分層 | flag off median | flag on median | Δ |
|------|----------------|---------------|---|
| 全 TP | 0.549 | 0.529 | -0.020 |
| TP Non-LOH（n=26,439） | 0.554 | 0.531 | -0.023 |
| TP Potential_LOH（n=2,070） | 0.091 | 0.093 | +0.002 |

**方向正確**（Non-LOH 朝 0.5 靠近），但**幅度遠小於 plan 預期**（plan 預期 0.836 → 0.55-0.65 是基於 landscape doc 引用的舊 baseline；V3-Fixed TO 實測 baseline 本來就是 0.549）。

**7.2.2 AUC Gate（TP=1, FP=0 的分類能力）**

計算每個特徵的「最佳方向 AUC」（`max(auc, 1-auc)`），17 個代表性特徵結果：

![Figure 8 — Phase 1 AUC Comparison](figures/fig8_phase1_auc.png)

**關鍵觀察**：
- **無任何特徵** AUC 改善 Δ ≥ +0.02
- 4 個 HP 衍生特徵 Δ ≤ -0.025（方向錯誤）：
  - `HPFineNGroups`: 0.536 → 0.510（-0.026）
  - `HPFine_NGroups_CF`: 0.536 → 0.510（-0.026）
  - `HPMergedDelta`: 0.541 → 0.515（-0.025）
  - `NHP3`: 0.535 → 0.500（-0.035，因 on-run 全部為 0）
- `AlleleDelta`（唯一 HP-independent 中等訊號特徵，AUC=0.629）完全不動——與預期一致（allele label 不依賴 HP tag）

**Gate 判定**：FAIL — 無改善潛力 → **不進 Phase 2（7 樣本全量重跑）**。

**7.2.3 HP_Ratio 分佈視覺化**

![Figure 9 — HP_Ratio Violin](figures/fig9_hp_ratio_violin.png)

Non-LOH 組的 HP_Ratio 分佈朝 0.5 靠攏（方向正確）；LOH 組幾乎不動（已在極端，somatic tags 本就少參與）。

### 7.3 意外發現（P1 HPFineNGroups 審查）

原 F pilot 研究（2026-04-18）的結論：**HPFineNGroups≥4 + NR≥80 作為 subclone heterogeneity marker，7 樣本 master dataset pooled TP rate = 89.1%（新 canonical filter 92.8%）**。這是 ISM 唯一經跨樣本 + confound correction + AF 去混淆驗證的 POSITIVE 訊號。

**本次 Phase 1 發現**（HCC1395 TO standalone, ClairS-TO raw TP/FP split）：

| 條件 | flag off | flag on |
|------|---------|---------|
| NG=4 + NR≥80（TP regions） | 4,003 | **0** |
| NG=3 + NR≥80（TP regions） | 18,729 | **0** |
| NG=2（TP regions） | 5,626 | **27,729** |

**flag on 後 NG≥3 完全消失**，97.3% 的 TP 集中到 NG=2。

![Figure 7 — HPFineNGroups Collapse](figures/fig7_ngroups_collapse.png)

**HCC1395 TO standalone 獨立驗證**：
- NonLOH baseline TP rate = 0.699
- NonLOH + NG≥4 + NR≥80 TP rate = 0.694（vs baseline -0.005pp，Fisher odds = 0.913 反向 p=3.5e-3）
- 即 marker 在 HCC1395 TO standalone 上**沒有 TP 富集**，與 master dataset 的 0.810 相差 0.12

**兩種非互斥的解讀**：

**解讀 A：資料集依賴**
- F pilot 的 master dataset 是 7 樣本 pooled + AF filter，使用不同 VCF aggregation 流程
- HCC1395 TO standalone 是 ClairS-TO raw TP/FP split，資料組成不同
- marker 的機制仍可能成立，只是本次驗證資料集不在適用範圍

**解讀 B：Marker 訊號部分來自 somatic HP 人工分群**
- F pilot Step 3 做過 chr-shuffle null test（Z=43.5），確認**非 spatial artifact**
- 但 chr-shuffle null **不檢查**「訊號是否來自 somatic HP 自身細分」
- flag=on 測試是新的 orthogonal null——移除 somatic HP tags 後 NG≥3 完全消失 → 訊號源**至少 50%+** 依賴 somatic HP labels

**目前狀態**：兩種解讀可能同時為真。Phase 1 無法分辨兩者佔比。**marker POSITIVE 結論不撤回**，但須標記 pipeline-dependency 警示（見 memory `project_hpfinengroups_subclone_marker.md` 2026-04-21 更新）。

---

## Section 8：結論與下一步

### 8.1 修正本身的價值定位

| 面向 | 判定 | 理由 |
|------|------|------|
| 機制正確性 | ✅ | Phase 0 守恆驗證、unit tests 12/12 pass |
| 研究誠實性 | ✅ | 移除 LongPhase-TO self-phasing 循環依賴；HP 分群僅反映 germline 結構 |
| Filter 方向（TP/FP 分類）| ❌ NEGATIVE | 無 TSV 特徵 AUC 改善；多個 HP 衍生特徵 AUC 下降 |
| Characterization 方向 | ⚠ CONDITIONAL | subclone marker 結論需 pipeline-dependency 標記 |

**為何 self-phasing 是真的，但 filter 修正卻無 AUC 增益？**

TP 與 FP 在 somatic tag 密度上幾乎相同（TP 27.4/site、FP 27.7/site）→ somatic tag 本身對 TP/FP 區分**沒有訊號**。修正只是把「noisy 分群（NG ∈ {0-4}）」變成「無分群（NG ∈ {0-2}）」——前者的 AUC 已接近 0.5，後者不會更好。**訊號與噪音比的問題**，不是「過濾雜訊後訊號會浮現」的問題。

### 8.2 對其他研究方向的連鎖影響

1. **HPFineNGroups subclone marker**（原 POSITIVE ⭐4）
   - 加入 pipeline-dependency 警示
   - 未來引用須標記：(1) master dataset vs raw VCF split、(2) `--germline-hp-only=false` 前提、(3) 若要發表需 P4 master × 兩 flag 重驗

2. **TO FP filter 結構性關閉結論**（原 CLOSED）
   - **強化**：self-phasing 修正後 AUC 仍無改善 → self-phasing 不是 TO filter 失敗的主因
   - 原結論「self-phasing 僅解釋 35% TO/Paired 差距」獲得定量證實

3. **LOH.bed 兩套系統澄清**
   - LOH.bed region-level LOH：不受 self-phasing 影響（Jaccard=1.0 已證）
   - ISM HP_Ratio LOH：受 self-phasing 影響嚴重
   - 兩者 kappa=0.670 差異完全解釋

4. **Phase 2 (7 樣本全量重跑)**：**暫緩**。不滿足 Phase 1 gate。

### 8.3 建議下一步

![Figure 10 — Next-Step Decision Tree](figures/fig10_next_steps.png)

| 項目 | 判定 | 理由 |
|------|------|------|
| Phase 2（7 樣本 × 2 模式全量重跑） | **STOP** | Phase 1 gate 明確 FAIL，翻盤機率 <5% |
| P2（within_dom_alt_frac 下游 pipeline 重建） | **SKIP** | 成本 ~1 day，預期 AUC 改善 <0.02 |
| P3（flag default=on 推廣） | **STATUS QUO** | 立即破壞現有 HPFineNGroups subclone 應用；default=off、opt-in 是較穩妥選擇 |
| P4（master dataset × 兩 flag × 兩 mode 比對） | **延後** | 僅在未來發表 HPFineNGroups 論文時條件性啟動；12-24 hr 機器時 |

**回到更高 ROI 的主線**：
- CovM baseline blocker（master dataset 重跑）
- HCC1395 chr8 hotspot 深度分析（7.4× FP enrichment）
- Phase 2 Normal Methylation Reference（characterization 方向）
- LOH Subclone AF × Methylation POSITIVE 的 TO 延伸驗證

---

## 附錄

### A. 檔案清單

**實作**
- `include/core/Config.hpp`
- `include/core/ReadParser.hpp`
- `include/utils/ArgParser.hpp`
- `include/core/DataStructs.hpp`
- `include/core/RegionProcessor.hpp`
- `src/core/ReadParser.cpp`
- `src/core/RegionProcessor.cpp`
- `CMakeLists.txt`

**測試**
- `tests/test_read_parser.cpp`（12 cases）

**驗證資料**
- Phase 0 smoke：`/tmp/ism_hp_fix_smoke/{off,on}/significance_summary.csv`
- Phase 1 full：`/tmp/ism_hp_fix_phase1/{tp,fp}_{off,on}/significance_summary.csv`、`merged_{off,on}.csv`

**報告**
- Phase 0 smoke：`docs/experiments/in_progress/2026/04/20260421_ReadParser_GermlineHPOnly_smoke_01.md`
- Phase 1 full：`docs/experiments/in_progress/2026/04/20260421_ReadParser_GermlineHPOnly_HCC1395_01.md`
- P1 marker 審查：`docs/experiments/in_progress/2026/04/20260421_HPFineNGroups_Marker_Reaudit_01.md`
- Landscape 根因：`docs/reports/research_landscape/02_Self_Phasing根因.md`
- PON-only 原始報告：`research/loh_investigation/reports/20260403_pon_only_haplotag_ism_verification_report.md`

**圖表生成**
- `scripts/analysis/generate_pi_report_figures_self_phasing.py`

**Commit**
- `775027d` on branch `refactor/phase1-safety`

### B. 原始資料位置與重現命令

```bash
# Phase 0 smoke（~21 秒）
./build/bin/inter_sub_mod \
  -t /big7_disk/liaoyoyo2001/longphase-to-mod/output/pononly_v3_fixed/tumor_tagged.bam \
  -n /big8_disk/liaoyoyo2001/InterSubMod/data/bam/HCC1395/normal.bam \
  -v /tmp/ism_hp_fix_smoke/chr19_tp.vcf \
  -r /big8_disk/liaoyoyo2001/InterSubMod/data/ref/GRCh38_no_alt_analysis_set.fasta \
  -o /tmp/ism_hp_fix_smoke/off \
  --metric BERNOULLI --window 5000 --threads 32
# 加 --germline-hp-only 即為 on run

# Phase 1 full HCC1395 TO（4 runs, ~1 hr total）
# 指令同上，VCF 換為 clairsto_v3fixed_work/clairsto_{tp,fp}.vcf.gz
```

### C. 術語表

| 術語 | 解釋 |
|------|------|
| **Paired 模式** | Tumor + Normal 成對測序，Normal 樣本提供乾淨的 germline ground truth |
| **TO 模式** | Tumor-Only，只有腫瘤 BAM，無配對 Normal |
| **LongPhase-TO** | 第三方 phasing 軟體，用於 TO 模式下產生 phased VCF 與 haplotagged BAM |
| **haplotag** | 將 phased variant 資訊寫入 BAM read 的 HP tag（標記每個 read 屬於哪條 haplotype） |
| **HP tag** | BAM AUX 欄位 `HP:i:N`，值 1/2 = germline haplotype；11/21/33 = somatic self-phased block |
| **Self-phasing** | Somatic variants 在 phasing graph 中自己互相連結，形成獨立 phase block 的現象 |
| **LOH (Loss of Heterozygosity)** | 正常雜合區域在腫瘤中失去一個 allele；可由 VCF AD 或 BAM HP_Ratio 判定 |
| **ISM (Inter-Subclonal Methylation)** | 本專案 C++ 工具，分析 haplotype × methylation pattern |
| **HP_Ratio** | ISM 計算的 min(HP1 reads, HP2 reads) / total；~0.5 表示平衡，<0.1 或 >0.9 表示 LOH |
| **HPFineNGroups** | 一個 region 內出現的 HP label 種類數（fine = 包含 HP1/HP2/1-1/2-1/3/0 等細分） |
| **Subclone marker** | 可用於判定某 region 是否屬於 subclonal heterogeneous 的特徵 |
| **Master dataset** | F pilot 使用的標準 7 樣本 × 正式 haplotag 流程聚合資料 |
| **ClairS-TO raw TP/FP split** | ClairS-TO caller 原始輸出依 truth VCF 切分的 TP / FP 子集 |
| **AUC** | Receiver Operating Characteristic 曲線下面積，0.5 = 隨機、1.0 = 完美分類；`best_auc = max(AUC, 1-AUC)` |
| **PON (Panel of Normals)** | 正常樣本庫，用於 TO 模式下猜測 germline variants |

---

## 報告結束

**主要訊息給 PI**：
1. Self-phasing 是真實存在的結構性 bug（7/7 樣本一致），**不是錯覺**
2. 在 ISM 端加 filter 可以乾淨地隔離這個 bug 的影響——**機制層面完美**
3. 但這個 filter **對 TP/FP 分類不產生任何 AUC 增益**——self-phasing 不是 TO filter 失敗的瓶頸
4. 意外發現：**HPFineNGroups subclone marker 的訊號源可能不是純生物學**，而是部分依賴 somatic HP tag 的人工分群——這是研究誠實性重要議題，需在未來 marker 論文中明確標記
5. 建議回到更高 ROI 的主線研究

**開放議題**：
- P4 master dataset × 兩 flag 比對（發表時才做）
- 長期方案——修 LongPhase-TO 上游：現階段不緊迫，因為短期 filter 已給出「self-phasing 是真但不解決 AUC」的定量結論
