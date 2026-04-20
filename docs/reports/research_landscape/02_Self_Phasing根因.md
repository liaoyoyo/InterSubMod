<!--
建立時間: 2026-04-04 21:20
目標: 完整說明 self-phasing circular dependency 的機制、影響與修正狀態
處理範圍: Self-phasing 因果鏈 S0-S5 + PON-only phasing 測試結果
關聯檔案:
  - research/loh_investigation/reports/20260403_pon_only_phasing_verification_report.md
  - research/loh_investigation/reports/20260403_pon_only_haplotag_ism_verification_report.md
  - docs/reports/research_landscape/04_暫停判定與重評估.md
-->

# Self-Phasing 根因：為什麼 TO 的 HP tag 不可信？

---

## 一句話說明

**LongPhase-TO 讓 somatic mutation 自己決定自己被分到哪個 haplotype，形成「球員兼裁判」的循環依賴。**

---

## 正確流程 vs 問題流程

```mermaid
graph TB
    subgraph Paired["✅ Paired 模式（正確）"]
        direction TB
        P1["Step 1: Normal 樣本<br/>呼叫 germline SNPs<br/>（百萬級、已驗證）"]
        P2["Step 2: Germline SNPs<br/>建立 phasing scaffold<br/>（兩條 haplotype 骨架）"]
        P3["Step 3: Somatic variants<br/>不參與 scaffold 建立<br/>只是被動接受分配"]
        P4["Step 4: Haplotag<br/>根據 germline scaffold<br/>分配 reads → HP1 或 HP2"]
        P5["結果: HP 分配<br/>反映真實 haplotype"]

        P1 --> P2 --> P3 --> P4 --> P5
    end

    subgraph TO["❌ TO 模式（循環依賴）"]
        direction TB
        T1["Step 1: 沒有 Normal 樣本<br/>用 PON 猜測 germline<br/>（+未標記的 somatic）"]
        T2["Step 2: Germline + Somatic<br/>混合建立 phasing scaffold<br/>⚠️ Somatic 也當 anchor"]
        T3["Step 3: Somatic variant 的<br/>ALT reads 互相連結<br/>→ 被拉向同一 HP"]
        T4["Step 4: Haplotag<br/>被 somatic-contaminated<br/>scaffold 驅動"]
        T5["結果: 94.6% somatic<br/>reads → HP1<br/>bias 17.3:1"]

        T1 --> T2 --> T3 --> T4 --> T5
    end

    style Paired fill:#e8f5e9
    style TO fill:#ffebee
```

---

## 用具體數字走過一個例子

假設有一個 **AF=0.3 的 somatic variant**（chr1 位置 X），腫瘤中 30% reads 帶 ALT，70% 帶 REF。

### Paired 模式

1. Normal 樣本沒有這個 variant → 確認為 somatic
2. Phasing scaffold 由 germline SNPs 建立（不含此 variant）
3. Haplotag 根據 scaffold 分配 reads
4. 30% ALT reads 隨機分到 HP1/HP2（取決於此 variant 實際位於哪條 haplotype）
5. **HP_Ratio ≈ 0.5**（如果 LOH 則會偏移）

### TO 模式

1. 此 variant 未被 PON 標記 → 被視為 germline 或 unknown
2. 在 phasing graph 中，此 variant 的 ALT reads 互相連結
3. 這些 ALT reads 上的其他 somatic variants 也連結進來
4. 形成一個巨大的 ALT-read cluster → **全部指向 HP1**
5. **HP_Ratio → 0.94**（偏向 HP1，產生假 LOH）

```mermaid
graph LR
    subgraph SomaticReads["Somatic Variant 的 ALT Reads"]
        R1["Read A<br/>帶 ALT @ pos X<br/>帶 ALT @ pos Y"]
        R2["Read B<br/>帶 ALT @ pos X<br/>帶 ALT @ pos Z"]
        R3["Read C<br/>帶 ALT @ pos Y<br/>帶 ALT @ pos Z"]
    end

    subgraph Graph["Phasing Graph"]
        G1["pos X ALT ─── pos Y ALT<br/>(via Read A)"]
        G2["pos X ALT ─── pos Z ALT<br/>(via Read B)"]
        G3["pos Y ALT ─── pos Z ALT<br/>(via Read C)"]
        G4["結果: X,Y,Z 的 ALT<br/>全部在同一 phase block<br/>→ 全指向 HP1"]
    end

    R1 --> G1
    R2 --> G2
    R3 --> G3
    G1 --> G4
    G2 --> G4
    G3 --> G4

    style SomaticReads fill:#ffcdd2
    style Graph fill:#fff9c4
```

---

## 量化證據總表（穩定度 4/5）

![Self-Phasing 影響量化](figures/03_self_phasing_impact.png)

| 證據項 | 數值 | 解讀 |
|--------|------|------|
| Somatic HP1:HP2 bias | **17.3:1** (614K vs 35K) | 94.6% somatic reads 偏向 HP1 |
| TO TP **ISM HP_Ratio LOH** 中 self-phasing 造成的比例 ⁽¹⁾ | **62%** 移除 self-phasing 後消失 | AF 0.1-0.8 近 100% |
| 全 TO **ISM HP_Ratio LOH** 中 self-phasing artifact ⁽¹⁾ | **31.2%** | 其餘 68.8% 為 structural ISM LOH |

> ⁽¹⁾ **LOH 層次說明（2026-04-19 P2-A 補註）**：此處 62%/31.2% 指 **ISM HP_Ratio LOH**（BAM HP tag 路徑），非 **LOH.bed region-level LOH**（VCF AF/VAF 路徑）。LOH.bed 在 PON-only 實驗中 Jaccard=1.0 完全不受 self-phasing 影響。
| 同位點 HP_Ratio 跨模式相關 | **r = 0.001**（288K pairs） | 完全不相關 |
| TO-only LOH 在 paired 下完全平衡 | **86.5%** | HP_Ratio 0.4-0.6 |
| Cohen's d (HP_Ratio 差異) | **-1.20** | 巨大效應量 |
| 跨樣本方向一致性 | **7/7 samples** | CV-2 全通過 |

> **7 個樣本全部觀察到相同方向的 self-phasing 效應，排除了樣本特異性的可能。**

---

## Self-Phasing 對 ISM 各輸出欄位的影響

### 影響程度分級

```mermaid
graph TD
    subgraph Severe["🔴 嚴重影響（結果不可信）"]
        S1["HP_Ratio → 假 LOH"]
        S2["Potential_LOH → 62% artifact"]
        S3["HPMergedDelta/Sig → 方向反轉"]
        S4["hp_assign_rate → 偏高"]
        S5["effective_hp_reads → 偏離"]
    end

    subgraph Moderate["🟡 中度影響（間接污染）"]
        M1["QualityScore → AUC 0.497<br/>(已移除 LOH penalty)"]
        M2["GlobalP → 取 HP/Allele 最小值<br/>HP 噪音可能偶然壓低"]
        M3["CramersV → 取 HP/Allele 最大值"]
        M4["VerificationClass → label_sig<br/>含 HP 成分"]
    end

    subgraph None["🟢 不受影響（結論穩固）"]
        N1["PairwiseMean/MedianDist<br/>全 reads 計算，不分 HP"]
        N2["AlleleDelta/AlleleP<br/>只用 allele label (ALT/REF)"]
        N3["Caller 特徵 (AF/GQ/DP/SB)<br/>來自 VCF"]
        N4["甲基化矩陣 (raw)<br/>來自 BAM MM/ML tag"]
        N5["CpG 座標、region_methyl_mean<br/>基因組固有特性"]
    end

    style Severe fill:#ffcdd2
    style Moderate fill:#fff9c4
    style None fill:#c8e6c9
```

### 完整特徵分類數量

| 分類 | 特徵數 | 佔比 | TO 結果狀態 |
|------|--------|------|------------|
| **A. 完全不依賴 HP** | ~42 | 55% | 結論全部穩固，不需重測 |
| **B. 直接依賴 HP** | ~29 | 38% | TO 結果不可信，需修正後重測 |
| **C. 間接依賴 HP** | ~14 | 7% | 大部分影響微弱或已程式碼移除 |

> 詳細分類表見 [03_ISM分析價值界定.md](03_ISM分析價值界定.md)

---

## PON-Only Phasing 修正測試

### 修改內容

在 `LongPhase-TO` 的 `PhasingProcess.cpp` 中，first-pass phasing 前呼叫 `convertNonGermlineToSomatic()`，使 somatic variants 以 reduced edge weight 處理。新增 `--pon-only-phasing` CLI flag。

### VCF 層面結果（✅ 全面正確）

| 指標 | Baseline | PON-only | 判定 |
|------|----------|----------|------|
| LOH.bed Jaccard | — | **1.0000** | 完全一致 |
| Somatic HP1:HP2 bias | 17.3:1 | **消除** | self-phasing 解除 |
| Phase block N50 | 4,061 | **8,109 (+99.7%)** | 品質翻倍 |
| Phased rate | 54.9% | **78.5% (+23.6pp)** | 大幅改善 |
| 執行時間 | 2,693s | **1,976s (1.36× 快)** | bonus |

### Haplotag + ISM 層面結果（❌ 發現新問題）

```mermaid
graph TD
    subgraph VCF["VCF 層面 ✅"]
        V1["Somatic bias 消除"]
        V2["N50 翻倍"]
        V3["LOH.bed 不變"]
    end

    subgraph Haplotag["Haplotag 層面 ⚠️"]
        H1["Somatic 位點<br/>GT=0|0, GT2=.|1"]
        H2["LongPhase haplotag 解析:<br/>refHaplotype = UNDEFINED"]
        H3["所有 reads 統一<br/>標記為 HP:i:21"]
        H4["6,485 個非 LOH 平衡位點<br/>100% 出現此問題"]
    end

    subgraph ISM["ISM 結果 ❌"]
        I1["HP_Ratio TP median<br/>0.5000(Paired)<br/>0.8358(Baseline)<br/>0.0000(PON-only)"]
        I2["ISM-only LOH excess<br/>15.4% → 54.8%"]
    end

    VCF --> Haplotag
    Haplotag --> ISM
    H1 --> H2 --> H3 --> H4

    style VCF fill:#c8e6c9
    style Haplotag fill:#fff9c4
    style ISM fill:#ffcdd2
```

### 根因與修復方向

**根因**：LongPhase-TO haplotag 的 GT 解析邏輯——對 PON-only phased VCF 中 somatic 位點的 `GT=0|0`（homozygous REF），`refHaplotype` 變成 `UNDEFINED`，導致所有覆蓋該位點的 reads 被統一標記為某一 HP。

**最低成本修復方案**：ISM 端的 ReadParser 忽略 somatic HP tags（HP:i:11/21/33），只使用 germline HP:i:1/HP:i:2。

**長期方案**：修正 LongPhase-TO haplotag 模組的 GT 解析邏輯，或實作二次 phasing。

### LOH.bed Jaccard=1.0 的意外發現

**一個重要問題**：如果 self-phasing 是 TO LOH 的主因（62% 消失），為什麼移除 self-phasing 後 LOH.bed 完全不變（Jaccard=1.0）？

```mermaid
graph TB
    Q["Self-phasing 修正後<br/>LOH.bed 為何不變？"]

    H1["假說 A:<br/>LOH.bed 基於 VCF allele depth<br/>（而非 BAM HP tag）"]
    H2["假說 B:<br/>LOH.bed 的計算路徑<br/>不經過 site-level HP assignment"]

    I1["推論:<br/>ISM 的 HP_Ratio LOH<br/>與 LOH.bed 用不同 LOH 定義"]
    I2["證據支持:<br/>ISM LOH vs LOH.bed<br/>kappa = 0.670（不完全一致）"]

    V["⚠️ 待驗證:<br/>確認 LOH.bed 生成機制<br/>是 VCF AD 還是 BAM HP tag"]

    Q --> H1
    Q --> H2
    H1 --> I1
    H2 --> I1
    I1 --> I2
    I2 --> V

    style Q fill:#e3f2fd
    style V fill:#fff9c4
```

> **這是 P0 優先級的待驗證事項**——它影響我們如何理解兩套 LOH 判定系統的關係。

---

## 本章小結

| 問題 | 答案 |
|------|------|
| Self-phasing 是什麼？ | Somatic variants 自己 phase 自己，造成 HP 分配偏向 HP1 |
| 影響有多大？ | 94.6% somatic reads → HP1；62% **ISM-level** LOH 是 artifact |
| 影響了哪些特徵？ | 29 個直接依賴 HP 的特徵（佔 38%）不可信 |
| 不影響什麼？ | 42 個不依賴 HP 的特徵（佔 55%）結論穩固 |
| 修正到什麼程度了？ | VCF 層面完美修正；Haplotag 層有新問題待修 |
| 修正後 TO 能接近 Paired 嗎？ | ❌ 不能。Self-phasing 只解釋 ~35% 問題，65% 是結構性差距 |

> **重要精確化（來自方法論審查）**：「62% TO LOH 是 artifact」中的 LOH 指的是 **ISM 計算的 HP_Ratio LOH**（HP_Ratio < 0.1 or > 0.9），而非 **LOH.bed region-level LOH**。PON-only phasing 實驗已證實 LOH.bed 完全不受 self-phasing 影響（Jaccard=1.0）。因此，self-phasing 的因果影響位於 **haplotag → ISM HP_Ratio** 這條路徑上，而非 LongPhase 的 LOH region detection。兩套 LOH 系統使用不同定義（kappa=0.670 的不一致性由此解釋）。
