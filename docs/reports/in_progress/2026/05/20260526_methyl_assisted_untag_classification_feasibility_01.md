<!--
build_date: 2026-05-26
agent: claude-opus-4-7 (literature-review + source-code drill, 無實驗)
status: in_progress (純可行性論證，不含實測)
report_class: feasibility-landscape (探索性研究，巨觀地圖式)
audience: 用戶自己未來 / PI / lab member
scope: longphase-S paired-mode 下 untag + HP3 reads 用甲基化二次分類 + confidence scoring 的可行性
tier: ⭐2 L4 (literature + 既有實作 + 既有 evidence 串接；無 wet-lab/實測)
pipeline_track: paired_full (longphase-S, NOT longphase-to-mod)
parent_inventory: 7 主題 A-G (見 §2)
inputs:
  - /big8_disk/liaoyoyo2001/Knowledge/05_tools/longphase-s.md (KB v2.2)
  - /big8_disk/liaoyoyo2001/Knowledge/05_tools/methyl-somatic-analysis.md (KB v1.2)
  - /big8_disk/liaoyoyo2001/longphase-s/src/somatic_haplotag/SomaticHaplotagProcess.cpp (654 lines, primary source)
  - /big8_disk/liaoyoyo2001/MethylSomaticAnalysis/tools/python/advanced_analysis.py (clustering_analysis L469)
  - /big8_disk/liaoyoyo2001/MethylSomaticAnalysis/tools/python/read_phylogeny_validator.py (read × CpG matrix L41-49)
  - InterSubMod/src/core/ReadParser.cpp:121-154 (HP:Z + HP:i mapping)
  - InterSubMod/research/paired_priority_bug_audit/00_audit_report.md (chr19 tag 分布)
  - InterSubMod/research/paired_priority_bug_audit/01_step_D_germline_absent_finding.md (HP3 子類)
  - InterSubMod/docs/experiments/INDEX.md 多 entry (Phase BCD / LOH×AF×Methyl / LOH-constrained / Phase 2 cycle ablation)
  - 7 個 paired tagged BAM (HCC1395/HCC1395_DORADO/HCC1937/HCC1954/COLO829/H1437/H2009) on big7
outputs:
  - 本檔 (.md, 教學版完整稿)
  - companion standalone HTML (待 P5 turn 2 產出)
verdict: feasibility-landscape 多軸地圖 (見 §6)；非 binary GO/NO-GO
last_verified: 2026-05-26
report_template: feasibility-landscape v1.0 (narrative-frame Pyramid + Pre-Mortem hybrid)
-->

# 甲基化輔助 untag / HP3 二次分類與信心評分

## 可行性論證報告（純文獻 + 既有實作分析）

> **本報告性質**：探索性研究**起跑前的可行性論證**，不含實測。目的是把「**用甲基化二次分類 longphase-S 漏標 reads 並輸出 confidence score**」這個 idea 拆解成可驗證單元，盤點既有工具/資料/證據/反例後，給出**巨觀地圖式** verdict (Feasibility Landscape) 而非 binary GO/NO-GO。
>
> **教學定位**：困難概念配「直白比喻」，所有 source code 段附完整解讀，讓不熟此專案的讀者也能掌握。

---

## 0. TL;DR / Executive Summary

### 核心 idea（用戶 2026-05-26 提出）

對 longphase-S paired-mode tagged BAM，用 **甲基化 clustering** 做三件事：
1. **Q1**：把**沒被標 tag 的 reads (untag)** 重新歸類到 HP1 / HP2 / somatic
2. **Q2**：對 **HP3 (ambiguous somatic)** 細分為兩子類並給不同 confidence
   - 子類 A：read 沒經過任何 germline 位點
   - 子類 B：read 經過 germline 位點但 votes 不一致
3. **Q3**：對既有 tagged reads (HP1/HP2/HP1-1/HP2-1) 做**二次驗證**並輸出 confidence

### 關鍵盤點結論

| 項目 | 結論 |
|---|---|
| **資料齊全度** | 7 樣本 paired tagged BAM **全在 big7**（修正用戶原認知「only HCC1395」）|
| **HP3 二次分類 novelty** | ⚠️ **longphase-S 原生已有 `inheritHaplotype()` 做 HP3 confidence-based re-classification**（用 somatic-variant derived HP, threshold 0.6）— 用戶 idea 差異化在「**用 methylation 而非 somatic-variant 做 inheritance**」（orthogonal signal）|
| **Untag 二次分類 novelty** | ✅ longphase-S **完全不處理 untag** — 真正的 gap，用戶 idea 集中novelty 在這 |
| **MSA infrastructure** | ✅ read × CpG matrix + scipy linkage/DBSCAN/silhouette 已實作於 `read_phylogeny_validator.py` — 可重用 base |
| **Subclone framework 串接** | ✅ Phase BCD 4-group framework（Normal Diploid / Epi-Het / LOH / Tumor-Specific）可作為 ground truth |
| **重要 警示** | 🔴 Phase 2 LR ablation 證 methylation 5 features 在 multi-axis F1 filter 為 5th-rank vestigial（drop 後 ΔF1 只損 +0.0007）— 但**本場景是 single-axis classifier 非 multi-axis filter，場景不同不可直接 transfer**（用戶 ack 5/26）|

### Feasibility Landscape verdict（巨觀地圖）

8 cell × 4 scenario × 7 sample → **見 §6 完整 landscape**。簡言：

- 🟢 **C2 LOH sanity / C4 HP3 germline-absent / C6 tag 二次驗證**：方法應用最可能成立
- 🟡 **C1 untag non-LOH / C3 多倍體 / C5 HP3 inconsistent**：partial 成立，需實測
- 🔴 **Sa ASM region**：必失敗（germline confound），列為 control

---

## 1. 研究方向與動機 (Problem Framing)

### 1.1 longphase-S 的設計與限制

LongPhase-S（CCU Bioinformatics Lab 2025，Ho et al. bioRxiv）是 tumor-normal paired ONT 長讀序列的 somatic haplotagging 工具。它根據 germline phasing 結果，把 tumor BAM 內的每條 read 標上 `HP:Z:` tag。

> **直白比喻**：phasing 好比把房間裡所有家具分成「來自爸爸的傢俱」和「來自媽媽的傢俱」兩堆。longphase-S 對每個傢俱（read），看它身上的線索（SNV），判斷它屬於哪一堆。

可能的 tag 結果有 6 種：

| Tag | 説明 |
|---|---|
| `HP:Z:1` | Germline haplotype 1（純爸爸來源的 read） |
| `HP:Z:2` | Germline haplotype 2（純媽媽來源的 read） |
| `HP:Z:1-1` | 支持 somatic ALT，可追溯到 germline hap 1（爸爸 hap 上有體細胞突變的 read） |
| `HP:Z:2-1` | 支持 somatic ALT，可追溯到 germline hap 2 |
| `HP:Z:3` | 支持 somatic ALT 但**無法**由既有 germline haplotype 唯一導出（讀不清是哪邊的）|
| **no tag** | longphase-S 完全沒標 — 線索不足（沒經過任何 informative SNV，或 vote 不夠 confidence）|

### 1.2 觀察到的問題

實證 HCC1395 paired chr19（出處：`InterSubMod/research/paired_priority_bug_audit/00_audit_report.md` §3.1）：

```
chr19 primary reads total: 582,566
chr19 tagged (HP:Z: present): 354,919 (60.9%)
chr19 untagged:               227,647 (39.1%)  ← 這塊完全沒分類
```

| HP:Z: | reads | % of tagged |
|---|---:|---:|
| HP:Z:2 | 183,309 | 51.6% |
| HP:Z:1 | 143,760 | 40.5% |
| HP:Z:2-1 (somatic on HP2) | 14,504 | 4.1% |
| HP:Z:1-1 (somatic on HP1) | 12,401 | 3.5% |
| **HP:Z:3 (ambiguous)** | **1,145** | **0.3%** |

**問題顯化**：
- **Untag (39.1%, 227K reads)** 完全沒分類 — 下游任何依賴 HP tag 的分析都丟棄這塊
- **HP:Z:3 (0.3%, 1.1K reads)** 標了但「不確定 hap」，**沒有可量化的 confidence score** — 下游 binary 信不信
- **Tagged HP1/HP2/HP1-1/HP2-1** 也沒 confidence score — 都當作 100% 信
- **用戶要的**：3 件事一起做 = untag 救回 + HP3 細分 + tagged 二次評分

### 1.3 為什麼用甲基化？

甲基化（CpG methylation）作為 orthogonal signal 的三個理由：

1. **獨立於 SNV vote**：甲基化是 read 上 epigenetic 標記，longphase-S 完全不用這資訊，所以**訊號獨立**
2. **長度不依賴 informative SNV 數**：即使 read 沒經過任何 germline het，仍可能有 CpG methylation 資訊
3. **同 hap reads 甲基 pattern 應相似**：物理上同一 hap 上 reads 共享 epigenetic 環境（cell-of-origin 同 sub-population），cluster 後應抱團

> **直白比喻**：把 read 比作學生。SNV vote 像看學生的「姓」（爸爸或媽媽姓），有時看不清。甲基化像看學生的「穿著風格」— 同班同學（同 sub-clone）穿著風格相似，即使姓被遮住也能猜出班級。

### 1.4 假設可驗證的 3 個 Q

| Q | 假設 | 期望結果 |
|---|---|---|
| Q1 | Untag reads 的甲基 pattern 與**鄰近** tagged HP1/HP2 reads 相似可分群 | Untag 30-50% 可分到 HP1 或 HP2 |
| Q2 | HP3 兩子類（germline-absent vs germline-inconsistent）的**甲基 confidence 分布不同** | A 子類 confidence 較 monotonic；B 子類 bimodal |
| Q3 | 既有 tagged reads 的甲基化與 tag 一致率可作 confidence score | 甲基矛盾 reads 多為 false tag / FP variant |

### 1.5 為什麼這份是「純可行性論證」而非實測

用戶 2026-05-26 明示：「**純 literature review + 可行性論證**」 — 因為：
1. 一旦進實測，前期 hypothesis register + cycle init + pre-reg + benchmark + ledger 等流程開銷大
2. 先做可行性論證可避免「**做完才發現 longphase-S 原生已實作類似機制**」這類重複造輪子
3. 探索性研究在「**見林**」階段（用戶語）需保留多軸選擇空間，不應過早收斂為 GO/NO-GO

---

## 2. 既有實作與資料盤點 (Inventory)

> 本章 7 個主題（Theme A-G），是後續可行性論證的底層 evidence base。每個主題列「來源」「核心結論」「對用戶 idea 的 implication」。

### 2.1 Theme A: LongPhase-S HP:Z 5 編碼（KB authoritative source）

**來源**：`/big8_disk/liaoyoyo2001/Knowledge/05_tools/longphase-s.md` v2.2 (2026-04-17), KB id `kb-05-tools-longphase-s`

**核心結論**：5 tag + 1 untag = 6 個 read 分類（已列於 §1.1）

**KB 判讀邊界**（明文）：
> 「`HP:z:3` 是 **ambiguous somatic haplotype**，不是 LOH 標記」
>
> 「Purity estimation 主要依賴 germline haplotype imbalance，**不是用 `HP:z:3` 的數量直接判定 LOH 或 purity**」

**對用戶 idea 的 implication**：
- HP3 是 ambig 不是 LOH — 用戶 idea「HP3 細分」是正確方向
- 工具 binary 名稱差異重要：**longphase-S (paired)** ≠ **longphase-to (TO mode)** — 不同 codebase，HP tag 編碼也不同（HP:Z: 字串 vs HP:i: 整數）

### 2.2 Theme B: HCC1395 chr19 paired tag 實證分布

**來源**：`InterSubMod/research/paired_priority_bug_audit/00_audit_report.md` §3.1 (2026-05-09)，由 5/9 paired priority bug audit 跑出

**核心數字**：（已列於 §1.2）

**對用戶 idea 的 implication**：
- **untag : HP3 = 200 : 1**（227K : 1.1K）— untag 是主戰場，HP3 細分是 polishing
- chr19 tagged 60.9% 不算特別高 — 全基因組類似比例的話，untag 救回潛力大

### 2.3 Theme C: HP3 子類細分（germline-absent vs germline-inconsistent）

**來源**：`InterSubMod/research/paired_priority_bug_audit/01_step_D_germline_absent_finding.md` (2026-05-09)，5/9 paired audit Step D 補充

**核心結論**：InterSubMod 已有「germline-absent」概念實作於 **event-level**（不是 read-level）：
- 從 T1.2 vote dump 篩 `cnt_HP1 + cnt_HP2 = 0` (沒有 germline vote) AND `somatic_total > 0` (但有 somatic vote)
- chr19 germline-absent events = **5,789**

**Step D 的 5,789 events cross-tab**（baseline vs V3F vs V5）：

| paired HP:Z: | events | baseline hp=11 | baseline hp=21 | V3F hp=33 | V5 hp=11 | V5 hp=21 |
|---|---:|---:|---:|---:|---:|---:|
| HP:Z:1-1 | 2,040 | **1,679** | 318 | 2,040 | 1,679 | 318 |
| HP:Z:2-1 | 1,588 | **1,291** | 295 | 1,588 | 1,291 | 295 |
| HP:Z:3 | 530 | 342 | 178 | 530 | 343 | 177 |
| **加總** | — | **3,312** | **791** | **2,158** | 3,313 | 790 |

**關鍵 ratio**：
```
baseline germline-absent  HP1 : HP2 = 3,312 : 791 = 4.19 : 1   ← priority bug 偏移
V3F      germline-absent  hp=33 全部                           ← 保守不選邊
V5       germline-absent  HP1 : HP2 = 3,313 : 790 = 4.19 : 1   ← 與 baseline 同
```

**對用戶 idea 的 implication**：
- 用戶 idea「HP3 子類 A (germline-absent)」**已有概念實作**，但是 **event-level**（不是 read-level）
- Read-level aggregation 需新加：
  - 子類 A：read 跨的**所有** SNV events 全 germline-absent
  - 子類 B：read 跨的 events **有些** 有 germline vote 但 inconsistent
- 不需新 C++，只要 Python 後處理 vote dump

### 2.4 Theme D: MSA architecture (Level 1/2/3 + AdvancedAnalyzer)

**來源**：`/big8_disk/liaoyoyo2001/Knowledge/05_tools/methyl-somatic-analysis.md` v1.2 (2026-04-01)

**核心結構**：

| Level | 粒度 | 內容 |
|---|---|---|
| Level 1 | per-CpG × per-read | 最細：CpG 位置 / strand / read ID / methylation 狀態 / quality |
| Level 2 | per-variant window | 中：variant position / 平均 methylation / CpG count / coverage |
| Level 3 | per-haplotype-group | 粗：每 hap group 的 methylation 分佈統計 |

**Python AdvancedAnalyzer** (`/big8_disk/liaoyoyo2001/MethylSomaticAnalysis/tools/python/advanced_analysis.py`) 含 4 個工具：
- `wavelet_analysis` (L81)：多尺度甲基化訊號
- `kde_analysis` (L224)：核密度估計
- `gpr_analysis` (L319)：高斯過程回歸
- `clustering_analysis` (L469)：分群（**用戶 idea 直接對應**）

**對用戶 idea 的 implication**：
- ✅ Level 1 已輸出 per-CpG × per-read methylation — 用戶要的 read-level 資料**已有**
- ⚠️ `clustering_analysis` (L469) input 是 **feature_matrix（variant-level）** 而非 read-level — 不能直接套
- ✅ 但 KB 明示「建議使用 haplotagged BAM 配 MSA」— MSA 設計初衷就是配 longphase-S

### 2.5 Theme E: InterSubMod 5 條既有 paired-mode evidence

**來源**：`InterSubMod/docs/experiments/INDEX.md` + 8 個 memory entries

| 研究 | 結論 | tier | 對用戶 idea implication |
|---|---|---|---|
| **Phase BCD Dual BAM Validation** (INDEX:207, 2026-04-13) | 4-group subclone（Normal Diploid 17.5% / Epi-Het 12.9% / LOH 2.6% / Tumor-Specific 67%）; LOH BED 94.1% concordance with ISM hp_ratio | ⭐3-4 | **可作為 ground truth framework** — untag classify 後 cross-tab 對 4-group |
| **LOH × AF × Methylation Paired POSITIVE** (INDEX:210, 2026-04-14) | ΔNG=+0.787 (7/7 p<10⁻⁶⁵), median \|r\|=0.755 | ⭐3 | Paired 比 TO 效應更強，methylation × AF 真實存在 |
| **LOH-constrained phasing discovery (TO)** (INDEX:263, 2026-04-22) | NG=2 Inner ≥93% same-hap (6/6 sample) | ⭐3-4 | 物理約束限制 untag 真實 hap — LOH 區 untag 預期甲基單群 |
| **ASM PERMANOVA POSITIVE** (INDEX:133, 2026-04-01) | ASM 32-66% but FP/TP 重疊大 | ⭐ characterization only | ASM 訊號真實但 noisy — 反例風險源頭 |
| **Phase 2 cycle 3 ablation: ISM vestigial** (INDEX:243, 2026-05-20) | methylation 5 features drop 後 ΔF1 損 +0.00065 (3% of total uplift) | ⭐4 | 🔴 warning：methylation 在 multi-axis filter 為 vestigial — **但 paired classifier 與 TO filter 場景不同，不直接 transfer** |

### 2.6 Theme F: longphase-s `inheritHaplotype()` 原生 confidence 設計

**來源**：`/big8_disk/liaoyoyo2001/longphase-s/src/somatic_haplotag/SomaticHaplotagProcess.cpp:371-527` (primary source code)

**呼叫點** (L371-374)：
```cpp
if(hpResult == ReadHP::H3){
    // Inherit haplotype information for H3 reads based on somatic variant derived haplotype
    hpResult = inheritHaplotype(deriveByHpSimilarity, params.percentageThreshold, 
                                somaticVarDeriveHP, hpCount, hpResult, aln);
}
```

**`inheritHaplotype()` 完整邏輯** (L461-527)：
```cpp
int SomaticHaplotagChrProcessor::inheritHaplotype(
    float &deriveByHpSimilarity,
    double percentageThreshold,
    std::map<int, std::pair<int , int>>& somaticVarDeriveHP,
    std::map<int, int>& hpCount,
    int &hpResult,
    const bam1_t &aln
){
    int deriveByH1 = 0;
    int deriveByH2 = 0;
    
    // Step 1: 對 read 上所有 somatic variants，看 deriveHP 是 H1 還是 H2
    for(auto& somaticVarIter : somaticVarDeriveHP){
        int baseHP = somaticVarIter.second.first;
        int deriveHP = somaticVarIter.second.second;
        if(baseHP == SnpHP::SOMATIC_H3){
            if(deriveHP == SnpHP::GERMLINE_H1)
                deriveByH1++;
            else if(deriveHP == SnpHP::GERMLINE_H2)
                deriveByH2++;
        }
    }
    
    // Step 2: 算 confidence score = max / (max + min)
    int max = (deriveByH1 > deriveByH2) ? deriveByH1 : deriveByH2;
    int min = (deriveByH1 > deriveByH2) ? deriveByH2 : deriveByH1;
    int maxHp = (deriveByH1 > deriveByH2) ? SnpHP::GERMLINE_H1 : SnpHP::GERMLINE_H2;
    
    deriveByHpSimilarity = (max == 0) ? 0.0 : ((float)max / ((float)max + (float)min));
    
    // Step 3: 若 confidence 達 threshold（預設 0.6）→ 把 HP3 升級為 H1_1 或 H2_1
    if(deriveByHpSimilarity >= percentageThreshold){
        switch(maxHp){
            case SnpHP::GERMLINE_H1: hpResult = ReadHP::H1_1; break;
            case SnpHP::GERMLINE_H2: hpResult = ReadHP::H2_1; break;
        }
    }
    return hpResult;
}
```

**逐步教學解讀**：

> 1. **Step 1**：read 上每個 somatic variant 都有個 `deriveHP`（這變異是從 hap 1 還是 hap 2 衍生出來的）。對 read 全部 variants 各自加票數到 `deriveByH1` 或 `deriveByH2`。
>
> 2. **Step 2**：算「**多數派比例**」= max 票 / 總票。例如 8 票全部到 H1 = similarity 1.0；4 票 H1 + 4 票 H2 = similarity 0.5；6 票 H1 + 2 票 H2 = similarity 0.75。
>
> 3. **Step 3**：若 similarity ≥ 0.6（預設 `-p` 參數）→ 把 HP3 升級為 H1_1（多數派 H1）或 H2_1（多數派 H2）。否則保留 HP3。

**對用戶 idea 的 implication**：

| 用戶 Q | 對應 longphase-S 原生 | 用戶 idea 差異 |
|---|---|---|
| Q1 untag 分類 | **不處理**（保留 unTag）| ✅ 真正的 gap — 用戶 idea 真正 novelty 在這 |
| Q2 HP3 二次分類 | **已實作**（inheritHaplotype, 用 somatic-variant derive）| ⚠️ 原生**用 somatic-variant**, 用戶要**用 methylation**（orthogonal signal）— **新角度**但不是新點子 |
| Q3 confidence score | **已實作**（deriveByHpSimilarity ∈ [0,1]）但只對 HP3，不對 HP1/HP2/HP1-1/HP2-1 | ✅ 擴展到所有 tagged reads 是 gap |

**為什麼 methylation 是真正的 orthogonal signal？**

> longphase-S 的 inheritance 邏輯有個 **inherent weakness**：所有 somatic variants 在 self-phasing 機制下 100% 共現於同一 hap → vote 必然偏向一邊。
>
> 這正是 **V5 Layer 1.5 在 germline-absent 區繼承 priority bug 4.19:1 偏 HP1** 的根因（已實證於 `paired_priority_bug_audit/01_step_D`）。
>
> 甲基化訊號**不受** somatic phasing 共現影響 — 它是獨立記錄 read 自己 hap 的 epigenetic 狀態。所以 orthogonal。

### 2.7 Theme G: MSA read-level methylation matrix 既有實作

**來源 1**：`/big8_disk/liaoyoyo2001/MethylSomaticAnalysis/tools/python/read_phylogeny_validator.py:41-49`

```python
def build_read_methylation_matrix(self) -> pd.DataFrame:
    """構建read-CpG甲基化矩陣"""
    # 創建 read-CpG 位點矩陣
    pivot_data = self.data.pivot_table(
        index='read_id',
        columns='methyl_pos', 
        values='meth_call',
        fill_value=-1  # -1 表示無數據
    )
    print(f"   矩陣大小: {pivot_data.shape[0]} reads × {pivot_data.shape[1]} CpG位點")
    self.read_matrix = pivot_data
    return pivot_data
```

> **教學解讀**：把所有 reads 攤平成「**reads × CpG 位點**」的矩陣。行 = 每條 read（用 read_id 索引），列 = 每個 CpG 位置。值 = 0（未甲基）/ 1（甲基）/ -1（read 沒 cover 該 CpG）。這就是「read 條碼」資料結構。

**來源 2**：`read_phylogeny_validator.py` 含的 scipy clustering 函數：
- `linkage` (hierarchical clustering)
- `fcluster` (從 linkage 切出 cluster labels)
- `DBSCAN` (density-based)
- `silhouette_score` (cluster validity metric)

**對用戶 idea 的 implication**：
- ✅ **Read × CpG matrix 構建已有**（不用重寫）
- ✅ **Clustering tooling 已有**（linkage / DBSCAN / silhouette）
- ⚠️ **缺**：把 untag reads 「project to nearest HP1/HP2 cluster centroid」的分類器邏輯需新增（估 +200-300 行 Python，1-2 day 開發）
- ⚠️ **缺**：confidence score 量化（仿 longphase-S `deriveByHpSimilarity` 公式 max/(max+min)，但用甲基相似度而非 vote 數）

---

## 3. 觀察方法設計（P2 v2 完整版）

> **本章對應用戶 5/26 指示**：「LOH 區域與一般區域的分法、多倍體區域的狀況、分析合理性、能否統計觀察確認、對甲基訊號是否可信、相信程度、適用範圍 — 都是可以觀察確認並統計驗證的」

### 3.1 區域分類 3 軸（為什麼這樣分）

| 軸 | 分類 | 為什麼影響甲基訊號可信度 |
|---|---|---|
| **CN (Copy Number)** | CN=1 (LOH 單 hap) / CN=2 (diploid 雙 hap) / CN≥3 (gain 或 cnLOH 多 hap) | 物理上有幾條 hap 決定**甲基 cluster 應該分幾群** — 與 longphase-S 永遠 tag 2 hap 的設計天生有落差 |
| **ASM** | Known ASM (germline imprinting / cis-mQTL) / Non-ASM | ASM 區兩 hap methylation **本來就不同**（germline-origin），不是 somatic 訊號 → confound 必須排除 |
| **CpG density** | Dense (≥20 CpG per 2kb window) / Sparse (<20) | Cluster algorithm 需要足夠 features 才能分群 — sparse 區條碼太短會誤分 |

→ 3×2×2 = 12 combinations，實際只觀察 6-8 個 informative cells。

### 3.2 🛡️ 甲基訊號可信度評估 framework

**Trust Score** = 4 個 binary 條件加總，per 2kb window 評估：

| 條件 | 加 1 分若... | 理由 |
|---|---|---|
| CpG density | dense (≥20 CpG) | 條碼長度足夠不撞名 |
| Coverage | high (>30x) | 多次測量降低 noise |
| ASM | non-ASM (germline VCF + normal BAM 無 imprinting signature) | 排除 germline-origin confound |
| CN | CN=2 diploid (兩 hap 物理上對稱) | 期待 cluster 平衡 |

| Trust 等級 | 分數 | 訊號可信度 | 適用範圍 |
|---|---|---|---|
| 🟢 **HIGH** | 4/4 | 甲基 cluster 反映真實 hap 結構 | **核心可推論區**，可做 untag classification |
| 🟡 **MID** | 2-3/4 | 部分可信，須加 caveat | 觀察用，結論加註限制條件 |
| 🔴 **LOW** | 0-1/4 | 不可信，cluster = noise | **排除**，不在 inference 之列 |

### 3.3 6 個 main cells + 2 個 sanity cells

#### Main cells（對應用戶 Q1-Q3）

| Cell | Read 類別 | CN | ASM | Trust | 對應 Q | 預期結果（假設成立） |
|---|---|---|---|---|---|---|
| **C1** | Untag | CN=2 | Non-ASM | 🟢 HIGH | Q1 | 甲基 cluster 對齊鄰近 tagged HP1/HP2 centroid；untag recall 30-50% |
| **C2** | Untag | **CN=1 LOH** | Non-ASM | 🟢 HIGH | Q1 + sanity | 甲基**單群**（物理 1 hap）；若分 2 群 → 揭露 **epi-mosaic subclone** signal（bonus） |
| **C3** | Untag | **CN≥3 多倍體** | Non-ASM | 🟢 HIGH | Q1 | 甲基**可能 ≥2 群**（subclonal hap）；longphase-S 仍標 HP1/HP2 — 此區甲基 cluster **比 longphase-S 更精細** |
| **C4** | HP3-A (germline-absent) | CN=2 | Non-ASM | 🟢 HIGH | Q2 | 甲基推導 confidence > longphase-S 原生（後者在 germline-absent 受 priority bug 4.19:1 偏移） |
| **C5** | HP3-B (germline-inconsistent) | CN=2 | Non-ASM | 🟢 HIGH | Q2 | 甲基解決 germline vote conflict；confidence 高 reads promote 為 H1_1/H2_1 |
| **C6** | Tagged HP1-1/HP2-1 | CN=2 | Non-ASM | 🟢 HIGH | Q3 | 甲基與 tag 一致率作 confidence；矛盾 reads 標 low-confidence |

#### Sanity cells（**failure detection**，故意挑訊號不可信區）

| Cell | Read 類別 | 條件 | 為什麼是 sanity |
|---|---|---|---|
| **Sa** | Tagged HP1 | CN=2, **Known ASM** (imprinting) | ASM 訊號是 **germline-origin**；若 cluster 仍信任 → false positive somatic 歸類 |
| **Sb** | Tagged HP1 | **CN=1 LOH** | LOH 只 1 hap；若 cluster 仍分 2 群 → 確證 noise / overfit / contamination |

### 3.4 觀察與統計指標（per cell）

| Metric | 算法 | 用途 |
|---|---|---|
| **N reads in cell** | count | sample size 評估（CpG dense + CN=2 + Non-ASM 區可能不多）|
| **Trust score 分布** | 4 條件加總，per 2kb window | 排除 LOW |
| **Silhouette score** | scikit-learn `silhouette_score` | cluster 真實性（>0.3 視為真實 cluster）|
| **Confidence distribution** | `methyl_similarity ∈ [0,1]`，仿 longphase-S `deriveByHpSimilarity = max/(max+min)` 但用甲基相似度 | per-read confidence |
| **Classification accuracy** | held-out tagged HP1/HP2 reads 用同算法 re-classify | 監督式精度驗證 |
| **Failure rate** | silhouette < 0.1 比例 | unphasable 物理區比例 |
| **ASM confound rate** | germline VCF + normal BAM 找 imprinting biased CpG 比例 | false positive risk 量化 |
| **CN-cluster 落差** | C3 區 methylation cluster 數 vs longphase-S tag 數 | polyploid 區精細度差距 |

### 3.5 跨樣本重複設計（subclone framework 串接）

對 **7 樣本**（HCC1395 / HCC1395_DORADO / HCC1937 / HCC1954 / COLO829 / H1437 / H2009）跑同 8 cells 觀察。

統計：
- Wilcoxon paired test 跨樣本 metric 一致性
- 對齊 **Phase BCD 4-group subclone framework**（Normal Diploid 17.5% / Epi-Het 12.9% / LOH 2.6% / Tumor-Specific 67%）

→ 若 untag classification 後比例 **enrich 同方向** = 串接 subclone framework 成功。

### 3.6 教學版直白比喻

#### 為什麼 LOH 區（C2/Sb）是 sanity？

> LOH 區好比一個房間只放了 **1 條毛巾**（一條 hap）。你拍照分析時不可能看到 2 條不同顏色的毛巾。如果甲基 cluster algorithm 在這分出 2 群 — **一定是雜訊** 或 mosaic contamination。LOH 是天然 negative control，若這裡 fail，整套方法不可信。

#### 為什麼 CN≥3 (C3) 比 CN=2 (C1) 更需要看？

> CN=2 是普通 diploid（兩條 hap），longphase-S tag 為 HP1/HP2 兩群正好對應。
>
> CN=3 cnLOH 是「**一條 hap 複製成兩份 + 另一條 hap**」— 物理上仍 2 種 hap，但 reads 數量不均（3:1 而非 1:1）。
>
> CN=4 duplication 是「**各 hap 都複製**」— 物理仍 2 種 hap，但每 hap reads 都加倍。
>
> 真正的「3 群以上」要看 **subclonal 結構**：tumor 內存在 sub-population，每個 sub 有不同 hap state — 此時甲基 cluster 可能 ≥3 群，但 longphase-S 仍標 2 群 → **甲基比 longphase-S 更精細**。

#### 為什麼 ASM region 是 confound（Sa）？

> ASM = Allele-Specific Methylation = 兩條 hap 甲基狀態**天生不同**。3 個來源：
>
> - **來源 1: germline imprinting**（如 H19/IGF2 區，遺傳設定的）— 與 somatic 無關
> - **來源 2: cis-mQTL**（germline SNP 影響附近甲基）— 也是 germline-origin
> - **來源 3: subclonal ASM**（tumor 內 sub-population 各自 methylation）— **這才是我們要的 somatic signal**
>
> **區分方法**：cross-ref with **normal BAM**（normal 也有 ASM = germline 來源；只 tumor 才有 = somatic-only）

#### 為什麼甲基 cluster 比 longphase-S `inheritHaplotype` 更穩？

> longphase-S 的 inheritance 是「**用 read 上其他 somatic variants 各自 derive 出的 HP**」做 majority vote。但 self-phasing 機制下：所有 somatic variants 100% 共現於同 hap → vote 必然偏向一邊。
>
> **這正是 V5 Layer 1.5 在 germline-absent 區繼承 priority bug 4.19:1 偏 HP1 的根因**（paired_priority_bug_audit Step D 已實證）。
>
> 甲基訊號 **不受** somatic phasing 共現影響 — 它是獨立記錄 read 自己 hap 的 epigenetic 狀態。所以 orthogonal signal source。

#### 為什麼 CpG-dense 比 sparse 重要？

> 把每個 read 想成一個用 CpG methylation 0/1 串成的條碼。**條碼越長（CpG 越多），不同 hap 越容易分辨**。
>
> Read 只有 3 個 CpG = 只有 3 個字符的條碼，很容易撞名；
> 20 個 CpG = 20 字符條碼，幾乎不會撞名。
>
> 數學上：若 CpG 為獨立 50/50 機率，3 CpG 有 2^3 = 8 種可能條碼，1000 reads 必有大量碰撞；20 CpG 有 2^20 ≈ 100 萬種，1000 reads 幾乎不碰。

---

### 3.7 通用 Algorithmic Core（5 步驟 pseudocode）

> 本節給 8 cell 共用的算法核心。Per-cell 特殊化在 §3.8。

#### Step 1: Read × CpG Methylation Matrix Builder

```python
def build_read_methylation_matrix(reads_in_window, cpg_positions):
    """
    Args:
        reads_in_window: list of (read_id, hp_tag, [(cpg_pos, meth_call ∈ {0,1}), ...])
        cpg_positions: sorted list of unique CpG positions in window
    Returns:
        matrix: pd.DataFrame, index=read_id, columns=cpg_pos
                values ∈ {0=unmethylated, 1=methylated, -1=missing}
        meta: dict{read_id: {'hp_tag': ..., 'n_covered_cpg': int}}
    """
    matrix = pd.DataFrame(-1, index=[r[0] for r in reads_in_window], columns=cpg_positions, dtype=int)
    meta = {}
    for read_id, hp_tag, cpg_calls in reads_in_window:
        for pos, meth in cpg_calls:
            if pos in matrix.columns:
                matrix.loc[read_id, pos] = meth
        meta[read_id] = {'hp_tag': hp_tag, 'n_covered_cpg': (matrix.loc[read_id] != -1).sum()}
    # Drop reads with < min_cpg covered (default 5)
    valid_reads = [r for r, m in meta.items() if m['n_covered_cpg'] >= 5]
    return matrix.loc[valid_reads], {r: meta[r] for r in valid_reads}
```

**直白解讀**：把窗口內的 reads 攤成矩陣（行 = read, 列 = CpG 位置, 值 = 0/1/missing）；少於 5 個 CpG 的 read 視為「條碼太短」直接排除。

#### Step 2: Pairwise Distance (Jaccard with missing-data handling)

```python
def jaccard_distance_with_missing(r1, r2):
    """
    r1, r2: numpy arrays with values in {0, 1, -1}; -1 = missing
    Only compare positions where both have data (intersection).
    """
    mask = (r1 != -1) & (r2 != -1)
    n_shared = mask.sum()
    if n_shared < 3:  # too few shared CpG → unreliable
        return np.nan
    # Jaccard on methylated CpGs
    both_meth = ((r1[mask] == 1) & (r2[mask] == 1)).sum()
    either_meth = ((r1[mask] == 1) | (r2[mask] == 1)).sum()
    if either_meth == 0:
        return 0.0  # both unmethylated → identical pattern
    return 1.0 - (both_meth / either_meth)
```

**為什麼選 Jaccard 不選 Hamming？**
- Hamming 距離把 (0,0) 視為 match — 但兩 read 都 unmethylated 不代表同 hap（背景甲基化都低）
- Jaccard 只看 methylated CpG 的重疊比例 — 更聚焦 differential signal
- Missing data 用 intersection mask 處理（避免 false neighbors）

#### Step 3: Hierarchical Clustering（Ward linkage）

```python
from scipy.spatial.distance import squareform
from scipy.cluster.hierarchy import linkage, fcluster
from sklearn.metrics import pairwise_distances

# Compute pairwise distance matrix
def cluster_reads(read_matrix, n_clusters=2):
    D = pairwise_distances(
        read_matrix.values, 
        metric=jaccard_distance_with_missing
    )
    # Replace NaN (insufficient shared CpG) with median distance
    D = np.nan_to_num(D, nan=np.nanmedian(D))
    np.fill_diagonal(D, 0)
    
    # Ward linkage on condensed distance matrix
    Z = linkage(squareform(D, checks=False), method='ward')
    
    # Cut to n_clusters
    labels = fcluster(Z, t=n_clusters, criterion='maxclust')
    return labels, D
```

**為什麼選 Ward linkage 不選 average / single？**
- Ward 最小化 within-cluster variance — 對 binary 條碼資料較穩
- Single linkage 易產生 chaining (string-like cluster)
- Average linkage 對 outlier 敏感

**為什麼 n_clusters=2 是 default？**
- 對應 longphase-S HP1/HP2 兩群結構
- C3 多倍體 cell 可改 n_clusters=3 探索 subclonal hap

#### Step 4: Silhouette Score（cluster validity）

```python
from sklearn.metrics import silhouette_score

def cluster_quality(D, labels):
    """
    D: precomputed distance matrix (n × n)
    labels: cluster assignments
    Returns:
        score ∈ [-1, 1]
            > 0.5:  strong cluster structure
            0.3-0.5: reasonable
            0.1-0.3: weak
            < 0.1:  no real cluster (noise)
    """
    if len(set(labels)) < 2:
        return -1  # only one cluster, undefined
    return silhouette_score(D, labels, metric='precomputed')
```

**閾值判斷**：
| Silhouette | 解讀 | Action |
|---|---|---|
| > 0.5 | 強 cluster | 高信心 classification |
| 0.3 - 0.5 | 合理 cluster | 中信心，含 caveat |
| 0.1 - 0.3 | 弱 cluster | 低信心，標 mark |
| < 0.1 | noise / no structure | 排除該 window inference |

#### Step 5: Confidence Score（per-read assignment + similarity-based confidence）

```python
def assign_untagged_read(untag_read_profile, cluster_centroids):
    """
    Project an untagged read to its nearest cluster centroid.
    Confidence = max_similarity / (sum of all similarities)
    
    Args:
        untag_read_profile: 1D numpy array (binary methylation + missing)
        cluster_centroids: dict{cluster_id: mean methylation profile}
    Returns:
        (assigned_cluster, confidence ∈ [0.5, 1.0])
    """
    similarities = {}
    for cid, centroid in cluster_centroids.items():
        # Cosine similarity with missing-aware mask
        mask = (untag_read_profile != -1)
        if mask.sum() < 3:
            return (None, 0.0)  # too few CpG to classify
        u = untag_read_profile[mask].astype(float)
        c = centroid[mask].astype(float)
        # Cosine similarity
        if np.linalg.norm(u) == 0 or np.linalg.norm(c) == 0:
            similarities[cid] = 0.0
        else:
            similarities[cid] = np.dot(u, c) / (np.linalg.norm(u) * np.linalg.norm(c))
    
    if not similarities or max(similarities.values()) <= 0:
        return (None, 0.0)
    
    max_cid = max(similarities, key=similarities.get)
    # Confidence: longphase-S analog = max / sum (in 2-cluster case, ∈ [0.5, 1.0])
    sum_sim = sum(max(s, 0) for s in similarities.values())
    confidence = similarities[max_cid] / sum_sim if sum_sim > 0 else 0.0
    return (max_cid, confidence)
```

**Confidence threshold**（仿 longphase-S 0.6）：
- `confidence ≥ 0.6` → assign 給 max-cluster
- `confidence < 0.6` → 保留 untag（誠實 unphasable）
- 用戶可調此 threshold 做 sensitivity analysis

**直白解讀**：把 untag read 與每群 centroid 算 cosine 相似度，最像的一群得票最高，confidence = 領先程度 / 總分。完全壓倒 = 1.0（最像 HP1）；五五波 = 0.5（無法分）。

---

### 3.8 Per-Cell Algorithmic Specification（8 cells 詳細 spec）

> 每 cell 定義：(1) 篩 read 條件 (2) 特殊算法 (3) primary metric (4) pass criterion (5) 統計檢定。

#### C1: Untag, CN=2, Non-ASM, Trust HIGH（主軸 Q1 untag recall）

```yaml
filter:
  hp_tag: == 'untag'  # longphase-S 未標
  cn: == 2           # diploid 區
  asm_flag: False    # cross-ref normal BAM + germline VCF imprinting
  trust_score: == 4  # full HIGH
  window_cpg_density: >= 20 per 2kb
  window_coverage: > 30x

algorithm:
  1. Build read × CpG matrix for window (含 untag + nearby tagged HP1/HP2 as references)
  2. Cluster reference reads only (HP1/HP2 tagged) → derive centroid_HP1, centroid_HP2
  3. For each untag read: assign_untagged_read() → (cluster, confidence)
  4. Apply confidence threshold 0.6 → keep classified untag

primary_metrics:
  - rescue_rate = N_classified_untag / N_total_untag
  - silhouette_reference = silhouette(reference HP1/HP2)
  - centroid_separation = distance(centroid_HP1, centroid_HP2)

pass_criteria:
  - silhouette_reference > 0.3 (cluster real)
  - rescue_rate >= 0.30 (within trust HIGH subset)
  - held-out tagged HP1/HP2 re-classify accuracy >= 0.85

statistical_tests:
  - Wilcoxon paired test cross-sample rescue_rate
  - Bootstrap CI (n=1000) on confidence distribution
```

#### C2: Untag, CN=1 LOH, Non-ASM, Trust HIGH（sanity Q1 LOH）

```yaml
filter:
  hp_tag: == 'untag'
  cn: == 1  # LOH 區
  asm_flag: False
  trust_score: == 4

algorithm:
  1. Build read × CpG matrix
  2. Cluster with n_clusters=2 (forced)
  3. Compute silhouette
  4. Inverse check: n_clusters=1 alternative log-likelihood

primary_metrics:
  - silhouette(n=2): 若 > 0.3 → sanity FAIL (應 < 0.1)
  - bimodality_coefficient: 若 > 0.55 → epi-mosaic candidate
  - n_clusters_optimal: from gap statistic; 應 == 1

pass_criteria:
  - silhouette(n=2) < 0.1 (no fake 2-group) OR
  - bimodality_coefficient > 0.55 (epi-mosaic bonus finding)

statistical_tests:
  - Hartigan's dip test for unimodality (p > 0.05 = unimodal expected)
  - Gap statistic to determine optimal cluster count
```

#### C3: Untag, CN≥3 多倍體, Non-ASM, Trust HIGH（探索 Q1 polyploid）

```yaml
filter:
  hp_tag: == 'untag'
  cn: >= 3  # gain or cnLOH
  asm_flag: False
  trust_score: == 4
  subtype: separate {cnLOH=3 / gain=3 / cnLOH=4+ / gain=4+}

algorithm:
  1. Build read × CpG matrix
  2. Cluster with n_clusters ∈ {2, 3, 4} (gap statistic to select)
  3. Compare longphase-S tag (HP1/HP2 only 2 groups) vs methyl cluster
  4. If methyl cluster > 2 → potential subclonal hap state

primary_metrics:
  - optimal_n_clusters (gap statistic)
  - silhouette per n_clusters
  - subclonal_excess = optimal_n_clusters - 2 (>0 = methyl > longphase-S 精細)

pass_criteria:
  - silhouette(optimal_n) > 0.3
  - For subclonal_excess > 0: cross-ref Phase BCD 4-group "Tumor-Specific" (67%) enrichment

statistical_tests:
  - Permutation test: shuffle read assignments 1000 times, check observed silhouette > 95th percentile
```

#### C4: HP3-A germline-absent, CN=2, Non-ASM, Trust HIGH（Q2 子類 A）

```yaml
filter:
  hp_tag: == 'HP:Z:3'
  germline_absent: True  # cnt_HP1 + cnt_HP2 == 0 across all SNV events on read
  cn: == 2
  asm_flag: False
  trust_score: == 4

algorithm:
  1. Compute methyl confidence (Step 5 上述) for each HP3-A read
  2. Compute longphase-S native deriveByHpSimilarity (from somatic variant derive)
  3. Cross-compare: methyl_confidence vs native_similarity
  4. Disagreement reads: methyl assigns H1_1 but native assigns H2_1 (or vice versa)

primary_metrics:
  - methyl_confidence_distribution
  - native_similarity_distribution
  - agreement_rate (% reads same assignment both methods)
  - methyl_advantage = (% reads where methyl confidence > threshold but native < threshold)

pass_criteria:
  - methyl_advantage > 10% (methylation rescues 10% reads native cannot classify)
  - agreement_rate > 60% on classifiable reads

statistical_tests:
  - McNemar test on agreement/disagreement table
  - Compare priority-bug 4.19:1 偏 HP1 ratio: methyl_HP1:HP2 should be more balanced
```

#### C5: HP3-B germline-inconsistent, CN=2, Non-ASM, Trust HIGH（Q2 子類 B）

```yaml
filter:
  hp_tag: == 'HP:Z:3'
  germline_inconsistent: True  # cnt_HP1 > 0 AND cnt_HP2 > 0 on same read
  cn: == 2
  asm_flag: False
  trust_score: == 4

algorithm:
  1. Identify vote conflict severity = min(cnt_HP1, cnt_HP2) / max(cnt_HP1, cnt_HP2)
     (1.0 = tied, 0.0 = one-sided which shouldn't be classified as B)
  2. For severe conflict (severity > 0.5), compute methyl confidence
  3. Methylation 應「打破平局」 → expect bimodal confidence (high or low, few mid)

primary_metrics:
  - confidence_bimodality_coefficient
  - resolved_rate = % reads where methyl confidence >= 0.6
  - sequencing_error_rate (per-CpG quality score below 20)

pass_criteria:
  - resolved_rate > 40%
  - bimodality_coefficient > 0.5 (true tie-breaker pattern)
  - sequencing_error_rate < 15% (low conflict-from-noise)

statistical_tests:
  - Compare methyl-resolved vs random assignment baseline
```

#### C6: Tagged HP1-1/HP2-1, CN=2, Non-ASM, Trust HIGH（Q3 tagged 二次驗證）

```yaml
filter:
  hp_tag: ∈ {'HP:Z:1-1', 'HP:Z:2-1'}
  cn: == 2
  asm_flag: False
  trust_score: == 4
  region: ⊂ {known FP-rich zone OR matched random control zone}

algorithm:
  1. Use methyl cluster from tagged HP1/HP2 reads as ground truth centroids
  2. Project each HP1-1/HP2-1 read to centroids
  3. Compute methyl-tag agreement:
     - methyl_HP1 vs longphase-S HP1-1: agreement
     - methyl_HP1 vs longphase-S HP2-1: disagreement (→ low-confidence flag)
  4. Compare disagreement rate FP-rich zone vs control zone

primary_metrics:
  - agreement_rate per zone
  - disagreement_fp_enrichment = disagree_in_fp_zone / disagree_in_control_zone
  - confidence distribution

pass_criteria:
  - agreement_rate > 0.85 in trust HIGH cells
  - disagreement_fp_enrichment > 2.0 (disagreement reads enriched for FP variants)

statistical_tests:
  - Fisher exact test on (disagreement, FP/TP) 2x2
  - ROC curve for using confidence to predict FP variant
```

#### Sa: Tagged HP1, CN=2, Known ASM（confound stress test）

```yaml
filter:
  hp_tag: == 'HP:Z:1'
  cn: == 2
  asm_flag: True  # known ASM region (germline imprinting from public DB e.g. GeneImprint, or cis-mQTL Hi-C / GTEx eQTL)
  
algorithm:
  1. Run same clustering as C1 but on ASM region
  2. Compute hap-cluster correlation with germline VCF imprinting expectation
  3. Cross-ref normal BAM: if normal also shows 2-cluster ASM → germline origin

primary_metrics:
  - tumor_silhouette
  - normal_silhouette (must compute on matched normal BAM)
  - germline_imprinting_concordance = tumor_cluster matches normal_cluster

pass_criteria (sanity test, FAIL is expected behavior):
  - germline_imprinting_concordance > 0.7 → confirm ASM is germline-origin (sanity OK)
  - if confounded as somatic (concordance < 0.3) → false positive risk identified → exclude all ASM regions in main analysis

statistical_tests:
  - Chi-square test on tumor-normal cluster concordance
```

#### Sb: Tagged HP1, CN=1 LOH（overfit stress test）

```yaml
filter:
  hp_tag: == 'HP:Z:1'
  cn: == 1  # LOH 區 only 1 hap physically

algorithm:
  1. Build read × CpG matrix using only tagged HP1 reads (should be ~all reads in LOH)
  2. Force n_clusters=2 (test if algorithm splits noise)
  3. Compute silhouette; compute bimodality

primary_metrics:
  - forced_2cluster_silhouette
  - bimodality_coefficient
  - n_optimal_clusters (gap statistic, should be 1)

pass_criteria (sanity test, PASS = single cluster):
  - n_optimal == 1 (gap statistic favors 1 cluster)
  - forced_2cluster_silhouette < 0.1 (no fake split)
  - OR bimodality > 0.55 → epi-mosaic bonus

failure_implication:
  - If forced_2cluster_silhouette > 0.3 with bimodality < 0.55 → algorithm overfit noise → method NO-GO
```

---

### 3.9 Threshold Sensitivity Analysis（必跑 sub-experiment）

對下列 5 個 threshold 各跑 ±20% sensitivity，確認結論不依賴 threshold 任意選擇：

| Threshold | Default | Sensitivity range |
|---|---|---|
| `confidence_threshold` (assign untag) | 0.6 | {0.5, 0.6, 0.7, 0.8} |
| `min_cpg_per_read` (drop short reads) | 5 | {3, 5, 10, 15} |
| `silhouette_real_cluster` | 0.3 | {0.2, 0.3, 0.4, 0.5} |
| `window_size` (matrix scope) | 2kb | {1kb, 2kb, 5kb} |
| `cpg_density_dense_min` (trust scoring) | 20 per 2kb | {10, 20, 30} |

→ 每 cell 在每 threshold 跑一次，繪 heatmap 看 metric 對 threshold 的 robustness。

**Robust 判斷**：
- 若 main metric (e.g. C1 rescue_rate) 在 threshold ±20% 內變化 < 25% → robust
- 若變化 > 50% → fragile，需報告 threshold caveat
- 若變化 100%+ → algorithm 過度依賴 threshold，可能 overfit

---

## 4. 可行性論證 — Feasibility Landscape（巨觀見林）

> **本章對應用戶 5/26 指示**：「巨觀的觀察，全面到可以見林的程度」 — 不用 binary GO/NO-GO 收斂；改用**多軸 outlook landscape** 列出所有可能性。

### 4.1 8 cells × 4 scenario outcome（per-cell 可能結果地圖）

對每個 cell 列出 4 種可能結局，這樣讀者能看到完整可能性而非單一預期。

#### Cell C1 (Untag, CN=2, Non-ASM, Trust HIGH) — Q1 主軸

| Scenario | Outcome | 機率（先驗）| 對方法整體 implication |
|---|---|---:|---|
| **S1.1 POSITIVE** | 甲基 cluster 對齊鄰近 HP1/HP2 centroid; recall 30-50%, silhouette ≥0.3 | 40% | 方法核心可行；untag 救回價值確認 |
| **S1.2 PARTIAL** | recall 10-30%, silhouette 0.1-0.3 | 30% | 方法部分可用，需提高 trust threshold |
| **S1.3 NEGATIVE** | recall <10%, silhouette <0.1 | 20% | 物理不可分 dominant，untag 救回 ROI 低 |
| **S1.4 INCONCLUSIVE** | CN=2 + Non-ASM + Trust HIGH 區 cells 太少 (N < 100)，統計 power 不足 | 10% | 樣本範圍需擴展 |

#### Cell C2 (Untag, CN=1 LOH, Non-ASM, Trust HIGH) — Q1 sanity

| Scenario | Outcome | 機率 | implication |
|---|---|---:|---|
| **S2.1 POSITIVE sanity pass** | 甲基**單群** (silhouette near 0 or N_cluster=1)；確認 LOH 區方法行為一致 | 60% | 方法物理正確性確認；可信度提升 |
| **S2.2 BONUS finding** | 甲基**分 2 群**（但物理只 1 hap）→ 揭露 epi-mosaic / CTC contamination / subclone | 20% | 出乎意料的新發現，需追深 |
| **S2.3 sanity fail** | 甲基**分 2 群** + 對齊鄰近 tagged reads → 證明 algorithm overfit 在純 noise 上分群 | 15% | 整套方法**不可信**，NO-GO |
| **S2.4 INCONCLUSIVE** | LOH 區 untag reads N < 50，無 statistical power | 5% | 樣本擴展或 binary cluster threshold 調整 |

#### Cell C3 (Untag, CN≥3 多倍體, Non-ASM, Trust HIGH) — Q1 polyploid

| Scenario | Outcome | 機率 | implication |
|---|---|---:|---|
| **S3.1 EXPECTED** | 甲基 cluster 數 = 2 (對齊 longphase-S HP1/HP2)；untag recall 與 C1 類似 | 35% | 方法在 polyploid 區與 diploid 行為一致 |
| **S3.2 SUBCLONE BONUS** | 甲基 cluster ≥3 群（揭露 subclonal hap state）→ **甲基比 longphase-S 更精細** | 30% | 方法**超越** longphase-S — 新角度 |
| **S3.3 CONFUSED** | 甲基 cluster 數量不穩定（DBSCAN epsilon sensitivity）| 25% | 算法參數待 tune；先 inconclusive |
| **S3.4 LOW SIGNAL** | silhouette <0.1，分不出 | 10% | polyploid 區甲基訊號真的弱 |

#### Cell C4 (HP3-A germline-absent, CN=2, Non-ASM, Trust HIGH) — Q2 子類 A

| Scenario | Outcome | 機率 | implication |
|---|---|---:|---|
| **S4.1 POSITIVE** | 甲基推導 confidence > longphase-S inheritHaplotype（後者在 germline-absent 受 priority bug 偏 HP1 4.19:1 影響）| 50% | 用戶 idea Q2 對 A 子類成立；orthogonal signal value 確認 |
| **S4.2 PARTIAL** | 甲基 confidence 與 longphase-S 接近，無顯著優勢 | 30% | orthogonal value 微弱；HP3-A 兩 method 均可 |
| **S4.3 WORSE** | 甲基 confidence 比 longphase-S 還差（noisy）| 15% | HP3-A 內 CpG 太稀少 |
| **S4.4 INCONCLUSIVE** | HP3-A reads N < 30，無 power | 5% | 需擴 sample 範圍 |

#### Cell C5 (HP3-B germline-inconsistent, CN=2, Non-ASM, Trust HIGH) — Q2 子類 B

| Scenario | Outcome | 機率 | implication |
|---|---|---:|---|
| **S5.1 POSITIVE** | 甲基解決 vote conflict; confidence 高 reads promote 為 H1_1/H2_1 | 40% | Q2 對 B 子類成立 |
| **S5.2 NOISE** | conflict 來自 sequencing error 不是 subclone；甲基也 noisy | 30% | B 子類大部分為 noise，方法救不回 |
| **S5.3 MIXED** | 一半 reads 可救一半不行 | 25% | 需精細 threshold |
| **S5.4 INCONCLUSIVE** | HP3-B reads N < 30 | 5% | 樣本擴展 |

#### Cell C6 (Tagged HP1-1/HP2-1, CN=2, Non-ASM, Trust HIGH) — Q3 二次驗證

| Scenario | Outcome | 機率 | implication |
|---|---|---:|---|
| **S6.1 POSITIVE** | 甲基與 tag 一致率高（>85%），不一致 reads 多為 known FP variant | 45% | Q3 二次驗證機制成立；可用為 FP filter |
| **S6.2 BACKGROUND** | 甲基與 tag 一致率隨機（~50-60%），無 filtering 用途 | 30% | 方法不適用於 tagged 二次驗證 |
| **S6.3 CONFUSED** | 一致率與 tag 強烈 dependent on region（chr-specific）| 20% | per-region threshold needed |
| **S6.4 INCONCLUSIVE** | tagged HP1-1/HP2-1 在 FP-rich zone 太少 | 5% | 需擴 zone |

#### Sanity Cell Sa (Tagged HP1, CN=2, Known ASM) — confound stress test

| Scenario | Outcome | 機率 | implication |
|---|---|---:|---|
| **Sa.1 EXPECTED FAIL** | 甲基 cluster 強分 2 群（imprinting germline 訊號），如不排除 → false positive 把 germline ASM 當 somatic | 70% | 確認 ASM 必排除；驗證 Trust framework 必要 |
| **Sa.2 LUCKY** | imprinting 訊號弱於 cluster resolution → false positive 率低 | 20% | 可少擔心 ASM confound |
| **Sa.3 INCONCLUSIVE** | known ASM regions 在 paired BAM 內 read 不足 | 10% | 需 reference ASM panel |

#### Sanity Cell Sb (Tagged HP1, CN=1 LOH) — overfit stress test

| Scenario | Outcome | 機率 | implication |
|---|---|---:|---|
| **Sb.1 PASS** | 甲基**單群**（與物理 1 hap 一致）| 65% | algorithm 沒 overfit；可信 |
| **Sb.2 FAIL** | 甲基**分 2 群**（noise overfit）| 20% | algorithm 不穩定，NO-GO |
| **Sb.3 MOSAIC** | 甲基分 2 群但是 real mosaic / CTC contamination | 15% | bonus 發現，需 follow-up |

### 4.2 機制可信度光譜（mechanism trust spectrum）

從「最確定可用」到「最不確定」依序排：

```
🟢 確定可用 (Trust ≥ 0.9)
├─ 1. MSA read × CpG matrix 構建 (Theme G 已實作)
├─ 2. scipy linkage / DBSCAN / silhouette tooling (Theme G 已有)
├─ 3. longphase-S deriveByHpSimilarity 公式 max/(max+min) 可直接套用甲基相似度 (Theme F 證實設計)
├─ 4. paired tagged BAM 7 樣本可得 (P0b 驗證)
└─ 5. Phase BCD 4-group subclone framework 可作 ground truth (Theme E)

🟡 部分確定 (Trust 0.5-0.8)
├─ 6. Read-level CpG matrix 在 CpG-dense region 可分群 (理論支持, 需實測 silhouette)
├─ 7. Untag reads 30-50% 可被分類 (基於 chr19 untag 39.1% 估算, 假設 60-70% non-LOH non-ASM 為 trust HIGH)
├─ 8. HP3-A 子類 confidence > longphase-S 原生 (基於 V5 priority bug 4.19:1 證據, 但 methylation 在 germline-absent 區仍未實證)
└─ 9. CN≥3 區甲基 cluster ≥2 群 (理論預期, 但 subclonal hap state 未直接量化)

🟠 不確定 (Trust 0.3-0.5)
├─ 10. HP3-B 子類能解決 vote conflict (vote conflict 可能 sequencing error dominant)
├─ 11. Tagged HP1-1/HP2-1 二次驗證一致率 >85% (沒既有 baseline)
├─ 12. Cross-sample (7 samples) 一致 (paired 模式只 HCC1395 既有深度觀察)
└─ 13. Subclone classifications 對齊 Phase BCD 4-group framework (理論可接, 但 untag → subclone mapping 沒既有)

🔴 警示 (Trust < 0.3)
├─ 14. methylation 在 multi-axis filter 場景仍 vestigial (Phase 2 教訓; **本場景不同但仍需驗**)
├─ 15. HCC1937 BRCA1 outlier 可能破壞 cross-sample (已知 outlier sample)
└─ 16. 物理 unphasable untag (沒 CpG / read 太短 / MQ 低) 救不回 (需先量化 unphasable 比例)
```

### 4.3 與 longphase-S 既有 `inheritHaplotype()` 的 8 種 scenario 對照

把所有可能 read 結局列成 2×2×2 = 8 scenario matrix：

| # | longphase-S 處理 | methylation 介入 | 結果 |
|---|---|---|---|
| 1 | HP3 → longphase-S inherit 為 H1_1（confidence ≥0.6）| methylation 同意 H1_1 | **conflict-free**：信度提升至甲基+原生 average |
| 2 | HP3 → longphase-S inherit 為 H1_1 | methylation 投 H2_1 | **conflict**：保留 HP3 或標 low-confidence (per Cell C5 邏輯) |
| 3 | HP3 → longphase-S 留 HP3（confidence <0.6）| methylation 投 H1_1（confidence ≥0.6）| **methylation rescue**：promote 為 H1_1（**用戶 idea 直接 use case**）|
| 4 | HP3 → longphase-S 留 HP3 | methylation 也分不出 | 保留 HP3（誠實 unphasable）|
| 5 | untag → longphase-S **不處理** | methylation cluster 對齊 HP1 | **untag rescue**：assign HP1（**用戶 idea Q1 直接 use case**）|
| 6 | untag → longphase-S 不處理 | methylation cluster 對齊 HP2 | 同上，assign HP2 |
| 7 | untag → longphase-S 不處理 | methylation cluster 不對齊任一 hap | 保留 untag（誠實 unphasable）|
| 8 | tagged HP1/HP2/HP1-1/HP2-1 → longphase-S **不再評**| methylation 與 tag 矛盾 | 標 low-confidence；**Q3 use case** |

### 4.4 multi-sample × multi-cell 觀察 grid（理論期望表）

對 7 樣本 × 8 cells 跑同樣 metric，理論 expected matrix（**未實測，純理論預期**）：

| Sample / Cell | C1 | C2 | C3 | C4 | C5 | C6 | Sa | Sb |
|---|---|---|---|---|---|---|---|---|
| HCC1395 | 🟢 expected | 🟢 sanity pass | 🟢 expected | 🟢 expected | 🟡 mixed | 🟢 expected | 🔴 expected fail | 🟢 pass |
| HCC1395_DORADO | 🟢 expected | 🟢 sanity pass | 🟢 expected | 🟢 expected | 🟡 mixed | 🟢 expected | 🔴 expected fail | 🟢 pass |
| HCC1937 (BRCA1 outlier) | 🟡 mixed | 🟢 sanity pass | 🟡 sample-specific | 🟡 mixed | 🟡 mixed | 🟡 mixed | 🔴 expected fail | 🟢 pass |
| HCC1954 | 🟢 expected | 🟢 sanity pass | 🟢 expected | 🟢 expected | 🟡 mixed | 🟢 expected | 🔴 expected fail | 🟢 pass |
| COLO829 (ONT R10 無 methyl?) | ⚠️ method 不適用 | ⚠️ | ⚠️ | ⚠️ | ⚠️ | ⚠️ | ⚠️ | ⚠️ |
| H1437 | 🟢 expected | 🟢 sanity pass | 🟢 expected | 🟢 expected | 🟡 mixed | 🟢 expected | 🔴 expected fail | 🟢 pass |
| H2009 | 🟢 expected | 🟢 sanity pass | 🟢 expected | 🟢 expected | 🟡 mixed | 🟢 expected | 🔴 expected fail | 🟢 pass |

**讀法**：
- 🟢 **expected** = 方法應該成立
- 🟡 **mixed / partial** = 方法部分成立或 outlier 行為
- 🔴 **expected fail** = 方法在此 cell 預期失敗（sanity / confound stress test）
- ⚠️ **method 不適用** = COLO829 用 ONT R10 chemistry **可能無 methylation calling**（需 P0c 驗證 BAM 是否有 MM/ML tag）

**巨觀觀察**：
- 跨樣本一致性最強的 cells：C1 / C2 / C4 / C6（5 樣本 expected）
- 跨樣本一致性弱的 cells：C3 / C5（理論本身有 sample-specific 維度）
- COLO829 caveat：需先驗證是否有 methylation calling，若無則樣本範圍降為 6

---

## 5. Pre-mortem 反例與失敗模式（P4）

> **本章不是「肯定成立」框架，而是反向假設「方法已失敗，what could have caused this?」** — 列出 top 5 failure modes 與既有 evidence。

### FM1: ASM 訊號來自 germline imprinting / cis-mQTL（與 somatic 無關）

| 項目 | 內容 |
|---|---|
| 失敗機制 | C3 sanity cell Sa 內，甲基化兩 hap 差異本是 germline 設定，不是 somatic 訊號 → cluster 誤把 germline ASM 當 somatic subclone |
| 既有 evidence | `project_snv_methylation_association`: ASM 32-66% SNV 位點，FP/TP 重疊大 |
| 反例比例估算 | 32-66% of 甲基 differential CpG 為 germline-origin |
| Mitigation | (1) cross-ref with normal BAM ASM baseline (2) germline VCF imprinting region 排除 (3) 限定 cells C1/C4/C5/C6 在 non-ASM region |

### FM2: Untag 物理不可分（read 太短 / MQ 低 / coverage 不足，甲基也救不回）

| 項目 | 內容 |
|---|---|
| 失敗機制 | Untag 內可能 40-60% 為「物理 unphasable」 — read 沒經過任何 informative SNV、長度太短、MQ < 20 — 甲基化也面臨相同限制（CpG 太少）|
| 既有 evidence | longphase-S `qualityThreshold` 預設 1（很寬鬆）、`percentageThreshold` 預設 0.6 — 嚴格設可能更多 untag |
| 反例比例估算 | 40-60% untag 預期 trust LOW |
| Mitigation | (1) Trust framework 排除 LOW (2) 對 untag 內 trust HIGH subset 算 net recall (3) 報告 "rescue rate among phasable reads" 非 "rescue rate among all reads" |

### FM3: HP:Z:3 太稀少（chr19 0.3%）導致 cluster 信號 noise dominate

| 項目 | 內容 |
|---|---|
| 失敗機制 | chr19 HP3 只 1,145 reads (0.3%)，全基因組推估 < 100K；per 2kb window 多 <10 reads — 統計 power 不足 |
| 既有 evidence | Theme B 實證 chr19 數字 |
| 反例比例估算 | 預估全基因組 80-90% windows 內 HP3 reads < 10 |
| Mitigation | (1) aggregate per 1Mb window 提升 statistics (2) 跨 sample pool HP3 (3) HP3 細分降級為 secondary objective，主軸聚焦 untag |

### FM4: HCC1937 BRCA1 outlier 行為破壞 cross-sample generalize

| 項目 | 內容 |
|---|---|
| 失敗機制 | HCC1937 是 BRCA1 mutant + CNV-driven germline het sample，已知 FP/TP=0.194 vs 其他樣本 0.01（per `project_hpfinengroups_subclone_marker`）— 此 sample 在 untag classification 可能 outlier |
| 既有 evidence | 5/20 5-Goal validation Goal 2 HPFineN cross-sample HCC1937 outlier |
| 反例比例估算 | 7 樣本中 1/7 outlier |
| Mitigation | (1) per-sample 報告，不強行 pool (2) cross-sample summary 用 median 不用 mean (3) HCC1937 列為 sensitivity analysis 子項 |

### FM5: methylation cluster centroid 偏向 majority hap，untag 偏向 over-classified to HP1

| 項目 | 內容 |
|---|---|
| 失敗機制 | HCC1395 chr19 HP1:HP2 = 1:1.275 已小幅不平衡；若 hap-imbalance 區域更嚴重（LOH-adjacent），cluster centroid 會 dominantly 反映 majority hap → untag 偏向被歸類為 majority |
| 既有 evidence | paired_priority_bug_audit chr19 數字 |
| 反例比例估算 | per-window 預估 20-30% windows 有顯著 hap imbalance |
| Mitigation | (1) cluster centroid 用 balanced subsample 計算 (2) Z-score normalize per region (3) 報告 per-hap recall 分別不混合 |

### 止證項（falsification criteria）

若以下任一條件達成，方法**整體 NO-GO**：

| 條件 | NO-GO 判斷 |
|---|---|
| Sb (LOH sanity) sanity fail 比例 >30% (Scenario Sb.2) | algorithm overfit dominant，不可信 |
| C1 untag silhouette 中位數 < 0.1 跨 5+ 樣本 | 物理 unphasable dominant，untag 救不回 |
| Sa (ASM confound) false positive 率 > 50%（Scenario Sa.1）且 mitigation 無效 | germline confound 不可控 |
| C2 + Sb 跨 5+ 樣本 mosaic finding (Scenario Sb.3 / S2.2) 一致 → 改 pivot 為「epi-mosaic detector」研究 | 主目標轉向 |

---

## 6. 結論與探索建議 — 巨觀地圖

### 6.1 Feasibility Landscape Verdict（非 binary）

按用戶要求「巨觀見林」，verdict 以**地圖式**呈現：

```
                     高 ROI / 高 Novelty
                              │
                              │
        ┌─────────────────────┼─────────────────────┐
        │                     │                     │
        │   C4 (HP3-A)        │   C1 (untag CN=2)  │
        │   methyl orthogonal │   主軸救回 untag    │
        │   vs longphase-S    │   30-50% recall     │
        │   priority bug      │                     │
        │   ★ recommended     │   ★ recommended     │
        │                     │                     │
高──────┼─────────────────────┼─────────────────────┼────── 高
信       │                     │                     │   不確定性
心       │   C2 (LOH sanity)   │   C3 (polyploid)   │
        │   方法物理正確性     │   subclonal hap     │
        │   ✓ must include    │   △ exploratory     │
        │                     │                     │
        │                     │                     │
        ├─────────────────────┼─────────────────────┤
        │                     │                     │
        │   C6 (tag 二次驗證)  │   C5 (HP3-B)        │
        │   FP filter 潛力     │   confidence mixed  │
        │   △ partial value   │   △ partial         │
        │                     │                     │
        │   Sa (ASM control)  │   Sb (overfit ctrl) │
        │   ⛔ control         │   ⛔ control         │
        │                     │                     │
        └─────────────────────┴─────────────────────┘
                              │
                     低 ROI / 低 Novelty
```

### 6.2 探索建議（按 ROI 排序）

**Phase A（推薦先做，最低成本驗證 trust framework）**：
- Tooling MVP：移植 `read_phylogeny_validator.py` 的 read × CpG matrix builder + scipy linkage → 純 Python `untag_classifier.py`
- 在 **HCC1395 chr19 only** 跑 C1 + C2 + Sa + Sb 4 cells（避開 chr 數 statistical noise）
- 工程量 < 1 day，無需新 C++
- Pass criterion：C2 sanity pass + Sa confound 確認 + C1 silhouette > 0.1

**Phase B（Phase A 成功後）**：
- 擴 C3 (polyploid) + C4 + C5 + C6
- HCC1395 全基因組
- 工程量 ~2-3 day

**Phase C（Phase A+B 成功後）**：
- Cross-sample 7 樣本（先檢查 COLO829 是否有 methylation calling）
- Wilcoxon paired test 跨樣本一致性
- 串接 Phase BCD 4-group subclone framework
- 工程量 ~3-5 day

**Phase D（NOT recommend immediately）**：
- 整合進 MSA 或 longphase-S binary（C++ 改動）
- 待 Phase A-C 全 POSITIVE 後再啟動
- 工程量 ~2-3 週

### 6.3 主要不確定性 + 開放問題（巨觀視野）

| 不確定性 | 為什麼重要 | 解法路徑 |
|---|---|---|
| **物理 unphasable untag 真實比例** | 決定 untag 救回的 ceiling | Phase A 跑完看 C1 silhouette < 0.1 比例 |
| **HCC1937 在新方法下行為** | BRCA1 outlier 是否會破壞 cross-sample | Phase C 必含 HCC1937 sensitivity |
| **COLO829 methylation calling** | 是否減少樣本 N 從 7 → 6 | Phase 0 額外驗證 BAM 有無 MM/ML tag |
| **CN≥3 多倍體區的甲基 cluster 行為** | 是否真比 longphase-S 更精細 | Phase B C3 cell 核心觀察 |
| **與 Phase 2 LR vestigial 教訓的關係** | 場景不同但需驗證 single-axis 是否真 transferable | Phase B 加 LR baseline comparison |

### 6.4 推薦不做的方向（pre-emptive scope control）

| 不推薦 | 理由 |
|---|---|
| 直接整合 C++ binary | Phase A-C 未驗證前風險高（可能像 Phase 2 cycle 1+0.022 一樣後來 LOSO -0.00012） |
| HP3 細分為**主軸** | longphase-S 原生已做（雖然不用 methylation），novelty 集中在 untag 才合理 |
| 跑全 7 樣本 immediately | Phase A 單樣本 chr19 沒過則 7 樣本浪費；按階段推進 |
| 投稿 paper 前先試 | Phase 2 LR vestigial 教訓 — multi-axis filter context 已踩雷，類似機制 transferability 需 strong evidence |

---

## 7. Glossary / 教學專欄

> 給不熟此專案的讀者用，每個 term 一句直白比喻 + InterSubMod 既有引用。

| Term | 直白比喻 | InterSubMod 引用 |
|---|---|---|
| **Haplotype (hap)** | 來自爸爸 / 媽媽各一條的染色體 | Theme A KB |
| **Phasing** | 把每個變異/read 標屬於爸爸 hap 還是媽媽 hap | Theme A KB |
| **longphase-S** | 工具：tumor-normal paired ONT phasing + somatic tagging（CCU Lab）| `/big8_disk/.../longphase-s/` |
| **HP:Z: tag** | longphase-S 寫在 BAM read 上的 hap 標記（字串：1 / 2 / 1-1 / 2-1 / 3）| `SomaticHaplotagProcess.cpp:533` |
| **HP:Z:3 (ambig)** | 有 somatic 訊號但分不出 hap 的 read | Theme F |
| **untag** | longphase-S 完全沒標的 read（資訊不足）| Theme B 39.1% |
| **inheritHaplotype** | longphase-S 對 HP3 用 somatic-variant derived HP 做二次判斷的內部函數 | `SomaticHaplotagProcess.cpp:461-527` |
| **deriveByHpSimilarity** | longphase-S 的 confidence score = max-vote / total-vote ∈ [0,1] | 同上 L503 |
| **percentageThreshold** | longphase-S 的 confidence threshold，預設 0.6 | KB longphase-s.md `-p` 參數 |
| **LOH (Loss of Heterozygosity)** | 某 hap 在 tumor 中遺失了 — 該區域物理上只剩 1 hap | `project_loh_constrained_phasing_discovery` |
| **cnLOH** | LOH 但 CN 不變（另一 hap 複製補上）| `project_loh_research_status` |
| **CN (Copy Number)** | 該區域 hap 的份數（CN=1 缺一條、CN=2 正常、CN=3+ 多了）| |
| **ASM (Allele-Specific Methylation)** | 兩條 hap 甲基化天生不同的區域 | `project_snv_methylation_association` |
| **cis-mQTL** | germline SNP 影響附近甲基化的現象（germline-origin ASM）| Theme E |
| **Imprinting** | 父母來源決定基因/甲基狀態（如 H19/IGF2）| Theme E |
| **CpG density** | 區域內 CpG 位點密度（dense ≥20 per 2kb / sparse < 20）| Theme G |
| **Silhouette score** | cluster 真實性 metric（>0.3 視為真實，<0.1 為 noise）| scikit-learn |
| **MSA (MethylSomaticAnalysis)** | 工具：methylation × somatic variant 關聯分析（C++ + Python）| `/big8_disk/.../MethylSomaticAnalysis/` |
| **Phase BCD 4-group subclone** | Normal Diploid / Epi-Het / LOH / Tumor-Specific 4 群 framework | INDEX:207 |
| **Trust score** | 本報告新定義：4 條件加總 (CpG dense / Coverage high / Non-ASM / CN=2) 評估甲基訊號可信度 | §3.2 |
| **Priority bug** | longphase-to / V5 早期 binary 在 vector ordered check 上把 somatic vote 蓋過 germline vote 的 bug | `paired_priority_bug_audit/01_step_D` 4.19:1 偏 HP1 |

---

## 8. References / Source Provenance

### 8.1 KB sources（authoritative source-of-truth）

- `/big8_disk/liaoyoyo2001/Knowledge/05_tools/longphase-s.md` v2.2 (2026-04-17), kb-05-tools-longphase-s
- `/big8_disk/liaoyoyo2001/Knowledge/05_tools/methyl-somatic-analysis.md` v1.2 (2026-04-01), kb-05-tools-methylsomaticanalysis
- `/big8_disk/liaoyoyo2001/Knowledge/05_tools/longphase.md`
- `/big8_disk/liaoyoyo2001/Knowledge/05_tools/longphase-to.md`

### 8.2 Primary source code

- `/big8_disk/liaoyoyo2001/longphase-s/src/somatic_haplotag/SomaticHaplotagProcess.cpp:310-527` (judgeHaplotype + inheritHaplotype)
- `/big8_disk/liaoyoyo2001/longphase-s/src/somatic_haplotag/SomaticHaplotag.cpp`
- `/big8_disk/liaoyoyo2001/MethylSomaticAnalysis/tools/python/advanced_analysis.py:469` (clustering_analysis)
- `/big8_disk/liaoyoyo2001/MethylSomaticAnalysis/tools/python/read_phylogeny_validator.py:41-49` (build_read_methylation_matrix)
- `InterSubMod/src/core/ReadParser.cpp:121-154` (HP:Z + HP:i mapping)

### 8.3 InterSubMod previous reports

- `InterSubMod/research/paired_priority_bug_audit/00_audit_report.md` (chr19 tag 分布, 2026-05-09)
- `InterSubMod/research/paired_priority_bug_audit/01_step_D_germline_absent_finding.md` (HP3 子類, 2026-05-09)
- `InterSubMod/docs/experiments/INDEX.md` (5 條 paired-mode evidence)
- `InterSubMod/docs/experiments/in_progress/2026/05/20260518_Phase2_Cycle1_Global_FP_Filter_01.md` (Phase 2 LR vestigial 教訓)

### 8.4 Paired tagged BAM 7 樣本（已驗證齊全於 big7, 2026-05-26）

```
/big7_disk/liaoyoyo2001/big7_disk_output/canonical/<sample>/paired_full/<dated>_complete_matrix/longphase_s/<sample>_tagged.bam
```

| Sample | BAM 路徑（縮寫） | 大小 |
|---|---|---|
| HCC1395 | `.../20260314_..._complete_matrix/longphase_s/HCC1395_tagged.bam` | 260 GB |
| HCC1395_DORADO | `.../20260315_..._complete_matrix/longphase_s/HCC1395_DORADO_tagged.bam` | (TBD) |
| HCC1937 | `.../20260315_..._complete_matrix/longphase_s/HCC1937_tagged.bam` | (TBD) |
| HCC1954 | `.../20260315_..._complete_matrix/longphase_s/HCC1954_tagged.bam` | (TBD) |
| COLO829 | `.../20260315_..._complete_matrix/longphase_s/COLO829_tagged.bam` ⚠️ ONT R10 may lack methylation | (TBD) |
| H1437 | `.../20260315_..._complete_matrix/longphase_s/H1437_tagged.bam` | (TBD) |
| H2009 | `.../20260315_..._complete_matrix/longphase_s/H2009_tagged.bam` | (TBD) |

### 8.5 InterSubMod memory entries（cross-ref）

- `project_loh_subclone_af_methylation_positive.md` (Paired ΔNG=+0.787 7/7 p<10⁻⁶⁵)
- `project_loh_constrained_phasing_discovery.md` (NG=2 Inner ≥93% same-hap)
- `project_methyl_filter_pilot_marginal_positive.md` (Phase 2 ΔF1=+0.00242 marginal)
- `project_snv_methylation_association.md` (ASM 32-66%)
- `project_loh_research_status.md` (LOH = paired 16-52×)
- `project_hp_integer_tag_fix.md`
- `project_pon_only_phasing_verification.md`
- `project_v3_fixed_haplotag_verification.md`

### 8.6 Companion HTML

- `InterSubMod/docs/reports/in_progress/2026/05/20260526_methyl_assisted_untag_classification_feasibility_01.standalone.html` (待 P5 turn 2 產出)

---

**報告版本**：v1.0 (2026-05-26)
**作者**：Claude Opus 4.7
**狀態**：in_progress（純可行性論證；下一步啟動 Phase A MVP 需用戶 explicit go）
**Tier**：⭐2 L4（literature + 既有實作 + 既有 evidence 串接；無實測）
