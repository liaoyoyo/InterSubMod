<!--
建立時間: 2026-03-24
目標: InterSubMod 簡報用概念說明與未來方向，以投影片版面格式呈現
處理範圍: ISM 核心概念、技術流程、實驗結論、文獻定位、突破方向
關聯檔案:
  - docs/references/20260324_InterSubMod研究方法與相關文獻突破方向全域分析_01.md
  - docs/methodology/20260324_方法學審查全域結論報告_01.md
  - docs/architecture/deep-research-report.md
-->

# InterSubMod — 簡報用概念與方向說明

> 本文件以「投影片 + 備忘錄」格式撰寫：
> - **投影片區塊**（`[SLIDE]`）：僅含重點文字，適合直接放入 PPT
> - **備忘錄區塊**（`[NOTES]`）：詳細說明與背景，供講者參考

---

## ═══════════════════════════════════════════
## SLIDE 1 — 封面
## ═══════════════════════════════════════════

```
┌─────────────────────────────────────────────────────┐
│                                                     │
│              InterSubMod                            │
│    Inter-Subclonal Methylation Analysis             │
│                                                     │
│    ─────────────────────────────────────             │
│                                                     │
│    以長讀序體細胞變異為錨點                           │
│    整合 read-level 甲基化、單倍型、拷貝數資訊         │
│    提供表觀遺傳層級的變異解讀                          │
│                                                     │
│    ─────────────────────────────────────             │
│                                                     │
│    2026-03                                          │
│                                                     │
└─────────────────────────────────────────────────────┘
```

### [NOTES] 備忘錄

- InterSubMod = Inter-Subclonal Methylation Analysis Module
- 核心：利用 Oxford Nanopore（ONT）長讀序定序同時攜帶突變資訊與甲基化資訊的特性
- 目標：在每個體細胞突變（somatic SNV）附近建立甲基化結構分析，判讀其表觀遺傳意義
- 定位：不是 variant caller（突變偵測器），而是 variant interpreter（突變解讀器）

---

## ═══════════════════════════════════════════
## SLIDE 2 — 為什麼需要 InterSubMod？
## ═══════════════════════════════════════════

```
┌─────────────────────────────────────────────────────┐
│                                                     │
│  ▎ 問題：腫瘤基因體分析的三個盲點                      │
│                                                     │
│  ┌───────────────────────────────────────┐           │
│  │  1. 突變真偽難辨                        │           │
│  │     Somatic caller 的 FP 中有些帶有      │           │
│  │     真實的甲基化差異（germline ASM）       │           │
│  └───────────────────────────────────────┘           │
│  ┌───────────────────────────────────────┐           │
│  │  2. 等位特異性甲基化 (ASM)               │           │
│  │     哪個等位基因被甲基化？                 │           │
│  │     → 影響基因表現與 second hit 判讀       │           │
│  └───────────────────────────────────────┘           │
│  ┌───────────────────────────────────────┐           │
│  │  3. 亞克隆表觀結構                      │           │
│  │     腫瘤內部的甲基化異質性                 │           │
│  │     → 演化分支、治療抗性線索               │           │
│  └───────────────────────────────────────┘           │
│                                                     │
│  ◆ ISM 的解法：把同一條 read 上的突變+甲基化+           │
│    單倍型三種訊號串起來                                │
│                                                     │
└─────────────────────────────────────────────────────┘
```

### [NOTES] 備忘錄

**盲點 1 — 突變真偽**：
現有 somatic variant caller（如 ClairS、DeepSomatic）已能達到 F1~0.85-0.98，但殘留的假陽性中有一部分來自 germline allele-specific methylation（胚系等位特異性甲基化）。這些位點在甲基化層面確實有差異，但不是腫瘤驅動的。ISM 可以標記這些區域，輔助判讀。

**盲點 2 — ASM 與 second hit**：
腫瘤抑制基因（TSG）的雙等位失活（biallelic inactivation）是癌症關鍵機制。經典模式：一個等位基因發生突變，另一個等位基因被甲基化沉默（promoter methylation + LOH）。長讀序可以在同一分子上同時看到突變和甲基化狀態，直接判斷 cis/trans 關係。

**盲點 3 — 亞克隆表觀結構**：
腫瘤並非均質。不同亞克隆可能有不同的甲基化圖譜。TRACERx 等大型研究顯示，腫瘤內甲基化異質性（ITMD）與基因體不穩定性高度相關。ISM 可在 read-level 觀察這些結構。

---

## ═══════════════════════════════════════════
## SLIDE 3 — ONT 長讀序的獨特優勢
## ═══════════════════════════════════════════

```
┌─────────────────────────────────────────────────────┐
│                                                     │
│  ▎ 一條 ONT Read 同時攜帶三種資訊                      │
│                                                     │
│   ┌──────────── 一條長讀序 (~10-100 kb) ────────────┐ │
│   │                                                │ │
│   │  ◆ 鹼基序列 ──→ SNV / Indel 偵測               │ │
│   │  ◆ MM/ML 標籤 ──→ CpG 甲基化機率               │ │
│   │  ◆ HP/PS 標籤 ──→ 單倍型歸屬（phasing）         │ │
│   │                                                │ │
│   └────────────────────────────────────────────────┘ │
│                                                     │
│       短讀序                    長讀序                │
│   ┌───────────┐           ┌───────────────┐         │
│   │ 單一資訊層  │           │ 三層資訊整合    │         │
│   │ 需分開實驗  │    vs     │ 單次定序完成    │         │
│   │ 無法確認    │           │ 可追溯同一分子  │         │
│   │ 同一分子    │           │ cis/trans 關係  │         │
│   └───────────┘           └───────────────┘         │
│                                                     │
│  ◆ 關鍵：三種訊號在同一分子上 → 可做因果推論            │
│                                                     │
└─────────────────────────────────────────────────────┘
```

### [NOTES] 備忘錄

**MM/ML 標籤**：
- MM（Modified base Markup）：記錄 CpG 位點相對於讀序起點的偏移量
- ML（Modified base Likelihood）：每個修飾位點的甲基化機率（0-255，mapping 到 0.0-1.0）
- 無需化學轉換（如 bisulfite），直接從原生電訊號推斷甲基化狀態

**HP/PS 標籤**：
- HP（Haplotype Phase）：此 read 歸屬哪個親代單倍型（1 或 2）
- PS（Phase Set）：同一相位集合的 ID，確保同一區段的 reads 相位一致
- 由 LongPhase / WhatsHap 等工具根據雜合子 SNP 分配

**為什麼重要**：
傳統短讀序定序（Illumina）需要分別做 WGS（突變）+ WGBS/EM-seq（甲基化）+ Hi-C/Trio（相位），且無法確認訊號來自同一 DNA 分子。ONT 長讀序一次解決，且可追溯到單一分子層級。

---

## ═══════════════════════════════════════════
## SLIDE 4 — ISM 核心概念：以 SNV 為錨點的甲基化聚類
## ═══════════════════════════════════════════

```
┌─────────────────────────────────────────────────────┐
│                                                     │
│  ▎ 核心假設                                          │
│                                                     │
│  「在每個 somatic SNV 的局部視窗內，                    │
│    read-level 甲基化圖樣應能形成穩定群集；               │
│    若該 SNV 是真實腫瘤事件，                            │
│    則群集應與生物標籤呈統計顯著關聯」                     │
│                                                     │
│  ─────────────────────────────────────────           │
│                                                     │
│  SNV 位置                                            │
│    ↓                                                │
│  ──────────[====== ±5000bp 視窗 ======]──────────    │
│                                                     │
│  在這個視窗內：                                        │
│  • 提取所有覆蓋此區域的 reads                           │
│  • 對每條 read 取出其 CpG 甲基化機率                    │
│  • 建立 Read × CpG 甲基化矩陣                         │
│  • 計算 Read-Read 距離 → 聚類                          │
│  • 檢定聚類結果 vs 生物標籤 的關聯性                     │
│                                                     │
│  生物標籤：ALT/REF（突變等位）、HP1/HP2（單倍型）         │
│                                                     │
└─────────────────────────────────────────────────────┘
```

### [NOTES] 備忘錄

**為什麼以 SNV 為錨點？**
- Somatic SNV 是已知的腫瘤事件標記
- 以 SNV 為中心開窗，可以自然地取到「帶有突變的 reads」與「不帶突變的 reads」
- 這提供了天然的分組標籤：ALT reads（帶突變）vs REF reads（不帶突變）
- 進一步結合 HP（phasing）標籤，可以判斷突變與甲基化是在同一等位基因（cis）還是不同等位基因（trans）

**視窗大小（±5000bp）**：
- 經測試，5000bp 視窗內平均可涵蓋 ~75 個 CpG 位點
- 太小（1000bp）→ CpG 數量不足（~15 個），統計檢力不足
- 太大 → 甲基化結構被遠端背景稀釋

**Read × CpG 矩陣**：
- 行 = 每條覆蓋此視窗的 read
- 列 = 每個 CpG 位點
- 值 = 甲基化機率（0.0 - 1.0）或二值化（0/1）
- 缺失值 = 該 read 未覆蓋此 CpG（以 NaN 表示）

---

## ═══════════════════════════════════════════
## SLIDE 5 — ISM 技術流程總覽
## ═══════════════════════════════════════════

```
┌─────────────────────────────────────────────────────┐
│                                                     │
│  ▎ InterSubMod 六層處理架構                            │
│                                                     │
│  ┌─────────────────────────────────────────────┐     │
│  │ Layer 1: 甲基化解析                           │     │
│  │  BAM MM/ML 標籤 → CpG 座標 + 甲基化機率       │     │
│  │  （含 CIGAR 校正、正反股處理）                  │     │
│  └──────────────────────┬──────────────────────┘     │
│                         ▼                            │
│  ┌─────────────────────────────────────────────┐     │
│  │ Layer 2: 矩陣建構                            │     │
│  │  Read × CpG 機率矩陣 + 二值矩陣              │     │
│  └──────────────────────┬──────────────────────┘     │
│                         ▼                            │
│  ┌─────────────────────────────────────────────┐     │
│  │ Layer 3: 距離計算                            │     │
│  │  6 種度量: NHD / L1 / L2 / Bernoulli /       │     │
│  │          Jaccard / Correlation                │     │
│  └──────────────────────┬──────────────────────┘     │
│                         ▼                            │
│  ┌─────────────────────────────────────────────┐     │
│  │ Layer 4: 階層式聚類                           │     │
│  │  UPGMA → cluster labels（自動最佳群數）        │     │
│  └──────────────────────┬──────────────────────┘     │
│                         ▼                            │
│  ┌─────────────────────────────────────────────┐     │
│  │ Layer 5: 雙向顯著性檢定                       │     │
│  │  Cluster-First: Fisher + Cramer's V           │     │
│  │  Label-First: PERMANOVA + Delta Distance      │     │
│  └──────────────────────┬──────────────────────┘     │
│                         ▼                            │
│  ┌─────────────────────────────────────────────┐     │
│  │ Layer 6: VerificationClass 判定               │     │
│  │  Strong / Weak / Subclone / Noise             │     │
│  │  + LOH 分層 + QualityScore                    │     │
│  └─────────────────────────────────────────────┘     │
│                                                     │
└─────────────────────────────────────────────────────┘
```

### [NOTES] 備忘錄

**Layer 1 — MethylationParser.cpp**：
- 從 BAM 檔的 MM/ML 輔助標籤解碼甲基化位點
- MM 使用 delta-encoding（相鄰修飾位點的跳躍距離）
- 需處理正股（+C）和反股（-C 互補為 G）的 CpG 辨識
- CIGAR alignment 校正：insertion/deletion/soft-clip 會影響 read 座標到 reference 座標的對應

**Layer 2 — RegionProcessor.cpp**：
- 對每個 SNV 開出 ±5000bp 視窗，收集所有 overlapping reads
- 建立兩種矩陣：raw probability matrix（連續值 0-1）與 binary matrix（≥0.5 → 1, <0.5 → 0）
- 過濾條件：最少 N 條 reads、最少 M 個 CpG 位點

**Layer 3 — DistanceCalculator.cpp**：
- **NHD**（Normalized Hamming Distance）：二值矩陣上，不同位點比例
- **L1/L2**：機率矩陣上的曼哈頓/歐式距離
- **Bernoulli**：ISM 特有設計，信心加權 w(p)=2|p-0.5|，高信心位點貢獻大
- **Jaccard**：集合相似度，只看「都甲基化」的部分
- **Correlation**：Pearson 相關距離

**Layer 4 — 聚類**：
- UPGMA（Unweighted Pair Group Method with Arithmetic Mean）
- 自動用 silhouette score 決定最佳群數（k=2~10）
- 輸出 cluster label 與 dendrogram

**Layer 5 — 顯著性分析**：
- **Cluster-First 路徑**：先聚類，再問「聚類結果與 ALT/REF/HP 標籤是否顯著關聯？」→ Fisher-Freeman-Halton exact test + Cramer's V 效應量
- **Label-First 路徑**：先用 HP/Allele 標籤分組，再問「組間距離是否顯著大於組內？」→ PERMANOVA + Delta Distance
- 兩個路徑獨立運作，結果交叉比對

**Layer 6 — VerificationClass**：
- Strong：兩路徑都顯著且一致 → 最高可信度
- Weak：標籤顯著但全域效應弱 → 有信號但不強
- Subclone：聚類顯著但標籤不顯著 → 非 HP/Allele 驅動的結構
- Noise：都不顯著 → 隨機噪訊

---

## ═══════════════════════════════════════════
## SLIDE 6 — 關鍵定義速查
## ═══════════════════════════════════════════

```
┌─────────────────────────────────────────────────────┐
│                                                     │
│  ▎ 關鍵術語與定義                                     │
│                                                     │
│  ┌─────────────┬───────────────────────────────┐     │
│  │ 術語         │ 定義                          │     │
│  ├─────────────┼───────────────────────────────┤     │
│  │ Somatic SNV │ 腫瘤特有的單核苷酸變異          │     │
│  │ ASM          │ 等位特異性甲基化               │     │
│  │             │ (Allele-Specific Methylation)  │     │
│  │ LOH          │ 雜合度喪失                    │     │
│  │             │ (Loss of Heterozygosity)       │     │
│  │ Second Hit  │ 雙等位失活的第二次打擊           │     │
│  │ Subclone    │ 腫瘤內部的亞克隆群體            │     │
│  │ TP / FP     │ 真陽性 / 假陽性                │     │
│  │ HP (tag)    │ 單倍型標籤 (Haplotype Phase)   │     │
│  │ MM/ML       │ ONT 甲基化位置/機率標籤         │     │
│  │ UPGMA       │ 算術平均非加權配對組法           │     │
│  │ Cramer's V  │ 關聯效應量 (0~1)               │     │
│  │ PERMANOVA   │ 多變量置換變異數分析             │     │
│  │ F1 Score    │ Precision 與 Recall 調和平均    │     │
│  │ CNV         │ 拷貝數變異                     │     │
│  │ Purity      │ 腫瘤純度（腫瘤細胞佔比）         │     │
│  └─────────────┴───────────────────────────────┘     │
│                                                     │
└─────────────────────────────────────────────────────┘
```

### [NOTES] 備忘錄

**Somatic SNV vs Germline SNV**：
- Somatic SNV：只存在於腫瘤細胞中的突變，不遺傳
- Germline SNV：存在於所有細胞中的變異（包括正常細胞），可遺傳
- ISM 的輸入來自 somatic variant caller（如 ClairS），只分析被認為是 somatic 的位點

**ASM（等位特異性甲基化）**：
- 同一個 CpG 位點，在兩個等位基因上的甲基化狀態不同
- 正常組織中存在 imprinting ASM（如 H19/IGF2）
- 癌症中可出現 tumor-specific ASM（如 BRCA1 promoter 被甲基化沉默）
- ISM 的挑戰：無法區分 germline ASM 與 somatic ASM（主要 FP 來源）

**LOH（雜合度喪失）**：
- 某區域的一個等位基因整段遺失（deletion）或被另一個等位基因取代（UPD）
- 在 LOH 區域，ALT reads ≈ 腫瘤來源，REF reads ≈ 正常污染
- ISM 用 HP ratio（HP1/(HP1+HP2)）偵測 LOH：ratio < 0.1 或 > 0.9

**Second Hit（雙等位失活）**：
- Knudson 假說：TSG 需要兩次打擊才失活
- 常見組合：mutation + LOH、mutation + promoter methylation、LOH + methylation
- ISM 的長期目標：在 read-level 判斷 second hit 的先後順序

---

## ═══════════════════════════════════════════
## SLIDE 7 — VerificationClass：四層分類系統
## ═══════════════════════════════════════════

```
┌─────────────────────────────────────────────────────┐
│                                                     │
│  ▎ VerificationClass — ISM 核心輸出                   │
│                                                     │
│  ┌─────────────────────────────────────────────┐     │
│  │                                             │     │
│  │  Cluster-First     Label-First              │     │
│  │  (先聚類再看標籤)   (先用標籤分組再看距離)      │     │
│  │       ↓                  ↓                  │     │
│  │       └──── 交叉比對 ────┘                   │     │
│  │              ↓                              │     │
│  │  ┌──────────────────────────────────┐       │     │
│  │  │                                  │       │     │
│  │  │  Strong   ← 兩路徑都顯著          │       │     │
│  │  │             Precision 81~100%    │       │     │
│  │  │                                  │       │     │
│  │  │  Weak     ← 標籤顯著，全域弱      │       │     │
│  │  │             Precision 90~99%     │       │     │
│  │  │                                  │       │     │
│  │  │  Subclone ← 聚類顯著，非標籤驅動   │       │     │
│  │  │             需分層：LOH vs Non-LOH │       │     │
│  │  │                                  │       │     │
│  │  │  Noise    ← 都不顯著             │       │     │
│  │  │             含弱信號 TP（33~67%）  │       │     │
│  │  │                                  │       │     │
│  │  └──────────────────────────────────┘       │     │
│  │                                             │     │
│  └─────────────────────────────────────────────┘     │
│                                                     │
│  推薦使用方式：VC != Noise（捕獲率 33~68%）             │
│  比 Significant=True（1~12%）大幅提升                  │
│                                                     │
└─────────────────────────────────────────────────────┘
```

### [NOTES] 備忘錄

**Cluster-First 路徑詳解**：
1. 對 Read × CpG 矩陣做距離計算 → 聚類
2. 取得 cluster labels（每條 read 屬於哪個 cluster）
3. 建立 cluster × label 列聯表（label = ALT/REF 或 HP1/HP2）
4. Fisher-Freeman-Halton exact test（p-value）+ Cramer's V（效應量）
5. 通過 gate（p ≤ 0.1 AND V ≥ 0.1）→ 進入下一層分析

**Label-First 路徑詳解**：
1. 直接用 HP 或 Allele 標籤把 reads 分成兩組
2. 計算組間距離 vs 組內距離
3. PERMANOVA 檢定組間差異是否顯著
4. Delta Distance = mean(Between) - mean(Within)

**Subclone 分層的重要性**：
- LOH_Subclone（HCC1395/HCC1937 佔 87-89%）：LOH 驅動，語意明確，可靠
- Non-LOH Subclone（COLO829/H1437 佔 80-90%）：多為 unreliable CramersV 觸發的噪訊
- 使用建議：`(VC == 'Subclone') AND (LOH_Subtype == 'LOH_Subclone')` 才可信

**與 Significant 的比較**：
- 舊指標 `Significant=True` 僅捕獲 1-12% 的真陽性（太嚴格）
- `VC != Noise` 可捕獲 33-68%，且 Precision 相近
- 這是方法學審查的重要結論之一

---

## ═══════════════════════════════════════════
## SLIDE 8 — 完整分析管線示意
## ═══════════════════════════════════════════

```
┌─────────────────────────────────────────────────────┐
│                                                     │
│  ▎ InterSubMod 在整體分析管線中的位置                   │
│                                                     │
│  ┌────────────────────────────────────────────────┐  │
│  │              上游 Variant Calling                │  │
│  │  ONT BAM ──→ ClairS / DeepSomatic ──→ VCF     │  │
│  │              (somatic SNV calls)                │  │
│  └──────────────────────┬─────────────────────────┘  │
│                         ▼                            │
│  ┌────────────────────────────────────────────────┐  │
│  │              上游 Phasing                       │  │
│  │  VCF + BAM ──→ LongPhase-S ──→ Tagged BAM     │  │
│  │              (HP/PS tags added)                 │  │
│  └──────────────────────┬─────────────────────────┘  │
│                         ▼                            │
│  ┌────────────────────────────────────────────────┐  │
│  │         ★ InterSubMod 核心分析 ★                │  │
│  │  Tagged BAM + VCF                               │  │
│  │  → per-SNV 甲基化矩陣                            │  │
│  │  → 距離計算 + 聚類                               │  │
│  │  → 雙向顯著性檢定                                │  │
│  │  → VerificationClass + annotations              │  │
│  └──────────────────────┬─────────────────────────┘  │
│                         ▼                            │
│  ┌────────────────────────────────────────────────┐  │
│  │              下游整合與判讀                       │  │
│  │  • caller-first → methylation-support           │  │
│  │  • 甲基化 annotation 輔助變異優先級排序            │  │
│  │  • 基因層級 second hit 整合（未來方向）             │  │
│  └────────────────────────────────────────────────┘  │
│                                                     │
│  確認的最佳架構：                                      │
│  caller-first → methylation-support → artifact triage│
│                                                     │
└─────────────────────────────────────────────────────┘
```

### [NOTES] 備忘錄

**三層架構已通過方法學審查確認**：
1. **Caller-first**：以 variant caller 的結果為主要判決（GQ、QUAL 等）
2. **Methylation-support**：ISM 的 VerificationClass、QualityScore 等作為輔助支持
3. **Artifact triage**：對可疑位點做甲基化層面的標記（不直接刪除）

**為什麼不是 methylation-first？**
- ISM 實驗（37+ 次）反覆驗證：甲基化特徵的 TP/FP 鑑別力（AUROC~0.5-0.62）遠不及 caller 自身的 GQ/QUAL（AUROC~0.92-0.97）
- 甲基化的價值在於「提供額外維度的觀察」，而非「替代 caller 做決策」
- 類比：甲基化 = 病理切片染色，variant caller = 基因檢測報告。兩者互補，不能相互替代

**實際效能數據**：
- HCC1395 5kHz paired：ClairS F1=0.8443 → LongPhase-S F1=0.8522 → ISM F1=0.8532（+0.0010）
- ISM 的增益小但穩定，主要貢獻在 annotation 而非 hard filter

---

## ═══════════════════════════════════════════
## SLIDE 9 — 實驗驗證結果總覽
## ═══════════════════════════════════════════

```
┌─────────────────────────────────────────────────────┐
│                                                     │
│  ▎ 方法學審查結果（2025-11 ~ 2026-03）                 │
│                                                     │
│  37+ 次實驗 ┃ 16 份方法學觀察文件 ┃ 7 個 canonical 樣本 │
│                                                     │
│  ── 已確認的正面發現 ──────────────────────────        │
│                                                     │
│  ✓ VerificationClass 框架整體合理                      │
│  ✓ LOH_Subclone 是可靠的 Subclone 子類                │
│  ✓ VC != Noise 顯著優於 Significant=True              │
│  ✓ Tumor-Only: GQ>=3 rescue 可達 F1 +0.0094         │
│                                                     │
│  ── 已確認的瓶頸 ──────────────────────────            │
│                                                     │
│  ✗ 甲基化 TP/FP 鑑別力 AUROC ~0.5-0.62               │
│  ✗ FP 主因：germline ASM（98.7% FP 靠 CramersV 通過） │
│  ✗ 跨樣本不穩定：PairwiseMedianDist 方向依賴            │
│  ✗ ISM 覆蓋 FN pool 僅 7%（rescue 上限有限）           │
│                                                     │
│  ── 已排除的方向 ──────────────────────────            │
│                                                     │
│  ✗ AlleleDelta 單獨 filter（AUROC 0.41-0.76）         │
│  ✗ LOH 整合 VerificationClass（無鑑別力）              │
│  ✗ PERMANOVA 作為 filter（AUROC 0.53-0.59）           │
│                                                     │
└─────────────────────────────────────────────────────┘
```

### [NOTES] 備忘錄

**37+ 次實驗的範圍**：
從最初的 MM/ML 解析驗證、距離計算優化、PERMANOVA 實作，到 F1 最佳化、purity-aware 策略、跨樣本驗證、TO rescue 規則探索等。完整列表見 `docs/experiments/INDEX.md`。

**7 個 canonical 樣本**：
HCC1395（乳癌）、HCC1395_DORADO（同細胞株不同定序平台）、COLO829（黑色素瘤）、H1437（肺腺癌）、H2009（肺腺癌）、HCC1937（乳癌）、HCC1954（乳癌）。

**FP 根因分析**：
方法學審查最重要的結論之一：98.7% 的假陽性是靠 unreliable Cramer's V 通過 gate 的。這些 FP 來自 germline ASM 區域 — 即正常組織中本來就存在的等位特異性甲基化。ISM 正確偵測到甲基化差異，但無法區分這是「胚系」還是「體細胞」造成的。解決需要 normal sample 甲基化參考。

**Tumor-Only (TO) 場景**：
不使用 matched normal 的分析場景。ISM 在此場景下的 rescue 規則：GQ>=3 可從 FN pool 中救回部分漏掉的 TP，delta F1=+0.0094。但 ISM 覆蓋的 FN pool 僅佔全部 FN 的 7%（773/11051），因此 rescue 的絕對上限有限。

---

## ═══════════════════════════════════════════
## SLIDE 10 — ISM 在文獻中的定位
## ═══════════════════════════════════════════

```
┌─────────────────────────────────────────────────────┐
│                                                     │
│  ▎ 12 篇核心文獻對照：ISM 的獨特性                     │
│                                                     │
│    NanoMethPhase                  MethylBERT         │
│    (全基因組 ASM)                 (ML 甲基化分類)      │
│         ↑                             ↑             │
│         │                             │             │
│    ─────┴──── 相關但不同 ────┬─────────┘             │
│                              │                      │
│              ┌───────────────┴───────────────┐       │
│              │                               │       │
│              │    ★ InterSubMod 獨有 ★        │       │
│              │                               │       │
│              │  以 somatic SNV 為錨點          │       │
│              │  + read-level 甲基化聚類        │       │
│              │  + 雙向顯著性驗證               │       │
│              │  → 完整框架尚無對標             │       │
│              │                               │       │
│              └───────────────┬───────────────┘       │
│                              │                      │
│    ─────┬──── 互補上游 ──────┴──────────┐            │
│         │                              │            │
│    Long-Read POG              TRACERx CAMDAC        │
│    (gene-level                (purity/CN-aware      │
│     evidence)                  methylation)         │
│                                                     │
│  結論：沒有人完全做過 ISM 想做的事                       │
│  但 ISM 的獨特性也是瓶頸所在                            │
│                                                     │
└─────────────────────────────────────────────────────┘
```

### [NOTES] 備忘錄

**12 篇核心文獻**：
1. Long-Read POG (Cell Genomics 2024) — 長讀序整合 CNV+ASM+ASE+somatic
2. TRACERx NSCLC (Nature Genetics 2025) — CAMDAC purity/CN-aware methylation
3. EVOFLUx (Nature 2025) — fluctuating CpGs 追蹤癌症演化
4. DeepSomatic (bioRxiv 2024) — DL-based somatic variant caller
5. LongPhase-S (bioRxiv 2025) — somatic haplotyping
6. t-nanoEM (Cell Reports Methods 2025) — targeted 長讀序甲基化
7. NanoMethPhase (Genome Biology 2021) — 全基因組 haplotype-resolved methylation
8. MethSig (Cancer Discovery 2021) — 甲基化驅動基因識別
9. MethPhaser (Nature Communications 2024) — 甲基化延伸 phasing
10. MethylBERT (Nature Communications 2025) — Transformer read-level 甲基化
11. PRISM (Bioinformatics 2019) — 甲基化 pattern subclone inference
12. MHB pan-cancer (Cell Reports 2025) — 甲基化單倍型區塊

**ISM 的獨特組合**：
「somatic SNV anchoring + read-level methylation matrix + hierarchical clustering + dual-path significance testing (cluster-first & label-first) + VerificationClass」— 這個完整組合在現有文獻中確認無對標。

**最接近的競爭者**：
- NanoMethPhase：同為 ONT + methylation + phasing，但做全基因組 ASM，不是 per-SNV clustering
- MethylBERT：同為 read-level methylation pattern，但用 supervised ML，目標是 tumor deconvolution
- PRISM：同為 methylation → subclone，但用 beta-binomial mixture on bulk，不是 read-level with variant anchor

---

## ═══════════════════════════════════════════
## SLIDE 11 — FP 根因與核心瓶頸
## ═══════════════════════════════════════════

```
┌─────────────────────────────────────────────────────┐
│                                                     │
│  ▎ 為什麼甲基化鑑別力有限？根因分析                      │
│                                                     │
│  ┌─────────────────────────────────────────────┐     │
│  │                                             │     │
│  │   Somatic ASM                               │     │
│  │   (腫瘤驅動甲基化差異)                        │     │
│  │        │                                    │     │
│  │        │   ← ISM 無法區分                    │     │
│  │        │                                    │     │
│  │   Germline ASM                              │     │
│  │   (胚系等位特異性甲基化)                       │     │
│  │                                             │     │
│  │   ┌─────────────────────────────────┐       │     │
│  │   │ 98.7% FP 靠 CramersV 通過 gate  │       │     │
│  │   │ FP AlleleDelta 中位 = 0.155     │       │     │
│  │   │ FP AlleleSig=True 比率 = 94%    │       │     │
│  │   │                                 │       │     │
│  │   │ → 這些 FP 有真實甲基化差異       │       │     │
│  │   │   但差異源自胚系，非腫瘤          │       │     │
│  │   └─────────────────────────────────┘       │     │
│  │                                             │     │
│  └─────────────────────────────────────────────┘     │
│                                                     │
│  解決方案：引入 Normal Sample 甲基化參考                 │
│  （下一階段研究方向）                                   │
│                                                     │
└─────────────────────────────────────────────────────┘
```

### [NOTES] 備忘錄

**Germline ASM 是什麼？**
正常人體中，某些基因組區域在兩個等位基因上本來就有不同的甲基化模式（例如 imprinting、X-inactivation、random monoallelic expression）。當 somatic variant caller 在這些區域錯誤呼叫一個 germline variant 為 somatic 時，ISM 會偵測到真實的甲基化差異 — 但這個差異與腫瘤無關。

**數據支持（HCC1395）**：
- FP 的 AlleleDelta（等位甲基化差異）中位數 = 0.155，是 TP 的 2.5 倍
- FP 中 94% 帶有 AlleleSig=True（等位甲基化顯著）
- 相較之下，TP 中只有 76.4% 帶 AlleleSig=True
- 這意味著 FP 在甲基化層面的「信號」比 TP 還強！

**為什麼 gate 抓不到？**
ISM 的 gate 條件是 p ≤ 0.1 AND CramersV ≥ 0.1。Germline ASM 區域確實有顯著的聚類-標籤關聯（因為甲基化差異是真實的），所以它們能通過 gate。問題不在 gate 設計，而在 ISM 缺乏「這個甲基化差異是胚系還是體細胞」的判斷能力。

**解決方向**：
- Paired 場景：用 normal BAM 的甲基化作為 baseline，ISM 聚類結果若與 normal 一致 → downgrade 為 germline ASM
- Tumor-only 場景：需要 Panel of Normals for Methylation（目前不存在）

---

## ═══════════════════════════════════════════
## SLIDE 12 — 未來突破方向（四階段策略）
## ═══════════════════════════════════════════

```
┌─────────────────────────────────────────────────────┐
│                                                     │
│  ▎ 漸進式突破策略（已確認優先序）                        │
│                                                     │
│  ┌──────────────────────────────────────────┐        │
│  │ Phase 1: ML Read Classification           │        │
│  │  • MethylBERT 啟發的 read-level 模式辨識    │        │
│  │  • 突破 clustering 天花板                  │        │
│  │  • 在現有框架內改進                        │        │
│  └─────────────────┬────────────────────────┘        │
│                    ▼                                 │
│  ┌──────────────────────────────────────────┐        │
│  │ Phase 2: Normal Ref + CN/Purity-aware    │        │
│  │  • 引入 normal 甲基化基線                  │        │
│  │  • CAMDAC 啟發的 purity/LOH/CNV 校正      │        │
│  │  • 解決 germline ASM 根因                 │        │
│  └─────────────────┬────────────────────────┘        │
│                    ▼                                 │
│  ┌──────────────────────────────────────────┐        │
│  │ Phase 3: Gene-level Integration           │        │
│  │  • 多 SNV 窗口 → 基因層級證據彙總           │        │
│  │  • Biallelic inactivation narrative        │        │
│  │  • Second hit 順序推斷                     │        │
│  └─────────────────┬────────────────────────┘        │
│                    ▼                                 │
│  ┌──────────────────────────────────────────┐        │
│  │ Phase 4: CpG Functional Stratification   │        │
│  │  • Regulatory / fCpG / MHB 分層           │        │
│  │  • 高關聯位點篩選與輸出                     │        │
│  │  • Subclone 演化時間軸                     │        │
│  └──────────────────────────────────────────┘        │
│                                                     │
│  最終目標：判斷 somatic 位點與甲基位點                   │
│           是否為癌症驅動 TP                             │
│                                                     │
└─────────────────────────────────────────────────────┘
```

### [NOTES] 備忘錄

**Phase 1 — ML Read Classification（方向 E）**：
- 靈感來源：MethylBERT（Nature Communications 2025）
- 目標：用 Transformer 或其他 ML 模型學習 read-level methylation pattern → tumor/normal classification
- 取代或補充 ISM 的 unsupervised clustering
- 預期收益：MethylBERT 已展示 read-level pattern 可超越 mean methylation 的判讀能力
- 風險：需要 training data（labeled tumor/normal reads）；泛化能力未經 subclone 驗證

**Phase 2 — Normal Ref + CN/Purity-aware（方向 A + D）**：
- 靈感來源：TRACERx CAMDAC（Nature Genetics 2025）+ 方法學審查結論
- 目標：引入 normal sample 甲基化作為 baseline + 用 LOH-SNV 做 purity/CN 校正
- 這是解決 germline ASM 根因的唯一已知路徑
- 預期收益：直接消除目前最大的 FP 來源（98.7%）
- 風險：Tumor-only 場景需要 methylation PoN（目前不存在）

**Phase 3 — Gene-level Integration（方向 B）**：
- 靈感來源：Long-Read POG IMPALA（Cell Genomics 2024）
- 目標：把同一基因內多個 SNV 窗口的 ISM 結果彙總，建立 gene-level evidence panel
- 回答「此基因是否有 biallelic inactivation 證據？」
- 定義 second-hit archetypes：mutation→LOH、LOH→mutation、promoter_methylation→LOH 等
- 預期收益：從 region-level significance 升級到 clinically actionable gene-level narrative

**Phase 4 — CpG Functional Stratification（方向 C）**：
- 靈感來源：EVOFLUx fCpGs（Nature 2025）+ TRACERx M_R/M_N + MHB atlas（Cell Reports 2025）
- 目標：將 CpG 依功能分層（regulatory / fCpG-like / MHB-core），不同目標用不同 CpG 子集
- 預期收益：直接提升 signal-to-noise ratio
- 輸出：annotated CpG catalogue + subclone evolution timeline

---

## ═══════════════════════════════════════════
## SLIDE 13 — 學界共識與 ISM 的機會
## ═══════════════════════════════════════════

```
┌─────────────────────────────────────────────────────┐
│                                                     │
│  ▎ 學界已有共識（ISM 可直接借鑑）                       │
│                                                     │
│  1. Purity/CNV 校正是必要的                           │
│     → ISM Phase 2 目標                               │
│                                                     │
│  2. 長讀序的多模態價值已被確認                           │
│     → ISM 的立基優勢                                  │
│                                                     │
│  3. LOH-SNV = 天然 tumor/normal 分離標籤              │
│     → ISM 已有 ALT/REF 標記，可直接升級                │
│                                                     │
│  4. Promoter methylation + LOH = second hit          │
│     → ISM Phase 3 目標                               │
│                                                     │
│  ▎ 尚未解決（ISM 的機會空間）                           │
│                                                     │
│  ✦ read-level per-SNV 甲基化 cluster 的生物意義         │
│  ✦ clustering vs ML 比較（無 head-to-head 研究）       │
│  ✦ germline vs somatic ASM 的區分方法                  │
│  ✦ read-level 甲基化特徵的跨平台穩定性                  │
│                                                     │
└─────────────────────────────────────────────────────┘
```

### [NOTES] 備忘錄

**共識的外部支持**：
1. Purity/CNV 校正：TRACERx CAMDAC、Long-Read POG、多項 deconvolution 研究一致
2. 多模態價值：Long-Read POG、NanoMethPhase、t-nanoEM、LongPhase-S 一致展示
3. LOH-SNV 策略：TRACERx 用此驗證 CAMDAC PDR（R>0.8）；Long-Read POG 用此確認 biallelic inactivation
4. Second hit：BRCA1/RAD51C/CDKN2A 等已有多項臨床驗證

**ISM 的獨特機會空間**：
- 目前學界在 read-level per-SNV 的甲基化結構分析上仍是空白
- MethylBERT 最接近但做的是 deconvolution 而非 per-variant context
- 如果 ISM 能在 Phase 1-2 解決技術限制，Phase 3-4 的 gene-level narrative + CpG stratification 將進入無人探索的領域
- 這正是 ISM 最有可能產出高影響力論文的定位

---

## ═══════════════════════════════════════════
## SLIDE 14 — 論文定位建議
## ═══════════════════════════════════════════

```
┌─────────────────────────────────────────────────────┐
│                                                     │
│  ▎ ISM 的學術定位                                     │
│                                                     │
│  ┌───────────────────────────────────────────────┐   │
│  │                                               │   │
│  │   ✗ 不是 "variant filter"                     │   │
│  │     （鑑別力不足以獨立做過濾）                   │   │
│  │                                               │   │
│  │   ✓ 是 "variant interpreter"                  │   │
│  │                                               │   │
│  │   「以長讀序 somatic variant 為錨點，              │   │
│  │     整合甲基化 pattern、haplotype、               │   │
│  │     copy-number 資訊，                            │   │
│  │     提供 read-level epigenetic context           │   │
│  │     for variant interpretation」                 │   │
│  │                                               │   │
│  └───────────────────────────────────────────────┘   │
│                                                     │
│  與現有研究的互補關係：                                  │
│  • Long-Read POG ── gene-level evidence              │
│  • TRACERx ──── purity/CN-aware methylation          │
│  • EVOFLUx ──── evolution barcode                    │
│  • t-nanoEM ──── clinical targeted methylation       │
│  → ISM 提供：per-variant read-level 的深度分析         │
│                                                     │
└─────────────────────────────────────────────────────┘
```

### [NOTES] 備忘錄

**為什麼不定位為 variant filter？**
- 37+ 次實驗反覆驗證：甲基化特徵 AUROC 0.5-0.62，無法獨立超越 caller 的 GQ/QUAL（AUROC 0.92-0.97）
- 強行定位為 filter 會面臨「為什麼不直接用 caller 的分數？」的致命問題
- 學術上也無法與 DeepSomatic、ClairS 等 caller 競爭

**為什麼 variant interpreter 是更好的定位？**
- ISM 提供的資訊是 caller 無法給出的：
  - 這個 SNV 附近的甲基化結構是什麼樣的？
  - 這個 SNV 的突變等位基因與甲基化差異是 cis 還是 trans？
  - 這個 SNV 是否位於 LOH 區域？如果是，甲基化模式是否支持 second hit？
  - 這個 SNV 附近是否有亞克隆結構？
- 這些問題目前沒有任何工具可以從 read-level 回答

**與現有研究的互補定位**：
- Long-Read POG：在 gene level 做整合，ISM 在 variant level 做深度分析 → 可作為 POG-style pipeline 的 per-variant module
- TRACERx CAMDAC：在 bulk level 做 deconvolution，ISM 在 read-level 做聚類 → CAMDAC 的 purity 估計可作為 ISM 的先驗
- EVOFLUx：找 fCpGs 做演化追蹤，ISM 的 per-SNV 結構可提供 variant-anchored 的演化觀點
- t-nanoEM：targeted panel 高深度，ISM 在 WGS 全基因組 → 不同應用場景

---

## ═══════════════════════════════════════════
## SLIDE 15 — 技術實作亮點
## ═══════════════════════════════════════════

```
┌─────────────────────────────────────────────────────┐
│                                                     │
│  ▎ C++17 高效能核心                                   │
│                                                     │
│  ┌───────────────────────────────────────────────┐   │
│  │ 效能                                          │   │
│  │  • 單 Region 平均 < 300ms                     │   │
│  │  • OpenMP 平行化（thread-local BAM/Fasta）     │   │
│  │  • 全基因組 30K+ regions 可在數小時內完成        │   │
│  ├───────────────────────────────────────────────┤   │
│  │ 距離度量                                      │   │
│  │  • 6 種可選: NHD / L1 / L2 / Bernoulli /     │   │
│  │    Jaccard / Correlation                      │   │
│  │  • Bernoulli: 信心加權 w(p)=2|p-0.5|          │   │
│  │    → 降低 ONT p≈0.5 雜訊位點的影響             │   │
│  ├───────────────────────────────────────────────┤   │
│  │ 輸出格式                                      │   │
│  │  • CSV: 甲基化矩陣、距離矩陣、顯著性摘要        │   │
│  │  • TSV: read 資訊、cluster label              │   │
│  │  • PNG: 距離熱圖、分群熱圖（自動 Python 生成）   │   │
│  │  • TXT: region metadata                       │   │
│  └───────────────────────────────────────────────┘   │
│                                                     │
└─────────────────────────────────────────────────────┘
```

### [NOTES] 備忘錄

**Bernoulli 距離 — ISM 的特色設計**：
ONT 的甲基化機率不是 0/1 二值，而是連續機率（0.0-1.0）。中間值（p ≈ 0.5）表示不確定。傳統的 Hamming 或 Jaccard 距離在二值化後會丟失這個不確定性資訊。

Bernoulli 距離的設計：
- 視每個 CpG 為 Bernoulli 隨機變數
- 兩條 reads 在某 CpG 的期望不一致 = p1(1-p2) + p2(1-p1)
- 再乘以信心權重 w = 2|p-0.5|
- 效果：p=0.05 vs p=0.95 → 高距離、高信心；p=0.45 vs p=0.55 → 低距離、低信心

**輸出結構（每個 SNV region）**：
```
region_chr1_12345/
├── metadata.txt           # 區域基本資訊
├── reads.tsv              # 每條 read 的 ALT/REF、HP、cluster 等
├── methylation.csv        # Read × CpG 甲基化矩陣
├── distance_matrix_NHD.csv # Read-Read 距離矩陣
├── significance_summary.csv # 顯著性分析結果
├── distance_heatmap.png   # 距離矩陣熱圖
└── cluster_heatmap.png    # 甲基化+聚類熱圖
```

---

## ═══════════════════════════════════════════
## SLIDE 16 — 總結與展望
## ═══════════════════════════════════════════

```
┌─────────────────────────────────────────────────────┐
│                                                     │
│  ▎ InterSubMod — 總結                                │
│                                                     │
│  ┌───────────────────────────────────────────────┐   │
│  │ 我們做了什麼                                    │   │
│  │  • 建立完整的 per-SNV read-level 甲基化          │   │
│  │    聚類與雙向顯著性驗證框架                       │   │
│  │  • 37+ 實驗驗證，16 份方法學審查文件              │   │
│  │  • 確認框架合理性與當前瓶頸                       │   │
│  ├───────────────────────────────────────────────┤   │
│  │ 我們發現什麼                                    │   │
│  │  • 甲基化是 variant interpretation 的有力補充     │   │
│  │  • 但不足以獨立做 TP/FP 過濾                     │   │
│  │  • 主要限制：germline ASM 無法區分               │   │
│  │  • 文獻中尚無完整對標研究                        │   │
│  ├───────────────────────────────────────────────┤   │
│  │ 我們接下來要做什麼                               │   │
│  │  • Phase 1: ML 模式辨識（突破 clustering 天花板）  │   │
│  │  • Phase 2: Normal + CN/Purity 校正（解決根因）   │   │
│  │  • Phase 3: Gene-level 整合（臨床 narrative）     │   │
│  │  • Phase 4: CpG 分層（演化推斷）                  │   │
│  └───────────────────────────────────────────────┘   │
│                                                     │
│  願景：提供 read-level epigenetic context              │
│       讓每個 somatic variant 都有甲基化層面的解讀       │
│                                                     │
└─────────────────────────────────────────────────────┘
```

### [NOTES] 備忘錄

**回顧與前瞻**：

InterSubMod 從 2025 年 11 月開始開發，歷經 4 個月密集的實驗與方法學審查。目前已建立完整的技術基礎設施（C++ 核心 + Python 視覺化 + 自動化流程），並通過 37+ 次實驗驗證了框架的合理性與瓶頸。

最重要的學習是：甲基化資訊的價值不在「過濾」，而在「解讀」。這個發現重新定義了 ISM 的學術定位，也指引了未來四個階段的研究方向。

**長期願景**：
每個 somatic variant 不再只是一個「座標 + 基因型」，而是附帶完整的表觀遺傳 context：
- 這個位點附近的甲基化結構是什麼？
- 與哪個等位基因、哪個單倍型關聯？
- 是否位於 LOH/CNV 區域？
- 所在基因是否有 biallelic inactivation 證據？
- 甲基化變化是驅動事件還是乘客事件？

這就是「variant interpretation」的完整圖景。

---

## ═══════════════════════════════════════════
## 附錄 A — 樣本資訊
## ═══════════════════════════════════════════

### [NOTES] 備忘錄（附錄不需投影片）

**7 個 canonical 測試樣本**：

| 樣本 | 癌症類型 | 平台 | 特點 |
|------|---------|------|------|
| HCC1395 5kHz | 乳腺導管癌 | ONT 5kHz | 主要驗證樣本，paired with HCC1395BL |
| HCC1395_DORADO | 乳腺導管癌 | ONT DORADO | 跨平台對照 |
| COLO829 | 惡性黑色素瘤 | ONT | 高 SV 負荷 |
| H1437 | 肺腺癌 | ONT | 非小細胞肺癌 |
| H2009 | 肺腺癌 | ONT | 顯著性率最高（3.8%） |
| HCC1937 | 乳腺癌（BRCA1-/-） | ONT | 距離最低（0.097） |
| HCC1954 | 乳腺導管癌 | ONT | 高突變負荷 |

---

## ═══════════════════════════════════════════
## 附錄 B — 參考文獻
## ═══════════════════════════════════════════

### [NOTES] 備忘錄

1. O'Neill et al., Cell Genomics 4(11), 2024 — Long-Read POG
2. TRACERx consortium, Nat Genet, 2025 — NSCLC methylation + CAMDAC
3. Gabbutt & Duran-Ferrer et al., Nature, 2025 — EVOFLUx
4. Park et al., bioRxiv, 2024 — DeepSomatic
5. Ho et al., bioRxiv, 2025 — LongPhase-S
6. Kunigo et al., Cell Reports Methods, 2025 — t-nanoEM
7. Akbari et al., Genome Biology 22, 2021 — NanoMethPhase
8. Landau et al., Cancer Discovery 11(9), 2021 — MethSig
9. Fu et al., Nat Commun 15:5327, 2024 — MethPhaser
10. Jeong et al., Nat Commun, 2025 — MethylBERT
11. Lee et al., Bioinformatics 35(14), 2019 — PRISM
12. Cell Reports, 2025 — MHB pan-cancer atlas

---

## ═══════════════════════════════════════════
## 附錄 C — 詳細數據表
## ═══════════════════════════════════════════

### [NOTES] 備忘錄

**VerificationClass 跨樣本 Precision**：

| 樣本 | Strong | Weak | Subclone | Noise(含TP%) |
|------|--------|------|----------|-------------|
| HCC1395 | 81-89% | 90-95% | 需分層 | 33-45% |
| COLO829 | 95-100% | 97-99% | 需分層 | 50-67% |
| H1437 | 90-96% | 95-98% | 需分層 | 40-55% |

**TO Rescue 規則效能（HCC1395 5kHz）**：

| 規則 | delta F1 | Precision | 備註 |
|------|---------|-----------|------|
| GQ >= 3 | +0.0094 | 74.1% | 最高 F1 |
| QS >= 50 | +0.0086 | 76.9% | 最高 Precision |
| GQ >= 3 + QS >= 50 | = H011 | 子集關係 | 不可正交組合 |

**統計方法效能**：

| 指標 | 全域 AUROC | 備註 |
|------|-----------|------|
| CramersV | 0.52 | 幾乎無效 |
| QualityScore | 0.62 | 有限改善 |
| PairwiseMedianDist | 0.53 | 方向不穩定 |
| PERMANOVA | 0.53-0.59 | 僅 annotation |
| VCF QUAL (caller) | 0.97 | 上游最佳 |
| VCF AF (caller) | 0.92 | 上游次佳 |
