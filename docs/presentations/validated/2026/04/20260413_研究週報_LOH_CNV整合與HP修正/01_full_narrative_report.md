<!--
建立時間: 2026-04-13 22:00
目標: PPTX 完整敘事稿（延續式故事線）
處理範圍: 2026-04-07 ~ 2026-04-13 LOH×CNV 整合觀察與 HP Tag 修正
關聯檔案:
  - docs/presentations/validated/2026/04/20260408_研究週報v3_完整推論鏈與獨立建模/pptx_config.json
  - docs/experiments/in_progress/2026/04/20260409_SEQC2_CNV分層觀察_01.md
  - docs/experiments/validated/2026/04/20260411_GC校正與甲基化CN驗證_01.md
  - docs/experiments/validated/2026/04/20260413_Phase_BCD_Dual_BAM_Validation_01.md
  - research/literature_validation/reports/20260412_文獻驗證綜合報告_01.md
  - docs/experiments/finalized/2026/04/20260409_beyond_auc_comprehensive_validation_01.md
-->

# 完整敘事稿：LOH×CNV 整合觀察與 HP Tag 修正

## 敘事主軸：延續 → 深化 → 修正 → 現況

上週確立 LOH.bed 可信（SEQC2 Jaccard=0.928）、LOH 內外特徵全面關閉、ISM 轉向 characterization。本週從上週的 LOH 分析出發，沿四條路線深化：(1) LOH 與 CNV 的交叉觀察揭示分割判斷的必要性；(2) HP Tag 上游修正從 Self-Phasing 到 V5；(3) 七種統計方法確認甲基化特徵空間耗盡；(4) Phase 2 Dual-BAM 驗證成功，確立轉向方向。

---

## Act I：從 LOH 到 LOH×CNV（~11 min, Slides 3-12）

### 上週回顧（5 頁精簡版，~5 min）

#### Slide 3：推論地圖

上週的完整推論鏈：

```
LOH 定義驗證 → ISM LOH 影響 → Filter FAIL → Non-LOH 無突破 → 特徵設計 → 獨立建模 → 轉向
```

七步推論鏈中，每一步都有數據支撐、每個判定都是確定性結論。這張地圖是本週所有延伸分析的起點。

#### Slide 4：SEQC2 + 四象限

- **SEQC2 金標準**：Jaccard=0.928，Sensitivity=0.961，Precision=0.964 — LOH.bed 與 FDA 多平台共識幾乎完全吻合
- **四象限分佈**：Q1 Both=26.7%, Q2 HP-only=15.2%（phasing 偏差非弱 LOH）, Q3 LOH.bed-only≈0%, Q4 Neither=58.1%

這兩個結果共同確立：LOH.bed 是可靠的分層基礎。

#### Slide 5：LOH 影響 + Filter FAIL

- **LOH 內甲基化全面失效**：PERMANOVA valid rate 僅 5-6%，7/7 樣本 AUC~0.50
- **10/10 filter 全 FAIL**：LOH FP rate=0.239 < Non-LOH 0.338 — LOH 區域反而是 **TP 富集區**
- 已實作 QS mode-aware：TO 模式停用 LOH penalty 與 verify bonus

#### Slide 6：Non-LOH + caller_af

- **Non-LOH 最高 AUC < 0.58**：HPFineNGroups 0.643 排除 confound 後失效
- **多特徵 Voting AUC=0.577**：cnLOH 0.587 是 Simpson's Paradox（per-sample=0.50）
- **caller_af AUC=0.654 唯一超過門檻** — 單一 caller 特徵超越全部 ISM 特徵
- **TO 獨立建模**：ISM 在 caller 之上僅增加 +0.003~0.030

#### Slide 7：上週結論 → 銜接

上週確定性結論：
1. LOH.bed 可信（J=0.928）但不能 filter（TP 富集）
2. Non-LOH 也無突破（AUC<0.58）
3. ISM 對 TO FP 過濾近乎無效，轉向 characterization

**帶出三個延伸問題**：
1. LOH 加上 CNV 維度，會看到什麼？— LOH+CNV 交叉是否揭示更精確的分佈模式？
2. HP Tag 本身可靠嗎？— ISM 依賴的 HP tags 來自 LongPhase-TO，上游有問題嗎？
3. 特徵真的全部耗盡了嗎？— 除了 AUC，有沒有其他統計方法能挖到信號？

---

### LOH × CNV 交叉觀察（~6 min）

#### Slide 8：Section Divider — 「LOH × CNV 交叉觀察」

#### Slide 9：TumorLens 方向佐證

**TumorLens**（medrxiv 2026.03.18）是第一個統一 long-read 分析框架：
- 單一 assay 同時偵測 SNVs、CNVs、LOH、CpG methylation
- 使用 tumor+normal paired 模式 + purity-aware CNV/LOH modeling
- 整合 LOH、CNV、甲基化三種資訊做聯合判斷

**意義**：競爭者也走向 LOH+CNV+甲基化整合路線，驗證了我們的觀察方向正確。我們的 Phase 2 Dual-BAM 架構與此高度一致。

#### Slide 10：雙面觀察 — TP 富集與 FP 熱點

上週已知 LOH FP rate=0.239，代表 LOH 區域是 TP 富集區。本週加入 SEQC2 CNV 維度後，看到更完整的圖景：

**觀察一：LOH-only zone TP enrichment**
- LOH-only TP rate=96.1%（最高）
- LOH FP rate=0.239 vs Non-LOH 0.338 — LOH 區域的 TP 比例異常高
- 生物學解釋：LOH 是腫瘤抑制基因失活的常見機制，區域內 variant 多為真陽性

**觀察二：Gain+LOH zone FP hotspot**
- Gain+LOH zone FP rate=10.2%（最高）
- CN=3（copy gain + LOH）FP rate 達 12.9% — FP 特別集中
- Allele imbalance 在 CN gain 環境下誤導 somatic caller

**兩者合一**：LOH 區域有 TP 富集（好事），但特定 LOH+CNV 組合也有 FP 集中（壞事）。這意味著 LOH 與 CNV 的分割判斷、或作為分析參數加入，是必要的。

#### Slide 11：CN=3 根因 + 跨樣本驗證

- **CN=3 是 FP 最高點**：FP rate 12.9%（所有 CN 值最高），FP rate 隨 CN 增加反而下降
- **根因**：CN=3+LOH 環境造成 allele imbalance，somatic caller 誤判
- **跨 7 樣本不一致**：HCC1395 Gain+LOH AUC=0.782（AlleleDelta），但跨樣本 mean AUC ≤ 0.641
- **Coverage_Multiple 可信**：與 SEQC2 真實 CN Pearson r=0.831，ISM read-depth proxy 可靠

#### Slide 12：CNV zone 結論

- **Zone-aware filter 不可行**：所有排除策略 trade-off 低於 break-even（排除 CN_Loss 移除 45% FP 但損失 11% TP）
- **Simpson's Paradox 排除**：CNV 非 Quality_Score pooling 問題根源（分層後 QS diff=-0.042）
- **但揭示了 FP 分佈的生物學根因**：不同 CNV zone 有截然不同的 TP/FP 比例

**結論**：CNV zone-aware filter 關閉。但 LOH+CNV 作為分析參數（而非 filter）是必要的 — 這正是 Phase 2 Dual-BAM 架構的設計方向。

---

## Act II：HP Tag 修正 + 甲基化特徵最終驗證（~13 min, Slides 13-23）

### Self-Phasing 修正與 V5（~6 min）

#### Slide 13：Section Divider — 「HP Tag 可靠性：修正與驗證」

#### Slide 14：Self-Phasing 問題

ISM 的所有 HP-dependent 特徵（佔信號的 ~95%）依賴 LongPhase-TO 輸出的 HP tags。但我們發現了根本性問題：

**Self-Phasing Circular Dependency**：
- Somatic variant 同時作為 phasing anchor 和被評估對象
- ALT reads 系統性偏向一個 haplotype → HP1/HP2 ratio 偏差 **17.3:1**
- 62.0% TO TP LOH 移除 self-phasing 後消失（AF 0.1-0.8: ~100%）
- 同位點 HP_Ratio 跨 TO/Paired 模式 r ≈ 0 — 完全不相關

**影響量化**：
- 31.2% TO LOH 是 self-phasing artifact
- 89.9% TP somatic 有 PS assignment，67.6% phase blocks 被汙染
- Cohen's d = -1.20（flip loci）

#### Slide 15：PON-only 修正

**解法**：從 phasing VCF 移除所有 somatic variants，僅使用 PON（Panel of Normals）germline 位點作為 phasing anchor。

**修正效果**：
- HP1/HP2 bias: 17.3:1 → **消除**
- Phase block N50: +99.7%（翻倍）
- Phased rate: +23.6 percentage points
- LOH.bed: **不變**（LOH.bed 使用 VCF AF/VAF，非 BAM HP tags）

**重要澄清**：LOH.bed 生成機制使用 VCF AF（`PhasingGraph.cpp:1817`，VAF≥0.8 → HOM），而 ISM hp_ratio 使用 BAM HP tags（`ReadParser.cpp:123`）。兩套系統使用不同數據源，解釋了 Jaccard=1.0（LOH.bed 不變）與 62% ISM LOH 消失的表面矛盾。

#### Slide 16：V5 版本演進

從 Baseline 到 V5 的 `getVote()` 修正歷程：

| 版本 | 修正內容 | SEQC2 F1 | AMB% |
|------|---------|---------|------|
| Baseline | Self-phasing (circular) | 0.7157 | N/A |
| V2b | PON-only, 移除 somatic bias | 0.7115 | ~0% |
| V3-Fixed | getVote 兩層分離 | 0.7153 | 17.5% |
| V4 | altHP UNDEFINED guard | 0.7153 | 17.5% |
| **V5** | **Layer 1.5 somatic fallback** | **0.7154** | **8.0%** |

V3-Fixed 矯枉過正：germline votes=0 時忽略 somatic directional evidence → AMB% 過高。V5 加入 Layer 1.5：用 somatic HP1_1 vs HP2_1 作為 fallback，而非一律返回 HP:i:33。

#### Slide 17：V5 = 接近 Paired

**V5 全基因組驗證**：
- 129,482 reads 從 HP:i:33 重新分配到 HP:i:11/21（-54%）
- Germline tagged reads **零變化**（17,435,063 不動）
- 新分配 accuracy: **95.0%**（267/281 in clean blocks）

**方法學推論鏈**：
- V5 somatic fallback concordance = **85.4%**
- Germline concordance = **84.8%**
- Somatic directional concordance = **84.8%**
- 三者幾乎相等 → V5 somatic fallback 的準確度等同 germline phasing，方法學合理

**F1=0.7154** 與 V3-Fixed 持平 — 修正 AMB% 沒有犧牲判別力。34 個 problem PS blocks（germline acc < 70%）是 TO phasing 的根本限制，非 V5 造成。

---

### ISM 甲基化特徵最終驗證（~7 min）

#### Slide 18：Section Divider — 「ISM 甲基化特徵最終驗證」

#### Slide 19：Beyond-AUC 七方法綜合驗證

上週用 AUC 判定特徵耗盡，本週用七種獨立統計方法交叉驗證：

1. **AUC-ROC**：25 features × 748K regions, pure methylation 全 ≤ 0.58
2. **Cohen's d**：甲基化特徵 |d| < 0.20（小效應量）
3. **Mutual Information**：MI < 0.005 bits（近乎零信息）
4. **Calibration 分析**：predicted probability 幾乎不隨 feature 變化
5. **Partial Dependence**：RF model 對甲基化特徵的 marginal effect 平坦
6. **SHAP**：甲基化特徵 importance 被 HP-dependent 特徵完全壓倒
7. **Permutation Importance**：pure methylation 特徵 shuffle 後 performance 不降

**結論**：七種方法全部一致 — 純甲基化特徵對 TP/FP 判別**沒有信號**。特徵空間確認耗盡。

#### Slide 20：文獻交叉驗證

四個文獻假說的系統性檢驗：

- **L1 Directional ASM**：TP 位點 directional ratio=49.4%（與隨機 50% 無異），FP 反而有方向性（p=1.16×10⁻¹⁴）→ 與文獻預期相反，TP 的 ASM 無一致方向
- **L2 PMD Proxy**：PMD proxy 與 hypomethylation 相關確認，但判別 AUC≤0.622，排除 AF 後 ≤0.58 → 確認現象但無用
- **L4 fCpG Density**：TP vs FP fCpG density p=0.77 → 完全無差異
- **L3 Normal Baseline** ✓：7 Normal BAM 全有 MM/ML tags，可用於建立 reference → **唯一正面結果**

#### Slide 21：特徵耗盡根因

**為什麼純甲基化無法區分 TP/FP？**

核心原因是 **Germline ASM >> Somatic passenger ASM**：
- Germline ASM 效應量是 somatic passenger 的 3-6 倍
- ISM 觀察到的甲基化差異被 germline ASM 主導
- Somatic passenger variants 的甲基化 context 與 germline 無法區分

這不是特徵設計問題、不是統計方法問題、不是 confound — 是**生物學限制**。Read-level 甲基化在 TP/FP 之間本質上沒有差異信號。

ISM 的價值在於 **characterization**（描述 somatic heterogeneity 模式）而非 **filtering**（過濾 FP）。

#### Slide 22：GC + CN + PON + H2009 排除

本週同時排除了四個「也許還有其他原因」的假說：

- **GC 校正不需要**：delta-r = -0.0002（門檻 ≥ 0.03），ONT 5kHz GC bias 極小
- **甲基化對 CN 無感**：HP-free 特徵 residualized |r| < 0.07，甲基化無法感知拷貝數
- **PON 跨樣本穩定**：7 樣本 mean 97.77%，refined > 98%，穩定度升級 3→4
- **H2009 負向非 ISM 問題**：Paired FP rate 僅 0.06%（86/132,994），改進天花板 +0.0004 F1 — caller 已近完美

這四項排除共同說明：不是 GC bias、不是 CN confound、不是 PON 不穩、不是某個樣本特殊 — **甲基化特徵耗盡是全域性、根本性的**。

#### Slide 23：正面觀察

在「耗盡」的大圖景中，仍有明確的正面發現：

- **HPFineNGroups 確認為 somatic heterogeneity marker**：N≥4 + NR≥80 → TP rate 89.1%，7/7 樣本一致，residualized AUC=0.617
- **ASM 32-66% SNV 位點有等位基因特異甲基化**：ISM PERMANOVA 是唯一能量化 haplotype 間 entropy imbalance 的工具
- **FP ASM 有方向性**（directional ratio 偏離 50%，p=1.16×10⁻¹⁴）：提示 FP 附近有系統性 methylation context

這些正面觀察支持 ISM 在 **epigenetic characterization** 上的獨特價值。

---

## Act III：Phase 2 轉向（~6 min, Slides 24-30）

#### Slide 24：Section Divider — 「Phase 2：Dual-BAM 架構驗證」

#### Slide 25：Phase B/C/D Dual-BAM 驗證

Phase 2 架構核心：同時處理 Tumor BAM + Normal BAM，利用 Normal 作為 baseline。

**HCC1395 全基因體驗證通過**：

| 模組 | 結果 | 判定 |
|------|------|------|
| Phase B: Sample ASM | 97.3% 位點顯著 | ✓ PASS |
| Phase C: Normal Baseline | 100% valid entries | ✓ PASS |
| Phase D: LOH Concordance | 98.4% 一致 | ✓ PASS |
| Phase D: Subclone Analysis | 4 groups detected | ✓ PASS |

Phase A（HP ASM）+ Phase B/C/D 的程式碼已完成（2026-04-13）。

#### Slide 26：Subclone 4-group 發現

Cross-region subclone analysis 在 HCC1395 上發現四個群組：

| 群組 | 比例 | 特徵 |
|------|------|------|
| Normal Diploid | 9.2% | LOH=0%, Normal baseline 一致 |
| Epi Heterogeneous | 8.3% | 高 entropy, 多等位 ASM |
| LOH Influenced | 5.4% | LOH region, 單 haplotype |
| Tumor-Specific | 77.1% | somatic 驅動的甲基化變化 |

這是 ISM 首次在全基因體尺度上自動偵測到有生物學意義的亞克隆分群。

#### Slide 27：L3 Normal Baseline 可行性

文獻假說 L3 是本週唯一正面結果：

- **7 Normal BAM 全有 MM/ML tags** — 甲基化數據可用
- FP ASM directional 1.7× > TP — Normal baseline 可望放大此差異
- **Phase 2A 路徑確認**：Normal methylation reference → 計算 Tumor-Normal 甲基化差異 → 利用差異作為新特徵

Phase B/C/D 已驗證 Dual-BAM 架構可行，L3 驗證 Normal baseline 數據可用，Phase 2A 的執行條件已齊全。

#### Slide 28：結論總表

本週 10 項判定一覽：

**確定性關閉（5 項）**：
1. Beyond-AUC 七方法 → 甲基化特徵空間 EXHAUSTED
2. SEQC2 CNV 分層 → CNV zone-aware filter 關閉
3. GC 校正 + 甲基化-CN → 全 NO-GO
4. 文獻 L1/L2/L4 → NEGATIVE
5. Fine-Pairwise + O9 FN + TO-pure + Option C → identifiability problem 匯聚

**正面成果（3 項）**：
1. Phase B/C/D Dual-BAM 全基因體驗證通過
2. L3 Normal Baseline 可行性確認
3. PON 跨樣本穩定度 3→4

**獨立診斷（2 項）**：
1. H2009 負向根因確認（caller 已完美）
2. LOH.bed 生成機制確認（VCF AF vs BAM HP 兩套系統）

#### Slide 29：下一步

```
V5 + PON-only  ──→  7 samples × 2 modes  ──→  Phase 2A
全量重跑                ISM 重跑              Normal Baseline
    │                      │                      │
    └──────────────────────┴──────────────────────┘
                           │
                           ▼
                   7 samples 全量驗證
                           │
                           ▼
                      論文撰寫
```

1. **短期**：V5 binary + PON-only phasing → haplotag + ISM 全量重跑（7 samples × paired + TO）
2. **中期**：Phase 2A Normal Methylation Reference baseline → 利用 Tumor-Normal 差異
3. **長期**：7 samples 全量驗證 → Platform normalization → 論文

#### Slide 30：附錄/路徑

所有報告路徑索引（供 PI 自行查閱）。

---

## 觀眾心理追蹤

| Act | 觀眾應該在想... | 我們回答了... |
|-----|-----------------|-------------|
| I-前半 | 「上週說 LOH 分析做完了，那現在呢？」 | 上週結論壓縮回顧，帶出三個延伸問題 |
| I-後半 | 「LOH+CNV 一起看會怎樣？」 | TP 富集 + FP 熱點雙面觀察 → 分割判斷必要 |
| II-前半 | 「ISM 的 HP tags 本身可靠嗎？」 | Self-Phasing 問題發現 → PON-only 修正 → V5 驗證 |
| II-後半 | 「特徵真的全部耗盡了？還有其他方法嗎？」 | 七方法+文獻+GC/CN/PON 全面排除 → 確認耗盡 |
| III | 「耗盡了那怎麼辦？」 | Phase 2 Dual-BAM 已有初步驗證 → 方向確認 |

---

## 時間分配

| Section | Slides | 時間 | 核心敘事 |
|---------|--------|------|---------|
| A 封面+導覽 | 1-2 | 1 min | 本週重點卡片 |
| B 上週回顧 | 3-7 | 5 min | 壓縮七步推論鏈 → 三個延伸問題 |
| C LOH×CNV | 8-12 | 6 min | TumorLens 佐證 → 雙面觀察 → zone 結論 |
| D Self-Phasing+V5 | 13-17 | 6 min | 問題→修正→V5 接近 Paired |
| E ISM 現況 | 18-23 | 7 min | 七方法耗盡 → 文獻 → 排除 → 正面觀察 |
| F Phase 2 | 24-27 | 4 min | Dual-BAM 驗證 → 4-group → Normal baseline |
| G 結論 | 28-30 | 1 min | 總表 → 下一步 → 路徑 |
| **Total** | **30** | **30 min** | |
