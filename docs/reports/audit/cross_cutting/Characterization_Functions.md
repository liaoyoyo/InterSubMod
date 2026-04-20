# Characterization Functions（主軸定位下的功能清單）

> **建立日期**: 2026-04-19

## 背景

**用戶定調（R-02，2026-04-20）**：
> 主軸是 **read-level epigenetic characterization** 定位；後續補齊 Phase 2 A+D 才有足夠故事線；同時也做 Subclone Marker Tool 的分析**作為其中的功能**。

**此文件目的**：
- 將 22 結論中「POSITIVE / CONDITIONAL POSITIVE」的結果歸類為 characterization 的**具體功能**
- 明確每個功能的 downstream 應用場景、數據依賴、輸出產物、成熟度
- 作為 Phase 2 論文故事線的功能清單基礎

**非此文件範圍**：
- 結論重新論證（詳見各 audit card）
- 論文投稿策略（詳見 R-03 決策）

---

## 核心功能樹（按 Phase 2 論文章節構想）

```
ISM read-level epigenetic characterization
├── F1: Self-Phasing 機制揭露（TO mode 診斷）
│   ├── F1.1 Phasing-based LOH 偽訊號識別（C09, C10, C21）
│   └── F1.2 LongPhase-TO bias 量化（17.3:1）
├── F2: Subclone Marker Tool（子功能）
│   ├── F2.1 HPFineNGroups somatic heterogeneity marker（C16）
│   └── F2.2 LOH Subclone AF × Methylation dual evidence（C17）
├── F3: Zone-Aware Characterization
│   ├── F3.1 Z1-Z5 confidence zone stratification（C22, bug-fix pending）
│   └── F3.2 HCC1954-like CNV reversal 診斷（C18）
├── F4: ASM Profiling（探索中）
│   ├── F4.1 ASM 32-66% 跨樣本分佈（C12）
│   └── F4.2 germline vs somatic ASM 區分（Phase 2 A+D 補齊）
└── F5: Variant Confidence Score（附加）
    ├── F5.1 LOH.bed clean Jaccard=1.0（C21）
    └── F5.2 Phase 1A ML filter（C11, F1=+0.0112 small gain）
```

---

## 功能明細

### F1: Self-Phasing 機制揭露

**描述**：ISM 最大原創貢獻 — 揭露 TO mode 的 phasing-based LOH circular dependency，並量化 LongPhase-TO 的 bias 機制。

**組成結論**：
- **C09** Self-Phasing causal chain CONFIRMED（⭐4，7/7 樣本）
- **C10** PON-only phasing 驗證（⭐4，HCC1395 only）
- **C21** LOH.bed 不受 self-phasing 汙染（⭐5，Jaccard=1.0000）

**下游應用**：
1. TO mode variant caller 開發警告：不可用 tumor-derived haplotype 做 filter
2. LOH.bed 作為「乾淨」的 LOH set（不讀 HP tag）供下游分析
3. 臨床：ONT TO mode 在 LOH-rich 癌型（HER2+、TNBC）的使用限制

**數據依賴**：VCF coordinate + LongPhase output（不需重跑）

**輸出產物**：
- `output/canonical/*/LOH.bed`（coordinate-only）
- `research/.../self_phasing_evidence_chain.md`
- `06_結論穩定性審查.md:62%` 的 LOH 消失 metric

**成熟度**：⭐⭐⭐⭐⭐（最成熟，可寫進論文 main finding）

**遺留修正**：
- **P2-A**：LOH.bed vs ISM HP_Ratio LOH 兩層概念在 CURRENT_FOCUS.md:53 / 07_LOH_CN_AF.md:102 需補註

---

### F2: Subclone Marker Tool（子功能）

**描述**：以 HPFineNGroups + LOH AF×Methylation 雙路徑識別 somatic subclone heterogeneity。

**組成結論**：
- **C16** HPFineNGroups NG=4+AF<0.4+NR≥80 filter（⭐4，7/7 樣本，F pilot 2026-04-18）
- **C17** LOH Subclone AF×NGroups r=+0.705（⭐4，7/7 p<10⁻³⁹）

**下游應用**：
1. Tumor heterogeneity profiling：per-region NG 分佈 + AF 聯合
2. Clinical subclone detection: Low-purity tumor 的 subclone emergence 追蹤
3. Longitudinal monitoring: treatment response 的 subclone shift 偵測

**數據依賴**：
- HP-tagged BAM + ISM baseline（已有）
- ⚠️ C17 受 CovM bug 影響，修正後需重算 r 值

**輸出產物**：
- `output/canonical/*/HPFineNGroups.tsv`
- `output/canonical/*/LOH_subclone_evidence.tsv`
- Methylation + AF joint plots

**成熟度**：⭐⭐⭐⭐（pending CovM bug fix + P0-B within-group OLS + P1-C bootstrap CI）

**遺留修正**：
- **P0-B**：within-group OLS 重算 r=+0.705（pooled OLS trap 可能保留分組信號）
- **P1-C**：1000× bootstrap CI for NGroups residualized AUC
- **P2-C**：NG=3 非單調現象（C-BIO-3）補充生物機制討論

---

### F3: Zone-Aware Characterization

**描述**：將 region 按 Coverage_Multiple + QS 分層為 Z1-Z5 五個 confidence zones，提供 variant 的可信區分層 characterization。

**組成結論**：
- **C22** Zone-Aware F1 NEGATIVE / Characterization CONFIRMED（⭐4/⭐3）
- **C18** HCC1954 CNV-driven reversal 雙路徑機制（⭐4）
- **C19** Z3 amplicon blacklist CONDITIONAL-NEGATIVE-for-canonical（⭐5）

**下游應用**：
1. Variant filter downstream：不用於 canonical FP removal（已證 NEGATIVE）
2. Characterization only：回答「此 variant 在哪個 confidence zone」
3. HER2+/pseudo-tetraploid 癌型診斷：HCC1954-like pattern 識別

**數據依賴**：🔴 **CovM baseline bug 修正後需重定義 Z1-Z5 邊界**

**輸出產物**：
- `research/zone_aware_framework/Z1-Z5_definitions.md`（bug-fix pending）
- Zone-stratified TP rate 表

**成熟度**：⭐⭐⭐（bug-fix gated；修正後升級）

**遺留修正**：
- **P0-A**：CovM bug 修正 + Zone 重定義
- **P1-A**：`08_Zone_Aware.md:49,73` 誤標 r=0.997 修正
- **P1-B**：Zone TP rate 差異 bootstrap + FDR
- **P2-B**：Characterization 具體功能清單化（即此文件）

---

### F4: ASM Profiling（探索中）

**描述**：ASM (Allele-Specific Methylation) 跨樣本分佈 32-66%，作為 read-level epigenetic heterogeneity 的量化 metric。

**組成結論**：
- **C12** ASM POSITIVE（⭐4，但 FP>>TP 重疊大）

**下游應用**：
1. Methylation haplotype profiling
2. Imprinting region 識別（germline ASM）
3. Tumor epigenetic evolution（somatic ASM，Phase 2 A+D 補齊）

**數據依賴**：
- HP-tagged BAM + CpG methylation
- Normal BAM（germline/somatic 區分，Phase 2）

**輸出產物**：
- `output/canonical/*/asm_*.tsv`
- 跨樣本 ASM 分佈直方圖

**成熟度**：⭐⭐⭐（germline vs somatic 未區分是最大缺口）

**遺留修正**：
- **P2-C**：germline (15-30%) vs somatic ASM 區分 + 文獻引用
- **Phase 2 A+D**：Normal BAM 整合提供 germline baseline

---

### F5: Variant Confidence Score（附加）

**描述**：基於 ISM 特徵的 variant-level confidence score，提供 Phase 1A ML filter 與 LOH.bed clean set。

**組成結論**：
- **C11** Phase 1A ML read classification（⭐3，paired F1=+0.0112）
- **C21** LOH.bed clean Jaccard=1.0（⭐5，SEQC2 F1=96.2%）

**下游應用**：
1. ClairS-TO 後端 filter（僅 paired mode，gain 微小）
2. LOH.bed 作為下游 CNV caller 的 priors
3. 不作為獨立 variant filter 推廣（F1 gain 過小）

**數據依賴**：已有 paired mode VCF + SEQC2 truth

**輸出產物**：
- `output/canonical/*/ml_filter_result.tsv`
- `output/canonical/*/LOH.bed`

**成熟度**：⭐⭐⭐（effect size 小是主要限制）

**遺留修正**：
- **P2-B**：Phase 1A F1 effect size 誠實揭露（per-sample CI overlap）
- 不進一步擴展（資源投入 Phase 2 A+D 優先）

---

## 功能間依賴關係

```
F1 Self-Phasing ──(揭露 bias)──> F5.1 LOH.bed clean
F3 Zone-Aware ──(CovM bug fix)──> F2.2 LOH Subclone (C17 step3)
F4 ASM ──(Phase 2 A+D)──> germline vs somatic 區分
F2 Subclone Marker <──(B.2 evidence)── F3.2 HCC1954 reversal (C18)
```

---

## Phase 2 論文故事線對接

| 章節構想 | 主要支撐功能 | 支援證據 |
|---------|-----------|---------|
| Intro: ONT TO mode 的 phasing 陷阱 | F1 | C09, C10, C21 |
| Methods: ISM read-level characterization framework | F3, F4 | C12, C22 |
| Results 1: Self-Phasing 揭露（最強貢獻） | F1 | C09, C21 |
| Results 2: Subclone marker（新觀察） | F2 | C16, C17 |
| Results 3: Zone-aware + CNV-driven reversal | F3 | C18, C19, C22 |
| Discussion: Characterization 定位 + F5 附屬 | F5 | C11 |
| Phase 2 A+D outlook: Normal BAM integration | F4 | ASM germline/somatic |

---

## 成熟度總表

| 功能 | 成熟度 | 主要 gate | 預期 ready 時機 |
|------|--------|----------|---------------|
| F1 Self-Phasing | ⭐⭐⭐⭐⭐ | 已完成 | 當下可寫 |
| F2 Subclone Marker | ⭐⭐⭐⭐ | P0-B + P1-C | 1-2 週 |
| F3 Zone-Aware | ⭐⭐⭐ | P0-A CovM bug fix | 2-3 週 |
| F4 ASM | ⭐⭐⭐ | Phase 2 A+D Normal BAM | 1-2 月 |
| F5 Variant Confidence | ⭐⭐⭐ | 不擴展 | 當下可寫 |

---

## 整體定位宣言

**ISM is not a variant filter. ISM is a read-level epigenetic characterization framework** that:
1. Reveals structural biases in TO mode variant calling (F1)
2. Identifies subclone heterogeneity signals (F2)
3. Stratifies variant confidence zones (F3)
4. Profiles allele-specific methylation (F4)
5. Provides complementary filters for specific use cases (F5)

**Primary contribution**: F1 self-phasing mechanism revelation.
**Secondary contribution**: F2 subclone marker tool (as a function, not the whole story).
**Tertiary contributions**: F3-F5 characterization features.
