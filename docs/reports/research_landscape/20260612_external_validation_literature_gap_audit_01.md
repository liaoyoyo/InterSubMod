<!--
建立時間: 2026-06-12
目標: 全面核對 /big7_disk/liaoyoyo2001/external_validation/ 的論文覆蓋，補查與 ISM 及目前論文主軸直接相關但尚未納入的來源
任務類型: B Comprehensive validation
處理範圍: external_validation 28 張 CONTEXT 卡 + InterSubMod KB/研究主軸文件 + Europe PMC 截至 2026-06-12 的主題式與精確查詢
關聯檔案:
  - /big7_disk/liaoyoyo2001/external_validation/REGISTRY.md
  - /big7_disk/liaoyoyo2001/external_validation/_landscape/
  - InterSubMod/docs/reports/research_landscape/20260611_subclonal_reconstruction_paper_foundation_01.md
  - InterSubMod/knowledge/11_external_literature/
限制: 新發現來源先完成一手 metadata/abstract 驗證；正式投稿引用前仍需逐篇全文親讀與 citation verification
-->

# external_validation 文獻缺口全面審查

用 SCQA + evidence-gap matrix：

**TL;DR：確認 external_validation 尚未完整覆蓋 ISM 論文主軸；目前 28 張卡品質高，但 axis3 為 0、axis6 僅 1，且至少 9 篇會改變核心 framing 的 P0 論文尚未建卡或未整合（影響：高，信心：高）。**

## 1. 最終裁決

### Situation

目前論文主軸是：

> **Subclonal reconstruction using somatic haplotagging and methylation profiles with Nanopore sequencing**

其中 somatic haplotagging 是重建骨幹，甲基是有界的 characterization / assist 訊號；內部結果已明確指出甲基在 germline-haplotype 層級強、within-haplotype subclone 層級弱，且不可當 variant TP/FP filter。

### Complication

`external_validation` 已由 24 張擴充到 **28 張 CONTEXT 卡**，但 coverage 與 SoT 尚未同步：

| 軸 | 現有卡數 | 判定 |
|---|---:|---|
| axis1 subclonal reconstruction | 9 | genetic/CCF 骨幹充足，但缺 methylation-based phylogeny canonical methods |
| axis2 long-read phasing/haplotag | 7 | cancer phasing 有覆蓋，但缺 NanoMethPhase、HapBridge、NANOME、base LongPhase 獨立卡 |
| axis3 cancer ASM | **0** | **結構性缺口**；核心 cis/ASM/CN-confound 文獻未建卡 |
| axis4 methylation tumor evolution | 7 | 已有 EVOFLUx/MethylTree/Epiclomal，但缺多個直接重建與 ONT epiallele 方法 |
| axis5 methylation clustering/distance | 4 | 已有最相近工具，但仍缺 DAMEfinder、pycoMeth、PoreMeth2 等定位來源 |
| axis6 somatic callers | **1** | **不足**；只有 DeepSomatic，缺實際 upstream ClairS-TO 與最新整合/benchmark |

額外 QA 問題：

1. `[L5-local]` `_landscape/05_context_cards_index.md` 標題仍寫「24」，正文卻已寫「28」。
2. `[L5-local]` `REGISTRY.md` 仍多處寫「24/24、24 卡、24 源」，已落後實際目錄。
3. `[L5-local]` 新卡中存在 `last_verified: 2026-06-13`，但本次審查日期為 **2026-06-12**；應校正日期或記錄時區/建立來源，否則 provenance 出現未來日期。
4. `[L5-local]` 「24 個源中 0 個真 CONFLICT」對原 24 源可成立，但對完整論文 framing 已不夠精確：新補來源未必是同口徑實驗衝突，卻形成明確的 **framing tension**。

### Answer

`external_validation` **不是不可靠，而是不完整且整合落後**。最需要修正的不是增加泛背景文獻，而是補齊會約束以下四個主張的來源：

1. **什麼才可稱 subclonal reconstruction**：需 clone 數、樹、branch、CN/purity deconvolution 或正交 ground truth。
2. **read-level / distance-based methylation 並非天然新穎**：已有多個直接方法。
3. **6 個 cell line 是 benchmark replication，不等於病人 cohort 泛化**。
4. **ISM 可防守 novelty** 應縮在 somatic haplotagging-conditioned、normal-anchored cis test、LOH/CN-aware 與 5mC/5hmC 有界分析，而非廣泛「用甲基重建 subclone」。

## 2. 驗證方法與範圍

### 2.1 本地盤點

輸入路徑：

- `/big7_disk/liaoyoyo2001/external_validation/`
- `InterSubMod/docs/CURRENT_FOCUS.md`
- `InterSubMod/docs/reports/research_landscape/20260611_subclonal_reconstruction_paper_foundation_01.md`
- `InterSubMod/knowledge/11_external_literature/`

執行命令：

```bash
find /big7_disk/liaoyoyo2001/external_validation -name CONTEXT.md -printf '%h\n' | sort
for d in /big7_disk/liaoyoyo2001/external_validation/axis*; do
  find "$d" -mindepth 1 -maxdepth 1 -type d | wc -l
done
rg -n "BitPhylogeny|CAMDAC|NanoMethPhase|methclone|ClairS-TO|..." \
  /big7_disk/liaoyoyo2001/external_validation InterSubMod/knowledge/11_external_literature
```

實際輸出摘要：

- CONTEXT 卡總數：**28**
- 各軸：`9 / 7 / 0 / 7 / 4 / 1`
- 沒有 context card 的論文層目錄：0
- 代表所有現存論文目錄都有卡，但**重要文獻根本尚未成為目錄/卡**。

### 2.2 網路搜尋

來源：Europe PMC REST API，一手 metadata + abstract；查詢截止日 **2026-06-12**。

主題式查詢軸：

1. `(nanopore OR long-read) AND methylation AND cancer/tumor AND subclone/evolution/phylogeny/lineage`
2. `methylation AND subclonal reconstruction/tumor phylogeny/lineage tracing`
3. `(nanopore OR long-read) AND methylation AND phasing/haplotype/allele-specific`
4. `(nanopore OR long-read) AND somatic AND benchmark/tumor-only/phasing`

另對 P0/P1 候選以題名或 DOI 精確查詢。新候選的存在、題名、日期、DOI/PMID/PMCID 與摘要主張已驗；尚未逐篇全文/code 親讀者不得當成 source-verified 細節引用。

### 2.3 證據標記

| 標記 | 意義 |
|---|---|
| L5-local | 本地全文/源碼/context card 或 repo SoT 已驗 |
| L4-primary | peer-reviewed 一手摘要與識別碼已驗 |
| L3-preprint | preprint 一手摘要與 DOI 已驗 |
| L2-review | review/secondary，只作導航，不承載核心 claim |

## 3. P0：論文 framing 凍結前必補

以下來源若不納入，reviewer 可合理質疑「reconstruction 定義、方法新穎性、ground truth、CN confound 或泛化」。

| 論文 | 證據 | 為何是缺口 | 對 ISM 主張的直接影響 | 建議軸 |
|---|---|---|---|---|
| [BitPhylogeny: a probabilistic framework for reconstructing intra-tumor phylogenies](https://doi.org/10.1186/s13059-015-0592-6) (2015) | L4-primary | KB 已知，但 external_validation 無卡 | full Bayesian 聯合推 clone 數、組成與樹；直接界定「reconstruction」高於 flat clustering | axis1/axis4 bridge |
| [Epigenetic evolution and lineage histories of CLL](https://europepmc.org/article/MED/31092926) (Gaiti et al., Nature 2019) | L4-primary | canonical single-cell methylation lineage reconstruction 未建卡 | genetic subclones 可映射到僅由 epimutation 推得的 clades；是 ISM 缺少的單細胞正交 ground truth | axis4 |
| [Sgootr: Single-cell methylation sequencing data reveal ... tumor progression models](https://doi.org/10.1101/gr.277608.122) (2023) | L4-primary | 網路搜尋確認，external_validation/KB 主索引未充分 foreground | **distance-based** methylation lineage tree，且用真腫瘤單細胞資料；直接縮限「距離式甲基重建」新穎性 | axis4/axis5 |
| [DNA methylation cooperates with genomic alterations during NSCLC evolution](https://europepmc.org/article/MED/40931149) (CAMDAC/TRACERx, Nat Genet 2025) | L4-primary | KB 已知，但 external_validation 無卡 | 217 tumor/normal regions、59 patients、顯式 methylation deconvolution 與 genomic evolution；要求 ISM 正面處理 CN/purity | axis4/axis3 |
| [Long-read sequencing of single cell-derived melanoma subclones](https://doi.org/10.1101/2025.08.28.672865) (2025) | L3-preprint | KB 已知，但 external_validation 無卡 | 23 單細胞衍生 subclone，聯合 SNV/SV/CNV/甲基並有 phylogeny-constrained calls；是最接近 ISM 的 orthogonal comparator | axis4/axis6 bridge |
| [PoreMeth2 for decoding the evolution of methylome alterations with nanopore sequencing](https://europepmc.org/article/MED/41115805) (2025) | L4-primary | 本次網路搜尋新確認，external_validation 無卡 | ONT epiallele diversity + methylation-frequency DMR；會壓縮「long-read methylation pattern/evolution」泛用 novelty | axis4/axis5 |
| [ClairS-TO](https://europepmc.org/article/MED/41173866) (Nat Commun 2025) | L4-primary | 實際 upstream caller，registry 有名但無卡；axis6 僅 DeepSomatic | HCC1395/COLO829、tumor-only、purity/VAF/coverage benchmark；必須交代 ISM input provenance 與單 pipeline 自參照 | axis6 |
| [Cancer genome standards for long-read sequencing using cancer cell line mixtures](https://europepmc.org/article/MED/41934171) (GigaScience 2026) | L4-primary | 本次網路搜尋新確認 | 10 purities、22 samples、60x；matched normal depth 對準確度關鍵。直接校準 ISM 跨樣本與 purity/depth 設計 | axis6 |
| [Substantial genomic and methylation variability between MCF-7 sublines](https://doi.org/10.64898/2026.02.17.706500) (2026) | L3-preprint | 本次網路搜尋新確認 | 證明同一 cancer cell line 的 subline 可有顯著 genomic/ASM/methylation 差異；限制「6 cell lines = biological generalization」 | axis3/limitations |

### P0 核心結論

`BitPhylogeny + Gaiti + Sgootr + CAMDAC + Lee2025` 形成完整壓力鏈：

```text
正式 reconstruction
= clone identity/number
+ lineage tree / evolutionary model
+ CN/purity handling
+ orthogonal or single-cell ground truth
```

ISM 目前已有 somatic haplotagging 骨幹與 read-level methylation structure，但尚缺上述完整鏈。除非 G-D 真正重建 demo 與 G-E 正交驗證完成，正文應優先使用：

> **characterize methylation-associated subclonal structure conditioned on somatic haplotags**

而非無限制地宣稱：

> **reconstruct subclones using methylation**

## 4. P1：Related Works / Methods 必補

這些來源多數已在 KB 被認識，但尚未落成 external_validation CONTEXT 卡；問題是「整合缺口」，不是完全沒搜尋過。

| 群組 | 應補來源 | 對 ISM 的用途 |
|---|---|---|
| methylation disorder / clonal dynamics | [methclone](https://europepmc.org/article/MED/25260792)、[Landau PDR](https://europepmc.org/article/MED/25490447) | 把 ISM structure 與既有 disorder/entropy/epiallele-shift 清楚分開 |
| methylation-assisted phasing | [NanoMethPhase](https://europepmc.org/article/MED/33618748)、[HapBridge](https://doi.org/10.1101/2025.11.07.687303)、[NANOME](https://doi.org/10.1101/2025.06.29.662079)、base [LongPhase](https://doi.org/10.1093/bioinformatics/btac058) | 證明 germline methyl-phasing 已擁擠；ISM novelty 只能落在 tumor/somatic/LOH extension 與有界評估 |
| cancer ASM / cis / CN confound | [DAMEfinder](https://europepmc.org/article/MED/32487212)、[Do et al. 2020](https://europepmc.org/article/MED/32594908)、[Martin-Trujillo et al. 2017](https://europepmc.org/article/MED/28883545) | axis3 必建基本盤；支撐 normal-anchored cis test，也防止把 CN/subclone drift 誤寫成 cis |
| latest evolution models | [Fluctuating DNA methylation sites encode colorectal tumour growth history](https://doi.org/10.64898/2026.06.04.730217) (2026-06-09) | 最新多區域、機制式 growth inference；再次說明「evolution」主張需模型與參數推論 |
| integrated long-read somatic workflows | [LRSomatic](https://doi.org/10.64898/2026.02.26.707772) (2026) | 同時支援 SNV/indel/SV/CN、paired/TO 與 epigenetic integration；可作 G-E 第二 pipeline 候選 |

## 5. P2：應追蹤，但不承載核心 novelty

| 來源 | 用途 | 為何不是 P0/P1 |
|---|---|---|
| [Read-level DNA methylation deconvolution enhances ctDNA detection](https://europepmc.org/article/MED/41115210) (Alpha, 2025) | read-level deconvolution 方法對照 | cfDNA/tumor-fraction，不是 tissue subclone tree |
| [MethylBench](https://doi.org/10.64898/2026.04.28.721268) (2026) | ONT/短讀/array 跨平台甲基 benchmark；提醒 coverage filtering | preprint，且偏平台比較 |
| [Clonal tracing with somatic epimutations reveals dynamics of blood ageing](https://europepmc.org/article/MED/40399669) (EPI-Clone, Nature 2025) | 甲基 epimutation 可作大規模 clonal barcode 的高階背景 | 非癌症重建主問題 |
| [Characterization of subclonal variants in HG002 GIAB](https://europepmc.org/article/MED/41421359) (2025/2026) | 85 個 orthogonally curated subclonal SNV，可作低 AF caller 驗證資源 | 非癌 cell line、非甲基 |
| [Benchmarking major somatic SV callers on HG008](https://europepmc.org/article/MED/42200198) (2026) | 多 caller ensemble 與 orthogonal SV benchmark | 與目前 small-variant haplotag 主線較遠 |
| [Methylation-based lineage tracing in cancer](https://europepmc.org/article/MED/41397346) (Blood 2026) | 最新 review，投稿前 bibliography completeness check | review，不應承載方法 claim |

## 6. 對目前論文主張的影響

### 6.1 「Subclonal reconstruction」

**裁決：需縮限或補實驗，不能只靠現有 flat clustering。**

- BitPhylogeny：推 clone number/composition/tree。
- Gaiti/Sgootr：單細胞甲基可直接重建 lineage tree。
- EVOFLUx/2026 CRC growth model：bulk 甲基也能推演化參數。
- CAMDAC：正式處理 CN/purity/deconvolution。
- Lee2025：長讀 + 單細胞衍生 subclone + phylogeny constrained ground truth。

因此，外部文獻不是直接 refute ISM 結果，而是 refute **過寬動詞**。

### 6.2 「read-level methylation pattern 是主要新穎點」

**裁決：不可單獨作 novelty。**

現有卡中的 ASMS、cvlr、MethylBERT，加上待補的 PoreMeth2、DAMEfinder、pycoMeth、NANOME，已涵蓋 read clustering、read-level pattern、ASM、epiallele diversity、haplotype-aware DMR 與 deconvolution。

可防守差異應寫成組合：

> somatic haplotag-conditioned read structure  
> + matched-normal/normal-anchored somatic cis test  
> + LOH/CN-aware interpretation  
> + bounded 5mC/5hmC characterization

### 6.3 「甲基輔助 tumor/LOH phasing 是白地」

**裁決：仍可防守，但必須 head-to-head。**

MethPhaser、NanoMethPhase、HapBridge、NANOME 都主要是 germline；目前搜尋仍未找到完成 tumor/LOH methyl-assisted phasing 的直接同軸論文。這是 ISM 最清楚的開放 niche。

但要成立為 paper claim，至少需：

1. 與 MethPhaser/HapBridge/LongHap 在相同 tumor 資料比較 switch error、block N50、read-tag recovery。
2. LOH 內有 orthogonal phase truth。
3. 多樣本、第二 caller/pipeline 驗證。

### 6.4 「6 cell lines × 3 cancers」

**裁決：可稱跨模型 reproducibility，不可稱 clinical/general population generalization。**

- Cancer genome standards 2026 支持 cell-line mixture、purity 與 depth benchmark 的合理性。
- MCF-7 subline 2026 顯示同一 cell line 仍可能因 subline/provenance 產生大幅 genomic 與 methylation divergence。

因此 Methods 必須記錄 cell-line source、subline、basecaller、library、passage/provenance；Discussion 必須把樣本數與病人 cohort 泛化分開。

### 6.5 「0 個真衝突」

建議改寫為：

> 尚未發現與 ISM 有界內部結果在同資料 regime 下直接相反的實驗結論；但存在多個高可信度 framing tensions，限制 reconstruction、novelty 與 generalization 的可宣稱範圍。

這比「0 conflict」更符合完整文獻地景。

## 7. external_validation 建議補卡與 SoT 修復順序

### 第一批：P0 context cards

1. axis1/axis4 bridge：BitPhylogeny、Gaiti2019、Sgootr2023、CAMDAC/TRACERx2025、Lee2025。
2. axis4/axis5：PoreMeth2。
3. axis6：ClairS-TO、Cancer genome standards 2026。
4. axis3/limitations：MCF-7 subline variability 2026。

### 第二批：補空軸與方法源流

1. axis3：Do2020、Martin-Trujillo2017、DAMEfinder、POG/O'Neill cross-link。
2. axis2：NanoMethPhase、HapBridge、NANOME、base LongPhase。
3. axis4/axis5：methclone、Landau、2026 colorectal growth-history model。
4. axis6：LRSomatic；之後再評估 HG002/HG008 benchmark cards。

### 同步修復

1. 更新 `REGISTRY.md` 的 24/24、repo/PDF/card 數。
2. 更新 `_landscape/01_verified_sources.md`、`02_consensus_contradiction_map.md`、`04_gaps_unverified_followups.md`、`05_context_cards_index.md`。
3. 將 `axis3_cancer_asm` 從空軸補成正式 evidence chain。
4. 校正本地 `2026-06-13` future-date provenance。
5. 每張新卡先標 `primary_abstract_verified`，全文/code 親讀後才升級。

## 8. 對研究 gap 的直接映射

| 內部 gap | 外部來源要求 | 優先動作 |
|---|---|---|
| G-A 多樣本驗證 | Cancer genome standards、MCF-7 subline variability | 報告 purity/depth/source/subline，避免把 cell-line 數量當 cohort 泛化 |
| G-B within-haplotype null | DAMEfinder、Do2020、Martin-Trujillo、CAMDAC | 建立 somatic-vs-baseline、CN/purity-aware null |
| G-C cis vs program | CAMDAC、POG/O'Neill、DAMEfinder | 分離 focal cis、subclone program、CN/LOH |
| G-D 真 reconstruction demo | BitPhylogeny、Gaiti、Sgootr、EVOFLUx、Lee2025 | 至少輸出 clone/群結構與可驗的 lineage/ground-truth 對照 |
| G-E orthogonal pipeline | ClairS-TO、DeepSomatic、LRSomatic、GIAB/HG008 | 第二 caller/第二 workflow，打破 longphase-S 自參照 |

## 9. 最終建議

### 可維持

- somatic haplotagging 是 genetic reconstruction backbone。
- 甲基可提供有界、read-level、haplotype-linked characterization。
- tumor/LOH methyl-assisted phasing 仍是合理開放 niche。
- 「甲基不是 variant filter」與外部來源無直接衝突。

### 必須調整

- 在完成 G-D/G-E 前，避免把現有 flat clustering 寫成已完成的 full subclonal reconstruction。
- 不把 read-level methylation、distance、clustering 或 methylation evolution 單獨宣稱為新穎。
- 不把 6 cell lines 寫成 patient/cohort generalization。
- 將「0 conflict」改成「無同 regime 直接衝突，但有高影響 framing tensions」。

### 信心

- **高信心**：現有 28 卡數量、各軸覆蓋、axis3/axis6 缺口、SoT 過期、P0 論文存在與摘要層關聯。
- **中高信心**：P0/P1 優先排序及其對 framing 的影響。
- **中信心**：新候選與 ISM 的精確方法差異；需全文/code 親讀後才能升為正式投稿級比較。

---

## 來源入口

- [Europe PMC](https://europepmc.org/)
- [BitPhylogeny](https://doi.org/10.1186/s13059-015-0592-6)
- [Gaiti et al. 2019](https://europepmc.org/article/MED/31092926)
- [Sgootr 2023](https://doi.org/10.1101/gr.277608.122)
- [CAMDAC/TRACERx 2025](https://europepmc.org/article/MED/40931149)
- [PoreMeth2 2025](https://europepmc.org/article/MED/41115805)
- [ClairS-TO 2025](https://europepmc.org/article/MED/41173866)
- [Cancer genome standards 2026](https://europepmc.org/article/MED/41934171)
- [2026 colorectal tumour growth-history preprint](https://doi.org/10.64898/2026.06.04.730217)
- [MCF-7 subline variability 2026](https://doi.org/10.64898/2026.02.17.706500)
