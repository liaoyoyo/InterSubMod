<!--
建立時間: 2026-06-15
任務: 對使用者提出的 HCC1395 chr2:18.07-18.11M 五群 subclone 與推論 1-6 做獨立資料複核
權威資料:
  - data/answer/SEQC2/*
  - /big8_disk/data/HCC1395/SEQC2/*
  - HKU tumor/normal tagged BAM
  - HCC1395_DORADO paired tagged tumor BAM + raw normal BAM
重算資產:
  - 20260615_chr2_18M_subclone_verification_assets/scripts/independent_subclone_audit.py
  - 20260615_chr2_18M_subclone_verification_assets/data/independent_audit.{json,md}
狀態: independent_audit_complete
-->

# HCC1395 chr2:18M Subclone 推論獨立複核

## 研究設計與範圍

- **任務類型**：B comprehensive validation；完整驗證使用者指定的 `chr2:18,066,480-18,110,828` 區域與 HCC1395 HKU／DORADO paired data，不外推為全基因組或跨病人結論。
- **服務目標**：G3 read-level epigenetic、G4 reproducibility、G5 external-verifiable evidence；LongPhase-S 是 genetic backbone，不重開 methylation-as-variant-filter DEAD 方向。
- **研究問題**：六個候選 sSNV、LongPhase-S haplotag 與 10 個候選 CpG 是否足以支持五群 subclone 與使用者提出的推論 1-6？
- **假設**：若是真實 subclone，應有可重現的同-read mutation linkage、matched-normal ALT 缺失、跨 basecaller 一致性，以及不能被 normal ASM／homopolymer artifact 解釋的 methylation coherence。
- **成功條件**：至少一組分支事件在 HKU 與 DORADO 方向一致；局部 allele-CpG association 通過 BH-FDR；LOH 有獨立 benchmark 支持；演化順序只在觀察到 nested genotype 時成立。
- **失敗條件**：群組只由單點低品質 call、homopolymer、單一 strand、normal ASM 或缺少跨距 read 支撐；此時不得稱 biological subclone 或確認演化順序。
- **主要指標**：unique primary read 的 BQ20 allele count、MAPQ≥20/BQ≥10 linkage `00/10/01/11`、Mann-Whitney/Fisher BH-FDR、HP imbalance、matched-normal ALT、SEQC2 truth/HC/LOH status。

### 主要輸入

- Reference：`/big8_disk/ref/GRCh38_no_alt_analysis_set.fasta`
- HKU tumor：`/big8_disk/liaoyoyo2001/data/bam/HCC1395_Tmode_tagged_ClairS_pileup_v040_woFilter.bam`
- HKU normal：`/big8_disk/liaoyoyo2001/data/bam/HCC1395BL_ONT_5khz_simplex_5mCG_5hmCG_tagged.bam`
- DORADO tagged tumor：`/big7_disk/liaoyoyo2001/big7_disk_output/canonical/HCC1395_DORADO/paired_full/20260315_HCC1395_DORADO_paired_full_full_complete_matrix/longphase_s/HCC1395_DORADO_tagged.bam`
- DORADO raw tumor／normal：`/big8_disk/data/HCC1395/ONT_Dorado/HCC1395.bam`、`/big8_disk/data/HCC1395/ONT_Dorado/HCC1395BL.bam`
- Truth／HC／LOH：SEQC2 truth VCF、High-Confidence BED 與 CNV/LOH BED。

### 輸出與驗證

- 重算腳本：`InterSubMod/docs/experiments/in_progress/2026/06/20260615_chr2_18M_subclone_verification_assets/scripts/independent_subclone_audit.py`
- 完整機器可讀結果：`InterSubMod/docs/experiments/in_progress/2026/06/20260615_chr2_18M_subclone_verification_assets/data/independent_audit.json`
- 人讀摘要：`InterSubMod/docs/experiments/in_progress/2026/06/20260615_chr2_18M_subclone_verification_assets/data/independent_audit.md`
- 實際執行結果：程式退出碼 `0`，輸出最後顯示兩個 `[done]` 路徑；腳本通過 `python3 -m py_compile` 與 `git diff --check`。

## 最終裁決

此區域可防守的結論是：

> **在 SEQC2-confirmed LOH 區中，存在可重現、由多個 somatic allele linkage 定義，且帶一致局部 5mC 結構的 regional subclonal states。資料支持至少兩個主要分支，以及 `(3) A → (5) C` 的巢狀子分支；甲基可作 lineage coherence 與區域 characterization，但不能把 pos4 的 G/T/DEL 證成三個 subclone，也不能單獨確認完整演化順序。**

原始五群模型需要收斂成：

- **可支持：至少 3 個 regional molecular states**
  - `alpha`: `(3)=A`, `(5)=G`
  - `alpha-1`: `(3)=A`, `(5)=C`
  - `beta-like`: `(3)=G`, `(4)=non-ref`, `(6)=G`
- **另有一個已觀察、但未與 beta-right 直接 genetic-link 的 beta-left block**
  - `(1)=G`, `(2)=C`
  - 只能以互斥與甲基 coherence 暫接到 broad beta-like program。
- **不可支持：5 個 biological subclones**
  - `(4)` 的 `G/T/DEL` 應合併成一個 **homopolymer-uncertain pos4-altered state**。

整體證據等級：

- **LOH、somatic allele existence、局部 read linkage：強**
- **regional subclonal structure：強支持，但仍屬 operational / regional**
- **完整 clone identity、clone fraction、完整 phylogeny：未證**
- **演化順序：只有 `(3) → (5)` 有直接 nesting 支持；其餘未定**
- **甲基造成突變：完全未證**

## 一、先修正三個資料口徑

### 1. `(4.1)` 座標誤植

圖中的 `(4.1)` 是 `chr2:18,096,341`，這是 CpG 的 G coordinate；正確 CpG C coordinate 是：

```text
chr2:18,096,340 C>G dinucleotide
```

文字中的 `chr2:18,096,041` 在 GRCh38 是 `A`，不是 CpG。既有腳本因此漏掉真正 `(4.1)`；本次已重算。

### 2. `(3) chr2:18,086,020` 不能標成 SEQC2 FP

- `(1)(2)(4)(5)(6)` 都存在於 SEQC2 truth VCF。
- `(3)` 位在 HC BED 的空隙：

```text
上一段 HC: chr2:18,085,979-18,085,984
下一段 HC: chr2:18,086,058-18,086,100
```

因此 `(3)` 是 **out-of-HC / truth-unevaluable candidate**，不是可由 SEQC2 判定的 FP。

它在兩套 tumor data 都為 `G>A PASS`，且 matched normal 沒有 A：

| Dataset | Tumor BQ20 A/G | Normal BQ20 A/G |
|---|---:|---:|
| HKU | 29 / 30 | 0 / 18 |
| DORADO | 21 / 27 | 0 / 27 |

最合理說法是：

> `(3)` 為強 tumor-specific somatic candidate，但因落在 SEQC2 HC gap，不能升成 truth-confirmed TP。

### 3. Normal 並非「每一個 call 都是 REF」

兩套 normal 都沒有六個預期 somatic ALT，但仍有少數 ONT error call，例如 DEL 或 T。

正確說法：

> **matched normal 未觀察到預期 somatic ALT；不是所有 normal read 的每個 base 都絕對等於 REF。**

## 二、LOH 與 haplotag 結論

### LOH：已確認

獨立證據一致：

1. SEQC2 CNV BED：`chr2:16,146,119-22,100,000 = loh`，完整涵蓋此區域。
2. HKU tagged tumor：HP2-parent `232`、HP1-parent `4`，tagged reads 中 HP2-parent 約 `98.3%`。
3. DORADO paired tagged tumor：HP2-parent `279`、HP1-parent `1`，tagged reads 中 HP2-parent 約 `99.6%`。

因此 `[推論 1] 發生 LOH` 可升為 **confirmed within available benchmark definitions**。

### LongPhase-S fine tag 不足以分出本區域 subclone：成立

HKU 為 `HP2=130 / HP2-1=102`；DORADO 為 `HP2=208 / HP2-1=71`。HKU 的六個主要 ALT 幾乎都落在 `HP2-1`；但 DORADO 的 `(3)A` 主要落在 `HP2`，而 `(1)G/(2)C/(5)C/(6)G` 主要落在 `HP2-1`。這表示 fine tag 不只無法完整拆開 alpha、alpha-1 與 beta-like，其細分類在兩套 basecalling 資料間也不穩定，不能直接當 subclone label。

但「HP2 突變到不見／沒有抽到 HP2 read」不成立：

- HP2 是主要 parent haplotype。
- HKU tumor 有 `17` 條 reads 在至少 2 個已覆蓋 SNV 上全為 REF，至少 3 點全 REF 也有 `4` 條。
- 沒有跨 ≥4 點的 all-REF read 是 read-length/coverage 限制造成，不能解讀成 ancestral HP2 消失。

## 三、read-level linkage 實際支持什麼

以下以 unique primary reads、MAPQ≥20、base quality≥10 重算。`11` 表示兩事件同時存在於同一 read。

| 事件組合 | HKU: 00 / 10 / 01 / 11 | DORADO: 00 / 10 / 01 / 11 | 判讀 |
|---|---|---|---|
| `(1)G` vs `(2)C` | 18 / 0 / 0 / **13** | 11 / 1 / 0 / **6** | `(1)(2)` 同支強支持 |
| `(1)G` vs `(3)A` | 3 / 3 / 7 / **0** | 無足夠跨距 read | 互斥 |
| `(2)C` vs `(3)A` | 3 / 5 / 9 / **0** | 0 / 1 / 0 / **0** | 互斥 |
| `(3)A` vs `(5)C` | 8 / 1 / 0 / **4** | 1 / 2 / 0 / **2** | `(5)C` 巢狀於 `(3)A` |
| `(3)A` vs `(4)non-ref` | 0 / 9 / 13 / **0** | 0 / 7 / 5 / **1** | 主要互斥；DORADO 1 discordant |
| `(3)A` vs `(6)G` | 0 / 3 / 3 / **0** | 無足夠跨距 read | 互斥 |
| `(4)non-ref` vs `(6)G` | 9 / 0 / 0 / **10** | 1 / 0 / 0 / **1** | 同 state；順序不可分 |
| `(5)C` vs `(6)G` | 3 / 7 / 11 / **0** | 1 / 7 / 0 / **0** | alpha-1 與 beta-like 互斥 |

### 可直接推出

1. `(1)(2)` 是同一 branch-defining event set。
2. `(3)` 與 `(1)(2)`、`(4)`、`(6)` 主要互斥，支持 alpha / beta-like 分叉。
3. `(5)C` 只在有 `(3)A` 的跨位點 reads 出現，支持 `(3)A → (5)C`。
4. `(4)non-ref` 與 `(6)G` 在 HKU 跨位點 reads 完全共現，但這只能說同 state，**不能排序誰先發生**。

### 不能直接推出

1. `(1)(2)` 與 `(4)(6)` 沒有足夠直接跨 36-40 kb 的 genetic linkage read。
2. 兩者可由互斥與甲基一致性推為同一 broad beta-like program，但仍需標為 **methylation-bridged inference**。
3. 所有六點的單一路徑完整 read 幾乎不存在，不能畫成已確認完整 phylogeny。

## 四、甲基關聯：哪些實際跨 basecaller 複製

下表是每個局部變異錨點的 ALT vs REF 5mC mean。q 為 10 個候選 CpG 內 Mann-Whitney BH-FDR。

| Variant anchor | CpG | HKU ALT / REF, q | DORADO ALT / REF, q | 結論 |
|---|---|---|---|---|
| `(2) G>C` | 2.1 | 0.141 / 0.719, `7.2e-4` | 0.012 / 0.858, `1.5e-5` | 複製 |
| `(2) G>C` | 2.2 | 0.112 / 0.843, `2.3e-4` | 0.241 / 1.000, `1.4e-3` | 複製 |
| `(3) G>A` | 3.1 | 0.222 / 0.917, `3.3e-4` | 0.074 / 0.998, `2.4e-5` | 複製 |
| `(3) G>A` | 3.2 | 0.970 / 0.144, `4.1e-5` | 1.000 / 0.006, `2.4e-5` | 複製 |
| `(3) G>A` | 3.3 | 0.416 / 0.844, `0.025` | 0.257 / 0.993, `0.0061` | 複製 |
| `(3) G>A` | 3.4 | 0.045 / 0.874, `4.1e-5` | 0.256 / 0.858, `0.060` | 方向複製，DORADO 未過 FDR |
| `(3) G>A` | 3.5 | 0.010 / 0.865, `0.0017` | 0.380 / 0.835, `0.084` | 方向複製，DORADO 未過 FDR |
| `(4) non-ref` | corrected 4.1 | 0.177 / 0.975, `2.2e-6` | 0.516 / 0.991, `8.7e-5` | 複製 |
| `(5) G>C` | 5.1 | 0.116 / 0.749, `0.0037` | 0.005 / 0.649, `0.0020` | 複製 |
| `(5) G>C` | 5.2 | 0.030 / 0.698, `0.0012` | 0.010 / 0.506, `0.0020` | 複製 |
| `(6) C>G` | 5.1 | 0.939 / 0.189, `0.0015` | 1.000 / 0.096, `0.0068` | 複製 |
| `(6) C>G` | 5.2 | 0.885 / 0.031, `0.0055` | 1.000 / 0.010, `0.0054` | 複製 |

這證明：

> **甲基不是只有視覺印象；多個局部 allele-CpG association 在 HKU 與 DORADO 方向一致並通過 FDR。**

但這只證明 association / lineage coherence，不證明甲基造成突變，也不證明甲基可獨立定義 clone。

## 五、使用者手標甲基表需要修正

### 明確修正

1. `(2.2)`：
   - 使用者表：alpha 低、beta 高。
   - 實測：alpha/reference 高、`(2)C` beta-like ALT 低。

2. `(5.1)(5.2)`：
   - 使用者表：`(5)C` 群高、beta-like `(5)G` 群低。
   - 實測：`(5)C` ALT 群低；`(6)G` beta-like 群高。

3. `(4.1)`：
   - 使用者文字座標錯誤。
   - 改用圖上的 `18,096,341` 後，`(4)non-ref` vs REF 的甲基差異可在兩 basecaller 複製。

因此後續圖與論文表格應直接由 audit JSON 產生，不再手填 H/L。

## 六、matched-normal 顯示的甲基 confound

HKU normal HP1 vs HP2 在下列 CpG 已有顯著 ASM：

| CpG | Normal HP1 mean | Normal HP2 mean | MW FDR |
|---|---:|---:|---:|
| 3.3 | 0.814 | 0.046 | 0.0070 |
| 3.4 | 0.033 | 0.850 | 0.0070 |
| 3.5 | 0.583 | 0.012 | 0.0496 |
| corrected 4.1 | 0.887 | 0.078 | 0.0070 |

這表示：

- 3.3-4.1 區段的 tumor lineage methylation 差異，部分可能受既存 haplotype ASM 影響。
- 不能把所有 tumor 甲基差異直接稱為 subclone 形成後的新 epimutation。

但 normal HP 在 `2.1/2.2`、`3.1/3.2`、`5.1/5.2` 沒有相同顯著差異，而 tumor allele association 在這些位置跨 basecaller 複製。這些位置是較乾淨的 **tumor-associated lineage methylation candidates**。

最嚴謹表述：

> 此區域同時包含 normal haplotype ASM 與 tumor-associated allele-methylation structure；ISM 的價值是將兩者拆開，而不是把所有甲基差異都歸因於 subclone。

### Tumor/normal differential methylation 也可跨資料重現

以所有 reads 直接比較 tumor 與 matched normal，HKU 與 DORADO 在 `2.1/2.2`、`3.1/3.2`、`5.1/5.2` 都有同方向且通過 FDR 的差異；例如：

| CpG | HKU tumor / normal, MW q | DORADO tumor / normal, MW q |
|---|---|---|
| 2.1 | 0.525 / 0.030, `7.7e-4` | 0.497 / 0.016, `3.5e-4` |
| 3.1 | 0.561 / 0.045, `1.2e-3` | 0.476 / 0.017, `1.0e-3` |
| 5.1 | 0.560 / 0.033, `7.7e-4` | 0.429 / 0.006, `6.1e-3` |

這支持使用者所述的「tumor/normal 有差異甲基位點」輸出目標，但 all-read tumor/normal 差異仍混合了 tumor composition、LOH 與 haplotype ASM；必須再用 allele-conditioned 與 normal-HP-conditioned 分析拆解。

## 七、為何 pos4 的 G/T/DEL 不能算三個 subclone

`chr2:18,096,269` 後方緊接約 20 bp poly-T：

```text
...AGTCT[C]TTTTTTTTTTTTTTTTTTTTGAGAC...
```

SEQC2 truth 是 `C>G MedConf`，但 ONT 在此位點產生 G/T/DEL 混合 call：

| Dataset | G | T | DEL | REF C |
|---|---:|---:|---:|---:|
| HKU BQ20 | 4 | 13 | 10 | 21 |
| DORADO BQ20 | 6 | 15 | 7 | 29 |

反對「三 subclone」的證據：

1. 兩套資料都出現相似 G/T/DEL 混合，符合平台/比對 homopolymer jitter。
2. HKU 中 G/T/DEL 與 `(6)G` 都落在同一 pos4-altered state。
3. G/T/DEL 在其他 CpG 大致共享 beta-like methylation；DORADO corrected 4.1 mean 也接近。
4. HKU 的 G call 只有 4 條且全在 reverse strand，不能視為乾淨獨立群。

因此 `[推論 6] pos4 多突變造成群三、四、五` 應 **撤回**。

## 八、推論 1-6 最終判定

| 推論 | 判定 | 修正版 |
|---|---|---|
| 1. 發生 LOH | **確認** | SEQC2 LOH BED + HKU/DORADO HP imbalance 一致 |
| 2. 發生 subclone | **強支持、有界** | 可稱 regional operational subclonal states；非完整 biological clone truth |
| 3. HP2 沒抽到／突變到不見 | **否定** | HP2 是主體；all-REF 短 linkage reads 存在；完整 read 缺失是 coverage |
| 4. `(1)(2)(3)(6)` 先突變 | **不成立／需重寫** | `(3)` 與 `(1)(2)(6)` 互斥，不能是同一 trunk；只可稱不同早期 branch-defining candidates |
| 5. `(5)` 將群一、群二分開 | **支持** | `(5)C` 對 `(3)A` 呈 perfect nesting；目前 direct support HKU 4、DORADO 2，仍屬 provisional |
| 6. `(4)` 多突變造成三群 | **否定** | 合併為 pos4-altered beta-like state；G/T/DEL 為 homopolymer-uncertain |

## 九、可防守的候選樹

```text
all-REF / unsampled common ancestor
  |
  +-- alpha: (3) G>A
  |     |
  |     +-- alpha-1: + (5) G>C
  |
  +-- beta-like program
        |
        +-- left block:  (1) C>G + (2) G>C
        |
        +-- right block: (4) C>G truth, observed as G/T/DEL jitter + (6) C>G
```

必要註記：

- alpha 與 beta-like 分叉有直接互斥支持。
- `(3) → (5)` 有 direct nesting 支持。
- beta left/right block 之間是以互斥 + methylation coherence 連接，缺直接 genetic bridge，應畫虛線。
- `(4)` 與 `(6)` 共現，不能排序先後。
- root 是 parsimony / contamination-aware 假說，不是被完整 spanning read 觀察到的 clone。

## 十、與 ISM 研究目標的關係

此案例確實支持 ISM 最有價值的研究定位：

1. 同一長讀上聯合觀測 somatic allele、haplotag 與 native 5mC。
2. 找出跨 basecaller 重現的 allele-CpG association。
3. 在 LongPhase-S fine tag 無法拆開的情況下，補充 regional molecular-state characterization。
4. 使用 matched normal 揭露哪些甲基差異可能是既存 ASM confound。
5. 用多點 read linkage 提出有不確定性標記的 regional subclone hypothesis。

但不能宣稱：

- 甲基造成突變。
- 甲基獨立重建完整 subclone tree。
- 五個 subclone 已被證實。
- 完整演化順序已確認。
- DORADO 是獨立病人／獨立 biological replicate；它是同細胞株的技術/資料重現。

## 十一、外部方法學對齊

外部文獻也支持上述 claim 邊界：

- [Tarabichi et al., Nature Methods 2021](https://doi.org/10.1038/s41592-020-01013-2)：正式 subclonal reconstruction 需處理 genetic backbone、CN/purity、clustering、phylogeny 與不確定性。
- [Liu/Goretsky et al. 2025 preprint](https://doi.org/10.1101/2025.08.28.672865)：最接近的 long-read subclone multi-omics 案例以 SNV 建樹，再把 methylation 掛到 lineage；支持「genetic backbone + methylation characterization」。
- [LongPhase-S preprint](https://doi.org/10.1101/2025.11.20.689492)：somatic haplotagging 是本案例的 genetic backbone；ISM 的差異化是加入有界、normal-aware 的 read-level methylation interpretation。

因此此案例適合論文用語是：

> **A regional, LOH-constrained, somatic-haplotag-conditioned subclonal structure with cross-basecaller methylation coherence.**

不適合：

> **Confirmed five-subclone evolutionary reconstruction.**

## 可重現命令

```bash
python3 \
  docs/experiments/in_progress/2026/06/20260615_chr2_18M_subclone_verification_assets/scripts/independent_subclone_audit.py \
  --json-out docs/experiments/in_progress/2026/06/20260615_chr2_18M_subclone_verification_assets/data/independent_audit.json \
  --md-out docs/experiments/in_progress/2026/06/20260615_chr2_18M_subclone_verification_assets/data/independent_audit.md
```
