<!--
建立時間: 2026-03-16 17:05
目標: 依原始 IGV template 的 panel 內容、格式、樣本與順序，整理 5 TP + 4 FP 的忠實版 case 導讀
處理範圍:
  - /big8_disk/liaoyoyo2001/IGV_session/template.xml
  - full-template 9 loci PNG 與 session XML
  - 整合 20260315 case study 與 20260316 phase4 機制分析
關聯檔案:
  - /big8_disk/liaoyoyo2001/IGV_session/template.xml
  - /big8_disk/liaoyoyo2001/InterSubMod/scripts/analysis/igv_snapshot_from_template.sh
  - /big8_disk/liaoyoyo2001/InterSubMod/docs/reports/validated/2026/03/20260315_subclone_phase4_casestudy_01.md
  - /big8_disk/liaoyoyo2001/InterSubMod/docs/reports/validated/2026/03/20260316_Phase4_甲基亞克隆視覺觀察與FP機制分析報告_01.md
  - /big8_disk/liaoyoyo2001/InterSubMod/docs/references/manual/20260316_IGV自動截圖操作手冊_01.md
-->

# IGV full template 忠實版案例導讀

## 開頭重點結論

1. 目前這批 9 張 full-template PNG 可視為 **忠實保留** `/big8_disk/liaoyoyo2001/IGV_session/template.xml` 的 panel 與樣本順序版本。
2. 這條路徑最適合做「跟你平常手動打開 IGV 時一樣」的 case review，尤其適合寫逐 panel 的觀察說明。
3. 真正最有效率的判讀順序不是把 11 個 panel 全部平均看，而是：
   - `DataPanel`
   - 主 `tumor BAM`
   - 主 `normal BAM`
   - `FeaturePanel`
   - 視需要再看 `LongPhase-S / Dorado / expanded` 面板
4. 5 個 TP 與 4 個 FP 的主要差異，不在「有沒有甲基差」，而在：
   - ALT reads 是否在 tumor 中形成可見群
   - 該群是由 allele 驅動，還是被 HP / 鄰近 INDEL / MNP 機制帶著走

## 背景與驗證

### 本文件用的是哪一批圖

- full template PNG：
  - `/big8_disk/liaoyoyo2001/InterSubMod/docs/reports/validated/2026/03/assets/20260316_subclone_phase4_igv_sessions/snapshots/`
- per-locus session XML：
  - `/big8_disk/liaoyoyo2001/InterSubMod/docs/reports/validated/2026/03/assets/20260316_subclone_phase4_igv_sessions/sessions/`
- 9 個位點都已存在：
  - `5 TP`
  - `4 FP`

### 為什麼說這是忠實版

根據 [igv_snapshot_from_template.sh](/big8_disk/liaoyoyo2001/InterSubMod/scripts/analysis/igv_snapshot_from_template.sh#L219)，產生 per-locus session 時只改：

1. `Session locus`
2. 所有 alignment track 的 `groupByPos`

因此：

1. panel 順序保留
2. track 順序保留
3. 樣本順序保留
4. coverage / alignment / feature 的上下層結構保留

這與 [IGV 自動截圖操作手冊](/big8_disk/liaoyoyo2001/InterSubMod/docs/references/manual/20260316_IGV自動截圖操作手冊_01.md#L110) 的結論一致。

### 整合哪些前一份文件

本文件主要整合兩份：

1. [20260315_subclone_phase4_casestudy_01.md](/big8_disk/liaoyoyo2001/InterSubMod/docs/reports/validated/2026/03/20260315_subclone_phase4_casestudy_01.md)
   - 原始 case 數值、heatmap、distance、methylartist
2. [20260316_Phase4_甲基亞克隆視覺觀察與FP機制分析報告_01.md](/big8_disk/liaoyoyo2001/InterSubMod/docs/reports/validated/2026/03/20260316_Phase4_甲基亞克隆視覺觀察與FP機制分析報告_01.md)
   - 後續的機制驗證、TP/FP 推論與修正口徑

## 原始 template 的 panel 地圖

### 建議主閱讀順序

| 順序 | 面板                                | 主要用途                                      | 何時最重要               |
| ---- | ----------------------------------- | --------------------------------------------- | ------------------------ |
| 1    | `DataPanel`                         | 先看這個位點在 keep/remove/TP/FP 的分類上下文 | 一開始定位 case          |
| 2    | 主 `tumor BAM`                      | 看 ALT/REF/HP 是否在 read-level 成群          | 最核心                   |
| 3    | 主 `normal BAM`                     | 看 normal 背景是否已有類似 pattern            | 判斷 tumor-specific      |
| 4    | `FeaturePanel`                      | 看 caller、SEQC2、鄰近 INDEL、annotation      | 做機制判讀               |
| 5    | `LongPhase-S tumor tagged BAM`      | 補看 somatic haplotag 版本的分群              | 想對照 haplotag 結果時   |
| 6    | `Dorado tumor / normal`             | 跨平台方向比較                                | 懷疑 5kHz 特異時         |
| 7    | `expanded tumor / normal alignment` | 看更細 read 細節                              | 懷疑 MNP/錯配/複雜事件時 |

### 每個面板看到什麼代表什麼

| 面板                 | 要看什麼現象                                                | 代表什麼意義                                  |
| -------------------- | ----------------------------------------------------------- | --------------------------------------------- |
| `DataPanel`          | 變異落在 kept / removed / tp / fp 哪個軌道                  | 先知道 pipeline 對此位點的原始分類            |
| 主 `tumor BAM`       | ALT read 是否集中、HP 顏色是否與 ALT 共線、甲基色塊是否分層 | 決定這是不是 read-level 真訊號                |
| 主 `normal BAM`      | 同一區域在 normal 是否也出現相同分層                        | 區分 tumor-specific 與 germline 背景          |
| `LongPhase-S tumor`  | 同一位點在另一套 somatic haplotag track 是否也有一致分層    | 驗證 HP/ALT 對應是否穩定                      |
| `Dorado` 對照        | Dorado 是否也看到相同方向的 read 結構                       | 判斷 5kHz 專屬或跨平台可重現                  |
| `expanded alignment` | 單 read 錯配、相鄰位點異常、局部複雜性                      | 判斷 INDEL adjacency、MNP、alignment artifact |
| `FeaturePanel`       | HC region、SEQC2 sSNV/sINDEL、caller output                 | 判斷 benchmark 是否支持或有 annotation gap    |

## 9 個案例總覽

| Case  | 類型 | 位點                 | full-template PNG                                                                       | 核心結論                      |
| ----- | ---- | -------------------- | --------------------------------------------------------------------------------------- | ----------------------------- |
| TP_01 | TP   | `chr6:145444893 G>A` | ![PNG](assets/20260316_subclone_phase4_igv_sessions/snapshots/TP_01_chr6_145444893.png) | allele-only，最乾淨的 ALT 群  |
| TP_02 | TP   | `chr4:70548355 G>A`  | ![PNG](assets/20260316_subclone_phase4_igv_sessions/snapshots/TP_02_chr4_70548355.png)  | 高 HP0 背景下仍有 allele 分離 |
| TP_03 | TP   | `chr5:153209947 C>A` | ![PNG](assets/20260316_subclone_phase4_igv_sessions/snapshots/TP_03_chr5_153209947.png) | 低甲基背景中的微弱差異        |
| TP_04 | TP   | `chr16:35118902 G>A` | ![PNG](assets/20260316_subclone_phase4_igv_sessions/snapshots/TP_04_chr16_35118902.png) | 最清楚的 ALT 低甲基雙峰       |
| TP_05 | TP   | `chr7:109185781 G>T` | ![PNG](assets/20260316_subclone_phase4_igv_sessions/snapshots/TP_05_chr7_109185781.png) | HP+Allele 共線                |
| FP_A1 | FP   | `chr8:93565727 C>T`  | ![PNG](assets/20260316_subclone_phase4_igv_sessions/snapshots/FP_A1_chr8_93565727.png)  | ALT 太少，不足成群            |
| FP_A2 | FP   | `chr9:137953060 T>C` | ![PNG](assets/20260316_subclone_phase4_igv_sessions/snapshots/FP_A2_chr9_137953060.png) | 高 CpG 假顯著                 |
| FP_B1 | FP   | `chr7:52087777 A>T`  | ![PNG](assets/20260316_subclone_phase4_igv_sessions/snapshots/FP_B1_chr7_52087777.png)  | HP-driven，且鄰近 SEQC2 INDEL |
| FP_B2 | FP   | `chr9:75383880 T>A`  | ![PNG](assets/20260316_subclone_phase4_igv_sessions/snapshots/FP_B2_chr9_75383880.png)  | MNP 機制，HP-driven           |

## TP 案例導讀

### TP_01 `chr6:145444893 G>A`

![TP_01 full template](/big8_disk/liaoyoyo2001/InterSubMod/docs/reports/validated/2026/03/assets/20260316_subclone_phase4_igv_sessions/snapshots/TP_01_chr6_145444893.png)

數值摘要：

| 指標          | 值                     |
| ------------- | ---------------------- |
| CramersV      | 0.722                  |
| ALT reads     | 19/58                  |
| ALT-REF delta | +0.039                 |
| 核心結論      | allele-only，HP 不分離 |

逐 panel 觀察：

| 先看哪個 panel  | 要看什麼現象                                         | 這個現象代表什麼意義                |
| --------------- | ---------------------------------------------------- | ----------------------------------- |
| `DataPanel`     | 該位點落在 TP / kept 軌道而非 fp_cluster             | 先建立它是正例背景                  |
| 主 `tumor BAM`  | ALT read 在同一層集中，且與一批 REF 有可見間隔       | 這是 allele-driven 群，不是雜訊散布 |
| 主 `normal BAM` | normal 沒有對應的 ALT 群                             | 支持 tumor-specific                 |
| `FeaturePanel`  | caller 與 HC region 都支持，但 HP 軌沒有額外分離訊號 | 表示分群主因是 allele，不是 HP      |

對照來源：

- ![原始 case study](/big8_disk/liaoyoyo2001/InterSubMod/docs/reports/validated/2026/03/20260315_subclone_phase4_casestudy_01.md#L51)
- [後續結論](/big8_disk/liaoyoyo2001/InterSubMod/docs/reports/validated/2026/03/20260316_Phase4_甲基亞克隆視覺觀察與FP機制分析報告_01.md#L56)

### TP_02 `chr4:70548355 G>A`

![TP_02 full template](/big8_disk/liaoyoyo2001/InterSubMod/docs/reports/validated/2026/03/assets/20260316_subclone_phase4_igv_sessions/snapshots/TP_02_chr4_70548355.png)

數值摘要：

| 指標      | 值                            |
| --------- | ----------------------------- |
| CramersV  | 0.502                         |
| ALT reads | 11/48                         |
| HP0 比例  | 32/48 = 67%                   |
| 核心結論  | 高 HP0 背景下仍有 allele 分離 |

逐 panel 觀察：

| 先看哪個 panel      | 要看什麼現象                                     | 這個現象代表什麼意義                           |
| ------------------- | ------------------------------------------------ | ---------------------------------------------- |
| 主 `tumor BAM`      | 雖然很多 read 沒有清楚 HP，但 ALT 仍偏在一個主群 | 甲基差異不完全依賴 phasing 才能看到            |
| 主 `normal BAM`     | normal 主要提供背景，不會出現同型 ALT 分群       | 不是 normal 先天分兩群造成                     |
| `LongPhase-S tumor` | HP0 很多，但非所有 HP0 都混成單一亂群            | HP 不足不等於整個位點不可讀                    |
| `FeaturePanel`      | HC region 內、caller 支持存在，但 phase 資訊不足 | 這類位點要寫成 allele-supported、phase-limited |

對照來源：

- [原始 case study](/big8_disk/liaoyoyo2001/InterSubMod/docs/reports/validated/2026/03/20260315_subclone_phase4_casestudy_01.md#L94)
- [後續結論](/big8_disk/liaoyoyo2001/InterSubMod/docs/reports/validated/2026/03/20260316_Phase4_甲基亞克隆視覺觀察與FP機制分析報告_01.md#L67)

### TP_03 `chr5:153209947 C>A`

![TP_03 full template](/big8_disk/liaoyoyo2001/InterSubMod/docs/reports/validated/2026/03/assets/20260316_subclone_phase4_igv_sessions/snapshots/TP_03_chr5_153209947.png)

數值摘要：

| 指標          | 值                         |
| ------------- | -------------------------- |
| CramersV      | 0.443                      |
| ALT reads     | 17/57                      |
| ALT-REF delta | +0.035                     |
| 核心結論      | 全域低甲基背景中的微弱差異 |

逐 panel 觀察：

| 先看哪個 panel       | 要看什麼現象                                           | 這個現象代表什麼意義                  |
| -------------------- | ------------------------------------------------------ | ------------------------------------- |
| 主 `tumor BAM`       | read 顏色整體很接近，沒有像 chr16 那樣的強雙峰         | 這是弱效應、不是肉眼超強分離型        |
| 主 `normal BAM`      | normal 也偏低甲基背景                                  | 整體區域本來就是低甲基場景            |
| `FeaturePanel`       | caller 支持存在，但僅靠 IGV 視覺不容易強判             | 這種 case 主要仍仰賴統計而非肉眼      |
| `expanded alignment` | 若看細節，重點不是單條 read 很特別，而是群體小幅 shift | 支持「微弱但一致」而非單讀段 artifact |

對照來源：

- [原始 case study](/big8_disk/liaoyoyo2001/InterSubMod/docs/reports/validated/2026/03/20260315_subclone_phase4_casestudy_01.md#L131)
- [後續結論](/big8_disk/liaoyoyo2001/InterSubMod/docs/reports/validated/2026/03/20260316_Phase4_甲基亞克隆視覺觀察與FP機制分析報告_01.md#L78)

### TP_04 `chr16:35118902 G>A`

![TP_04 full template](/big8_disk/liaoyoyo2001/InterSubMod/docs/reports/validated/2026/03/assets/20260316_subclone_phase4_igv_sessions/snapshots/TP_04_chr16_35118902.png)

數值摘要：

| 指標          | 值                      |
| ------------- | ----------------------- |
| CramersV      | 0.442                   |
| ALT reads     | 25/91                   |
| ALT-REF delta | −0.198                  |
| 核心結論      | 最清楚的 ALT 低甲基雙峰 |

逐 panel 觀察：

| 先看哪個 panel      | 要看什麼現象                                         | 這個現象代表什麼意義                                   |
| ------------------- | ---------------------------------------------------- | ------------------------------------------------------ |
| 主 `tumor BAM`      | 可以看到最明顯的高甲基群 vs 低甲基群對照             | 這是最像真正 subclone methylation reprogramming 的案例 |
| 主 `normal BAM`     | normal 偏向高甲基背景，沒有對應的低甲基 ALT 群       | 支持 tumor-specific hypomethylated clone               |
| `LongPhase-S tumor` | HP=3 / HP=1-1 的 ALT 都偏到低甲基側                  | ALT 群即使 phase 複雜，仍有一致方向                    |
| `FeaturePanel`      | caller 與 allele context 支持，非單純鄰近 INDEL 問題 | 可作為 prototype TP                                    |

對照來源：

- [原始 case study](/big8_disk/liaoyoyo2001/InterSubMod/docs/reports/validated/2026/03/20260315_subclone_phase4_casestudy_01.md#L169)
- [後續結論](/big8_disk/liaoyoyo2001/InterSubMod/docs/reports/validated/2026/03/20260316_Phase4_甲基亞克隆視覺觀察與FP機制分析報告_01.md#L90)

### TP_05 `chr7:109185781 G>T`

![TP_05 full template](/big8_disk/liaoyoyo2001/InterSubMod/docs/reports/validated/2026/03/assets/20260316_subclone_phase4_igv_sessions/snapshots/TP_05_chr7_109185781.png)

數值摘要：

| 指標          | 值             |
| ------------- | -------------- |
| CramersV      | 0.429          |
| ALT reads     | 28/127         |
| ALT-REF delta | −0.025         |
| 核心結論      | HP+Allele 共線 |

逐 panel 觀察：

| 先看哪個 panel  | 要看什麼現象                                     | 這個現象代表什麼意義              |
| --------------- | ------------------------------------------------ | --------------------------------- |
| 主 `tumor BAM`  | ALT 幾乎跟著某一側 HP 背景走，而不是獨立切成新層 | 這類位點要先懷疑 HP confounding   |
| 主 `normal BAM` | normal 也可能已有 HP1/HP2 背景差異               | 支持 germline HP 對甲基模式有影響 |
| `Dorado` 對照   | 若方向一致，代表不是單純 5kHz 顯示偏差           | 用來判斷是否跨平台共線            |
| `FeaturePanel`  | 雖然 allele 與 HP 都顯著，但 delta 很小          | 統計顯著不等於 allele 是主因      |

對照來源：

- [原始 case study](/big8_disk/liaoyoyo2001/InterSubMod/docs/reports/validated/2026/03/20260315_subclone_phase4_casestudy_01.md#L208)
- [後續結論](/big8_disk/liaoyoyo2001/InterSubMod/docs/reports/validated/2026/03/20260316_Phase4_甲基亞克隆視覺觀察與FP機制分析報告_01.md#L104)

## FP 案例導讀

### FP_A1 `chr8:93565727 C>T`

![FP_A1 full template](/big8_disk/liaoyoyo2001/InterSubMod/docs/reports/validated/2026/03/assets/20260316_subclone_phase4_igv_sessions/snapshots/FP_A1_chr8_93565727.png)

數值摘要：

| 指標               | 值                     |
| ------------------ | ---------------------- |
| CramersV           | 0.000                  |
| ALT reads          | 2/68                   |
| PairwiseMedianDist | 0.359                  |
| 核心結論           | VAF 太低，ALT 無法成群 |

逐 panel 觀察：

| 先看哪個 panel       | 要看什麼現象                         | 這個現象代表什麼意義                 |
| -------------------- | ------------------------------------ | ------------------------------------ |
| 主 `tumor BAM`       | ALT 只有零星幾條，且不在同一群       | 低 VAF 到不足形成穩定 visual cluster |
| 主 `normal BAM`      | 沒有可支持 tumor-specific 結構的對照 | 更像雜訊背景                         |
| `expanded alignment` | 高異質、局部亂                       | repeat / multi-mapping 風險高        |
| `FeaturePanel`       | 若缺乏 truth 支持且 ALT 極少         | 這類通常直接降權為 FP                |

對照來源：

- [FP 機制分析](/big8_disk/liaoyoyo2001/InterSubMod/docs/reports/validated/2026/03/20260316_Phase4_甲基亞克隆視覺觀察與FP機制分析報告_01.md#L158)

### FP_A2 `chr9:137953060 T>C`

![FP_A2 full template](/big8_disk/liaoyoyo2001/InterSubMod/docs/reports/validated/2026/03/assets/20260316_subclone_phase4_igv_sessions/snapshots/FP_A2_chr9_137953060.png)

數值摘要：

| 指標          | 值            |
| ------------- | ------------- |
| CramersV      | 0.000         |
| ALT reads     | 5/58          |
| ALT-REF delta | −0.003        |
| 核心結論      | 高 CpG 假顯著 |

逐 panel 觀察：

| 先看哪個 panel  | 要看什麼現象                                | 這個現象代表什麼意義                       |
| --------------- | ------------------------------------------- | ------------------------------------------ |
| 主 `tumor BAM`  | ALT 不集中成群，視覺上與 REF 幾乎無差       | read-level 觀察不支持真實 allele 分層      |
| 主 `normal BAM` | normal 只提供背景，無明顯對應異常           | 不像 tumor-only 的新群                     |
| `FeaturePanel`  | 雖然某些統計 p-value 可顯著，但效應量接近 0 | 這種位點要優先相信效應量而不是高維 p-value |

對照來源：

- [FP 機制分析](/big8_disk/liaoyoyo2001/InterSubMod/docs/reports/validated/2026/03/20260316_Phase4_甲基亞克隆視覺觀察與FP機制分析報告_01.md#L172)

### FP_B1 `chr7:52087777 A>T`

![FP_B1 full template](/big8_disk/liaoyoyo2001/InterSubMod/docs/reports/validated/2026/03/assets/20260316_subclone_phase4_igv_sessions/snapshots/FP_B1_chr7_52087777.png)

數值摘要：

| 指標          | 值                          |
| ------------- | --------------------------- |
| CramersV      | 1.000                       |
| ALT reads     | 9/63                        |
| ALT-REF delta | +0.058                      |
| 核心結論      | HP-driven + INDEL adjacency |

逐 panel 觀察：

| 先看哪個 panel       | 要看什麼現象                                       | 這個現象代表什麼意義                             |
| -------------------- | -------------------------------------------------- | ------------------------------------------------ |
| 主 `tumor BAM`       | 分群非常漂亮，但分的是 HP1 vs HP2，不是 ALT vs REF | 高 CramersV 也可能只是 germline HP 差            |
| 主 `normal BAM`      | normal 若也保有 HP 分層背景                        | 更支持是 germline imprinting / HP effect         |
| `FeaturePanel`       | 鄰近有 SEQC2 sINDEL truth                          | 這個 SNV 標成 FP 可能是 benchmark annotation gap |
| `expanded alignment` | 看 52087776-52087777 鄰近錯配 / shift              | 支持 INDEL-adjacent artifact 機制                |

對照來源：

- [FP 機制分析](/big8_disk/liaoyoyo2001/InterSubMod/docs/reports/validated/2026/03/20260316_Phase4_甲基亞克隆視覺觀察與FP機制分析報告_01.md#L185)

### FP_B2 `chr9:75383880 T>A`

![FP_B2 full template](/big8_disk/liaoyoyo2001/InterSubMod/docs/reports/validated/2026/03/assets/20260316_subclone_phase4_igv_sessions/snapshots/FP_B2_chr9_75383880.png)

數值摘要：

| 指標          | 值                                   |
| ------------- | ------------------------------------ |
| CramersV      | 0.736                                |
| ALT reads     | 14/62                                |
| ALT-REF delta | +0.005                               |
| 核心結論      | MNP 機制，且甲基分群主要是 HP-driven |

逐 panel 觀察：

| 先看哪個 panel       | 要看什麼現象                        | 這個現象代表什麼意義                    |
| -------------------- | ----------------------------------- | --------------------------------------- |
| 主 `tumor BAM`       | ALT read 在 target 前一格也一致改變 | 先懷疑不是單一獨立 SNV                  |
| `expanded alignment` | 75383879 與 75383880 兩格一起異常   | 支持 MNP 而非 indel                     |
| 主 `normal BAM`      | normal 沒有相同 A 型態              | 表示這是 tumor-specific haplotype 事件  |
| `FeaturePanel`       | truth set 沒有對應 MNP / SNV 標註   | 這類會在 SNV-only benchmark 中被標成 FP |

對照來源：

- [FP 機制分析](/big8_disk/liaoyoyo2001/InterSubMod/docs/reports/validated/2026/03/20260316_Phase4_甲基亞克隆視覺觀察與FP機制分析報告_01.md#L227)

## 如何把 full template PNG 寫成固定格式的觀察文字

建議每個 case 都照這個順序寫：

1. `DataPanel`：
   - 這個位點在 pipeline 的原始分類是什麼
2. 主 `tumor BAM`：
   - ALT / REF / HP / 甲基 read-level 現象是什麼
3. 主 `normal BAM`：
   - normal 是否也有類似現象
4. `FeaturePanel`：
   - truth / caller / annotation 是否支持
5. 結論：
   - 是 allele-driven TP
   - HP-driven confound
   - INDEL adjacency
   - MNP
   - 高 HP0 / 低 VAF / 高維假顯著

## 最後結論

1. full template 忠實版 PNG 確實適合做「逐 panel、逐樣本」的教學式觀察。
2. 它最大的價值不是比 slim template 更穩，而是 **更貼近你平常打開 template.xml 時的閱讀習慣**。
3. 若目標是大量自動截圖，優先用 slim template。
4. 若目標是寫 case 導讀、做固定順序的觀察說明與驗證，full template 更合適。

## 參考文件索引

1. [原始 template](/big8_disk/liaoyoyo2001/IGV_session/template.xml)
2. [IGV 自動截圖操作手冊](/big8_disk/liaoyoyo2001/InterSubMod/docs/references/manual/20260316_IGV自動截圖操作手冊_01.md)
3. [Phase4 case study 原始文件](/big8_disk/liaoyoyo2001/InterSubMod/docs/reports/validated/2026/03/20260315_subclone_phase4_casestudy_01.md)
4. [Phase4 TP/FP 機制分析](/big8_disk/liaoyoyo2001/InterSubMod/docs/reports/validated/2026/03/20260316_Phase4_甲基亞克隆視覺觀察與FP機制分析報告_01.md)
5. [full-template 截圖資產目錄](/big8_disk/liaoyoyo2001/InterSubMod/docs/reports/validated/2026/03/assets/20260316_subclone_phase4_igv_sessions)
