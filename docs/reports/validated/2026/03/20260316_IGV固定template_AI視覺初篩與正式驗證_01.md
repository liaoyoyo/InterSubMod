<!--
建立時間: 2026-03-16 19:10
目標: 以固定 IGV template 的標準化 PNG，實測 AI 視覺初篩是否能辨識 ALT 成群、HP/allele 驅動與相鄰異常，並用 BAM/VCF/統計值正式驗證
處理範圍:
  - full-template 9 loci PNG
  - HCC1395 5kHz paired tumor/normal BAM
  - SEQC2 sSNV / sINDEL truth
  - analysis/feature_matrix_20260315.tsv
關聯檔案:
  - /big8_disk/liaoyoyo2001/IGV_session/template.xml
  - assets/20260316_subclone_phase4_igv_sessions/snapshots/
  - assets/20260316_igv_case_validation_01.tsv
  - assets/20260316_fp_cases_for_mnp_screen_01.vcf
  - assets/20260316_fp_mnp_screen_result_01.tsv
  - ../../../../scripts/analysis/validate_igv_case_loci.py
  - ../../../../scripts/analysis/screen_mnp_adjacent_fp.py
  - 20260315_subclone_phase4_casestudy_01.md
  - 20260316_Phase4_甲基亞克隆視覺觀察與FP機制分析報告_01.md
-->

# IGV 固定 template AI 視覺初篩與正式驗證

## 一句結論

固定 template 的標準化 IGV PNG，**很適合做第一輪 AI 視覺初篩**；其中 `TP_01`、`TP_04`、`FP_A1`、`FP_B1` 的方向幾乎可直接從圖上先抓到。  
但 `TP_03`、`FP_A2` 這類弱效應或高維假顯著案例，**只靠圖片不夠**；而 `FP_B2` 更顯示出一個重要修正：先前看起來像 `MNP` 的現象，在主 tumor BAM 的嚴格 read-level 驗證下，較符合 **`target SNV + -1 deletion-like gap`**，不是乾淨的相鄰雙 SNV。

## 1. 本次研究問題與驗證設計

### 1.1 研究問題

1. 固定 template 批次產出的標準化 PNG，是否足夠讓 AI 先做穩定的視覺初篩？
2. AI 能否先抓出以下 5 類現象：
   - ALT 是否成群
   - normal 是否也有相同結構
   - 是否有相鄰異常柱
   - 是否像 HP-driven 而不是 allele-driven
   - 是否有明顯 coverage / alignment anomaly
3. 哪些現象必須回到 BAM / VCF / 統計值才能正式定論？

### 1.2 本次假設

1. `強視覺分離` 的 case，AI 可以先抓方向，再由 BAM/VCF 補證。
2. `弱效應` 或 `高維假顯著` 的 case，AI 會不穩，必須依賴統計量。
3. `相鄰異常柱` 只靠 IGV 畫面容易誤判，必須用 read-level 對齊資訊確認。

### 1.3 成功條件

1. 9 個 loci 都有固定格式的 AI 初篩結論。
2. 9 個 loci 都有對應的 tumor/normal pileup 與 truth 查核。
3. 能清楚區分：
   - `AI 可以安全先判的現象`
   - `必須回 BAM / truth / feature matrix 才能定論的現象`

## 2. 資料與方法

### 2.1 固定 PNG 與對照來源

- full-template PNG：[/big8_disk/liaoyoyo2001/InterSubMod/docs/reports/validated/2026/03/assets/20260316_subclone_phase4_igv_sessions/snapshots](/big8_disk/liaoyoyo2001/InterSubMod/docs/reports/validated/2026/03/assets/20260316_subclone_phase4_igv_sessions/snapshots)
- 原始 template：[/big8_disk/liaoyoyo2001/IGV_session/template.xml](/big8_disk/liaoyoyo2001/IGV_session/template.xml)
- AI 初篩對照的既有 case 文件：
  - [20260315_subclone_phase4_casestudy_01.md](/big8_disk/liaoyoyo2001/InterSubMod/docs/reports/validated/2026/03/20260315_subclone_phase4_casestudy_01.md)
  - [20260316_Phase4_甲基亞克隆視覺觀察與FP機制分析報告_01.md](/big8_disk/liaoyoyo2001/InterSubMod/docs/reports/validated/2026/03/20260316_Phase4_甲基亞克隆視覺觀察與FP機制分析報告_01.md)

### 2.2 正式驗證資料

- tumor BAM：`/big8_disk/liaoyoyo2001/data/bam/HCC1395_Tmode_tagged_ClairS_pileup_v040_woFilter.bam`
- normal BAM：`/big8_disk/liaoyoyo2001/data/bam/HCC1395BL_ONT_5khz_simplex_5mCG_5hmCG_tagged.bam`
- reference：`/big8_disk/ref/GRCh38_no_alt_analysis_set.fasta`
- truth SNV：`/big8_disk/data/HCC1395/SEQC2/high-confidence_sSNV_in_HC_regions_v1.2.1.vcf.gz`
- truth INDEL：`/big8_disk/data/HCC1395/SEQC2/high-confidence_sINDEL_in_HC_regions_v1.2.1.vcf.gz`
- feature matrix：`analysis/feature_matrix_20260315.tsv`

### 2.3 本輪新增的正式驗證輸出

- case 驗證表：[20260316_igv_case_validation_01.tsv](/big8_disk/liaoyoyo2001/InterSubMod/docs/reports/validated/2026/03/assets/20260316_igv_case_validation_01.tsv)
- 4 個 FP 的 MNP 嚴格檢查輸入：[20260316_fp_cases_for_mnp_screen_01.vcf](/big8_disk/liaoyoyo2001/InterSubMod/docs/reports/validated/2026/03/assets/20260316_fp_cases_for_mnp_screen_01.vcf)
- 4 個 FP 的 MNP 嚴格檢查結果：[20260316_fp_mnp_screen_result_01.tsv](/big8_disk/liaoyoyo2001/InterSubMod/docs/reports/validated/2026/03/assets/20260316_fp_mnp_screen_result_01.tsv)

### 2.4 實際使用的命令

```bash
python3 scripts/analysis/validate_igv_case_loci.py \
  --output-tsv docs/reports/validated/2026/03/assets/20260316_igv_case_validation_01.tsv

python3 scripts/analysis/screen_mnp_adjacent_fp.py \
  --fp-vcf docs/reports/validated/2026/03/assets/20260316_fp_cases_for_mnp_screen_01.vcf \
  --tumor-bam /big8_disk/liaoyoyo2001/data/bam/HCC1395_Tmode_tagged_ClairS_pileup_v040_woFilter.bam \
  --normal-bam /big8_disk/liaoyoyo2001/data/bam/HCC1395BL_ONT_5khz_simplex_5mCG_5hmCG_tagged.bam \
  --output-tsv docs/reports/validated/2026/03/assets/20260316_fp_mnp_screen_result_01.tsv
```

## 3. 9 個案例總覽

| Case | 位點 | 圖檔 | AI ALT 成群 | AI normal 是否同構 | AI 相鄰異常柱 | AI 初判驅動 | AI coverage / alignment 異常 | 正式驗證重點 | 本輪結論 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| TP_01 | `chr6:145444893 G>A` | [圖](assets/20260316_subclone_phase4_igv_sessions/snapshots/TP_01_chr6_145444893.png) | 清楚 | 否 | 否 | allele | 否 | exact truth SNV；tumor ALT=`32/92`；ALT HP=`2-1:32` | AI 初篩可靠，屬乾淨 allele-driven TP |
| TP_02 | `chr4:70548355 G>A` | [圖](assets/20260316_subclone_phase4_igv_sessions/snapshots/TP_02_chr4_70548355.png) | 中等 | 否 | 否 | allele | 中等，phase 資訊不足 | exact truth SNV；tumor ALT=`24/81`；ALT HP=`1-1:24`；REF 側大量 `NA` HP | AI 可抓方向，但需寫成 allele-supported、phase-limited |
| TP_03 | `chr5:153209947 C>A` | [圖](assets/20260316_subclone_phase4_igv_sessions/snapshots/TP_03_chr5_153209947.png) | 弱 | 部分相似 | 否 | 弱 allele | 否 | exact truth SNV；tumor ALT=`30/92`；`AlleleDelta=0.0148` | 主要靠統計，不適合只憑圖定論 |
| TP_04 | `chr16:35118902 G>A` | [圖](assets/20260316_subclone_phase4_igv_sessions/snapshots/TP_04_chr16_35118902.png) | 清楚 | 否 | 否 | allele | phase 複雜，但非明顯對齊異常 | exact truth SNV；tumor ALT=`35/132`；ALT HP=`1-1:15,3:20` | AI 初篩可靠，是最像真實亞克隆重程式化的 TP |
| TP_05 | `chr7:109185781 G>T` | [圖](assets/20260316_subclone_phase4_igv_sessions/snapshots/TP_05_chr7_109185781.png) | 中等 | 部分 | 否 | HP > allele | 否 | exact truth SNV；`DominantLabel=hp`；`LabelHPP=1e-3` | AI 可先懷疑 HP confounding，正式驗證支持此判讀 |
| FP_A1 | `chr8:93565727 C>T` | [圖](assets/20260316_subclone_phase4_igv_sessions/snapshots/FP_A1_chr8_93565727.png) | 否 | 否 | 否 | none | 有，低 ALT 支持 + 高未相位味道 | 無 truth；tumor ALT=`5/112`；`CramersV=0` | AI 初篩可靠，屬低 VAF / ALT 不成群型 FP |
| FP_A2 | `chr9:137953060 T>C` | [圖](assets/20260316_subclone_phase4_igv_sessions/snapshots/FP_A2_chr9_137953060.png) | 否 | 近似無差 | 否 | 視覺上不像 allele | 否明顯 | 無 truth；tumor ALT=`9/84`；`LabelAlleleP=0.01` 但 `CramersV=0` | 典型高維假顯著，必須靠效應量而不是看圖 |
| FP_B1 | `chr7:52087777 A>T` | [圖](assets/20260316_subclone_phase4_igv_sessions/snapshots/FP_B1_chr7_52087777.png) | 有，但跟 HP 走 | 部分 | 是 | hp | 有，鄰近 truth INDEL | 無 exact truth SNV；`chr7:52087776:TA>T` nearby truth INDEL；`DominantLabel=hp` | HP-driven + INDEL-adjacent artifact 的證據充足 |
| FP_B2 | `chr9:75383880 T>A` | [圖](assets/20260316_subclone_phase4_igv_sessions/snapshots/FP_B2_chr9_75383880.png) | 有，小群 | 否 | 是 | hp + local alignment anomaly | 有 | 無 truth；tumor ALT=`15/106`；`adj_minus1_alt_del_frac=1.000`；MNP screen `0 candidates` | **本輪修正：更像 `target SNV + -1 deletion-like gap`，不是乾淨 MNP** |

## 4. 單張圖片觀察

### 4.1 TP_01 `chr6:145444893 G>A`

![TP_01 full template](assets/20260316_subclone_phase4_igv_sessions/snapshots/TP_01_chr6_145444893.png)

觀察重點：

1. 主 `tumor BAM` 可見 target 支持 reads 形成可辨識小群。
2. `normal BAM` 沒有對應 ALT 結構，因此偏向 tumor-specific。
3. 這張圖適合先判成 `allele-driven`。

### 4.2 TP_02 `chr4:70548355 G>A`

![TP_02 full template](assets/20260316_subclone_phase4_igv_sessions/snapshots/TP_02_chr4_70548355.png)

觀察重點：

1. 圖上有一群較一致的 target 支持 read，但邊界不如 `TP_01` 清楚。
2. 問題更像 phase 背景不足，不是完全沒有訊號。
3. 適合先標成 `partial ALT clustering / phase-limited`。

### 4.3 TP_03 `chr5:153209947 C>A`

![TP_03 full template](assets/20260316_subclone_phase4_igv_sessions/snapshots/TP_03_chr5_153209947.png)

觀察重點：

1. 肉眼很難抓到強烈 ALT 小群。
2. 這類型更像低甲基背景上的小幅群體位移。
3. 單靠圖不適合直接下 TP 結論。

### 4.4 TP_04 `chr16:35118902 G>A`

![TP_04 full template](assets/20260316_subclone_phase4_igv_sessions/snapshots/TP_04_chr16_35118902.png)

觀察重點：

1. 這張圖最容易看出兩層 read 結構。
2. `normal BAM` 沒有同樣的小群，因此很像 tumor-specific。
3. 這是最適合做 `AI 先判 + 正式驗證再確認` 的代表 TP。

### 4.5 TP_05 `chr7:109185781 G>T`

![TP_05 full template](assets/20260316_subclone_phase4_igv_sessions/snapshots/TP_05_chr7_109185781.png)

觀察重點：

1. target 支持 read 比較像跟某個 HP 背景一起移動。
2. 因此這張圖適合先標成 `HP confounding suspect`。
3. 正式驗證的 `DominantLabel=hp` 與這個方向一致。

### 4.6 FP_A1 `chr8:93565727 C>T`

![FP_A1 full template](assets/20260316_subclone_phase4_igv_sessions/snapshots/FP_A1_chr8_93565727.png)

觀察重點：

1. 圖上看不出可信 ALT 小群。
2. 這類型常是低 VAF / 低功效，而不是穩定亞群。
3. AI 在這型 FP 的辨識通常可靠。

### 4.7 FP_A2 `chr9:137953060 T>C`

![FP_A2 full template](assets/20260316_subclone_phase4_igv_sessions/snapshots/FP_A2_chr9_137953060.png)

觀察重點：

1. 視覺上沒有明確 ALT 群。
2. 這種 case 優先看效應量，不要先看 p-value。
3. 本輪正式驗證把它歸到高維假顯著型 FP。

### 4.8 FP_B1 `chr7:52087777 A>T`

![FP_B1 full template](assets/20260316_subclone_phase4_igv_sessions/snapshots/FP_B1_chr7_52087777.png)

觀察重點：

1. 分層雖然明顯，但更像 HP 背景，不像 ALT 自己形成新層。
2. 再看鄰近 annotation，可發現它緊貼 truth INDEL。
3. 這是 `HP-driven + INDEL-adjacent` 的標準反例。

### 4.9 FP_B2 `chr9:75383880 T>A`

![FP_B2 full template](assets/20260316_subclone_phase4_igv_sessions/snapshots/FP_B2_chr9_75383880.png)

觀察重點：

1. 圖上確實有「target 附近不只一格怪」的感覺。
2. 但這種圖像只能先標成 `adjacent anomaly suspect`，不能直接寫成 MNP。
3. 正式驗證顯示這條主 BAM 路徑更支持 `-1 deletion-like gap`。

## 5. 重要驗證結論

### 5.1 AI 真的可以先抓到的現象

1. **ALT 是否成群**
   - `TP_01`、`TP_04`、`FP_B1`、`FP_B2` 都能先從圖上看出 target 支持 reads 不是完全散開的。
   - `FP_A1` 則能先看出 ALT 太少，不足以自成可信群。
2. **HP-driven vs allele-driven**
   - `TP_05`、`FP_B1` 只看圖就會讓人先懷疑 HP 共線，正式驗證也支持這個方向。
   - `TP_01`、`TP_04` 看起來則更接近 allele-driven，正式驗證同樣支持。
3. **強烈的局部異常**
   - `FP_B1` 的鄰近 truth INDEL 問題
   - `FP_B2` 的 `-1` 格 deletion-like 訊號
   這兩型都適合先由 AI 從圖上標示為「需要 read-level 深查」。

### 5.2 AI 不適合單獨定論的現象

1. **弱效應**
   - `TP_03` 的肉眼分離度很低，但 exact truth SNV 與 `LabelAlleleP=0.003` 仍支持它是 TP。
2. **高維假顯著**
   - `FP_A2` 視覺上缺乏明確 ALT 群，正式驗證顯示 `CramersV=0`、`AlleleDelta=-0.0029`，只有 label-aware p-value 還顯著。
   - 這種 case 若只看 p-value，容易被過度解讀。

### 5.3 本輪最重要的修正：`FP_B2` 不是乾淨 MNP

先前的工作假設是：`chr9:75383880` 可能是 `75383879-75383880` 的相鄰雙 SNV / MNP。  
但本輪用主 tumor BAM 做嚴格 read-level 驗證後，得到的是：

1. target `75383880 T>A` 在 tumor ALT reads 上明確存在：
   - `tumor ALT = 15/106`
   - `tumor ALT HP = 1-1:15`
2. ALT reads 在 `75383879` 的主訊號不是穩定的非 reference base，而是：
   - `adj_minus1_alt_del_frac = 1.000`
3. 用 [screen_mnp_adjacent_fp.py](/big8_disk/liaoyoyo2001/InterSubMod/scripts/analysis/screen_mnp_adjacent_fp.py) 對 4 個 FP 做嚴格 MNP 篩查：
   - [20260316_fp_mnp_screen_result_01.tsv](/big8_disk/liaoyoyo2001/InterSubMod/docs/reports/validated/2026/03/assets/20260316_fp_mnp_screen_result_01.tsv)
   - 結果是 **`0 candidates`**

**因此在這條主 BAM 路徑下，`FP_B2` 應先改寫成：**

- `target SNV + -1 deletion-like gap`
- `HP-driven / local alignment anomaly`

而不是直接寫成「已確認為 MNP」。

### 5.4 `FP_B1` 的機制在這輪被再次強化

`FP_B1` 的正式驗證與先前判讀一致，而且更穩：

1. 無 exact truth SNV
2. `±2 bp` 內有 nearby truth INDEL：
   - `chr7:52087776:TA>T`
3. feature matrix 上：
   - `DominantLabel=hp`
   - `CramersV=1.0`

這使它非常適合作為 `HP-driven + INDEL-adjacent artifact` 的標準反例。

## 6. 對未來 IGV 自動截圖流程的操作建議

### 6.1 建議的固定流程

1. 先用固定 template 批次出 full-template PNG。
2. 用固定欄位做 AI 視覺初篩：
   - ALT 是否成群
   - normal 是否同構
   - 是否有相鄰異常柱
   - HP 還是 allele 驅動
   - 是否有 coverage / alignment anomaly
3. 對 AI 標成 `強訊號` 或 `局部異常` 的位點，立刻跑：
   - [validate_igv_case_loci.py](/big8_disk/liaoyoyo2001/InterSubMod/scripts/analysis/validate_igv_case_loci.py)
4. 若懷疑相鄰 MNP / indel-like 事件，再補：
   - [screen_mnp_adjacent_fp.py](/big8_disk/liaoyoyo2001/InterSubMod/scripts/analysis/screen_mnp_adjacent_fp.py)

### 6.2 AI 適合扮演的角色

AI 最適合做：

1. **快速 triage**
   - 幫你先把 `清楚 allele-driven`、`可能 HP-driven`、`相鄰異常需要深查` 分開
2. **把看圖語言標準化**
   - 例如把「有點像一群但不太乾淨」固定翻成 `partial ALT clustering`
3. **提醒哪裡不能只靠圖**
   - 弱效應 case
   - 高維假顯著 case
   - 相鄰柱異常但未經 aligned-pairs / pileup 驗證的 case

### 6.3 不該只靠 AI 看圖就下結論的項目

1. 真正的 TP / FP 定性
2. MNP 與 indel adjacency 的正式分類
3. LOH / CNV 這類區域層級機制
4. 弱效應位點的生物意義

## 7. 本輪結論整理

1. **workflow 本身可行**：固定 template PNG → AI 視覺初篩 → BAM/VCF/統計值正式驗證，這條鏈已能落地。
2. **AI 不是最終證據**：它能先抓出 `TP_01`、`TP_04`、`FP_A1`、`FP_B1` 這類強圖像訊號，但不能單獨判定 `TP_03`、`FP_A2`。
3. **正式驗證能真正修正圖像直覺**：`FP_B2` 就是最好的例子。IGV 視覺上像相鄰雙變異，但主 tumor BAM 的嚴格檢查更支持 `-1 deletion-like gap`。
4. **因此之後的標準寫法應是**：
   - `圖片觀察` 寫圖片觀察
   - `BAM/VCF/統計值驗證` 另列為正式驗證
   - 兩者不要混寫成同一層結論
