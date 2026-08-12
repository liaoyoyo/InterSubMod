<!--
建立時間: 2026-08-10
狀態: in_progress
任務類型: B comprehensive validation（HCC1395 全基因組 chr1-22）
目標: 量化 HCC1395 subclonal reconstruction 結果受 CNV 影響的程度；以 SEQC2 為結對真值、以 SAVANA 為流程驗證；設計矯正方案；裁決既有結論是否仍成立
處理範圍: HCC1395 canonical 2026-07-24 frozen cohort（chr1-22）。不重算 CCF、不重跑樹枚舉、不修改任何 canonical 輸出
data_sources: ../data/cn_annotation_summary.json,../data/af_vs_cn_expectation.json,../data/savana_vs_seqc2_adjudication.json,../data/savana_refit_grid.json,../data/purity_selfconsistency_audit.json,../data/resolution_vs_cn_stratified.json,../data/mechanism_and_intervals.json,../data/mechanism_v2_distinctness.json,../data/robustness_checks.json,../data/consolidated_findings.json
驗證方式: 全部數字由 scripts/ 下腳本產出至 data/*.json 後讀回；分母與 20260801 authority_manifest / denominator_registry 對帳
-->

# CNV 混淆稽核：SEQC2 結對驗證與 SAVANA 流程驗證（HCC1395 全基因組）

## 0. 一句話結論

**共現與拓撲骨幹沒有被 CNV 抬高，可以保留；但「read-AF 挑出唯一最佳樹」這一層被 CNV 顯著抬高，
其中約 15–24% 的「唯一樹」是拷貝數的產物而非譜系訊號。**
機制已釐清：CN 不是讓 read-AF 排錯順序，而是決定它**有沒有東西可排**。
既有結論不需要撤回，但「單一拓撲率」這類數字必須改標為 CN-conditional，
且不得在未整合 allele-specific CN 的前提下升格為生物學宣稱 —— 這與 2026-08-01 規格
`§5.2` 已規定的 `RAW_AF_UNIQUE_REPRESENTATIVE` 上限一致。

---

## 1. 這次為什麼要重做（前一輪的三個問題）

2026-07-10 曾做過一輪 CN confound 檢查，得到「共現區 94.2% CN-altered、巢狀 subclone 區 98%」。
本輪重做，因為那一輪有三個具體問題：

| # | 前一輪的做法 | 問題 | 本輪的做法 |
|---|---|---|---|
| 1 | 骨幹用 23,810 sSNV 的 `mlhp_part_*.json` | 該骨幹已被 20260801 權威登錄列為 **historical**；canonical 已改為 469,849 dataset-records / 98,955 units | 直接讀 canonical frozen cohort `HCC1395.topology.jsonl` |
| 2 | CN 只用 categorical 的 `ngs_benchmark_cnvs_gain_loss_loh.bed`（660 段，只有 gain/loss/loh 標籤） | 同目錄下**一直有整數 CN**（`ngs_benchmark_cnv_gain_cn.bed` 1,502 段、`ngs_benchmark_cnv_loss_cn.bed` 63 段），當時未使用，導致「要真 CCF 需整數 CN」被誤認為原料缺失 | 使用整數 CN，並把 LOH 分離成獨立軸 |
| 3 | 用「區域中點」的單一 CN 標籤 | categorical 檔有 110 處重疊（gain 與 loh 共存＝GAINLOH），中點法會丟失共存資訊 | 以 bp 加權涵蓋率指派，跨斷點者標記為 mixed |

補充：第 3 點在本輪實測影響很小 —— 9,624 個可定位 unit 中只有 **21 個**跨 CN 斷點，
因為 unit 的基因體跨度遠小於 CN 段。**這一項前一輪其實沒有實質錯誤**，據實記錄。

---

## 2. 分析流程（可重跑）

```
canonical frozen cohort (2026-07-24)
  HCC1395.topology.jsonl   11,590 units  ← chrom + active_positions + af_coverage(alt/ref/exact fraction)
  HCC1395.census.jsonl      9,130 units  ← resolution_class / topology signature
        │
        ├── cn_annotate_units.py ──────────► 每 unit / 每 sSNV 的 CN 註釋
        │      SEQC2 整數 CN（gain 1,502 段 + loss 63 段，其餘視為 CN=2）
        │      SEQC2 LOH（獨立軸）+ exclusion.bed
        │      SAVANA WGS CN（獨立第二來源）
        │
        ├── af_vs_cn_expectation.py ───────► read-AF 落在 m/c 格點的程度（含隨機命中基線）
        ├── savana_vs_seqc2_adjudication.py ► SAVANA vs SEQC2 段層一致性 + 回歸 + 逐染色體偏差
        ├── savana_refit_grid.py ──────────► 用 SAVANA 自己的 log2r 做 purity×ploidy grid search
        ├── purity_selfconsistency_audit.py ► 逐段 BAF 上限違反率（可套用到無真值的樣本）
        ├── resolution_vs_cn_stratified.py ─► 唯一性拆兩層 + 候選樹數/k 分層 CMH
        ├── mechanism_v2_distinctness.py ──► 機制檢驗（AF 相異性）
        └── robustness_checks.py ──────────► 染色體混淆檢驗 + 套套邏輯檢驗
```

分母對帳：`ranked 9,130 + zero_denominator 227 + no_active_alt 1,966 + family_incomplete 267 = 11,590` ✓
其中 1,966 個 `no_active_alt` 無 active position 故無法定位，可定位者 **9,624**。

> **分母警語**：`denominator_registry.tsv` 的 98,955 / 71,955 / 63,506 等是**全七樣本**合計。
> 本報告所有比例的分母都是 **HCC1395 單樣本**，不可與登錄表數字直接相除。

---

## 3. 問題一：有多少坐在 CNV 區？

| 層級 | CN-altered | 分母 | 比例 |
|---|---|---|---|
| unit（區域） | 9,185 | 9,624 | **95.44%** |
| sSNV site | 18,833 | 20,060 | **93.88%** |

分類明細：

| CN 類別 | unit 數 | unit % | site 數 | site % |
|---|---:|---:|---:|---:|
| neutral（CN=2 且無 LOH） | 439 | **4.56%** | 1,227 | **6.12%** |
| gain（CN>2，無 LOH） | 5,970 | 62.03% | 10,617 | 52.93% |
| gain + LOH | 1,683 | 17.49% | 4,687 | 23.36% |
| copy-neutral LOH | 1,459 | 15.16% | 3,375 | 16.82% |
| loss / loss+LOH | 73 | 0.76% | 154 | 0.77% |

sSNV 落點的整數 CN 以 **CN=3（5,877 個）與 CN=4（5,184 個）為主**，CN=2 僅 4,602 個，
最高到 CN=9.5。HCC1395 是高度非整倍體樣本，這個結果本身不意外。

### 分布極度不均（這是後面推論的關鍵限制）

| 染色體 | unit 數 | neutral | neutral % |
|---|---:|---:|---:|
| **chr21** | 114 | 104 | **91.23%** |
| chr2 | 805 | 159 | 19.75% |
| chr6 | 384 | 71 | 18.49% |
| chr16 | 293 | 51 | 17.41% |
| chr4 | 776 | 29 | 3.74% |
| **chr7 / chr10 / chr12 / chr14 / chr18 / chr19** | 2,553 | **0** | **0.00%** |

chr8 的 613 個 unit 中有 522 個是 gain+LOH（85.15%），與既知的 chr8 全 LOH 一致。

**限制**：可用於比較的 CN-neutral 樣本高度集中 —— 在 ranked 且有 ≥2 個 AF 值的 351 個
neutral unit 中，chr2(43.59%)＋chr21(29.63%)＋chr6(8.55%)＋chr4(7.69%) 就佔 **89.46%**。
因此「neutral vs altered」的對比天生帶有染色體混淆，§6 對此做了專門檢驗。

---

## 4. 問題二：拓撲／共現骨幹受影響嗎？——**沒有**

把「唯一性」拆成兩層是本次分析的關鍵：

| 層 | 定義 | 是否受 CN 影響 |
|---|---|---|
| **結構唯一** | `total_tree_count == 1`，枚舉本身就只有一棵樹，read-AF 沒有參與 | 理論上 CN-robust |
| **AF 破 tie** | `total_tree_count > 1` 且 `best_tree_unique`，多棵結構上合法的樹由 read-AF 選出贏家 | 正是 CN 能污染的一層 |

結果：

| | 全部 ranked | CN-neutral | CN-altered |
|---|---:|---:|---:|
| n | 9,130 | 351 | 8,779 |
| **結構唯一率** | 49.11% | **56.41%** | **48.82%** |
| unique_best_tree 率 | 77.19% | 67.52% | 77.57% |

**結構唯一率在 CN-neutral 反而更高**（56.41% > 48.82%）。也就是說，共現連鎖本身＋最小樹枚舉
這一層沒有因為 CNV 而變得「更容易得到唯一解」，方向甚至相反。
**這支持既有的核心立場：sSNV 單分子共現骨幹是 CN-robust 的，拓撲軸存活。**

---

## 5. 問題三：read-AF 這一層受影響嗎？——**顯著受影響**

只看 read-AF 真正有作用空間的 unit（候選樹 > 1，共 4,646 個）：

| CN 類別 | n | read-AF 成功破 tie |
|---|---:|---:|
| **neutral** | 153 | **25.49%**（95% CI 19.24–32.94） |
| loss_loh | 40 | 12.5% |
| loss | 17 | 23.53% |
| cnloh | 1,001 | 45.75% |
| gain | 2,391 | 51.32% |
| **gain + LOH** | 1,044 | **79.5%** |
| 合計 CN-altered | 4,493 | **56.18%**（95% CI 54.72–57.62） |

差距 **+30.69 個百分點**，Fisher OR = 3.747（p = 6.5e-14）。

### 分層控制後仍成立

局部突變密度是本專案已知的必要控制變項（2026-07-26 甲基那次就是分層後訊號全滅），
因此逐項控制：

| 分層變項 | CMH 合併 OR | p |
|---|---:|---|
| 候選樹數 | 3.5307 | 8.2e-12 |
| active bit count k | 3.7282 | 1.6e-13 |
| 候選樹數 × k | 3.5317 | 8.1e-12 |
| **染色體**（見 §6） | **2.1612** | 2.7e-04 |

前三項幾乎不衰減（3.747 → 3.5307），**與甲基那次的崩解模式相反，這次是真的**。
但按染色體分層後縮小到 2.16 —— 效應存活，量級減半，據實記錄於 §6。

---

## 6. 問題四：機制是什麼？（含一次被對抗審查抓出的計算錯誤）

### 6.0 先說修正：初版把不進分數的位點算了進去

初版在算「unit 內 AF 的離散度與相異性」時，遍歷了 `af_coverage` 的**全部**列。
但 read-AF 分數只看 `active_positions`；`af_coverage` 另外帶有非 active 位點的全參考列
（`fraction = 0/1`），那些列永遠不會進入 `s(p,c)`。9,130 個 ranked unit 中有 4,734 個帶這種
多餘列（最多 +11 列），於是「分數上根本不可能被解開」的 unit 被誤標為「AF 相異」，
而且兩個 CN 組受影響的程度不同。

此錯誤由 2026-08-10 的對抗審查抓出（49 項發現中唯一被確認為真者）。
修正＝把相異性與離散度都限制在 active 位點上。**修正後兩項檢驗的結論都改變，且方向一致。**

| 檢驗 | 初版（含非 active 位點） | 修正後（僅 active 位點） |
|---|---|---|
| AF 離散度 altered vs neutral | 0.5223 < 0.7826，p = 0.9997 → **推翻** | **0.2265 > 0.0，p = 4.0e-10 → 成立** |
| AF 全同比例 neutral vs altered | 20.80% vs 11.78%，OR 1.9669 | **58.76% vs 25.23%，OR 4.2224，p = 2.2e-20** |
| 中介檢驗 CMH | 2.819，p = 3.9e-06 → 部分中介 | **1.2761，p = 6.7e-01 → 完全中介** |

**未受影響**：§3 CN 佔比、§4 結構層、§5 主要關聯與其分層、§7 超額估計、§8–§9 SAVANA。
這些都不使用相異性旗標，已逐項確認。

### 6.1 機制：CN 決定 AF 值在算術上是否可區分

read-AF 以 exact rational 比較，tie 發生在候選樹**分數相等**時。
分數 `s(p,c) = Σ_{i∈p}(AF_i − AF_j(c))`：若一個 unit 的所有 active AF 都等於同一個值 v，
每一項都是 0、所有候選樹分數皆為 0 → **tie 是數學必然**。
實測驗證：1,569 個 active AF 全同的 unit，`best_score_fraction` 全部是 `0/1`，零例外。

在純腫瘤細胞株的 CN-neutral 區，單 copy haplotype 上的 clonal 突變讀出恰好 1/1，
同一 unit 內多個位點共享同值；CN gain/LOH 下同一個 clone 的突變落在不同的 m/c 格點，值就分開了。

| CN 類別 | unit 內 AF 中位離散度 |
|---|---:|
| neutral | **0.0** |
| gain | 0.1833 |
| gain + LOH | 0.2629 |
| copy-neutral LOH | 0.3058 |

neutral 的中位離散度為 0，正是因為過半 unit 的 active AF 完全相同。

### 6.2 完全中介：CN 的效應全部經由「算術可解性」

在 contested unit（候選樹 > 1）中，**算術上不可能被 read-AF 解開**的比例：

| | 比例（95% CI） | n |
|---|---|---:|
| CN-neutral | **67.97%**（60.22–74.85） | 153 |
| CN-altered | **32.54%**（31.18–33.92） | 4,493 |

一旦排除這些算術不可解的 unit，剩下的破 tie 率：

| | 比例（95% CI） | n |
|---|---|---:|
| CN-altered | 83.27%（81.9–84.56） | 3,031 |
| CN-neutral | 79.59%（66.36–88.52） | 49 |

Fisher OR = 1.2765，p = 0.4459；控制相異性與候選樹數的 CMH OR = 1.2761，p = 6.7e-01。
**在算術可解性之外，偵測不到殘餘的 CN 效應。**

這比初版的「部分中介、機制不明」乾淨得多，而且是**強化**主結論而非削弱：
CN 不是讓 read-AF「排錯順序」，而是決定 read-AF **有沒有東西可排**。
CN-neutral 區有三分之二的 unit 在算術上本來就無解，CN-altered 區把這個比例砍半 ——
多出來的那些「可解」單元，其可解性來自拷貝數結構，不是來自譜系。

⚠ 非退化子集的 n_neutral 僅 49，區間寬（66.36–88.52），OR = 1 並未被證明。
只能說「無法宣稱存在殘餘效應」，不能說「已證明沒有殘餘效應」。

### 6.3 AF 有多少能被 CN 單獨解釋

以 m/c 格點命中率衡量（容忍度 0.05，並扣除格點密度造成的隨機命中）：

| CN 類別 | 命中率 | 隨機基線 | 超額 |
|---|---:|---:|---:|
| neutral | 54.34% | 5.0% | +49.34% |
| cnloh | 63.66% | 15.0% | +48.66% |
| gain | 58.82% | 16.53% | +42.29% |
| gain+LOH | 56.28% | 33.56% | +22.72% |

各類別都遠超隨機。注意 neutral 的格點只有 {1/1} 一點，其 54.34% 命中即「AF≈1.0 的 clonal 位點」，
反過來說 **neutral 區約 45.7% 的位點 AF 明顯低於 1.0，那才是可乾淨解讀的 subclonal 候選訊號**。

---

## 7. 問題五：影響量化 —— 有多少「唯一樹」是 CNV 的產物？

若 CN-altered 區的破 tie 率應與 CN-neutral 區相同：

| | 數量 |
|---|---:|
| CN-altered 實際破 tie | 2,524 |
| 以 neutral 率預期 | 1,145.3 |
| **超額（CN 可歸因）** | **1,378.7**（範圍 1,044.0 – 1,659.5） |
| 佔 7,047 個 unique_best_tree unit | **19.56%**（範圍 14.81% – 23.55%） |

**約每 5 棵「唯一最佳樹」就有 1 棵，其唯一性來自拷貝數而非譜系。**

⚠ 此估計的三個前提必須同時聲明：
1. 假設 CN-neutral 的破 tie 率就是「無 CN 基線」，但 neutral 區仍可能含未報告的 subclonal CN；
2. neutral contested 樣本只有 153 個，區間寬；
3. 區間僅傳播 neutral 率的 Wilson 區間，未涵蓋其他不確定性。
   若改採染色體分層後的 OR 2.16，超額比例會低於上述點估計 —— 上表應視為**上緣估計**。

---

## 7A. 非 CN 區還剩多少資訊？——「乾淨」比 4.56% 還要少

「乾淨區佔 4.56%」是**區域**的比例，不是**可用資訊**的比例。要建樹至少需要 2 個 sSNV：

| | ranked unit | 其中 k=1（單一 sSNV，**數學上不可能建樹**） | k≥2（真正有樹結構資訊） |
|---|---:|---:|---:|
| CN-neutral | 351 | **174**（49.57%） | **177**（50.43%） |
| CN-altered | 8,779 | 2,972 | 5,807（66.15%） |

→ **真正能提供樹結構資訊的乾淨 unit 只佔 177 / 5,984 = 2.96%。**
且 CN-neutral 的 k 分布明顯偏小（k=1 佔一半，altered 只佔三分之一），
因為乾淨區多半是散落的單點，不是密集共現群。

**染色體代表性**：這 177 個分布在 16 條染色體，但 chr2(64) + chr21(49) 就佔 **63.84%**。
不能當成全基因組的隨機樣本。

**統計檢定力**：以現有 n（contested neutral 153 vs altered 4,493）、基線 25.49%、α=0.05、power=0.8，
能偵測到的最小效應是 **36.25%（OR 1.662）**。**小於 OR 1.66 的 CN 效應，這個樣本量看不到** ——
所以 §6.2「排除算術不可解後偵測不到殘餘效應」必須理解為「偵測不到」，不是「不存在」。

### 但有一項乾淨區的資訊是不可替代的

CN-neutral 非 LOH、純腫瘤細胞株、單 copy haplotype 上的 clonal 突變，read-AF 應讀出 1/1。
因此**低於 1 的 AF 在這裡是唯一沒有拷貝數解釋的 subclonal 讀數**：

| CN-neutral 的 1,227 個 sSNV | 數量 | 比例（95% CI） |
|---|---:|---|
| clonal（read-AF ≥ 0.95） | 787 | 64.14%（61.42–66.78） |
| **subclonal（read-AF < 0.95）** | **440** | **35.86%（33.22–38.58）** |

**這 440 個位點是全樣本中唯一能在「無拷貝數解釋」前提下讀出 subclonal 分率的證據。**
數量不大，但它們不是殘渣 —— 是目前唯一乾淨的頻率軸來源。

---

## 7B. CN 區的特殊現象（三項，都不是「只是被污染」）

### 現象一：read-AF 隨拷貝數單調下降，且下限精確落在 1/CN

| total CN | n（sSNV） | read-AF 中位 | 模型預測最低格點 1/CN | AF 落在 m/c 格點 ±0.05 內 |
|---:|---:|---:|---:|---:|
| 1 | 151 | 1.0 | 1.0 | 86.09% |
| 2 | 3,375 | 0.931 | 0.5 | 64.56% |
| 3 | 5,877 | 0.8788 | 0.3333 | 67.19% |
| 4 | 5,184 | 0.5769 | 0.25 | 65.7% |
| 5 | 2,422 | 0.4634 | 0.2 | 65.94% |
| 6 | 1,175 | 0.365 | 0.1667 | 71.06% |
| 7 | 420 | 0.1918 | 0.1429 | 75.71% |
| 8 | 205 | 0.2169 | 0.125 | 85.37% |

中位 AF 從 CN=1 的 1.0 一路降到 CN=7 的 0.19，**而且格點吻合率在高 CN 端反而上升**
（CN=8 時 85.37%）。這是 mutation multiplicity 的直接讀數：拷貝數越高，
AF 越被鎖在 m/c 的離散格點上、越不可能是連續的細胞分率。
**在 CN≥4 的區域把低 AF 解讀成「小 subclone」，幾乎一定是誤讀。**

### 現象二：LOH 有雙重且方向相反的效應

| CN 類別 | 結構唯一率 | contested 中 AF 破 tie 率 |
|---|---:|---:|
| gain（無 LOH） | **58.48%** | 51.32% |
| neutral | 56.41% | 25.49% |
| loss | 37.04% | 23.53% |
| gain + LOH | 31.54% | **79.5%** |
| copy-neutral LOH | **29.71%** | 45.75% |
| loss + LOH | **9.09%** | 12.5% |

**帶 LOH 的區段結構唯一率全面偏低（9–32%），沒有 LOH 的反而高（56–58%）。**
也就是 LOH 讓**結構**更難定（候選樹更多），卻同時讓 **read-AF** 更容易挑出唯一解 ——
gain+LOH 兩者兼具，於是成為「候選樹多、又幾乎總能被 AF 選掉」的最危險組合（79.5%）。
機制上合理：LOH 使一條 haplotype 消失，該 HP 家族的 read 支持變薄、partial pattern 變多（結構更糊），
同時多個 copy 集中在單一 haplotype 上、AF 落在不同 m/c（AF 更能分辨）。

### 現象三：CN **不**改變樹的形狀 —— 控制 k 之後差異消失

分支（Sister 類）比例，限 k≥2：

| | n | 分支比例（95% CI） |
|---|---:|---|
| CN-neutral | 177 | 15.25%（10.7–21.28） |
| CN-altered | 5,807 | 20.51%（19.49–21.57） |

**信賴區間重疊**；逐 k 分層後更清楚：

| k | CN-neutral | CN-altered |
|---:|---|---|
| 2 | 11.19%（6.9–17.65） | 13.59%（12.57–14.68） |
| 3 | 29.03%（16.1–46.59） | 30.43%（27.92–33.08） |
| 4 | 37.5%（13.68–69.43） | 47.11%（42.14–52.13） |

每一層都重疊。**拷貝數不改變重建出來的樹長什麼樣，只改變我們能不能從候選中挑出唯一一棵。**
未控制 k 時看起來有差（neutral 的 `Single-only` 49.57% 遠高於 altered 33.85%），
那純粹是 k=1 佔比不同造成的假象。

> **三層總結**：形狀 CN-robust（現象三）→ 結構唯一性 CN-robust（§4）→ **AF 挑唯一樹被 CN 抬高（§5）**。
> 受影響的只有最後一層。

---

## 7C. 扣掉 CN 區之後，還夠不夠確認樹的演化？——**不夠，而且比表面數字更不夠**

要判斷「夠不夠」，`unit 數`不是正確的單位。一個 unit 要真的貢獻演化資訊，必須同時滿足兩件事：

1. **k ≥ 2** —— k=1 只有一個突變，樹必然是 `ROOT→A`，**沒有任何先後順序資訊**
2. **樹被定出來** —— 分成 `結構唯一`（枚舉即唯一，read-AF 沒參與）與
   `AF 定出`（多候選由 read-AF 選一棵；**在 CN 中性區這個選擇沒有拷貝數混淆，因此可用**）

真正的計量單位是 **deep edge**：parent 不是 ROOT 的邊，也就是真正表達「兩個突變誰先誰後」的關係。

### 逐層拆解（CN 中性）

| | unit 數 | 結構唯一 | AF 定出 | 仍 tied | **可定樹合計** | 總邊 | **deep edge** | 染色體 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| 全部 ranked | 351 | 198 | 39 | 114 | 237（67.52%） | 318 | **72** | 16 |
| **k ≥ 2** | **177** | 24 | 39 | 114 | **63（35.59%）** | 144 | **72** | 16 |
| **k ≥ 3** | **43** | 1 | 13 | 29 | **14（32.56%）** | 46 | **29** | 5 |

逐 k 看更清楚：

| k | unit | 可定樹率 | deep edge |
|---:|---:|---:|---:|
| 1 | 174 | 100% | **0** ← 全部是空的 |
| 2 | 134 | 36.57% | 43 |
| 3 | 31 | 35.48% | 19 |
| ≥4 | 12 | 25.0% | 10 |

**「乾淨區可定樹率 67.52%」是假象** —— 那個數字被 174 個 k=1 撐起來，而 k=1 貢獻 **0 條** deep edge。
扣掉之後可定樹率掉到 **35.59%**。

### 與 CN 變異區對照

| | unit（k≥2） | 可定樹率 | deep edge |
|---|---:|---:|---:|
| CN 中性 | 177 | **35.59%**（28.91–42.88） | **72** |
| CN 變異 | 5,807 | **66.09%**（64.86–67.3） | **4,547** |

→ **CN 中性區只佔全部 deep edge 的 1.56%。**

### 🔑 一個反直覺但關鍵的發現：乾淨區反而更定不出樹

CN 中性區的 `仍 tied` 比例是 **64.41%**，CN 變異區只有 33.91%。原因正是 §6 的機制：
乾淨區的 read-AF 大量等於 1/1、彼此相同 → 分數相同 → **數學上無法破 tie**。
CN 反而「提供」了打破平手所需的數值差異 —— 只是那個差異來自拷貝數，不是時間。

**所以這是一個取捨，不是單純的好壞：**

| | 可辨識性 | 可信度 |
|---|---|---|
| CN 變異區 | 高（66.09% 可定樹） | 低（AF 差異可能是 multiplicity） |
| CN 中性區 | 低（35.59%） | 高（AF 選擇無拷貝數混淆） |

### 結論：夠做什麼、不夠做什麼

| 用途 | 判定 | 理由 |
|---|---|---|
| 證明方法在乾淨地上可運作（proof-of-concept） | ✅ **夠** | 63 個可定樹 unit、72 條無爭議 deep edge |
| 作為 CN-free 的驗證錨點 / 陰性對照 | ✅ **夠** | 這 72 條是全樣本唯一不需拷貝數假設的譜系關係 |
| 描述這個腫瘤的演化史 | ❌ **不夠** | 能看多步譜系的只有 14 個 unit、29 條 deep edge，且集中在 5 條染色體（chr2 12／chr21 9／chr16 8 就佔 29/43） |
| 統計推論（偵測中小效應） | ❌ **不夠** | 檢定力只到 OR 1.662 |
| 取代 CN 變異區重做全基因組分析 | ❌ **不可能** | 98.44% 的 deep edge 在 CN 變異區 |

**實務結論：不能靠「扣掉 CN 區」來做研究。** 正確做法是**保留全部、標註 CN 狀態、分層報告** ——
CN 變異區照樣輸出候選樹與共識骨幹（形狀與結構層都是 CN-robust 的），
但其 read-AF 挑出的「唯一樹」一律停在 `RAW_AF_UNIQUE_REPRESENTATIVE`；
CN 中性的那 72 條 deep edge 則作為可對外宣稱的乾淨錨點。

---

## 7D. 跨樣本累積乾淨 deep edge（H1437／H2009／HCC1954）

### 先驗證「用 SAVANA 判乾淨區」準不準

擴到無真值樣本前，必須先在 HCC1395 上量測這件事本身。判準是
`total CN ≈ 2（±0.5）且段 mean BAF < 0.75`（1+1 才會讓單 copy haplotype 的 clonal 突變讀出 1/1）：

| 變體 | 判為 neutral 的 bin | Precision | Recall |
|---|---:|---:|---:|
| SAVANA 發布 CN，無 BAF gate | 9,174 | 10.83% | 34.23% |
| SAVANA 發布 CN + BAF gate | 349 | 25.21% | 3.03% |
| 校正 CN（purity 1.0／ploidy 2.95），無 BAF gate | 10,414 | 4.9% | 21.78% |
| **校正 CN + BAF gate** | 206 | **90.29%** | 7.94% |

**只有「校正 CN + BAF gate」可用** —— precision 90.29%，但 recall 僅 7.94%。
對建立錨點而言這正是要的性質：**寧可少而乾淨**。

### 三樣本結果（每個 unit 的所有 active 位點都須落在乾淨區）

| 樣本 | ranked | 判為乾淨 | 乾淨 k≥2 | 可定樹 | **deep edge** | k≥3 deep | 嚴格版 deep | 染色體 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| H1437 | 13,740 | 186 | 105 | 8.6% | **12** | 5 | 12 | 3 |
| H2009 | 23,128 | 8,700 | 6,595 | 5.9% | **496** | 221 | 483 | 17 |
| HCC1954 | 5,647 | 743 | 379 | 64.4% | **222** | 55 | 106 | 11 |

（嚴格版另要求 minor-allele CN ≥ 0.5；HCC1954 因 purity 僅 0.66，LOH gate 較弱，嚴格版砍半到 106。）

### 累積結果

| | k≥2 unit | k≥2 deep edge | k≥3 deep edge |
|---|---:|---:|---:|
| HCC1395（SEQC2 真值） | 177 | 72 | 29 |
| + H1437 / H2009 / HCC1954（SAVANA + gate） | 7,256 | **802** | **310** |

**乾淨 deep edge 從 72 增至 802（11 倍）**，HCC1395 只佔 8.98%。

### 🔑 跨樣本獨立重現了 §6 的機制

| 樣本 | 乾淨區可定樹率 | CN 變異區可定樹率 |
|---|---:|---:|
| HCC1395 | 35.59% | 66.09% |
| H1437 | 8.6% | 23.7% |
| H2009 | 5.9% | 27.1% |
| HCC1954 | 64.4% | 84.9% |

**4/4 樣本都是「乾淨區反而更定不出樹」。** 這在四個獨立樣本上重現了
「AF 值相同 → 分數相同 → 數學上無法破 tie」的機制，比單樣本結論強得多。

⚠ **這 730 條非 HCC1395 的 deep edge，其可信度是「期望」而非「量測」**：
三樣本用 SAVANA 發布的 purity/ploidy（通過 BAF 自洽性檢驗，違反率 5.62–9.81%），
但無外部 CN 真值，故 90.29% 的 precision 是從 HCC1395 遷移過來的推論。

---

## 7E. 與 SEQC2 外部 clonality 分析對照 —— 尺度不同，不能直接比

SEQC2 主論文（Fang et al., *Nat Biotechnol* 2021）對 HCC1395 做過 clonality 分析：

| SEQC2 的做法 | 我們的做法 |
|---|---|
| **S1 = MRCA，S2–S10 = 9 個 subclone**（S2 的 cancer cell fraction 60%） | local per-PS×HP block 的 mutation-state 樹 |
| 10x Single Cell CNV：**1,270 個 tumor 細胞** + 638 個 normal 細胞 | 單一 ONT bulk，無單細胞 |
| SNV clonality = VAF **經 local copy number 校正** | **做不到**（缺 allele-specific CN，見 §10 L3） |
| 21 個 WGS replicates 的 hierarchical clustering | 單一資料集 |
| subclone 主要由**染色體 gain/loss 事件**定義 | 由 **read-level sSNV 共現**定義 |

**兩者是不同尺度的物件，不是同一棵樹的兩個版本：**
SEQC2 的 S1–S10 是**全基因組的細胞群演化樹**；我們的是**局部突變狀態樹**，
而且 2026-08-01 規格 §7 明文禁止跨 block 合併（無 bridging read 時不可宣稱同一 cellular lineage）。

### 可以對照的兩點（都吻合）

1. **「大部分突變是 clonal」** —— SEQC2 稱大部分 driver 點突變位於 MRCA。
   我們在 CN-neutral 區（depth≥20，n=514）觀察到 **48.05% 的位點 read-AF ≈ 1.0**，
   形成清楚的 clonal 主峰。方向一致。
2. **「高度異質、大量 CNA 在 sub-population」** —— SEQC2 由 1,270 個單細胞確認。
   我們獨立觀察到 **95.44% 的 unit 坐在 CN-altered 區**。方向一致，
   且這正好解釋了為什麼我們的乾淨區如此稀少 —— 不是方法問題，是這顆腫瘤的性質。

### 不能對照、且我們做不到的部分

CN-neutral 的 267 個 subclonal 位點（AF < 0.95）其 AF 呈**連續分布**（中位 0.4583），
落在 0.55–0.65（對照 SEQC2 S2 的 CCF 60%）的僅 40 個（14.98%），
**沒有對應 S1–S10 離散 CCF 的清晰峰**。

原因不必歸咎於方法：n 太小、位點分散於不同 local block（可能屬於不同 subclone）、
且 depth 20–50 時 AF 的抽樣誤差已達 ±0.1。
**結論是「解析度不足以重建 S1–S10」，不是「與 SEQC2 矛盾」。**

---

## 7F. 結構層能不能反過來「限制」VAF-based 的 subclone 推論？——**可以，而且在 CN 區最有用**

### 兩類方法的資訊落差

VAF-based 重建（PyClone／SciClone／DPClust／CITUP 等）只有**邊際 VAF**。
它看不到哪些突變在同一條分子上，只能從頻率反推共存關係（pigeonhole／sum rule）：
`f_A ≥ f_B` 就允許 A 為 B 的祖先，`f` 相近就歸為同一 clone。

長讀長直接提供那個缺失的觀測量。對一組共享 read 的突變對，觀測到的 haplotype pattern
**不需要任何頻率論證**就能定案：

| 觀測到的 pattern | 結論 |
|---|---|
| `{RR, AR, AA}` | 巢狀，A 先於 B |
| `{RR, RA, AA}` | 巢狀，B 先於 A |
| `{RR, AR, RA}` | **互斥 —— 不同分支，不存在祖先關係** |
| 四種都出現 | 違反完美系統發生（需 recurrence 或 loss） |

### 量測：VAF-only 會在哪裡判錯

取 HCC1395 中 `k=2` 的 unit（每個恰好隔離出一組成對關係，兩位點深度皆 ≥10），共 3,926 對：

| 關係 | 比例 |
|---|---:|
| 巢狀（可定先後） | 86.19% |
| **互斥（不同分支）** | **13.81%** |
| 四配子違反 | 0（ranked unit 已排除 recurrence screen 的 45 個） |

**關鍵數字**：那 542 對互斥關係中，**26.01%（141 對）的兩個突變 AF 差距 ≥ 0.20**。
VAF-only 方法看到這種差距會判為祖先→後代，**但分子證據說它們在不同分支、不可能有祖先關係**。

反之，巢狀對違反 AF 單調性的只有 **0.09%（3/3,384）** —— 在真正巢狀的情況下 VAF 大致可靠。
**問題不在 VAF 排錯順序，而在 VAF 無法察覺「根本不該排序」。**

### 🔑 這個落差正好集中在 CN 變異區

| | 互斥對中 AF 差距 ≥0.2 的比例 | 互斥對的 AF 差距中位 |
|---|---:|---:|
| CN-neutral | **7.14%**（1/14） | **0.0** |
| CN-altered | **26.52%**（140/528） | 0.0965 |

在拷貝數乾淨的地方，互斥的兩個突變 AF **中位差距為 0** —— VAF 方法看到相同的頻率，
不會硬排順序，判斷不會錯得離譜。
在 CN 變異區，multiplicity 給了互斥突變一個**假的頻率梯度**，
VAF 方法就會把它讀成祖先關係。

**於是出現一個漂亮的互補性：CN 破壞 AF 層的地方，正是 linkage 約束最有價值的地方。**

### 我們能提供多少硬約束（全 7 樣本）

| 樣本 | k≥2 unit | 順序約束（deep edge） | 互斥宣告（分支 unit） |
|---|---:|---:|---:|
| HCC1395 | 5,984 | 7,554 | 1,167 |
| H1437 | 9,646 | 15,902 | 1,757 |
| H2009 | 18,612 | 37,712 | 5,138 |
| HCC1954 | 3,003 | 3,177 | 843 |
| HCC1395_DORADO | 3,395 | 3,605 | 603 |
| COLO829 | 6,158 | 7,288 | 1,013 |
| HCC1937 | 3,022 | 3,417 | 756 |
| **合計** | **49,820** | **78,655** | **11,277** |

### 定位結論（⚠ 已於 §7G 修正，見該節）

我們的結構層**不是**要取代 VAF-based 重建，而是提供它拿不到的兩類硬約束：

1. **排除約束**（11,277 個 unit）：「這兩個突變不在同一條分子上」——
   VAF 無法產生此判斷。**但 §7G 證明此約束只在 CN-neutral 區能推論到細胞層**，
   在 CN 區僅為分子層陳述。
2. **順序約束**（78,655 條）：「A 必定先於 B」——
   由分子直接觀測，不依賴頻率、不需 CN 校正。**§7G 確認此約束 CN-free by construction。**

**誠實的邊界**：這些約束都是 **local**（同一 PS×HP block 內、有共享 read 的位點對），
覆蓋率遠低於全基因組 VAF 方法。它們能**約束與否證**一棵候選演化樹，
但**不能單獨建構**全樣本的 clone 樹（規格 §7 禁止跨 block 合併）。
正確用法是：**VAF 方法產生候選 → 我們的約束篩掉不相容者**。

---

## 7G. 約束方法是否成立？能否排除 CN 影響？——**兩類約束地位不同，須分開回答**

§7F 把兩類約束並列，但它們的效力**不對等**。逐一稽核。

### 約束一：順序（`{RR, AR, AA}` → A 先於 B）—— ✅ **成立，且 CN-free by construction**

**論證**：B 只出現在已帶 A 的分子上。在 infinite-sites 模型下突變沿分子譜系只增不減，
因此 A+B 的分子必然衍生自 A 的分子。
**這個推論從頭到尾沒有引用 copy number、ploidy 或 purity** —— 它是關於**分子繼承**的陳述，
不是關於頻率的陳述。CN 改變的是「同一條 haplotype 有幾份」，不改變「B 出現在 A 的背景上」。

| | 條數 |
|---|---:|
| 全部順序約束（HCC1395） | 7,554 |
| 其中 CN-neutral | 214 |
| **其中 CN-altered** | **7,340（全部有效）** |

⚠ **層級限定**：這是**分子譜系**的順序，不是**細胞譜系**的順序。可以說
「B 這個突變發生在已帶 A 的 DNA 分子上」，不能直接說「帶 B 的細胞是帶 A 細胞的後代」。

⚠ **真正的殘餘威脅是抽樣，不是 CN**：宣稱「沒看到 `RA`」可能只是漏抽。以連結深度量化
（若某 pattern 真實比例為 p，在 D 條分子中漏掉的機率是 (1−p)^D）：

| 連結深度 | unit 數 | 佔比 | 能以 95% 信心排除的 pattern 比例 |
|---|---:|---:|---|
| 3–9 | 468 | 7.82% | 僅能排除 >28.31–63.16%（**很弱**） |
| 10–19 | 902 | 15.07% | >14.59–25.89% |
| 20–49 | 3,314 | 55.38% | >5.93–13.91% |
| ≥50 | 1,300 | 21.72% | **>1.49–5.82%（強）** |

中位連結深度 33，中位可排除比例 **8.68%**。
→ **順序約束對「低比例的第三種 pattern」不敏感**；depth<20 的 22.89% unit 應標為低信心。

### 約束二：互斥（`{RR, AR, RA}` → 不同分支）—— ⚠ **只在 CN-neutral 區能推到細胞層**

**在 CN-neutral（1+1）**：每個細胞對該 haplotype 只貢獻一份分子。
「沒有分子同時帶 A 和 B」⇒「沒有細胞同時帶 A 和 B」⇒ **真的是不同細胞分支** ✅

**在 CN gain**：該 haplotype 每個細胞有 c > 1 份。
**A 可以在 copy 1、B 在 copy 2，屬於同一個細胞** —— 沒有任何單一分子同時帶 A+B，
但這個細胞兩個突變都有。
**分子層的觀測完全不變，細胞層的推論卻整個崩掉。** ❌

| | 數量 | 佔比 |
|---|---:|---:|
| 互斥 unit（HCC1395） | 1,167 | 100% |
| **能推論到細胞層（CN-neutral）** | **25** | **2.14%**（CI 1.46–3.14） |
| 只能停在分子層（CN-altered） | 1,142 | 97.86% |

### 🔴 對 §7F 的修正

§7F 把「排除約束」講成兩者中較有價值的一項，**那個框架是錯的**。正確的是相反：

| 約束 | 能否排除 CN 影響 | 全基因組可用性 |
|---|---|---|
| **順序** | ✅ **可以，by construction** | 7,554 條全部可用（分子譜系層） |
| **互斥** | ❌ **不能** | 細胞層僅 25 條（2.14%）可用 |

**原因很簡單**：順序約束問的是「B 是否出現在 A 的背景上」，這與有幾份 copy 無關；
互斥約束問的是「有沒有細胞同時帶兩者」，而那正是 copy number 直接介入的地方。

### 修正後的對外措辭

可以寫：

> Long-read linkage yields 7,554 ordering constraints of the form "mutation B was
> acquired on a molecule already carrying A". This is a statement about molecular
> descent under infinite sites and is independent of copy number, ploidy and purity.
> Mutual-exclusivity observations are additionally reported, but are promoted to a
> cellular-branch claim only on copy-neutral segments (25 units here), since under
> copy gain two mutations on different copies of the same haplotype in the same cell
> produce no doubly-mutant molecule.

不可以寫：

- 「互斥的兩個突變屬於不同 subclone」（在 CN 區不成立）
- 「A 的細胞是 B 的細胞的祖先」（順序約束只到分子層）
- 「未觀測到某 pattern 即該 pattern 不存在」（depth<20 時只能排除 >15–26% 的比例）

---

## 8. 問題六：SAVANA 作為流程驗證 —— 訊號好，校準錯

HCC1395 是唯一有外部 CN 真值的樣本，因此是唯一能驗證 SAVANA 流程本身的地方。

| | purity | ploidy | 狀態一致率 | 整數 CN 吻合 | 平均絕對誤差 |
|---|---:|---:|---:|---:|---:|
| SAVANA 發布值 | 0.76 | 1.83 | **28.55%** | **3.64%** | 1.1206 |
| grid search 最佳 | **1.0** | **2.95** | **89.91%** | **83.3%** | 0.4834 |
| 整數吻合最佳 | 1.0 | 2.8 | 89.78% | 83.64% | 0.4171 |

- SEQC2 自身隱含的常染色體 ploidy（bp 加權）＝ **2.9104**，與 grid 最佳的 2.95 幾乎重合。
- 逐染色體偏差幾乎一致為 **−0.8 ～ −1.2 個 copy**（僅 chr16 +1.6387、chr6 +0.401 例外）。
- 段層回歸：全部 bin 的 R² 僅 0.414，但**非中性區 R² 達 0.8693**（slope 0.8228）。
  全體 R² 偏低是 CN=2 常數值稀釋造成的統計 artifact，**非 segmentation 品質不良**。

**裁決：SAVANA 的 log2r 訊號與 segmentation 是可用的；錯的是 purity/ploidy 擬合。**
這是校準失敗（原則上可矯正），不是訊號失敗（不可救）。
HCC1395 真值為純腫瘤細胞株 purity≈1.0、ploidy≈2.9104，SAVANA 挑了退化解。

---

## 9. 問題七：回頭裁決其他樣本

用一個**不需外部真值**的內部檢驗：對每段以其自身 CN 計算 BAF 理論上限
`BAF_max = (ρn + (1−ρ)) / (ρn + 2(1−ρ))`（完全 LOH 是最極端情形），統計超出上限的段比例。

**先用 HCC1395 校準此檢驗**：

| HCC1395 | 違反率 |
|---|---:|
| 在發布的 purity 0.76 下 | **54.65%** |
| 在校正後 purity 1.0 下 | **0.00%** |

→ 這個純內部檢驗**確實能抓到** mis-fit，且其判定被 SEQC2 外部真值背書。套用到其他樣本：

| 樣本 | purity | ploidy | 可用段數 | 違反率 | 裁決 |
|---|---:|---:|---:|---:|---|
| HCC1395 | 0.76 | 1.83 | 666 | **54.65%** | mis-fit（已由 SEQC2 證實） |
| HCC1937 | 0.62 | 3.13 | 613 | **50.41%** | mis-fit，維持 cn=unknown |
| H2009 | 0.95 | 2.22 | 693 | 9.81% | 通過 |
| H1437 | 0.95 | 2.9 | 907 | 7.17% | 通過 |
| HCC1954 | 0.66 | 3.21 | 818 | 5.62% | 通過（但最差超額 0.414，有少數極端段） |

**這次獨立檢驗支持而非推翻既有判定**：H1437/H2009/HCC1954 可用、HCC1937 不可用。
COLO829 有三個 fit 版本（原始／refit／overrule），本輪未納入比較，屬缺口。

---

## 10. 可以矯正嗎？—— 方案設計（本輪不執行重算）

對照 2026-08-01 規格 §8.4 已列出的 CN/LOH 欄位，缺口與可行性：

| 矯正層級 | 需要什麼 | HCC1395 | 其他樣本 | 評估 |
|---|---|---|---|---|
| **L1 CN 註釋層** | total CN + LOH | ✅ SEQC2 整數 CN 已可用（本輪已建） | SAVANA 校準後可得 | **現在就可做** |
| **L2 期望 AF 帶** | total CN + purity + LOH | ✅ purity=1.0、CN 已知 | H1437/H2009/HCC1954 可近似 | **可做**，但需假設 allele 分配 |
| **L3 真 multiplicity 校正** | **allele-specific CN**（major/minor 分別） | ❌ SEQC2 只給 total CN + LOH 標籤 | ❌ | **做不到**，這是真正的阻斷點 |
| **L4 HP↔allele 對應** | 哪個 germline HP 對應 major allele | ❌ 無此對應 | ❌ | **做不到** |

**建議的可行做法（不需 L3/L4）＝ fail-closed gate，而非數值校正：**

1. 為每個 unit 標註 `cn_loh_status`（規格 §8.4 已定義欄位名），值域：
   `CN_NEUTRAL_NO_LOH` / `CN_ALTERED` / `CN_UNKNOWN`。
2. **只有 `CN_NEUTRAL_NO_LOH` 的 unit，其 read-AF 排序才可進入生物學解讀**；
   其餘一律停在規格 §5.2 已規定的 `RAW_AF_UNIQUE_REPRESENTATIVE`，不得升格。
3. 對 `CN_ALTERED` 的 unit，額外標 `af_grid_consistent`（AF 是否落在 m/c 格點容忍區內）
   作為警示旗標 —— 落在格點上者，其 AF 差異優先以 multiplicity 解釋。
4. 所有對外數字（單一拓撲率等）**同時報 CN-neutral 子集版本**。

這是把「矯正」實作成**適用性閘門**而非數值修正。理由：L3 缺 allele-specific CN，
任何數值校正都必須假設 allele 分配，會把假設偷渡成結果 —— 那正是本專案要避免的模式。

**代價要說清楚**：套用此 gate 後，HCC1395 可用於 AF 生物學解讀的 unit 只剩 **439 個（4.56%）**，
且其中 89.46% 集中在 4 條染色體。這是誠實的可用範圍，不是失敗。

---

## 11. 之前的 subclone 結論還站得住嗎？

| 既有結論 | 裁決 | 理由 |
|---|---|---|
| sSNV 單分子共現是非循環骨幹 | ✅ **維持** | 共現是直接觀測，結構唯一率在 neutral 反而更高（56.41% vs 48.82%），無 CN 抬高跡象 |
| 拓撲軸 CN-robust | ✅ **維持** | 同上，本輪以 canonical 分母重新驗證 |
| 甲基化為 bounded-auxiliary、不進 likelihood | ✅ **不受影響** | 本輪未觸及；CN 問題出在 read-AF 層 |
| read-AF 僅作 model-conditional ordering，不可作生物確認 | ✅ **強化** | 本輪給出量化：CN 可歸因超額 14.81–23.55% |
| 「88.26% 單一拓撲」「55.1% 唯一最佳樹」 | ⚠ **必須加註** | 這些數字是 CN-conditional；未標註即等同宣稱生物學拓撲盛行率，而規格 §10 已明文禁止 |
| VAF/頻率軸嚴重 CN-limited | ✅ **維持並量化** | 95.44% unit 坐在 CN-altered 區 |
| 「純 parsimony 單一拓撲率下界 64.89%」 | ⚠ **仍無出處** | 本輪未能補上，維持既有警告 |

**總評：不需要撤回任何已發表結論，但需要為 AF 層的數字補上 CN-conditional 標註與可用範圍聲明。**
本輪的發現與 2026-08-01 規格的方向一致 —— 該規格已寫明 allele-specific CN/LOH 未整合、
AF-unique 只能停在 `RAW_AF_UNIQUE_REPRESENTATIVE`。本輪把那個「未整合」從定性警告
變成了帶區間的量化，並補上規格 §9 第 6 項所要求的 independent CN annotation。

---

## 12. 限制（不因本輪完成而消失）

1. **單一樣本**。HCC1395 是唯一有 CN 真值者，且高度非整倍體，CN-neutral 對照組僅 4.56%。
   本輪結論不可外推為「所有樣本的 CN 影響都是這個量級」。
2. **CN-neutral 組染色體高度集中**（89.46% 在 4 條染色體），按染色體分層後 OR 由 3.53 降至 2.16。
   效應存活但量級減半，§7 的超額估計應視為上緣。
3. **無 allele-specific CN**。SEQC2 只提供 total CN 與 LOH 標籤，
   §6.2 的 m/c 格點分析依賴「LOH 時 c = total CN、非 LOH 時均分」的假設，是 model-conditional。
4. **機制在校正後與「完全經由算術可解性中介」一致**（非退化子集 OR 1.2765，p = 0.4459），
   但該子集 n_neutral 僅 49，尚不能排除達 OR ≈ 2.6 的部分中介。
   初版曾宣稱存在顯著殘餘（OR 2.82），該宣稱源自 §6.0 的位點集錯誤，已撤回。
5. **未重算**。依任務界定，本輪不重跑 CCF、樹枚舉或 subclone 判定，
   因此「矯正後結論如何改變」只有設計與估計，沒有實跑結果。
6. **COLO829 未納入** SAVANA 自洽性比較（有三個 fit 版本，需先決定用哪個）。
7. **SEQC2 CN 本身**是 NGS benchmark 推導的共識，非細胞層真值，且含 half-integer（如 3.5）
   代表 subclonal CN —— 本輪一律當作段層代表值處理。

---

## 13. 產物

| 檔案 | 內容 |
|---|---|
| `data/cn_annotation_summary.json` | unit / site 層 CN 分類、逐染色體分布、AF×CN |
| `data/hcc1395_unit_cn_annotation.jsonl` | 9,624 個 unit 的逐項 CN 註釋 |
| `data/hcc1395_site_cn_annotation.jsonl` | 20,060 個 sSNV 的逐點 CN 與 read-AF |
| `data/savana_vs_seqc2_adjudication.json` | 段層一致性、回歸、逐染色體偏差 |
| `data/savana_refit_grid.json` | purity×ploidy grid search（1,586 點） |
| `data/purity_selfconsistency_audit.json` | 五樣本 BAF 違反率 + HCC1395 校準 |
| `data/resolution_vs_cn_stratified.json` | 兩層唯一性拆解 + CMH |
| `data/mechanism_v2_distinctness.json` | 機制檢驗（相異性） |
| `data/robustness_checks.json` | 染色體混淆檢驗 + 套套邏輯檢驗 |
| `data/consolidated_findings.json` | 報告注入用彙總 |

### 具體案例（見樹）

**`chr8|PS=126637291|HP=1`** — chr8:127,307,102–127,323,260（16 kb），4 個 sSNV，
SEQC2 **total CN 7.31、LOH 涵蓋率 1.0**，12 棵候選樹經 read-AF 選出唯一贏家，
標為 `UNIQUE_TREE`。在一條 haplotype 帶約 7 個 copy 的區段上，
AF 只能取 m/7，這些差異幾乎確定來自 multiplicity 而非發生順序。
**這是 §7 那 1,379 個超額唯一樹的典型形態。**

對照組 **`chr2|PS=196625788|HP=2`** — CN=2、LOH=0，單一 sSNV，結構本身即唯一樹（read-AF 未參與）。
