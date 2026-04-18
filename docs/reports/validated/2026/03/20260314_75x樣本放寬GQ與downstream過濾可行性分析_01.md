<!--
建立時間: 2026-03-14 22:40
目標: 分析在 paired HCC1395 樣本中，是否值得放寬 GQ，並依賴 LongPhase-S / InterSubMod / 甲基特徵 / LOH-like / CNV proxy 在後段壓回 FP
處理範圍:
  - HCC1395 5kHz paired（主樣本）
  - HCC1395_DORADO paired（交叉檢查）
關聯檔案:
  - /big8_disk/liaoyoyo2001/InterSubMod/scripts/analysis/analyze_relaxed_gq_downstream_strategy.py
  - /big8_disk/liaoyoyo2001/InterSubMod_runs/output/big8_disk_output/research_rounds/20260314_relaxed_gq_downstream_strategy_audit/analysis_summary.md
  - /big8_disk/liaoyoyo2001/InterSubMod_runs/output/big8_disk_output/research_rounds/20260314_relaxed_gq_downstream_strategy_audit/rule_strategy_summary.tsv
  - /big8_disk/liaoyoyo2001/InterSubMod_runs/output/big8_disk_output/research_rounds/20260314_relaxed_gq_downstream_strategy_audit/rule_distribution_summary.tsv
  - /big8_disk/liaoyoyo2001/knowledge/05_tools/longphase_s.md
  - /big8_disk/liaoyoyo2001/knowledge/06_workflows/phasing_workflow.md
  - /big8_disk/liaoyoyo2001/knowledge/06_workflows/benchmark_workflow.md
-->


<!-- PATH WARNING (2026-04-11 驗證):
本報告中部分絕對路徑已失效。類別: Knowledge 路徑大小寫差異
詳細路徑清單見 docs/data_specs/20260411_path_inventory.tsv
報告結論仍有效，但引用的外部路徑可能需要更新。
-->

# 75x 級 paired 樣本放寬 GQ 與 downstream 過濾可行性分析

## 一句結論

在目前證據下，**不建議把 `HCC1395 5kHz paired` 的 primary GQ gate 放寬到 `gq>=12` 以下，更不建議放到 `gq>=10` 再期待後段 LongPhase-S / InterSubMod / 甲基 / 弱 LOH/CNV proxy 幫忙把 FP 壓回來**。  
主因不是 downstream 完全沒訊號，而是它們的淨效果**不足以抵消 relaxed gate 帶進來的額外 FP**。

更精確地說：

1. `5kHz paired` 最接近可接受的 relaxed 策略是 `gq>=12 + Quality_Score>=60`，但它相對 `LongPhase-S` 幾乎只是回到 baseline，**沒有超過**現行 `gq>=15 / 18 / 20`。
2. `gq>=10` 在 `5kHz paired` 會額外帶入 `616` 個 FP；目前所有已驗證的 downstream 支持、artifact veto、弱 LOH/CNV proxy 都**救不回來**。
3. 這個 75x 級樣本裡，`DP` 本身不是主分離器；真正把 `TP` 與 `FP` 拉開的，是 `GQ/QUAL`、`ref/alt 結構`、以及 `LowQual` error 狀態。
4. `haplotype / 弱 LOH/CNV proxy` 訊號在主樣本裡不是安全 veto。某些高 `GQ` 的剩餘 `FP` 反而**更常**帶有 `has_h_flag` 或 `Potential_LOH` / 極端 `AF` proxy。

## 背景

根據 [Knowledge/05_tools/longphase_s.md](/big8_disk/liaoyoyo2001/knowledge/05_tools/longphase_s.md) 與 [Knowledge/06_workflows/phasing_workflow.md](/big8_disk/liaoyoyo2001/knowledge/06_workflows/phasing_workflow.md)，paired 流程的高層分工是：

1. `ClairS` 產生 paired somatic call set。
2. `LongPhase-S` 進行 somatic haplotagging 與 caller-side 收斂。
3. `InterSubMod` 使用 tagged BAM 與甲基特徵做後段判讀。

根據 [Knowledge/06_workflows/benchmark_workflow.md](/big8_disk/liaoyoyo2001/knowledge/06_workflows/benchmark_workflow.md)，本輪比較口徑固定對齊 `SEQC2` truth，並同時列出 `TP / FP / FN / precision / recall / F1`。

這一題真正要回答的是：

1. 若前段把 `GQ` 放寬，能否靠後段 `InterSubMod` 特徵把新增 `FP` 壓回去？
2. 這件事在 read 數、倍體 / LOH-like proxy、haplotype 狀況上，是否有足夠機制證據支持？

## 方法

本輪新增可重跑腳本：

- [analyze_relaxed_gq_downstream_strategy.py](/big8_disk/liaoyoyo2001/InterSubMod/scripts/analysis/analyze_relaxed_gq_downstream_strategy.py)

主輸出目錄：

- [/big8_disk/liaoyoyo2001/InterSubMod_runs/output/big8_disk_output/research_rounds/20260314_relaxed_gq_downstream_strategy_audit](/big8_disk/liaoyoyo2001/InterSubMod_runs/output/big8_disk_output/research_rounds/20260314_relaxed_gq_downstream_strategy_audit)

關鍵輸出：

1. [analysis_summary.md](/big8_disk/liaoyoyo2001/InterSubMod_runs/output/big8_disk_output/research_rounds/20260314_relaxed_gq_downstream_strategy_audit/analysis_summary.md)
2. [rule_strategy_summary.tsv](/big8_disk/liaoyoyo2001/InterSubMod_runs/output/big8_disk_output/research_rounds/20260314_relaxed_gq_downstream_strategy_audit/rule_strategy_summary.tsv)
3. [feature_availability.tsv](/big8_disk/liaoyoyo2001/InterSubMod_runs/output/big8_disk_output/research_rounds/20260314_relaxed_gq_downstream_strategy_audit/feature_availability.tsv)
4. [rule_distribution_summary.tsv](/big8_disk/liaoyoyo2001/InterSubMod_runs/output/big8_disk_output/research_rounds/20260314_relaxed_gq_downstream_strategy_audit/rule_distribution_summary.tsv)
5. [rule_category_summary.tsv](/big8_disk/liaoyoyo2001/InterSubMod_runs/output/big8_disk_output/research_rounds/20260314_relaxed_gq_downstream_strategy_audit/rule_category_summary.tsv)

分析邏輯分兩層：

1. 規則層：比較 `gq>=10/12/15/18/20` 與幾種 downstream 支持/過濾組合的 `delta F1`。
2. 機制層：對被選中的 `TP / FP` 比較 `DP / AF / AD_ref / AD_alt / PairwiseMedianDist / Quality_Score / HP / LOH-like proxy` 分布。

## 主要結果

### 1. 主樣本 `HCC1395 5kHz paired`：放寬到 `gq>=10` 明確不划算

基線 benchmark：

| Method | TP | FP | FN | F1 |
| --- | ---: | ---: | ---: | ---: |
| `LongPhase-S` | `29754` | `627` | `9693` | `0.8522` |
| `InterSubMod` | `29752` | `544` | `9695` | `0.8532` |

主要策略比較：

| 策略 | rescued TP | reintroduced FP | delta F1 vs LongPhase-S | delta F1 vs InterSubMod |
| --- | ---: | ---: | ---: | ---: |
| `gq>=10` | `183` | `616` | `-0.004451` | `-0.004486` |
| `gq>=12` | `142` | `274` | `-0.000996` | `-0.001021` |
| `gq>=12 + Quality_Score>=60` | `68` | `91` | `+0.000015` | `-0.000006` |
| `gq>=15` | `106` | `75` | `+0.000833` | `+0.000813` |
| `gq>=18` | `73` | `21` | `+0.000951` | `+0.000932` |
| `gq>=20` | `59` | `8` | `+0.000880` | `+0.000861` |

解讀：

1. `gq>=10` 帶入的 FP 太多，downstream 目前無法回收。
2. `gq>=12` 雖然比 `gq>=10` 好很多，但仍然輸給 `gq>=15`。
3. `gq>=12 + Quality_Score>=60` 是最接近可用的 relaxed 組合，但本質上只是把集合重新縮回一個更嚴格的高品質子集；它**沒有帶來比 `gq>=15` 更好的最終 F1**。
4. 如果目標是這個樣本的最高 F1，`gq>=18` 與 `gq>=20` 其實比 relaxed + downstream 更乾淨。

### 2. read 數量不是主分離器，75x 樣本中的 FP 也很有 read 支持

`5kHz paired` 幾個關鍵策略下的 `DP / AF` 中位數：

| 策略 | subset | DP median | AF median |
| --- | --- | ---: | ---: |
| `gq>=10` | TP | `73` | `0.386` |
| `gq>=10` | FP | `71` | `0.130` |
| `gq>=12` | TP | `63` | `0.495` |
| `gq>=12` | FP | `73` | `0.165` |
| `gq>=15` | TP | `56` | `0.544` |
| `gq>=15` | FP | `88` | `0.206` |
| `gq>=20` | TP | `50` | `0.622` |
| `gq>=20` | FP | `96.5` | `0.264` |

這裡最重要的觀察是：

1. `FP` 並不因為 read 少就被排掉。很多 `FP` 的 `DP` 其實不低，甚至比 `TP` 還高。
2. 真正差異在 `AF` 與 `ref/alt` 結構：
   - `TP` 更常呈現較高 `AF` 與明顯 ALT 主導
   - 但 relaxed gate 放進來的 `FP` 也可能有大量 reads，只是仍帶著較重 ref burden 或較差 caller 品質
3. 這代表在 75x 級樣本中，單純「coverage 很高」不等於可以安全放寬 `GQ`。

### 3. 真正的主分離器是 caller-side 結構，不是後段甲基訊號

`5kHz paired` 在 `gq>=15` 下：

| 特徵 | rescued TP median | rescued FP median |
| --- | ---: | ---: |
| `GQ` | `20` | `16` |
| `QUAL` | `18.47` | `0.00` |
| `DP` | `56` | `88` |
| `AF` | `0.5436` | `0.2063` |
| `AD_ref` | `0` | `63` |
| `AD_alt` | `30` | `16` |

這組差異很明確地說明：

1. `TP` 比較像 caller 邊界附近但仍合理的真變異。
2. `FP` 則常是 read 很多、但仍保留明顯 ref 支持，或直接還是 `LowQual` 的位點。

也就是說，paired 主樣本裡最穩定的分離力還是在 caller 這一層，而不是把 gate 放鬆後再期待後段甲基特徵幫忙補救。

### 4. `Quality_Score` 與 `PairwiseMedianDist` 能篩掉一部分 FP，但不夠

`gq>=12 + Quality_Score>=60` 的性質很有代表性：

| subset | n | DP median | AF median | PASS fraction | Strong/Subclone fraction |
| --- | ---: | ---: | ---: | ---: | ---: |
| TP | `68` | `63` | `0.544` | `100%` | `30.88%` |
| FP | `91` | `68` | `0.170` | `100%` | `51.65%` |

解讀：

1. 這組規則會把 `LowQual` 幾乎全清掉，所以表面上看起來變乾淨。
2. 但它同時也保留了不少 `PASS` 且 `Quality_Score` 很高的 `FP`。
3. 更關鍵的是，這些 `FP` 並不一定落在「很弱的甲基類別」，反而有相當比例進到 `Strong/Subclone`。

因此這裡不能把 `Quality_Score` 解讀成一個足以接手 primary gating 的特徵；它比較像是後段 ranking/support，而不是彌補 relaxed GQ 的主防線。

### 5. haplotype 與弱 LOH/CNV proxy 在主樣本中不是安全 veto

`5kHz paired` 沒有這輪完整的 paired candidate-specific `Potential_LOH / Coverage_Category`，所以本輪只能用弱 proxy。這裡的 proxy 只表示 haplotype-imbalance / LOH-CNV 線索，**不是正式 LOH call**：

1. `has_h_flag_variant`
2. `extreme_AF_proxy`：`AF>=0.9 and AD_ref=0`

把兩者合成 `merged_loh_proxy` 後，結果並不支持它作為 relaxed gate 的補救條件。

例子：

| 策略 | subset | has_h_flag fraction | merged_loh_proxy fraction |
| --- | --- | ---: | ---: |
| `gq>=15` | TP | `44.34%` | `57.55%` |
| `gq>=15` | FP | `57.33%` | `60.00%` |
| `gq>=20` | TP | `61.02%` | `79.66%` |
| `gq>=20` | FP | `100.00%` | `100.00%` |

這代表：

1. 弱 LOH/CNV proxy 樣態不是 TP 專屬。
2. 在高 `GQ` 殘留 `FP` 裡，這種訊號反而可能更強。
3. 也就是說，若把 relaxed gate 的候選交給弱 LOH/CNV proxy 決定去留，風險很高。

### 6. `DORADO paired` 交叉檢查：downstream 可以幫一點，但仍不值得用來合理化 relaxed gate

`DORADO paired` 最好的 relaxed downstream 組合是：

- `gq>=10 + Quality_Score>=60 + hp_assign_rate>=0.99`
- `50 TP / 21 FP`
- `delta F1 vs LongPhase-S = +0.000536`

但同一個樣本更簡單的 caller-first 規則仍更好：

- `gq>=18`: `67 TP / 24 FP`, `+0.000777`
- `gq>=20`: `51 TP / 7 FP`, `+0.000725`

另外 `DORADO paired` 有可用的 `Potential_LOH / Coverage_Category`，但：

1. `gq>=10 + Potential_LOH` 是負增益。
2. `gq>=15 + Potential_LOH` 也沒有比單純 `gq>=15` 更好。

這說明 even 在有較完整 LOH proxy 的 paired dataset，中後段 LOH/CNV annotation 目前也比較像背景解釋，而不是 enough-to-relax-GQ 的主規則。

## 結論

### 1. 是否應放寬 GQ？

就目前證據，**不應把 `HCC1395 5kHz paired` 的 primary gate 放寬到 `gq>=12` 以下**。  
如果目標是最高 F1，現有證據反而支持維持 caller-first，甚至偏向更保守的 `gq>=18/20`。

### 2. 為什麼 75x 樣本下 still 不值得？

因為在這個 depth 等級裡：

1. `FP` 已經常有足夠 read 支持，不能指望 read 數自動把它們淘汰。
2. `haplotype / 弱 LOH/CNV proxy` 在 `FP` 也很常見，不能當安全 veto。
3. 真正穩定的分離器還是 caller-side 的 `GQ / QUAL / ref-alt 結構 / LowQual` 狀態。

### 3. downstream 最合理的位置

目前最合理的分工仍是：

1. 第一層：`caller-first`
2. 第二層：`InterSubMod / methylation support / ranking`
3. annotation / QC：`PairwiseMedianDist`、`Quality_Score`、`hp_assign_rate`、弱 LOH/CNV proxy

也就是說，downstream 特徵現在更適合**在安全 gate 之後幫助理解與排序**，而不是拿來交換 primary GQ 的鬆綁。

## 限制

1. `HCC1395 5kHz paired` 本輪沒有 paired candidate-specific 的正式 `Potential_LOH / Coverage_Category`，因此 copy-number / LOH 判讀只能先用弱 proxy。
2. 本輪是 post hoc rescue 模擬，不是把 relaxed caller VCF 重新餵回完整 pipeline 的 full rerun。
3. `DORADO paired` 只能當交叉檢查，不是主樣本的直接替代。

## 下一步建議

1. 若要正式宣稱 relaxed GQ 可行，下一輪必須做 **full rerun**：
   - 重寫 caller input VCF
   - 重新跑 `LongPhase-S`
   - 重新跑 `InterSubMod`
   - 再 benchmark
2. 若短期目標是穩定提升 paired F1，現階段更值得投資的是：
   - 在 `gq>=15/18` 之後做更細的 artifact triage
   - 補齊 `5kHz paired` 的正式 LOH/CNV proxy
   - 對高 `GQ` 但 `LowQual / ref-heavy` 的殘留 FP 做 read-level 回查
