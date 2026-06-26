<!--
建立時間: 2026-06-27
類型: L4 ISM 甲基 corroboration — 單變量（clone/subclone 整合報告）
狀態: in_progress（HCC1395 ⭐3；🔴 既有甲基輸出覆蓋不足，genome-wide 待重抽）
data_sources: data/sm_methyl_corroboration.json
-->

# L4 — ISM 甲基 corroboration（單變量；既有輸出可用性 + 小樣本結果）

> 問題：① 既有 ISM 甲基輸出能否直接用於全基因組結構區的甲基 corroboration？② 甲基對遺傳定義群的輔助程度？
> 🔴 紅線：甲基 = **corroborate 非 detect**；只在 **genotype-anchored**（用 sSNV 定好的群，非 clustering）才非循環。
> 圖：`figures/05_methylation_corroboration.png`。資料：`data/sm_methyl_corroboration.json`。

## §1 決定性發現 — 既有甲基輸出覆蓋不足（不能直接全基因組用）

- 既有 ISM 甲基窗 = **1,139 個 cis-candidate anchor 窗**（records_v6 候選位置，**±5kb 寬**）。
- 這組位置與我的 TP∪FP 連鎖 sSNV **是不同集合** → 窗不對齊我的區域。
- 實測：CN-clean 結構區（span≤8kb）**323 個中，只有 9 個**有 covering candidate 窗（其餘 314 缺）。

→ **既有甲基輸出無法直接做全基因組區域甲基 corroboration**。🔴 **誠實分母（修對抗稽核）**：9/323 = 2.8% 是「CN-clean + span≤8kb」這個最有利子集的覆蓋率；對**全部 4,678 個有結構區域**，覆蓋率僅 **9/4,678 = 0.19%**。**全基因組需從 BAM（MM/ML tag）重抽甲基、對齊連鎖區域** = pipeline 整合步驟（後補；對應你「缺的後續整合成流程」）。

## §2 可測的 8 區小樣本結果（genotype-anchored，非循環）

方法：re-pileup 區域 somatic sSNV → per-read 基因型群（genotype-anchored）→ join methylation.csv（by read_name）→ 兩大 population 間 per-CpG Mann-Whitney U + BH-FDR（sig = q<0.05 且 |Δβ|≥0.2）。cis-ASM proxy：sig CpG 是否同時在 HP1 vs HP2 區分（= germline cis 非 subclone）。

| 指標 | 值 |
|---|---|
| 可測區（有窗+≥2 population）| 8 |
| 甲基 corroborate（≥1 sig CpG）| **4（50%）** |
| median sig CpG / 測試區 | 1.5 |
| cis-ASM 可解釋的 corroboration | **0 / 4（0%）= 這 4 區是 subclone-specific（非 germline cis）** |

🔴 **n=8 過小，僅 anecdotal**，不可當 genome-wide 比例。方向與 chr17 canonical 一致（甲基能區分遺傳群、且非 cis），但**統計力不足以下定論**。

## §3 chr17 canonical（已驗證，established）
chr17:48360161（有窗）先前完整分析：甲基 BERNOULLI/UPGMA **只切得出「有無 somatic」（L0/γ vs L1+L2），切不出 L1↔L2 fine subclone**；per-CpG 可歸因（23 CpG→α 祖先 ASM / 6→L1-L2）。→ 甲基 = **corroborate（刻畫已驗群）非 detect**，且細分 subclone 需遺傳錨。

## §4 data-supported 結論（誠實）
1. **既有甲基輸出 = 稀疏 anchor 窗，覆蓋 ~3% 結構區** → genome-wide 甲基 corroboration **尚不可用既有輸出**，需重抽（後補）。
2. 小樣本（4/8）+ chr17 canonical **暗示**甲基可 corroborate 遺傳群且部分 subclone-specific，但 **n 不足、不下定論**。
3. 🔴 甲基 corroborate 非 detect；細分 subclone 必須靠遺傳（sSNV 連鎖 + HP）。

## §5 後補流程（待執行，已規劃）
- 從 tumor BAM 解析 MM/ML（5mC）tag → 對**每個連鎖區域**重抽 read×CpG 甲基矩陣（對齊我的 region，非 ISM anchor 窗）。
- 全基因組 genotype-anchored per_cpg_diff + NACT cis-control（normal）→ 真 genome-wide 甲基輔助程度。
- 這是論文「methylation profiles = characterize」一章的完整數據來源。
