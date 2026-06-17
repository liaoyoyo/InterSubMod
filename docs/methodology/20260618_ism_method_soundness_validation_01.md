<!--
建立: 2026-06-18
報告類型: 方法健全性驗證 — ISM 核心範式(read×read 甲基距離→階層分群→PERMANOVA)是否合理 + 套用距離矩陣分群驗證概念是否正當 + guardrails
任務類型: D handoff（供論文 Methods「validity / limitations」段；給 reviewer 防守）
狀態: 文獻驗證(L3 citation) + ISM source 親讀交叉確認(check_dispersion / min_reads grep 實證)
data_sources:
  - workflow wt6tsvyla（3-angle method-soundness survey, researcher WebFetch）
  - ISM source: include/core/StructureTest.hpp + src/core/{StructureTest,SignificanceAnalyzer}.cpp（check_dispersion / min_reads_for_permanova grep 實證 2026-06-18）
  - external_validation/_landscape/06 §D/§F-1（longphase-s somatic_haplotag 零甲基 源碼層證實 → HP-label 非循環）
  - 配套方法選單: docs/method_comparison/20260618_distance_matrix_cluster_validation_methods_for_methylation_01.md
provenance_note: 不產新數字；citation 為文獻層(L3)，引用前過 /citation-verification 補一手。ISM 現況 claim 有 source file 佐證(grep)。
-->

# ISM 方法健全性驗證 — read×read 甲基距離 → 分群 → PERMANOVA 是否合理

> **裁決：✅ SOUND_WITH_CAVEATS** — 範式是公認成熟的多變量框架（非 ad-hoc），設計非循環，且最關鍵的 dispersion 防線 ISM C++ 已實作；但有 7 個真實 validity 風險須處理，否則 characterization 結論會被高估。

## 1. 為何合理（established 根基）
| 點 | citation | 適用 ISM |
|---|---|---|
| 「任意 dissimilarity → 幾何分解 → permutation 推論」是數十年成熟框架 | Legendre & Legendre《Numerical Ecology》; Gower 1966 (PCoA, Biometrika 53:325) | ISM read×read 距離矩陣 = 此框架標準輸入 |
| PERMANOVA = canonical（唯一硬假設 = null 下 exchangeable，不需常態/變異同質）| Anderson 2001 (Austral Ecology 26:32-46) | pseudo-F + N-perm + HP 單因子 = 教科書級正確用法 |
| 非歐/半度量距離進 PERMANOVA 被明確支援（additive constant 處理負特徵值，不需先 Euclidean embed）| McArdle & Anderson 2001 (Ecology 82:290) | NHD/L1/Bernoulli/Jaccard 進 PERMANOVA 合法 |
| 微生物組 beta-diversity 大量同構移植先例（事實標準）| vegan::adonis2; QIIME2 adonis | 「CpG 狀態 × read」≈「物種 × 樣本」合理同構 |
| read-level/epiallele 是公認概念 | Landau 2014 PDR; Li 2014 methclone; Scherer 2020 qFDRP(NAR 48:e46) | ISM NHD=qFDRP 核心；Bernoulli+window=cvlr 同做法 |
| UPGMA/Ward 為奠基標準分群 | Sokal&Michener 1958; Ward 1963 (JASA 58:236) | 分群法選擇本身無爭議（Ward 用法見 §3 caveat）|
| **🟢 HP-label 非循環（最強項，源碼證實）**：HP 來自 longphase-s `somatic_haplotag` 路徑，該路徑**零甲基**（external_validation/06 §D/§F-1 grep 證實）| McArdle&Anderson 2001（外生因子合法）| PERMANOVA(甲基距離 ~ HP)=「甲基結構 vs **獨立**基因標籤」對齊檢定，非循環。**前提**：phasing 完全不用甲基當 tie-break |

## 2. 套用「過去概念」（silhouette/clusterboot/consensus）是否正當？→ 合理，3 硬條件
- 這些工具本就為「任意 dissimilarity + 小樣本 + 非常態」設計 → 移植正當。
- **條件①** 匹配 metric：非歐只用 **UPGMA/average-linkage + silhouette/clusterboot**（任意 dissimilarity 一致）；**Ward/k-means/部分 gap-stat 假設 Euclidean → 非歐下語意失效**。
- **條件②** 借 PERMANOVA **必配 PERMDISP**（Anderson & Walsh 2013：兩者是一組）。
- **條件③** **不換檢定對象**：silhouette/consensus 描述「分群品質/穩定度」可以，但不能取代「結構 vs 獨立 HP」對齊檢定，也不能用甲基定義的群當 ground truth 自證（循環）。

## 3. 必守 guardrails（7 風險 + ISM 現況）
| # | 風險 | citation | ISM 現況 / 行動 |
|---|---|---|---|
| 1 | 🔴 **PERMANOVA 混淆 location vs dispersion**（顯著可能只是「甲基較離散」非「型態不同」）| Warton 2012 (MEE 3:89); Anderson&Walsh 2013 (Ecol Monogr 83:557) | ✅ **C++ 已實作 `StructureTest::check_dispersion()`**（distance-to-centroid + ANOVA F + α=0.05 warning；`SignificanceAnalyzer` 呼叫）→ **行動：確認 `enable_dispersion` 開 + warning 隨 pseudo-F 一併報告（勿默默丟；連動 06-16 FN 審查）** |
| 2 | 🔴 **不平衡設計放大偏誤**（小組高 dispersion→偽陽性）；HP1/HP2 幾乎總不平衡(LOH skew) | Anderson&Walsh 2013 | ⚠ **行動：每 region 報 HP1/HP2 read 數 + dispersion 方向；不平衡且 dispersion 顯著的結論降級** |
| 3 | 🔴 小樣本 permutation p 有硬下限 1/#perm；**未顯著≠無結構** | permutation 理論; Anderson 2001 | ✅ `min_reads_for_permanova=5` + `C_min≥3`；⚠ 行動：報 permutation 數、小組標 low-power、勿報精細 p |
| 4 | 🔴 **Ward × 非歐**：Ward 目標函數假設 squared-Euclidean，非歐下變異分解語意失效 | Murtagh&Legendre 2014 (J Classif 31:274); SciPy docs | ⚠ **行動：非歐 metric 分群改 UPGMA/average-linkage；Ward 結果僅標視覺啟發** |
| 5 | Jaccard 非對稱（排除雙非甲基 M00 共享資訊）| Legendre&Legendre double-zero | ⚠ 行動：primary 用對稱 NHD/L1/Bernoulli；Jaccard 限「只關心共甲基」並明示假設 |
| 6 | 🔴 **per-CpG Fisher 忽略 over-dispersion**（同 read 多 CpG 相關 + ONT 噪音 → 低估 p）| NanoMethPhase(Genome Biol 22:68); DSS | 🔴 **KB 既列「必修」**：至少 FDR 校正，理想改 beta-binomial/DSS |
| 7 | read exchangeability（重複/同分子 read 違反獨立性）| Anderson 2001 | ⚠ 行動：dedup amplicon/重複 read（ONT simplex 無 PCR 風險較低）；CpG 間相關性本身非問題 |

## 4. ISM 紅線一致 + 對外措辭
- 範式定位 **characterization**（甲基結構是否與獨立 haplotype 對齊），與紅線一致——**甲基非重建驅動、filter 方向 DEAD**，不可用此回頭宣稱判別 TP/FP variant。
- 🔴 **對外勿宣稱「ISM 獨有顯著性檢定／對手缺檢定」**（cvlr Table1 randomization P、ASMS pval.rs 兩個 permutation 皆有，且皆 ONT-capable；external_validation/06 §A-1/§F-1 源碼證實）。ISM 真正增量三件 = **① 對 read×read 距離矩陣「結構」做 PERMANOVA（非 cluster 間甲基差量）② normal-baseline somatic cis-test ③ somatic-subclone 目標**——是「檢定對象/錨點」而非「有無檢定/平台」。

## 5. 一句結論
用 ISM 這方法做 characterization **合理且站得住**（成熟框架 Legendre/Anderson + 非循環設計 + read-level 有 cvlr/qFDRP 先例 + dispersion 防線已在 C++）；套用過去驗證概念也合理，前提守「非歐用 UPGMA、PERMANOVA 配 PERMDISP、不換檢定對象」。最高優先行動：確認 `check_dispersion` 結果有被報告（接 06-16 FN 審查）+ 修 per-CpG Fisher over-dispersion（KB 必修）。配套方法選單見 `InterSubMod/docs/method_comparison/20260618_distance_matrix_cluster_validation_methods_for_methylation_01.md`。
