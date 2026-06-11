---
title: ISM 輔助 tag — 圖片觀察 × 判斷依據 互動審查
date: 2026-06-01
sample: HCC1395 paired_full + TO（單樣本）
status: in_progress
audience: PI / self-review
companion_html: 20260601_ISM_aux_tag_figure_review_interactive_01.standalone.html
commit: 859c55a
partial_scope: single-sample
---

# ISM 輔助 tag — 圖片觀察 × 判斷依據 互動審查

> **grep-able SoT；視覺/互動版見 companion `.standalone.html`**（11 圖 base64 內嵌、分頁圖庫、三態 localStorage 審查鈕、燈箱放大）。
> 用途：把整個調查的圖片 + 判斷 + 依據分層攤開，供逐項檢視確認，並記錄「什麼證據會推翻它」供未來校對。

## L0 裁決

ISM 能正確量測 read-level 甲基結構，但「somatic 驅動 ASM/子克隆甲基聚集」不成立——兩方法兩軸 + validated 尺 + TO + 6-agent 紅隊 + 嚴格 patch（FDR/jackknife/coverage/null）全同向，且 correction-robust + CI-quantified + hotspot-robust。BRCA2 真 cluster 反映 germline Layer A。唯一未關閉：chr16:34707014 + 高覆蓋 residual。

## L1 八個判斷（供確認）

| ID | 判斷 | 依據 | tier |
|----|------|------|------|
| J1 | ISM 機器量測 SOUND | BRCA2 blind-ARI 0.79 + 尺 | L1 |
| J2 | 分群品質尺 VALIDATED | PC1 0.557/NC1 −0.005/imprinted 1.0 | L1 |
| J3 | 核心 NEGATIVE（somatic 非特異） | TP 0.135<<0.758；FN>TP；TP−null CI 跨 0 | L1 |
| J4 | strong-ASM anti-discriminative ~3× | OR 5.16→2.84（chr8 61%） | L1 |
| J5 | NEGATIVE correction-robust + CI | Bonf+FDR 都 null≥TP | L1 |
| J6 | BRCA2 真 cluster = germline allelic | b2_cluster MDS + per-cpg | L1/L2 |
| J7 | 唯一 open：chr16 + 高覆蓋 residual | q3 map + p4 crossover | L3 open |
| J8 | ASM = 有效且足夠的 supplement null | 整體嚴格度 | L2 |

## L3 判斷依據 × 校對記錄（未來新證據驗證點）

- **J3 推翻條件**：跨樣本（COLO829）somatic median ARI > imprinted 尺，或 genome-wide matched null 顯示 TP 顯著 > null（CI 下界>0）。
- **J4 推翻條件**：用真測 AF/DP 重算後 FP-enrichment 消失，或 chr8 外 jackknife OR≈1。
- **J5 翻正條件**：genome-wide matched null strong-rate 明顯低於 TP（兩校正一致）。
- **J7 驗證**：BAM 核 chr16 甲基分裂對齊真 phase-block（真 ASM）vs phasing 失敗；H5 full-sample logistic coverage×|Δβ| 交互顯著 + OOF AUC>0.55。
- **J8**：若決策改 standalone-positive 目標（需 C1/C2/C3），充分性門檻提高，J7 須先關閉。

## 遺漏 / 待補

單樣本（cross-sample 未驗）· 5mC/5hmC 未分離（Level1 已刪，需 MSA 重跑，off-goal）· AF/DP placeholder（FP 機制 asserted 非 measured）· obs24_H2 圖 CJK tofu 未嵌（數據已在文字）· null n=1,334 高覆蓋 bin 小（33-163）。

## 最有價值 3 處

1. validated 尺（可重用工程資產）。2. 「FN 比 TP 更聚類」（NEGATIVE 最有力單一證據）。3. chr16:34707014（唯一 VAF-robust allele-orthogonal 候選，待 BAM）。

## Provenance

圖：`research/tsg_promoter_asm_reviewer/figures/{b2_method_validity_proof, q3_perlocus_map, brca2_per_cpg_delta, obs26_H5_roc, obs25_H3_cn_stratification}.png` + `b2_clusters/chr13` + `human_eye/chr13` + `credibility_patch/p{1-4}`。
數據：`genome_survey_v2/{credibility_patch_results, b2_broad_scan_results, enriched_perlocus, obs2*_stats}.json`。
ledger entry 56-65。所有數字本 session 親自重算（L1）。
