---
title: ISM 輔助 tag — 甲基 × haplotype × somatic 觀察推論漏斗
date: 2026-05-31
sample: HCC1395 paired_full（單樣本）+ TO 對比
status: in_progress
audience: PI / advisor
companion_html: 20260531_ISM_aux_tag_observation_funnel_01.standalone.html
cycle_id: 20260531_clustering_source_attribution_Q1Q2Q3
commit: 859c55a
partial_scope: single-sample（cross-sample 未跑）
---

# ISM 輔助 tag — 甲基 × haplotype × somatic 觀察推論漏斗

> **本檔為 grep-able SoT；視覺版見 companion `.standalone.html`（圖片 base64 內嵌、sticky TOC）。**
> 整份調查（2026-05-29 → 05-31）把「somatic 事件能否驅動 read-level 甲基聚集」拆 3 層、用兩方法兩軸測，並先驗一把分群品質尺。

## TL;DR

ISM 能正確**量測**每個 somatic 位點的甲基化結構，但「somatic 驅動 ASM／子克隆甲基聚集」**不成立**（兩方法 B1+B2、兩軸 HP+ALLELE、validated ruler、TO 對比全 NEGATIVE）。BRCA2 等個別位點真實強 cluster，反映 pre-existing 等位甲基化（germline Layer A），非 somatic 驅動。

## 觀察推論漏斗（全樣本 → 篩除 → 保留）— 全數字 L1 親測

| 級 | 保留 | 篩除方法 |
|----|------|---------|
| 0 universe | **39,447** SEQC2 高信度 somatic TP | — |
| 1 進入 ASM | **30,511** variants（51,171 records；valid 51,091） | 篩除 8,936（23%）無 CpG／覆蓋不足 → 無可測軸 |
| 2 名目顯著 | **11,830** records（HP 4,694 / ALLELE 7,136） | Wilcoxon paired p<0.05 |
| 3 Bonferroni | **313**（HP 98 / ALLELE 215） | Bonferroni α=9.79×10⁻⁷ |
| 4 strong-ASM | **172**（HP 57 / ALLELE 115；佔 0.34%） | + effect-size gate \|Δβ\|≥0.1 |

- 軸分佈：HP-axis（somatic-controlled）23,840 ｜ ALLELE-axis（baseline-confounded）27,331；LOH 20,564 / nonLOH 30,607。
- 172 strong 現象：hypo 76 / hyper 96（44.2%，**無方向**）；HP-strong median |Δβ|=0.173、n_cpg=108；BRCA2 HP-axis Δβ=−0.122。
- **B-discrimination 反向**：TP strong 0.337% vs FP strong 1.71% → Fisher OR=0.194, p=1.8×10⁻²⁸（strong-ASM 在 FP 富集 5×）。FP-strong = 低覆蓋+極端 baseline 回歸極值 artifact（單一 regime）。

## 各層結論

- **B1 平均**：門檻 artifact（p<0.05 OR=1.79 → Bonferroni OR=0.81 反轉）；effect 隨 coverage 下降 ρ≈−0.245 → NEGATIVE。
- **H2/H3/H5 深查**（credible regime TP n=198）：穩定但非方向（signed Δ=0.0003 跨 0，sign-test p=0.943，chr-split KS p=0.276）；不能預測 TP（AUC 0.505/credible 0.570 落 null）；LOH-ASM 非 CNV-dose（Spearman ρ=−0.087）。
- **分群品質尺（先驗尺）**：PC1 模擬真 cluster ARI=0.557 偵測 ✓ / NC1 隨機 −0.005 拒絕 ✓ / NORMAL imprinted median 0.758（GNAS/RB1=1.0）→ 尺 VALIDATED。
- **B2 廣掃**：TP 0.135（11/45）/ FP 0.099（9/36）/ FN 0.200（17/45）/ het-NULL 0.177（9/50）；TP vs het MW p=0.922。somatic median << imprinted 尺 0.758 → **NEGATIVE（非特異）**。
- **TO 對比**：TP 0.069 / FP 0.214（18/45）→ FP/TP=3.11×（paired 0.73×）→ 強化 TO germline-FP NO-GO。
- **per-locus 歸因（Q1/Q2/Q3，186 loci）**：9/9 clustered het 全 HP-aligned；ahl 退化 het 1.00→TP 0.83→FP 0.60→FN 0.34（MW p=6.2×10⁻⁹）；非 error/非 subclone（TP 高-ARI 4%<germline 18%）。

## 唯一未關閉 lead（L3，須 BAM follow-up）

allele-only somatic 5/85（vs het 0/50），~2 count-robust：**chr16:34707014**（FN, VAF=0.66, blind_allele=0.945 / blind_hp=0.165, MAPQ60）+ chr14:86405691。從 JSON 無法分辨「真 somatic-ASM」vs「germline phasing 失敗」。

## 取樣條件與限制

1. 單樣本 HCC1395（cross-sample 未驗）。2. 漏斗只覆蓋 77% 有可測軸的 TP（23% 因覆蓋排除）。3. pipeline VCF 的 AF/DP/AD 是 constant placeholder（全 AF=0.2）→ FP allelic-imbalance 機制 asserted 非 measured，降為 hypothesis。4. LOH=whole-region NGS（粗）。5. CN class categorical（無 dose-response）。6. B2 廣掃 curated ~45/class（非 genome-wide 率）。7. 5mC+5hmC max-collapse（未分離）。8. chr16/chr14 lead 須 BAM 定論。

## Provenance

- artifacts：`research/tsg_promoter_asm_reviewer/`（scripts 18-31、output/{cluster_quality_eval_framework, clustering_source_attribution, phase2_dataverify_synthesis}.md、genome_survey_v2/）
- 圖：`figures/{b2_method_validity_proof, brca2_per_cpg_delta, q3_perlocus_map}.png`（q3 本次重繪：圖例移軸外+縮點）
- ledger：`research/autoresearch/evidence_ledger.jsonl` entry 56-64（本調查）
- 所有數字本 session 親自從 JSON/TSV re-derive（L1）。
