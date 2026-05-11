---
title: ClairS-TO Verdict Characterization Pilot — Index
date: 2026-04-20
status: in_progress
sample: HCC1395 subsample t20_n30 (purity ≈ 0.40)
scope: Characterization only · no C++ change · no pipeline filter integration
---

# ClairS-TO Verdict Characterization Pilot

## 觸發

- 2026-04-19 Z3 × HCC1954 amplicon blacklist pilot 確認跨樣本 region-level filter 不可行
- 文獻調查發現 ClairS-TO v0.2.0+ 已內建 Verdict 模組（ASCAT CN segmentation → purity/ploidy → 二項式檢定 → per-variant germline/somatic/subclonal 分類）
- 本專案 v0.3.0 主 pipeline 輸出 VCF 中 Verdict flags 全為空（6 樣本缺 loci / HCC1395 purity>0.8 被跳過）
- 唯獨 HCC1395 t20_n30 subsample（purity=0.40）實際產出 14,875 tags → 現成測試材料

## 假說

| ID | 假說 | Null 門檻 |
|----|------|----------|
| H-V1 | Verdict_Germline 富集於 SEQC2 FP（≥70%） | <55% |
| H-V2 | Verdict_Somatic 與 SEQC2 TP 一致率 ≥85% | >20% FP |
| H-V3 | 若 H-V1 POSITIVE，可推廣至其他 6 樣本 | — |

## 預 flight 觀察（腳本執行前已檢）

ClairS-TO v0.3.0 輸出的 Verdict × FILTER 分佈：

| Verdict class | n | ClairS-TO FILTER |
|---------------|---|------------------|
| Verdict_Germline | 4,633 | **全部 LowQual** |
| Verdict_Somatic | 9,602 | 全部 PASS |
| Verdict_SubclonalSomatic | 640 | 全部 PASS |

- **重要發現**：Verdict_Germline 不會出現在 PASS 集合中 → 啟用 Verdict 對 production PASS set 的直接 FP 移除效益 = 0
- 後續分析仍需驗證：(a) LowQual ∩ Verdict_Germline 是否真為 germline（校準正確性）；(b) PASS ∩ Verdict_Somatic vs PASS ∩ no_Verdict 的 TP rate 差異

## 執行步驟

| Step | 腳本 | 狀態 |
|------|------|------|
| 1 | `scripts/step1_verdict_vs_seqc2.py` — Verdict × SEQC2 交叉 + Wilson CI | pending |
| 2 | `scripts/step2_reference_f1.py` — 假設性過濾 ΔF1 | pending |
| 3 | `scripts/step3_verdict_zone_crosstab.py`（條件） | conditional |
| 4 | Wakhan / SAVANA 文獻（僅 NEGATIVE） | conditional |

## 輸入（唯讀）

- ClairS-TO VCF: `/big8_disk/data/HCC1395/ONT/subsample/t20_n30/ClairS_TO_v0_3_0/snv.vcf.gz`
- SEQC2 sSNV truth: `/big8_disk/data/HCC1395/SEQC2/high-confidence_sSNV_in_HC_regions_v1.2.1.vcf.gz`
- SEQC2 HighConf BED: `/big8_disk/data/HCC1395/SEQC2/High-Confidence_Regions_v1.2.bed`

## 產出

- `data/step1_verdict_vs_seqc2.tsv` — per-variant annotation
- `data/step1_confusion_matrix.tsv` — Verdict × {TP,FP} + Wilson CI
- `data/step2_reference_f1.tsv` — 假設性 ΔF1
- `figures/step1_verdict_fp_rate_per_class.png`
- `figures/step2_reference_f1_barchart.png`

## 最終報告

`docs/experiments/in_progress/2026/04/20260420_ClairS_TO_Verdict_Characterization_Pilot_01.md`
