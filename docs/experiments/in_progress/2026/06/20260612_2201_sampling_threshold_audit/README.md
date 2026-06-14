---
title: ISM 取樣門檻全量觀察（C_min / min_site_coverage）— 實驗 #1
date: 2026-06-12
type: experiment
status: DONE
data_sources:
  - baseline/threshold_summary.json
  - baseline/cmin_common_cpg_scan.tsv
  - baseline/site_coverage_distribution.tsv
  - raw_summaries/{tp,fp,fn}.run.log
---

# 取樣門檻全量觀察 — 實驗 #1（參數實驗框架 dogfood）

## L0 結論
- ✅ **C_min=3 數值合理保留**（文檔改 3 對齊實跑）
- 🔴 **真正待改 = MAX_DIST→SKIP 策略**：C_min=3 排除 12-20% read-pair，全被設 1.0 當「最遠」= 潛在 clustering artifact
- ◽ **min_site_coverage=5 死參數**：啟用只砍 <0.32% column → wire 清死碼或刪（低優先）
- 💡 **中間帶丟失量化**：二值化 (0.2,0.8) 當 missing，額外造成 2.3% read-pair 被 C_min 排除

## 數據（全量 35435 region / 2.04 億 read-pair）

| C_min=3 排除率 | raw(methylation.csv) | binary(實際距離) | avg 共同 CpG |
|---|---|---|---|
| TP (30490 region) | 9.7% | 12.0% | 12.97 |
| FP (4842 region) | 11.4% | 13.8% | 12.97 |
| FN (109 region) | 17.8% | 20.5% | 13.01 |

C_min 掃描（raw, TP）：1→3.6% / 3→9.7% / 5→19.5% / 10→45.6%（漸進無拐點，median 共同 CpG=10）。
min_site_coverage 掃描：覆蓋<5 砍 0.27-0.32% column（median 覆蓋 77 reads）。

## 產物分類（你的 #3 規範）
- **baseline/**（重複比對主要）：`cmin_common_cpg_scan.tsv`, `site_coverage_distribution.tsv`, `threshold_summary.json`
- **figures/**（重複觀察）：fig1 共同CpG分佈 / fig2 C_min掃描 / fig3 column覆蓋（png，依 no-binary 政策不入 git，留 worktree）
- **raw_summaries/**（這次全量副產品）：TP/FP/FN significance_summary + run.log（FP/FN 比舊 20260420 run 全）
- **record**：manifest.yaml, README.md
- **ephemeral（已刪）**：9.6G per-region 矩陣

## 成本（供未來預估）
全量 TP+FP+FN 開矩陣 = **1762s (29min) / 9.6GB 臨時**（寫 ~60萬檔 I/O 重）。純彙總（不開矩陣）會快很多（用戶基準 3-5min）。

## 下一步
1. **MAX_DIST vs SKIP 對比實驗**（改 `--nan-distance-strategy SKIP` 重跑，比 clustering/significance 差異）
2. C_min 文檔改 3 對齊（走 cpp-change/doc 修正）
3. 位點排除原因記錄（C++ 加 RejectReason，見 governance framework A）
