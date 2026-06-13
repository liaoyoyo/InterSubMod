---
title: "#3 SNV 早取 — tumor SNV-point fetch"
date: 2026-06-13
type: experiment / perf
status: DONE（結果一致 + 加速，已 commit）
data_sources:
  - baseline/benchmark.jsonl（bench_ism.sh）
  - HCC1395 paired chr1 ±5000 BERNOULLI -j16
---

# #3 SNV 早取：tumor SNV-point fetch（減少讀取量）

## L0
✅ **結果完全一致 + 加速 16%**——tumor fetch 改 SNV point（只取覆蓋 SNV 的 read），結果不變、wall 68.7→57.6s、RSS 略降。

## 改動（理論上不變複雜，只改提取 read）
`RegionProcessor.cpp:761`：tumor fetch region `(region_start, region_end)` 窗範圍 → `(snv.pos, snv.pos+1)` SNV 點。
- 只取覆蓋 SNV 的 read（省 fetch + alt_support ~40% 不過 SNV 的 read）。
- 甲基讀取仍用 window ref_seq（covering read 長、span 窗）。
- **normal 維持 window fetch**（甲基 baseline 需要，不過 SNV 也保留）。
- **alt_support filter 保留為 safety net**（SNV 在 covering read 的 deletion 等 edge case）。
- 下游（甲基/距離/clustering）全不變。

## 驗證（chr1 ±5000 BERNOULLI -j16，2624 位點）
| 指標 | before（window）| after（SNV-point）| 一致？ |
|---|---|---|---|
| NumReads_med | 118 | 118 | ✅ |
| NumCpGs_med | 78 | 78 | ✅ |
| Significant | 91 (3.47%) | 91 (3.47%) | ✅ |
| **wall** | 1:08.74 | **0:57.60 (−16%)** | ⚡ |
| max RSS | 3.58GB | 3.41GB | ⚡ |
| ctest | — | 218/218 pass | ✅ |

## 結論
**結果一致（無漏 read）+ 16% 加速**。tumor SNV-point fetch 是低風險、有效的優化：省 htslib fetch + alt_support 處理 40.7% 不過 SNV 的 tumor read，不影響任何結論。

> ⚠ benchmark.jsonl 的 commit 欄記 git HEAD（commit 時為 a6ded02，#3 尚未 commit）——數值是 #3 改後 binary 實測，僅 commit 標記滯後。
