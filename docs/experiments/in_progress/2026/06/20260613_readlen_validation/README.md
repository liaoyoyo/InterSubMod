---
title: "#2 read-len 驗證 — min_read_length=1000 vs 不限"
date: 2026-06-13
type: experiment
status: DONE
data_sources:
  - baseline/benchmark.jsonl（bench_ism.sh 記錄 wall/RSS/disk/commit）
  - HCC1395 paired chr1 ±5000 BERNOULLI
---

# #2 read 長度限制是否需要？（min_read_length=1000 vs 1）

## L0
✅ **read len≥1000 必要、不該移除**——不限讓保留率 3.47%→1.14%、Strong 降、Noise 升（短 read 稀釋分群品質）；效能完全不省。

## 數據（chr1 ±5000 BERNOULLI，2624 位點）
| min_read_length | 保留率 | Strong | Noise | Weak | NumReads_med | NumCpGs_med |
|---|---|---|---|---|---|---|
| **1000（現狀）** | **3.47%**(91) | 389 | 262 | 1961 | 118 | 78 |
| 1（幾乎不限）| 1.14%(30) | 261 | 286 | 2069 | 126 | 78 |

機制：不限→短 read 進入（NumReads 118→126）→短 read 共同 CpG 少→marginal pair 多→距離雜訊→**分群品質降（Strong 389→261、Noise 262→286、保留 3.47%→1.14%）**。

## 效能（benchmark.jsonl）
| | wall | max RSS | disk |
|---|---|---|---|
| len1000 | 1:08.74 | 3.58GB | 536MB |
| len1 | 1:07.76 | 3.57GB | 556MB |

→ **不限 read 長度不省時間**（短 read 仍需 fetch + 處理；fetch 是 by-region 不論長度）。

## 結論
**保持 min_read_length=1000**——驗證確認是有數據依據的必要設定（不是任意門檻）：移除/降低會降分群品質且不省效能。
