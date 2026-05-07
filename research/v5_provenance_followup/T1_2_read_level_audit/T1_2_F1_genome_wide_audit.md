# T1.2 Read-Level Vote Audit — 4 路徑驗證結果（全基因組 (chr1-22 + chrX/Y)）

**Input**: 3 vote dumps from `/big7_disk/liaoyoyo2001/InterSubMod/research/v5_provenance_followup/T1_2_read_level_audit` (suffix=`genome`)
**Region**: 全基因組 (chr1-22 + chrX/Y) (HCC1395 5kHz, V2b PON-only phased VCF as input)
**Binaries**:
  - baseline = `8b8c1fd` (V2b PON-only flag, getVote priority bug 仍在)
  - v3f      = `380e8d2` (V3F two-layer + INDEL guard)
  - v5       = HEAD `938f0df` (Layer 1.5 + ploidy fix + threshold 0.9)

## Summary

| 指標 | 值 |
|------|---:|
| baseline rows | 29,973,253 |
| v3f rows | 29,973,253 |
| v5 rows | 29,973,253 |
| 3-way merged | 59,061,835 |
| 雙向矛盾 reads (germline_maj ≠ somatic_maj, both >0) | 34,855 |
| **Priority bug confirmed victims** (baseline 跟 somatic) | **34,855** |
| V3F 修正比例 (改向 germline_maj) | **100.00%** |
| V5 修正比例 (改向 germline_maj) | **100.00%** |

## 4 路徑驗證

### ① 個案 trace（≥10 條）
總數 34,855 條 — **PASS**（門檻 ≥10）

前 10 個案例（read_name + pos + countMap + hpResult 三版對比）：

| read_name (前 12 chars) | chr | pos | HP1/HP2 | HP1_1/HP2_1 | germline_maj | somatic_maj | hp_b → hp_v3f → hp_v5 |
|---|---|---:|---|---|:---:|:---:|---|
| bba1a3e0-b31 | chr1 | 1,187,904 | 0/25 | 1/0 | 2 | 1 | 11 → 21 → 21 |
| 5d33bc6d-468 | chr1 | 2,325,085 | 0/1 | 1/0 | 2 | 1 | 11 → 21 → 21 |
| 5d33bc6d-468 | chr1 | 2,325,085 | 0/1 | 1/0 | 2 | 1 | 11 → 21 → 21 |
| 5d33bc6d-468 | chr1 | 2,325,085 | 0/1 | 1/0 | 2 | 1 | 11 → 21 → 21 |
| 5d33bc6d-468 | chr1 | 2,325,085 | 0/1 | 1/0 | 2 | 1 | 11 → 21 → 21 |
| 5d33bc6d-468 | chr1 | 2,325,085 | 0/1 | 1/0 | 2 | 1 | 11 → 21 → 21 |
| 5d33bc6d-468 | chr1 | 2,325,085 | 0/1 | 1/0 | 2 | 1 | 11 → 21 → 21 |
| 5d33bc6d-468 | chr1 | 2,325,085 | 0/1 | 1/0 | 2 | 1 | 11 → 21 → 21 |
| 5d33bc6d-468 | chr1 | 2,325,085 | 0/1 | 1/0 | 2 | 1 | 11 → 21 → 21 |
| 357e4fa6-6da | chr1 | 2,771,103 | 1/0 | 0/1 | 1 | 2 | 21 → 11 → 11 |


### ② 區域聚集 (1Mb windows，total ≥100 過濾)
Top 10 priority bug 比例最高的 windows：

| window | bug confirmed | total reads | enrichment |
|---|---:|---:|---:|
| chr19:30 | 215 | 24,732 | 0.87% |
| chr20:29 | 177 | 23,415 | 0.76% |
| chr15:20 | 162 | 22,915 | 0.71% |
| chr21:10 | 239 | 34,160 | 0.70% |
| chrX:62 | 70 | 10,252 | 0.68% |
| chr7:66 | 177 | 26,753 | 0.66% |
| chr2:242 | 41 | 6,285 | 0.65% |
| chr19:27 | 133 | 23,639 | 0.56% |
| chr3:75 | 119 | 21,555 | 0.55% |
| chr6:64 | 71 | 13,626 | 0.52% |


### ③ Somatic density 共變
分組對比（high somatic vote vs low）：

| 群體 | bug confirmed N | V3F 修正率 |
|---|---:|---:|
| **high** somatic vote ≥5 | 26 | 100.00% |
| **low**  somatic vote <5 | 34,829 | 100.00% |

→ **BORDERLINE**：High somatic density reads 不是 priority bug 主要受害者

### ④ 修正後消失
- V3F 修正率 = 100.00% — **PASS (≥80%)**
- V5 修正率  = 100.00% — **PASS (≥80%)**

## ⑤ Per-chromosome priority bug 分布

| chr | bug victims | total reads | enrichment ‰ | rank |
|---|---:|---:|---:|---:|
| chr1 | 2,674 | 4,942,520 | 0.541 | 3 |
| chr2 | 2,792 | 4,525,538 | 0.617 | 2 |
| chr3 | 1,595 | 4,102,109 | 0.389 | 11 |
| chr4 | 2,451 | 3,737,586 | 0.656 | 5 |
| chr5 | 2,194 | 3,506,133 | 0.626 | 6 |
| chr6 | 1,824 | 3,319,118 | 0.550 | 8 |
| chr7 | 3,508 | 4,852,872 | 0.723 | 1 |
| chr8 | 666 | 3,324,020 | 0.200 | 21 |
| chr9 | 1,696 | 2,392,136 | 0.709 | 9 |
| chr10 | 1,010 | 2,452,467 | 0.412 | 15 |
| chr11 | 776 | 2,213,598 | 0.351 | 18 |
| chr12 | 1,024 | 2,406,447 | 0.426 | 14 |
| chr13 | 819 | 1,497,933 | 0.547 | 16 |
| chr14 | 1,653 | 2,266,658 | 0.729 | 10 |
| chr15 | 1,205 | 1,791,787 | 0.673 | 12 |
| chr16 | 2,584 | 2,267,135 | 1.140 | 4 |
| chr17 | 285 | 1,602,470 | 0.178 | 23 |
| chr18 | 1,084 | 1,700,337 | 0.638 | 13 |
| chr19 | 752 | 1,069,832 | 0.703 | 19 |
| chr20 | 2,101 | 1,609,083 | 1.306 | 7 |
| chr21 | 792 | 619,219 | 1.279 | 17 |
| chr22 | 673 | 1,098,683 | 0.613 | 20 |
| chrX | 630 | 1,719,017 | 0.366 | 22 |
| chrY | 67 | 45,137 | 1.484 | 24 |

**chr8 hotspot 驗證**：chr8 enrichment = 0.200‰ vs genome-wide 0.590‰ → 0.34× (低於或等於 genome avg)

## ⑥ Layer 1.5 觸發偵測（genome only）

| 指標 | 值 |
|---|---:|
| germline_vote=0 reads | 21,765,669 |
| V3F tagged 數 | 0 |
| V5 tagged 數 | 560,881 |
| **Layer 1.5 額外觸發**（V5 tag 而 V3F 沒 tag） | **560,881** |
| V5 - V3F 差異 | +560,881 |

→ **Layer 1.5 確實觸發**

## 機制因果結論

如果 4 路徑 ≥3 通過 → priority bug 機制因果**確立**；V5 修對是真實。
如果 ≤2 通過 → 機制詮釋需降級。
