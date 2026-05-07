# T1.2 Read-Level Vote Audit — 4 路徑驗證結果

**Input**: 3 vote dumps from `research/v5_provenance_followup/T1_2_read_level_audit`
**Region**: chr19 (HCC1395 5kHz, V2b PON-only phased VCF as input)
**Binaries**:
  - baseline = `8b8c1fd` (V2b PON-only flag, getVote priority bug 仍在)
  - v3f      = `380e8d2` (V3F two-layer + INDEL guard)
  - v5       = HEAD `938f0df` (Layer 1.5 + ploidy fix + threshold 0.9)

## Summary

| 指標 | 值 |
|------|---:|
| baseline rows | 549,206 |
| v3f rows | 549,206 |
| v5 rows | 549,206 |
| 3-way merged | 1,069,832 |
| 雙向矛盾 reads (germline_maj ≠ somatic_maj, both >0) | 752 |
| **Priority bug confirmed victims** (baseline 跟 somatic) | **752** |
| V3F 修正比例 (改向 germline_maj) | **100.00%** |
| V5 修正比例 (改向 germline_maj) | **100.00%** |

## 4 路徑驗證

### ① 個案 trace（≥10 條）
總數 752 條 — **PASS**（門檻 ≥10）

前 10 個案例（read_name + pos + countMap + hpResult 三版對比）：

| read_name (前 12 chars) | pos | HP1/HP2 | HP1_1/HP2_1 | germline_maj | somatic_maj | hp_b → hp_v3f → hp_v5 |
|---|---:|---|---|:---:|:---:|---|
| 1c50034a-f0f | 201,417 | 1/3 | 1/0 | 2 | 1 | 11 → 21 → 21 |
| afb8e89b-893 | 585,252 | 1/2 | 1/0 | 2 | 1 | 11 → 21 → 21 |
| 35c7e166-ec3 | 824,360 | 0/1 | 1/0 | 2 | 1 | 11 → 21 → 21 |
| cd9ed883-f97 | 854,138 | 0/1 | 1/0 | 2 | 1 | 11 → 21 → 21 |
| 096ab9a7-030 | 1,574,442 | 0/3 | 1/0 | 2 | 1 | 11 → 21 → 21 |
| 120f85f6-6f8 | 2,107,550 | 0/1 | 1/0 | 2 | 1 | 11 → 21 → 21 |
| 7e23e9cc-26d | 2,117,352 | 0/1 | 1/0 | 2 | 1 | 11 → 21 → 21 |
| ccc8185d-f9b | 2,558,240 | 0/1 | 2/0 | 2 | 1 | 11 → 21 → 21 |
| 303ae34f-1ce | 2,560,802 | 0/1 | 2/0 | 2 | 1 | 11 → 21 → 21 |
| 71bcb0c9-2dd | 3,744,892 | 0/1 | 1/0 | 2 | 1 | 11 → 21 → 21 |


### ② 區域聚集 (chr19 1Mb windows)
Top 10 priority bug 比例最高的 windows：

| window 1Mb | bug confirmed | total reads | enrichment |
|---:|---:|---:|---:|
| chr19:30M | 215 | 24,732 | 0.87% |
| chr19:27M | 133 | 23,639 | 0.56% |
| chr19:16M | 41 | 22,078 | 0.19% |
| chr19:14M | 37 | 23,845 | 0.16% |
| chr19:38M | 23 | 18,057 | 0.13% |
| chr19:21M | 27 | 22,114 | 0.12% |
| chr19:31M | 27 | 24,871 | 0.11% |
| chr19:56M | 20 | 18,814 | 0.11% |
| chr19:18M | 22 | 23,371 | 0.09% |
| chr19:29M | 21 | 24,773 | 0.08% |


### ③ Somatic density 共變
分組對比（high somatic vote vs low）：

| 群體 | bug confirmed N | V3F 修正率 |
|---|---:|---:|
| **high** somatic vote ≥5 | 0 | nan% |
| **low**  somatic vote <5 | 752 | 100.00% |

→ **BORDERLINE**：High somatic density reads 不是 priority bug 主要受害者

### ④ 修正後消失
- V3F 修正率 = 100.00% — **PASS (≥80%)**
- V5 修正率  = 100.00% — **PASS (≥80%)**

## 機制因果結論

如果 4 路徑 ≥3 通過 → priority bug 機制因果**確立**；V5 修對是真實。
如果 ≤2 通過 → 機制詮釋需降級。
