<!--
建立時間: 2026-07-18
目標: 獨立重算 positional singleton 分母、甲基多群比例與 claim ceiling
處理範圍: 7 datasets / 6 biological samples / chr1-22 / 469,849 dataset-sites
關聯檔案:
  - InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/results/independent_m2_gate_recount.v3.json
  - InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/claim-contract-v5.md
  - /big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/positional_singleton_methyl_multigroup_audit_v1_source_attested/positional_singleton_audit_summary.json
-->

# Positional singleton 數值與 claim 獨立 agent 稽核

## 審查資訊

- Reviewer agent: `019f72b4-530b-7440-887d-b09a41b2a198`
- 審查方式: 不呼叫 production gate 函式，直接讀取 site TSV、stable assignments、
  independent M2 receipt 與 claim contract 後獨立重算。
- Verdict: **數值與母數 APPROVE；若外推 `30/48` 或稱為 subclone / linear
  evolution，敘述必須 FIX。**

## 獨立重算結果

`positional singleton` 定義為同一 dataset、同一 chromosome，以相鄰距離
`<=50,000 bp` 建立 transitive connected component 後 `component_size=1`。

| 指標 | 分子 / 分母 | 比例 |
|---|---:|---:|
| positional singleton / 全 dataset-sites | 50,432 / 469,849 | 10.733661% |
| M1-evaluable / singleton | 48,347 / 50,432 | 95.865720% |
| M1 stable multigroup / singleton | 5,961 / 50,432 | 11.819876% |
| M1 stable multigroup / M1-evaluable | 5,961 / 48,347 | 12.329617% |
| M2-evaluable / M1 multigroup | 48 / 5,961 | 0.805234% |
| M2 PASS / singleton | 30 / 50,432 | 0.059486% |
| M2 PASS / selected M2-evaluable | 30 / 48 | 62.5% |

守恆重算：

```text
30 PASS + 18 FAIL + 5,913 NOT_EVALUABLE + 44,471 NOT_RUN = 50,432
```

M1 不可評估 `2,085 = 2,069` 個 matrix-joined focal-ALT reads 不足，加上
`16` 個 distance 資訊不完整。singleton 的最小有限 nearest gap 為
`50,003 bp`，與 50 kb contract 一致。

## 分層結果

| Truth stratum | Singleton | M1-evaluable | M1 multigroup | M2 P / F / NE / NR |
|---|---:|---:|---:|---:|
| TP | 45,193 | 44,171 | 5,494 | 30 / 16 / 5,448 / 39,699 |
| FP | 1,084 | 514 | 52 | 0 / 1 / 51 / 1,032 |
| UNASSESSED | 4,155 | 3,662 | 415 | 0 / 1 / 414 / 3,740 |

FP 的 M2-evaluable 只有 `n=1`，因此不能估計 specificity，也不能宣稱 M2
已能區分真實與偽陽性位點。

## 技術資料集比較

HCC1395 與 HCC1395_DORADO 是同一 biological sample，不是兩個獨立生物重現。

| 層級 | Intersection | Union | Jaccard |
|---|---:|---:|---:|
| singleton loci | 7,484 | 9,116 | 82.10% |
| M1 multigroup loci | 407 | 1,289 | 31.57% |
| M2 PASS loci | 0 | 4 | 0% |

兩邊各有 2 個 M2 PASS，但沒有共同 PASS locus，因此目前沒有 locus-level
technical replication。

## Claim ceiling

1. `48,347 / 50,432` 支持「技術上可完成 focal-ALT 甲基 screen」。
2. `5,961 / 50,432` 支持「有 operational stable methyl-multigroup signal」。
3. `30 / 50,432` 支持「八個 measured axes 下仍有 residual read-level
   epigenetic partition」。
4. 目前可確認的 cellular subclone 數為 `0`；可確認的 linear ancestry
   順序也為 `0`。

共同 ancestral ALT 被兩個後代 clone 攜帶、再由甲基狀態分群，是與觀察相容的
模型；但同一觀察也可由 epigenetic state、cis-ASM、CN、purity、read geometry
或未測因素產生。單一 focal ALT 缺乏 multi-marker co-segregation，不能識別
clone 數或 linear / branching order。

## 必須保留的限制

- `50,432` 是 positional singleton，不可改稱所有可能定義下的
  read-sharing degree-zero。
- `30/48` 只描述高度選擇後的 conditional subset，不可外推為 population rate。
- 7 datasets 只代表 6 biological samples。
- G1、G2、formal R1 在這 50,432 個 singleton 中全部 `NOT_RUN`。
