<!--
建立時間: 2026-07-27
目標: 診斷 HCC1395 chr9 S4–S5 為何未形成 HP2 direct read edge
處理範圍: 單一樣本、單一區域 exploratory diagnostic（PARTIAL，不可外推全樣本）
關聯檔案: raw_bridge_audit_v2.json
-->

# HCC1395 chr9 S4–S5 raw bridge audit

> **PARTIAL / exploratory diagnostic**：只涵蓋 HCC1395 chr9:5,074,198–5,095,768。

## 結論

S4–S5 並非沒有 alignment 跨越。BAM 有 23 個 continuous-reference-span
alignment rows、21 個 distinct QNAME；分到 exact PS=3625170、HP2 後有 8 個
canonical molecules，但只有 2 個在兩端皆為 fixed R/A，低於 strict edge
門檻 3。其餘 6 個為 `LA=5`、`LD=1`，主要是 S4 BQ<20。

HP2 仍可透過 S4–S3（17 molecules）及 S3–S5（3 molecules）形成
threshold-3 transitive component；這不是 PS 或 k≤12 partition cut 所造成。

## 輸入

- BAM: `/big7_disk/liaoyoyo2001/data/HCC1395_5kHz/HCC1395.bam`
- canonical sparse calls:
  `/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260722_exact_ps_k12_hcc1395_pilot/hcc1395_chr1_22_direct_big7_v2/chromosomes/chr9/extraction/HCC1395.chr9.molecule_sparse_calls.tsv.gz`
- MLHP:
  `/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260724_exact_ps_cpp_topology_af_all_samples/all7_strict_guard1000_v1/samples/HCC1395/HCC1395.exact_ps_mlhp.json`

## 執行命令

```bash
python research/20260727_hcc1395_chr9_s4_s5_raw_bridge_audit/scripts/audit_raw_bridge.py \
  --bam /big7_disk/liaoyoyo2001/data/HCC1395_5kHz/HCC1395.bam \
  --sparse-calls /big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260722_exact_ps_k12_hcc1395_pilot/hcc1395_chr1_22_direct_big7_v2/chromosomes/chr9/extraction/HCC1395.chr9.molecule_sparse_calls.tsv.gz \
  --mlhp /big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260724_exact_ps_cpp_topology_af_all_samples/all7_strict_guard1000_v1/samples/HCC1395/HCC1395.exact_ps_mlhp.json \
  --output research/20260727_hcc1395_chr9_s4_s5_raw_bridge_audit/results/raw_bridge_audit_v2.json
```

## 驗證結果

```text
raw span rows = 23
distinct raw QNAME touching both = 21
raw fixed-R/A rows = 12
canonical HP1 = RR 9, RL 2, LR 1, LL 1
canonical HP2 = AA 2, LA 5, LD 1
HP2 threshold-3 path = S4(78) --17-- S3(77) --3-- S5(79)
```

`raw_bridge_audit_v2.json` 為有效輸出。較早的 `raw_bridge_audit.json`
只在 transitive-path 搜尋範圍誤從 S4 起算，會錯失 S3；不得用於
transitive-path 判定，其 raw/canonical pair counts 不受影響。
