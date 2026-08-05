<!--
建立時間: 2026-07-18
目標: 以 exact integer product 修正 legacy Decimal log1p sensitivity 比較
處理範圍: HCC1395 chr22/chr6 真實 probe；chr1-22 執行介面與 gate
關聯檔案: InterSubMod/research/20260718_k_gt8_read_supported_segmentation/scripts/remediate_exact_log_sensitivity.py
-->

# Exact-log sensitivity remediation 驗證紀錄

## TL;DR

獨立、不可覆寫、fail-closed 的 exact-product 重算已完成，服務 G4（可重現性）與 G5（外部可審核性）。它只修正 `log1p` sensitivity cut 的比較，不會改寫 primary raw-molecule 分割，也不會修改 frozen worker 或 live full outputs。HCC1395 chr22 的 6/6 components 與 chr6 的 83/83 components，legacy log cut 均與 exact-product cut 相同；這兩個真實 probe 沒有發生 cut-level 改判。

## 為何需要 exact product

對每個 read-pattern 的正整數 molecule support `m`：

```text
argmax Σ log(m+1) = argmax log(Π(m+1)) = argmax Π(m+1)
```

`Π(m+1)` 是任意精度整數，可以完全避免 `Decimal.ln()` 近似值與加總順序在 near-tie 時翻轉比較。後續 tie-break 保持不變：

1. retained pattern count 較多；
2. block 較少；
3. cut genomic gap sum 較大；
4. cut indices lexicographic 較小。

Primary cut 仍是 raw molecule sum；exact product 只取代 sensitivity 的 legacy log 比較。

## 重要 schema 發現

只讀 `legacy_components.tsv.gz` 與 `cut_constraints.tsv.gz` 無法在一般情況下重建完整 DP：

- component table 只有 `start/end/k/positions_sha256`；
- constraint table 只含有 read evidence 的位點；
- inactive sites 的實際座標與 local index 不在上述兩檔；
- 因而無法重建合法 block、constraint span 與 genomic-gap tie-break。

因此工具以這兩檔為主要輸入，並強制驗證同 partition 的 `site_membership.tsv.gz` 作為 coordinate witness。缺檔、receipt/SHA 不符、positions semantic SHA 不符、constraint count/mass 不守恆時，會在建立 output directory 前拒絕執行。

## 實作與驗證

新增：

- `InterSubMod/research/20260718_k_gt8_read_supported_segmentation/scripts/remediate_exact_log_sensitivity.py`
- `InterSubMod/research/20260718_k_gt8_read_supported_segmentation/tests/test_remediate_exact_log_sensitivity.py`

未修改：

- `InterSubMod/research/20260718_k_gt8_read_supported_segmentation/scripts/build_k_gt8_partitions.py`
- `InterSubMod/research/20260718_k_gt8_read_supported_segmentation/scripts/read_support_partition.py`
- live full output 內任何檔案

測試命令：

```bash
PYTHONDONTWRITEBYTECODE=1 \
/bip7_disk/liaoyoyo2001/miniconda3/envs/cnvtools/bin/python \
  -m pytest -q -p no:cacheprovider \
  InterSubMod/research/20260718_k_gt8_read_supported_segmentation/tests/test_remediate_exact_log_sensitivity.py
```

實際輸出：

```text
.........
9 passed in 11.37s
Elapsed: 11.90 s
Maximum RSS: 90,388 KiB
Exit status: 0
```

覆蓋：

- red-team：兩組均為 raw support 11、3 patterns、exact product 72，確認用後續 tie-break 決定，不讓近似 log 誤差翻轉；
- 250 randomized cases × 3 objectives 的 brute-force oracle；
- output directory 已存在時拒絕覆寫；
- coordinate witness 被竄改時，在 output 建立前 fail closed；
- full root 沒有 authenticated `_SUCCESS` 時，在 output 建立前 fail closed；
- HCC1395 chr22 與 chr6 真實 partition。

## 真實 probe 結果

| Scope | components | sites | patterns | legacy log = exact product | corrected weight-stable | stability 改判 | wall | max RSS |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| HCC1395 chr22 | 6 | 98 | 371 | 6/6 | 4/6 | 0 | 0.25 s | 22,864 KiB |
| HCC1395 chr6 | 83 | 25,657 | 1,920 | 83/83 | 82/83 | 0 | 9.50 s | 78,896 KiB |

解讀：

- 兩個 probe 的 legacy log cut 都通過 exact-product 驗證；
- 沒有證據顯示 primary raw cut 受 near-tie 影響；
- chr6 的 1 個 weight-sensitive component 仍是 raw/equal/exact 三者切法不同，不是 Decimal artifact；
- exact remediation 是防止未來/其他 chromosome near-tie 的 correctness guard，不應把「probe 沒改判」解讀為 remediation 無必要。

輸出：

- `/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260718_k_gt8_read_supported_segmentation/sensitivity/exact_log_remediation_chr22_v2/`
- `/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260718_k_gt8_read_supported_segmentation/sensitivity/exact_log_remediation_chr6_v2/`

每個 output 含：

- `exact_log_sensitivity.tsv.gz`
- `summary.json` 與 `summary.json.sha256`
- `receipt.json` 與 `receipt.json.sha256`
- `_SUCCESS`

## 完整 chr1–22 執行 gate

工具的 `--full-root` 模式只在來源 root 存在可驗證的 comprehensive `_SUCCESS` 後執行，並強制確認 canonical totals 為 408 components / 47,570 sites。2026-07-18 此紀錄建立時，來源 full run 尚未完成，因此沒有越過 gate 執行全量 exact remediation。

來源完成後的命令：

```bash
/usr/bin/time -v \
/bip7_disk/liaoyoyo2001/miniconda3/envs/cnvtools/bin/python \
  InterSubMod/research/20260718_k_gt8_read_supported_segmentation/scripts/remediate_exact_log_sensitivity.py \
  --full-root /big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260718_k_gt8_read_supported_segmentation/full/HCC1395_chr1_22_v1 \
  --output-dir /big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260718_k_gt8_read_supported_segmentation/sensitivity/HCC1395_exact_log_remediation_chr1_22_v1 \
  --expected-components 408
```

## Claim ceiling

本結果只證明：

- 在既有 ordered contiguous `k≤8` objective 下，可用 exact integer product 正確取代 log sensitivity 比較；
- raw/equal objective 可從 frozen outputs 獨立重現；
- 目前 chr22/chr6 的 legacy log cuts 沒有 cut-level 誤判。

本結果不證明：

- 每個 block 有唯一或真實演化樹；
- VAF 候選樹排序正確；
- clone/subclone 生物學 truth；
- HCC1395 與 HCC1395_DORADO 跨技術一致性。
