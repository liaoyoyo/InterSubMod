<!--
建立時間: 2026-07-18 08:20
目標: 獨立驗證 k>8 分割後各 HP/PS observed constraint unit 的證據保留率
處理範圍: HCC1395；先 chr6+chr22 probe，canonical _SUCCESS 後才執行 chr1-22
關聯檔案:
  - InterSubMod/research/20260718_k_gt8_read_supported_segmentation/scripts/audit_hp_ps_unit_retention.py
  - InterSubMod/research/20260718_k_gt8_read_supported_segmentation/tests/test_audit_hp_ps_unit_retention.py
-->

# HP/PS unit retention 獨立紅隊稽核

## 1. 任務分類與服務目標

- Task Type：**B — Comprehensive validation**
- 當前狀態：probe 已完成；chr1–22 full 尚待 canonical `_SUCCESS`
- 服務目標：**G4 多樣本一致性與 reproducibility**、**G5 可被外部驗證的證據契約**
- 本稽核不修改 full runner、extractor、partitioner、read-support solver 或其既有輸出。

研究啟動 5 問：

1. Thread D read-level epigenetic：否；本輪只有 read-linked SNV constraint retention。
2. Thread B 撤回範圍：否。
3. KDE-corrected：不適用。
4. VCF caller AF：不需要；VAF/AF 不參與切割保留率。
5. 長計算/C++/搬移/NO-GO：無 C++、無搬移；full audit 只在上游長計算完成並寫出 `_SUCCESS` 後啟動。

## 2. 研究問題

全體 component 的 aggregate read-support retention 可能掩蓋 HP1、HP2 或特定 PS unit 的低保留情形。本稽核要回答：

1. 每個 `dataset × chromosome × legacy component × HP family × exact PS` observed unit 保留多少 exact pattern 與 molecule weight？
2. HP1、HP2 在同一 legacy component 的 retention 是否大致相同？
3. 哪些低 retention 是切割造成，哪些是 `k≤8` 下本來就無法完整容納的 unavoidable pattern？

## 3. 假設、成功與失敗條件

### 假設

1. `cut_constraints.tsv.gz` 已是 upstream solver 的 exact-pattern sufficient statistic。
2. `molecule_weight` 的正確 aggregate 單位是 **molecule × component incidence**，不是 unique read/molecule。
3. 不同 HP family 與 exact phase set 不得先混池。

### 成功條件

1. 每個 observed unit 均滿足：
   - `total = retained + cut_lost + unavoidable`
   - pattern rows 與 molecule-incidence weight 各自守恆
2. source receipt、SHA-256、scope、HP/PS contract 與 upstream count 全部一致。
3. 未觀察到 constraint 的 component 不被回填為 0% 或 100% retention。
4. full mode 只接受 canonical HCC1395 chr1–22 `_SUCCESS`。

### 失敗條件

- 未知 disposition、HP/PS 混池、空 phase set、output hash drift、count/mass 不守恆、scope 不完整，任一項即 fail closed。

## 4. 統計口徑

### 最小 unit

`dataset × chrom × legacy_component_id × hp_family × exact known phase_set`

### 分母

- `total_pattern_rows`：該 observed unit 的 exact pattern rows。
- `total_molecule_component_incidence_weight`：各 pattern 的 molecule weight 加總。
- `retention_ratio = retained_weight / total_weight`。

### disposition

- `retained`：完整落在單一新 block。
- `cut_lost`：`disposition == cut`。
- `unavoidable`：`disposition` 以 `unavoidable_` 開頭。
- `nonretained = cut_lost + unavoidable`。

### 支持度與 headline

- disjoint strata：`1–4`、`5–19`、`20–49`、`≥50`
- headline eligible：`total_weight ≥ 20 AND total_pattern_rows ≥ 5`

### 無 observed constraint 的 component

本版不從 sparse calls 反推「應存在但沒有 ≥2 fixed R/A observations」的 HP/PS unit，因為 partition output 沒有提供獨立且無歧義的 zero-opportunity denominator。這些 component 只計入：

`components_without_observed_constraint_units`

不輸出合成 unit row，也不給 0%/100% retention。

## 5. Step → Verify

1. 驗證 source receipt/SHA 與 upstream identities  
   → 驗證：任何 file drift 均在讀取統計前 fail closed。
2. 逐 component 驗證 membership 與 constraint  
   → 驗證：site membership、active-site union、pattern row、molecule mass 與 receipt 全一致。
3. 先按 exact HP/PS unit 聚合  
   → 驗證：`unit_id == HP{family}|PS{phase_set}`，HP 僅 1/2，PS 非空。
4. 同 component 內分別加總 HP1/HP2 後配對  
   → 驗證：輸出 pair table 與 signed/absolute delta。
5. full mode 驗證 canonical `_SUCCESS`  
   → 驗證：22 autosomes、79,687 sSNV、408 k>8 components、47,570 sites、21 completed + 1 zero-target skip。

## 6. chr6 + chr22 probe

### 輸入

- chr6：
  `/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260718_k_gt8_read_supported_segmentation/full/HCC1395_chr1_22_v1/chromosomes/chr6/partition`
- chr22：
  `/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260718_k_gt8_read_supported_segmentation/probes/HCC1395_chr22/partition_v2`

### 執行命令

```bash
/usr/bin/time -v \
  /bip7_disk/liaoyoyo2001/miniconda3/envs/cnvtools/bin/python \
  /big7_disk/liaoyoyo2001/InterSubMod/research/20260718_k_gt8_read_supported_segmentation/scripts/audit_hp_ps_unit_retention.py \
  --mode probe \
  --partition-dir /big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260718_k_gt8_read_supported_segmentation/full/HCC1395_chr1_22_v1/chromosomes/chr6/partition \
  --partition-dir /big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260718_k_gt8_read_supported_segmentation/probes/HCC1395_chr22/partition_v2 \
  --output-dir /big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260718_k_gt8_read_supported_segmentation/validation/hp_ps_unit_retention_probe_chr6_chr22_v5
```

### 輸出

`/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260718_k_gt8_read_supported_segmentation/validation/hp_ps_unit_retention_probe_chr6_chr22_v5/`

主要檔案：

- `HCC1395.hp_ps_observed_units.tsv.gz`
- `HCC1395.hp1_hp2_paired_components.tsv.gz`
- `summary.json`
- `summary.tsv`
- `receipt.json`
- `receipt.json.sha256`

### 實際輸出片段

```text
all_pass: true
observed_constraint_units: 57
eligible_headline_units: 44
hp1_hp2_paired_components: 17
weighted_retention_ratio: 0.718342942558
receipt_sha256: 818f45ddb61874035a53b163066db770e918249c4e1b40754818a73d252e79e2
```

Runtime：

- wall：1.39 秒
- max RSS：33,052 KiB
- exit status：0

## 7. Probe 結果

### Overall

| 指標 | 數值 |
|---|---:|
| k>8 components | 89 |
| 有 observed HP/PS constraint unit | 19 |
| 無 observed HP/PS constraint unit | 70 |
| observed units | 57 |
| headline eligible units | 44 |
| molecule×component incidence total | 4,683 |
| retained | 3,364 |
| cut lost | 227 |
| unavoidable | 1,092 |
| weighted retention | 71.8343% |

Component HP coverage：

| 類型 | 數量 |
|---|---:|
| HP1 and HP2 都有 observed unit | 17 |
| HP1 only | 1 |
| HP2 only | 1 |
| 無 observed HP1/HP2 unit | 70 |

按 chromosome 分開後差異很大：

| chromosome | observed units | eligible | incidence total | retained | weighted retention | unit median |
|---|---:|---:|---:|---:|---:|---:|
| chr6 | 41 | 29 | 3,532 | 2,327 | 65.8834% | 67.5676% |
| chr22 | 16 | 15 | 1,151 | 1,037 | 90.0956% | 98.9904% |

### Unit distribution

| subset | n | median | `<0.5` | `=1` |
|---|---:|---:|---:|---:|
| all observed units | 57 | 94.4444% | 13 | 19 |
| headline eligible | 44 | 94.2868% | 12 | 12 |

解讀：中位 unit 保留率高，但 aggregate weighted retention 僅 71.83%，是因 chr6 的大型 legacy chain 集中大量低 retention / unavoidable incidence。不能只報一個 aggregate 比例。

### HP1/HP2 paired component

| 指標 | 數值 |
|---|---:|
| paired components | 17 |
| pair 兩側皆 headline eligible | 13 |
| absolute delta median | 1.25 percentage points |
| absolute delta `<5 pp` | 12 / 17 |
| absolute delta `<10 pp` | 13 / 17 |
| absolute delta `≥25 pp` | 2 / 17 |
| maximum absolute delta | 38.26 pp |

這支持「多數 paired component 的 HP1/HP2 retention 接近」，但不支持宣稱所有 HP 都公平保留；大型或局部 evidence-sparse component 仍有明顯 HP-specific 差異。

## 8. 測試

輸入：

`InterSubMod/research/20260718_k_gt8_read_supported_segmentation/tests/test_audit_hp_ps_unit_retention.py`

命令：

```bash
PYTHONDONTWRITEBYTECODE=1 \
/bip7_disk/liaoyoyo2001/miniconda3/envs/cnvtools/bin/python \
  -m pytest -q -p no:cacheprovider \
  InterSubMod/research/20260718_k_gt8_read_supported_segmentation/tests/test_audit_hp_ps_unit_retention.py
```

實際輸出：

```text
.........                                                                [100%]
9 passed in 1.52s
```

另以不載入 Python audit 模組的 `gzip + awk` 重算原始 constraint tables：

```text
observed_units                  57
pattern_total                 2291
pattern_retained              1035
pattern_cut                    192
pattern_unavoidable           1064
weight_total                  4683
weight_retained               3364
weight_cut                     227
weight_unavoidable            1092
hp1_hp2_paired_components       17
unknown_dispositions              0
```

與 `summary.json` 完全一致。

完整相關測試矩陣（segmentation、span cap、summarizer、report builder、HP/PS audit）：

```text
........................................................................ [ 79%]
...................                                                      [100%]
91 passed in 21.00s
```

## 9. 報告整合建議

不修改 `build_report_artifact.py`，但正式 HTML 建議增加獨立章節：

1. KPI：
   - observed HP/PS units
   - headline eligible units
   - weighted retention
   - components without observed units
2. 分佈圖：
   - eligible unit retention histogram / ECDF
   - 支持度 strata × retention
3. HP fairness：
   - paired component HP1 vs HP2 scatter，附 `y=x`
   - absolute delta 分佈與 top outliers
4. 必備註解：
   - aggregate 是 molecule×component incidence
   - zero-opportunity unit 未重建
   - retention 只衡量工程切割後 constraint 保存，不等於 topology 正確率
   - VAF 未參與本稽核；VAF 只能在 downstream candidate-tree ranking 使用

## 10. 現階段結論

Probe 支持這個 audit 方法可實作、守恆且可追溯。它補出 aggregate component 統計看不到的 HP/PS unit 尾部風險；但目前仍是 chr6+chr22 partial evidence，必須等 full runner 的 canonical `_SUCCESS` 後才能下 HCC1395 chr1–22 結論。
