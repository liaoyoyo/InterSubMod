<!--
建立時間: 2026-03-25 14:25
目標: 針對 Phase 1A 進行更嚴格的 context-only vs methyl+context incremental test，並擴到更大 shard 與其他 biological samples，確認甲基特徵是否真有穩定增益
處理範圍:
  - sample400 mode-mixed 大 shard benchmark
  - paired-only multi-bio external validation
  - paired comparison / McNemar / region bootstrap CI
  - dataset-level 增益拆解
關聯檔案:
  - /big7_disk/liaoyoyo2001/InterSubMod/scripts/analysis/run_phase1a_read_classifier_benchmark.py
  - /big7_disk/liaoyoyo2001/InterSubMod/scripts/analysis/run_phase1a_incremental_test.py
  - /big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260325_phase1_manifest_shard_export_sample400
  - /big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260325_phase1a_read_classifier_benchmark_sample400_v1
  - /big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260325_phase1a_incremental_test_sample400_v1
  - /big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260325_phase1_training_manifest_paired_multibio_v1
  - /big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260325_phase1_manifest_shard_export_paired_multibio_sample637
  - /big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260325_phase1a_read_classifier_benchmark_paired_multibio_sample637_v1
  - /big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260325_phase1a_incremental_test_paired_multibio_sample637_v1
-->

# Phase 1A incremental test 與 multi-bio validation round 2

## 1. 本輪要回答的問題

round 1 已證明 `Phase 1A` 可行，但還不能回答更關鍵的問題：

1. `methyl + context` 相對於 `context-only` 到底是穩定增益，還是只是小樣本波動？
2. 這個差異在更大 shard 下會變得更清楚，還是會反轉？
3. 增益能否跨 `biological samples` 成立，而不是只停留在 `HCC1395`？
4. 若甲基特徵不是全域有效，邊界條件是什麼？

## 2. 本輪設計

### 2.1 資料

mode-mixed 大 shard：

- [/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260325_phase1_manifest_shard_export_sample400](/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260325_phase1_manifest_shard_export_sample400)
- `400` 個 regions
- `58,795` 筆 reads
- discovery: `HCC1395 5kHz paired + TO`
- validation: `HCC1395 DORADO paired + TO`

paired-only multi-bio shard：

- [/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260325_phase1_training_manifest_paired_multibio_v1](/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260325_phase1_training_manifest_paired_multibio_v1)
- [/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260325_phase1_manifest_shard_export_paired_multibio_sample637](/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260325_phase1_manifest_shard_export_paired_multibio_sample637)
- `637` 個 regions
- `86,521` 筆 reads
- discovery: `HCC1395 5kHz paired`
- external validation:
  - `HCC1395 DORADO paired`
  - `COLO829 PAO paired`
  - `H1437 Google paired`
  - `H2009 Google paired`
  - `HCC1937 Google paired`
  - `HCC1954 Google paired`

### 2.2 比較方法

benchmark：

- [run_phase1a_read_classifier_benchmark.py](/big7_disk/liaoyoyo2001/InterSubMod/scripts/analysis/run_phase1a_read_classifier_benchmark.py)
- 新增 `benchmark_dataset_metrics.tsv` 與 `read_id` 級 prediction 輸出，讓 dataset-level 解讀更直接

incremental test：

- [run_phase1a_incremental_test.py](/big7_disk/liaoyoyo2001/InterSubMod/scripts/analysis/run_phase1a_incremental_test.py)
- 固定比較：
  - `logistic_context_only`
  - `logistic_methyl_context`
- 評估方式：
  - `McNemar paired test`
  - 以 `region_key` 為單位做 `500` 次 bootstrap
  - 輸出 `delta_F1`、`delta_accuracy`、`95% CI`、`positive_fraction`

## 3. 主要輸出

mode-mixed sample400：

- benchmark：
  [/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260325_phase1a_read_classifier_benchmark_sample400_v1](/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260325_phase1a_read_classifier_benchmark_sample400_v1)
- incremental：
  [/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260325_phase1a_incremental_test_sample400_v1](/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260325_phase1a_incremental_test_sample400_v1)

paired-only multi-bio sample637：

- benchmark：
  [/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260325_phase1a_read_classifier_benchmark_paired_multibio_sample637_v1](/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260325_phase1a_read_classifier_benchmark_paired_multibio_sample637_v1)
- incremental：
  [/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260325_phase1a_incremental_test_paired_multibio_sample637_v1](/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260325_phase1a_incremental_test_paired_multibio_sample637_v1)
- model-error compare：
  [/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260325_phase1a_model_error_compare_paired_multibio_sample637_v1](/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260325_phase1a_model_error_compare_paired_multibio_sample637_v1)

sample400 `to-pure` failure diagnosis：

- `context-only`：
  [/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260325_phase1a_sample400_to_pure_failure_context_only_v1](/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260325_phase1a_sample400_to_pure_failure_context_only_v1)
- `methyl+context`：
  [/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260325_phase1a_sample400_to_pure_failure_methyl_context_v1](/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260325_phase1a_sample400_to_pure_failure_methyl_context_v1)

## 4. 結果

### 4.1 sample400：更大 mode-mixed shard 下，甲基增益不成立

整體結果：

| split | context-only F1 | methyl+context F1 | delta F1 | bootstrap 95% CI | McNemar p-value |
| --- | ---: | ---: | ---: | --- | ---: |
| discovery_holdout | `0.8331` | `0.8244` | `-0.0088` | `[-0.0227, -0.0005]` | `7.71e-09` |
| external_validation | `0.9312` | `0.9106` | `-0.0206` | `[-0.0352, -0.0095]` | `1.30e-118` |

dataset-level external validation：

| dataset | delta F1 | bootstrap 95% CI | 判讀 |
| --- | ---: | --- | --- |
| `HCC1395 DORADO paired` | `0.0000` | `[0.0000, 0.0000]` | 無差異 |
| `HCC1395 DORADO TO` | `-0.0328` | `[-0.0537, -0.0145]` | 明確負增益 |

補充 benchmark：

- `logistic_context_only`
  - discovery holdout `F1=0.8331`
  - external validation `F1=0.9312`
- `logistic_methyl_context`
  - discovery holdout `F1=0.8244`
  - external validation `F1=0.9106`
- 在 `ONT_Dorado|to-pure`：
  - `context-only F1=0.9079`
  - `methyl+context F1=0.8751`

結論：

1. 在 `paired + to-pure` 混合任務下，甲基特徵沒有帶來增益。
2. 主要負面來源不是 `paired-pure`，而是 `to-pure`。
3. 因此不能把 `methyl+context` 當成所有 `Phase 1A` 任務的預設更優方案。

### 4.1.1 sample400 `to-pure` 大 shard failure diagnosis

`ONT_Dorado|to-pure` 的整體 error rate：

- `context-only`: `9.78%`
- `methyl+context`: `13.64%`

主要惡化 bucket：

| bucket | context error rate | methyl error rate | delta |
| --- | ---: | ---: | ---: |
| `TP + Noise + PassedGating=False + REF` | `0.3853` | `0.8543` | `+0.4690` |
| `TP + Strong + PassedGating=True + REF` | `0.1237` | `0.1706` | `+0.0469` |
| `FP + Subclone + PassedGating=True + ALT` | `0.4014` | `0.4965` | `+0.0950` |
| `FP + Strong + PassedGating=True + REF` | `0.1932` | `0.2074` | `+0.0143` |

這表示 `to-pure` 的負增益不只是少數邊界樣本，而是會同時：

1. 增加一部分 `TP REF` 誤判
2. 放大 `FP Subclone / Strong` 的 ALT 偏誤
3. 在 `Noise / Weak` 子群再引入額外錯誤

### 4.2 paired-only multi-bio sample637：跨 biological samples 出現穩定小幅正增益

整體結果：

| split | context-only F1 | methyl+context F1 | delta F1 | bootstrap 95% CI | McNemar p-value | 是否支持增益 |
| --- | ---: | ---: | ---: | --- | ---: | --- |
| discovery_holdout | `0.9195` | `0.9269` | `+0.0073` | `[-0.0078, +0.0323]` | `2.91e-02` | 否 |
| external_validation | `0.8611` | `0.8722` | `+0.0112` | `[+0.0044, +0.0188]` | `1.44e-51` | 是 |

補充 benchmark：

- `context-only`
  - discovery holdout `F1=0.9195`
  - external validation `F1=0.8611`
- `methyl+context`
  - discovery holdout `F1=0.9269`
  - external validation `F1=0.8722`
- `harmonization_group` 外部驗證：
  - `ONT_Dorado|paired-pure`: `0.9597 -> 0.9652`
  - `ONT_Google|paired-pure`: `0.8530 -> 0.8618`
  - `ONT_PAO|paired-pure`: `0.7532 -> 0.8024`

### 4.3 dataset-level 拆解：增益不是全域，但確實跨多個 sample

external validation dataset-level incremental：

| dataset | delta F1 | bootstrap 95% CI | positive_fraction | 判讀 |
| --- | ---: | --- | ---: | --- |
| `HCC1954 Google paired` | `+0.0496` | `[+0.0165, +0.0992]` | `1.000` | 明確支持增益 |
| `COLO829 PAO paired` | `+0.0491` | `[+0.0101, +0.0976]` | `0.998` | 明確支持增益 |
| `H1437 Google paired` | `+0.0227` | `[-0.0153, +0.0660]` | `0.868` | 方向正，但證據未收斂 |
| `HCC1395 DORADO paired` | `+0.0055` | `[+0.0009, +0.0128]` | `0.986` | 小幅但支持增益 |
| `HCC1937 Google paired` | `+0.0001` | `[0.0000, +0.0004]` | `0.606` | 幾乎無差異 |
| `H2009 Google paired` | `-0.0039` | `[-0.0128, +0.0054]` | `0.186` | 輕微負向，未收斂 |

這代表：

1. `methyl+context` 不是只在單一 biological sample 有效。
2. 但也不是「所有 samples 都穩定提升」。
3. 真正正確的說法應該是：
   - 在 `paired-pure multi-bio` 任務中，甲基特徵有**小幅但可重現的整體增益**
   - 在 sample 層級，增益具有明顯異質性

### 4.4 paired-only multi-bio bucket-level diagnostics

以 `context-only -> methyl+context` 的 error delta 來看，最明顯的改善 bucket 包括：

| dataset | bucket | delta error rate |
| --- | --- | ---: |
| `HCC1954` | `TP + Subclone + PassedGating=True` | `-0.4047` |
| `HCC1954` | `FP + Strong + PassedGating=True` | `-0.2564` |
| `COLO829` | `FP + Weak + PassedGating=False` | `-0.1315` |
| `COLO829` | `TP + Weak + PassedGating=False` | `-0.1039` |
| `H1437` | `TP + Strong + PassedGating=True` | `-0.0432` |
| `HCC1395_DORADO paired` | `FP + Weak + PassedGating=False` | `-0.0289` |

最明顯的惡化 bucket 則集中在 `H2009`：

| dataset | bucket | delta error rate |
| --- | --- | ---: |
| `H2009` | `FP + Weak + PassedGating=True` | `+0.3041` |
| `H2009` | `FP + Weak + PassedGating=False` | `+0.1404` |
| `H2009` | `TP + Strong + PassedGating=True` | `+0.0029` |

這說明：

1. `H2009` 的負向表現主要不是全域崩潰，而是集中在 `FP + Weak`。
2. `HCC1954` 與 `COLO829` 的正增益則來自幾個很明確的 verification buckets。
3. 下一輪 diagnostics 應以 `Weak / Subclone / Strong` bucket 為主，而不是只看 sample 整體 F1。

## 5. 綜合解讀

本輪之後，`Phase 1A` 可以明確拆成兩條完全不同的子任務：

### 5.1 paired-pure 路線

結論：

1. `context-only` 已經很強。
2. `methyl+context` 在 `paired-only multi-bio external validation` 上有 `+0.0112 F1` 的穩定整體增益。
3. 這個增益不是假象，因為 bootstrap CI 已完全在零以上。
4. 但增益主要集中在部分 samples，不能直接宣稱「所有樣本都提升」。

### 5.2 to-pure 路線

結論：

1. `to-pure` 仍是主要困難子任務。
2. 在 `sample400` mode-mixed 測試中，負增益主要由 `to-pure` 驅動。
3. 所以目前不應把 paired-pure 的正面結論直接外推到 `to-pure`。

## 6. 目前可下的最保守且正確結論

1. `Phase 1A` read-level ALT/REF classification 已確立為可行研究主線。
2. `context-rich` 模型穩定優於純甲基 baseline，這點已成立。
3. `methyl+context` **不是全域有效**。
4. 但在 `paired-pure multi-bio external validation` 上，`methyl+context` 已出現**小幅但有統計支持的整體增益**。
5. `methyl+context` 對 `to-pure` 目前則沒有增益，甚至會造成明確退步。

## 7. 下一步

1. 將 `paired-pure multi-bio` 正式升級為 `Phase 1A` 的主要驗證軸。
2. `to-pure` 與 `paired-pure` 分開建模，不再共用同一個「是否加甲基」結論。
3. 針對 `H2009 / HCC1937 / H1437` 做 dataset-specific diagnostics，確認為何增益弱或不穩定。
4. 補 `platform/modcall-aware normalization`：
   - `5mCG-only`
   - `5mCG + 5hmCG`
5. 優先追 `H2009 FP + Weak` 與 `to-pure TP REF` 誤判來源，確認是否需要 mode-aware calibration。
6. 若資源允許，將 paired-only multi-bio 從 sampled shard 擴到更大區域數，確認 `+0.0112 F1` 是否仍穩定。
