<!--
建立時間: 2026-03-25 22:05
目標: 針對 Phase 1A（within-tumor alt-support read classification）建立第一版可重跑 benchmark，確認 read-level baseline 是否已具體可行
處理範圍:
  - sample80 與 sample200 的 read-level benchmark
  - discovery holdout 與 external validation
  - pure methyl baseline vs methyl+context logistic 比較
  - 失敗點與下一步
關聯檔案:
  - /big7_disk/liaoyoyo2001/InterSubMod/docs/experiments/in_progress/2026/03/20260325_Phase1_label_schema與harmonization規格_01.md
  - /big7_disk/liaoyoyo2001/InterSubMod/scripts/analysis/run_phase1a_read_classifier_benchmark.py
  - /big7_disk/liaoyoyo2001/InterSubMod/scripts/analysis/export_phase1_manifest_shard.py
  - /big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260325_phase1_manifest_shard_export_sample80
  - /big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260325_phase1_manifest_shard_export_sample200
  - /big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260325_phase1a_read_classifier_benchmark_sample80_v1
  - /big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260325_phase1a_read_classifier_benchmark_sample200_v1
-->

# Phase 1A read classifier benchmark round 1

## 1. 本次研究問題

在已完成 `per-read training/export layer`、`Phase 1A split manifest` 與 `head-to-head baseline table` 之後，現在要回答的問題是：

1. `Phase 1A` 的 ALT/REF read classification 是否已經具備可重跑 benchmark 的條件？
2. 單純依賴 methylation summary 的 baseline，與加入 InterSubMod context 的 read-level model，相差多少？
3. 這個訊號在 `5kHz discovery` 與 `DORADO external validation` 間是否能維持？

## 2. 本次假設

1. `methyl_mean` 單欄位會有弱訊號，但不足以穩定跨平台泛化。
2. `methylation summary + region context + mode / hp / gating` 的 logistic baseline，應明顯優於純甲基 threshold 與純甲基 logistic。
3. `paired-pure` 應明顯比 `to-pure` 容易。

## 3. 成功條件

1. `logistic_methyl_context` 在 discovery holdout 與 external validation 都明顯優於：
   - `majority_ref`
   - `methyl_mean_threshold`
   - `logistic_methyl_only`
2. `sample80` 與 `sample200` 的結論方向一致。
3. 能指出困難子族群，而不是只給整體 F1。

## 4. 失敗條件

1. 若 external validation 完全崩潰，表示目前特徵只是在記住 discovery sample。
2. 若 `sample80` 與 `sample200` 方向相反，表示 benchmark 還不穩定。
3. 若 `paired-pure` 與 `to-pure` 無法分開解讀，則本輪 benchmark 不能作為後續主線依據。

## 5. 主要資料與限制

### 5.1 使用資料

sample80：

- [/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260325_phase1_manifest_shard_export_sample80](/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260325_phase1_manifest_shard_export_sample80)
- `80` 個 regions
- `14,165` 筆 reads

sample200：

- [/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260325_phase1_manifest_shard_export_sample200](/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260325_phase1_manifest_shard_export_sample200)
- `200` 個 regions
- `31,969` 筆 reads

benchmark 輸出：

- sample80：
  [/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260325_phase1a_read_classifier_benchmark_sample80_v1](/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260325_phase1a_read_classifier_benchmark_sample80_v1)
- sample200：
  [/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260325_phase1a_read_classifier_benchmark_sample200_v1](/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260325_phase1a_read_classifier_benchmark_sample200_v1)

### 5.2 benchmark 設定

1. 正樣本定義：`ALT`
2. discovery holdout：
   - 僅使用 `dataset_role=discovery`
   - 以 `region_key` 做 group split，避免 region leakage
3. external validation：
   - 用全部 discovery rows 訓練
   - 用全部 validation rows 測試

### 5.3 比較模型

1. `majority_ref`
2. `methyl_mean_threshold`
3. `logistic_methyl_only`
4. `logistic_context_only`
5. `logistic_methyl_context`

### 5.4 主要限制

1. 目前仍是 selected-region sampled shard，不是 full read universe。
2. `Phase 1B` 仍不可做，因為 normal reads 尚未進入 exporter universe。
3. `sample80 / sample200` 都偏向高品質與可 resolve regions，因此不代表 full manifest 全域難度。

## 6. 本次改動

新增 benchmark 腳本：

- [run_phase1a_read_classifier_benchmark.py](/big7_disk/liaoyoyo2001/InterSubMod/scripts/analysis/run_phase1a_read_classifier_benchmark.py)

功能：

1. 自動切 discovery holdout / external validation
2. 比較 4 種 baseline
3. 輸出：
   - `benchmark_metrics.tsv`
   - `benchmark_group_metrics.tsv`
   - `benchmark_top_coefficients.tsv`

## 7. Benchmark 指標

依研究手冊固定列出：

- `truth_total`
- `calls_total`
- `TP`
- `FP`
- `FN`
- `precision`
- `recall`
- `F1`

此外補充：

- `accuracy`
- `balanced_accuracy`

## 8. 主要結果

### 8.1 sample80 結果

| 模型 | 評估 | truth_total | TP | FP | FN | precision | recall | F1 | accuracy |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| majority_ref | discovery_holdout | 1422 | 0 | 0 | 575 | 0.0000 | 0.0000 | 0.0000 | 0.5956 |
| methyl_mean_threshold | discovery_holdout | 1422 | 532 | 606 | 43 | 0.4675 | 0.9252 | 0.6211 | 0.5436 |
| logistic_methyl_only | discovery_holdout | 1422 | 450 | 331 | 125 | 0.5762 | 0.7826 | 0.6637 | 0.6793 |
| logistic_context_only | discovery_holdout | 1422 | 545 | 191 | 30 | 0.7405 | 0.9478 | 0.8314 | 0.8446 |
| logistic_methyl_context | discovery_holdout | 1422 | 509 | 151 | 66 | 0.7712 | 0.8852 | 0.8243 | 0.8474 |
| majority_ref | external_validation | 7745 | 0 | 0 | 3174 | 0.0000 | 0.0000 | 0.0000 | 0.5902 |
| methyl_mean_threshold | external_validation | 7745 | 2819 | 4143 | 355 | 0.4049 | 0.8882 | 0.5562 | 0.4192 |
| logistic_methyl_only | external_validation | 7745 | 2802 | 4100 | 372 | 0.4060 | 0.8828 | 0.5562 | 0.4226 |
| logistic_context_only | external_validation | 7745 | 2938 | 762 | 236 | 0.7941 | 0.9256 | 0.8548 | 0.8711 |
| logistic_methyl_context | external_validation | 7745 | 3108 | 741 | 66 | 0.8075 | 0.9792 | 0.8851 | 0.8958 |

### 8.2 sample200 結果

| 模型 | 評估 | truth_total | TP | FP | FN | precision | recall | F1 | accuracy |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| majority_ref | discovery_holdout | 2778 | 0 | 0 | 1110 | 0.0000 | 0.0000 | 0.0000 | 0.6004 |
| methyl_mean_threshold | discovery_holdout | 2778 | 1013 | 1277 | 97 | 0.4424 | 0.9126 | 0.5959 | 0.5054 |
| logistic_methyl_only | discovery_holdout | 2778 | 767 | 734 | 343 | 0.5110 | 0.6910 | 0.5875 | 0.6123 |
| logistic_context_only | discovery_holdout | 2778 | 807 | 148 | 303 | 0.8450 | 0.7270 | 0.7816 | 0.8377 |
| logistic_methyl_context | discovery_holdout | 2778 | 826 | 160 | 284 | 0.8377 | 0.7441 | 0.7882 | 0.8402 |
| majority_ref | external_validation | 17159 | 0 | 0 | 7326 | 0.0000 | 0.0000 | 0.0000 | 0.5731 |
| methyl_mean_threshold | external_validation | 17159 | 5896 | 8654 | 1430 | 0.4052 | 0.8048 | 0.5390 | 0.4123 |
| logistic_methyl_only | external_validation | 17159 | 5439 | 7645 | 1887 | 0.4157 | 0.7424 | 0.5330 | 0.4445 |
| logistic_context_only | external_validation | 17159 | 6918 | 1113 | 408 | 0.8614 | 0.9444 | 0.9010 | 0.9114 |
| logistic_methyl_context | external_validation | 17159 | 6882 | 1243 | 444 | 0.8470 | 0.9394 | 0.8908 | 0.9017 |

### 8.3 sample80 vs sample200 穩定性

關鍵比較：

- `logistic_context_only`
  - discovery holdout：`0.8314 -> 0.7816`
  - external validation：`0.8548 -> 0.9010`
- `logistic_methyl_context`
  - discovery holdout：`0.8243 -> 0.7882`
  - external validation：`0.8851 -> 0.8908`

這代表：

1. 結論方向沒有翻盤
2. discovery holdout 在更大樣本下略降，但 context-rich 模型仍遠高於純甲基 baseline
3. external validation 甚至略升，代表跨平台 generalization 不是 sample80 偶然值
4. `methyl + context` 並沒有穩定壓倒 `context-only`

## 9. 分群觀察

以 sample200 的 `logistic_methyl_context` 為主：

| 評估 | harmonization_group | rows_total | precision | recall | F1 | accuracy |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| discovery_holdout | ONT_5kHz\|paired-pure | 1553 | 0.8835 | 0.9504 | 0.9157 | 0.9524 |
| discovery_holdout | ONT_5kHz\|to-pure | 1225 | 0.7985 | 0.6172 | 0.6962 | 0.6980 |
| external_validation | ONT_Dorado\|paired-pure | 8065 | 0.9720 | 0.9980 | 0.9848 | 0.9887 |
| external_validation | ONT_Dorado\|to-pure | 9094 | 0.7723 | 0.8997 | 0.8311 | 0.8245 |

這一輪最清楚的結論不是「平台差異導致崩潰」，而是：

1. `paired-pure` 明顯比 `to-pure` 容易
2. `to-pure` 才是目前真正的困難子任務
3. 在有 context 特徵時，`5kHz -> DORADO` 的跨平台遷移並沒有崩掉

### 9.2 sample200 的 to-pure failure diagnosis（external validation）

正式輸出：

- [/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260325_phase1a_to_pure_failure_diagnosis_v1](/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260325_phase1a_to_pure_failure_diagnosis_v1)

以 `logistic_methyl_context` 為例，`ONT_Dorado|to-pure`：

- 整體 error rate：`17.55%`
- 依 truth_status：
  - `FP`: `25.21%`
  - `TP`: `11.09%`
- 依 VerificationClass：
  - `Subclone`: `99.32%`
  - `Weak`: `55.33%`
  - `Strong`: `16.08%`
  - `Noise`: `1.71%`

這表示 `to-pure` 的主要困難不是所有資料都難，而是集中在：

1. `FP + Subclone`
2. `FP + Weak`
3. 一部分 `FP + Strong`

## 10. 係數觀察與初步解讀

sample200 的 `logistic_context_only` / `logistic_methyl_context` 高權重特徵都包含：

- `hp`
- `mode_to-pure`
- `CramersV`
- `PassedGating_False`

而 `methyl_mean` 只有在 `logistic_methyl_context` 版本中才進入高權重清單。

這表示模型並不是只靠單一甲基均值在判斷，而是同時吃：

1. read-level methylation pattern
2. haplotype / read assignment 訊息
3. region-level context
4. mode 差異

但這裡也要保守：

1. `hp` 的高權重不能直接解讀成生物機制
2. 它也可能反映目前 exporter universe 中的 allele / assignment 結構
3. 因此下一輪應加做 ablation，而不是直接把係數大小當結論

## 11. 初步結論

本輪已能給出 5 個明確結論：

1. `Phase 1A` 已經不只是 schema，而是有可重跑 benchmark 的任務。
2. 純甲基 baseline 不夠：
   - `methyl_mean_threshold`
   - `logistic_methyl_only`
   在 sample80 / sample200 都顯著落後。
3. `context-rich` 模型已顯示明確可行性：
   - sample200 `logistic_context_only` discovery holdout `F1 = 0.7816`
   - sample200 `logistic_context_only` external validation `F1 = 0.9010`
   - sample200 `logistic_methyl_context` discovery holdout `F1 = 0.7882`
   - sample200 `logistic_methyl_context` external validation `F1 = 0.8908`
4. `methyl + context` 相對於 `context-only` 的增益目前尚未被穩定證明。
5. 現階段真正的難點不是 cross-platform，而是 `to-pure`。

## 12. 後續進展

本節原先列出的下一步已在同日完成第二輪：

- [20260325_Phase1A_incremental_test與multi_bio_validation_round2_01.md](/big7_disk/liaoyoyo2001/InterSubMod/docs/experiments/in_progress/2026/03/20260325_Phase1A_incremental_test與multi_bio_validation_round2_01.md)

round 2 已補上：

1. 更嚴格的 `context-only` vs `methyl+context` paired comparison
2. `sample400` 更大 shard 驗證
3. paired-only multi-bio external validation

最新結論以 round 2 為準：

1. `methyl + context` 不是全域有效
2. 在 `paired-pure multi-bio` 有小幅但有統計支持的整體增益
3. 在 `to-pure` 目前沒有增益，甚至會退步
