<!--
建立時間: 2026-07-12T00:48:00+08:00
目標: 獨立查核 clean layered-v3 最新 terminal 狀態、辨識最新完整可用 topology 資料版本，並界定 HCC1395 技術配對與 cancer/drug annotation 的 claim ceiling
處理範圍: chr1-22、7 dataset rows / 6 biological samples；只讀 producer/canonical/sensitivity/verification 與 historical layered-v2 證據
關聯檔案:
  - InterSubMod/research/20260710_layered_reconstruction_v2/00_INDEX.md
  - InterSubMod/scripts/complete_raw_all_layered_v3_validation.sh
  - InterSubMod/scripts/layered_v3_lifecycle.py
  - InterSubMod/research/20260711_read_group_C_tree_T_topology_report/data/exact_topology_catalog.json
  - InterSubMod/research/20260712_hcc1395_pair_coarse_topology_gene_drug_validation/agents/topology_pair_analysis.md
狀態: point-in-time independent audit
-->

# 最新資料版本與 HCC1395 pair claim ceiling 獨立稽核

## Final-freeze addendum — 2026-07-12 02:23:35 +08:00

> **此段更新下方 00:48 point-in-time 的 live 狀態，但不改變資料來源與 claim ceiling。** Raw-all producer 已於 01:31 完成 **7/7 PASS**，並建立 producer `verification_summary.json` 與 aggregate `_SUCCESS`；然而 receipt closeout 最新重試於 01:57 以 `E_SIDECAR_VALIDATION` 失敗，continuation execution 亦為 `_FAILED`，canonical／sensitivity clean layered-v3 roots仍不存在。因此 historical layered-v2 仍是唯一完整 topology 定量來源，報告仍為 `PARTIAL / SCIENTIFIC NO-GO`。

- Producer：7/7 PASS；`run_status.tsv` SHA-256 `064e66c1...be02`；verification SHA-256 `bc3e0380...53e4`。
- Latest closeout failure：`InterSubMod/research/20260710_layered_reconstruction_v2/diagnostics/20260712_redundant_identity_contract_v1/raw_all_receipt_closeout.failure.attempt2.json`；`E_SIDECAR_VALIDATION: native validation reports duplicates/conflicts/key loss`。
- Continuation：`InterSubMod/research/20260710_layered_reconstruction_v2/execution/20260711_raw_all_to_layered_v3_completion_v1/_FAILED`，exit code 4。
- Canonical／sensitivity：兩個 target roots皆 absent；沒有 clean tree 數字可替換 historical-v2。

> **用 SCQA + PREP：clean layered-v3 截至 2026-07-12 00:48 +08:00 仍只有 producer 6/7 PASS、HCC1954 active；aggregate `_SUCCESS`、canonical tree、sensitivity tree 與 verification receipt 全部尚不存在。因此目前最新完整 topology 數字只能使用 historical layered-v2 engineering snapshot，並維持 `PARTIAL / SCIENTIFIC NO-GO`。HCC1395 與 HCC1395_DORADO 只能支持同一生物樣本跨 basecalling／處理流程的技術重現性，不能當兩個獨立生物樣本，也不能配合 cancer/drug overlap 證明樹為真或方法具臨床效力。**

## 1. 任務分類與研究邊界

- Task type：**(B) Comprehensive validation**。
- 服務目標：**G4** 多樣本一致性／reproducibility；**G5** 可外部驗證的 provenance 與 claim boundary。
- Thread D：否；本輪查核 regional mutation-state topology，甲基化不是主要證據。
- Thread B 撤回方向：否。
- KDE-corrected：不適用於本輪 topology lifecycle 判定。
- VCF caller AF：本輪沒有拿 caller AF 作外部真值；既有多樹排序使用 family-specific read-AF，仍受 CN／coverage confounding。
- 長計算／C++／檔案搬移／NO-GO gate：只讀觀察既有長計算，**未啟停 producer、未修改 canonical output**；科學 claim gate 為 NO-GO。

## 2. 一句 verdict

| 問題 | 判定 | 信心 | 理由 |
|---|---|---:|---|
| clean raw-all producer 是否 7/7 terminal PASS？ | **否；6 PASS / 1 active** | 高 | 6 個 sample `_SUCCESS`；HCC1954 只有 START，host process 仍在執行 |
| producer aggregate `_SUCCESS`／verification 是否存在？ | **否** | 高 | root `_SUCCESS`、`verification_summary.json`、`_RAW_ALL_RECEIPTS_SUCCESS` 全缺 |
| clean canonical LongPhase-S PASS tree 是否存在？ | **否** | 高 | target run root與 v3 manifest皆缺 |
| ClairS PASS sensitivity tree 是否存在？ | **否** | 高 | target run root與 v3 manifest皆缺 |
| 是否已有 canonical-vs-sensitivity comparison／BAM immutability closeout？ | **否** | 高 | comparison、post-BAM receipt、execution `_SUCCESS` 皆缺 |
| clean-v3 是否可取代 historical layered-v2？ | **不可** | 高 | 尚未跨過 producer、receipt、runner、verification 任一全量終點 |
| 目前最新完整可重算 topology 版本？ | **historical layered-v2 engineering snapshot** | 高 | 7/7 internal verification PASS，但 `VALIDATION_SCOPE.md` 明確禁止 comprehensive scientific claim |
| HCC pair 是否可證明方法有效？ | **不可；只支持 technical reproducibility** | 高 | 同一 biological ID、共享上游、無 clone truth；strict tree-set concordance遠低於 coarse agreement |
| cancer/drug overlap 是否補足真值？ | **不可；只能作 exploratory face validity** | 高 | gene/drug annotation不含 region tree truth、drug response或本樣本功能驗證 |

## 3. clean layered-v3 point-in-time 狀態

### 3.1 Exact roots

| Stage | Exact root／target | 截止狀態 |
|---|---|---|
| raw-all producer | `/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260711_longphase_s_raw_all_production_sidecars_v2` | 6/7 sample PASS；HCC1954 active；無 aggregate `_SUCCESS` |
| canonical tree | `/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260711_layered_reconstruction_v3_raw_all_lps_pass_v1` | **root absent** |
| sensitivity tree | `/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260711_layered_reconstruction_v3_raw_all_clairs_pass_sensitivity_v1` | **root absent** |
| continuation execution | `InterSubMod/research/20260710_layered_reconstruction_v2/execution/20260711_raw_all_to_layered_v3_completion_v1` | wait loop active；無 `_SUCCESS`／`_FAILED` |

`InterSubMod/scripts/complete_raw_all_layered_v3_validation.sh` 明確規定：producer aggregate `_SUCCESS` 出現後，才依序 finalise receipts、建立兩份 v3 manifest、執行 canonical tree、執行 sensitivity tree、比較 backbone、驗證 canonical BAM 未變，最後才建立 execution `_SUCCESS`。目前仍停在第一個 wait gate。

### 3.2 Producer 逐 dataset 證據

| Dataset row | biological_id | terminal state | marker mtime (+08) | sidecar rows | rescued / removed | zero-read warnings | sample `_SUCCESS` SHA-256 |
|---|---|---|---|---:|---:|---:|---|
| HCC1395 | HCC1395 | PASS | 2026-07-11 11:33:39 | 40,859,727 | 4,592 / 5,528 | 2 | `8d9c458184952edba99ca016dca4c2f3ba6dd85fae2d7b11e29dfb7ed364a2c6` |
| HCC1395_DORADO | HCC1395 | PASS | 2026-07-11 13:45:37 | 40,033,094 | 5,162 / 4,404 | 0 | `49c8f934c6c58fcc8a7f2acd1eb74606acb47087f479753b66d720e0b9d9ef98` |
| COLO829 | COLO829 | PASS | 2026-07-11 14:49:40 | 8,255,461 | 2,002 / 740 | 0 | `eb2c939d9694e6eb01e863974141c39498625bb7dd586700fa073c03e4d830cf` |
| H1437 | H1437 | PASS | 2026-07-11 17:05:45 | 14,434,968 | 5,168 / 1,091 | 0 | `139cdc54017f77c3f02a5af56c55f429d49fd693966ca412fbaa83c86d7cd89c` |
| H2009 | H2009 | PASS | 2026-07-11 19:57:09 | 21,690,297 | 6,285 / 2,095 | 0 | `8ad0e0a526571b109318e4e3df1af6fea68f22d94f2b8fec5a820d35d71ec770` |
| HCC1937 | HCC1937 | PASS | 2026-07-11 23:27:14 | 23,460,628 | 5,616 / 3,049 | 0 | `2ce69170dabbe443d9974813b8e3fe6f231a8a83b36a48c02f4ba061b2697930` |
| HCC1954 | HCC1954 | **START / active** | START 2026-07-11 23:27:22 | 尚未 terminal | 尚未 terminal | 尚未 terminal | **absent** |

即時 host process 證據（約 00:44 +08）：HCC1954 的 LongPhase-S PID `1677324`、sidecar capture PID `1677321` 仍存活；producer continuation shell PID `1174422` 仍在等待。00:47 producer log已進到 HCC1954 `chr7`，而 `run_status.tsv` 尚未追加 PASS。process PID 只作即時輔證；terminal marker才是 authority。

Frozen producer authority：

- `input_manifest.json` SHA-256 `a62aefd5961bdc9cdd292c44e67b062d58c57c3c39d24b17d2774e7124b67468`
- `code.sha256` file SHA-256 `ba7e9eadd88e7e8b78413ff487da68a8bf1ecfca98b31430b0335e961a33f5a6`
- `params.json` SHA-256 `7973a764a68df9d4f9e3a3e10caacfbe3ec1c4b39354506c9ab8ea290f82501d`
- `run_status.tsv` 截止 23:27:22 的 SHA-256 `e6dfca7ab58f4fc76a77f29b4db2ec87d74b7e18450186dea986c952e2c81f61`

### 3.3 尚缺的 terminal artifacts

截至 cutoff，下列全部 **absent**：

- producer root：`_SUCCESS`、`_FAILED`、`verification_summary.json`、`_RAW_ALL_RECEIPTS_SUCCESS`。
- manifests：`layered_input_manifest_v3_raw_all_lps_pass_v1.json`、`layered_input_manifest_v3_raw_all_clairs_pass_sensitivity_v1.json`。
- canonical/sensitivity run roots與各自 `_SUCCESS`、`verification_summary.json`。
- `backbone_sensitivity_v3_raw_all_v1.json`。
- `20260711_canonical_longphase_tagged_bam_post_full_validation_v1.json`。
- continuation execution root的 `_SUCCESS`／`_FAILED`。

因此，**sample-level producer PASS 不等於 clean tree output PASS**；更不能把兩個 HCC producer PASS 誤當已存在 clean-v3 HCC topology。

## 4. 最新完整可用資料版本

### 4.1 可以計算，但只能作 engineering baseline

目前唯一涵蓋 7 dataset rows、chr1–22 並有完整 region/tree output 的 root：

`/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260710_232501_layered_reconstruction_v2`

其證據為：

- `run_status.tsv` 最後一列：`ALL verify PASS 7/7 datasets`。
- `verification_summary.json`：`dataset_count=7`、`biological_sample_count=6`、`n_pass=7`、`n_fail=0`、`all_pass=true`。
- `verification_summary.json` SHA-256：`b44fe9390cb2bfafce89917c4cc7c476a1bbbcb402a47991030c7832b5d5fc4d`。
- `input_manifest.json` SHA-256：`fd1a9aec8514e602e7ae1407b6e388735488a2679e80a0ac17b41527a59415dc`。

但同 root 的 `VALIDATION_SCOPE.md`（SHA-256 `8a8c8909f00e7489942071f851fd32563ed57e1795ec50904fee73b17c7025a8`）明確指出：6/7 historical LongPhase-S commands 使用 `--truth-bed`，HP tagging並非 genome-wide；只允許 upstream-mismatched engineering baseline／old-vs-clean sensitivity，不可作 canonical comprehensive validation或外部 handoff。

### 4.2 最新 topology derivative 仍繼承上述限制

目前 2026-07-12 pair/coarse-topology 報告的母表仍來自 historical root：

- HCC1395 layered SHA-256 `6b0398c4f32cdc1f8380e675440c2bffb5a83e159768b9111ce9596cf3280b60`。
- HCC1395_DORADO layered SHA-256 `16f38a178952f1ae5ef1c5ac3d2bd9fd1db93c162bb945b2a4b2946d37aa7e7f`。
- `exact_topology_catalog.json` SHA-256 `e0f29a01f3ddffd4a76083f922eed713f0fd1fbf90b2aaf4df497548c4f82c58`。
- `coarse_topology_all_regions.tsv` SHA-256 `90fea7488281d800bab9276bf3cf66a82516c797bd1c40788d2d0ae0719e338c`。

所以應在 HTML 寫成：

> 「以下為目前最新完整可重算的 historical layered-v2 engineering snapshot；clean layered-v3 producer／canonical／sensitivity closeout尚未完成，數值不可升格為 final scientific result。」

## 5. HCC1395 pair 的 provenance 與安全 wording

### 5.1 可直接證明的 provenance

clean producer `input_manifest.json` 對兩列明示：

- HCC1395：`biological_id="HCC1395"`；tumor BAM位於 `.../ONT_5khz_simplex_5mCG_5hmCG/HCC1395.bam`。
- HCC1395_DORADO：`biological_id="HCC1395"`；tumor BAM位於 `.../ONT_Dorado/HCC1395.bam`。
- 兩列皆使用 `ClairS paired`、相同 GRCh38 reference、相同 LongPhase-S/ISM方法家族；DORADO與主列共用 HCC1395 SEQC2 truth。
- repo KB稱 DORADO 為「同細胞株、不同 basecall 方法」及「技術重複」，並明示不可當獨立樣本。

### 5.2 建議使用的精確句子

> **HCC1395 與 HCC1395_DORADO 是同一 biological cell line/specimen 的兩個 basecalling／資料處理 dataset rows；它們可測試 pipeline 對處理版本的 repeatability／robustness，但不是兩個獨立 biological replicates。**

> **在 historical engineering snapshot 中，兩列於大量相同區間呈現高於染色體內 permutation null 的粗拓撲一致性，支持部分 region-level technical reproducibility；這不等同 biological topology accuracy。**

### 5.3 禁止使用的句子

- 「兩個 HCC1395 樣本都證明相同生物克隆樹。」
- 「跨兩個獨立樣本重現。」
- 「因兩者都命中癌症基因／藥物，所以樹為真。」
- 「DORADO replicate驗證方法對不同病人／不同細胞株有效。」
- 「VAF選出的最可能樹就是已確認真樹。」

## 6. 目前 pair 一致性能支持到哪裡

既有獨立 pair analysis 的 exact-coordinate complete pairs為 5,720：

- 五類 raw agreement：69.39%（95% CI 68.18%–70.57%）。
- Cohen's κ：0.497；chromosome-preserving null mean 39.51%，permutation `p=0.00020`。
- aggregate composition不相同：五類 total-variation distance 20.03%；Topo>1為 50.46% vs 63.79%，直系 only為 29.58% vs 13.82%。
- dominant `Topo>1/Topo>1` 同格移除後，剩餘 agreement降為 45.50%。
- exact candidate-tree-set digest agreement：ordered 19.74%；允許 HP1↔HP2 phase swap後 35.21%。

判讀：

1. **有訊號**：粗類別 concordance明顯高於保留 chromosome composition 的 null，不能說完全隨機。
2. **不可互換**：marginal distribution shift與 strict tree-set agreement很低，兩列不是相同結果的複本。
3. **只到 reproducibility，不到 accuracy**：沒有 single-cell／multi-region／known lineage ground truth可判哪一列正確。
4. **粗化會放大表面一致性**：五類會把許多不同 exact trees壓成同一類，因此 69.39%不是 exact-tree accuracy。

## 7. cancer gene／drug annotation 為何不能變成真值

### 7.1 可以使用的角色

- region→GENCODE gene body／promoter：定位與可讀性。
- COSMIC CGC v104：known cancer-gene context／face validity。
- DGIdb local snapshot：exploratory gene–drug interaction annotation；最強展示子集可限制為 approved + antineoplastic。
- 可問：「技術 concordance在已知癌症／可藥物基因區是否與背景不同？」

### 7.2 不能使用的角色

- CGC／DGIdb沒有 HCC1395 的 region tree、clone ancestry或正確 topology label。
- gene–drug interaction不等於該 HCC1395 mutation可治療、藥物有效、或拓撲正確。
- COSMIC resistance records要求 exact variant allele＋drug；目前 coarse region artifact沒有可安全支援的 exact resistance join，gene-only join已被 source audit阻擋。

### 7.3 本輪數據揭示的 circularity

Exact-coordinate matched pair 的 gene/drug flag由座標決定，因此同一 pair兩邊共享 annotation是**配對定義的結果**，不是第二份獨立證據。annotation-stratified agreement只能比較 strata，不能把「兩邊都命中同一基因」當 replicate confirmation。

描述性結果：

- 全體 exact complete agreement 69.39%。
- CGC gene-body present 72.39% vs absent 69.24%，差 +3.15pp。
- approved antineoplastic DGIdb gene-body present 72.25% vs absent 69.14%，差 +3.11pp。

這些差距尚未控制 gene length、region length、SNV count、coverage、CN、VAF、mappability與多重比較；因此只能標 `descriptive / hypothesis-generating`，不能寫成 enrichment被確認，更不能寫成 efficacy proof。

## 8. 五個主要偏差來源與風險

| 風險 | 具體機制 | 若忽略會誤寫成 | 嚴重度 |
|---|---|---|---|
| technical ≠ biological replicate | 同一 biological_id、同細胞株、共用 truth；DORADO是 basecall variant | 跨樣本 biological replication | High |
| shared upstream circularity | 同 caller family、same pipeline、same reference、same truth benchmark；可能同源 raw signal | 獨立方法互相確認 | High |
| region boundary／selection bias | 只分析 exact/overlap matched且 complete的交集會排除 boundary drift與 failed/incomplete regions | 全基因組 accuracy | High |
| CN／coverage／VAF confounding | read-state群、MINREAD、Topo ambiguity及VAF排序都受 depth、allelic CN、multiplicity與basecall影響 | 生物 ancestry真值 | High |
| gene/drug annotation非 ground truth | gene density與著名基因／藥物資料庫 coverage造成廣泛 overlap；interaction≠response | biological/clinical validation | High |

補充：HCC pair exact complete coverage約各 complete regions的 82.4%／84.7%，仍有未配對／不完整區域；任何 headline必須同時報 coverage，不能只報 matched subset agreement。

## 9. Claim ceiling

### 現在可宣稱

1. historical engineering snapshot通過自身 schema/funnel/V1–V7 consistency，適合做 classifier與報告管線的工程檢查。
2. 同一 HCC1395 biological sample 的兩個 processing rows在 coarse topology上有高於 chromosome-preserving null 的 region-level technical concordance。
3. 方法對一部分區域訊號具 processing robustness，但結果存在明顯 basecall/read-support敏感性，不能互換。
4. cancer/drug annotation提供可探索的生物背景與候選區域導覽。

### 現在不可宣稱

1. clean layered-v3已全量完成或可取代 historical v2。
2. 兩個獨立 biological samples重現相同 topology。
3. inferred Steiner tree是已確認 clone tree／cell lineage。
4. cancer-gene或drug overlap驗證 topology accuracy。
5. 方法已被證明「有效」於 biological truth、臨床療效或藥物選擇。

**最高安全結論**：`ENGINEERING REPRODUCIBILITY SUPPORTED / BIOLOGICAL VALIDITY NOT ESTABLISHED / SCIENTIFIC PROOF NO-GO`。

## 10. 升級為方法有效證據所需的最小門檻

1. clean producer exact 7/7 + aggregate `_SUCCESS` + receipts + verification PASS。
2. canonical及sensitivity run各有 frozen launch receipt、`verification_summary.json`、hash-bound `_SUCCESS`，並完成 backbone comparison與 post-run BAM immutability receipt。
3. 用相同五類 classifier重生 clean-v3全 7 rows，完整報 matched與unmatched／incomplete denominators。
4. HCC pair按 CN-neutral/altered、coverage decile、VAF、SNV count、C、Topo1/Topon分層；排除winner read-AF自我驗證。
5. region matching需報 exact、reciprocal overlap與boundary perturbation，並用未配對區作 failure census。
6. cancer/drug analysis採 chromosome／region length／SNV count／gene length／coverage／CN／mappability matched null，預註冊 primary endpoint與多重校正。
7. 至少新增獨立 biological replicates；若要驗證 tree accuracy，需 single-cell DNA、multi-region/longitudinal lineage、orthogonal molecule phasing或帶已知 tree的 simulation/spike-in truth。
8. 藥物 claim需 variant-specific evidence、腫瘤型別／biomarker context及功能或臨床 response；gene-level DGIdb interaction不夠。

## 11. 建議自動化 gate

- producer：exact 7 sample START→PASS＋唯一最後 `ALL verify PASS`；root `_SUCCESS`與closeout hash互綁。
- layered：canonical/sensitivity run root各須有 `_SUCCESS`且其 `verification_sha256`可重算。
- report：資料 source root若不是 clean-v3 verified root，HTML自動加 `HISTORICAL / PARTIAL / SCIENTIFIC NO-GO` ribbon。
- replicate：dataset manifest若 `biological_id`重複，aggregate顯示 dataset-n與bio-n，統計不得當獨立 n。
- annotation：exact-coordinate pairs的annotation一致性不得列為獨立 replicate metric；drug table預設顯示 source/version與 `interaction ≠ efficacy`。

## 12. 主要證據 ledger

| Evidence | Path | 截止觀察／hash |
|---|---|---|
| producer status | `/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260711_longphase_s_raw_all_production_sidecars_v2/run_status.tsv` | 6 PASS + HCC1954 START；SHA `e6dfca...1f61` |
| producer input authority | 同 root `input_manifest.json` | SHA `a62aefd...7468`；HCC pair `biological_id`相同 |
| continuation log | `InterSubMod/research/20260710_layered_reconstruction_v2/execution/20260711_raw_all_to_layered_v3_completion_v1/execution.log` | 00:47仍為 WAIT producer |
| lifecycle source | `InterSubMod/scripts/layered_v3_lifecycle.py` | `_SUCCESS`只在 verification/provenance readback後建立 |
| historical verification | historical root `verification_summary.json` | 7/7 internal PASS；SHA `b44fe939...4d` |
| historical scope override | historical root `VALIDATION_SCOPE.md` | 禁止 comprehensive scientific use；SHA `8a8c890...5a8` |
| topology pair audit | `InterSubMod/research/20260712_hcc1395_pair_coarse_topology_gene_drug_validation/agents/topology_pair_analysis.md` | exact complete n=5,720；69.39%；κ=0.497；strict digest 19.74% |
| gene/drug source inventory | 同 topic `agents/source_inventory.tsv` | GENCODE v46、CGC v104、DGIdb heterogeneous local snapshot；all non-truth |

## 13. 執行命令、輸入、輸出與實際片段

### 輸入

- `InterSubMod/research/20260710_layered_reconstruction_v2/{00_INDEX.md,execution/,audit_receipts/,launch_plans/}`
- `InterSubMod/scripts/{complete_raw_all_layered_v3_validation.sh,layered_v3_lifecycle.py,verify_layered_v3.py}`
- clean producer root與 historical layered-v2 root。
- `InterSubMod/research/autoresearch/evidence_ledger.jsonl` lines 112, 117–118；只用於歷史 provenance，不拿舊 `production_state`當即時狀態。

### 關鍵只讀命令

```bash
find output/synthesis/research_rounds/20260711_longphase_s_raw_all_production_sidecars_v2 \
  -maxdepth 3 -type f \( -name '_SUCCESS' -o -name '_FAILED' -o -name 'filter_transition_audit.json' \)
sha256sum output/synthesis/research_rounds/20260711_longphase_s_raw_all_production_sidecars_v2/samples/*/_SUCCESS
ps -eo pid,ppid,user,stat,etime,lstart,pcpu,pmem,cmd | rg '20260711_longphase_s_raw_all_production_sidecars_v2|complete_raw_all_layered_v3_validation'
find output/synthesis/research_rounds/20260710_232501_layered_reconstruction_v2 -maxdepth 2 -type f
jq '.' output/synthesis/research_rounds/20260710_232501_layered_reconstruction_v2/verification_summary.json
```

### 實際輸出片段

```text
HCC1395 ... PASS
HCC1395_DORADO ... PASS
COLO829 ... PASS
H1437 ... PASS
H2009 ... PASS
HCC1937 ... PASS
HCC1954 production_tagging START

ABSENT producer/_SUCCESS
ABSENT producer/verification_summary.json
ABSENT canonical run root
ABSENT sensitivity run root

historical verification: dataset_count=7; biological_sample_count=6; n_pass=7; n_fail=0; all_pass=true
```

### 輸出

- 本 audit：`InterSubMod/research/20260712_hcc1395_pair_coarse_topology_gene_drug_validation/agents/latest_data_claim_audit.md`
- point-in-time machine-readable status：`InterSubMod/research/20260712_hcc1395_pair_coarse_topology_gene_drug_validation/data/latest_pipeline_status.json`（SHA-256 `de182ef54c4165e32b07b9fcd0f5b1437d0703b7702a1ed7b4c3158dc81d97e4`）

## 14. 最終判定

**資料版本 verdict：PARTIAL / BLOCKED FOR CLEAN SCIENTIFIC RATES。** clean-v3 producer仍 active，canonical/sensitivity/verification尚未開始；不得取代 historical v2。

**方法 claim verdict：SUPPORTED FOR TECHNICAL REPRODUCIBILITY ONLY。** HCC pair顯示部分可重現的 coarse regional signal，但不是 biological replicate，也沒有 topology truth。

**gene/drug claim verdict：EXPLORATORY CONTEXT ONLY。** 可做可追溯的基因圖與候選導覽，不可作 tree ground truth、藥物有效性或方法有效性的證明。
