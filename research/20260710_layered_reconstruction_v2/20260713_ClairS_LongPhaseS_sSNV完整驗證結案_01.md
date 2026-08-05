<!--
建立時間: 2026-07-13T10:43:16+08:00
目標: 結案 ClairS raw/all 位點進入 LongPhase-S、sSNV read-tag join、演化樹 backbone 與 canonical BAM 不覆寫的全樣本驗證
處理範圍: 7 個資料集、全基因組 ClairS raw/all VCF、LongPhase-S recalibration、chr1-22 regional sSNV reconstruction、canonical/sensitivity 雙 backbone
關聯檔案: InterSubMod/research/20260710_layered_reconstruction_v2/backbone_sensitivity_v3_raw_all_v6.json；InterSubMod/research/20260710_layered_reconstruction_v2/audit_receipts/20260713_canonical_longphase_tagged_bam_post_full_validation_v6.json
-->

# ClairS → LongPhase-S → sSNV 完整驗證結案

## 1. 結論先行

本次為 **Task Type B：Comprehensive validation**，範圍是 7 個資料集與全基因組 ClairS raw/all records，服務 G2、G3、G4、G5。所有指定目標已完成，沒有等待使用者補充才能成立的結論。

1. **ClairS raw/all 位點確實完整交給 LongPhase-S**：輸入 638,259 records，LongPhase-S 輸出 ledger 仍為 638,259 records，沒有位點遺失。
2. **sSNV 分析使用正確且完整的 read tag**：canonical 35 個 parts 共 11,513,224 alignment-group exposures，與 sidecar 11,513,224 筆逐筆 exact match；missing、conflict、extra、malformed、allele conflict 均為 0。
3. **主演化樹輸入應採同一 run 的 LongPhase-S recalibrated `FILTER=PASS` sSNV**：這是 canonical 定義；原始 ClairS `FILTER=PASS` 僅作 sensitivity backbone。
4. **canonical LongPhase-S tagged BAM 沒有被輸出覆寫**：7/7 storage identity 完全相同；producer 使用 named FIFO，沒有建立持久化 BAM，producer root 中 regular BAM 數量為 0。
5. **已合理放大機器資源**：使用 runner 允許上限 `--parallel-samples 4`；主機 48 CPU、約 516 GB 可用記憶體、約 969 GB 可用磁碟均通過 gate。工作負載主要受 BAM I/O 限制，因此不採 5 個以上樣本併行。

## 2. 完整資料流程

| 階段 | 實際資料 | 驗證結果 |
|---|---|---|
| ClairS caller 輸入 | 7 個 raw/all VCF，638,259 records | 全部納入 LongPhase-S input ledger |
| LongPhase-S somatic haplotag | raw/all VCF + tumor BAM，含 supplementary，MAPQ gate 20 | input/output records 都是 638,259 |
| LongPhase-S variant recalibration | `--output-somatic-vcf` | rescue 32,184；remove 17,444；總 records 不變 |
| canonical tree backbone | LongPhase-S recalibrated `FILTER=PASS` | 582,820 records |
| canonical sSNV scope | chr1-22、biallelic sSNV | 469,849 records |
| regional selection | singleton、MAX_SNV 與 scope 分支均留 ledger | retained 194,149 sSNVs |
| read-tag join | primary + supplementary alignment groups | 11,513,224/11,513,224 exact match |
| regional reconstruction | retained sSNV + HP/PS read evidence | 51,815 groups = 51,815 regions |
| sensitivity backbone | 原始 ClairS `FILTER=PASS` | 568,080 records；48,960 regions |

### 2.1 位點守恆與分支

canonical 的 638,259 個 raw/all records 在流程中皆有可追蹤去向：

- LongPhase-S recalibrated `FILTER=PASS`：582,820。
- LongPhase-S filter 排除：55,439。
- canonical PASS 中非 chr1-22 scope：112,971。
- chr1-22 biallelic sSNV：469,849。
- positional singleton：50,432。
- MAX_SNV cap 排除：225,268。
- 最終 retained：194,149。
- read unsupported：0。

「所有位點都被 LongPhase-S 讀取」不等於「所有位點都應進入樹」。前者是 transport/record conservation；後者必須套用同一 run 的 recalibration、染色體與 variant-type scope、regional determinacy 與 read evidence。被排除的 records 仍保留在 site ledger，可追溯原因。

## 3. 為何樹應使用 LongPhase-S sSNV

### 3.1 Canonical 決策

主樹使用 **同一個 LongPhase-S run 產出的 recalibrated `FILTER=PASS` sSNV**。理由如下：

1. LongPhase-S 不只附加 read-level HP/PS；也對 somatic variants 做雙向品質重校正。
2. 本次實測有 32,184 個 non-PASS/LowQual → PASS rescue，也有 17,444 個 PASS → non-PASS/LowQual remove。
3. 若用原始 ClairS PASS 建樹，variant backbone 與後續 LongPhase-S read-tag evidence 的定義不一致。
4. `PS` 是 phase context/QC 欄位，不應單獨當作 clone 或 tree label；樹節點必須由 sSNV mutation state 與 read support 建立。

知識庫依據：`/big8_disk/Knowledge/05_tools/longphase-s.md`，其工具定義包含 somatic haplotagging、read-level HP tags，以及經 `--output-somatic-vcf` 輸出的 PASS/LowQual 雙向 recalibration。

### 3.2 Sensitivity 驗證

| 指標 | LongPhase-S PASS canonical | ClairS PASS sensitivity |
|---|---:|---:|
| tree-input records | 582,820 | 568,080 |
| chr1-22 biallelic sSNV | 469,849 | 455,210 |
| retained sSNV | 194,149 | 182,400 |
| groups/regions | 51,815 | 48,960 |
| read-tag missing/conflict | 0/0 | 0/0 |

跨樣本 comparison verdict 為 `backbone_sensitive`：

- aggregate metric 最大絕對差：2.056768 percentage points。
- 最低 retained-position Jaccard：0.5772569444。
- 最低 primary-unit-key Jaccard：0.4740272187。
- 最低 shared-topology digest concordance：0.9361099879。

因此，aggregate 比例雖只移動約 2.06 pp，實際 site/unit membership 已有實質差異；backbone 不能視為可任意互換。ClairS PASS 結果應保留作 sensitivity，不取代 canonical LongPhase-S PASS。

## 4. Read Tag 完整性

### 4.1 Canonical

- parts：35。
- alignment-group exposures：11,513,224。
- exact sidecar matches：11,513,224。
- primary：11,066,022。
- supplementary：447,202。
- sidecar missing：0。
- sidecar conflict：0。
- sidecar extra：0。
- malformed：0。
- allele conflict：0。

### 4.2 ClairS PASS Sensitivity

- parts：35。
- alignment-group exposures：10,631,171。
- exact sidecar matches：10,631,171。
- primary：10,231,396。
- supplementary：399,775。
- missing、conflict、extra、malformed、allele conflict：全部 0。

這裡驗證的是 read alignment group 與 sidecar 的逐筆 join，而不只是檔案存在。producer 全部使用 `--tagSupplementary`，並以 `collapse_redundant_rows_with_identical_HP_PS` 合併 12,681,029 個 HP/PS 完全相同的冗餘 alignment identities；151,572,508 個 unique coordinate identities 中 conflict 為 0。

## 5. BAM 沒有被覆寫

### 5.1 Producer 行為

- 7 個資料集全部 `persisted_bam=false`。
- transport 全部為 named FIFO。
- producer root regular BAM 數量：0。
- LongPhase-S tagged output 在串流中直接轉成 sidecar，不落地取代 canonical BAM。

### 5.2 Post-validation identity audit

7/7 canonical LongPhase-S tagged BAM 均 `match=true`，且 `differing_fields=[]`。比對欄位包括：

- requested path 與 realpath。
- size、device/inode。
- mtime、ctime。
- first/middle/last 1 MiB SHA-256 chunks。

結論是本次流程沒有輸出並覆蓋原本 tagged BAM。限制是大型 BAM 採 `storage_identity_v1` metadata 加三段 sampled SHA-256，而非對整個 BAM 做 full-content SHA-256；但 inode、size、timestamps 與三段內容雜湊均一致，且 producer 沒有 regular BAM 寫入路徑，兩條獨立證據一致。

範例 canonical BAM：

`/big7_disk/liaoyoyo2001/big7_disk_output/canonical/HCC1395/paired_full/20260314_HCC1395_paired_full_full_complete_matrix/longphase_s/HCC1395_tagged.bam`

## 6. HCC1937 回歸確認

canonical HCC1937 產出 3,612 groups = 3,612 regions。先前缺失的 zero-population key `chr9:61184408-61184992` 已存在：

- 8 個 sSNVs。
- 5 個 full-coverage reads。
- `hp_multiplicity=0`。
- primary lineages 為空。
- `region_determinacy=no_primary_lineage`。
- 沒有捏造演化樹。

這證明「保留區域紀錄」與「強迫產生可解釋樹」已被正確區分。

## 7. 發現、質疑與修正

### 7.1 發現的 contract bug

第一次 ClairS PASS sensitivity v5 在 COLO829 最後 site-ledger contract 失敗。錯誤條件要求 LongPhase-S raw/all input 43,014 records 必須等於 ClairS PASS tree input 38,196 records；兩者本來就代表不同集合。

### 7.2 修正

修改：

`InterSubMod/docs/methodology/_assets/20260627_subclone_4axis_teaching/scripts/build_ssnv_site_ledger.py`

新語意：

- `clairs_PASS_input + clairs_raw_all` 檢查 `caller_raw_PASS_equals_tree_input`。
- 保留 legacy pass-only contract，不改變舊流程行為。

測試：

- 完整 suite：103/103 PASS。
- real COLO829 replay：43,014 rows，7/7 checks true。
- sensitivity v6：7/7 datasets PASS。

失敗的 v5 與原始 v15 execution 都保留，沒有刪除或覆寫，供 audit 重現問題。

## 8. 資源使用判斷

主機有 48 CPUs，執行前 gate 顯示約 516 GB available memory 與約 969 GB free disk。本次以 4 個樣本並行，是 runner/verifier 明確允許的上限；1 到 4 會接受，5 以上會 fail closed。

合理性：每個 sample worker 主要是單核心與大型 BAM 順序/區域讀取，瓶頸是儲存 I/O，不是 CPU 總數。執行期間觀察到其他廣域唯讀掃描使 I/O wait 上升到約 33–35%；停止無關掃描後降到約 0–4%。因此增加到 5 個以上並不會縮短關鍵路徑，反而提高 I/O contention 與結果延遲風險。

結論：**可以且已經更大規模使用機器資源，但合理上限是 4-way bounded sample concurrency，不是無限制把 48 CPU 全部塞滿。**

## 9. 實際命令、輸入與輸出

### 9.1 Canonical + 初次 sensitivity execution

輸入：

- ClairS raw/all VCF 與 canonical tumor BAM：由 7-sample manifest 解析。
- producer root：`/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260711_longphase_s_raw_all_production_sidecars_v2`

命令：

```bash
env PYTHONDONTWRITEBYTECODE=1 \
  INTERSUBMOD_EXEC_ROOT=/big7_disk/liaoyoyo2001/InterSubMod/research/20260710_layered_reconstruction_v2/execution/20260713_raw_all_to_layered_v3_completion_v15 \
  bash scripts/complete_raw_all_layered_v3_validation.sh
```

輸出：

- canonical：`/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260713_layered_reconstruction_v3_raw_all_lps_pass_v5`
- 初次 sensitivity（contract failed，保留）：`/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260713_layered_reconstruction_v3_raw_all_clairs_pass_sensitivity_v5`

實際輸出片段：

```text
canonical: all_pass=true, n_pass=7, n_fail=0
sensitivity v5: E_CHILD_FAILED at COLO829 site-ledger contract
false check: longphase_input_equals_tree_input
```

### 9.2 Fail-closed sensitivity resume

輸入：

- 已成功且 immutable 的 canonical v5。
- 修正後的 ledger builder。
- 同一批 raw/all producer sidecars。

命令：

```bash
env PYTHONDONTWRITEBYTECODE=1 \
  INTERSUBMOD_EXEC_ROOT=/big7_disk/liaoyoyo2001/InterSubMod/research/20260710_layered_reconstruction_v2/execution/20260713_raw_all_to_layered_v3_sensitivity_resume_v16 \
  bash scripts/complete_raw_all_layered_v3_sensitivity_resume.sh
```

輸出：

- sensitivity v6：`/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260713_layered_reconstruction_v3_raw_all_clairs_pass_sensitivity_v6`
- execution root：`InterSubMod/research/20260710_layered_reconstruction_v2/execution/20260713_raw_all_to_layered_v3_sensitivity_resume_v16`

實際輸出片段：

```text
status=SUCCESS
timestamp=2026-07-13T02:40:39+00:00
artifacts_sha256=07a19f51885a4ef16814b247b93ab11ef2a2c537cb237b3d16ffd413d0d5383c
canonical verification: all_pass=true, 7/7
sensitivity verification: all_pass=true, 7/7
comparison aggregate_verdict: backbone_sensitive
BAM post-validation identity: all_match=true, 7/7
```

## 10. 證據索引

- Canonical summary：`/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260713_layered_reconstruction_v3_raw_all_lps_pass_v5/verification_summary.json`
- Sensitivity summary：`/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260713_layered_reconstruction_v3_raw_all_clairs_pass_sensitivity_v6/verification_summary.json`
- Backbone comparison：`InterSubMod/research/20260710_layered_reconstruction_v2/backbone_sensitivity_v3_raw_all_v6.json`
- Comparison table：`InterSubMod/research/20260710_layered_reconstruction_v2/backbone_sensitivity_v3_raw_all_v6.tsv`
- raw/all receipt：`/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260711_longphase_s_raw_all_production_sidecars_v2/raw_all_receipt_closeout.json`
- BAM identity receipt：`InterSubMod/research/20260710_layered_reconstruction_v2/audit_receipts/20260713_canonical_longphase_tagged_bam_post_full_validation_v6.json`
- Fix resolution：`InterSubMod/research/20260710_layered_reconstruction_v2/diagnostics/20260713_clairs_pass_raw_all_ledger_contract_v1/clairs_pass_raw_all_ledger_resolution.json`
- 103-test log：`InterSubMod/research/20260710_layered_reconstruction_v2/diagnostics/20260713_clairs_pass_raw_all_ledger_contract_v1/raw_all_contract_tests.103.log`
- COLO real-data replay：`InterSubMod/research/20260710_layered_reconstruction_v2/diagnostics/20260713_clairs_pass_raw_all_ledger_contract_v1/real_data_replay/COLO829/ssnv_site_ledger_COLO829.summary.json`
- Frozen artifact manifest：`InterSubMod/research/20260710_layered_reconstruction_v2/execution/20260713_raw_all_to_layered_v3_sensitivity_resume_v16/artifacts.sha256`

## 11. 尚需完成或釐清事項

就本次使用者要求的驗證範圍，**沒有未完成目標，也沒有阻塞結論的待釐清問題**。

只有下列三項屬於未來研究政策，不是本次缺漏：

1. 若要把 chrX/chrY 納入正式樹，需要另定性染色體 ploidy、性別與 PAR 規則；目前它們完整保留在 ledger，但 canonical tree scope 是 chr1-22。
2. 若要改 `MAX_SNV` 或 `TIER_R`，需要預先定義新成功標準並重新做跨樣本 sensitivity；不能直接沿用本次數字。
3. 本次輸出是 read-supported regional mutation-state trees；若要宣稱 cellular clones 或病人層級演化史，仍需 purity/CNV、跨區域 clone linking 與獨立 truth evidence，不能只靠 HP/PS 命名。

## 12. 最終判定

本次證據鏈支持以下正式口徑：

> 全部 ClairS raw/all records 已完整提供給 LongPhase-S；主 sSNV regional tree 使用同一 LongPhase-S run recalibrated `FILTER=PASS` 的 chr1-22 biallelic sSNV，並以完整匹配的 primary/supplementary HP/PS read evidence 重建。原始 ClairS PASS 僅作 sensitivity。流程沒有覆寫任何 canonical LongPhase-S tagged BAM。
