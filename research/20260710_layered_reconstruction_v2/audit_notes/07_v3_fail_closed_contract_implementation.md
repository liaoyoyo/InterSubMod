<!--
建立時間: 2026-07-11
目標: 紀錄 layered v3 strict source contract、frozen lock 與 adversarial tests 的實作證據。
處理範圍: 新增 v3 preparer/validator/tests/schemas；未修改 v2/live producer/consumer，未讀取 production BAM/sidecar 全量內容。
關聯檔案: InterSubMod/research/20260710_layered_reconstruction_v2/audit_notes/05_fail_closed_wiring_audit.md；InterSubMod/research/20260710_layered_reconstruction_v2/probes/20260711_COLO829_chr1_4386684_4388348_coordinate_join_v1/equivalence_probe.json
任務類型: (B) Comprehensive validation；contract 固定 7 datasets / 6 biological samples / chr1-22。
研究目標: G4 reproducibility；G5 external verifiability。
-->

# Layered v3 fail-closed contract 實作紀錄

## 1. 結論

新增版本化 v3 source manifest 與 v1 frozen lock；shared tiny adversarial suite **17/17 PASS**，post-run receipt normalizer suite **8/8 PASS**。此 contract 不把現行 9 欄 sidecar 偽稱為 payload identity v2，而固定宣告：

- `identity_schema=coordinate_join_v1`
- `assurance=bounded_coordinate_equivalence_not_global_payload_identity`
- `claim_limit=join_method_only_not_per_sample_global_payload_identity`

COLO829 bounded real-data receipt與既有 3-case synthetic receipt僅證明 join 方法；每個 production dataset仍須各自通過 artifact/producer/global-duplicate contract。

## 2. 新增檔案

- `InterSubMod/docs/methodology/_assets/20260627_subclone_4axis_teaching/scripts/prepare_clean_layered_manifest_v3.py`
- `InterSubMod/docs/methodology/_assets/20260627_subclone_4axis_teaching/scripts/validate_layered_v3_inputs.py`
- `InterSubMod/docs/methodology/_assets/20260627_subclone_4axis_teaching/scripts/test_layered_v3_contract.py`
- `InterSubMod/docs/methodology/_assets/20260627_subclone_4axis_teaching/scripts/build_longphase_production_capture_receipt_v3.py`
- `InterSubMod/docs/methodology/_assets/20260627_subclone_4axis_teaching/scripts/test_build_longphase_production_capture_receipt_v3.py`
- `InterSubMod/docs/methodology/_assets/20260627_subclone_4axis_teaching/schemas/layered_input_manifest_v3.schema.json`
- `InterSubMod/docs/methodology/_assets/20260627_subclone_4axis_teaching/schemas/layered_input_lock_v1.schema.json`
- `InterSubMod/docs/methodology/_assets/20260627_subclone_4axis_teaching/schemas/longphase_production_capture_receipt_v1.schema.json`
- `InterSubMod/research/20260710_layered_reconstruction_v2/probes/20260711_COLO829_chr1_4386684_4388348_coordinate_join_v1/layered_v3_contract_test_receipt.json`

## 3. Load-bearing gates

1. Exact canonical scope：7 dataset IDs、6 biological IDs、ordered chr1–22；zero/duplicate/unsafe/unknown field皆拒絕。
2. Production semantics：`production_default_filter`，LongPhase internal filter enabled；threads 12、MAPQ 20、supplementary tagging、somatic VCF output皆 explicit；purity/read-assignment 0.6標記 default；truth flags與`--disableFilter` token-level absent。
3. Artifact identity：sidecar/index/validation與4種角色VCF使用 full SHA-256；BAM/normal/reference使用明示為**非 full content hash**的 `storage_identity_v1`（metadata、first/middle/last chunks、index full hash）。
4. Producer binding：獨立 post-run normalized receipt綁 germline phased VCF+index、normal BAM+BAI、tumor BAM+BAI、LongPhase input ClairS PASS VCF+index、reference+FAI+dict、LongPhase binary/version、runner/capture/validator source、frozen manifest/params/input/output hash manifests及normalizer environment；caller raw VCF只列lossless ledger，不偽稱LongPhase輸入。
5. Tree role：canonical tree backbone固定為LongPhase-S recalibrated FILTER=PASS；ClairS PASS只作LongPhase input，recalibrated-all保留lossless LongPhase ledger。
6. Sidecar binding：重算9欄 gzip sidecar全域 row/key multiplicity，duplicate/conflict必須0；subject同時綁 sidecar/index/native validation/normalized producer receipt、BAM storage identity、producer argv/options/input digest與4角色VCF hashes。
7. Method evidence：real-data與synthetic receipts以 full hash綁入 top-level contract；validator另重算 receipt引用的 bounded sidecar/index與test/consumer source hashes。
8. Atomic lifecycle：錯誤只原子寫 failure report，`valid_lock_written=false`；合法 lock先在candidate完成schema/digest readback後才 publish，publish durability失敗則 quarantine，不留下正式 lock。
9. CN contract：measured branch完整保留source/semantics/coordinate system/unlisted-position semantics/allowed states/overlap policy與artifact hashes；unavailable branch要求non-empty reason且禁止任何measured path。Frozen lock不再只留`availability`。

## 4. 實際驗證

輸入：三個新增 Python檔、兩份新增 schema、全部由 temporary tiny artifacts產生，不連 production資料。

執行命令：

```bash
PYTHONDONTWRITEBYTECODE=1 python3 \
  docs/methodology/_assets/20260627_subclone_4axis_teaching/scripts/test_layered_v3_contract.py \
  --receipt research/20260710_layered_reconstruction_v2/probes/20260711_COLO829_chr1_4386684_4388348_coordinate_join_v1/layered_v3_contract_test_receipt.json
```

實際輸出片段：

```text
Ran 17 tests in 32.860s
OK
```

退出碼：`0`。Receipt machine gate：`pass=true, tests_run=17, failures=0, errors=0`。

Post-run normalizer測試：

```bash
PYTHONDONTWRITEBYTECODE=1 python3 \
  docs/methodology/_assets/20260627_subclone_4axis_teaching/scripts/test_build_longphase_production_capture_receipt_v3.py
```

實際：`Ran 8 tests in 0.375s`、`OK`、exit `0`。涵蓋valid、truth token、input swap、native validation缺欄、output hash mutation、incomplete sample、output不可覆寫，以及HCC1395_DORADO型 symlink BAM inventory。後者要求producer `stat -c` 的logical-path `lstat` metadata與normalizer一致，並另由`storage_identity_v1`綁target realpath/metadata/chunks及BAI full hash。

修正後SHA-256：builder `9f3b7f11393091fec8dc4b0ab2b1b8ca2c5bfd9376443cedc0ae060a1b8ccd3a`；test `64ea2174003689cf411e0b4ae39a97770465731c5ed5c89e7549fb5569259510`。Schema未改動。

另以 exact receipt paths執行 method check，輸出：

```text
ACTUAL_METHOD_RECEIPTS_PASS 2/2
```

## 5. 尚未閉合，不可 launch

1. Active production sidecar validation是native run evidence，不能被改寫；run完成後須以新增normalizer另產7份 `producer_capture_receipt_v3.json`。目前active HCC1395仍缺 `output.sha256`，實測normalizer exit `3`、`E_REQUIRED_ARTIFACT`，且未產receipt。
2. v3 runner現已由獨立版本化entry point接上frozen lock/source bundle/atomic lifecycle與verifier（runner 12/12、lifecycle 13/13、verifier 16/16）；但active producer仍0/7且300秒實機resource gate未跑，因此仍不可launch。詳見`InterSubMod/research/20260710_layered_reconstruction_v2/audit_notes/08_layered_v3_runner_integration.md`。
3. `coordinate_join_v1`不含MAPQ/SEQ/QUAL；方法 receipt僅降低 join implementation風險，不能宣稱7 dataset global payload exactness。對外/L1要求時須升級sidecar payload identity v2並重產。
