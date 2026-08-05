<!--
建立時間: 2026-07-11
目標: 為 layered reconstruction v2 設計 fail-closed input validator、frozen manifest/provenance、exact sidecar join、atomic lifecycle 與 adversarial tests。
處理範圍: run_layered_7samples_newbb.sh、validate_layered_v2_inputs.py、verify_layered_v2.py、manifest schema，以及 read-tag sidecar 消費端/producer contract 的必要接線。
關聯檔案: InterSubMod/docs/methodology/_assets/20260627_subclone_4axis_teaching/scripts/；InterSubMod/docs/methodology/_assets/20260627_subclone_4axis_teaching/data/；InterSubMod/research/20260710_layered_reconstruction_v2/audit_notes/01_upstream_core_audit.md；InterSubMod/research/20260710_layered_reconstruction_v2/audit_notes/03_verification_risk_audit.md
任務類型: (B) Comprehensive validation；完整 7 datasets / 6 biological samples；不可 subset。
研究目標: G4 多樣本一致性與 reproducibility；G5 可被外部驗證的業界級貢獻。
限制: 本輪僅 read-only audit 與文件；未修改 pipeline/source/manifest，未啟動新 run，未執行 production BAM/VCF 全量驗證。
-->

# Layered reconstruction v2 fail-closed wiring 與 frozen provenance 稽核

> **用 SCQA + FMEA：完成 fail-closed wiring 設計 — 目前 live runner 不具備安全 launch 條件；至少有 4 個 P0 blocker，必須先修 strict schema、sidecar-to-BAM cryptographic binding、frozen source/environment、以及 atomic lifecycle，才可開始新的正式 7-dataset run（影響：高；信心：高）。**

## 0. 結論先行

### 0.1 Launch verdict

**目前狀態：NO-LAUNCH（技術性 launch gate，並非研究方向 NO-GO）。**

理由不是「程式無法 parse」；4 個 audited scripts 目前均可通過 shell/Python syntax parse。真正問題是 trust boundary 仍然 fail-open：validator、runner、worker、verifier 並未消費同一份不可變 contract，也沒有證據證明 sidecar、consumer BAM、producer command、source code 與 environment 是同一個 run 所宣稱的組合。

### 0.2 必須先解的 blockers

| 等級 | Blocker | 目前直接證據 | 不修的結果 |
|---|---|---|---|
| **P0-1** | Manifest schema 與 runner wiring 不一致 | 現行 manifest 是 `2.0`，7/7 samples 都缺 6 個 runner 已讀取/要求的 sidecar 與 caller 欄位；runner line 24 仍接受 2.0，validator line 123 又把 sidecar 當 optional | Preflight 可先通過，直到 worker line 79 才以空路徑失敗；不是 fail-closed |
| **P0-2** | 0-dataset vacuous PASS | runner 的 jq gate 接受 `dataset_count=0,samples=[]`；validator line 181 與 verifier line 153 都使用 `all([])` | 空 manifest 可報 `all_pass=true`，甚至產生 0/0「成功」摘要 |
| **P0-3** | Sidecar validation 沒綁定 exact bytes/producer/consumer BAM | validator 只 hash sidecar index 與 validation JSON，**沒有 hash sidecar 本體**；`pysam.TabixFile(str(sidecar))` 沒用 manifest 指定 index；validation JSON 也沒 subject hashes | 可把一份 PASS JSON 複製到另一份 sidecar，或提供另一個 index/BAM，而 preflight 無法辨識 |
| **P0-4** | 舊 tagged BAM 不符合 clean haplotag scope | 7/7 BAM 的 LongPhase-S `@PG CL` 都含 `--truth-vcf`，6/7 亦含 `--truth-bed` | 舊 BAM 的 HP/PS 不可作 clean genome-wide evidence；若只靠欄位名稱或 `pass:true`，會混用 truth-restricted tags |
| **P1-1** | Mutable manifest/source 在 validation 後仍被直接讀取 | runner line 38–42、206、210–211 仍讀 `$MANIFEST`；line 121、147、157、164 從 mutable repo 執行 source | 驗證後修改 manifest/source 會改變實際執行內容（TOCTOU） |
| **P1-2** | 只有 hash list，沒有 source bundle/end attestation | runner line 184–189 未 hash runner 本身，也未保存 source；既有 run 的 8 個 code entries 目前 5 個 mismatch | 完成結果無法從 run root 還原 exact executable source，也不能判斷修改發生於 run 中或 run 後 |
| **P1-3** | Existing exact join key 不足以保證 allele payload 相同 | consumer key 只有 QNAME/chrom/start/end/FLAG/CIGAR digest；MAPQ、SEQ、QUAL 未進 key，且 duplicate identical identity 可被 dict 合併 | 兩個 BAM 可在座標/CIGAR 相同時具有不同 allele/base quality，卻誤認為 exact join |
| **P1-4** | Lifecycle 非 atomic、無 trap/heartbeat | final run root 先建立；status 由多 process 直接 append；JSON/TSV 多為 direct write；`xargs -P` 無 tracked child registry；無 `_SUCCESS` | 中斷、部分 write、worker orphan、stale run 都可能留下外觀近似完成的目錄 |
| **P1-5** | Verifier 信任 producer boolean 且仍讀 mutable manifest | verifier line 96–99 只相信 output 內的 `check_exact_sidecar_join` / site-ledger `pass`；line 144–146 讀原 manifest | producer 自己寫錯 boolean、sample dir 被替換、manifest 被換掉時，final verifier 無獨立證據 |

### 0.3 最小安全修補範圍

新正式 run 前的**最小必要範圍**不是只改一個 jq expression，而是以下完整 chain：

1. 新增 strict JSON Schemas：source manifest v3、frozen lock v1、run receipt/state v1。
2. 修改 `validate_layered_v2_inputs.py`：strict schema、artifact/producer/sidecar identity、atomic frozen lock；任何不確定都 non-zero。
3. 修改 `run_layered_7samples_newbb.sh`：先 preflight、後 publish run root；worker 只讀 frozen lock/source bundle；加入 trap/state/heartbeat/atomic output gates。
4. 修改 `verify_layered_v2.py`：只讀 frozen lock + launch receipt；獨立驗證 provenance、dataset set、output hashes/state，拒絕空集合。
5. 修改 sidecar producer/capture/validator **或拒絕現行 sidecar schema**：必須產生 alignment identity v2、subject hashes 與 global stream digest。
6. 修改 `sm_linkage_genomewide.py`：explicit sidecar index、identity v2、duplicate/malformed/missing fail hard、不得 fallback；輸出 lock/sidecar digest。
7. 加入 tiny synthetic fixtures 與 adversarial tests；測試通過前不得 full launch。

若第 5–6 點不在本次可改範圍，唯一 fail-closed 行為是 validator 對現行 sidecar schema 回傳 `E_SIDECAR_CONTRACT_UNSUPPORTED`，**不可把目前的 9-column sidecar 宣稱為 cryptographically exact join**。

---

## 1. 稽核邊界、證據語彙與 snapshot

### 1.1 證據標記

- `[O-L4]`：本輪以 exact path + hash/header/程式行號直接觀察。
- `[F-L3]`：程式碼的可重現控制流事實。
- `[I-L2]`：由多個直接證據導出的推論，另列限制。
- `[U-L1]`：本輪沒有 exact artifact/path，因此不得當結論。

### 1.2 Audited source snapshot

以下 line references **只對這組 SHA-256 有效**；live files 在本輪前曾有 concurrent edits，因此後續實作前必須先重新 hash/diff。

| Alias | Repo path | Lines | SHA-256 | mtime (Asia/Taipei) |
|---|---|---:|---|---|
| `RUNNER` | `InterSubMod/docs/methodology/_assets/20260627_subclone_4axis_teaching/scripts/run_layered_7samples_newbb.sh` | 217 | `d8a8037f71f9185333ab9e4c311fdb5ee532099b616815c4364619478e63ef2a` | 2026-07-11 04:06:34 |
| `PREFLIGHT` | `InterSubMod/docs/methodology/_assets/20260627_subclone_4axis_teaching/scripts/validate_layered_v2_inputs.py` | 194 | `6e7ee45b0a5861761f678d67050649c9e4f139dfa955e80725caee34ae000500` | 2026-07-11 04:06:34 |
| `VERIFY` | `InterSubMod/docs/methodology/_assets/20260627_subclone_4axis_teaching/scripts/verify_layered_v2.py` | 178 | `1e2cfc4326ba45df21805618e370498c556555c7a9d88fb74a6a2922ccb94f86` | 2026-07-11 04:04:59 |
| `CONSUMER` | `InterSubMod/docs/methodology/_assets/20260627_subclone_4axis_teaching/scripts/sm_linkage_genomewide.py` | 371 | `b01be8e3bba276012fa28bb798062353c3ee00978d58081ebf5577d58555557a` | 2026-07-11 03:58:39 |

`bash -n RUNNER` exit 0；`ast.parse(PREFLIGHT, VERIFY, CONSUMER)` 為 3/3。這只代表 syntax 可載入，**不代表 contract、provenance 或科學語意正確**。[O-L4]

> **Post-write drift notice（2026-07-11 04:23:50 +08:00）：**外部 hardening session 在本文件完成後再次修改 live `RUNNER` 與 `PREFLIGHT`；readback hash/lines 分別變成 `29a3a7602cb3561dc3f9d47170d86ccef9ec24e181274fafc178da8bdde9d2a5` / 233 lines，以及 `95453fd4412d1a60113909c2228e329e04bb69fa26f8b71cdf329a2fa37ce3e6` / 217 lines。本文件沒有修改兩個 live files，也沒有把新版本混入舊稽核；§4/§6 的 line-level patch points仍只適用上表 04:06 snapshot。實作者必須對新 SHA重新逐項確認 blocker是否已消除，不能以 line number機械套 patch。

### 1.3 Current manifest mismatch

Current manifest：

`InterSubMod/docs/methodology/_assets/20260627_subclone_4axis_teaching/data/layered_v2_input_manifest.json`

- SHA-256：`fd1a9aec8514e602e7ae1407b6e388735488a2679e80a0ac17b41527a59415dc`
- `schema_version=2.0`、`dataset_count=7`、`biological_sample_count=6`
- 7/7 sample objects 都沒有：`read_tag_sidecar`、`read_tag_sidecar_index`、`read_tag_validation`、`caller_raw_vcf`、`longphase_recalibrated_all_vcf`、`somatic_vcf_role`。[O-L4]

然而 `RUNNER` line 61–65、79–82 將前五個 path/validation 當作必須存在；line 24 卻接受 2.0；`PREFLIGHT` line 123–137 將 sidecar 當 optional。[F-L3] 因此 schema、preflight、worker 三者目前不是同一 contract。

### 1.4 Existing completed run 不能當 frozen-source 證明

已存在的：

`/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260710_232501_layered_reconstruction_v2/`

其 `run_status.tsv` 最後為 `ALL verify PASS 7/7 datasets`，`verification_summary.json` 為 `all_pass=true,n_pass=7,n_fail=0`。[O-L4] 這證明該 run 的**內部 verifier checks**通過，但不能補足 provenance：

- `sha256sum -c code.sha256` 對 current tree 為 3 OK / 5 FAILED。
- FAILED：`sm_linkage_genomewide.py`、`sm_multilocus_combinations.py`、`layered_tree_reconstruction.py`、`build_region_view.py`、`verify_layered_v2.py`。
- run root 只有 hash list，沒有 source bundle；hash 指向 repo 的 mutable absolute paths。
- 沒有 end-of-run source digest，不能由 run artifact 判定 edits 是「執行中」還是「執行後」發生。[I-L2]

因此這個既有 7/7 PASS 可作 historical observation，**不可宣稱可由 run root exact reproduce**。

### 1.5 Historical BAM producer scope

對 current manifest 的 7 個 `tumor_bam` 只讀 BAM header：

| Dataset | `SO:coordinate` | `@RG SM` | LongPhase-S `@PG` | `--truth-vcf` | `--truth-bed` | `-q 20` | `--tagSupplementary` |
|---|---|---|---|---|---|---|---|
| HCC1395 | true | empty | present | true | true | true | true |
| HCC1395_DORADO | true | empty | present | true | true | true | true |
| COLO829 | true | empty | present | true | false | true | true |
| H1437 | true | empty | present | true | true | true | true |
| H2009 | true | empty | present | true | true | true | true |
| HCC1937 | true | empty | present | true | true | true | true |
| HCC1954 | true | empty | present | true | true | true | true |

結論：這些 BAM 可以在明確政策下作為 **alignment payload**，但其 embedded HP/PS 必須標記 `ignored`；不可當 clean sidecar 的 tag authority。[O-L4] `@RG SM` 空白本身不等於 sample 錯誤，但 schema 必須明確選擇：要求非空一致，或以 reviewed manifest + artifact hash 作 sample binding；不得默默推論。

### 1.6 本輪未驗證事項

- `[U-L1]` 外部 session 正在產生的新 sidecar run root 未提供 exact path；本輪沒有掃描磁碟尋找它，也沒有稽核其實際 output hashes/counts。
- `[U-L1]` 本輪沒有重新執行全量 validator、BAM content scan 或 pipeline；以下 adversarial tests 是**設計規格，不是已執行結果**。
- `[U-L1]` 尚未決定正式政策要對大型 BAM 使用 `full_sha256` 或 lower-assurance `storage_identity_v1`；兩者不得混稱。

---

## 2. Fail-closed trust model

### 2.1 唯一允許的 trust chain

```text
reviewed source manifest v3
        │ strict schema + semantic cross-field validation
        ▼
observed artifact identities + producer command/input/scope validation
        │ exact sidecar ↔ alignment-payload global join proof
        ▼
atomic frozen input lock v1
        ├── source bundle digest + stored source bytes
        ├── environment lock digest
        └── canonical lock digest
        │ launch receipt binds all three digests
        ▼
workers read ONLY frozen lock + extracted source bundle
        │ atomic sample outputs + state events + heartbeat
        ▼
post-input identity recheck + independent final verifier
        │ output hashes + provenance graph + terminal state
        ▼
_SUCCESS
```

任何 missing、unknown、ambiguous、mixed mode、digest mismatch、parser exception、stale index、duplicate identity、unexpected environment variable 都是 error；不得轉成 empty dict、skip、warning-only 或 default neutral。

### 2.2 Authority 分離

| Data | 唯一 authority | 必須忽略/拒絕的替代來源 |
|---|---|---|
| Alignment bases/qualities/coordinates | `alignment_payload_bam` frozen identity | 不能由 sidecar補 bases；BAM 在 run 中改變即 fail |
| HP/PS | `external_exact_sidecar` | consumer BAM embedded HP/PS 必須 `ignored`; sidecar missing 不得 fallback |
| Somatic backbone | ClairS paired PASS artifact 明確 role | 不得把 raw/recalibrated VCF 靜默替代 backbone |
| CN | manifest 宣告的 measured source 或 explicit unavailable | missing CN 不得當 neutral；unlisted position semantics 必須明示 |
| Source | frozen source bundle | mutable working tree 不再是 executable authority |
| Manifest | frozen lock | source manifest 只作 provenance input；worker 不重讀它 |
| Success | final verifier + `_SUCCESS` | 目錄存在、單一 `pass:true`、status 最後一行都不等於 success |

---

## 3. Manifest / frozen lock schema

### 3.1 Schema 分層

建議不要再把 source intent、observed identity、run state 混在同一個 `schema_version=2.x`：

| Artifact | 建議 schema name/version | 誰寫 | 誰讀 |
|---|---|---|---|
| Source manifest | `intersubmod.layered_input_manifest` / `3.0.0` | reviewed builder/human | strict validator only |
| Frozen lock | `intersubmod.layered_input_lock` / `1.0.0` | validator atomic write | runner/workers/verifier |
| Environment lock | `intersubmod.layered_environment_lock` / `1.0.0` | validator/provenance helper | runner/verifier |
| Source bundle manifest | `intersubmod.layered_source_bundle` / `1.0.0` | bundler | runner/verifier |
| Launch receipt/state | `intersubmod.layered_run_receipt` / `1.0.0` | runner | verifier |
| Scientific output | existing layered output `2.0` | workers | verifier/report |

每個 JSON Schema 必須：`additionalProperties:false`、所有 load-bearing fields 放入 `required`、enum 精確、integer range 精確、URI/path 不接受 control character；version mismatch 不做 best-effort migration。

### 3.2 Source manifest top-level required fields

| Field | Type / constraint | 判斷目的 |
|---|---|---|
| `$schema` | repo-relative schema path/string const | 防止用錯 validator/schema |
| `schema_name` | const `intersubmod.layered_input_manifest` | 不只靠裸版本號 |
| `schema_version` | const `3.0.0` | sidecar mode 是 breaking contract，不接受 2.0/2.1 |
| `manifest_id` | non-empty, immutable identifier | provenance node ID |
| `created_at_utc` | RFC3339 UTC | audit ordering；不是 content identity |
| `dataset_count` | const 7, minimum 1 | 防空集合與漏 sample |
| `biological_sample_count` | const 6, minimum 1 | 防生物樣本計數錯 |
| `analysis_contract.scope` | const `whole_autosomes_chr1_22` + exact contig array | 禁止 subset 混入 comprehensive run |
| `analysis_contract.read_tag_mode` | const `external_exact_sidecar` | 禁止 mixed embedded/sidecar mode |
| `analysis_contract.embedded_tag_policy` | const `ignore` | 歷史 BAM tag 不得成 authority |
| `analysis_contract.require_exact_join` | const true | schema-level gate |
| `analysis_contract.launch_profile` | const `comprehensive_validation` | 綁定 `VERIFY_EVERY=1`,`ANALYSIS_TREE_CAP=0` |
| `samples` | array, `minItems=maxItems=7` | 每個 dataset 一筆；sample/biological cardinality另作 semantic check |

### 3.3 Per-sample required fields

```json
{
  "sample": "HCC1395",
  "biological_id": "HCC1395",
  "platform": "ONT_PROMETHION",
  "replicate_role": "canonical_or_platform_replica",
  "alignment_payload": {
    "role": "tumor_alignment_payload",
    "path": "/absolute/path/tumor.bam",
    "index": {"path": "/absolute/path/tumor.bam.bai"},
    "embedded_hp_ps_policy": "ignore",
    "identity_policy": "full_sha256",
    "expected_identity": {}
  },
  "somatic": {
    "backbone_role": "clairs_paired_filter_pass",
    "caller_pass_vcf": {"path": "...", "index": {"path": "..."}, "expected_identity": {}},
    "caller_raw_vcf": {"path": "...", "index": {"path": "..."}, "expected_identity": {}},
    "longphase_recalibrated_all_vcf": {"path": "...", "index": {"path": "..."}, "expected_identity": {}}
  },
  "read_tags": {
    "authority": "external_exact_sidecar",
    "sidecar": {"path": "...bgz", "expected_identity": {}},
    "index": {"path": "...bgz.tbi", "expected_identity": {}},
    "validation": {"path": "...json", "expected_identity": {}},
    "identity_schema": "alignment_identity_v2",
    "require_unique_identity": true,
    "require_global_payload_match": true,
    "fallback_policy": "forbidden",
    "producer": {}
  },
  "copy_number": {
    "availability": "measured",
    "source_id": "reviewed_source_name",
    "segments": {"path": "...", "expected_identity": {}},
    "integer_gain": {"path": "...", "expected_identity": {}},
    "integer_loss": {"path": "...", "expected_identity": {}},
    "coordinate_system": "0_based_half_open",
    "allowed_states": ["gain", "loss", "loh", "neutral"],
    "unlisted_position_semantics": "neutral",
    "overlap_policy": "forbid",
    "integration": {"path": "...json", "expected_identity": {}}
  }
}
```

若 CN unavailable，必須使用另一個 `oneOf` branch：`availability:"unavailable"`、`reason` 必填，所有 CN path 必須 absent，`unlisted_position_semantics:"unavailable"`；不得 `cn_source=unavailable` 但仍給 measured files，或 path 空白又默認 neutral。

### 3.4 Producer command / clean haplotag scope

`read_tags.producer` 至少要求：

| Field | Required rule |
|---|---|
| `program` / `version` | exact allowlist，例如 `longphase-s` 與 reviewed version |
| `pg_id` / `pp_chain_sha256` | 綁定 BAM `@PG` graph；若 producer output 不保留 BAM，則由 signed execution record提供 |
| `command_argv` | JSON string array；不得只存 shell-rendered string |
| `command_argv_sha256` | 對 canonical JSON array 的 SHA-256 |
| `effective_options` | 每個 option 含 value 與 origin (`explicit`/`default`)；至少 MAPQ、supplementary、threshold、output mode |
| `haplotag_scope` | const `genome_wide_clean` |
| `truth_filtering` | const `forbidden` |
| `forbidden_options_absent` | 必含 `--truth-vcf`,`--truth-bed` 且 result true |
| `input_bindings` | germline phased VCF、normal BAM、raw tumor BAM、somatic VCF、reference 各自 artifact identity |
| `execution_log` | full hash；parsed truth fields只是附加證據，不能取代 argv |
| `source_bundle_sha256` / `environment_lock_sha256` | producer 本身也可重建 |

**判斷規則：**argv token exact match 優先；不以 substring 猜參數。`--truth-vcf=`、縮寫、response file、quoted value 都要由同一 parser canonicalize。發現 forbidden option，即使 path 空白也 fail；clean contract 不依賴「log 看起來沒有移除」這種弱推論。

### 3.5 Full identity 與 storage identity

所有 small artifacts（manifest、VCF、index、CN、integration、sidecar、validation、source、environment）一律 `full_sha256`。大型 BAM 提供兩個**明確不同 assurance** 的 policy：

#### `full_sha256`（正式/對外交付預設）

- `realpath`, `lstat` symlink metadata（若允許 symlink）
- `size_bytes`, `sha256`
- BAM header SHA、reference dictionary SHA
- BAI full SHA
- `samtools quickcheck -v` exit 0、EOF block valid、coordinate sort/header/index reference一致
- pre-launch 與 post-run full hash相同

#### `storage_identity_v1`（只有顯式核准才可用；不得標成 content hash）

- `realpath`, `st_dev`, `st_ino`, `size_bytes`, `mtime_ns`, `ctime_ns`
- first/middle/last fixed chunk SHA-256（offset/length寫入 lock）
- BAM header SHA、reference dictionary SHA、BAI full SHA
- `quickcheck/index/sort` checks
- filesystem/storage contract ID
- pre-launch、每 worker dispatch 前、post-run全數一致

storage identity 仍無法排除精心維持 metadata/chunks 的中段 mutation；因此 comprehensive/final policy 應要求 full hash，或在 report 明確降級 assurance。validator 不可在 full hash失敗/太慢時自動 fallback。

### 3.6 Frozen lock contents

Frozen lock 必須包含：

- source manifest 的原始 byte hash與 canonical JSON hash；最好嵌入 canonical snapshot。
- validator schema/version/source SHA、exact `command_argv`。
- 每個 artifact 的 `expected_identity`、`observed_identity`、comparison result。
- producer command/input/scope 的 observed evidence與 error-free結果。
- sidecar global identity-stream proof（見 §5）。
- source bundle manifest/hash、environment lock/hash。
- generated timestamp、host ID、user/uid、working tree commit/tree/dirty diff bundle identity。
- canonical digest algorithm：UTF-8、sorted keys、compact separators、禁止 NaN、計算時排除 self-digest field。
- `lock_sha256`，寫入 launch receipt與每個 sample output。

Frozen lock 使用 same-directory temp、flush、`fsync(file)`、`os.replace()`、`fsync(parent dir)`；權限在 publish 後 read-only。任何 write failure 不得留下正式檔名。

---

## 4. `validate_layered_v2_inputs.py` fail-closed 設計

### 4.1 現行判斷的精確缺口

| Current lines | 現況 | Fail-closed replacement |
|---|---|---|
| 13–15 | 只有 autosomes/CN label constants | 增加 schema name/version、allowed samples/contigs、identity policy、stable error code constants |
| 17–27 | full SHA helper只有一般 file | 提供 canonical JSON、atomic write、lstat/realpath、full/storage identity、command argv digest |
| 30–69 | 非 biallelic/malformed records line 42–43 被 silent `continue`；只統計 SNV FILTER | 所有 record 分類並守恆；malformed/unknown contig/index mismatch直接 error；raw/PASS/recal record-key contract明確 |
| 72–86 | BAM 只記 size/mtime/header/ref-dict/BAI SHA | 加 full/storage identity、quickcheck、EOF、sort/index、`@PG` graph、RG sample policy、pre/post recheck |
| 89–100 | CN 少欄位直接 skip，未檢 overlap/sorted/coord/unlisted semantics | strict parse；空/壞 interval、負座標、`start>=end`、overlap、unknown state、source/path矛盾都 fail |
| 103–107 | CLI 只有 manifest/output | 加 `--schema --freeze --identity-policy --source-bundle --environment-lock --require-clean-sidecar` |
| 108–110 | 直接 `manifest["samples"]`；無 schema/cardinality/unique | 先 strict JSON Schema，再 semantic cross-field checks；必須 exact 7/6，unique safe IDs |
| 112–122 | 只檢幾個 path存在 | 所有 artifact/index/integration必須符合 role、type與 expected identity |
| 123–137 | sidecar optional | v3 comprehensive branch中 sidecar/validation/index必填；mixed mode fail |
| 138–155 | validation只看 `pass`/truth boolean；sidecar不 hash；指定 index沒被使用；chr22空不會 fail | full subject hash binding、explicit index、all 22 contigs、global identity stream、producer argv/input binding；任何空/exception fail |
| 156–161 | VCF/CN shallow check | 完整 invariants + conservation；不接受 silent skip |
| 162–181 | 由 results 推 counts；`all([])` | counts 必須先 match manifest contract且 `len==7`；zero/duplicate/unknown sample fail |
| 184–185 | direct JSON write | atomic frozen lock write；若任何 fail只寫 separate failure report，不產生 valid lock |

### 4.2 Validator algorithm（固定順序）

1. **Parse bytes**：UTF-8 strict、duplicate JSON keys fail、NaN/Infinity fail。
2. **Schema gate**：schema name/version exact；`additionalProperties:false`；7/6 cardinality。
3. **Semantic gate**：sample ID unique且 path-safe；biological IDs six；HCC1395 platform pair關係明示；scope exact chr1–22；profile是 comprehensive。
4. **Path gate**：`resolve(strict=True)`；記 requested path/realpath/lstat；symlink policy明示；file type regular；index path exact。
5. **Artifact identity gate**：對每個 artifact計 expected vs observed；不存在「只有存在就通過」。
6. **BAM gate**：quickcheck/header/sort/ref dictionary/index；producer `@PG`/input binding；embedded tag policy。
7. **VCF gate**：header/caller/paired role/index/record parse/contig/filter/key conservation；禁止 silent record drop。
8. **CN/integration gate**：availability branch、coord/state/overlap/coverage semantics、cross-file relationship。
9. **Sidecar producer gate**：clean command argv、forbidden flags、producer inputs、sidecar/index/validation hashes。
10. **Exact global join gate**：alignment identity v2/count/digest與 consumer BAM exact match。
11. **Source/environment gate**：source bundle與 environment lock完整且 hash一致；未知 `SM_*` fail。
12. **Re-observe critical identities**：避免 validation期間 TOCTOU。
13. **Atomic freeze**：只有 1–12 全 PASS才建立 frozen lock；failure report與 lock檔名分離。

### 4.3 Stable error taxonomy

| Exit | Error family | Example codes |
|---:|---|---|
| 2 | schema/CLI | `E_SCHEMA_VERSION`, `E_SCHEMA_EMPTY`, `E_SAMPLE_DUPLICATE`, `E_UNKNOWN_FIELD` |
| 3 | artifact identity | `E_HASH_MISMATCH`, `E_INDEX_MISMATCH`, `E_BAM_QUICKCHECK`, `E_PATH_RETARGETED` |
| 4 | producer/scope | `E_PRODUCER_AMBIGUOUS`, `E_TRUTH_FLAG_PRESENT`, `E_INPUT_BINDING_MISMATCH` |
| 5 | sidecar exact join | `E_SIDECAR_SUBJECT_MISMATCH`, `E_JOIN_MISSING`, `E_JOIN_EXTRA`, `E_JOIN_DUPLICATE`, `E_JOIN_PAYLOAD_MISMATCH` |
| 6 | source/environment | `E_SOURCE_BUNDLE_MISMATCH`, `E_ENVIRONMENT_MISMATCH`, `E_UNKNOWN_SM_ENV` |
| 7 | lifecycle/output | `E_STATE_INVALID`, `E_HEARTBEAT_STALE`, `E_OUTPUT_HASH_MISMATCH`, `E_DATASET_SET_MISMATCH` |
| 70 | unexpected internal error | stack trace另存 failure report；不得轉成普通 FAIL/PASS混合 |

Machine-readable failure 每筆包含 `code,stage,sample,artifact,expected,observed,message`；human message 不能取代 stable code。

---

## 5. Exact sidecar join contract

### 5.1 現有 consumer 的優點與限制

正面：`CONSUMER` line 162–175 在 sidecar configured 時沒有 fallback 到 BAM HP/PS；`RUNNER` line 133–138 又要求 analyzed alignment exposures 全部 exact matched、missing/conflict為 0。[F-L3]

但這仍不足：

- line 91–94 的 key 不含 MAPQ、query sequence、query qualities；實際 allele call在 line 159–160讀 sequence，filter在 line 156讀 MAPQ，pileup也受 qualities影響。
- sidecar row雖含 MAPQ（line 114），consumer把它讀成 `_mapq` 後忽略。
- line 103 未傳指定 index；manifest 的 `read_tag_sidecar_index` 可能不是實際被 tabix使用的 index。
- line 108–109 將 fetch error轉成空 dict；只在 expensive worker後由 count gate發現。
- identical duplicate keys不算 conflict，dict會 collapse；multiplicity不守恆。
- exact check只覆蓋 sSNV group exposure，不是 consumer BAM 與 sidecar的全域 alignment payload。

### 5.2 Alignment identity v2

建議 sidecar v2 每個 mapped alignment的 authority key為：

```text
sample_id
+ QNAME_BYTES_SHA256
+ RNAME
+ START0
+ END0
+ FLAG
+ MAPQ
+ BLAKE2b-128(CIGAR string bytes)
+ BLAKE2b-128(SEQ bytes; "*" canonicalized)
+ BLAKE2b-128(QUAL bytes; missing canonicalized)
```

並要求 key **globally unique**。若真實資料允許完全重複 alignment，則 schema 必須加入 deterministic `occurrence_index` 並對 multiplicity做 multiset digest；不能讓 Python dict自行 collapse。

選擇 SEQ/QUAL不是過度設計：本流程的 REF/ALT與 base-quality判斷直接依賴它們。若只匹配 coordinates/CIGAR，所謂 exact join並不保證分析讀到同一 allele payload。

### 5.3 Global binding proof

Sidecar capture summary/validation至少輸出並互相綁定：

- `sidecar.sha256`, `sidecar_index.sha256`, `validation_subject_schema`。
- `alignment_payload_bam.expected/observed identity`。
- `mapped_alignment_count`, `identity_unique_count`, `duplicate_count=0`。
- 對 BAM mapped alignment stream與 sidecar identity columns分別 canonicalize後的 `identity_stream_sha256`；兩者 count/digest exact equal。
- `missing=0, extra=0, duplicate=0, malformed=0, conflict=0`。
- producer input VCF/reference/BAM hashes、command argv hash、source/env hash。
- validation JSON自身不得 self-assert未被驗證的 subject；downstream validator必須重算 subject hashes。

### 5.4 Consumer runtime rules

1. frozen lock指定 `sidecar_path` 與 **exact `sidecar_index_path`**；`pysam.TabixFile(path,index=index_path)`。
2. consumer startup先驗 sidecar/index SHA與 lock digest；不是每個 worker自行相信 path。
3. malformed row、unknown HP、duplicate key、fetch exception立刻 raise；不回 `{}`。
4. sidecar configured時 embedded HP/PS永遠 ignore；missing sidecar tag或 key立刻 fail，不 fallback。
5. 每個 part output寫 `input_lock_sha256,sidecar_sha256,sidecar_index_sha256,identity_schema,join counts`。
6. final verifier要求 5 parts的 digests一致、exposures守恆，且等於 frozen global proof的 subject。

### 5.5 Sidecar producer 的不可省略 patch points

Current capture/validation paths：

- `capture_longphase_tagged_bam_sidecar.py` line 31、49、59–60、62–75：升級 columns/key v2、global digest/count、subject hashes、atomic write。
- `validate_streamed_longphase_sidecar.py` line 47–81：duplicates任何一筆都 fail；重算 v2 multiset/stream digest。
- 同檔 line 102–133：validation output加入 sidecar/index/input/output/producer/source/env subject hashes；atomic write。

若新 sidecar仍是 schema 1.0、只有 `CHROM START END QNAME FLAG MAPQ CIGAR_B2 HP PS`，downstream只能標為 `coordinate_join_v1`，不得升級宣稱 `exact_alignment_payload_join`。

---

## 6. `run_layered_7samples_newbb.sh` precise patch points

以下以 snapshot SHA `d8a803...` 為基準：

| Lines | Required change | Verify |
|---|---|---|
| 5–20 | 只解析 bootstrap args；禁止在 strict validation前用 jq導出 counts。加入 `--manifest --schema --run-id --run-parent`；validate env numeric/ranges；reject unknown `SM_*` | 壞/空 manifest在建立 published root前 exit 2 |
| 22–24 | 取代 shallow jq gate為 validator schema gate；只接受 v3 comprehensive | 2.0/2.1、0/0、unknown field均 fail |
| 25–31 | parent directory `flock`；建立 `.RUN_ID.staging.PID`；manifest snapshot/lock atomic；完成 launch receipt後才 rename成 final root | 同 run ID並發只能一個；preflight fail不留下看似正式 root |
| 33–36 | `run_state.json` atomic state machine；event append加 `flock`或每事件 atomic file；detail sanitize tabs/newlines | 2 parallel workers不交錯/破壞 TSV/JSON |
| 38–42 | `manifest_value` 改讀 frozen lock；更佳為 validator輸出每 sample canonical JSON，worker不動態 jq source manifest | launch後修改 source manifest不影響 worker |
| 44–51 | 不只 `-f`；呼叫 frozen identity recheck（realpath/type/hash/storage token/index） | same path被替換/retarget立即 fail |
| 53–86 | 移除 worker-time首次 preflight；只接受 lock中 `sample.pass=true`且 exact fields；safe sample ID/dir containment；sidecar v2必填 | 不會到 line 79才發現空 sidecar |
| 88–107 | `input_files.tsv/input.sha256`由 frozen lock生成；加入 BAM identity、sidecar本體、integration、producer/env/source digests；禁止重新觀察後不比較 | manifest record與實際 consumed artifact一一對應 |
| 115–122 | 從 extracted source bundle執行；傳 explicit sidecar index/lock digest；env只用 allowlist；consumer啟動即 verify digest | working tree在run中修改不改執行 code |
| 125–130 | tracked PID/process group；wait任何 child non-zero即取消同 sample siblings並標 FAILED | split fail不留下其他 orphan |
| 132–139 | 保留 exposure exact gate，但另要求 identity v2、lock/sidecar/index digest一致、duplicate/malformed/extra=0 | 不能只由 producer boolean過關 |
| 143–175 | 所有 JSON/TSV/GZ先寫 `.partial`，producer成功+schema/hash驗證後 atomic rename；`output.sha256`寫後立即 `sha256sum -c` | SIGKILL/direct-write不會留下正式半檔 |
| 179–189 | 改為 source bundle：包含 runner本身、validator、verifier、全部 imported scripts/schema；保存 bytes、manifest、tar hash、git commit/tree/diff/untracked inventory | run root可在repo改變後重建 source |
| 190–191 | validator必須在 publish/dispatch前產生 frozen lock；runner使用 validator exit taxonomy；failure不得繼續 | worker啟動前所有7 samples已完整 PASS |
| 192–201 | params atomic且綁 lock/env/source digest；強制 `VERIFY_EVERY=1`,`ANALYSIS_TREE_CAP=0`; numeric range；記所有 effective defaults/origins | comprehensive profile不可被 env靜默降級 |
| 203–208 | 以明確 PID pool取代 opaque `xargs -P`，或建立可追蹤 process groups；安裝 trap與 heartbeat；samples只來自 frozen lock | SIGTERM可收齊 children；dataset set exact 7 |
| 210–217 | verifier只讀 frozen lock/receipt；先 post-input recheck；verification hashes readback；state `SUCCEEDED`後 atomic `_SUCCESS` | 沒 `_SUCCESS`就不是成功；empty/extra/swap均 fail |

### 6.1 Source bundle minimum

至少 bundle：runner、validator、verifier、`sm_linkage_genomewide.py`、`sm_multilocus_combinations.py`、tree solver、layered reconstruction、region view、site ledger、所有 imported local modules、JSON Schemas、provenance utility。Manifest逐檔記 path/mode/size/SHA；bundle本體另 SHA。若 Git dirty，要保存 `git diff --binary` 與 untracked source bytes；只記 `HEAD` 不夠。

**執行 authority必須是 bundle extract directory**，不是 bundle建立後仍跑 `$SCD/*.py`。

### 6.2 Environment lock minimum

- Python executable realpath/SHA/version；`pysam`與所有 Python distributions exact version/lock。
- `jq,samtools,bcftools,bgzip,tabix,xargs,/usr/bin/time` realpath/version/hash。
- OS/kernel/glibc、hostname/storage ID、locale、timezone。
- `LC_ALL=C`,`TZ=UTC`,`PYTHONHASHSEED=0`等 effective determinism settings。
- 完整 allowlisted `SM_*` value + origin；任何 unknown `SM_*` reject。
- environment lock SHA寫進 launch receipt與 outputs。

---

## 7. Atomic lifecycle、trap、heartbeat

### 7.1 State machine

```text
CREATING -> PREFLIGHT -> READY -> RUNNING -> VERIFYING -> SUCCEEDED
     \          \          \         \           \
      +----------+----------+---------+-----------> FAILED
                              signal ------------> ABORTED
```

- 只能按合法 edge前進；state JSON含 sequence number、previous digest、timestamp、stage、sample、launcher PID、reason/error code。
- `SUCCEEDED` 只能由 verifier exit 0 + post-input identity pass + output hash readback建立。
- `_SUCCESS` 是最後一個 atomic marker；內容至少含 run receipt/frozen lock/source/env/verification hashes。
- `FAILED/ABORTED` 保留完整 staging/artifacts供稽核，不刪除；必要時 rename到明確 failed namespace。

### 7.2 Trap contract

Runner一開始安裝 `EXIT ERR INT TERM HUP` handlers：

- 記錄 exit code/signal/current stage/current sample。
- 停止 heartbeat；向 tracked process groups送 TERM，grace period後才 KILL；逐一 `wait`。
- atomic寫 terminal state；若尚未 SUCCEEDED，絕不建立 `_SUCCESS`。
- trap本身 idempotent，避免 EXIT再次覆寫原始 failure。
- 不靠 `xargs`猜 grandchildren；每個 child/process group都登記 PID、start time、command digest。

### 7.3 Heartbeat

每 60 秒 atomic更新 `heartbeat.json`：

`seq,wall_time_utc,monotonic_seconds,host,launcher_pid,state,stage,active_samples,child_pid_start_times,last_event_seq,frozen_lock_sha256,source_bundle_sha256`。

監控/驗證規則：RUNNING/VERIFYING若超過 180 秒未更新，標 `E_HEARTBEAT_STALE`；PID存在但 start time不同也視為 PID reuse/stale。Heartbeat只代表 liveness，不代表 progress或success。

### 7.4 Atomic file rules

所有 JSON/TSV/checksum/summary：same-filesystem temp → flush/fsync → validate parse/schema → `os.replace`/`mv` → fsync directory。BGZF/Tabix需先完成 data、再index、再重開fetch smoke、再一次 publish成 artifact set；不可 data正式、index partial。

---

## 8. `verify_layered_v2.py` precise patch points

| Current lines | Required change | Independent check |
|---|---|---|
| 11–16 | 共用 canonical digest/atomic read-write/provenance helper | verifier自己的 source也由bundle hash固定 |
| before 28 | 新增 run receipt/frozen lock/source/env/state驗證 | 任何 digest mismatch在讀scientific outputs前fail |
| 28–40 | sample ID path containment、symlink escape、exact expected filenames；extra/missing sample dirs fail | sample swap/path traversal不能混入 |
| 42–50 | 每個 output先 schema parse、hash readback、embedded sample/lock digest match | 不只因檔案存在就信任 |
| 51–95 | 保留現有 scientific recomputation；補 numeric type/range/NaN/duplicate unit checks | malformed truthy值不能通過 |
| 96–99 | 不只信任 `check_exact_sidecar_join`/site-ledger `pass`；重算 part sums、digest consistency、ledger row/hash/conservation | producer寫假 boolean仍 fail |
| 101–107 | CN semantics從 frozen schema判斷，不用一個 string推 missing | unavailable/neutral不混淆 |
| 108–132 | region/detail keys唯一、sample identity、input/output lock digests；必要時cross-file unit IDs exact set | two JSON來自不同run會fail |
| 137–146 | CLI改成 `--run-root --frozen-lock --launch-receipt --output`；拒讀 mutable source manifest；assert exact 7/6與sample set | empty `samples=[]`、extra dir、duplicate ID fail |
| 147–157 | summary加入 receipt/lock/source/env/output-manifest hashes與provenance verdict | summary可追到exact execution |
| 158–169 | JSON/TSV atomic；TSV寫完後readback count/schema/hash | 中斷不留下正式summary |
| 170–174 | exit 0前驗 post-input identity、state transition與verification hash | verifier pass才可由runner建立 `_SUCCESS` |

Verifier不應重新使用 source manifest作 expected truth；expected dataset/artifacts必須來自 frozen lock。Source manifest可以被保存作 provenance，但 launch後修改它不應改變驗證集合。

---

## 9. Launch gate（任何 worker 前全部為 true）

### 9.1 Schema/scope

- [ ] schema name/version exact v3；沒有 unknown fields。
- [ ] dataset_count=7、samples length=7、unique IDs=7；biological IDs unique=6。
- [ ] sample IDs符合 `[A-Za-z0-9][A-Za-z0-9_.-]*`，無 `/`, `..`, control chars。
- [ ] scope exact chr1–22，沒有 pilot/subset flag。
- [ ] params profile exact：`VERIFY_EVERY=1`,`ANALYSIS_TREE_CAP=0`；所有 integer/range合法。

### 9.2 Inputs/producers

- [ ] 每個 artifact path/type/index/identity exact；沒有 silent symlink retarget。
- [ ] BAM header/ref/index/quickcheck/sort通過；sample binding policy明示。
- [ ] Somatic VCF role/ClairS paired producer/filter/index/record parsing符合contract。
- [ ] CN availability/source/coords/states/overlap/unlisted semantics一致。
- [ ] 7/7 read-tag mode均是 external sidecar，無 mixed mode。
- [ ] clean producer argv/input binding exact；`--truth-vcf/--truth-bed` absent。
- [ ] sidecar/index/validation full hashes互綁；actual Tabix index就是manifest index。
- [ ] global alignment identity v2 count/digest exact；missing/extra/duplicate/conflict/malformed皆 0。

### 9.3 Provenance/lifecycle

- [ ] source bundle完整、包含runner、hash readback；執行路徑指向bundle。
- [ ] environment lock完整；無 unknown `SM_*`。
- [ ] frozen lock atomic完成並readback；source manifest再次hash一致。
- [ ] run-parent lock取得；run ID未存在；staging與final在同filesystem。
- [ ] trap、state、heartbeat、child registry已初始化。
- [ ] launch receipt綁定 manifest/lock/source/env/params hashes，state從 PREFLIGHT到READY。
- [ ] **以上任一 false：worker count必須仍為 0。**

---

## 10. Adversarial test matrix（設計；本輪未執行）

測試應使用 `tempfile` + tiny synthetic BAM/VCF/BED/BGZF sidecar；禁止以 production 7-sample run作 unit test。每個 case同時 assert：exit code、stable error code、worker marker不存在、frozen lock/_SUCCESS不存在。

### 10.1 Manifest/schema

| ID | 故意破壞 | Expected |
|---|---|---|
| M01 | `dataset_count=0,biological=0,samples=[]` | exit 2 `E_SCHEMA_EMPTY`; no lock/worker |
| M02 | 兩筆同 `sample` | exit 2 `E_SAMPLE_DUPLICATE` |
| M03 | sample=`../escape`或含 newline/tab | exit 2 `E_SAMPLE_ID_UNSAFE` |
| M04 | declared counts與array/unique biological不符 | exit 2 `E_COUNT_MISMATCH` |
| M05 | typo field `read_tag_sidecr` | exit 2 `E_UNKNOWN_FIELD`；不得被additionalProperties吞掉 |
| M06 | 6 sidecar + 1 embedded-BAM mode | exit 2 `E_MIXED_TAG_MODE` |
| M07 | scope少chr22或加chrX | exit 2 `E_SCOPE_NOT_COMPREHENSIVE` |
| M08 | `cn_source=unavailable`但仍有CN path，或measured卻無path | exit 2 `E_CN_CONTRACT` |
| M09 | integration/integer-CN欄位缺失但下游會使用 | exit 2 `E_REQUIRED_ARTIFACT` |
| M10 | schema 2.0/2.1送入v3 runner | exit 2 `E_SCHEMA_VERSION`; 不做best-effort |

### 10.2 BAM / producer / storage identity

| ID | 故意破壞 | Expected |
|---|---|---|
| B01 | `@PG CL`加入 `--truth-vcf` | exit 4 `E_TRUTH_FLAG_PRESENT` |
| B02 | 加 `--truth-bed`或response-file等價參數 | exit 4 `E_TRUTH_FLAG_PRESENT` |
| B03 | LongPhase `@PG`缺失、兩個ambiguous leaf、PP chain斷裂 | exit 4 `E_PRODUCER_AMBIGUOUS` |
| B04 | version/MAPQ/tagSupplementary與manifest effective options不同 | exit 4 `E_PRODUCER_OPTION_MISMATCH` |
| B05 | producer command綁另一份 somatic/germline/normal/tumor/reference | exit 4 `E_INPUT_BINDING_MISMATCH` |
| B06 | BAM `SO`非coordinate、reference dictionary不同、BAI屬於另一BAM | exit 3 `E_BAM_SORT_OR_INDEX` |
| B07 | full-hash mode改一 byte | exit 3 `E_HASH_MISMATCH` |
| B08 | storage mode以同size/mtime替換inode或retarget symlink | exit 3 `E_PATH_RETARGETED` / `E_STORAGE_IDENTITY` |
| B09 | 在 preflight後、dispatch前改BAM | re-observe fail；worker marker仍不存在 |
| B10 | `@RG SM`空但schema要求match；或值屬另一sample | exit 3 `E_SAMPLE_BINDING`；若空白例外須manifest明示 |

### 10.3 VCF/CN

| ID | 故意破壞 | Expected |
|---|---|---|
| V01 | 插入 malformed/non-biallelic/non-SNV record；舊parser會skip | strict分類/守恆；按role fail `E_VCF_RECORD_CONTRACT`，不得沉默 |
| V02 | stale `.tbi/.csi`、manifest index不是actual opened index | exit 3 `E_INDEX_MISMATCH` |
| V03 | PASS backbone含非PASS FILTER或duplicate/unsorted keys | exit 3 `E_VCF_FILTER_OR_ORDER` |
| V04 | ClairS header存在但不是paired producer/normal input不符 | exit 4 `E_CALLER_ROLE_MISMATCH` |
| C01 | CN負座標、start>=end、overlap、unknown state | exit 2/3 `E_CN_PARSE` |
| C02 | measured CN漏coverage semantics，或unlisted默認在檔間不一致 | exit 2 `E_CN_SEMANTICS` |

### 10.4 Sidecar exact join

| ID | 故意破壞 | Expected |
|---|---|---|
| S01 | manifest指定 index A，但旁邊預設 `.tbi` 是 B | explicit index hash/open fail `E_INDEX_MISMATCH` |
| S02 | 把 PASS validation JSON複製到另一sidecar | subject sidecar SHA mismatch `E_SIDECAR_SUBJECT_MISMATCH` |
| S03 | sidecar改一行HP但validation不改 | sidecar SHA/validation subject mismatch |
| S04 | 移除一個BAM alignment key | `E_JOIN_MISSING` |
| S05 | 增加不存在BAM的key | `E_JOIN_EXTRA` |
| S06 | 同key不同HP/PS | `E_JOIN_CONFLICT` |
| S07 | 同key相同HP/PS重複一行 | `E_JOIN_DUPLICATE`；不能只檢conflict |
| S08 | coordinates/CIGAR不變但BAM sequence或quality改變 | identity v2/global digest mismatch `E_JOIN_PAYLOAD_MISMATCH` |
| S09 | MAPQ改變而其他key不變 | identity v2 mismatch |
| S10 | malformed 8/10-column row、unknown HP | `E_SIDECAR_PARSE`; 不得skip |
| S11 | sidecar缺chr22/任一contig，fetch回empty | `E_SIDECAR_SCOPE`; 不得只記 `chr22_fetch_nonempty:false`後仍pass |
| S12 | runtime sidecar path消失或fetch exception | worker startup fail；不得fallback BAM HP/PS |

### 10.5 Frozen provenance / lifecycle / verifier

| ID | 故意破壞 | Expected |
|---|---|---|
| P01 | freeze後修改 source manifest | workers仍讀原lock；source recheck記mismatch但不改dataset set；policy可abort |
| P02 | source bundle改一byte/少一import | launch前 exit 6 `E_SOURCE_BUNDLE_MISMATCH` |
| P03 | 改pysam/Python executable或注入unknown `SM_*` | launch前 exit 6 `E_ENVIRONMENT_MISMATCH/E_UNKNOWN_SM_ENV` |
| P04 | run中修改working-tree source | bundled execution結果不受影響；end bundle readback仍pass |
| L01 | writer在rename前被SIGKILL | 正式artifact不存在；只有`.partial`/FAILED state |
| L02 | SIGTERM dummy worker | state ABORTED、children全waited、無 `_SUCCESS` |
| L03 | 一個split exit non-zero | siblings取消、sample/run FAILED、final verifier不啟動 |
| L04 | heartbeat停止>180s | `E_HEARTBEAT_STALE`，監控/verify拒絕success |
| L05 | 兩個launcher同RUN_ID | 只有一個取得flock，另一個在worker前fail |
| R01 | verifier收到empty frozen sample set | exit 7 `E_DATASET_SET_MISMATCH`，不是 `all([])` |
| R02 | swap兩個sample directories或嵌入sample ID不符 | exit 7 `E_SAMPLE_OUTPUT_BINDING` |
| R03 | 增加第8個sample dir/漏一個dir | exit 7 `E_DATASET_SET_MISMATCH` |
| R04 | 修改output但保留producer `pass:true` | output hash/recomputation fail |
| R05 | site ledger/layered/view來自不同lock/run | embedded lock digest mismatch |
| R06 | verifier summary direct-write中斷 | 正式summary/TSV不存在；無 `_SUCCESS` |

### 10.6 Proposed test files

- `InterSubMod/tests/test_layered_v2_fail_closed.py`：schema、identity、producer、sidecar、verifier tiny fixtures。
- `InterSubMod/tests/test_layered_v2_lifecycle.py`：dummy worker、signal、heartbeat、atomic publish、flock。
- `InterSubMod/tests/fixtures/layered_v3/`：2–3-read tiny coordinate-sorted BAM/BAI、small VCF/index、CN、sidecar v2/validation；以factory生成亦可。
- `InterSubMod/scripts/test_layered_v2_launch_gate.sh`：只呼叫 fixture profile，不接 production paths。

測試本身要用 temporary run parent，結束後依 repo政策保留/由test framework清理其tempdir；不可對 production output執行 destructive cleanup。

---

## 11. Minimum implementation sequence（Step → Verify）

1. **固定 v3 schema 與 error taxonomy**  
   → 驗證：M01–M10 全部在 worker marker前以 expected code失敗；valid 7/6 fixture產生 canonical intent hash。

2. **升級 sidecar identity/producer validation**  
   → 驗證：S01–S12、B01–B05 全部fail-closed；valid fixture BAM與sidecar global count/digest exact equal。

3. **重構 PREFLIGHT 成 atomic frozen-lock producer**  
   → 驗證：全 artifact expected=observed；failure不產生 lock；valid lock readback SHA一致。

4. **建立 source bundle + environment lock**  
   → 驗證：P02/P03失敗；修改working tree不改bundle execution；bundle包含runner本身與所有imports。

5. **重接 RUNNER lifecycle**  
   → 驗證：L01–L05；state edges合法、trap收children、heartbeat可監測、`_SUCCESS`只在最後出現。

6. **重構 VERIFY 只吃 frozen artifacts**  
   → 驗證：R01–R06；空集合/extra/swap/false boolean/output mutation全fail。

7. **小型 end-to-end fixture smoke**  
   → 驗證：exit 0、exact expected outputs、所有hash readback、terminal SUCCEEDED、`_SUCCESS`內容可重算。

8. **先 review audit evidence，再決定 production launch**  
   → 驗證：7/6 full scope launch checklist逐項有machine-readable PASS；若使用storage identity，報告明確標低於full hash的assurance。

本文件不執行第 1–8 步的 source change或run；它們是下一輪 implementation/verification acceptance criteria。

---

## 12. 本輪實際執行紀錄

### 12.1 Inputs

- `InterSubMod/docs/methodology/_assets/20260627_subclone_4axis_teaching/scripts/run_layered_7samples_newbb.sh`
- `InterSubMod/docs/methodology/_assets/20260627_subclone_4axis_teaching/scripts/validate_layered_v2_inputs.py`
- `InterSubMod/docs/methodology/_assets/20260627_subclone_4axis_teaching/scripts/verify_layered_v2.py`
- `InterSubMod/docs/methodology/_assets/20260627_subclone_4axis_teaching/scripts/sm_linkage_genomewide.py`
- `InterSubMod/docs/methodology/_assets/20260627_subclone_4axis_teaching/scripts/sm_multilocus_combinations.py`
- `InterSubMod/docs/methodology/_assets/20260627_subclone_4axis_teaching/scripts/capture_longphase_tagged_bam_sidecar.py`
- `InterSubMod/docs/methodology/_assets/20260627_subclone_4axis_teaching/scripts/validate_streamed_longphase_sidecar.py`
- `InterSubMod/docs/methodology/_assets/20260627_subclone_4axis_teaching/data/layered_v2_input_manifest.json`
- `InterSubMod/docs/methodology/_assets/20260627_subclone_4axis_teaching/data/longphase_production_sidecar_manifest.json`
- Exact historical run root：`/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260710_232501_layered_reconstruction_v2/`

### 12.2 Read-only commands與關鍵輸出

1. Source snapshot：`sha256sum`, `stat`, `wc -l`, `git status --short -- <exact files>`  
   → 實際：四個指定檔案 hash/mtime/line count如 §1.2；runner/consumer modified，validator/verifier untracked。

2. Manifest key audit：

   ```bash
   jq '{schema_version,dataset_count,missing_sidecar_fields:[.samples[]|{sample,missing:(["read_tag_sidecar","read_tag_sidecar_index","read_tag_validation","caller_raw_vcf","longphase_recalibrated_all_vcf","somatic_vcf_role"]-keys)}]}' \
     InterSubMod/docs/methodology/_assets/20260627_subclone_4axis_teaching/data/layered_v2_input_manifest.json
   ```

   → 實際：7/7 samples皆缺上述六欄。

3. Empty-manifest expression audit：

   ```bash
   jq -e '(.schema_version == "2.0" or .schema_version == "2.1") and .dataset_count == (.samples | length)' \
     <<< '{"schema_version":"2.0","dataset_count":0,"biological_sample_count":0,"samples":[]}'
   python3 -c 'print(all([]))'
   ```

   → 實際：`runner_empty_manifest_gate_exit=0`；`python_all_empty=True`。

4. Existing run source readback：`(cd <exact run root> && sha256sum -c code.sha256)`  
   → 實際：3 OK / 5 FAILED；exit 1。

5. BAM header audit：由 manifest exact 7 BAM paths逐一 `samtools view -H`，只解析 `@HD/@RG/@PG`  
   → 實際：7/7 `--truth-vcf`、6/7 `--truth-bed`、7/7 MAPQ20 + tagSupplementary、7/7 RG SM empty。

6. Syntax-only：`bash -n RUNNER`；Python `ast.parse`三檔  
   → 實際：`bash_n_exit=0`、`ast_parse=3/3`。

### 12.3 Output

本輪唯一新增檔案：

`InterSubMod/research/20260710_layered_reconstruction_v2/audit_notes/05_fail_closed_wiring_audit.md`

未修改 source/manifest，未啟動新 run，未寫入 production output。

---

## 13. 最終 acceptance definition

只有同時滿足以下條件，才可以把新 7-dataset result稱為「可重現、fail-closed、clean sidecar」：

1. strict v3 schema exact 7/6且全部artifact/role完整；0/0、mixed mode、unknown field必定fail。
2. 7/7 producer command與inputs可追溯，truth flags token-level absent；舊 truth-restricted BAM tags明確 ignored。
3. sidecar本體/index/validation/producer/BAM互相以full hashes與global alignment identity v2綁定；missing/extra/duplicate/conflict/malformed皆0。
4. source bytes與environment皆保存在run root的content-addressed bundle/lock；worker不讀mutable repo/manifest。
5. lifecycle具atomic writes、flock、tracked children、trap、heartbeat、合法state transitions；失敗/中斷無 `_SUCCESS`。
6. verifier拒絕空集合，獨立重算provenance/output invariants，exact sample set 7/6，post-input identity不變。
7. §10 adversarial matrix全通過，且 tiny end-to-end fixture有可重算的 frozen lock/source/env/output hashes。
8. 正式 run 完成後，run root單獨即可回答：用了哪些 bytes、哪個command/env、每個sample哪個input、sidecar如何對上BAM、何時/為何成功。

在這些條件以前，現有的 `7/7 PASS` 應表述為「historical internal invariant pass」，不能升級為「frozen provenance pass」或「clean sidecar exact-join pass」。
