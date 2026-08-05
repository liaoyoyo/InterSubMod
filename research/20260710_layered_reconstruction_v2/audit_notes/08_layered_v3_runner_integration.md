<!--
建立時間: 2026-07-11
目標: 紀錄 layered-v3 runner、fixed full-scope lock、frozen-source worker adapter 與 tiny integration tests。
處理範圍: InterSubMod/scripts/run_layered_v3.py、InterSubMod/scripts/test_run_layered_v3.py；不修改 v2/live producer/consumer。
關聯檔案: InterSubMod/scripts/layered_v3_lifecycle.py；InterSubMod/scripts/verify_layered_v3.py；InterSubMod/research/20260710_layered_reconstruction_v2/audit_notes/05_fail_closed_wiring_audit.md。
任務類型: (B) Comprehensive validation launch-hardening；tiny fixtures only，未啟動 production/full computation。
研究目標: G4 reproducibility；G5 外部可稽核 provenance。
-->

# Layered-v3 runner integration 與 full-scope exclusivity

> **用 SCQA：v3 runner integration 完成 — bundled validator/worker/verifier、7/6 chr1–22 exact roles、fixed global flock、process observation、CN frozen semantics 與最後 `_SUCCESS` 已接線並通過 tiny/cross suites；production launch 仍必須等待有效 7/7 frozen manifest，本文沒有啟動全量計算（影響：高；信心：高）。**

## 1. Outcome

新增：

- `InterSubMod/scripts/run_layered_v3.py`
- `InterSubMod/scripts/test_run_layered_v3.py`

未修改：

- `InterSubMod/docs/methodology/_assets/20260627_subclone_4axis_teaching/scripts/run_layered_7samples_newbb.sh`
- live producer / sidecar consumer / v2 runner

Runner 的 production authority chain：

```text
fixed full-scope flock
  -> RunLifecycle parent run-id flock
  -> PREFLIGHT + exact /proc conflict observation
  -> source bundle (runner/validator/verifier/all local imports/all schemas)
  -> environment lock
  -> bundled validator -> strict candidate lock
  -> relocation-only validator/schema path normalization + lock digest recompute
  -> lifecycle.write_frozen_lock
  -> publish READY
  -> bundled worker adapter (sample parallel <=2; part parallel=1)
  -> output_manifests.json exact 7 path/hash entries
  -> VERIFYING + bundled verifier
  -> lifecycle.succeed -> SUCCEEDED -> `_SUCCESS`
```

## 2. Exact scientific roles

Worker only resolves sample inputs from `frozen_input_lock.json`：

| Operational role | Frozen lock role | Rule |
|---|---|---|
| LongPhase input | `somatic.longphase_input_vcf` | ClairS FILTER=PASS |
| Tree backbone | `somatic.tree_vcf` | distinct LongPhase recalibrated PASS artifact；runner streams every record and rejects non-PASS/empty |
| Site ledger universe | `somatic.caller_raw_vcf` | raw ClairS |
| Recalibration ledger | `somatic.longphase_recalibrated_all_vcf` | LongPhase all input keys / PASS+excluded |
| Read tags | `read_tags.sidecar` | coordinate_join_v1；BAM embedded HP/PS ignored |

Tree VCF is never generated dynamically from `recalibrated_all` and never replaced by ClairS PASS.

## 3. Frozen provenance details

### 3.1 Reviewed contract hashes pinned by the runner

| File | SHA-256 |
|---|---|
| `layered_input_manifest_v3.schema.json` | `2d64b544a3caed8ab01761061786a226f1e20a7800680b3a648f4715db2e81ba` |
| `layered_input_lock_v1.schema.json` | `1a40107f696e19e375836247587c728fec781a1f6c69c34ad60cbec8d72a4fdb` |
| `longphase_production_capture_receipt_v1.schema.json` | `36f9e44286b946d7abd6bff6d963f121203492d4159bbeca01d5f65f6aac55b5` |
| `validate_layered_v3_inputs.py` | `14b3bf26a1a28c31a6fcf546f0872ba11f1491666088a1263a3cff3d73f39618` |
| `verify_layered_v3.py` | `a5f384c1f65030b5c72136d5f5e90acaa29b7f7b5c4b222be5cf35cd6c7076d7` |

任一 reviewed contract drift 會在 bundle/PREFLIGHT 前以 `E_FROZEN_CONTRACT_DRIFT` 拒絕。

### 3.2 Source execution

- Validator 由 bundle core bytes 載入；schemas 從 bundle bytes 複製到 validator runtime resource directory。
- Worker 子程序執行 bundled runner；local imports 由 bundle manifest exact path finder 載入。
- Working-tree source 在 bundle 後改變，不影響 worker bytes；tiny mutation test 將原始 MLHP source 改成 exit 99，bundled worker 仍產出 frozen-v1 結果。
- Environment lock 強制 `LC_ALL=C`、`TZ=UTC`、`PYTHONHASHSEED=0`，unknown `SM_*` 仍由 lifecycle fail-closed。

### 3.3 Output binding

每個 sample 的 `output_manifest.json` 精確記錄：

- frozen lock / launch receipt / environment lock / source bundle manifest/content / input-set SHA；
- somatic four-role hashes；
- canonical frozen `copy_number_contract`；
- 五個 MLHP parts、layered、region、site ledger、site ledger summary 的 path/size/SHA。

Adapter 在 MLHP output 增加 sample/part/chromosome、somatic roles、provenance 與 frozen sidecar identities；若 consumer 已輸出衝突 diagnostic，adapter 拒絕覆寫。

## 4. Full-scope exclusivity

Per-run lock 之外，runner 從 preflight 至 success/failure 全程持有：

`<run-parent>/.layered_chr1_22_7dataset_full.lock`

所以不同 run IDs 也不能重複啟動同一個 7-dataset chr1–22 full run。PREFLIGHT 另建立 `process_observation.json`：

- observer PID + `/proc/<pid>/stat` start time；
- fixed lock path / held / run ID；
- bounded exact process basename policy（含 layered workers 與 `longphase-s`、production sidecar runner/capture/validator）；
- matching PID/start-time/full argv/argv SHA；
- 300 秒 baseline 開始與結束兩份 snapshot timestamp/SHA；兩時點 union 後判斷 conflict，關閉 start-only TOCTOU；
- conflict count 與 PASS。

非空 conflict 以 `E_CONFLICTING_FULL_RUN` 拒絕；receipt 綁 process observation path/hash/zero count；verifier 在 lock 尚持有時重讀 owner。

同一份 hash-bound observation 也包含 resource gate：

| Resource | Production threshold |
|---|---:|
| logical CPUs | >= 8 |
| `/proc/meminfo` MemAvailable | >= 128 GiB |
| run-parent free disk | >= 500 GiB |
| 1-minute load / logical CPUs | <= 1.25 |
| aggregate `/proc/stat` iowait（300 s sample） | <= 20% |
| `/proc/self/mountstats` `/big8_disk` NFS READ bytes_recv delta | < 80 decimal MB/s over same 300 s |

NFS evidence 保存 exact mount path/device/fstype、counter source、start/end counter、delta、UTC 起訖、實際秒數、decimal MB/s 與 threshold。mount stanza/READ row 缺失、counter reset或 mount identity drift 均以 `E_NFS_BASELINE_UNAVAILABLE` fail-closed，並在 staging 保存含 error 與 SHA 的 failure evidence。任一 resource check false 則以 `E_RESOURCE_GATE` 在 bundle/worker 前拒絕。tiny tests顯式使用 0.01 秒 baseline，另以 impossible `max_nfs_read_mbps=0` 證明 NFS rate 在 bundle 前拒絕；production inventory 會拒絕任何低於上述正式門檻的 CLI/env override。

## 5. CN semantics

Runner 不推測 CN：

- measured 必須保存 reviewed `source`、`semantics`、0-based half-open、unlisted=neutral、allowed states、overlap=forbid、CN BED identity；
- unavailable 必須有固定 reason，四個 measured artifacts 全為 null，且不可被解讀成 neutral；
- frozen branch 原樣複製到 output manifest 與 region view。

Tiny full-run fixture 同時覆蓋一個 measured sample 與六個 unavailable samples。

## 6. Step -> Verify 與執行紀錄

1. 建立 v3 runner + bundled adapter
   -> 驗證：AST parse、CLI help、source inventory/pinned hash gate PASS。
2. Invalid/empty/extra/resource/NFS-evidence fail before publish
   -> 驗證：正式 run root不存在，scientific sample directory不存在；stable error code。
3. Source mutation isolation
   -> 驗證：原 worktree stub改為 exit 99 後，7/7 bundled outputs仍含 `source_token=frozen-v1`。
4. Fixed global lock
   -> 驗證：第一把 lock 使用 `first-full-run`，第二個不同 run ID 得 `E_GLOBAL_RUN_LOCKED`；無 staging/root。
5. Producer/layered TOCTOU exclusivity
   -> 驗證：baseline start 空集合、end 才出現 production-sidecar runner，仍以 `E_CONFLICTING_FULL_RUN` 在 bundle 前拒絕；目前 read-only `/proc` 實測亦正確辨識 active runner/capture/LongPhase PIDs。
6. Frozen CN output binding
   -> 驗證：output manifest與region copy-number contract byte-equivalent；measured source及 unavailable reason均保留。
7. Independent verifier
   -> 驗證：final verifier cross suite 16/16 PASS；contract suite 17/17 PASS；lifecycle suite 13/13 PASS。

### Commands

```bash
cd /big7_disk/liaoyoyo2001/InterSubMod
env LC_ALL=C TZ=UTC PYTHONHASHSEED=0 PYTHONDONTWRITEBYTECODE=1 \
  python3 scripts/test_run_layered_v3.py

env LC_ALL=C TZ=UTC PYTHONHASHSEED=0 PYTHONDONTWRITEBYTECODE=1 \
  python3 scripts/test_layered_v3_lifecycle.py

env LC_ALL=C TZ=UTC PYTHONHASHSEED=0 PYTHONDONTWRITEBYTECODE=1 \
  python3 scripts/test_verify_layered_v3.py

env LC_ALL=C TZ=UTC PYTHONHASHSEED=0 PYTHONDONTWRITEBYTECODE=1 \
  python3 docs/methodology/_assets/20260627_subclone_4axis_teaching/scripts/test_layered_v3_contract.py
```

### Actual output excerpts

```text
Ran 12 tests in 9.708s
OK

Ran 13 tests in 1.152s
OK

Ran 16 tests in 57.283s
OK

Ran 17 tests in 32.703s
OK
```

## 7. Remaining launch gate

Component verdict：**runner integration PASS**。

Production verdict：**仍不可由本文單獨 launch**。正常 full invocation 必須先取得 upstream 已完成的 7/7 production sidecar + recalibrated all/PASS VCF、strict v3 manifest，且 production preflight 實際跑完 validator、PASS-only tree scan、global process observation。本文只使用 temporary tiny fixtures，沒有碰 production BAM/VCF、沒有啟動 full run。
