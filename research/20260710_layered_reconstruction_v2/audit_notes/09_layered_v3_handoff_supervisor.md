<!--
建立時間: 2026-07-11
目標: 記錄 LongPhase production 7/7 完成後，fail-closed closeout → 7 receipts → v3 manifest → full runner 的 one-shot handoff supervisor。
處理範圍: Task Type B；7 datasets / 6 biological samples / chr1-22；本文件只記錄 synthetic/mock 驗證，未執行正式 handoff 或 full run。
關聯檔案: InterSubMod/scripts/continue_layered_v3_after_producer.py; InterSubMod/scripts/test_continue_layered_v3_after_producer.py; InterSubMod/research/20260710_layered_reconstruction_v2/launch_plans/20260711_layered_v3_handoff_launch_plan_unreviewed_v1.json
-->

# Layered-v3 production handoff supervisor audit

> **TL;DR — one-shot handoff supervisor 與 execute=false launch-plan template 已完成，28/28 synthetic tests PASS；正式 supervisor 未執行、active producer artifacts 未被修改、full run 未啟動（影響：高，程式信心：高，真資料執行信心：待 7/7 gate）。**

## 1. 任務分類與服務目標

- Task Type：**B — Comprehensive validation**。
- Scope：固定 7 datasets、6 biological samples、chr1–22；不可 subset。
- 服務目標：G4 reproducibility、G5 externally verifiable provenance。
- 本次授權邊界：只實作、mock 測試與建立 `execute=false` template；**禁止**執行正式 supervisor、禁止啟動 full run。

## 2. Authority 與固定路徑

| 角色 | Authority |
|---|---|
| Production root | `/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260711_longphase_s_production_sidecars_v1` |
| Reviewed producer manifest | `InterSubMod/docs/methodology/_assets/20260627_subclone_4axis_teaching/data/longphase_production_sidecar_manifest.json` |
| Frozen producer snapshot | production root 的 `input_manifest.json`；必須與 reviewed manifest byte-identical |
| Base manifest | `InterSubMod/docs/methodology/_assets/20260627_subclone_4axis_teaching/data/layered_v2_input_manifest.json` |
| Real method receipt | `InterSubMod/research/20260710_layered_reconstruction_v2/probes/20260711_COLO829_chr1_4386684_4388348_coordinate_join_v1/equivalence_probe.json` |
| Synthetic method receipt | 同目錄 `synthetic_contract_receipt.json` |
| Candidate run ID | `20260711_layered_reconstruction_v3_full_no_truth_v1` |
| Candidate run parent | `/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds` |
| Handoff lock | production root 的 sibling `.20260711_longphase_s_production_sidecars_v1.handoff.lock` |

重要判斷：normalizer 的 `--production-manifest` 使用 **run-root snapshot**，因 receipt 要記錄真正 producer runtime authority；finalizer 與 v3 preparer 的 source-set authority 使用 **reviewed method manifest**，避免拿 run-root copy 自我證明。

## 3. Step → Verify 流程

1. 讀 launch plan 並取得 handoff flock
   → 驗證：plan schema `1.0.0`、固定 paths、5 inputs + 25 executable/source/schema/tool pins 全部 hash/readlink 相符；第二個 supervisor 得 `E_HANDOFF_LOCKED`。
2. 等待 producer（僅明示 `--wait-for-producer`）
   → 驗證：每 30 秒只重讀 exact `run_status.tsv`、`verification_summary.json`、`_SUCCESS` authority 與 `ps` argv；不做磁碟遞迴掃描。
3. 判斷 7/7 terminal state
   → 驗證：每 sample 恰一個 START 後恰一個 PASS；`ALL verify PASS` 唯一且最後；verification 恰 7/7；wrapper、LongPhase、capture 全退出。
4. 若尚未 closeout，前景呼叫 frozen finalizer
   → 驗證：finalizer exit 0；`closeout/production_closeout.json`、`artifacts.final.sha256`、`_SUCCESS` 三者互相 hash-bind；重新計算所有 sha256 manifest entries。
5. 依固定 sample order 執行 7 次 v3 normalizer
   → 驗證：每次執行前仍為 production PASS、closeout 未 drift、所有 pins 未 drift；每 sample 只能產生一個 valid receipt，failure JSON 不得被當 receipt。
6. 建立 v3 source manifest
   → 驗證：preparer output/failure 都是新路徑；output 為 schema `3.0.0`、固定 manifest ID、7/6、固定 sample order。
7. 啟動前重跑所有 gates
   → 驗證：producer/closeout/receipts/pins/source manifest/run-id/process overlap 全部 readback；任何差異立即非零退出。
8. 前景 `Popen` v3 runner並持續持有 handoff flock
   → 驗證：記錄 exact argv、controlled env/PATH、tool resolution、PID、log、exit code；runner 非零不 retry，亦不建立 handoff success。
9. runner exit 0 後讀回 layered `_SUCCESS`
   → 驗證：run ID、mode=`full`、scope=`chr1-22`、7 datasets、verification path/hash 皆一致，才建立 `handoff.success.json`。

## 4. WAIT 與 fatal 分界

只有下列狀況可視為 `E_PRODUCER_PENDING` 並繼續等待：

- status 是合法 append-only progress（未開始、START、START→PASS）；
- 對應 producer wrapper / LongPhase / capture process 仍存在；或
- 已 7/7 verified，但上述 process 尚未完全退出。

下列一律 fatal，不等待、不重試：

- 任一 `FAIL`、未知 sample/stage/status、PASS-before-START、duplicate START/PASS、ALL PASS 非唯一/非最後；
- wait 過程 `run_status.tsv` 被 truncate/rewrite（非 append-only）；
- producer status 未完成但已無 producer process；
- reviewed/snapshot/source/schema/binary/tool/code hash drift；
- v2 或另一個 v3 launcher overlap；
- partial closeout、partial receipts、run ID/workspace 已存在；
- normalizer、preparer、runner 任一非零；
- post-closeout artifact drift 或 `_SUCCESS` binding 不成立。

## 5. Execution gate

預設沒有 mode flag時仍是 dry-run；正式執行必須同時滿足三層：

1. CLI `--execute-reviewed-plan`；
2. JSON `execution_authorization.execute=true` 且 `reviewed_by/reviewed_at_utc` 非空；
3. CLI `--expected-launch-plan-sha256 <out-of-band SHA>` 與 plan bytes 完全一致。

目前 template 明確為：

```json
{"execute": false, "reviewed_by": null, "reviewed_at_utc": null}
```

因此 template **不可執行正式 handoff**。reviewer 必須逐段審核後，另以 `apply_patch` 建立 versioned `execute=true` plan，再從該新檔計算 out-of-band SHA；不可原地偷改 template 後宣稱已 review。

## 6. Truth flag 與命令判斷

命令 policy 是 token-level：只拒絕 `--truth-vcf`、`--truth-bed` 及其 `--option=value` 形式。不能用 substring `truth`，否則合法 run ID `...no_truth_v1` 會被誤判。

相同 policy 由 `build_commands()` 統一套用於 dry-run 與 execute，避免兩條路徑判斷分叉。producer closeout、receipt normalizer、manifest preparer、runner 的 argv 都不含 truth flag。

## 7. Frozen provenance 與 PATH

Launch plan 固定 5 個 inputs 與 25 個會執行/載入的 artifacts：

- supervisor、Python、LongPhase binary；
- finalizer、normalizer、preparer、v3 runner/lifecycle/validator/verifier；
- 3 schemas、6 scientific modules；
- `ps`、`stat`、`samtools`、`bcftools`、`bgzip`、`tabix`。

實際 child env 移除所有 inherited `SM_*`，固定：

```text
LC_ALL=C
TZ=UTC
PYTHONHASHSEED=0
PYTHONDONTWRITEBYTECODE=1
PATH=/usr/local/bin:/usr/bin
```

啟動前用 `shutil.which` 證明 PATH 精確解析至 pinned tools：samtools/bgzip/tabix → `/usr/local/bin`；stat/bcftools → `/usr/bin`。`ps` 直接用 pinned absolute path。runner receipt 保存 PATH 與解析表。

完整 30 pins 與 full SHA-256 在：

`InterSubMod/research/20260710_layered_reconstruction_v2/launch_plans/20260711_layered_v3_handoff_launch_plan_unreviewed_v1.json`

Template readback：5 inputs + 25 executables = 30 pins；plan SHA-256：

`9fa47571f0a116a6fd0994df360f4e1b5d1592b0b8a66cde20f98527c8aecb95`

## 8. No-overwrite / no-delete / no-retry

- Supervisor JSON 使用 final path `O_CREAT|O_EXCL`，寫入後 `fsync(file)` + `fsync(parent)`；已存在直接 `E_OUTPUT_EXISTS`，沒有 rename-overwrite race。
- 不含 `unlink`、`remove`、`rmtree`，不刪除任何檔案。
- command stage 只執行一次；非零立即停止，不 retry。
- execute 前 workspace、receipt/failure、manifest/failure、run root 都必須不存在；任何 partial state 留作稽核，不自動清理。
- dry-run 不建 handoff workspace或 scientific output，但為了全程 flock，會建立/開啟 sibling lock inode。

## 9. Synthetic test evidence

### 輸入

- `InterSubMod/scripts/continue_layered_v3_after_producer.py`
- `InterSubMod/scripts/test_continue_layered_v3_after_producer.py`
- 全部 fixtures 在 `tempfile`；所有 subprocess/foreground runner 使用 mock。

### 命令

```bash
cd /big7_disk/liaoyoyo2001/InterSubMod
python3 -m py_compile scripts/continue_layered_v3_after_producer.py scripts/test_continue_layered_v3_after_producer.py
python3 scripts/test_continue_layered_v3_after_producer.py
```

### 實際輸出

```text
Ran 28 tests in 10.192s
OK
```

28 tests 覆蓋：happy dry-run、authority/command exactness、`no_truth` regression、happy mocked execute、wait START→7/7→quiescent、FAIL/hash drift/status regression、active producer、v2 overlap、第二 supervisor flock、START-only abandoned、PASS-before-START、非 7/7、code/source/post-closeout drift、partial closeout/receipts、existing run ID、三層 execution gate、controlled PATH、no-replace、normalizer/preparer/runner failures與不偽造 success。

### 輸出與 hashes

| Artifact | SHA-256 |
|---|---|
| `InterSubMod/scripts/continue_layered_v3_after_producer.py` | `cc91239368bb79732956f2e08fb608de9dea6dcb70d4d50bf724a1e9cf95e49f` |
| `InterSubMod/scripts/test_continue_layered_v3_after_producer.py` | `9ecd79bc2e126f77964aa9ed8d50052c23647efeb24f24a53d74e9b58bd6da5a` |
| unreviewed launch plan template | `9fa47571f0a116a6fd0994df360f4e1b5d1592b0b8a66cde20f98527c8aecb95` |

## 10. Reviewed plan 與 active wait supervisor（2026-07-11 06:16 +08:00）

另一個同父 Codex reviewer只把unreviewed template的authorization三欄改成`execute=true`、reviewer與UTC時間；移除authorization後兩份JSON的sorted diff為空。Reviewed plan：

`InterSubMod/research/20260710_layered_reconstruction_v2/launch_plans/20260711_layered_v3_handoff_launch_plan_reviewed_v1.json`

Byte SHA-256=`7e0ab871ee2fa15e772d17da75ebe0e836cc7dcc6c2726e7b136a3b305e3da6a`；主代理已將reviewed/template permissions改為`0444`（bytes/hash不變），並獨立重算5 inputs、25 executables與5個PATH tool resolutions全部PASS。

實際等待命令（PID `929656`；start `2026-07-11 06:16:43 +0800`）：

```bash
cd /big7_disk/liaoyoyo2001/InterSubMod
python3 \
  scripts/continue_layered_v3_after_producer.py \
  --launch-plan research/20260710_layered_reconstruction_v2/launch_plans/20260711_layered_v3_handoff_launch_plan_reviewed_v1.json \
  --wait-for-producer \
  --execute-reviewed-plan \
  --expected-launch-plan-sha256 7e0ab871ee2fa15e772d17da75ebe0e836cc7dcc6c2726e7b136a3b305e3da6a
```

`/proc/929656/exe`=`/bip7_disk/liaoyoyo2001/miniconda3/bin/python3.9`；process持有唯一handoff flock。啟動後主代理看到handoff workspace與full-run root皆不存在，表示目前只在exact 30-second polling，不是layered worker已啟動。任何第二個supervisor的dry-run皆以`E_HANDOFF_LOCKED`拒絕。

預期輸出路徑：

- Handoff workspace：`/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/.layered_v3_handoff/20260711_layered_v3_handoff_no_truth_v1/`
- v3 manifest：上述 workspace 的 `layered_input_manifest_v3.json`
- Full run：`/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260711_layered_reconstruction_v3_full_no_truth_v1/`

## 11. 未解風險與 claim limit

1. 本次只有 synthetic/mock execution；真實 7 份大型 sidecar receipt normalization、finalizer I/O、resource gate 尚未實跑。
2. Plan pin 直接執行的 source/schema/binaries/tools，但未把 libc、Python stdlib、pysam shared objects等完整 transitive dependency逐 byte pin；v3 runner仍會另建 environment/source bundle lock。
3. `ps` 判斷依 exact argv basename/root；設計偏 fail-closed，可能把同名檢查程序視為 overlap而暫停。
4. `--wait-for-producer` 無 timeout，符合「等到 terminal gate」語意；操作端仍需監看 producer 是否出現 fatal/abandoned。
5. 每步 hash readback 與實際 `exec` 間仍有極短 TOCTOU；supervisor重查、controlled PATH與runner source bundle降低風險，但不是 kernel-level fd-exec attestation。
6. one-shot 中途失敗會留下 partial receipts/workspace，依政策不 retry、不刪除；後續需人工稽核並建立新 versioned remediation plan。

**結論：程式與 reviewed plan 已通過人工逐段readback；supervisor目前只等待producer。Producer 7/7與runner 300秒resource gate未通過前，layered scientific workers仍維持 NO-LAUNCH。**
