<!--
建立時間: 2026-07-11
目標: 紀錄 layered v3 atomic lifecycle / frozen provenance 元件的可重用契約與 tiny test 證據。
處理範圍: InterSubMod/scripts/layered_v3_lifecycle.py、InterSubMod/scripts/test_layered_v3_lifecycle.py。
關聯檔案: InterSubMod/research/20260710_layered_reconstruction_v2/audit_notes/05_fail_closed_wiring_audit.md。
任務類型: (B) Comprehensive validation 的 launch-hardening 子元件；本文件不是全量資料驗證證據。
研究目標: G4 reproducibility；G5 可外部稽核的 production provenance。
-->

# Layered v3 atomic lifecycle / source-bundle contract

> **用 PREP：v3 lifecycle 子元件完成且 13/13 tiny tests PASS（影響：高，信心：高）；但尚未接線至 production runner，因此不能據此啟動全量 run。**

## 1. 結論與落點

- 元件：`InterSubMod/scripts/layered_v3_lifecycle.py`
- 測試：`InterSubMod/scripts/test_layered_v3_lifecycle.py`
- live v2 runner / validator / verifier：**未修改**
- production 或長計算：**未啟動**

此元件實作 `05_fail_closed_wiring_audit.md` §6–7、§10.5 的 lifecycle/provenance 邊界：parent flock、staging publish、合法 state edge、atomic JSON、heartbeat、tracked process group、source bytes bundle、environment allowlist、signal trap，以及最後才建立 `_SUCCESS`。

## 2. 唯一合法整合順序

```text
RunLifecycle(run_parent, run_id)
  -> begin_preflight()
  -> write_frozen_lock(...)
  -> build_source_bundle(runner, validator, verifier, all imported local files)
  -> capture_environment_lock(...)
  -> publish_ready()             # 此時正式 run root 才第一次出現
  -> begin_running()
  -> tracked children + heartbeat
  -> begin_verifying()
  -> succeed(verification.json)  # SUCCEEDED state 後，_SUCCESS 最後原子建立
```

任何 preflight exception、digest mismatch、unknown `SM_*`、source bundle 少 import、非法 state edge，均不得 publish 正式 run root。Preflight staging 保留供稽核，不刪除。

## 3. Fail-closed invariants

| 邊界 | 強制判斷 | 失敗結果 |
|---|---|---|
| 同 run ID 並發 | run parent 內 per-run `flock(LOCK_NB)` | 第二 launcher `E_RUN_LOCKED`，不建 staging |
| Publish | frozen lock + environment lock + source manifest/content 4 digests 全存在且 readback 相同 | staging state `FAILED`；無正式 root |
| State | 只接受 `CREATING→PREFLIGHT→READY→RUNNING→VERIFYING→SUCCEEDED` 或往 `FAILED/ABORTED` | `E_STATE_TRANSITION` |
| Source authority | bundle 保存 runner / validator / verifier / imported files 的 bytes、mode、size、SHA-256 | 一 byte drift / extra / missing file 均 `E_SOURCE_BUNDLE_MISMATCH` |
| Environment | 只序列化 allowlisted `SM_*`；任何其他 `SM_*` 拒絕 | `E_UNKNOWN_SM_ENV` |
| Child | 每個 child 是獨立 process group，記 PID start time + command digest | 任一 non-zero 取消並 wait siblings |
| Signal | `INT/TERM/HUP` 進同一 idempotent abort path | children 全回收、state `ABORTED`、無 `_SUCCESS` |
| Heartbeat | atomic payload 含 PID identity、state、active samples、children、lock/source digests | >180 秒或 PID start mismatch 為 `E_HEARTBEAT_STALE` |
| Success | verification artifact、receipt、lock/env/source 全部 readback 後才轉 `SUCCEEDED` | `_SUCCESS` 缺席即不算成功 |

## 4. Caller 的不可省略義務

1. `imported_scripts` 必須列出 runner 的**全部** local imports與 JSON schemas；元件只能驗證 caller 提交的集合，不能猜測 dynamic import。
2. `capture_environment_lock()` 必須傳 production 的 Python distributions 與 external tools 清單；空清單只允許 tiny fixture，不能作 production lock。
3. Worker 必須從 `SourceBundleResult.role_paths` 執行 frozen bytes，不可回頭執行 mutable repo path。
4. Frozen lock 必須由 strict v3 input validator 產生；此元件只保護 lifecycle/provenance，不取代 BAM/VCF/sidecar 科學 contract validation。
5. 正式 verifier 必須重算 output binding 與全 7 datasets；`succeed()` 只接受其已寫出的 verification artifact，不替代 verifier。

## 5. Step → Verify 與實際結果

1. 新增 standalone lifecycle component
   → 驗證：兩個目標 scripts 均為新檔；live v2 files 未修改。
2. 使用 tiny TemporaryDirectory fixtures 覆蓋 lifecycle adversarial cases
   → 驗證：parent flock、preflight no-publish、source tamper、unknown env、heartbeat stale、child fail-fast、direct handler 與真實 SIGTERM 均有 test。
3. 重跑完整 unit suite
   → 驗證：exit code 0；`Ran 13 tests in 1.117s`；`OK`。

輸入路徑：

- `InterSubMod/scripts/layered_v3_lifecycle.py`
- `InterSubMod/scripts/test_layered_v3_lifecycle.py`
- runtime-only tiny files：Python `tempfile.TemporaryDirectory`，不使用 production input。

執行命令：

```bash
cd /big7_disk/liaoyoyo2001/InterSubMod
PYTHONDONTWRITEBYTECODE=1 python3 scripts/test_layered_v3_lifecycle.py
```

實際輸出片段：

```text
test_installed_sigterm_trap_aborts_real_subprocess ... ok
test_publish_heartbeat_verify_and_success_marker ... ok
test_source_bundle_copies_bytes_and_detects_tamper ... ok
test_unknown_sm_environment_is_rejected_before_publish ... ok
Ran 13 tests in 1.117s
OK
```

輸出路徑與 SHA-256：

| 路徑 | SHA-256 |
|---|---|
| `InterSubMod/scripts/layered_v3_lifecycle.py` | `90fbd296c7c2c724af542e99b793baf54cbac14e9ec90567b8fb51d5f38aac4f` |
| `InterSubMod/scripts/test_layered_v3_lifecycle.py` | `4a2f5b99d024cf4a4ee986ed6f2b2df5d34bc4fdaa19e90c616ddb0511526b99` |

## 6. Launch gate

**元件 verdict：PASS。Production launch verdict：仍為 NO-LAUNCH。**

必須等 v3 manifest/validator/verifier、sidecar exact binding、完整 source/import inventory 與 production environment/tool allowlist 全部接至本元件，並通過 7-dataset preflight，才能將正式 full run root publish。
