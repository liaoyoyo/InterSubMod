<!--
建立時間: 2026-07-19T12:47:00+08:00
目標: 記錄第一次 v23 formal canonical pytest 因 bound-FD testcase inventory digest 尚未鎖定而 fail-closed
處理範圍: 738 個 topic tests；未發布 canonical JUnit、test-run manifest 或 detached signature
關聯檔案: InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/audit_notes/attested_release_evidence_v2.py；InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/audit_notes/run_attested_canonical_pytest_v2.py
-->

# v23 bound-FD inventory digest 未鎖定失敗紀錄

## 判定

- 狀態：`FAILED_BEFORE_FORMAL_PUBLICATION / EXCLUDED`
- pytest return code：`0`
- JUnit：`738 tests / 0 failures / 0 errors / 0 skipped`
- 失敗原因：正式 runner 透過 `/proc/self/fd/*` 執行每個 test source，JUnit 的 `classname` 為空；先前 direct-path probe 的 `classname` 含 repo path，因此 testcase inventory digest 不同。
- direct-path probe inventory SHA-256：`d9995207f1859c1a45315f49a3317b42d7258f7807540c6e742e5aec3beac2e0`
- bound-FD formal inventory SHA-256：`5206e2ed1a097cce773563cc3e23595631887d91ca2919e4df838a0c7747e5d6`
- 結論：兩者測試名稱與數量相同，但執行入口語意不同；正式 authority 必須鎖定 bound-FD digest。

## 執行

輸入：

- `InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/tests/*.py`
- `InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/audit_notes/run_attested_canonical_pytest_v2.py`

命令：

```bash
exec 3<InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/audit_notes/run_attested_canonical_pytest_v2.py
/usr/bin/env -i PATH=/usr/bin:/bin HOME=/tmp LANG=C.UTF-8 LC_ALL=C.UTF-8 \
  PYTHONHASHSEED=0 PYTHONNOUSERSITE=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 BLIS_NUM_THREADS=1 \
  PYTEST_DISABLE_PLUGIN_AUTOLOAD=1 \
  /bip7_disk/liaoyoyo2001/miniconda3/envs/cnvtools/bin/python -I /proc/self/fd/3
```

實際錯誤片段：

```text
TestRunError: Canonical pytest failed: returncode=0,
counts={'skipped': 0, 'failures': 0, 'errors': 0, 'tests': 738}
```

## 封存輸出

- raw JUnit：`pytest_full_pre_authority_v23_command_parity_raw_evidence/pytest.xml`
  - size=`78,935`
  - SHA-256=`d2d5009d70b46c2e902d18c55a55902fbc14dd075f20075b6d88b0dfa407c520`
- captured stdout：`pytest_full_pre_authority_v23_command_parity_attested_canonical.stdout.txt`
  - size=`1,974`
  - SHA-256=`117af248e783d6d97502b8b771d83825875e65e092d40e9ce7bae86b304c5866`
- captured stderr：空檔，SHA-256=`e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855`
- bound test FD directory：`.pytest_full_pre_authority_v23_command_parity_bound_test_fds/`

所有上述資料已搬入本目錄；原 formal output slots 已清空。此次沒有建立：

- `pytest_full_pre_authority_v23_command_parity_attested_canonical.xml`
- `pytest_full_pre_authority_v23_command_parity_test_run_manifest.v1.json`
- detached Ed25519 signature

因此不得將本次 raw evidence 引用為正式 source-authority 測試證據。

## 修正與再驗證

1. 將 v2 evidence validator 的 pinned testcase inventory 改為 bound-FD digest
   `5206e2ed...e5d6`。
2. 重算 v2 evidence validator SHA-256，更新 v6 release consumer 的 support allowlist。
3. 重新執行同一 v23 no-clobber runner；只有 canonical JUnit、manifest 與 detached signature
   全部成功時才可進外部 Claude A/B review。
