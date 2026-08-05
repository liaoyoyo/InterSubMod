<!--
建立時間: 2026-07-18
目標: 複核 positional singleton 全量稽核程式的來源鎖定、全量重算與發布安全
處理範圍: supplemental auditor、對應 tests、真實 469,849-row 執行
關聯檔案:
  - InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/scripts/audit_positional_singleton_methyl_multigroup.py
  - InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/tests/test_audit_positional_singleton_methyl_multigroup.py
  - /big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/positional_singleton_methyl_multigroup_audit_v1_source_attested/positional_singleton_audit_summary.json
-->

# Positional singleton 稽核程式獨立 agent 複核

## 審查資訊

- Reviewer agent: `019f72b0-3957-7ed3-bae8-0659ce2652d2`
- Auditor SHA-256:
  `5b259c7bbfc9e43e858cb01fae0efce71030075562e3f92bcb25c7f6bb23695f`
- Tests SHA-256:
  `dd81ed2b4bd16212bcfc816dc14de6d6d9ebabef089888c2ab36490ece2945aa`
- Verdict: **APPROVE，限上述精確 source/test SHA。**

## 已關閉 findings

| Finding | 狀態 | 實際封閉方式 |
|---|---|---|
| v4 authority 與 M2 source 未可信錨定 | CLOSED | 固定 authority、approval、signature、public key SHA；驗 Ed25519；重算 23-role source-set；驗 M2 與 primary-v3 cross-binding |
| resolved identity 與 parse path 分離 | CLOSED | module、兩輪 site TSV 與 assignments 全部使用已解析且鎖定的 input path |
| 全量 schema / enum / nonfinite 未驗 | CLOSED | 469,849 rows 全部解析並逐列重算 M1 stable membership |
| atomic publication 可能覆寫競爭目標 | CLOSED | Linux `renameat2(RENAME_NOREPLACE)`；目標已存在時 fail closed |
| E2E negative coverage 不足 | CLOSED 至核准門檻 | authority drift、strict-R1、nonfinite、existing output、concurrent target 均有拒絕測試 |

## 實際驗證

定向測試：

```bash
env -u PYTHONPATH PATH=/usr/bin:/bin PYTHONNOUSERSITE=1 \
  PYTHONHASHSEED=0 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
  MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 BLIS_NUM_THREADS=1 \
  /bip7_disk/liaoyoyo2001/miniconda3/envs/cnvtools/bin/python \
  -m pytest -q \
  research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/tests/test_audit_positional_singleton_methyl_multigroup.py
```

Reviewer 實際結果為 `15 passed`；真實全鏈 exit code `0`、`pass=true`。

正式輸出：

```text
/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/
20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/
positional_singleton_methyl_multigroup_audit_v1_source_attested/
```

- 目錄 mode `0555`
- 所有輸出 mode `0444`
- 469,849 全量 site rows 與 102,842 stable assignment keys 交叉驗證
- strict-R1 mismatches: `{}`
- singleton site audit: 50,432 rows
- M2 PASS case table: 30 rows
- 39 / 39 checks PASS

## 剩餘風險與處理

- Supplemental auditor 不在既有 v4 23-role source set；正式 singleton summary 已固定
  auditor SHA，後續 source-link receipt 必再綁定 final signed dataset/report receipts。
- 現有 before/after guard 是 path-based identity，不是全程持有同一 FD；惡意
  replace-and-restore 屬理論低風險。
- `renameat2` 是 Linux release scope；其他平台 fail closed。
- Signed data-release manifest 仍由後續 final dataset/report detached signatures 提供，
  本 supplemental receipt 不自行升級為新的 release authority。
