<!--
建立時間: 2026-07-17T14:18:00+08:00
目標: 獨立檢查 tumor-REF retrospective source attestation schema 1.2 與 final release builder 的契約一致性
處理範圍: v2 verifier、completion runner、final dataset builder、真實 snapshot / run manifest / pre-release probe
關聯檔案:
  - InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/scripts/verify_retrospective_running_source_identity_v2.py
  - InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/scripts/run_m2v5_recovered_completion_chain.sh
  - InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/scripts/build_all_ssnv_final_report_dataset.py
  - InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/results/tumor_ref_source_attestation_strict_repo_relative_preprobe.v1.json
-->

# Tumor-REF source attestation 與 release 接縫獨立稽核

## 稽核資訊

- 有效 reviewer：subagent `019f6eac-6ea8-7fd1-8271-87beb582dd72`（Fermat），唯讀。
- 稽核目標：確認執行中 snapshot creator、producer manifest、post-run verifier、relative script token，以及 completion runner / final builder 是否使用同一契約。
- 原始判定：無 High；2 Medium、2 Low。真實 probe 內部一致，但正式 receipt 尚未建立，且 builder 接受面需要收緊。
- 修正後判定：兩項程式契約問題已封閉；正式 receipt 仍依設計等待 completion runner 建立，不能以 pre-probe 取代。

## Findings 與處理

| 等級 | Finding | 影響 | 處理狀態 |
|---|---|---|---|
| Medium | final builder 只檢查 receipt checks 與 verifier 宣告身份，未逐欄重算 `command_binding` | standalone builder 可接受比 verifier 1.2 更寬的替代 receipt | **CLOSED**：builder 現在讀回 hash-bound snapshot，重驗 live command、manifest command、analyzer token與完整 binding |
| Medium | 正式 `post_run_source_identity.receipt.json` 尚不存在 | 尚不能稱 source-attested final release | **EXPECTED PENDING**：正式路徑只由 completion runner 原子建立；pre-probe 明確不可替代 |
| Low | verifier 接受任意 source-path suffix，例如 basename-only | 一般化契約可能誤認同名檔案 | **CLOSED**：relative token 必須精確等於 analyzer 相對 repo root 的完整路徑；mode=`repo_relative_exact` |
| Low | 舊 `/tmp` probe 記錄 verifier mode 0664，現在為 0444 | 舊 probe 的 stat metadata 過時 | **CLOSED / SUPERSEDED**：fresh probe 記錄 mode 0444、current SHA；舊 probe不作release evidence |

## 修正後驗證

輸入：

- `/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/tumor_ref_recovery_source_identity_v1/observed_during_execution.snapshot.json`
- `/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/all_ssnv_tumor_ref_controls_v2_prefix_recovered_seed_parallel/run_manifest.json`

命令：

```bash
/bip7_disk/liaoyoyo2001/miniconda3/envs/cnvtools/bin/python -B \
  research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/scripts/verify_retrospective_running_source_identity_v2.py verify \
  --snapshot /big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/tumor_ref_recovery_source_identity_v1/observed_during_execution.snapshot.json \
  --run-manifest /big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/all_ssnv_tumor_ref_controls_v2_prefix_recovered_seed_parallel/run_manifest.json \
  --output research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/results/tumor_ref_source_attestation_strict_repo_relative_preprobe.v1.json
```

實際結果：

- `pass=true`、schema `1.2.0`。
- `manifest_script_token_mode=repo_relative_exact`。
- manifest token精確為`research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/scripts/recover_all_ssnv_tumor_ref_controls.py`。
- v2 verifier SHA-256=`dc785d129b3eb3c34bf8bce12907d57f7069be5cf71fccbf6ae9b0307c9733ea`、mode=`0444`。
- verifier + final builder定向測試：`51 passed in 7.41s`。
- 故障注入已涵蓋錯 verifier、竄改 command binding、basename-only／短 suffix、`.`／`..`與source drift。

## Release 判定

目前可稱「strict repo-relative bounded retrospective attestation pre-probe PASS」，不可稱 final release PASS。只有 cooccurrence v5 完成後，由`run_m2v5_recovered_completion_chain.sh`在正式路徑建立receipt、再由final builder獨立重驗，才完成source-attested release gate。

## Bounded re-review

原 reviewer 於修補後重新執行四項契約 replay，結果為 High=0、Medium=0、Low=0，判定 **GO**：

- 完整repo-relative與absolute exact token通過；basename、短`scripts/<basename>`、`./`、`../`全部拒絕。
- Fresh pre-probe經builder實際讀回為`release_gate_pass=true`、`publishable_task_b_release=true`。
- Wrong verifier、receipt binding drift、manifest command drift、snapshot SHA drift均由production loader拒絕。
- 精確token定向pytest=`1 passed`；五個指定檔案在審查前後SHA-256完全相同，reviewer未修改檔案。

此GO只涵蓋bounded source-attestation code/release接縫；不代表尚在執行的cooccurrence或其後科學分析已完成。
