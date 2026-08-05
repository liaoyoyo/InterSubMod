<!--
建立時間: 2026-07-19T11:56:20+08:00
目標: 保存 cooccurrence v7 正式 runner 與 analyzer canonical command 不一致的 fail-closed 證據
處理範圍: 正式 producer 啟動階段；尚未建立任何科學結果或 release receipt
關聯檔案:
  - InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/scripts/run_cooccurrence_v6_source_locked.sh
  - InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/scripts/analyze_methyl_ssnv_cooccurrence.py
  - InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/results/cooccurrence_task_contract_preflight.v8_summary_hotfix_full_runtime.json
-->

# Cooccurrence v7 canonical command / pycache 參數不一致失敗紀錄

## 判定

- 狀態：**ABORTED / EXCLUDED / FAIL-CLOSED**。
- 失敗階段：analyzer 建立 output directory 或任何科學結果之前。
- 程序 exit code：`1`。
- 科學輸出目錄：不存在。
- `run_receipt.json`、`summary.json`、site/pair TSV、`release_receipt.json`：全部不存在。
- 已完成的 preflight v8 不受影響；其 receipt 維持 mode `0444`、`pass=true`。

## 輸入

- Manifest：
  `InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/results/all_ssnv_input_manifest.json`
- Stable assignments：
  `/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/all_ssnv_focal_alt_multigroup_v10_source_locked_thread_pinned_recovered_full/all_ssnv_stable_multigroup_read_assignments.jsonl.gz`
- Site results：
  `/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/all_ssnv_focal_alt_multigroup_v10_source_locked_thread_pinned_recovered_full/all_ssnv_site_results.tsv.gz`
- Preflight v8：
  `InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/results/cooccurrence_task_contract_preflight.v8_summary_hotfix_full_runtime.json`
  - SHA-256：`c09861197ec57c6c00845f5c1a684e8f98d43df21a12cdbe261c05d1eb5e61c7`
  - tasks：`102,842 / 102,842`
  - pass：`true`

## 執行命令

```bash
/usr/bin/bash \
  /big7_disk/liaoyoyo2001/InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/scripts/run_cooccurrence_v6_source_locked.sh
```

Runner 實際啟動 producer 的命令前綴：

```text
/bip7_disk/liaoyoyo2001/miniconda3/envs/cnvtools/bin/python
-I
-X
pycache_prefix=/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/.python_cache_cooccurrence_v7_bound_bootstrap
.../analyze_methyl_ssnv_cooccurrence.py
```

Analyzer `canonical_producer_command()` 允許的前綴：

```text
/bip7_disk/liaoyoyo2001/miniconda3/envs/cnvtools/bin/python
-I
.../analyze_methyl_ssnv_cooccurrence.py
```

## 實際輸出片段

```text
ValueError: Producer process command must match the canonical Task-B command exactly
FAILED exit=1 line=265 command="${PYTHON}" -I -X "pycache_prefix=${PYTHON_CACHE_ROOT}" ...
```

完整 log：

`InterSubMod/research/20260715_all_ssnv_focal_alt_multigroup_cooccurrence_validation/audit_notes/failed_formal_runs/20260719_v7_canonical_command_pycache_mismatch/methyl_ssnv_cooccurrence_v7_m2v5_raw_identity_contract_source_locked_summary_hotfix.log`

- Log size：`4,872` bytes
- Log SHA-256：`54dbbdd35d2b0130df64fb79b1856e68d3c87039362b13b58b5741d05f845cf2`
- Python cache：同一 archive directory 下
  `.python_cache_cooccurrence_v7_bound_bootstrap/`，約 `16 MiB`
- 原 workspace log/cache 路徑已清空，未刪除或覆寫任何證據。

## 根因

受保護 runner 強制所有正式 Python 呼叫使用 isolated mode 與獨立
`pycache_prefix`；受保護 analyzer 則直接讀取 `/proc/self/cmdline`，要求實際
argv 與 `canonical_producer_command()` 完全相等。兩邊各自有單元測試，但缺少
跨檔整合測試，沒有在 source authorization 前發現 `-X pycache_prefix=...`
token 差異。

## 影響

- 沒有任何 cooccurrence 科學數字由此 attempt 產生或引用。
- M1/M2、positional-singleton、preflight v8 與 LongPhase-S 輸入結論不變。
- Cooccurrence L4/L5、final dataset、HTML 與 result/report signatures 仍為
  `NOT_COMPLETE`。

## 修復 gate

1. Runner 實際 argv 與 analyzer canonical builder 必須由同一契約產生或由
   cross-file test 精確比較。
2. 保留 `-I`、受控 pycache 與 `/proc/self/cmdline` exact-match 的安全目的。
3. 任何受保護 source 修改後，舊 source authority 不得沿用；需重新跑 canonical
   tests、內外部 review、source authority 與受影響的 preflight。
4. 不允許直接執行 analyzer 後手動拼接 receipt 來宣稱 runner formal PASS。

