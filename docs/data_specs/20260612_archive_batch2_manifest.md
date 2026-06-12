<!--
建立時間: 2026-06-12
報告類型: Archive batch 2 manifest（舊日期 plain canonical run，供 AI session 理解確認）
任務類型: 輸出清理 batch 2（grep-check 後候選清單 + verdict）
狀態: ✅ 完成 — 11 SAFE 已 archive；6 B3 依賴用戶決定 KEEP（2026-06-12）
data_sources: 主回合 2026-06-12 ls/grep-check（存在性+mtime+大BAM+token邊界引用檢查）
build_branch: research/subclonal-reconstruction-202606
驗證方式: ls -dl mtime（不 du）；grep -rIlE token 邊界（避免誤匹配 _complete_matrix）scoped docs/+scripts/
-->

# Archive Batch 2 Manifest — 舊日期 plain canonical run

> 組織規則見 `InterSubMod/docs/data_specs/20260612_output_organization_guide_subclonal_01.md` §6；接續 batch 1（`20260612_archive_batch1_manifest.md`）。
> **狀態：PENDING — 已 grep-check，等用戶確認才 `mv`（不刪）。**

## L0 結論

- batch 2 = **17 個 plain run**（非-complete_matrix、無 _1/_2/test 後綴）。**全部 0 大 BAM**（省 inode 非大空間）。
- grep-check（token 邊界，避免誤匹配正典 `_complete_matrix`）：**11 SAFE archive** + **6 KEEP**（B3_paired_obs18 腳本硬依賴，與 batch 1 一致）。

## 11 SAFE archive（只有軟引用 = stale catalog / 歷史 doc，非腳本功能依賴）

| run | mtime | 軟引用 |
|---|---|---|
| `20260211_HCC1395_paired_full_full` | 03-14 | path_inventory.tsv（stale catalog）|
| `20260307_HCC1395_paired_full_full` | 03-14 | path_inventory + 03-14 標準 line 50（列表）+ 手冊 line 71（範例指令）|
| `20260314_HCC1395_paired_full_full`（plain）| 03-14 | 無 |
| `20260420_HCC1395_paired_full_full`（plain）| 04-20 | 無 |
| `20260212_HCC1937_paired_full_full` | 03-14 | path_inventory |
| `20260212_HCC1954_paired_full_full` | 03-14 | path_inventory |
| `20260212_H1437_paired_full_full` | 03-14 | path_inventory |
| `20260212_H2009_paired_full_full` | 03-14 | path_inventory |
| `20260212_COLO829_paired_full_full` | 03-14 | path_inventory |
| `20260211_HCC1395_DORADO_paired_full_full` | 03-14 | path_inventory |
| `20260307_HCC1395_DORADO_paired_full_full` | 03-14 | path_inventory |

> 軟引用非功能依賴：`path_inventory.tsv`(04-11) 已 stale/危險（memory `project_disk_cleanup_2026_06_01`，不維護）；0307 標準/手冊是 pre-complete_matrix 歷史列表+範例，archive 不破壞執行（current canonical = complete_matrix）。

## 6 KEEP（B3_paired_obs18.py 硬編碼輸入 line 77-82；與 batch 1 同一 phasing 對照）

`20260421_HCC1937` · `20260421_HCC1954` · `20260421_H1437` · `20260421_H2009` · `20260421_COLO829` · `20260420_HCC1395_DORADO`
> 理由同 batch 1 KEEP：B3 = obs18 LOH phasing paired 對照（降支撐但仍在新論文脊柱）；無大 BAM、留著成本低；HD-1 走 positive 可能重用。

## Archive 機制（同 batch 1）
`mv` 到 `_ARCHIVE_pending_cleanup_202606/canonical/{sample}/paired_full/`（同 fs 瞬間、保留結構可還原、不刪）；刪除留 ≥1 週 + 用戶確認（Hard Gate）。

## 執行紀錄

**2026-06-12 已 mv 11 SAFE run** → `_ARCHIVE_pending_cleanup_202606/canonical/{sample}/paired_full/`：
HCC1395 (0211/0307/0314-plain/0420-plain) · HCC1937/HCC1954/H1437/H2009/COLO829 (各 0212) · HCC1395_DORADO (0211/0307)。
**6 B3 依賴用戶決定 KEEP**（0420/0421_* + DORADO 0420）。

## Batch 2 結論

✅ **完成**：11 archived + 6 KEEP。**canonical paired_full 徹底清爽 — 每樣本剩 2 run：`complete_matrix`(正典) + B3-KEEP**。
archive 累計 16 run（batch1 5 + batch2 11，全 0 大 BAM）。
**還原**：`mv _ARCHIVE_pending_cleanup_202606/canonical/{sample}/paired_full/{run} canonical/.../`。
**刪除**：留 ≥1 週 + 用戶確認（Hard Gate）。
**剩 batch 3**：頂層 legacy pilot（見組織指南 §6）+ 各樣本 `paired_pileup/`/`to_pileup/` 模式（如需另起）。
