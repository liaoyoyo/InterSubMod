<!--
建立時間: 2026-06-12
報告類型: Archive batch 1 manifest（冗餘 canonical run 搬移紀錄，供其他 AI session 理解確認）
任務類型: 輸出清理 batch 1（grep-check 後候選清單 + verdict）
狀態: 5 安全 run 已 archive（2026-06-12 mv 完成）；HOLD 1 個待 B3 決策
data_sources: 主回合 2026-06-12 ls/grep-check（存在性+mtime+大BAM數+引用檢查）
build_branch: research/subclonal-reconstruction-202606
驗證方式: ls -dl mtime（不 du）；grep -rIl scoped docs/+scripts/（§12）；manifest 比對
-->

# Archive Batch 1 Manifest — 冗餘 canonical run（最安全批）

> 組織規則見 `InterSubMod/docs/data_specs/20260612_output_organization_guide_subclonal_01.md` §6。
> **狀態：PENDING — 已 grep-check，等用戶確認才 `mv`（不刪）。**

## L0 結論

- batch 1 = **6 個最明確的測試/重試 run**（`smokecheck`/`parallel_test`/`_1`/`_2`）。
- **全部 0 個大 BAM** → 搬移省的是 ISM 小檔/inode，**不省大空間**（1.7TB 正典 BAM 不在此批）。
- grep-check：**5 個安全**（無功能引用）；**1 個 hold**（`20260420_HCC1395_paired_full_full_2` 被 B3_obs18 腳本硬編碼當輸入）。
- OUTPUT_INDEX / master_run_manifest **未列**這些 run → 搬移不影響索引。

## 候選清單 + verdict（grep-check 驗證）

| run（canonical/.../paired_full/）| mtime | 大BAM | 引用檢查 | verdict |
|---|---|:---:|---|---|
| `20260314_HCC1395_paired_full_full_1` | 03-14 | 0 | 僅命名 catalog（list，非功能）| ✅ **archive** |
| `20260314_HCC1395_paired_full_full_2` | 03-14 | 0 | 無功能引用 | ✅ **archive** |
| `20260420_HCC1395_paired_full_full_1` | 04-20 | 0 | 無功能引用 | ✅ **archive** |
| `20260314_COLO829_paired_full_full_smokecheck` | 03-14 | 0 | 無（測試 run）| ✅ **archive** |
| `20260315_COLO829_paired_full_full_parallel_test` | 03-15 | 0 | 無（測試 run）| ✅ **archive** |
| `20260420_HCC1395_paired_full_full_2` | 04-20 | 0 | 🔴 **B3_obs18 腳本 line 76 硬編碼輸入** + 報告 line 41 | ⏸ **HOLD** |

**HOLD 理由**：`scripts/analysis/20260423_B3_paired_obs18.py:76` = `"HCC1395": "20260420_HCC1395_paired_full_full_2"`；`docs/experiments/in_progress/2026/04/20260423_B3_Paired_obs18_Control_01.md:41` 引用。B3_obs18 = 舊軸 Paired obs18 Control 實驗。**搬移前需確認 B3_obs18 已 concluded 且不重跑**（concluded → 可連同標記後一起 archive；若可能重跑 → 留）。

## Archive 機制（確認後執行）

1. 目標：`/big7_disk/liaoyoyo2001/big7_disk_output/_ARCHIVE_pending_cleanup_202606/`
2. `mv`（**同 fs 瞬間、不佔空間、不 du**）5 個安全 run。
3. 本 manifest 記錄搬移時間 + 來源路徑（可追溯/可還原）。
4. **刪除 = Hard Gate**：archive 後留存 ≥1 週，確認新主軸不需 → 用戶確認才 `rm`。

## 後續批次（見組織指南 §6）

- 批 2：舊日期非-complete_matrix run（20260211/0212/0307）— 需逐一 grep-check。
- 批 3：頂層 legacy pilot（先確認 `hcc1395_normal_pilot*` 非 V10 來源）。

## 執行紀錄

**2026-06-12 已 mv 5 個安全 run** → `big7_disk_output/_ARCHIVE_pending_cleanup_202606/canonical/{sample}/paired_full/`（保留結構供還原；同 fs 瞬間）：
| run | 來源 → 目標 | 狀態 |
|---|---|---|
| `20260314_HCC1395_paired_full_full_1` | canonical → _ARCHIVE_pending_cleanup_202606 | ✅ moved |
| `20260314_HCC1395_paired_full_full_2` | 同上 | ✅ moved |
| `20260420_HCC1395_paired_full_full_1` | 同上 | ✅ moved |
| `20260314_COLO829_paired_full_full_smokecheck` | 同上 | ✅ moved |
| `20260315_COLO829_paired_full_full_parallel_test` | 同上 | ✅ moved |

**還原法**：`mv _ARCHIVE_pending_cleanup_202606/canonical/{sample}/paired_full/{run} canonical/{sample}/paired_full/`。
**刪除**：留 ≥1 週確認新主軸不需 → 用戶確認才 `rm`（Hard Gate）。

## HOLD 決策依據（B3 查證，2026-06-12）

`20260420_HCC1395_paired_full_full_2` 被 `scripts/analysis/20260423_B3_paired_obs18.py:76` 硬編碼當輸入。查 `experiments/INDEX.md`：
- B3_paired_obs18 = **obs18 LOH-NG2 phasing 分析的 paired 對照**（INDEX line 251/296）。obs18 = G6 LOH-constrained phasing 主軸工作（NG=2 same-hap、Wilcoxon p=0.0078），**已降支撐但仍是新論文的 phasing 脊柱材料**。
- 該 run **無大 BAM**（0）→ 留著幾乎不佔空間。
- **建議：KEEP**（留著成本低 + 是仍在論文裡的 phasing 脊柱的 live 依賴；若 HD-1 走 phasing positive 路徑可能重用）。用戶最終決定。
