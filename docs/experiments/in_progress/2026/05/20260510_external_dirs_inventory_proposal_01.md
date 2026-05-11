---
title: 外部資料目錄盤整 — /big7_disk/liaoyoyo2001/data/ 與 big7_disk_output/
date: 2026-05-10
status: pending_user_review
type: inventory_proposal
classification: housekeeping_external_dirs
generated_by: AI scan (no deletion / move executed)
related:
  - InterSubMod/docs/experiments/in_progress/2026/05/20260510_cleanup_proposal_obsolete_files_01.md (InterSubMod/ repo 內 cleanup)
  - InterSubMod/docs/data_specs/20260414_output資料信任度與生命週期_01.md (信任度規範)
  - big7_disk_output/OUTPUT_INDEX.md (output 索引)
---

# 外部資料目錄盤整 — Inventory Proposal

> **Bottom line**：盤點 `/big7_disk/liaoyoyo2001/data/`（**412 GB**，輸入資料）與 `/big7_disk/liaoyoyo2001/big7_disk_output/`（11 個頂層目錄 + 兩個索引檔，總大小待用戶離線 du 補完）。**未執行任何刪除/搬移**。  
> 整體結論：(a) `data/` 三大子目錄（HCC1395_5kHz / PON / ref）皆 active 在用 — **不動**；只有空目錄 `data/bam/` 可清；  
> (b) `big7_disk_output/` 已有 `OUTPUT_INDEX.md` + 信任度 4 級規範，多個 SUPERSEDED archive 目錄可依 INDEX 判斷封存；  
> (c) **9 項候選需 user review**（含 2 個 SUPERSEDED archive 共 44 dirs、2 個 March pilot/runbook、2 個 KDE rerun 過程檔等）。

---

## §1 `/big7_disk/liaoyoyo2001/data/` — 412 GB 輸入資料

| 子目錄 | 大小 | mtime | 內容 | 引用方 | 建議 |
|---|---|---|---|---|---|
| `HCC1395_5kHz/` | **408 GB** | 2026-04-13 | `HCC1395.bam` (292 GB tumor) + `HCC1395BL_5kHz.tagged.bam` (145 GB normal) + .bai + `snv.vcf.gz` | `scripts/analysis/v5_*.py`、`scripts/analysis/v3f_ablation_caller_f1.sh`、`.claude/settings.local.json` | 🟢 **保留**（active 主資料） |
| `PON/` | 1.4 GB | 2026-04-03 | 4 個 PON VCF（1000G + CoLoRSdb + dbsnp + gnomad）+ tbi | reference data | 🟢 **保留**（reference） |
| `ref/` | 3.0 GB | 2026-04-03 | `GRCh38_no_alt_analysis_set.fasta` + .fai | reference data | 🟢 **保留**（reference） |
| `bam/` | **4.0 KB（空目錄）** | 2026-04-13 | 無檔案 | 無引用 | 🟡 **建議封存或刪除**（空目錄） |

**InterSubMod 對 data/ 的引用點**（4 處 — 證實 active 在用）：
- `InterSubMod/scripts/analysis/v5_whole_genome_per_site_audit.py`
- `InterSubMod/scripts/analysis/v5_purity06_somatic_alt_ratio.py`
- `InterSubMod/scripts/analysis/v3f_ablation_caller_f1.sh`
- `InterSubMod/.claude/settings.local.json`

**data/ 結論**：3/4 子目錄 active；只有 `bam/` 空目錄可清。

---

## §2 `/big7_disk/liaoyoyo2001/big7_disk_output/` — 結構盤點

> 本目錄已有完整 `OUTPUT_INDEX.md`（2026-04-14 建立，15 KB）+ `README.md`（1.8 KB）+ 「信任度 4 級規範」（CURRENT / PRE-FIX / SUPERSEDED / DEPRECATED）。  
> **本盤整文件不取代 OUTPUT_INDEX.md**；只新增從 v1.8 觀點的 housekeeping 評估。

### 2.1 頂層 11 個項目（無 size — 須離線 du 補；§5 提供 helper script）

| 項目 | mtime | 信任度 (per OUTPUT_INDEX) | 是否被引用 | 建議 |
|---|---|---|---|---|
| `canonical/` | 2026-04-22 | **PRE-FIX** | 7 樣本 × 3 模式 ISM baseline 19 runs | 🟢 **保留**（CURRENT_FOCUS 主軸） |
| `synthesis/` | 2026-04-20 | 混合 | 多 analysis script 引用 | 🟢 **保留**（含子目錄細審 §2.2） |
| `bip8_output_archive/` | 2026-04-20 | **SUPERSEDED** | 39 dirs 歷史資料 | 🟡 **可依 INDEX 判斷封存**（OUTPUT_INDEX §3 明寫「可依 INDEX 判斷刪除」） |
| `big8_output_archive/` | 2026-04-20 | **SUPERSEDED** | 5 dirs 歷史資料 | 🟡 **可依 INDEX 判斷封存**（同上） |
| `multilayer_hp_benchmark/` | 2026-04-03 | **PRE-FIX** | 仍被 self-phasing 驗證引用 | 🟢 **保留**（OUTPUT_INDEX §3 明寫「仍被引用」） |
| `v5_provenance_followup/` | 2026-05-06 | CURRENT | 最近 v5 audit；commit `388a437` E5/D4=DONE | 🟡 **若 v5 已 finalize → review 中間檔**（已在 cleanup proposal §3.3 列） |
| `InterSubMod_big7_runbook/` | 2026-03-19 | **CURRENT**（標籤）但 1.5 月未動 | 4 檔（遷移手冊 + 狀態 + scripts + manifests） | 🟡 **review 是否仍需 active**（migration 完成後可降為 SUPERSEDED） |
| `hcc1395_normal_pilot/` | 2026-04-19 | 未在 OUTPUT_INDEX 列 | 早期 normal pilot | 🟡 **review** |
| `hcc1395_normal_pilot_global/` | 2026-04-21 | 未在 OUTPUT_INDEX 列 | 同上 global 版 | 🟡 **review**（與 hcc1395_normal_pilot/ 是否重複/取代關係） |
| `kde_smoke_test/` | 2026-04-24 | 未在 OUTPUT_INDEX 列 | KDE smoke test | 🟡 **若 KDE smoke 已併入主分析 → 可封存** |
| `20260307_hcc1395_to_pilot_1` (symlink) | 2026-03-16 | symlink only | → `synthesis/research_rounds/legacy_partials/...` | 🟢 **保留**（symlink 不佔空間） |

### 2.2 `synthesis/` 子目錄（11 項）

| 子目錄 | mtime | 內容性質 | 建議 |
|---|---|---|---|
| `archive_synthesis_manifest.tsv` (file) | — | manifest | 🟢 保留 |
| `batch_runs/` | 2026-03-17 | 早期批次跑 | 🟡 **review**（1.5 月舊） |
| `concluded/` | 2026-04-14 | 已 concluded 研究 | 🟢 保留（concluded 收尾資料） |
| `final_closeout/` | 2026-03-16 | closeout 收尾 | 🟡 **review** |
| `final_closeout_debug/` | 2026-03-15 | debug 殘留 | 🟡 **建議封存**（debug 過程檔，1.5 月未動） |
| `kde_rerun_B_14combos/` | 2026-04-21 | KDE B 14 組合 | 🟡 **若 KDE 已 finalize → review** |
| `kde_rerun_pilot/` | 2026-04-19 | KDE pilot | 🟡 **若 KDE 已 finalize → review** |
| `master_experiment_matrix.tsv` (file) | — | master matrix | 🟢 保留 |
| `master_report.md` (file) | — | master report | 🟢 保留 |
| `master_run_manifest.tsv` (file) | — | master manifest | 🟢 保留 |
| `observation_workspaces/` | 2026-04-11 | observation O-series | 🟢 保留（O-series workspace） |
| `research_rounds/` | 2026-04-23 | research round 結果 | 🟢 保留（含 legacy_partials 是 symlink target） |

---

## §3 候選封存清單（彙整）

| Tier | 候選 | 預估價值 | 風險 |
|---|---|---|---|
| **Tier 1 — 確認可清** (small) | `data/bam/`（空目錄）| 釋放 inode | 🟢 無風險 |
| **Tier 2 — OUTPUT_INDEX 已標 SUPERSEDED**（大）| `bip8_output_archive/`（39 dirs）+ `big8_output_archive/`（5 dirs）| 大量空間 | 🟡 中（屬 audit history；建議**先離線 du 確認 size** 與 cross-check `bip8_output_archive/INDEX.md` 與 `big8_output_archive/INDEX.md` 內紀錄） |
| **Tier 3 — Workflow 收尾類，1.5 月未動**（中）| `synthesis/final_closeout_debug/`（debug 殘留）| 中等空間 | 🟡 中（debug 過程檔，但可能仍需 audit） |
| **Tier 4 — 需用戶判斷 finalize 狀態** | `synthesis/batch_runs/`（Mar 17 舊）/ `synthesis/final_closeout/`（Mar 16 舊）/ `synthesis/kde_rerun_*`（Apr 19-21）/ `kde_smoke_test/`（Apr 24）/ `hcc1395_normal_pilot{,_global}/`（Apr 19-21）/ `InterSubMod_big7_runbook/`（Mar 19）/ `v5_provenance_followup/T1_2_read_level_audit/vote_dump_*.tsv.gz` | 視 finalize 狀態 | 🟠 視 task 鏈 |
| **Tier 5 — 不動** | `data/HCC1395_5kHz/` + `data/PON/` + `data/ref/` + `canonical/` + `multilayer_hp_benchmark/`（被引用）+ `synthesis/concluded/`（concluded archive）+ `synthesis/observation_workspaces/`（O-series active）+ `synthesis/research_rounds/`（含 symlink target） | — | 🟢 active reference 必保 |

---

## §4 待用戶確認的 task 關聯（決定 Tier 2-4 取捨）

| # | 問題 | 決定哪些封存 |
|---|---|---|
| 1 | bip8_output_archive 39 dirs 是否仍需查閱？(per OUTPUT_INDEX 標 SUPERSEDED「可依 INDEX 判斷刪除」) | 若 INDEX 內紀錄完備 → 可整個搬到外部冷儲存 |
| 2 | big8_output_archive 同上問題 | 同上 |
| 3 | KDE 修補（plan v1.8 P-08 涵蓋）已 finalize 嗎？kde_rerun_B_14combos / kde_rerun_pilot / kde_smoke_test 是否已併入主結論？ | 若併入 → 全部 3 個可封存到冷儲存 |
| 4 | `InterSubMod_big7_runbook/` 遷移手冊（big7 從 big8 遷移）是否已完成？(Mar 19 後未動，OUTPUT_INDEX 標 CURRENT 但似已過完成期) | 若遷移已 done → 可降為 SUPERSEDED 並搬 archive |
| 5 | `hcc1395_normal_pilot{,_global}/`（4 月中旬）是否被 v1.8 之後的 cycle 取代？ | 若有 cycle 取代 → 封存 |
| 6 | v5 provenance（commit 388a437「E5/D4 狀態 DONE」）已 finalize 嗎？ | 若 done → vote_dump_*.tsv.gz 可封存（與 cleanup proposal §3.3 一致） |
| 7 | `synthesis/final_closeout_debug/` 是否仍需 audit？ | 若僅 debug 殘留 → 直接封存 |
| 8 | `synthesis/batch_runs/`（Mar 17 舊批次）是否被新 canonical/ 取代？ | 若被取代 → 封存 |
| 9 | `synthesis/final_closeout/`（Mar 16）內容是否與 INDEX 重複？ | 若重複 → 封存；否則保留 |

---

## §5 Helper Scripts（用戶確認後執行）

### 5.1 先跑 du 拿真實 size 做最終決定（離線可在背景跑）

```bash
# 跑離線 du；time-out 後可分次跑
nohup bash -c 'du -sh /big7_disk/liaoyoyo2001/big7_disk_output/*/ 2>/dev/null > /tmp/big7_output_sizes.txt' &
# 等 jobs done 後 cat /tmp/big7_output_sizes.txt 拿到完整 size
```

預估時間：含 bip8_output_archive（39 dirs）可能 5-30 分鐘（視 inode 數量）。

### 5.2 Tier 1 一定可以做（無風險）

```bash
# data/bam 空目錄可直接 rmdir（不需建 archive）
rmdir /big7_disk/liaoyoyo2001/data/bam
```

### 5.3 Tier 2 SUPERSEDED archive 範例（**du 拿到 size 後決定**）

```bash
# 假設用戶決定把整個 archive 搬到冷儲存（外部硬碟 / NAS / 壓縮）
ARCHIVE_DEST="/big7_disk/liaoyoyo2001/big7_disk_output_cold_archive"
mkdir -p "$ARCHIVE_DEST"
mv /big7_disk/liaoyoyo2001/big7_disk_output/bip8_output_archive "$ARCHIVE_DEST/"
mv /big7_disk/liaoyoyo2001/big7_disk_output/big8_output_archive "$ARCHIVE_DEST/"
# 同步更新 OUTPUT_INDEX.md 把這兩項標為 ARCHIVED_TO_COLD
```

### 5.4 Tier 3-4 須個別評估，不寫範本（避免誤殺）

---

## §6 與既有規範的整合

本盤整與 InterSubMod 既有 3 個規範交叉對齊：

| 既有規範 | 本盤整用法 |
|---|---|
| `InterSubMod/docs/data_specs/20260414_output資料信任度與生命週期_01.md` | 信任度 4 級 → 本盤整套用作為 Tier 分類基礎 |
| `big7_disk_output/OUTPUT_INDEX.md` | 已有完整 INDEX 與信任度標註 → 本盤整在其上新增「mtime 觀察 + 引用 cross-ref」 |
| `InterSubMod/docs/experiments/in_progress/2026/05/20260510_cleanup_proposal_obsolete_files_01.md`（v1.8 後續） | InterSubMod/ repo 內的 housekeeping → 與本外部目錄盤整互補不重疊 |

**建議**：用戶先跑 §5.1 du → 拿到 size 後針對 §3 Tier 2-4 候選逐項決定 → 執行 §5.2-5.4 對應 helper。

---

## §7 結論摘要

> 1. **`data/`** 412 GB：3/4 子目錄是 active 主資料（HCC1395 408GB BAM + PON 1.4G + ref 3G）— **保留**。空目錄 `data/bam/` 可 rmdir。
> 2. **`big7_disk_output/`** 11 頂層項目：已有完整 `OUTPUT_INDEX.md` 與信任度規範；本盤整新增 housekeeping 視角。
> 3. **大可釋放空間** = bip8_output_archive (39 dirs SUPERSEDED) + big8_output_archive (5 dirs SUPERSEDED) — OUTPUT_INDEX 已標可依 INDEX 判斷刪除。
> 4. **9 項候選 (§4)** 待用戶確認 task finalize 狀態後決定 — 涵蓋 KDE rerun、normal pilot、runbook、v5 vote_dump 等。
> 5. **未執行任何操作** — 用戶 review 後可參考 §5 helper script 執行。

---

## §8 References

- v1.8 retro：`InterSubMod/docs/experiments/in_progress/2026/05/20260510_v1.8_implementation_retro_01.md`
- InterSubMod cleanup proposal：`InterSubMod/docs/experiments/in_progress/2026/05/20260510_cleanup_proposal_obsolete_files_01.md`
- 資料信任度規範：`InterSubMod/docs/data_specs/20260414_output資料信任度與生命週期_01.md`
- 既有 OUTPUT_INDEX：`big7_disk_output/OUTPUT_INDEX.md`（不在 InterSubMod/ repo 內）
