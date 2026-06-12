<!--
建立時間: 2026-06-12
報告類型: 輸出組織與生命週期指南 — Subclonal Reconstruction 新主軸
任務類型: 輸出層整理規劃（讓後續 AI 清楚找到/確認/使用輸出，並識別 archive 候選）
狀態: current（所有結構/大小/run 數經 2026-06-12 安全 stat 驗證，未 du 遞迴）
data_sources: 主回合 2026-06-12 ls -d/ls -dl/ls -la 安全 stat（run 數+mtime+BAM 大小）+ docs/data_specs/20260414_output資料信任度與生命週期_01.md + docs/data_specs/20260601_cleanup_audit/
build_branch: research/subclonal-reconstruction-202606
驗證方式: 禁用 du/find 遞迴（ISM region >1萬子目錄會卡）；改 ls 淺掃 + mtime + ls 大檔 stat
-->

# 輸出組織與生命週期指南 — Subclonal Reconstruction 新主軸

> **本檔職責**：(1) 確認專案輸出特性；(2) 定義新主軸「最終主輸出 + 每次比較測試輸出」放哪；(3) 提供可尋性 + 重用舊檔 provenance；(4) 列 archive/清理候選 + **防卡安全檢查法**。
> **建在既有標準上不重複**：信任度 4 級 → `20260414_output資料信任度與生命週期_01.md`；已做的清理 → `20260601_cleanup_audit/`（回收 623 GiB）。

---

## L0 一眼結論

1. **輸出分三層放**：① **正典輸入**（`big7_disk_output/canonical/{sample}/{mode}/{date}_*_complete_matrix/`，7 個 tagged BAM 共 ~1.7 TB，**必留唯讀**）② **新主軸最終主輸出**（論文 figures/tables → `research/subclonal_reconstruction/final/`）③ **每次比較測試輸出**（per-cycle benchmark → `research/subclonal_reconstruction/cycles/{date}_{topic}/`）。
2. 🔴 **防卡鐵則**：ISM 輸出（`intersubmod_tp/filtered_snv_tp/`）有**上萬 region 子目錄** → **禁 `du`/`find .` 遞迴**（會卡死）。只用 `ls -d` 淺掃 / `ls -dl` 看 mtime / `ls -la` stat 大檔 / 讀 manifest。
3. **archive 候選（無大 BAM，省 inode 非大空間）**：**23 個冗餘 canonical run**（`_1`/`_2`/`smokecheck`/`parallel_test`/舊日期非-complete_matrix）+ **6 個頂層 legacy pilot**。**刪除=Hard Gate 須你確認**；本檔只列候選 + 提議搬 archive。

---

## §1 專案輸出特性（confirmed 2026-06-12）

**輸出根**：`/big7_disk/liaoyoyo2001/big7_disk_output/`（**非** repo 內；repo `InterSubMod/output/` 是 symlink）。big7 **98% 滿**（剩 1.1 TB）。

**結構文法**：`canonical/{sample}/{mode}/{date}_{sample}_{mode}_{suffix}/`
- `{sample}`：HCC1395 / HCC1937 / HCC1954 / H1437 / H2009 / COLO829（+ HCC1395_DORADO basecaller 對照）
- `{mode}`：`paired_full`（新主軸用）/ `paired_pileup` / `to_pileup`
- `{suffix}`：`complete_matrix`=**正典**（每樣本剛好 1 個）；其餘 = 冗餘/測試

**每個 complete_matrix run 的輸出單元**（= 一次完整 ISM 分析）：
| 物件 | 內容 | 大小/特性 |
|---|---|---|
| `longphase_s/{sample}_tagged.bam` | somatic-haplotag BAM（**大檔**）| 88–416 GB |
| `intersubmod_tp/filtered_snv_tp/{24 chr}/{region}/` | ISM TP 逐位點輸出 | 🔴 **上萬 region 子目錄**（卡住源）|
| `intersubmod_fp/` | ISM FP 輸出（判別對照另一半）| 同上 |
| `metrics.json` · `plots/` · `benchmark_comparison.{md,tsv}` · `round_summary.md` · `run_context.json` · `*.log` | 該 run 的彙總/圖/log | 小檔 |

**7 個正典 tagged BAM（= 新主軸 ISM 輸入，必留）**：
| 樣本 | tagged BAM | 大小 |
|---|---|---|
| HCC1937 | 20260315_*_complete_matrix | 416.5 GB |
| H2009 | 20260315 | 289.4 GB |
| HCC1395 | **20260314** | 259.2 GB |
| HCC1954 | 20260315 | 224.3 GB |
| HCC1395_DORADO | 20260315 | 222.1 GB |
| H1437 | 20260315 | 214.3 GB |
| COLO829 | 20260315 | 88.8 GB |
| **合計** | | **~1,714 GB（1.7 TB）** |

> **意涵**：big7 98% 主因 = 這 1.7 TB 正典 BAM（必留）+ 冗餘 run 的 ISM 小檔（高 inode）+ legacy。**冗餘 run 無大 BAM**（B 驗證：BAM 只出現在 complete_matrix）→ 清冗餘省的是檔案數/inode 與中量空間，非大空間。

---

## §2 信任度分類（verified；建在 04-14 4 級標準上）

| 等級 | 內容（verified by mtime+naming，未 du）| 動作 |
|---|---|---|
| **CURRENT 正典** | 7 個 `*_complete_matrix` run（mtime 2026-03-16）+ 其 tagged BAM | **KEEP 唯讀**（新主軸輸入）|
| **SUPERSEDED 冗餘** | **23 個** 非-complete_matrix run（見下清單）| archive 候選 |
| **legacy pilot** | 頂層 6 個舊 pilot/test 目錄 | archive 候選 |
| **已 archive** | `_ARCHIVE_2026_06_01/` · `big8_output_archive/` · `bip8_output_archive/` | 已歸檔（06-01）|

**23 個冗餘 canonical run 清單**（by sample；正典 complete_matrix 已排除）：
- **HCC1395 (8)**：20260211 · 20260307 · 20260314 · 20260314_1 · 20260314_2 · 20260420 · 20260420_1 · 20260420_2
- **HCC1937 (2)**：20260212 · 20260421
- **HCC1954 (2)**：20260212 · 20260421
- **H1437 (2)**：20260212 · 20260421
- **H2009 (2)**：20260212 · 20260421
- **COLO829 (4)**：20260212 · **20260314_smokecheck** · **20260315_parallel_test** · 20260421
- **HCC1395_DORADO (3)**：20260211 · 20260307 · 20260420

**頂層 legacy pilot（6，archive 候選）**：`20260307_hcc1395_to_pilot_1/`（03-16）· `multilayer_hp_benchmark/`（04-03）· `hcc1395_normal_pilot/`（04-19）· `hcc1395_normal_pilot_global/`（04-21）· `kde_smoke_test/`（04-24）· `v5_provenance_followup/`（05-06）。
> ⚠ 搬移前 SAFE check：`grep -rl "<run 名>" docs/ research/ 2>/dev/null`（scope 到輕目錄）確認無引用；尤其 `hcc1395_normal_pilot*` 可能與 V10 normal 分析相關，需先確認。

---

## §3 新主軸輸出組織 scheme（核心 — 解「最終主輸出 + 每次測試輸出 放哪不混亂」）

**放哪決策表**：
| 輸出類型 | 放哪 | 命名 | 範例 |
|---|---|---|---|
| **正典輸入**（ISM run）| `big7_disk_output/canonical/{sample}/paired_full/{date}_*_complete_matrix/` | 既有 | 唯讀，勿動 |
| **最終主輸出**（論文 figure/table/result）| `research/subclonal_reconstruction/final/{figures,tables,results}/` | `Fig{N}_{desc}` `Tab{N}_{desc}` | `final/figures/Fig1_reconstruction_demo.png` |
| **每次比較測試輸出**（per-cycle benchmark）| `research/subclonal_reconstruction/cycles/{YYYYMMDD}_{topic}/` | 內 `00_README.md`+`data/`+`figures/`+`results.json` | `cycles/20260615_GA_cross5sample/` |
| **真值/可重跑腳本**（已存）| `docs/experiments/in_progress/2026/05/20260531_methyl_phasing_A0_assets/` | 既有 | VERIFIED_RESULTS=SoT |
| **per-cycle 索引**（找得到）| `research/subclonal_reconstruction/00_INDEX.md` | L0+L1 分層 | 每 cycle 一行 hook + verdict + path |

**規則**（避免混亂）：
1. **一 cycle 一資料夾**（`cycles/{date}_{topic}/`），內含 `00_README.md`（目的/輸入/輸出/verdict）+ `data/` + `figures/` + `results.json`。**禁散落**。
2. **最終 vs 測試分離**：升上論文的圖表才進 `final/`；探索/比較留 `cycles/`。
3. **每 cycle 完成更新 `00_INDEX.md`**（L0 一行：date · topic · verdict · tier · path）→ 後續 AI 一眼知道有哪些、哪個能用。
4. **大中間檔（BAM/大 TSV）寫 `cycles/.../data/` 並 gitignore**；結論數字落 `results.json` + README。
5. **數字誠信**：results.json 為真值來源，README 引用前先 Read（§13.0）。

---

## §4 可尋性（後續 AI 怎麼找）

| 要找什麼 | 去哪 |
|---|---|
| 新主軸現況/下一步 | `docs/CURRENT_FOCUS.md` pinned 區塊 |
| 新主軸資料/知識/流程 | `docs/reports/research_landscape/20260611_subclonal_reconstruction_LAUNCH_READINESS_01.md` |
| 正典輸入路徑 | 本檔 §1 + `20260612_external_data_dependencies_01.md`（normal BAM）|
| 每次 cycle 輸出 | `research/subclonal_reconstruction/00_INDEX.md`（新建）|
| 真值數字 | A0_assets `VERIFIED_RESULTS.md` |
| 舊輸出信任度/已清理 | `20260414_output資料信任度與生命週期_01.md` + `20260601_cleanup_audit/MASTER_CLEANUP_INDEX.md` |
> ⚠ `big7_disk_output/OUTPUT_INDEX.md`（04-14）**已 stale 危險**（memory `project_disk_cleanup_2026_06_01`：早期 pilot 誤標 SAFE_DELETE 實為 LIVE）→ **勿信其 SAFE_DELETE 標記**，以本檔 §2 + 06-01 ledger 為準。

---

## §5 重用舊檔 provenance（新主軸會用、已驗證可信）

| 舊檔/輸出 | 新主軸用途 | 驗證狀態 |
|---|---|---|
| 7 個 `*_complete_matrix/` tagged BAM + ISM TP/FP | G-A 跨樣本 V10/V11c 輸入 | ✅ 2026-06-12 ls 驗存在（§1）|
| A0_assets（VERIFIED_RESULTS + ~150 腳本）| V1-V12 真值 + 可重跑（改 --bam 跑他樣本）| ✅ SoT 完整 |
| big8 normal BAM（5/6 有甲基）| V10 copy-clean 對照 | ✅ 契約 `20260612_external_data_dependencies_01.md` |
| `docs/method_comparison/20260609_ism_vs_external.../` | ISM vs 外部法定位 | ✅ Phase B 待核准 |
| `synthesis/master_run_manifest.tsv` | run 清單對照 | 🟡 2 月舊，參考非權威 |
> 其餘 `research/` 38 個舊 dir 多為**舊軸探索**（filter-DEAD/LOH/HPFine 等）→ 新主軸**不直接用**，留 provenance（見 MEMORY Concluded 區）。

---

## §6 archive / 清理計畫（候選 + 安全法；🔴 刪除須你確認）

**提議 archive 目標**：`big7_disk_output/_ARCHIVE_pending_cleanup_202606/`（同 fs `mv` = 瞬間、不佔空間、不 du）。

**分批（每批先 SAFE grep-check 無引用）**：
- **批 1（最安全）**：COLO829 `smokecheck`+`parallel_test`（明確測試）+ 所有 `_1`/`_2` 後綴（明確重試）。
- **批 2**：舊日期非-complete_matrix run（20260211/0212/0307）。
- **批 3**：頂層 legacy pilot（先確認 `hcc1395_normal_pilot*` 非 V10 來源）。

**安全法（given 上萬 region 卡住風險）**：
1. **判舊**：`ls -dl --time-style=long-iso <dir>`（mtime，不走子目錄）。
2. **判引用**：`grep -rl "<run名>" docs/ research/ scripts/`（scope 輕目錄，非 repo root）。
3. **搬移**：`mv <run> _ARCHIVE_pending_cleanup_202606/`（同 fs 瞬間；**禁先 du**）。
4. **記錄**：搬移前 append `{run, mtime, 理由}` 到 archive manifest（可追溯）。
5. **刪除**：archive 後**留存 ≥1 週**，確認新主軸不需 → **你確認**才 `rm`（Hard Gate）。
> 與 06-01 audit 關係：06-01 已清頂層大 BAM（回收 623 GiB）；本批是**canonical 內冗餘 run**（06-01 未涵蓋 canonical），多為 ISM 小檔（省 inode + 中量空間）。

---

## §7 防卡安全檢查法（鐵則，因 ISM region >1萬子目錄）

| 禁用（會卡）| 改用（安全）|
|---|---|
| `du -sh <ISM run>` | `ls -la <run>/longphase_s/*.bam`（只 stat 大檔）|
| `find <run> -name ...` | `ls -d <run>/*/`（淺掃一層）|
| `du */` / `find .` | `ls -dl --time-style=long-iso`（mtime 判舊）|
| 走 `intersubmod_tp/filtered_snv_tp/` | 讀 `metrics.json` / `round_summary.md`（彙總）|
> 需要 region 級細節時，先 `ls intersubmod_tp/filtered_snv_tp/ | head` 確認層級，再針對**單一 chr/region** 精準 ls，**永不全掃**。

---

> **使用方式**：新主軸跑分析前讀 §1（正典路徑）+ §3（輸出放哪）；找舊輸出讀 §4-§5；清理讀 §6-§7。archive/刪除一律先 SAFE check + 你確認。
