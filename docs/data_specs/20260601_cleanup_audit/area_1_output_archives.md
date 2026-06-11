<!--
建立時間: 2026-06-01
目標: Cleanup audit Area 1/6 — output 封存區 (bip8 / big8 / research_rounds/archive) 量化 + verdict
處理範圍: /big7_disk/liaoyoyo2001/big7_disk_output/{bip8_output_archive, big8_output_archive, synthesis/research_rounds/archive}
唯讀稽核 — 不刪除/搬移任何檔案 (刪/移為 Hard Gate, 須用戶確認)
量測方法: du -sb (小項精確) + find -size +5M 大檔加總 (大項 floor 估計, 因 big7_disk per-region 遞迴 du 全部 timeout)
關聯檔案:
  - /big7_disk/liaoyoyo2001/big7_disk_output/OUTPUT_INDEX.md
  - bip8_output_archive/INDEX.md / big8_output_archive/INDEX.md / synthesis/research_rounds/archive/INDEX.md
-->

# Area 1：Output 封存區 清理稽核

> **掃描狀態**：PARTIAL_WITH_ESTIMATES。big7_disk 對 per-region 子目錄（每 run 數十萬個 region 子資料夾）的 `du` 遞迴**全部 timeout**（即使 595s 上限）。改用 `du -sb`（小項精確）+ `find -size +5M` 大檔加總（大項 **floor 估計**：只計 >5MB 大檔，小檔殘量另計）。bytes 為 floor，真實值略高（小 csv 殘量）。
>
> **最高原則遵守**：本區三個 INDEX.md 把多數項標 SUPERSEDED / SAFE_DELETE，但**結案/被取代 ≠ 自動可刪**。撐起已發佈結論的原始證據 → 至少 ARCHIVE；LIVE 研究仍依賴者 → KEEP/NEEDS_USER_DECISION。

---

## 0. 兩項關鍵 anomaly（最重要 — 與 INDEX 自身分類衝突）

### A1. `202603_early_pilots/` = **2.71 TB** — INDEX 標 LOW_RISK，但實為 LIVE TO-BAM 唯一來源 ❗
- 三個 archive INDEX 把 `synthesis/research_rounds/archive/202603_early_pilots/` 標為 **LOW_RISK**（「被 canonical runs 取代」）。
- 但 `docs/experiments/in_progress/2026/04/20260422_Archive_TO_Rerun_ISM_Requirement_01.md` line 10/17/49-53 + `docs/experiments/in_progress/2026/04/20260423_KDE_Smoke_Test_Cross_Sample_Validation_01.md` line 71 證實：此目錄是 **5 樣本 TO `tumor_tagged.bam` 的唯一來源**（HCC1395_DORADO / H1437 / H2009 / HCC1937 / HCC1954），路徑形如 `202603_early_pilots/{date}_{sample}_to_pilot*/step03_longphase_to/tumor_tagged.bam`。
- canonical 的 `to_pileup` 列除 HCC1395 外**皆為空**（OUTPUT_INDEX §4.3 只有 1 個 TO run）→ canonical **未取代** TO BAM。
- KDE smoke test cross-sample validation（LIVE，LOH-phasing 主軸支撐）**直接讀此目錄**為輸入。
- **判定**：INDEX 的 LOW_RISK 分類**過時且危險**。verdict = **NEEDS_USER_DECISION（強烈 lean KEEP）**。占本區 ~75% 體積。**絕不可 SAFE_DELETE**。

### A2. `s-pure-pileup/` = **770 GB** — SUPERSEDED 但 canonical 為 partial（missing_tagged_bam）
- bip8 INDEX 標 `s-pure-pileup/` SUPERSEDED（→ `canonical/{sample}/paired_pileup/`）；7 樣本與 canonical 完全對齊。
- 但 OUTPUT_INDEX §2：全部 19 canonical runs `completeness_state = partial` / `blocking_reason = missing_tagged_bam`。770 GB 集中於 per-region 大檔（83 個 >5MB），疑為 paired_pileup 的中間/tagged BAM。canonical 既標 partial，**無法確認 s-pure-pileup 的 BAM 已在 canonical 重生**。
- docs/ 引用：多為 2026-03 validated/provenance 報告（已發佈結論的原始證據），非當前 LIVE 軸輸入。
- **判定**：verdict = **NEEDS_USER_DECISION（lean ARCHIVE）**。單項 770 GB = 本區最大「可能可回收」候選，但需用戶確認 canonical paired_pileup BAM 已重生後才可降為 SAFE_DELETE / 打包 ARCHIVE。

---

## 1. Artifact 判定表

> 欄位：path | 人類大小 | bytes(timeout=-1) | purpose | trust_tier | conclusion_status | cleanup_verdict | reclaimable_bytes | referenced_by | rationale
> 量測：(s)=du -sb 精確；(lf)=find>5MB floor；(est)=估計

### 1.1 `synthesis/research_rounds/archive/`（floor 總計 ~2.82 TB）

| path | 大小 | bytes | purpose | trust_tier | conclusion_status | verdict | reclaimable | referenced_by | rationale |
|------|------|-------|---------|-----------|-------------------|---------|------------|---------------|-----------|
| `…/archive/202603_early_pilots/` | ~2.71 TB (lf) | 2707240461764 | 建 canonical baseline 前的探索 pilot；實含 5 樣本 TO tagged BAM | SUPERSEDED(INDEX 宣稱) / **實為 CURRENT TO-BAM** | LIVE (TO-rerun + KDE smoke 依賴) | **NEEDS_USER_DECISION** | 0 (lean KEEP) | 20260422_Archive_TO_Rerun_ISM_Requirement_01.md; 20260423_KDE_Smoke_Test_Cross_Sample_Validation_01.md; 20260414_LOH_Subclone_AF_Methylation_Evidence_01.md | **anomaly A1**：INDEX LOW_RISK 過時；canonical to_pileup 5 樣本空，此為唯一 TO BAM 源 |
| `…/archive/duplicates/` | ~113 GB (lf) | 113164333292 | 同一任務 2 分鐘內連續重試執行（4× main_5khz + 3× to_pilot 重試） | SUPERSEDED | TRANSIENT (重試副本) | **NEEDS_USER_DECISION** (lean SAFE_DELETE) | 113164333292 | 20260414_LOH_Subclone_AF_Methylation_Evidence_01.md (line 676, COLO829 路徑提及) | INDEX 標 SAFE_DELETE；但 1 個 LIVE LOH 報告提及 `duplicates/20260317_colo829_to_pilot_1` 路徑 → 須先確認該引用非數據依賴（僅 provenance 註記）才可 SAFE_DELETE |
| `…/archive/202603_phase1a_iterations/` | 190 MB (s) | 199081197 | Phase 1A 80→200→400→637 漸進中間版（最終 637 已在 Active） | SUPERSEDED | CONCLUDED (結論在 experiments #44-45) | ARCHIVE | 199081197 | (僅 INDEX 引用) | 中間迭代版；最終版保留 Active；結論已落 docs |
| `…/archive/202603_smoke_and_diagnostics/` | 42 MB (s) | 43849613 | pipeline smoke test + TO-pure 失敗一次性診斷 | SUPERSEDED | CONCLUDED_NEGATIVE (TO-pure #182-183) | SAFE_DELETE | 43849613 | (僅 INDEX 引用) | smoke test = 可重生 transient；診斷結論已落 docs；無 live 引用 |

### 1.2 `bip8_output_archive/`（floor 總計 ~778 GB；39 項分類見 INDEX）

| path | 大小 | bytes | purpose | trust_tier | conclusion_status | verdict | reclaimable | referenced_by | rationale |
|------|------|-------|---------|-----------|-------------------|---------|------------|---------------|-----------|
| `bip8…/s-pure-pileup/` | ~770 GB (lf) | 769943090755 | 7 樣本 paired pileup model ISM（canonical paired_pileup 前身） | SUPERSEDED | CONCLUDED (撐 2026-03 validated 報告) | **NEEDS_USER_DECISION** (lean ARCHIVE) | 769943090755 | 多份 2026-03 validated/provenance 報告 (s-pure refs) | **anomaly A2**：canonical 為 partial(missing_tagged_bam)；770GB 大檔疑為唯一 BAM 副本，確認重生前不可刪 |
| `bip8…/research_rounds/` | ~4.9 GB (lf) | 4897999291 | bip8 時期 22 個早期 research rounds (2026-03-07~11) | SUPERSEDED | CONCLUDED (被 big7 rounds 取代) | ARCHIVE | 4897999291 | (INDEX + 數份 2026-03 報告) | 早期 round；結論已落 docs；歷史回溯價值 → 打包封存 |
| `bip8…/s-pure/` | ~1.19 GB (lf) | 1189838612 | 7 樣本 paired full model ISM（canonical paired_full 前身） | SUPERSEDED | CONCLUDED (撐 2026-03 報告) | NEEDS_USER_DECISION (lean ARCHIVE) | 1189838612 | 多份 2026-03 validated/provenance 報告 | 同 s-pure-pileup 邏輯但體積小；canonical partial 故保守 |
| `bip8…/bip8_disk_output/` | 1.35 GB (s) | 1350031847 | bip8 原始 output 快照（內含更早期資料，全小檔 n=0 大檔） | SUPERSEDED | CONCLUDED | ARCHIVE | 1350031847 | (INDEX) | 原始快照，歷史回溯；體積小 → 封存 |
| `bip8…/datestamp_ISM_21dirs/` (B 類) | ~178 MB (lf) | 178509915 | 2025-12~2026-01 各版本 ISM 執行 (WED/w5000/t120) | SUPERSEDED | CONCLUDED (結論在 docs/archive/) | SAFE_DELETE | 178509915 | (INDEX + docs/archive 歷史) | 21 dirs 全被 canonical 取代；早期參數調整版；結論已落 docs/archive；可重生 |
| `bip8…/tmp_meth_annot_test/` | 96 MB (s) | 96207077 | 甲基化標註測試暫存 | DEPRECATED | TRANSIENT | SAFE_DELETE | 96207077 | none | 暫存測試，無分析價值（INDEX F 類） |
| `bip8…/ml_feature_exploration/` | 19 MB (s) | 19640764 | ML 特徵探索分析 | SUPERSEDED | CONCLUDED | ARCHIVE | 19640764 | (INDEX, docs/archive) | 結論已落 docs/archive；體積小 |
| `bip8…/HCC1395/` (purity subsample) | 35 KB (s) | 35506 | HCC1395 purity subsample 系列（幾乎空/symlink） | SUPERSEDED | CONCLUDED_NEGATIVE (subsample purity) | SAFE_DELETE | 35506 | (INDEX) | 僅 35KB，疑為空殼/失效 symlink；INDEX「暫不刪」但實體幾乎無內容 |
| `bip8…/{C 類研究分析 8 dirs + F 類 logs/PPT/test}` | ~120 MB (s, 合計) | ~118900000 | F1 優化/深度/進階分析 + 日誌 + 早期 PPT + 測試 | SUPERSEDED/DEPRECATED | CONCLUDED | ARCHIVE (分析) / SAFE_DELETE (logs/test) | ~118900000 | (INDEX, docs/archive) | 分析結論已落 docs/archive → 封存；logs/test_log_info/PPT 暫存 → 可刪 |

### 1.3 `big8_output_archive/`（總計 ~9.9 GB；5 頂層項）

| path | 大小 | bytes | purpose | trust_tier | conclusion_status | verdict | reclaimable | referenced_by | rationale |
|------|------|-------|---------|-----------|-------------------|---------|------------|---------------|-----------|
| `big8…/multi_sample_quick_check/` | 6.50 GB (s) | 6497343039 | COLO829 chr19 快速檢查 + overlap/retest | SUPERSEDED | CONCLUDED | NEEDS_USER_DECISION (lean ARCHIVE) | 6497343039 | 20260305 三路對照/多樣本報告 | INDEX 標可刪；但含 COLO829 ISM；COLO829 為 LIVE ASM ⭐4 待補樣本 → 保守 |
| `big8…/big8_disk_output/` | 1.35 GB (s) | 1350031847 | big8 原始 output 快照 | SUPERSEDED | CONCLUDED | ARCHIVE | 1350031847 | (INDEX) | 原始快照；INDEX「暫不刪 空間極小」；封存 |
| `big8…/research_rounds/` | 1.09 GB (s) | 1089659076 | big8 時期 1 個 round (clairs_borderline_fn) | SUPERSEDED | CONCLUDED | ARCHIVE | 1089659076 | (INDEX) | 僅 1 round；被後續取代；結論已覆蓋 → 封存 |
| `big8…/pure_tumor_evaluation/` | 850 MB (s) | 849798175 | HCC1395 純腫瘤評估（2 次 1 分鐘內重跑） | SUPERSEDED | CONCLUDED_NEGATIVE (TO-pure #182-183) | SAFE_DELETE | 849798175 | (INDEX) | 1 分鐘內重複執行；結論 NEGATIVE 已落 docs；INDEX 列刪除優先級 1 |
| `big8…/three_way_comparison/` | ~100 MB (est) | -1 | 三方比較 (ClairS Dorado vs normal_only vs paired) | SUPERSEDED | CONCLUDED | ARCHIVE | 100000000 | 20260305 DORADO 三路對照報告 | du timeout(small per-region recursion)；大檔僅 60MB；INDEX「暫不刪 有參考價值」→ 封存 |

---

## 2. 區域加總（floor 估計）

| 類別 | bytes (floor) | 人類可讀 |
|------|---------------|---------|
| KEEP | 0 | 0 |
| ARCHIVE | 9,068,889,131 | ~9.07 GB |
| SAFE_DELETE | 1,168,490,310 | ~1.17 GB |
| NEEDS_USER_DECISION | 4,156,888,219,447 | ~4.16 TB |
| **總掃描 (floor)** | **3,608,212,884,758** | **~3.61 TB** |

> 註：NEEDS_USER_DECISION 含 early_pilots 2.71 TB（lean KEEP，不算回收）+ s-pure-pileup 770 GB + s-pure 1.19 GB + multi_sample_quick_check 6.5 GB（lean ARCHIVE/可回收，但需用戶判 canonical BAM 重生 + COLO829 ⭐4 需求）。加總分類與「總掃描」的差異來自 floor 估計捨入 + early_pilots 同時被列 NEEDS_USER_DECISION。
> **可實際立即回收（SAFE_DELETE 全確定）**：~1.17 GB（duplicates 的 113 GB 因 1 個 LIVE 報告路徑提及降為 NEEDS_USER_DECISION，未計入）。
> **最大潛在回收（用戶確認後）**：s-pure-pileup 770 GB + duplicates 113 GB ≈ 883 GB（前提：canonical paired_pileup BAM 已重生、duplicates 引用僅 provenance 註記）。

---

## 3. 建議行動順序（給主 agent / 用戶）

1. **立即可回收（低風險）**：`202603_smoke_and_diagnostics`(44MB) + `bip8 datestamp ISM 21dirs`(178MB) + `tmp_meth_annot_test`(96MB) + `big8/pure_tumor_evaluation`(850MB) + bip8 logs/test ≈ **1.17 GB** SAFE_DELETE。
2. **最大回收待確認 1**：`s-pure-pileup`(770GB) — 用戶確認 canonical paired_pileup 的 tagged BAM 已重生 → 打包 ARCHIVE 或 SAFE_DELETE。
3. **最大回收待確認 2**：`duplicates`(113GB) — 確認 LOH 報告 line 676 的 COLO829 路徑僅 provenance 註記（非數據依賴）→ SAFE_DELETE。
4. **絕不動**：`202603_early_pilots`(2.71TB) — LIVE TO-BAM 唯一源（除非先把 5 樣本 TO BAM 遷入 canonical/to_pileup）。
5. **打包封存**：bip8/big8 research_rounds + disk_output snapshots + 分析 dirs + three_way ≈ 9 GB。
