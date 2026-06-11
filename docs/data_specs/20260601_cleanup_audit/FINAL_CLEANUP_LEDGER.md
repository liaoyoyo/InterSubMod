<!--
title: Final Cleanup Ledger — 最終清理執行清單（搬移/封存/刪除 + 檢視保留）
date: 2026-06-01
type: cleanup_execution_ledger
companion: CLEANUP_PROPOSAL.md / MASTER_CLEANUP_INDEX.md
status: PLAN OF RECORD — 尚未執行任何刪除/搬移；所有刪/移為 Hard Gate，待用戶最終確認
byte_policy: 所有 bytes 來自本 session 實測 ls/du/find（§13 無編造）；未實測者標「待 du」
-->

# 最終清理執行清單（Final Cleanup Ledger）

> **🔴 狀態：尚未執行任何刪除/搬移。** 本檔是「將執行的清單 + 檢視過但保留的清單 + 各自原因」。待用戶最終確認後，主 agent 在互動模式逐項執行，並回填 §6 實際結果。
> **基於**：cleanup-audit workflow（8 agents）+ 對抗驗證（NEEDS_WORK 已修）+ 本 session 對每個 TiB 級候選的逐一實測（ls/du/find/grep）。

---

## §0 TL;DR — 驗證後的真實可回收（誠實版）

| 區塊 | big7 可回收 | 風險 | 一句話 |
|------|-----------:|------|--------|
| **A. 可重生 TO BAM（只刪 .bam，保留衍生 LOH.bed/vcf）** | **~623 GiB** | 低 | BAM 可由 big8 源 + longphase-TO 重生；論文用的 LOH.bed/VCF 全留 |
| **B. 純 transient / 棄用** | ~330 MB（big7）+ 13 GiB（root /tmp）| 零 | 跑飛備份、bytecode、棄用首跑、trash |
| **C. s-pure-pileup（待你定奪）** | 待 du（疑 ~770 GB）| 中 | ISM 衍生 CSV，被多份 validated/PI 報告引用 → 預設 KEEP，除非你接受「可重跑 ISM 重生」trade-off |
| **KEEP（檢視過不動）** | — | — | canonical 14 BAM 3.37 TiB + early_pilots 2.46 TiB + 現役 COLO829 + 全 LOH.bed + 狀態機 |

> **與你問卷選擇的差異（驗證後修正，見 §1）**：你勾的「duplicates 113G 可清」其中含 **COLO829 LOH.bed（LIVE n7 論文引用）** → 改為「只刪 BAM 留 LOH.bed」；「canonical 14 BAM」你選 KEEP（已遵照）；s-pure 你選可回收，但驗證發現它被多份 validated 報告當 provenance 引用 → 降為 §1 待定奪。

---

## §1 對問卷選擇的驗證修正（必讀）

| 你的選擇 | 驗證發現 | 修正後做法 |
|---------|---------|-----------|
| **duplicates 113G 可清** | `duplicates/20260317_colo829_to_pilot_1/` 的 `tumor_phased_LOH.bed` 是 **LIVE 論文 `20260414_LOH_Subclone_AF_Methylation_Evidence_01.md` §P2.2 引用的 COLO829 LOH 真值**（7 樣本 LOH.bed 來源之一）。且該 COLO829 與 `20260423_colo829_to_pilot` 幾乎位元相同（BAM 差 1 byte）。 | **不整刪**。只刪 `tumor_tagged.bam`(95.4G) + `tumor_tagged.partial_*.bam`(16.8G) + `.bai`；**保留 LOH.bed / GE.bed / vcf / logs**（provenance）。 |
| **s-pure 770G + 3 TO BAM 652G 可回收** | 3 TO BAM 確認可重生（big8 源在），但 **s-pure-pileup 是 ISM 衍生 CSV（非 BAM）**，被 `20260305_全樣本放寬標準...`、`20260422_truth_sets`、`20260422_Self_Phasing_PI` 等 validated/PI 報告引用。刪它＝刪報告引用的衍生數據（只能重跑 ISM 重生，成本高）。 | TO BAM → §2A 執行；**s-pure-pileup → §5 待你定奪**（預設 KEEP）。 |
| **canonical 14 BAM (3.37 TiB) KEEP** | 已實測確認 14 檔實存（HCC1937 兩檔各 ~450GB 最大；manifest 標 tagged_bam_ready=false 是過時錯誤）。 | 遵照你選擇 **KEEP**（§4）。 |
| **特別保留 `bip8_output_archive/20260119_all-with-w5000_1/`** | 實測：含 `analysis/{plots,tables}`（plots 有 correlation/distribution/evaluation/heatmaps 4 類子目錄）+ filtered_snv_tp/fp。⚠ `.png` 實測 0 張（depth≤4）→ 圖可能是其他格式或在更深層。 | **KEEP**（§4，遵照你指令）；附註：請確認 plots 子目錄是否含你預期的圖。 |

---

## §2 將執行：刪除清單（待最終確認）

### §2A 可重生 TO tagged BAM — 只刪 `.bam`+`.bai`，保留衍生輸出（big7 回收 ~623 GiB）

> **原理**：每個 TO 目錄 = 大 BAM（可重生）+ 小衍生（LOH.bed/GE.bed/vcf/logs，撐論文與 provenance）。只刪 BAM，衍生全留。
> **重生路徑**：`longphase --somaticMode`（即 LongPhase-S/TO）對 big8 源 tumor BAM（`data/bam/HCC1395/tumor.bam` → big8）重跑 → 復原 `tumor_tagged.bam`。

| # | 要刪的檔 | 大小 | 時間 | 保留的衍生（provenance） | 原因 |
|---|---------|-----:|------|------------------------|------|
| 1 | `synthesis/research_rounds/legacy_partials/20260307_hcc1395_to_pilot_1_tagged_bam_only/step03_longphase_to/tumor_tagged.bam`(+.bai) | 278.42 GB | Mar 11 | （此目錄**只有** BAM+bai，無衍生）→ 整目錄可刪 | HCC1395 TO **legacy**，已被 `20260315_hcc1395_to_pilot`（有完整衍生）取代；名稱即 "tagged_bam_only" |
| 2 | `synthesis/research_rounds/20260315_hcc1395_to_pilot/step03_longphase_to/tumor_tagged.bam`(+.bai) | 278.42 GB | Mar 15 | **保留** LOH.bed(74KB)/GE.bed/vcf(655MB)/logs | HCC1395 TO 現役衍生（LIVE 論文用 LOH.bed）；BAM 本身可重生，論文不需常駐 BAM |
| 3 | `synthesis/research_rounds/archive/duplicates/20260317_colo829_to_pilot_1/step03_longphase_to/tumor_tagged.bam` + `tumor_tagged.partial_20260317_212953.bam`(+.bai) | 112.23 GB | Mar 17 | **保留** LOH.bed(12.8KB,論文引用)/GE.bed/vcf/logs | COLO829 Mar17 與 Apr23 幾乎位元相同（重複）；partial 為中斷殘檔；LOH.bed 留給論文 |

**§2A 小計（big7 回收）= 669,074,617,720 B ≈ 623.1 GiB**（278.42+278.42+112.23 GB）。
> ⚠ COLO829 **Apr23**（現役，95.4G）**不刪**（ASM ⭐4 待補樣本可能需重讀；見 §4）。HCC1395 legacy(#1) symlink `output/20260307_hcc1395_to_pilot_1` 與數份 3 月報告引用此路徑 → 刪後我會在原目錄留 `_DELETED_README.txt` tombstone（記重生法），避免斷指標。

### §2B 純 transient / 棄用（零風險）

| 要刪 | 大小 | 位置 | 原因 |
|------|-----:|------|------|
| `/tmp/claude-1067.bak/` | 13 GiB | root `/`（非 big7）| Apr 29 過時備份，內含跑飛 per-base depth dump；緩解 root /tmp |
| `big8 .../pure_tumor_evaluation/pure_tumor_eval_20260305_143318/` | 17 MB | big8 | 棄用首跑（正式版 `_143418` 保留，見 §4 對抗驗證 #1）|
| `bip8_output_archive/datestamp_ISM 早期 21 dirs` | ~178 MB | big7 | 2025-12~2026-01 早期 ISM，全被 canonical 取代 |
| `bip8_output_archive/tmp_meth_annot_test/` | ~96 MB | big7 | 甲基標註測試暫存 |
| `synthesis/research_rounds/archive/202603_smoke_and_diagnostics/` | ~42 MB | big7 | smoke test + TO-pure 一次性診斷（結論已落 docs）|
| `docs/trash/to_pipeline_staging_v1/` | 8.27 MB | big7 | README 自述棄用待刪；v2 存於 research/ |
| `docs/**/__pycache__`、`tsg.../msa.log`、`colo829_to_kde_full/` stub log、`state/compact_snapshots/`（留最新 1-2）| <1 MB | big7 | bytecode / build log / 孤兒 stub / 可重生快照 |

**§2B 小計**：big7 ~324 MB + root /tmp 13 GiB + big8 17 MB。

---

## §3 將執行：搬移/封存（tar → `big7_disk_output/_ARCHIVE_2026_06_01/`）

> 你選「tar 留 big7 _ARCHIVE_2026_06_01/」。流程：tar → `tar -tzf` 驗證 → file count 比對 → **確認後**刪原檔（§CLEANUP_PROPOSAL §5.2 範本）。
> ⚠ 封存只回收「原 dir 變 tar」的差（壓縮比），不跨 disk；對可重生中間檔，封存主要是「整理＋留可回溯」而非省大空間。

| 封存項 | 大小 | 原因 | caveat |
|--------|-----:|------|--------|
| `bip8_output_archive` + `big8_output_archive`（除 §4 保留的 20260119 與 pure_tumor _143418）| ~9 GiB | 早期 round/F1 分析快照；結論已落 docs/archive；SUPERSEDED 但留歷史回溯 | 排除 20260119（§4 KEEP）|
| `research/v5_provenance_followup/.../vote_dump_{baseline,v3f,v5}_genome.tsv.gz` | 2.07 GiB | 全基因組 vote dump，可由 3 binaries 重生 | ⚠ 是 validated `20260513_V6_Attribution_Errata` 表的源 → 封存可、**絕不可刪**（對抗驗證 #2）|
| `research/.../ng_kde_rescaling/data/`（含 2 個 ledger-cited JSON）| ~35 MB | 可重生中間 TSV | ⚠ 須逐檔保留 `B5_colo829_s1_fold_detail.summary.json`+`b7_loh_noise_summary.json`（ledger human_decision=keep，對抗驗證 #3）|
| `research/methyl_augmented_filter_phase2/cycle2/data/` | 77 MB | per-sample master_augmented 中間檔；可重生 | ⚠ **不可動** cycle4 LOSO + phase2_completeness_audit（G6 論文 §3 證據）|
| `hcc1395_normal_pilot/HCC1395_normal_subset.bam` | 1.84 GiB | normal-read FP pilot subset（CONDITIONAL NO-GO）；可由 full normal BAM 重 subset | — |

**§3 小計（封存後可選刪原檔）≈ 13 GiB**（多為留可回溯，非主要空間回收）。

---

## §4 檢視過但不動（KEEP）— 含標註原因

| 保留項 | 大小 | 為何不動（原因） |
|--------|-----:|----------------|
| `canonical/*/paired_*/*/longphase_s/*_tagged.bam`（14 檔）| **3.37 TiB** | **你選 KEEP**；論文 baseline；14 檔實測確存（HCC1937 ~450GB×2 最大）|
| `synthesis/research_rounds/archive/202603_early_pilots/`（5 樣本 TO）| **2.46 TiB（待 du 校）** | **LIVE n7 論文 §P2.2 的 DORADO/H1437/H2009/HCC1937/HCC1954 LOH.bed 唯一來源** + KDE smoke 依賴；INDEX 標 LOW_RISK 是過時危險分類 |
| `20260423_colo829_to_pilot/.../tumor_tagged.bam` | 95.4 GB | COLO829 **現役** TO tagged BAM；ASM ⭐4 待補樣本可能需重讀；非重複 |
| 全部 `tumor_phased_LOH.bed` / GE.bed / vcf（各 TO 目錄）| KB–MB | LIVE 論文真值與 provenance；§2A 刪 BAM 時全保留 |
| `bip8_output_archive/20260119_all-with-w5000_1/` | 待 du | **你指定保留**（較完整基準 + analysis/plots/tables）|
| `bip8 s-pure-pileup/` | 待 du（疑 ~770G）| 被多份 validated/PI 報告引用 → §5 待定奪（預設 KEEP）|
| `big8 .../pure_tumor_eval_20260305_143418/` | 794 MB | validated `20260305_全樣本放寬標準...` §2.1-2.3 引用源（對抗驗證 #1 救回）|
| `research/{autoresearch,data_registry,_template}/` | — | 活躍狀態機 / SoT，**絕不動**（CLAUDE.md）|
| LIVE 母專案（tpfp_loh_af_kde / germline_asm / v6_bam / paired_priority_bug_audit / hku_collaboration / loh_* 等）| — | 論文主軸 / characterization / handoff 紀錄 |
| `research/methyl_augmented_filter_phase2/cycle4/loso_validation/` + `phase2_completeness_audit/` | — | G6 論文 §3 LOSO-NEGATIVE 證據 + PI signoff 引用（甲基 filter 雖 DEAD，此證據是論文 methods）|
| 各結案專案小 tracked 結論檔（manifest/00_INDEX/summary/report）| KB 級 | 撐 INDEX/MEMORY/evidence_ledger |
| validated 終版 docs / `docs/archive/` / references 手冊 | — | 終版 / 封存區 / 啟動上下文 |
| `data/bam,ref,answer/*` symlink | 0 byte | 全 symlink→big8，big7 不佔空間 |

---

## §5 待你最終定奪

1. **s-pure-pileup（疑 ~770 GB，big7 最大單一候選）**：是 PRE-FIX 的 ClairS 7-樣本 paired-pileup ISM 衍生 CSV，已被 canonical/paired_pileup（含 tagged BAM，實測確存）取代；但被 `20260305`/`20260422_truth_sets`/`20260422_Self_Phasing_PI` 等 validated/PI 報告當 provenance 引用。
   - **選項 A（預設）KEEP**：保 provenance，big7 不回收這 770G。
   - **選項 B 刪**：接受「報告引用的是已算好的結論，原始 ISM CSV 可重跑重生」，刪前我先在那幾份報告加 `regenerate-on-demand` 註記。回收 ~770G。
2. **§2A 三個 TO BAM 刪除**（623 GiB）— 確認執行？（LOH.bed/vcf 全留，BAM 可重生）
3. **執行模式**：要我（a）只做 §2B 零風險（~330MB + /tmp 13G），還是（b）含 §2A 的 623 GiB BAM，還是（c）全部含 §5 s-pure？

---

## §6 執行後回填（2026-06-01，用戶確認「§2B+§2A+§3」後執行）

> 用戶選定「回收 ~623 GiB BAM + 封存」。執行模式：互動模式逐項，守衛式（先驗證衍生輸出存在才刪 BAM），每步 df 對比。

### ✅ §2A 已執行 — 3 個可重生 TO BAM（衍生輸出全保留）

| # | 刪除的檔 | 釋出 | 保留的衍生（provenance）| 狀態 |
|---|---------|-----:|----------------------|------|
| 1 | `legacy_partials/20260307_hcc1395_to_pilot_1_tagged_bam_only/.../tumor_tagged.bam`+.bai | 278.42 GB | （此目錄原本只有 BAM）→ 留 `_DELETED_README.txt` tombstone（記重生法）| ✓ DELETED |
| 2 | `20260315_hcc1395_to_pilot/step03_longphase_to/tumor_tagged.bam`+.bai | 278.42 GB | ✓ 保留 LOH.bed / GE.bed / vcf(655MB) / logs | ✓ DELETED |
| 3 | `archive/duplicates/20260317_colo829_to_pilot_1/.../tumor_tagged.bam` + partial + .bai | 112.23 GB | ✓ 保留 LOH.bed(論文引用) / GE.bed / vcf / logs | ✓ DELETED |

**§2A 實際釋出 = 669.07 GB ≈ 623.1 GiB。df big7 實測：39T used / 1.4T free / 97% → 38T used / 2.0T free / 96%（≈ +600 GiB free）。**
> COLO829 **Apr23**（現役，95.4G）未動（KEEP，ASM ⭐4）。所有 LOH.bed/vcf 衍生完整保留 → LIVE 論文 provenance 無損。

### ✅ §2B 已執行 — transient / 棄用

| 刪除 | 釋出 | 狀態 |
|------|-----:|------|
| `/tmp/claude-1067.bak/`（跑飛 depth dump 備份）| 13 GiB（root `/`）| ✓ DELETED（root 701G→714G free）|
| `docs/trash/to_pipeline_staging_v1/`（v2 在 research/）| 8.27 MB | ✓ DELETED |
| `synthesis/.../20260423_colo829_to_kde_full/`（孤兒 pipeline.log stub）| 13.7 KB | ✓ DELETED |
| `bip8_output_archive/tmp_meth_annot_test/`（甲基測試暫存）| 92 MB | ✓ DELETED |
| `research/tsg_promoter_asm_reviewer/msa.log` | 470 KB | ✓ DELETED |
| `state/compact_snapshots/` 舊 7 檔（留最新 2）| ~70 KB | ✓ DELETED |

### ⚠ 未能執行 — big8 唯讀
- `big8 .../pure_tumor_eval_20260305_143318`（17M 棄用首跑）：**`/big8_disk` 是唯讀 NFS 掛載**（`140.123.104.246:/big8_disk`）→ **無法從本機刪除**。big8 上所有清理皆不可行。維持不動（且 `_143418` 正式版本就 KEEP）。

### §3 封存
- 建立封存資料夾 `big7_disk_output/_ARCHIVE_2026_06_01/`（+ README manifest）。
- ⚠ **誠實工程現實**：big8 唯讀（無冷儲存目標）+ 同 disk tar 對已壓縮/可重生資料回收 ~0 → §3「tar 留 big7」不回收空間，只是整理/隔離。故大型 SUPERSEDED archive（bip8/big8、s-pure）**改列「需另做刪除決策才能回收」**（見下方後續建議），不做無謂的同 disk tar。
- 隔離進封存區的項目見 `_ARCHIVE_2026_06_01/ARCHIVE_MANIFEST.md`。

### 📊 本次總回收
- **big7：≈ 623 GiB（97% → 96%，1.4T → 2.0T free）** — 全部來自 §2A 可重生 TO BAM（衍生/provenance 完整保留）。
- **root /tmp：13 GiB**（過時備份）。
- big7 上 §2B 小檔 ~100 MB。

### 後續可選的更大回收槓桿（需你再決策）
1. **s-pure-pileup（疑 ~770G）**：§5 trade-off（接受 ISM 可重跑 → 刪，先在報告加 regenerate 註記）。
2. **SUPERSEDED bip8/big8 早期 archive**：同 disk「封存」不回收 → 要回收需**刪除**（非封存）。它們被 3 月歷史報告描述性引用 → 刪需加 tombstone。這是一個獨立 Hard Gate 決策。
3. **canonical 14 tagged BAM（3.37 TiB）→ 經 2026-06-01 深查證實應 KEEP（非槓桿）**：原以為是「regenerate-on-demand 最大槓桿」，但 grep 反證 —— HCC1395 canonical tagged BAM 被 **LIVE ASM/methyl-phasing 分析（2026-05-31 的 germline_het_null.py / shuffle_control.py / seqc2_cn_methyl.py + 20260527 TSG ASM reviewer response）直接讀取**（HP:Z + MM/ML tags）；其他 6 樣本的 tagged BAM 正是 **cross-sample ASM ⭐4 升級的 gating 資源**（專案自述他樣本 tagged BAM 稀缺）。→ **這 3.37 TiB 是 ASM 軸的 live working set，非可回收 junk。刪除會直接傷害本週進行中的 ASM 研究。維持 KEEP。**

---

*稽核：cleanup-audit workflow + 主 agent 逐項實測 | 2026-06-01 | 唯讀提案，未執行任何刪除/搬移 | bytes 全來自本 session 實測*
