<!--
title: Cleanup Proposal — 分層清理提案（SAFE_DELETE / ARCHIVE / NEEDS_USER_DECISION / KEEP）
date: 2026-06-01
type: cleanup_audit_proposal
generated_by: cleanup-audit 綜合 agent（整合 6 區唯讀稽核）
companion: MASTER_CLEANUP_INDEX.md
status: PROPOSAL ONLY — 無任何刪除被執行；所有刪/移為 Hard Gate，須用戶逐項確認
byte_policy: 可回收小計只加總實測 bytes；timeout(large) 項標「需實測」不納入數字
-->

# Cleanup Proposal — 分層清理提案

> **🔴 本檔僅提案，未執行任何刪除/搬移。** 刪除/搬移檔案是 Hard Gate（CLAUDE.md §1，永遠 🔴），須用戶**逐項確認**後由主 agent 執行。
> **誠信聲明（§13）**：所有可回收 bytes 均為 6 區實測 `du -sb`/`stat`/`find` 加總；per-region ISM 樹遞迴 du 全 timeout → 標「需實測」，**未編造**。

> **🛡️ 對抗式驗證修正（2026-06-01，verdict=NEEDS_WORK → 已套用 3 項）**：fresh-context 驗證 agent 逐項 grep 4 類 live 來源（validated reports + evidence_ledger + 論文/HKU 草稿 + CURRENT_FOCUS）後，攔到 1 高風險誤刪 + 2 ARCHIVE caveat，本檔已修正：
> 1. **[高風險，已改]** `big8 pure_tumor_evaluation/` 原列 §1 SAFE_DELETE priority-1 → **改判**：其子目錄 `_143418/`(794M) 是 validated 報告 `20260305_全樣本放寬標準拆解與TPFPF1比較報告_01.md` §2.1–2.3 引用的資料源（per_variant_decision.tsv.gz + 2 README）→ 移入 §2 ARCHIVE 保留；§1 只留棄用首跑 `_143318/`(17M)。「結論 NEGATIVE」≠「報告引用的證據可丟」。
> 2. **[caveat]** §2 `v5 vote_dump_*_genome.tsv.gz` 是 validated `20260513_V6_Attribution_Errata_01.md` 已發佈表（42.76% hp=11 / 17,404 victim，L216-223）的資料源 → ARCHIVE 可，但封存前須確認可由 3 binaries 重生，**絕不可降 SAFE_DELETE**。
> 3. **[caveat]** §2 `ng_kde data/` 含 evidence_ledger L24-25（human_decision=keep）的 2 個 summary JSON → ARCHIVE 須逐檔保留，不可丟。
> **已確認安全**：methyl LOSO（G6 §3 證據）、G6 三軸（phasing/ASM/LOSO-NEGATIVE）證據全未列入任何刪除/封存清單。

---

## 0. 四層可回收空間總表

| 層級 | 定義 | 實測可回收 bytes | 人類可讀 | 另有 timeout(需實測) |
|------|------|----------------:|---------|---------------------|
| **SAFE_DELETE** | 確定可重生 / 已棄用 / 純 transient | 14,977,900,582 | **13.95 GiB** | 無 |
| **ARCHIVE** | 撐已發佈結論的可重生中間檔（打包→封存→確認後刪） | 12,279,921,537 | **11.44 GiB** | kde_rerun_B、O09/fine_pairwise(per-region) |
| **NEEDS_USER_DECISION** | 結論狀態不明 / 可重生但成本高 / 需確認非依賴 | 4,813,703,793,893 | **4.38 TiB** | paired phaseC、s-pure-pileup 部分 per-region |
| **KEEP** | 結論證據 / LIVE infra / baseline 本體 / 狀態機 | 3,709,334,840,333（+ docs/research 小報告 KEEP bulk） | **3.37 TiB+** | canonical/obs per-region |

### 可回收三檔位（給用戶選 risk 等級）

| 情境 | 可回收 | 說明 |
|------|-------:|------|
| **A. 立即（只清 SAFE_DELETE）** | **13.95 GiB** | 零風險；13.8GiB 是 /tmp depth dump，清完即緩解 /tmp 壓力 |
| **B. + ARCHIVE 打包後** | **25.39 GiB** | 打包封存 11.44GiB 中間檔後刪原檔 |
| **C. 全部用戶確認後（含 NEEDS_DECISION）** | **理論上限 4.40 TiB** | 但其中 **2.46 TiB 是 early_pilots（lean KEEP，TO-BAM 唯一源，實際不該回收）** |
| **C′. 扣 early_pilots 的現實上限** | **~1.92 TiB** | s-pure-pileup 770GB + 3 TO BAM 651.9GB + duplicates 113GB + v5 dumps 2.07GB + multi_sample 6.5GB + … |

> ⚠ **B 與 C 的最大現實回收**：用戶確認 canonical paired_pileup 的 tagged BAM 已重生 → `s-pure-pileup`(770GB) 可降 ARCHIVE/SAFE_DELETE；確認 longphase-TO 可重生且不回查 → 3 個 synthesis TO BAM(651.9GB) 可清。這兩塊（~1.42 TiB）是真正「大空間但需驗證」的核心決策。

---

## 1. SAFE_DELETE（13.95 GiB；零/極低風險）

每項：path | 大小 | reason | 引用衝擊。

| path | 大小 | reason | 引用衝擊 |
|------|------|--------|---------|
| `/tmp/claude-1067.bak/.../tasks/b29591ml3.output` | 13.7 GiB | 4 月跑飛的 per-base depth dump（chr1-4 `pos\tdepth`），在過時 .bak 備份內 | none（緩解 /tmp 災情，呼應 feedback_tmp_disk_full） |
| `/tmp/claude-1067.bak/` 其餘 | 數百 MB | active /tmp/claude-1067 的 Apr 29 過時備份副本 | none |
| `bip8 datestamp_ISM_21dirs/` | 178 MB | 2025-12~2026-01 各版本 ISM，全被 canonical 取代 | docs/archive 歷史（結論已落 docs） |
| `bip8 tmp_meth_annot_test/` | 96 MB | 甲基化標註測試暫存 | none |
| `big8 pure_tumor_evaluation/pure_tumor_eval_20260305_143318/` | 17 MB | 棄用首跑（同任務正式版是 `_143418`）；僅此棄用子目錄可刪 | none（_143418 已移 §2 ARCHIVE 保留，見頂部對抗驗證修正 #1）|
| `archive/202603_smoke_and_diagnostics/` | 42 MB | smoke test + TO-pure 失敗一次性診斷 | none（結論在 #182-183） |
| `bip8 logs/test_log_info/PPT`（F 類） | ~部分 of 120MB | 日誌 + 早期 PPT + 測試暫存 | none |
| `bip8 HCC1395 purity subsample` | 35 KB | 幾乎空殼/失效 symlink | none |
| `docs/trash/to_pipeline_staging_v1/` | 8.27 MB | README 自述「已棄用待刪除」；錯 VCF（F1=0.649）；v2 已確認存在 | self README only；v2 = research/to_pipeline_staging（EXISTS✓） |
| `docs/**/__pycache__`（2 dirs） | 120 KB | report-build 腳本 bytecode | none |
| `tsg_promoter_asm_reviewer/msa.log` | 470 KB | MSA build 過程 log | none |
| `synthesis/research_rounds/20260423_colo829_to_kde_full/`（stub log） | 13.7 KB | 孤兒 pipeline.log；實際輸出在 sibling colo829_to_pilot | none |
| `state/compact_snapshots/`（保留最新 1-2） | ~80 KB | precompact hook 自動快照，可重生 | none |

**SAFE_DELETE 小計（原）= 14,977,900,582 B（13.95 GiB）**；**對抗驗證修正後 ≈ 13.18 GiB**（pure_tumor 850M 中移出 _143418 794M 到 ARCHIVE，僅留 _143318 17M）。
> 來源加總：A6 transient 13,800,521,944 + A1 safe_delete 1,168,490,310 + A5 trash+pycache 8,393,279 + A4 msa.log 480,969 + A2 stub log 14,080。修正：A1 safe_delete 扣 pure_tumor _143418（794M）→ 改入 §2 ARCHIVE。

---

## 2. ARCHIVE（11.44 GiB 實測 + timeout 項；打包→封存→確認後刪）

| path | 大小 | reason | 引用衝擊 |
|------|------|--------|---------|
| `bip8/big8_output_archive`（rounds + disk_output snapshots + 分析 dirs + three_way） | ~9.07 GiB | 早期 round/快照/F1 分析，結論已落 docs/archive；歷史回溯價值 | INDEX + 數份 2026-03 報告（描述性） |
| `big8 pure_tumor_evaluation/pure_tumor_eval_20260305_143418/` ⚠ | 794 MB | **對抗驗證 #1 改判**：validated `20260305_全樣本放寬標準拆解與TPFPF1比較報告_01.md` §2.1–2.3 引用源（per_variant_decision.tsv.gz + global_relaxed_standard/README + standard_decomposition/README）；NEGATIVE 結論但報告 provenance 須留 → ARCHIVE 不刪 | 20260305 報告（load-bearing）|
| `hcc1395_normal_pilot/HCC1395_normal_subset.bam` | 1.84 GiB | normal-read FP pilot subset；CONDITIONAL NO-GO；可由 full normal BAM 重 subset | memory read_level_germline_fp（非 load-bearing 證據） |
| `v5_provenance_followup/.../vote_dump_*_genome.tsv.gz`（3 檔） ⚠ | 2.07 GiB | **對抗驗證 #2 caveat**：是 validated `20260513_V6_Attribution_Errata_01.md` 已發佈表（42.76% hp=11 / 17,404 victim，L216-223）的資料源；ARCHIVE 前須確認可由 3 binaries 重生，**絕不可降 SAFE_DELETE** | T1_2 §10 + V6 Errata（load-bearing）|
| `methyl_augmented_filter_phase2/cycle2/data/` | 77 MB | per-sample master_augmented 中間特徵檔；可重生 | none（結論已分離進小 summary） |
| `feature_layered_observation/data + figures` | 406+33 MB | G1-G10 中間特徵 + 圖；LIVE skill 首次示範 | none（建議用戶確認是否留示範） |
| `fp_provenance/.../*_master.tsv.gz` | 57 MB | 2.7M-row FP provenance master（pre-fix 口徑） | none |
| `feature_layered/seqc2/ng_kde/to_staging/z3 各 data/` ⚠ | 合計 ~110 MB | 可重生中間 TSV | **對抗驗證 #3 caveat**：ng_kde data/ 含 evidence_ledger L24-25（human_decision=keep）的 `B5_colo829_s1_fold_detail.summary.json` + `b7_loh_noise_summary.json` → 封存須逐檔保留；ng_kde/seqc2 整體仍 §3 NEEDS_DECISION（無 verdict）|
| `tsg_promoter_asm_reviewer/data/{hg38.gtf.gz,cpgIslandExt.txt.gz}` | 40.6 MB | 公開 reference annotation，可重新下載；移共用 reference | scripts/17（改指向共用路徑即可） |
| `v6_bam_tpfp_hp_loh_cn/step4/.../file_lists/*.txt`（16 檔） | 83 MB | 純路徑列舉中間檔；on/off 成對重複；可重生 | step4_findings（引用 master 非 file_lists） |
| `synthesis/kde_rerun_B_14combos/` | timeout(large) | KDE fix 全量重跑證據本體；aggregate_summary 是 SoT | validated KDE acceptance（留 aggregate，per-region 可封存）|
| `observation_workspaces` NEGATIVE 結案 O 系列（germline_fp 200M / to_fp_post 113M / O09 / fine_pairwise） | 200M+113M+timeout | NO-GO/NEGATIVE 結論證據本體；per-region 可封存 | memory（撐 concluded 結論） |

**ARCHIVE 小計（實測部分）= 12,279,921,537 B（11.44 GiB）**。`kde_rerun_B` / O09 / fine_pairwise per-region 樹 = **需實測**（du timeout）。
> 來源加總：A1 archive 9,068,889,131 + A6 normal subset 1,840,301,808 + A3 archive 711,233,910 + A4 (v6 filelists 86,891,610 + tsg ref 42,605,078) + A2 archive 可量部分 ~530,000,000。v5 dumps 2.07GB 在 A4 標 NEEDS_DECISION→ARCHIVE，已計入 NEEDS_DECISION 避免雙算；若用戶批 ARCHIVE 則從 §3 移入此層。

---

## 3. NEEDS_USER_DECISION（4.38 TiB 實測 + timeout；須用戶判定才能降層）

> 每項標「決策問句」— 用戶答覆後可降 ARCHIVE/SAFE_DELETE 或升 KEEP。

| path | 大小 | reason / 決策問句 | 引用衝擊 |
|------|------|------------------|---------|
| `archive/202603_early_pilots/` | **2.46 TiB** | **lean KEEP**。是 5 樣本 TO tagged BAM 唯一來源 + KDE smoke 依賴。**問句：5 樣本 TO BAM 是否已遷入 canonical/to_pileup？** 否 → 絕不可動 | 20260422_Archive_TO_Rerun / 20260423_KDE_Smoke（LIVE） |
| `bip8 s-pure-pileup/` | 770 GB | lean ARCHIVE。canonical paired_pileup 全 partial(missing_tagged_bam)。**問句：canonical paired_pileup 的 tagged BAM 已重生？** 是 → ARCHIVE/SAFE_DELETE | 多份 2026-03 validated（描述性） |
| 3× `synthesis/.../tumor_tagged.bam`（HCC1395×2+COLO829） | 651.9 GB | 可由 longphase-TO 重生但成本高。**問句：source BAM 仍在 big8 + longphase-TO 環境可用 + 不再回查歷史 tagged？** 最安全先清 #1 legacy partial(278GB，已被 #2 取代) | LIVE phasing n7 Grade A（COLO829）+ 多份 3 月 TO docs |
| `archive/duplicates/` | 113 GB | lean SAFE_DELETE。同任務 2 分鐘內重試副本。**問句：LOH 報告 line 676 的 COLO829 路徑僅 provenance 註記（非數據依賴）？** 是 → SAFE_DELETE | 20260414_LOH_Subclone（1 處路徑提及） |
| `big8 multi_sample_quick_check/` | 6.50 GB | lean ARCHIVE。含 COLO829 ISM。**問句：COLO829 ASM ⭐4 補樣本是否需此 chr19 檢查？** | 20260305 三路對照 |
| `clairs_to_verdict_pilot/data/step1_verdict_vs_seqc2.tsv` | 284 MB | research 最大單檔；可重生。**問句：結案 pilot 原始 join 表是否需保留？** 否 → SAFE_DELETE | none（結論在 confusion_matrix 1.2KB） |
| `seqc2_cnv_stratification/`（無 README） | 31 MB | **docs 10 引用但無 verdict 檔**。**問句：結案狀態？建議補 README** | docs 10 處（撐 CovM/CN-zone 結論） |
| `ng_kde_rescaling/`（無 README） | 35 MB | **docs 13 引用但 conclusion UNKNOWN**。**問句：結案狀態？** | docs 13 處 |
| `coverage_multiple_validation/` | 1.7 MB | **status: blocked 非 concluded**。**問句：待 master rerun，是否續做？** 勿當可清 | docs 4 處 |
| `paired_priority_bug_audit/phaseC_genome_three_way*` | **timeout(需實測)** | 本區最大回收標的（估 tens of GB）。**問句：phaseC vs _with_significance 哪個 superseded？** | validated .md（引用 summary 非 per-region） |
| `tsg_promoter_asm_reviewer/genome_survey/`(v1) | 3.0 MB | memory 暗示 v1=buggy 前版但 script 17 仍引用 v1。**問句：v1 是否已被 v2 取代？** | scripts/17 |
| `presentations/.../20260429_教授報告/`(跨版本 dup figs) | dedup ~7-11 MB | validated 交付物。**問句：base/v2 哪版 canonical？** 確認後 dedup | INDEX/週報 |
| docs untracked `*.standalone.html`（有 .md，5 檔） | ~4.84 MB | 可由 .md 重生但屬 LIVE 研究線。**問句：可改 regenerate-on-demand？** | 各自 .md |
| docs PI `*.standalone.html`（無 .md，3 檔） | ~3.0 MB | **不可重生（HTML 即本體）→ KEEP-leaning**。**問句：確認無 .md 來源？** | PI 消費 |
| `InterSubMod_big7_runbook/` | timeout | 遷移完成 1.5 月未動。**問句：遷移 finalize？** 是 → 降封存 | OUTPUT_INDEX §7 |
| `kde_smoke_test/` + `kde_rerun_pilot/` | timeout | 被 kde_rerun_B 全量取代的過程檔。**問句：可重生不回查？** | validated acceptance（過渡證據） |
| `hcc1395_normal_pilot_global/` | timeout | normal pilot global 版。**問句：與 hcc1395_normal_pilot 重複/取代關係？** | memory read_level_germline_fp |
| observation `*before_hp_fix*` 對照（620M+113M+104M smoke） | ~837 MB | fix 前後對照基準。**問句：post-fix 是否足夠，可清 before？** | none |
| `/tmp/claude-scholar` | 9.6 MB | plugin 安裝物（非 workflow 輸出）。交工具盤點 | n/a |

**NEEDS_USER_DECISION 小計（實測部分）= 4,813,703,793,893 B（4.38 TiB）**，其中 2.46 TiB 是 early_pilots（lean KEEP，不應計入實際回收）。`paired phaseC` 等 per-region 樹 = **需實測**。
> 扣 early_pilots 後的現實可回收上限 ≈ **1.92 TiB**（s-pure-pileup 770GB + 3 TO BAM 651.9GB + duplicates 113GB + multi_sample 6.5GB + v5 dumps 2.07GB + 其餘 < 1GB 級）。

---

## 4. KEEP（不列回收；含絕不動清單）

| 類別 | 代表 path | 理由 |
|------|----------|------|
| **baseline 本體** | 14× `canonical/.../longphase_s/*_tagged.bam`（3.37 TiB） | paper LIVE baseline；重生需重跑 longphase-s + source BAM |
| **狀態機 / SoT（絕不動）** | `research/{autoresearch,data_registry,_template}/` | 活躍狀態機 + 不可 git rm 清單（CLAUDE.md） |
| **LIVE 母專案** | `tpfp_loh_af_kde_discrimination`、`loh_investigation`、`germline_asm_analysis`、`v6_bam_tpfp_hp_loh_cn`、`paired_priority_bug_audit`、`hku_collaboration`、`v5_provenance_followup`(results) 等 15+ | 論文主軸 / LIVE characterization / handoff 紀錄 |
| **結論證據（小 tracked）** | 各結案專案 `manifest.yaml / 00_INDEX.md / *_report.md / summary.tsv`（KB 級） | 撐 INDEX/MEMORY/evidence_ledger |
| **收尾分析** | `synthesis/concluded/`（202M）、`multilayer_hp_benchmark`、`observation loh_round1`（master dataset 入口） | 撐 self-phasing/LOH 已發佈結論 |
| **docs 交付物** | validated 終版報告、`docs/archive/`、references 手冊、tracked dashboard HTML | 終版/封存區/啟動上下文 |
| **data symlink** | `data/bam,ref,answer/*`（0 byte big7 實佔） | 全 symlink→big8，不佔空間 |

---

## 5. archive-then-delete 指令範本（須用戶逐項確認）

> 🔴 以下為**範本**，非自動執行。每條須用戶確認後由主 agent 在互動模式逐項跑。所有路徑用實際絕對路徑替換。

### 5.1 SAFE_DELETE 範本（先 du 確認大小 → 刪 → 確認消失）

```bash
# Step 1 確認 target 大小與內容（不刪）
TARGET="/tmp/claude-1067.bak"
du -sh "$TARGET"
ls -la "$TARGET"

# Step 2（須用戶 ack）刪除 — Hard Gate，僅在用戶確認後執行
# rm -rf "$TARGET"

# Step 3 驗證已回收
df -h /tmp
```

### 5.2 ARCHIVE 範本（tar → 驗證 → 比對 count → 確認後刪原檔）

```bash
# 變數
SRC="/big7_disk/liaoyoyo2001/research/methyl_augmented_filter_phase2/cycle2/data"
ARCHIVE_DIR="/big7_disk/liaoyoyo2001/big7_disk_output/_ARCHIVE_2026_06_01/research_intermediate"
NAME="methyl_filter_phase2_cycle2_data"

# Step 1 建封存資料夾 + 記錄來源大小與 file count
mkdir -p "$ARCHIVE_DIR"
SRC_BYTES=$(du -sb "$SRC" | cut -f1)
SRC_FILES=$(find "$SRC" -type f | wc -l)
echo "source: $SRC_BYTES bytes, $SRC_FILES files"

# Step 2 打包（gzip）
tar -czf "$ARCHIVE_DIR/${NAME}.tar.gz" -C "$(dirname "$SRC")" "$(basename "$SRC")"

# Step 3 驗證 tar 完整 + file count 一致（不刪原檔前必過）
tar -tzf "$ARCHIVE_DIR/${NAME}.tar.gz" >/dev/null && echo "tar OK"
TAR_FILES=$(tar -tzf "$ARCHIVE_DIR/${NAME}.tar.gz" | grep -c -v '/$')
echo "tar files: $TAR_FILES (expect ~$SRC_FILES)"

# Step 4（須用戶 ack）file count 相符 + tar OK → 才刪原檔（Hard Gate）
# [ "$TAR_FILES" -eq "$SRC_FILES" ] && rm -rf "$SRC"

# Step 5 紀錄到封存 manifest
echo "$NAME | $SRC | $SRC_BYTES bytes | $SRC_FILES files | archived 2026-06-01" >> "$ARCHIVE_DIR/ARCHIVE_MANIFEST.txt"
```

### 5.3 大型 BAM（s-pure-pileup / TO tagged BAM）前置確認範本

```bash
# 770GB / 651.9GB 級不打 tar（無壓縮空間且耗時）；改：先確認可重生 → 直接 SAFE_DELETE 或移外部 disk
# 前置（用戶確認 canonical tagged BAM 已重生才能進行）：
CANON="/big7_disk/liaoyoyo2001/big7_disk_output/canonical"
find "$CANON" -maxdepth 6 -name '*_tagged.bam' -path '*paired_pileup*' -printf '%p  %s bytes\n'
# ↑ 若列出對應 sample 的 paired_pileup tagged BAM（非空）→ 確認 s-pure-pileup 可降 SAFE_DELETE
# 移外部 disk（若保留歷史）：
# rsync -a --info=progress2 "$SRC_BAM" /big8_disk/_COLD_ARCHIVE/ && <驗證 md5> && rm "$SRC_BAM"
```

---

## 6. 建議執行順序（風險遞增）

1. **零風險立即**：清 SAFE_DELETE 13.95 GiB（先 `/tmp/.bak` 13.7GiB 解 /tmp 壓力，再 docs trash + pycache + stub）。
2. **低風險打包**：ARCHIVE 11.44 GiB research/output 中間檔 → `_ARCHIVE_2026_06_01/`（tar+驗證+確認後刪）。
3. **中風險確認**：`duplicates`(113GB) 確認 LOH 引用僅 provenance → SAFE_DELETE。
4. **最大回收待驗證**：`s-pure-pileup`(770GB) 確認 canonical 重生 + 3 TO BAM(651.9GB) 確認 longphase-TO 可重生 → 清/移。
5. **絕不動**：`early_pilots`(2.46TiB TO-BAM 源)、`autoresearch/data_registry/_template`、14 canonical BAM、`coverage_multiple_validation`(blocked)、2 個無 verdict 高引用專案（先補 README）。

---

*稽核者：cleanup-audit 綜合 agent | 2026-06-01 | 唯讀提案，未執行任何刪除/搬移 | 可回收 bytes 來自 6 區實測；timeout 項標「需實測」未編造*
