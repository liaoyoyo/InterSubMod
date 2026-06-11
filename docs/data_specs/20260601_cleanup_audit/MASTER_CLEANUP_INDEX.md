<!--
title: Master Cleanup Index — 全專案大型輸出/資料夾索引與說明
date: 2026-06-01
type: cleanup_audit_master_index
generated_by: cleanup-audit 綜合 agent（整合 6 區唯讀稽核；READ-ONLY，未刪除/搬移任何檔案）
scope:
  - big7_disk_output（= InterSubMod/output symlink target；canonical + synthesis + 頂層 standalone + bip8/big8 archive）
  - research/（已結案 + LIVE 專案）
  - docs/
  - data/ + transient（/tmp、state/compact_snapshots）
sources: area_1_output_archives.md, area_2_output_live.md, area_3_research_concluded.md,
         area_4_research_live.md, area_5_docs.md, area_6_bam_transient.md
byte_policy: bytes 來自各 area 實測 du -sb / stat / find -maxdepth；per-region ISM 樹遞迴 du 全部 timeout → 標 unknown，不編造
companion: CLEANUP_PROPOSAL.md（分層清理提案 + archive-then-delete 指令範本）
-->

# Master Cleanup Index — 全專案大型輸出/資料夾索引與說明

> **唯讀稽核成果整合。本檔不執行任何刪除/搬移**（刪/移為 Hard Gate，須用戶逐項確認後由主 agent 執行）。
> **誠信聲明（§13）**：所有 byte 數字均來自 6 區 area 報告的實測 `du -sb` / `stat` / `find -maxdepth`。per-region ISM CSV 子目錄（每 run 數十萬個 region 資料夾）的遞迴 `du` **全部 timeout**，該等項目大小標 **unknown / timeout(large)**，**未編造任何數字**。可回收總計僅加總「實測 bytes」；timeout 項另列「需實測」。

---

## 0. TL;DR — 四句話

1. **磁碟絕大多數體積在兩處**：`big7_disk_output/synthesis/research_rounds/archive/202603_early_pilots/`（**2.46 TiB**，但它是 5 樣本 TO BAM 唯一來源 → **lean KEEP**）與 `canonical/*/paired_*/*_complete_matrix/longphase_s/*_tagged.bam`（**14 檔 3.37 TiB**，paper baseline → **KEEP**）。這兩塊合計 ~5.8 TiB 幾乎都不該刪。
2. **立即可清（SAFE_DELETE，確定可重生/已棄用）= ~13.95 GiB**，其中 **13.8 GiB 是 `/tmp/claude-1067.bak` 一個跑飛的 depth dump** + docs trash v1（8.27 MB）+ 一堆 stub log / pycache。
3. **真正大的「可回收但需確認」候選**：`s-pure-pileup`（770 GB，須確認 canonical paired_pileup tagged BAM 已重生）、3 個 synthesis TO tagged BAM（651.9 GB，可由 longphase-TO 重生）、`duplicates/`（113 GB，須確認 LOH 報告引用僅 provenance）。全列 **NEEDS_USER_DECISION**。
4. **兩個跨 area 重大校正**：(a) Area 6 漏計 canonical 14 個 3.37 TB tagged BAM，本檔以 Area 2 實測為準；(b) 兩個高 docs 引用 research 專案（`ng_kde_rescaling` 13 引用 / `seqc2_cnv_stratification` 10 引用）**無 README/verdict 檔**，結論狀態無法自證 → 勿輕動。

> 註：本檔 byte 單位以 **二進位（GiB/TiB，÷1024）** 呈現人類可讀值（與 area 報告 `du -sb` 一致）；area 報告偶用十進位 GB/TB 描述，數值來源相同。

---

## 1. 四大區域地圖（大小 + 子結構）

```
big7_disk/liaoyoyo2001/
├── big7_disk_output/                       ← InterSubMod/output symlink target（最大宗）
│   ├── canonical/                          ← [Area2] 19 runs baseline；★含 14 個 tagged BAM = 3.37 TiB（KEEP）
│   ├── synthesis/                          ← [Area2] 分析證據庫（concluded/observation_workspaces/research_rounds/kde_rerun_*）
│   │   └── research_rounds/
│   │       ├── (active rounds)             ← [Area2] LIVE LOH/phase1a 分析產物（KEEP）
│   │       ├── *_to_pilot/.../tumor_tagged.bam ← [Area6] 3 個 TO tagged BAM = 651.9 GB（NEEDS_DECISION）
│   │       └── archive/                     ← [Area1] 34 dirs 歷史 round（early_pilots 2.46TiB / duplicates 113GB …）
│   ├── bip8_output_archive/                ← [Area1] 39 dirs；s-pure-pileup 770GB SUPERSEDED（NEEDS_DECISION）
│   ├── big8_output_archive/                ← [Area1] 5 dirs ~9.9GB（多 ARCHIVE）
│   └── (頂層 standalone)                    ← [Area2] multilayer_hp_benchmark / kde_smoke_test / normal_pilot / runbook
├── research/                               ← [Area3 結案 + Area4 LIVE]
│   ├── (LIVE 母專案 15+)                    ← KEEP（autoresearch/data_registry/_template 絕不動）
│   │   └── */data|figures/                  ← 中間特徵 TSV（.gitignored，可重生）→ ARCHIVE
│   └── (結案 26 專案)                        ← KEEP 報告 + ARCHIVE data/
├── InterSubMod/data/                       ← [Area6] 幾乎全 symlink→big8（big7 實佔 ~0）
├── InterSubMod/docs/                       ← [Area5] 106 MiB（presentations 52M / reports 32M / experiments 12M / trash 8M）
└── /tmp + state/compact_snapshots          ← [Area6] transient（.bak 13.7GB depth dump = 最大單一 SAFE_DELETE）
```

### 1.1 區域大小一覽（實測；timeout 標明）

| 區域 | 主要 scope | 實測大小 | timeout 項 | 主結論 |
|------|-----------|---------|-----------|--------|
| **output — archive 區**（Area 1） | bip8/big8/research_rounds archive | floor ~3.61 TiB（only >5MB 大檔加總） | per-region 全 timeout | 多 SUPERSEDED 但含 2 個 anomaly（early_pilots LIVE / s-pure-pileup canonical-partial） |
| **output — 現役區**（Area 2） | canonical + synthesis 非 archive + standalone | ~3.38 TiB（可量部分） | canonical/obs/kde per-region timeout | 14 canonical BAM 3.37 TiB = KEEP 主體；分析產物撐結論 |
| **research — 結案**（Area 3） | 26 結案/DEAD 專案 | ~1.03 GiB | 部分 per-region timeout | 結論在小 tracked 報告（KEEP）；data/ 中間檔可 ARCHIVE |
| **research — LIVE**（Area 4） | 15+ LIVE 母專案 | ~2.64 GiB（可量部分） | paired_priority_bug_audit timeout | 母專案全 KEEP；內部 v5 dumps 2.26GB + phaseC（timeout）可回收 |
| **docs**（Area 5） | docs/ 20 subdirs | 106.4 MiB | 無 | trash v1 8.27MB SAFE_DELETE；其餘多 KEEP/小 NEEDS_DECISION |
| **data + transient**（Area 6） | data/ + BAM 清單 + /tmp | data ~0；transient 13.8GB | archive 深掃 timeout | data 全 symlink；3 TO BAM 651.9GB + transient 13.8GB |

---

## 2. 主要 artifact 一覽表（purpose + trust_tier + conclusion_status + verdict）

> 只列「大型 / 高引用 / 有決策意義」的 artifact。完整逐項見各 area_*.md。
> trust_tier 詞彙：CURRENT / PRE-FIX（HP/KDE/self-phasing 修正前）/ SUPERSEDED / DEPRECATED / TRANSIENT / NA(reference)。

### 2.1 output 區（最大宗）

| artifact | 大小 | purpose（為何產生） | trust_tier | conclusion | verdict | area |
|----------|------|---------------------|-----------|-----------|---------|------|
| `synthesis/research_rounds/archive/202603_early_pilots/` | 2.46 TiB | 建 canonical 前探索 pilot；**實含 5 樣本 TO tagged BAM 唯一來源** | INDEX 宣稱 SUPERSEDED / **實為 CURRENT TO-BAM** | LIVE（TO-rerun + KDE smoke 依賴） | **NEEDS_DECISION（lean KEEP）** | 1 |
| `canonical/.../longphase_s/*_tagged.bam`（14 檔） | 3.37 TiB | longphase-s paired tagged BAM；HP 欄位原始載體 | PRE-FIX(HP) | LIVE | **KEEP**（baseline 本體） | 2 |
| `bip8_output_archive/s-pure-pileup/` | 770 GB | 7 樣本 paired_pileup ISM（canonical 前身） | SUPERSEDED | CONCLUDED（撐 2026-03 validated） | **NEEDS_DECISION（lean ARCHIVE）** | 1 |
| 3× `synthesis/.../tumor_tagged.bam`（HCC1395×2 + COLO829） | 651.9 GB | TO longphase tagged BAM | PRE-FIX/SUPERSEDED | CONCLUDED + LIVE-adjacent | **NEEDS_DECISION**（可 longphase-TO 重生） | 6 |
| `synthesis/research_rounds/archive/duplicates/` | 113 GB | 同任務 2 分鐘內連續重試副本（4×main + 3×to） | SUPERSEDED | TRANSIENT（重試副本） | **NEEDS_DECISION（lean SAFE_DELETE）** | 1 |
| `big8_output_archive/multi_sample_quick_check/` | 6.50 GB | COLO829 chr19 快速檢查 + overlap/retest | SUPERSEDED | CONCLUDED | NEEDS_DECISION（lean ARCHIVE；COLO829 ⭐4） | 1 |
| `bip8/big8_output_archive`（其餘 rounds/snapshots/分析） | ~9.07 GiB | 早期 round + 原始快照 + F1 分析 | SUPERSEDED | CONCLUDED | **ARCHIVE**（封存歷史） | 1 |
| `synthesis/observation_workspaces/*loh_round1*`（3×） | 各 ~621 MB | LOH master dataset 入口（含 all_region_rows.tsv.gz） | CURRENT/PRE-FIX | LIVE | **KEEP**（主軸命脈） | 2 |
| `synthesis/kde_rerun_B_14combos/` | timeout(large) | KDE fix 全量重跑 14 組合 | CURRENT | CONCLUDED_POSITIVE | **ARCHIVE**（撐 validated KDE acceptance） | 2 |
| `synthesis/concluded/`（6 子目錄） | 202 MB | chr19/clairsto/loh 已收尾分析 | CURRENT | CONCLUDED | **KEEP**（收尾證據） | 2 |
| `multilayer_hp_benchmark/` | timeout(large) | self-phasing 修正前後 HP 多層比較 | PRE-FIX | LIVE（仍被引用） | **KEEP** | 2 |
| `InterSubMod_big7_runbook/` | timeout(large) | big7 遷移執行紀錄 | CURRENT/1.5 月未動 | CONCLUDED（遷移完成） | **NEEDS_DECISION**（遷移完成可降封存） | 2 |
| `big8_output_archive/pure_tumor_evaluation/` | 850 MB | HCC1395 純腫瘤評估（1 分鐘內 2 次重跑） | SUPERSEDED | CONCLUDED_NEGATIVE(TO-pure) | **SAFE_DELETE** | 1 |

### 2.2 research 區

| artifact | 大小 | purpose | trust_tier | conclusion | verdict | area |
|----------|------|---------|-----------|-----------|---------|------|
| `methyl_augmented_filter_phase2/`（全） | 87 MB | 甲基當 FP filter Phase2（LOSO）；**docs 37 引用最高** | DEPRECATED(⭐2 L4) | CONCLUDED_NEGATIVE | **NEEDS_DECISION**（G6 論文 methods 證據 — 整體勿刪） | 3 |
| └ `cycle2/data/` master_augmented TSV | 77 MB | per-sample 特徵注入中間檔 | DEPRECATED | TRANSIENT（可重生） | **ARCHIVE** | 3 |
| `v5_provenance_followup/T1_2_read_level_audit/vote_dump_*_genome.tsv.gz`（3 檔） | 2.07 GiB | baseline/V5/V3F 全基因組 read-level vote dump | PRE-FIX/CURRENT | CONCLUDED_POSITIVE | **NEEDS_DECISION→ARCHIVE**（可由 3 binaries 重生） | 4 |
| `paired_priority_bug_audit/phaseC_genome_three_way*/` | timeout(large) | baseline/V3F/V5/V6 × on/off × tp/fp 全基因組 per-region ISM | 混 | CONCLUDED_POSITIVE | **NEEDS_DECISION→ARCHIVE**（本區最大回收標的；需實測） | 4 |
| `clairs_to_verdict_pilot/data/step1_verdict_vs_seqc2.tsv` | 284 MB | 單檔 per-variant Verdict×SEQC2 join 中間表 | SUPERSEDED | CONCLUDED_NEGATIVE | **NEEDS_DECISION**（最大單檔 research 可回收；可重生） | 3 |
| `feature_layered_observation/data + figures` | 406+33 MB | G1-G10 分層觀察中間特徵 + 圖（LIVE skill 首次示範） | PRE-FIX | TRANSIENT（可重生） | **ARCHIVE**（建議用戶確認是否留示範） | 3 |
| `fp_provenance/.../*_master.tsv.gz` | 57 MB | 2.7M-row FP provenance master（pre-fix 口徑） | PRE-FIX | CONCLUDED_NEGATIVE | **ARCHIVE** | 3 |
| `seqc2_cnv_stratification/`（無 README） | 31 MB | SEQC2 CNV truth × zone/CN 分層；**docs 10 引用** | PRE-FIX | CONCLUDED_NEGATIVE（推得） | **NEEDS_DECISION**（無 verdict 檔 → 勿輕動） | 3 |
| `ng_kde_rescaling/`（無 README） | 35 MB | NG×KDE CN-tier rescaling；**docs 13 引用** | PRE-FIX | UNKNOWN（無 verdict 檔） | **NEEDS_DECISION**（結論狀態不明） | 3 |
| `coverage_multiple_validation/` | 1.7 MB | Coverage_Multiple 作 CN proxy | PRE-FIX | **blocked**（非 concluded） | **NEEDS_DECISION**（未結案，勿當可清） | 3 |
| `tsg_promoter_asm_reviewer/data/hg38.ncbiRefSeq.gtf.gz`+cpg | 40.6 MB | 公開 reference annotation（可重下載） | NA(reference) | NA | **ARCHIVE**（移共用 reference 或外部） | 4 |
| `tsg_promoter_asm_reviewer/msa.log` | 470 KB | MSA build log | NA | TRANSIENT | **SAFE_DELETE** | 4 |
| `autoresearch / data_registry / _template` | 小 | 活躍狀態機 + SoT + 範本 | CURRENT | LIVE | **KEEP（絕不動 — 不可 git rm 清單）** | 4 |

### 2.3 docs 區

| artifact | 大小 | purpose | trust_tier | conclusion | verdict | area |
|----------|------|---------|-----------|-----------|---------|------|
| `docs/trash/to_pipeline_staging_v1/` | 8.27 MB | v1 多階段 TO（**錯 VCF** F1=0.649）；README「已棄用待刪除」 | DEPRECATED | DEAD | **SAFE_DELETE**（v2 已確認存在於 research/） | 5 |
| `presentations/.../20260429_教授報告/`（base+v2+v3） | 34 MB | 4/29 教授報告，跨版本重複 figures | CURRENT(validated) | LIVE | **NEEDS_DECISION**（dedup ~7-11MB；validated 不自動刪） | 5 |
| untracked `*.standalone.html`（有 .md 伴生，5 檔） | ~4.84 MB | in_progress preview（可由 .md 重生） | in_progress | LIVE | **NEEDS_DECISION**（LIVE 研究線，保守） | 5 |
| PI-report `*.standalone.html`（**無 .md**，3 檔） | ~3.0 MB | PI 交付物，HTML 即本體（不可重生） | CURRENT(pi) | LIVE | **NEEDS_DECISION（KEEP-leaning）** | 5 |
| `docs/**/__pycache__`（2 dirs） | 120 KB | report-build 腳本 bytecode | TRANSIENT | DEAD | **SAFE_DELETE** | 5 |
| `docs/archive/` | 3.5 MB | 封存區本體 | SUPERSEDED | DEAD-but-archived | **KEEP**（它就是 archive） | 5 |

### 2.4 transient 區

| artifact | 大小 | purpose | trust_tier | conclusion | verdict | area |
|----------|------|---------|-----------|-----------|---------|------|
| `/tmp/claude-1067.bak/.../b29591ml3.output` | 13.7 GB | 跑飛的 per-base depth dump（chr1-4 `pos\tdepth`） | TRANSIENT | DEAD | **SAFE_DELETE**（最大單一立即回收） | 6 |
| `/tmp/glsdash` | 884 KB | goal-landscape dashboard 建構腳本（產物已 commit） | TRANSIENT | DEAD | **SAFE_DELETE** | 6 |
| `state/compact_snapshots/`（9 份） | 80 KB | precompact hook 自動快照 | TRANSIENT | DEAD | **SAFE_DELETE（保留最新 1-2）** | 6 |
| `/tmp/claude-scholar` | 9.6 MB | plugin git checkout（工具安裝物） | NA | NA | **NEEDS_DECISION**（非 workflow 輸出，交工具盤點） | 6 |
| `data/bam/*` + `data/ref/*` + `data/answer/*` | 0 byte | 全 symlink → big8 | NA | NA | **KEEP**（不佔 big7 空間） | 6 |

---

## 3. 跨 area 校正與邊界聲明（主 agent / 用戶必讀）

整合 6 區時發現 area 之間 3 處衝突，本檔以下列裁決為準，避免雙重計與假陰性：

1. **🔴 canonical tagged BAM 歸屬（Area 2 vs Area 6）**：Area 6 line 73/101 稱「canonical 全 ISM metadata、0 BAM >100MB」。**實測錯誤** — Area 2 用 `stat` 確認 `canonical/*/paired_*/*_complete_matrix/longphase_s/*_tagged.bam` 共 **14 檔 = 3,709,322,840,333 B（3.37 TiB）**，單檔 95-453 GB。Area 6 淺掃 maxdepth 4-6 未觸及 depth-6 longphase_s 或 find 提前 timeout → 假陰性。**本檔以 Area 2 為準，14 BAM 歸 output 現役區 KEEP**。manifest `tagged_bam_ready=false` 與磁碟現實脱鈎（BAM 後跑出未回寫 manifest）。
2. **3 個 synthesis TO BAM 歸 Area 6**：651.9 GB（278.3+278.3+95.4）由 Area 6 列管 NEEDS_DECISION；Area 2 只計其非 BAM 分析產物（避免雙重計）。本檔遵循此邊界。
3. **archive/ 歸 Area 1**：Area 2 的 `synthesis/research_rounds/archive/`（34 dirs）明確 defer 給 archive 區；Area 1 才是 bip8/big8/research_rounds/archive 的權威稽核者。本檔 archive 數字一律取 Area 1。

其他需注意 anomaly（不影響加總，但影響刪除安全性）：
- **2 個高引用 research 專案無 verdict 檔**：`ng_kde_rescaling`（13 引用）/ `seqc2_cnv_stratification`（10 引用）只有 scripts+data，無 README/manifest/00_。conclusion_status 無法自證 → **勿輕動，建議補 README 或用戶確認**。
- **`coverage_multiple_validation` 是 `status: blocked` 非 concluded** — 先驗清單誤列「暫緩」，實際未結案 → 不可當可清。
- **`legacy_partials` BAM 與 20260315 BAM 非 byte-identical**（差 1.47M）→ 兩次獨立 longphase run，非純副本 → legacy 標 NEEDS_DECISION 而非 duplicate-SAFE_DELETE（保守）。
- **OUTPUT_INDEX stale**：§5.6 `final_closeout/`、§5.4 `batch_runs/` 實測已不存在；kde_rerun_* / 2 個 colo829 round 未進 OUTPUT_INDEX。索引需更新（非清理範疇）。
- **`v5 vote_dump` .gitignore 註解 + 2026-05-10 proposal 嚴重低估**：兩處寫 ~26-40MB，實測 **2.07 GiB**（chr19 時代估值未更新）。

---

## 4. 「如何整理各種資料夾」操作建議

### 4.1 掃描紀律（避免 subdir 爆炸 / 卡住）— 對齊 CLAUDE.md §12

> **本次稽核最大教訓**：per-region ISM 輸出每 run 有**數十萬個 region 子資料夾**，`du` / `find` 遞迴對 root 或對 archive 樹**必 timeout**。

| 要做的事 | 正確指令 | 禁止 |
|---------|---------|------|
| 量某目錄總大小 | `du -sb <dir>`（單一）；大樹改 `du --max-depth=1`（仍可能 timeout → 接受 floor） | `du -sh */`（對 root 遞迴） |
| 找大檔 | `find <dir> -maxdepth N -size +100M -type f`（必加 maxdepth + scope） | `find . -size +100M`（無 maxdepth，掃 GB 級 data/output） |
| 大樹估體積 | `find <dir> -size +5M -printf '%s\n' \| awk '{s+=$1}END{print s}'`（只計大檔，floor 估計） | 對 archive 整樹遞迴 du |
| 廣搜引用 | **Grep/Glob 工具**（自動跳過 .gitignore 9 重目錄） | `grep -r . `（掃全 repo） |

### 4.2 三種整理動作的判斷

**(a) 直接刪（SAFE_DELETE）— 確定可重生 / 已棄用 / 純 transient**
- 條件：(i) README 自述「待刪除/棄用」；或 (ii) 純 build/log/pycache/depth-dump transient；或 (iii) 1 分鐘內重複執行的 NEGATIVE 結案重跑副本。
- 典型：`/tmp/*.bak` depth dump、`docs/trash/`、`__pycache__`、stub `pipeline.log`-only 目錄、`pure_tumor_evaluation`、`bip8 datestamp ISM 21dirs`、`tsg msa.log`。

**(b) 打包封存（ARCHIVE）— 撐已發佈結論的可重生中間檔**
- 條件：data/figures 是 `.gitignored` 可重生中間特徵檔，**但結論依賴它作為原始計算輸入**（結論本身已蒸餾進小 tracked 報告）。
- 動作：`tar` 打包 → 移到封存資料夾（見 §4.3）→ 確認後刪原檔。**不直接刪**（保留歷史回溯）。
- 典型：research `*/data/*.tsv`、`*/figures/`、`fp_provenance master.tsv.gz`、bip8/big8 rounds+snapshots、`kde_rerun_B`、tsg public reference。

**(c) 保留（KEEP）— 結論證據 / LIVE infra / baseline 本體**
- 條件：(i) 小 tracked 報告/manifest/summary（撐 INDEX/MEMORY/ledger）；或 (ii) LIVE 母專案 + 狀態機（autoresearch/data_registry）；或 (iii) paper baseline（14 canonical BAM）；或 (iv) PI 交付物（無 .md 的 standalone HTML）。
- **絕不動**：`autoresearch/ data_registry/ _template/`（不可 git rm 清單）、`data/` symlink、validated 終版報告、`docs/archive/`。

### 4.3 封存資料夾路徑建議

| 用途 | 建議路徑 | 理由 |
|------|---------|------|
| big7_disk_output 內 ARCHIVE | `big7_disk_output/_ARCHIVE_2026_06_01/` | 與 output 同 disk，tar 後同卷移動快；命名帶日期防混淆 |
| research/ 內中間檔 ARCHIVE | `big7_disk_output/_ARCHIVE_2026_06_01/research_intermediate/` | research data/ 已 .gitignore，移走不動 git；集中管理 |
| 跨 disk 冷儲存（最大宗，如 s-pure-pileup 770GB） | 外部 disk（如 `/big8_disk/_COLD_ARCHIVE/`）或用戶指定 | 770GB 級不宜佔 big7 主卷；若確認 canonical 已重生則可直接 SAFE_DELETE 免封存 |

> ⚠ 封存前務必 `du -sb` 確認來源大小 + tar 後 `tar -tzf` 驗證完整 + 比對 file count，**確認後才刪原檔**（archive-then-delete，逐項用戶確認）。具體指令範本見 `CLEANUP_PROPOSAL.md §5`。

---

## 5. 索引交叉引用

- 分層清理提案（SAFE_DELETE / ARCHIVE / NEEDS_DECISION / KEEP + 可回收小計 + 指令範本）→ **`InterSubMod/docs/data_specs/20260601_cleanup_audit/CLEANUP_PROPOSAL.md`**
- 6 區逐項詳細表 → `InterSubMod/docs/data_specs/20260601_cleanup_audit/area_{1..6}_*.md`
- output 唯一入口索引（需更新）→ `big7_disk_output/OUTPUT_INDEX.md`
- archive 三 INDEX → `big7_disk_output/{bip8_output_archive,big8_output_archive,synthesis/research_rounds/archive}/INDEX.md`

---

*稽核者：cleanup-audit 綜合 agent | 2026-06-01 | 唯讀整合，未刪除/搬移任何檔案 | byte 來自 6 區實測 du/stat/find；timeout 項標 unknown 未編造*
