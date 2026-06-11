<!--
建立時間: 2026-06-01
目標: Cleanup 稽核 Area [6/6] — BAM 清單 + transient 工作檔 + 跨 disk 大輸出（用戶核心問題：可移除的 BAM）
處理範圍:
  - InterSubMod/data/{bam,answer,vcf,ref}（symlink vs 真實大檔判定）
  - /big7_disk/liaoyoyo2001/big7_disk_output（= InterSubMod/output symlink target）真實 BAM/CRAM >100MB
  - InterSubMod/research 真實 BAM/CRAM >100MB
  - transient：/tmp/glsdash、/tmp/claude-*/.../tasks/*.output、state/compact_snapshots
方法: 唯讀 du / ls -la / scoped find(-maxdepth, timeout, ALLOW_FULL_SCAN) / scoped grep。絕不刪除/搬移。
-->

# Area [6/6]：BAM 清單 + Transient 工作檔 + 跨 disk 大輸出

> **scan_status**: PARTIAL_COMPLETE — 所有真實 BAM 已透過 targeted find 完整定位；兩個 archive 樹（bip8/big8）只能 maxdepth=5 淺掃（深掃 per-region ISM CSV 子目錄會 timeout），但淺掃 + manifest `tagged_bam_ready=false` 雙重佐證 archive 內**無大 BAM**。
>
> **用戶核心問題答案（可移除的 BAM）**：big7 上的真實大 BAM 全部集中在 `big7_disk_output/synthesis/research_rounds/.../step03_longphase_to/tumor_tagged.bam`（longphase-TO tagged BAM，**可由 source BAM 重跑重生**）。共 **3 個巨型 tagged BAM = 651.9 GB** + 1 個 normal subset 1.84 GB。`data/` 幾乎全是 symlink→big8（big7 上 0 byte），可回收量趨近 0。最大單一 transient = `/tmp/claude-1067.bak` 內一個 **13.7 GB 跑飛的 per-base depth dump**。

---

## 0. 關鍵發現摘要（先給結論）

| # | artifact | 大小 | verdict | 一句理由 |
|---|----------|------|---------|---------|
| 1 | `synthesis/research_rounds/legacy_partials/20260307_hcc1395_to_pilot_1_tagged_bam_only/step03_longphase_to/tumor_tagged.bam` | **278.3 GB** | NEEDS_USER_DECISION | HCC1395 TO legacy partial tagged BAM；被 20260315 complete run 取代但屬獨立歷史 tagged 狀態，且被多份 3 月 TO pilot 文件引用 |
| 2 | `synthesis/research_rounds/20260315_hcc1395_to_pilot/step03_longphase_to/tumor_tagged.bam` | **278.3 GB** | NEEDS_USER_DECISION | HCC1395 TO「完整」run 的 tagged BAM（legacy README 指此為權威來源）；canonical run 全 `tagged_bam_ready=false`，此為少數真 tagged BAM |
| 3 | `synthesis/research_rounds/20260423_colo829_to_pilot/step03_longphase_to/tumor_tagged.bam` | **95.4 GB** | NEEDS_USER_DECISION | COLO829 TO pilot tagged BAM；被 LIVE 的 LOH-phasing n7 Grade A + COLO829 ASM ⭐4 工作引用 |
| 4 | `hcc1395_normal_pilot/HCC1395_normal_subset.bam` | 1.84 GB | ARCHIVE | normal-read FP pilot 的 subset BAM；pilot 已 CONDITIONAL NO-GO，可由 full normal BAM(symlink→big8) 重 subset |
| 5 | `/tmp/claude-1067.bak/.../tasks/b29591ml3.output` | **13.7 GB** | SAFE_DELETE | 4 月一個跑飛的 per-base depth dump（chr1..chr4 `pos\tdepth`），在 `.bak` 備份目錄內，純 transient |
| 6 | `/tmp/glsdash` | 884 KB | SAFE_DELETE | 上次 session 的 goal-landscape dashboard 建構腳本 + 翻譯 chunk，純 /tmp transient，無 docs 引用 |
| 7 | `state/compact_snapshots` | 80 KB | SAFE_DELETE（保留最新 1-2） | `precompact_autosave.sh` hook 自動產的壓縮前快照；可重生、無 docs 引用 |

**3 個巨型 tagged BAM 全部 verdict = NEEDS_USER_DECISION 而非 SAFE_DELETE**，理由：(a) 雖可由 longphase-TO 重跑重生，但重跑成本高（每個數小時 + 需 source BAM 在 big8）；(b) canonical runs 全部 `tagged_bam_ready=false`，這 3 個是磁碟上少數真實 tagged BAM；(c) 皆被現行 docs 引用（含 LIVE 的 phasing 主軸）。**符合最高原則「結案 NEGATIVE ≠ 自動可刪；撐結論的原始證據至少 ARCHIVE」**。若用戶確認可重生且不再回查，#1（legacy partial，已被 #2 取代）是最安全的單一刪除目標 → 立即回收 278 GB。

---

## 1. `data/` — 幾乎全 symlink，可回收量≈0

`data/README.md` 明示「只放最小測試素材，大型 BAM/VCF/ref 用外部路徑」。實測：

| path | 大小(big7 實佔) | 性質 | verdict | 說明 |
|------|------|------|---------|------|
| `data/bam/HCC1395/{normal,tumor}.bam(.bai)` | 0 byte | **symlink → /big8_disk/.../*.bam** | KEEP | symlink，不佔 big7 空間 |
| `data/bam/H1437/` | 空目錄 | — | KEEP | 空 |
| `data/ref/*.fasta(.fai), hg38.fa*` | 0 byte | **symlink → /big8_disk/ref/** | KEEP | 4 個 ref symlink |
| `data/answer/{ClairS_output,ClairS_pileup,ClairS_pileup_longphase-s,phase_pos,SEQC2}/*` | 100 KB（全 symlink 條目） | **全 symlink → /big8_disk/** | KEEP | TP/FP/truth VCF + SEQC2 truth，皆 symlink |
| `data/vcf/HCC1395/pileup` | 0 byte | symlink → big8 | KEEP | |
| `data/vcf/COLO829_HKU_writable/*.tbi` | 732 KB | **真實檔（僅 2 個 .tbi index）** | KEEP | HKU handoff 的 TP/FP tbi index；極小、LIVE handoff 相關 |
| `data/test_snvs_32.tsv` | 68 byte | 真實 fixture | KEEP | 測試 fixture |

**結論**：`data/` 下無任何真實大檔；`du data/answer`=100K、`du data/vcf`=748K（其中 732K 是 COLO829 HKU tbi）。可回收量 ≈ 0。`data/bam` 全 symlink → 用戶問的「可移除 BAM」**不在這裡**。

---

## 2. 真實 BAM/CRAM >100MB 完整清單（用戶核心問題）

### 2.1 `/big7_disk/liaoyoyo2001/big7_disk_output`（= `InterSubMod/output` symlink target）

> 註：`InterSubMod/output` → `big7_disk_output` 為 symlink，兩 scope 同源，不重複計。

| path | bytes | 大小 | trust_tier | conclusion | verdict | reclaimable | referenced_by | rationale |
|------|------:|------|-----------|-----------|---------|------------:|---------------|-----------|
| `synthesis/research_rounds/legacy_partials/20260307_hcc1395_to_pilot_1_tagged_bam_only/step03_longphase_to/tumor_tagged.bam` | 278306777912 | 278.3 GB | SUPERSEDED | CONCLUDED (legacy partial) | NEEDS_USER_DECISION | 278.3 GB | `docs/.../20260307_HCC1395_5kHz_TO_pilot啟動與執行紀錄`, `20260308_*`(多份), `docs/standards/20260314_big7_canonical輸出規範`, `path_inventory.tsv` | legacy README 自述「不是完整可驗證 bundle，不應當主結果入口」；被 20260315 complete run 取代；可由 longphase-TO 重生但屬獨立歷史 tagged 狀態 |
| `synthesis/research_rounds/20260315_hcc1395_to_pilot/step03_longphase_to/tumor_tagged.bam` | 278305302826 | 278.3 GB | PRE-FIX | CONCLUDED (foundational TO run) | NEEDS_USER_DECISION | 278.3 GB | legacy README 指此為「完整完成版 TO run」權威來源；INDEX | HCC1395 TO 主 pilot 的 tagged BAM；canonical 全 `tagged_bam_ready=false`，此為真 tagged 來源；longphase-TO 可重生 |
| `synthesis/research_rounds/20260423_colo829_to_pilot/step03_longphase_to/tumor_tagged.bam` | 95388340944 | 95.4 GB | PRE-FIX | LIVE-adjacent | NEEDS_USER_DECISION | 95.4 GB | `docs/experiments/in_progress/2026/05/20260530_LOH_phasing_n7_cross_sample_Grade_A.md`(**LIVE 主軸**), `20260423_COLO829_TO_Append_Plan`, `20260414_LOH_Subclone_AF_Methylation_Evidence`(validated), INDEX | COLO829 TO tagged BAM；撐 LIVE LOH-phasing n7 Grade A + COLO829 ASM ⭐4；longphase-TO 可重生但成本高 |
| `hcc1395_normal_pilot/HCC1395_normal_subset.bam` | 1840301808 | 1.84 GB | PRE-FIX | CONCLUDED_NEGATIVE (normal-read FP pilot) | ARCHIVE | 1.84 GB | `docs/experiments/in_progress/2026/05/20260510_external_dirs_inventory_proposal` | normal-read germline-FP pilot 的 subset；對應 memory「Read-Level Germline FP CONDITIONAL NO-GO」；可由 full normal BAM(symlink→big8) 重 subset；非 load-bearing 結論證據 |

**Archive 樹（bip8_output_archive 39 dirs / big8_output_archive 5 dirs）**：
- INDEX 兩份皆標「全部 SUPERSEDED，已被 big7 canonical runs 取代，僅供極端歷史回溯」。
- `synthesis/master_run_manifest.tsv` 所有 archived run 皆 `tagged_bam_ready=false / blocking_reason=missing_tagged_bam` → archive 內**無 tagged BAM**。
- 實測：maxdepth=5 淺掃兩 archive 樹 → **0 個 BAM/CRAM >100MB**（exit 0 clean）。深掃 timeout（per-region ISM CSV 子目錄數十萬 inode）。
- 結論：archive 大在「inode/CSV 數量」非 BAM；BAM 回收量 = 0。archive 整體刪除屬其他 area 的 ISM-CSV 範疇，不在本 area BAM 清單。

**其他 big7_disk_output 子目錄（淺掃 maxdepth 4-6，皆 0 BAM）**：`canonical/`（19 runs，全 ISM metadata，`tagged_bam_ready=false`）、`multilayer_hp_benchmark/`、`v5_provenance_followup/`、`kde_smoke_test/`、`hcc1395_normal_pilot_global/`、`synthesis/concluded|kde_rerun_*` 均無 BAM >100MB。

### 2.2 `InterSubMod/research`（per 頂層子目錄 scoped find）

| 範圍 | 真實 BAM >100MB | verdict |
|------|----------------|---------|
| `research/*/`（全 43 個頂層子目錄 scoped find）| **0 個** | KEEP | research/ 內無真實大 BAM（BAM 皆 symlink 或不存在；研究產物為 .tsv/.csv/.md/.json/.png） |

→ 用戶問的「可移除 BAM」**完全不在 research/**。

---

## 3. Transient 工作檔（純可重生 → SAFE_DELETE）

| path | bytes | 大小 | 性質 | verdict | referenced_by | rationale |
|------|------:|------|------|---------|---------------|-----------|
| `/tmp/claude-1067.bak/-big7-disk-liaoyoyo2001-InterSubMod/a649aed0-.../tasks/b29591ml3.output` | 13718519808 | **13.7 GB** | 跑飛的 per-base depth dump（chr1..chr4 `pos\tdepth`） | SAFE_DELETE | none | 4 月一次 `samtools depth` 類指令 stdout 被整個截獲；在 `.bak`（舊 tmp 備份，Apr 29 快照）內；純 transient，無價值 |
| `/tmp/claude-1067.bak/` 其餘 | ~（13G 總量扣上項約幾百 MB） | — | 舊 claude tmp 備份（subagent jsonl symlink + task output） | SAFE_DELETE | none | 整個 `.bak` 是 active `/tmp/claude-1067` 的 Apr 29 過時備份副本 |
| `/tmp/glsdash` | ~905216 | 884 KB | goal-landscape dashboard 建構腳本 + gathered.json + 翻譯 miss_*.json chunk | SAFE_DELETE | none | 上次 session 一次性 /tmp 建構物；產物已落地 commit；無 docs 引用 |
| `state/compact_snapshots/` | ~81920 | 80 KB | `precompact_autosave.sh` hook 自動產的 9 份壓縮前快照 | SAFE_DELETE（建議保留最新 1-2） | none（scoped grep docs/scripts 無引用） | 可重生、自動產；保留最新 1-2 份做交接保險即可 |
| `/tmp/claude-scholar` | ~10066329 | 9.6 MB | claude-scholar plugin 的 git checkout（agents/commands/hooks/.git） | NEEDS_USER_DECISION | n/a | 非 workflow 輸出 transient，是工具安裝物；非本 area 範疇，標記供 area-config/工具盤點處理 |
| `/tmp/claude-1067/.../tasks/*.output`（active） | 多為 67MB cap | — | 當前 session 的 task output（含本次 audit 的 find/du）| KEEP（active session） | n/a | active session 工作檔，session 結束後自然清理；勿手動刪 |

---

## 4. Anomalies / 需在主 agent 注意的衝突

1. **`tumor_tagged.bam` 兩份 HCC1395 TO 非 byte-identical 重複**：legacy(Mar 11, 278306777912) vs 20260315(Mar 15, 278305302826)，差 1,475,086 byte → 兩次獨立 longphase-TO run，**不是純副本**。故 legacy 雖被取代，仍標 NEEDS_USER_DECISION 而非「duplicate → SAFE_DELETE」（保守）。
2. **canonical runs 全部 `tagged_bam_ready=false`**：磁碟上真實 tagged BAM 只有 synthesis 這 3 個。刪除前需確認重跑能力（source BAM 仍在 big8 + longphase-TO 環境可用），否則刪掉等於失去唯一 tagged 狀態。
3. **Archive 樹深掃不可行**：bip8/big8 archive 含數十萬 per-region ISM CSV 子目錄，`du --max-depth=1` 與遞迴 find 皆 timeout。本 area 已用 maxdepth=5 淺掃 + manifest 雙重佐證「archive 無大 BAM」；archive 的 ISM-CSV 整體去留屬其他 area 判定，非本 area。
4. **`/tmp/claude-1067.bak` 13 GB 在 /tmp**：呼應 memory `feedback_tmp_disk_full_pipeline_pitfall`（/tmp 寫滿災情）。此 13.7 GB depth dump 是 /tmp 膨脹源之一，SAFE_DELETE 直接緩解 /tmp 壓力。
5. **跨 disk 註記**：本 area 只掃 big7。3 個 tagged BAM 的 source（HCC1395/COLO829 tumor BAM）在 big8（symlink 指向）。刪 big7 tagged BAM **不影響** big8 source，重生路徑保留。

---

## 5. 本 area 可回收量加總

| 類別 | bytes | 大小 | verdict 分布 |
|------|------:|------|-------------|
| 3 巨型 tagged BAM | 652000421682 | 651.9 GB | NEEDS_USER_DECISION（可重生但需確認） |
| normal subset BAM | 1840301808 | 1.84 GB | ARCHIVE |
| transient（.bak 13.7G + glsdash + compact_snap） | 13800521944 | 13.8 GB | SAFE_DELETE |
| data/ | ~870000 | ~850 KB | KEEP（全 symlink/fixture） |

- **safe_delete（立即可清）**：13.8 GB（transient）
- **needs_user_decision（最大宗，須確認可重生）**：651.9 GB（3 tagged BAM）
- **archive**：1.84 GB（normal subset）
- **keep**：data/ ~0

> 若用戶確認 longphase-TO 可重生且不再回查歷史 tagged 狀態：3 BAM 全清 = **回收 651.9 GB**（big7 上最大單筆機會）。最安全先刪 = #1 legacy partial（已被 #2 取代）= 278 GB。
