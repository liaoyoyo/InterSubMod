<!--
建立時間: 2026-06-12
報告類型: 數據準確度 + 檔案/資料可尋性 改進稽核（盤點 + 討論）
任務類型: 整理盤點 + 流程改進討論（非研究啟動；HD-1 仍為唯一研究 gate）
狀態: 討論稿（待用戶確認要實作哪些改進）
data_sources: workflow wf_a3a61cc0-5a1 結構化盤點輸出 + 主回合 fresh ls/git 驗證(2026-06-12) + docs/reports/research_landscape/20260611_subclonal_reconstruction_paper_foundation_01.md
驗證方式: 所有 status=current 的關鍵事實經 2026-06-12 主回合 ls/git/wc/grep fresh 驗證（非憑記憶）
-->

# 數據準確度 × 檔案可尋性 改進稽核（盤點 + 討論稿）

> **框架**：盤點地圖（Inventory-Map）→ 缺口轉行動（Gap-to-Action）。L0 一眼 → L1 重點邏輯 → L2 細節表 → L3 溯源。讀到夠就停。
> **產出方式**：4 平行唯讀盤點員 + 1 完整性批判（workflow wf_a3a61cc0-5a1）→ 主回合 fresh 驗證爭議事實 → 本稿只寫**驗證過**的真值。

---

## L0 一眼結論

1. **本次盤點本身重現了要改進的問題**：盤點員（與批判）多處「憑記憶寫 status，沒去 ls/grep」→ 7 處事實與真值不符（含批判自己 1 處）。這**就是**「數據準確度 + 可尋性」最該補的流程缺口的活案例 → 對應改進 ①②（純紀律、low-effort、high-impact）。
2. **數據資產本身齊且可尋**（6 樣本 paired_full tagged BAM + somatic VCF + 甲基 + ISM TP/FP 輸出全部 fresh 驗證存在），但有 **3 個結構性風險**：(a) big8 normal BAM 是**他人帳號檔案的單點失效**；(b) 5 個 worktree 讓「append-only 權威 SoT」(evidence_ledger) **物理分叉成 99/15/15 行**；(c) stale 索引網把新 session 導向已 park 的舊主軸。
3. **記錄/格式標準分三層**：文件規範（完整）／機械 hard-gate（4 個已 wired、有效）／自動化採用（缺口——`data_sources` 全 repo 只 1 檔填、`fill_report.py` 零部署、口徑標籤完全沒有）。

---

## L1 重點邏輯（為什麼這些是真問題）

- **目標 A（準確度）的根因**：過去捏造/口徑混用/feature 誤解，共同根因 = 「**寫進去的東西沒回去檔案 grep**」。§13 三層防捏造已堵住「報告數字」這一類，但**盤點/狀態/索引類**事實沒有同等約束 → 本次再犯。
- **目標 B（可尋性）的根因**：資產在哪其實一次 ls 就有；真正讓人翻找/誤判的是 **(i) 不知道自己在哪個 branch/worktree**（5 並行副本）、**(ii) 多個 stale 索引互相引用成網**、**(iii) 同名 SoT 檔在不同 worktree 內容不同**。單修 CURRENT_FOCUS 解不了，因為新鮮層是孤島。

---

## L2-A 數據資產地圖（will_use 第一線 · fresh 驗證 2026-06-12）

**路徑文法（6 樣本一致）**：
```
/big7_disk/liaoyoyo2001/big7_disk_output/canonical/{SAMPLE}/paired_full/{DATE}_{SAMPLE}_paired_full_full_complete_matrix/
    ├── longphase_s/{SAMPLE}_tagged.bam        ← somatic-haplotag BAM（HP1/HP2/HP1-1/HP2-1/HP3 + MM/ML 甲基）
    ├── longphase_s/somatic_pass.vcf.gz        ← somatic callset（ClairS-TO）
    ├── intersubmod_tp/                         ← ISM TP 輸出
    │     ├── filtered_snv_tp/{24 chr 目錄}/    ← ⚠ region 在這層（非直接在 intersubmod_tp 下）
    │     ├── significance_summary.csv          ← 全域彙總（直放此層）
    │     └── debug/
    └── intersubmod_fp/                         ← ISM FP 輸出（判別性分析的另一半，勿漏）
```
- **DATE**：HCC1395=`20260314`（**僅 1 個 run，無版本混亂**）；其餘 5 樣本=`20260315`。
- **6 樣本**：HCC1395 / HCC1937 / HCC1954（乳腺）· H1437 / H2009（肺）· COLO829（黑色素瘤）。+`HCC1395_DORADO`（basecaller 對照，know_exists）。
- **真值 SoT**：`InterSubMod/docs/experiments/in_progress/2026/05/20260531_methyl_phasing_A0_assets/VERIFIED_RESULTS.md`（V1-V12 全數字）+ 同目錄 ~150 可重跑腳本（改 `--bam` 跑其他樣本）。

**🔴 SPOF — big8 normal BAM**：
- `/big8_disk/data/HCC1395/ONT_5khz_simplex_5mCG_5hmCG/HCC1395BL.bam`：**owner=`zhenyu112`（非本帳號）**、149 GB、mtime **2025-02-25**。V10 copy-clean 陰性對照 + cis-test 整條依賴它。對方清理/改權限 → 斷鏈且**我方無寫權挽救**。其他 5 樣本 normal：big8 下目錄**都在**（COLO829/H1437/H2009/HCC1937/HCC1954），但各自的對應 basecaller normal BAM 仍需逐一 ls 確認檔名。

## L2-B 過去數據準確度問題 → 防線 → 殘留缺口

| 問題類 | 現有防線（已 wired？） | 殘留缺口 |
|---|---|---|
| **報告數字捏造**（06-01 BRCA2 寫 0.572 真值 0.866 方向反） | §13 三層：A `fill_report.py`(by-construction)、B `number_provenance_check.sh`(✅ wired, validated/pi exit-2)、C audit 工具 | A 零部署、C 未融入 workflow；**非數字類事實（盤點/狀態）無同等約束 → 本次再犯** |
| **口徑混用**（Tmode/Paired/TO、read-instance/unique、merged AF≠caller_af、max-collapse vs 5mC-only） | 文件（data_specs 信任度 4 級 + pitfalls P-12） | **frontmatter 無 build_mode/build_date 欄 → 無機械標籤**；二次紀錄易 drift |
| **feature 字面誤解**（HPFineNGroups=phasing 非 methylation） | pitfalls P-10 + 規則「先讀 C++ 定義」 | plan.preconditions 未強制列 source_code_refs |
| **confound**（L2 collider / pooled OLS / n_reads / 空間自相關 / CN / ALLELE-axis germline 基線） | pitfalls P-01~P-09 + methodology-audit + scientific-rigor L4 | confound disclosure 靠人工 invoke，未自動稽核 |
| **stale-binary 下游**（expected_coverage hardcoded 75.0；pre-fix master 殘留） | check-staleness P2 PRECHECK（build_commit ≥ fix_commit）✅ | 部分 HPFineNGroups-依賴舊報告未回溯重標 |
| **tier 升級缺第 4 軌正交證據** | `pre_tier_upgrade_check.sh`(✅ wired exit-2) | **active.json tier 欄位髒資料**（見下）會讓 grep tier 漏判 |

## L2-C 導航/可尋性層（哪個查什麼 · stale 清單）

**新鮮 will_use 主入口**：
- `InterSubMod/docs/CURRENT_FOCUS.md`（live 主軸，06-11 更新）· `InterSubMod/docs/experiments/INDEX.md`（實驗 SoT，06-11）
- `InterSubMod/docs/reports/research_landscape/20260611_subclonal_reconstruction_{LAUNCH_READINESS,paper_foundation}_01.md`（✅ 當前 branch 都在）+ `InterSubMod/docs/concepts/2026/06/20260611_..._整合篇章_01.md`
- `InterSubMod/research/autoresearch/evidence_ledger.jsonl`（why-a-conclusion；當前 worktree **99 行**）
- `InterSubMod/docs/references/20260611_master_workflow_architecture_01.md`（how-to-execute + git 決策表）

**🟡 stale 索引網（know_exists；會把新 session 導向已 park 的 G6/G1）**：
- `state/active.json`（06-02 凍；仍列 G6/G1 為 active，**tier 欄位髒**：active=`"3"` vs concluded=`"⭐⭐ L4"`）
- `knowledge/10_research_status/{01..05}`（05-18，25 天，G6-era，標 2-週 TTL 已過期）· `knowledge/09_conclusions/`（05-18）
- `docs/reports/INDEX.md`（04-11，未含 06 新檔）· `research/autoresearch/research_direction.md`（指向已廢 Thread D）
- `MEMORY.md`（**溢出 34KB > 24.4KB 上限**，索引層本身在裂）

## L2-D 格式/記錄標準（已 wired vs 只有文件 vs gap）

| 標準 | 狀態 |
|---|---|
| 報告命名 `YYYYMMDD_描述_NN` / in_progress→validated 流程 / INDEX SoT 同步 | 文件化 ✅，**無機械強制**（orphan 偵測靠人工，05-29 曾出事） |
| §13 三層：`number_provenance_check.sh` | ✅ wired hard-gate（最堅實層） |
| `pre_tier_upgrade_check.sh` / `creation_guard.sh` / `skill_registry_sync.sh` / `health_drift_advisor.sh` | ✅ wired（後 3 advisory） |
| `data_sources:` frontmatter 慣例 | **gap**：validated/ 全 repo 只 **1 檔**填 |
| `fill_report.py`(Layer A) / `number_provenance.py audit`(Layer C) | **gap**：代碼齊全、**零報告部署** |
| build_mode/口徑標籤 | **gap**：完全缺欄位 |
| VERIFIED_RESULTS.md / DATA_PROVENANCE_LEDGER.md 模式 | 甲基案例✅ 示範，**無通用模板**、未見他處複製 |
| output 信任度 4 級 / 116 欄字典 / workspace 命名 / Mermaid 流程圖 | 文件化 ✅（參考手冊，know_exists） |

---

## L1 改進建議（ranked effort×impact；分 A 準確度 / B 可尋性）

| # | 改進 | 類 | effort | impact | 解什麼 |
|---|---|---|---|---|---|
| ① | 每份盤點/報告**首行標 `branch + HEAD + worktree`** | B | low | **high** | 消滅 cross-branch 幻覺（本次 80% 誤判源）；可進 doc-standards frontmatter `build_branch` |
| ② | **status=current 必附「驗證指令 + 一行輸出」否則自動降 unverified** | A | low | **high** | 堵「憑記憶寫 ≠ 去 grep」（§13.0 同源，擴到非數字事實） |
| ③ | **SoT 類檔跨 worktree 對帳**（ledger 99/15/15 分叉）→ 收斂單一權威副本 + 並行 session 不在各自 branch 維護 SoT | B | medium | **high** | 防 append-only 真值物理分叉（比 stale 嚴重一級） |
| ④ | big8 normal BAM 建「**外部依賴契約**」條目（owner/path/size/last-seen/我方無寫權 + checksum） | A | low | medium | 防 paper 關鍵證據 SPOF（zhenyu112 一年前檔） |
| ⑤ | **修 active.json tier 髒資料**（`"3"` vs `"⭐⭐ L4"` 正規化）+ 30 行 JSON schema lint 接 harness_health | A | low | medium | 防研究誠信 hard-gate 被髒格式繞過 |
| ⑥ | CURRENT_FOCUS pivot 時對 stale 索引頂部**自動插 redirect banner**（不重寫內容） | B | medium | medium | 把「誤導」降級為「明示過時+指路」 |
| ⑦ | frontmatter 補強制欄 `build_mode/build_date/data_sources` + 通用 template + `creation_guard` lint | A+B | medium | medium | 解口徑混用 + `data_sources` 採用率近零 |

---

## will_use vs know_exists 總表（給你的核心分類）

**🟢 WILL USE（研究啟動就會開、擺第一線）**
- 6 樣本 `{SAMPLE}_tagged.bam` + `somatic_pass.vcf.gz`（上方路徑文法）
- ISM `intersubmod_tp/filtered_snv_tp/` + `intersubmod_fp/`（判別分析兩半）
- `VERIFIED_RESULTS.md` + A0 ~150 腳本（真值 SoT + 可重跑）
- HCC1395 normal BAM（SPOF，V10/cis 對照）
- CURRENT_FOCUS / experiments-INDEX / foundation+LAUNCH_READINESS / 整合篇章（導航）
- evidence_ledger.jsonl（why）· master_workflow_architecture（how）
- §13 + number_provenance_check.sh（誠信）· known-pitfalls P-01~P-14 + confound catalog（準確度防線）
- data_specs 信任度 4 級 + Region 結構 spec（格式查表）

**🔵 KNOW EXISTS（知道有，後續再查）**
- HCC1395_DORADO（basecaller 對照）· path_inventory.tsv（stale，**勿信 SAFE_DELETE**）
- reports/INDEX.md（04-11 stale）· KB 10_research_status / 09_conclusions（05-18 G6-era stale）
- INDEX_DETAIL_ARCHIVE（深歷史）· postmortems（出事時查）
- master dataset 116 欄字典 / workspace 命名表 / Mermaid 流程圖（參考手冊）
- `fill_report.py` / `number_provenance.py audit`（工具在，未部署 → 改進 ⑦ 才啟用）

---

## 實作進度（2026-06-12，用戶選全做）

- **③ ledger 對帳結果（已驗證，免 merge）**：MAIN(99 行) 對 4 個 side-worktree（49/15/43/15）跑 `grep -Fxvf`：**side-only 行數全 = 0**，MAIN 反多 50/84/56/84 行 → **MAIN 是乾淨 superset，無內容分叉**（side-worktree 只是較舊快照）。批判的「永久分叉」是過度推論。**動作 = 政策非 merge**：並行 session 不在自己 branch 維護 ledger SoT；merge-back 一律取 MAIN ledger。**未覆寫 ledger**（Hard Gate）。
- **④ 外部依賴契約**：→ `InterSubMod/docs/data_specs/20260612_external_data_dependencies_01.md`。**研究級發現**：6 normal BAM 全 `zhenyu112` 擁有、4/6 是 symlink 指到 `/big8_disk/Google_somatic_data/bams/`（SPOF）。**甲基狀態 2026-06-12 samtools MM-tag 實測矯正**：先前憑 dir-name 推論「只有 HCC1395 有甲基→5/6 卡」**是錯的（P-17 活案例）**；實測 **5/6 normal 有甲基**（HCC1395 5mC+5hmC · HCC1937/HCC1954/H1437/H2009 5mC），**只 COLO829 R10 normal 無 MM** → G-A 對 **1/6 卡（僅 COLO829）**，5 樣本（乳腺3+肺2）可直接跑衝 ⭐4。foundation §5 G-A 應更新此實測前置。
- **①②⑤⑥⑦**：見下方各自落點（doc-standards / known-pitfalls P-17 / active.json+harness_health / redirect banner / frontmatter 欄位）。

## 待你確認

- 改進 ①②③④⑤⑥⑦ 哪些**現在做**、哪些**只記錄待後續**？
- 這份盤點稿要不要 commit 到當前 canonical 線（`docs/method-comparison-ism-external-202606`）？
- HD-1（研究啟動 gate）仍是你的決定，與本次「資料/流程整理」正交。
