# InterSubMod Knowledge Base — Changelog

本檔案記錄 KB 架構、規範、重大內容變動。內容級別修訂追蹤於個別文件的 `last_verified` 欄位。

---

## [v0.8] — 2026-07-12（全方法／流程外部驗證與 producer closeout）

### Added / changed

- 新增`11_external_literature/14_20260712_full_method_process_external_validation.md`：整理raw-all→LongPhase-S→partial-read regional tree全流程、外部一致/差異/衝突、方法空白與P0–P2驗證路線。
- raw-all producer/receipt現況由pending更新為7/7 PASS；638,259 records、582,820 `_sc PASS`、32,184 rescue、17,444 remove、164,253,537 sidecar rows。
- 明確新增`partial-read-constrained, model-determined regional mutation-state tree`口徑：model-determined不保證full-span read，也不等於observed clone或唯一biological history。
- 補TreeClone、PhyClone、CNAqc、Nomura 2026、Lee 2026等最新/直接prior art；ASMS升2026 VOR，Blood review升fulltext-verified。
- fresh canonical downstream仍未7/7；sensitivity/comparison/final verifier/immutable receipt與biological truth均pending，scientific release維持NO-GO。

## [v0.7] — 2026-07-11（claim/search memory、最新文獻與 raw-all anti-regression）

### Added / changed

- 新增 `11_external_literature/12_20260711_paper_claim_search_and_memory_router.md`：以 stable `ISM-*` claim ID、section、evidence role、relation 與 citation readiness 路由外部 evidence。
- 新增 `11_external_literature/13_20260711_latest_research_theory_contract_delta.md`：收錄 Rosenski/McAuley/PCDH/Blood/Foltz 最新查證、rooted three-gamete 術語與 raw-all LongPhase-S contract。
- 現行 input contract 更新為 normalized paired ClairS raw-all → LongPhase-S `_sc.vcf` PASS；新增 `ISM-X04` 防止 PASS-only producer 回流。
- PCDH 由 abstract-only 升 fulltext-verified；Blood legal OA、environment retrieval、read depth 分開記錄；新增 Rosenski 2026 fulltext+source card。
- 外部 evidence semantic mapping 已全部逐卡 curated：97 cards / 20 claims / 201 links / 118 artifacts，strict audit 0 errors/0 warnings；禁止以外部來源直接 `supports` 內部 quantitative/engineering result。
- 修正 11–13 三份新 KB 文件的 frontmatter enum，統一為既有 `active` / `reference-summary` schema；補 `related_ids` 對稱連結，三個 KB validators 全數可重跑。

### Current gate

- HCC1954 chr22 與 patched HCC1395 whole-chrX raw-all probes通過；7/7 fail-closed producer/closeout、U0–U7 與 fresh downstream run未完成，`ISM-M03/R03`保持 pending。

## [v0.6] — 2026-06-05（新頂層群組：外部文獻地景）

### Added — `11_external_literature/`（10 份 .md + 1 standalone HTML companion）
- **7 角度交叉驗證外部文獻地景**：把「我們的研究目標/做法/困難/觀察」對應到外部文獻是否佐證/反駁/不同看法，每篇標真實 PMID/DOI + relation(CONFIRM/REFUTE/DIFFERENT-VIEW) + caliber 差異 + credibility；內部 L1-L2 / 外部 L3 清楚分離。
  - `00_index.md`（群組索引 + 夢的狀態死/活/開放）｜ `01`-`07`（Q1 長短讀方法時間軸 / Q2 找甲基差異位點 / Q3 ISM 方法對比 / Q4 subclone / Q5 分 TP/FP / Q6 甲基輔助 phasing / Q7 ASM-cis 癌症影響）｜ `08_cross_question_correlation.md`（理想情境鏈逐環 REAL/BLOCKED/POSSIBLE + 可信度衝突矩陣 + reviewer 防守）｜ `09_references_and_provenance.md`（~38 篇 + 識別碼核對 + 未覆蓋工具）。
  - `00_index.standalone.html`：PI/協作者視覺版（sticky TOC + 依賴鏈紅綠燈 + 7 角度折疊卡 + 衝突色標）。
- **來源**：workflow `wf_37b2cc97-663`（16 agents = 7 scout + 7 verify〔WebFetch 實證〕+ cross-synthesis + critic；2026-06-05）。內部數字沿用既有 validated SoT（master_draft / 研究故事 §7 / tsg snapshot），未跑新分析（§13 撰寫與查證分批隔離）。
- README 速查表 +1 列（外部文獻地景）。
- ⚠ 投稿前須 `/citation-verification` 核 `09 §B` 三個識別碼衝突（POG 三 DOI / epihet 作者-PMID / HiCancer venue）+ 補 `09 §C` 未覆蓋工具（WhatsHap/LongPhase/methylKit/Gaiti2019/amrfinder/MLML/COLO829 benchmark）。

---

## [v0.2] — 2026-04-23

### Added（2026-04-23 S5 PPT 製作過程發現資料陷阱）
- **[05_data_formats/06_merged_dataset_pitfalls.md](05_data_formats/06_merged_dataset_pitfalls.md)** ⚠ 新增
  - 陷阱 1：`merged_7samples_paired_full_plus_hcc1395_to.tsv.gz` 的 `AF` 欄位**非 caller_af**（p75 <0.1，正確 caller_af median 0.44-0.55）
  - 陷阱 2：HCC1395 `kde_status='phase1_new'` 的 LOH annotation 殘缺（Inner 5.7% vs archive 58.8%）
  - 含偵測腳本、根因、正確替代路徑（stale master archive）、修正範例
- 相關更新：
  - `05_data_formats/00_index.md`：新增 06 陷阱文件連結
  - `10_research_status/04_blockers_and_risks.md`：新增 B0 陷阱 blocker
  - `02_samples/01_hcc1395.md`：新增「phase1_new LOH 殘缺」警告區塊

### Context
發現過程：2026-04-23 S5 PPT 圖製作時，先發現 HCC1395 phase1_new 的 LOH Inner n=2,303 異常 → 改用 archive 後發現其他 5 樣本 AF 分佈全堆 <0.1 底部 → 追 merged 檔發現 `AF` 欄位非 `caller_af` → 統一所有 6 樣本從 stale master archive 讀 + 真 caller_af 後圖恢復正常。

---

## [v0.1] — 2026-04-22

### Added（首次建立）
- 建立 `knowledge/` 目錄樹：11 個頂層分組 + `scripts/`
- **Phase 1 骨架與規範**（9 份）：
  - `README.md`（5 秒查詢決策樹 + 快速導航）
  - `AGENT.md`（AI agent 使用 SOP）
  - `CHANGELOG.md`（本檔案）
  - `00_governance/`：00_index、01_frontmatter_schema、02_naming_conventions、03_cross_reference_rules、04_freshness_policy、05_query_protocol、06_update_workflow
- 初始 frontmatter schema（對齊 `/big8_disk/liaoyoyo2001/knowledge/` 規範）
- 5 種典型查詢情境 SOP（情境 A-E）

### Scope
- 本 KB 採 git 版控於 `/big7_disk/liaoyoyo2001/InterSubMod/knowledge/`
- 不接入 `/big8_disk/liaoyoyo2001/knowledge/` MCP server（靜態文件 + scripts 驗證）
- 策略：核心規格（pipeline/參數/truth）自撰 + 研究結論（landscape/experiments）用索引連結

---

## [v0.2] — 2026-04-22（Phase 2-6 全部完成）

### Added
- **Phase 2**（28 份）：03_pipelines（5）+ 04_parameters（6）+ 05_data_formats（6）+ 06_workflows（7）+ 08_truth_and_benchmark（4）
- **Phase 3**（14 份）：02_samples（8）+ 07_derived_features（6）
- **Phase 4**（17 份）：01_project_overview（5）+ 09_conclusions（6）+ 10_research_status（6）
- **Phase 5** 維運腳本（5 py 檔）：
  - `scripts/validate_frontmatter.py`
  - `scripts/check_related_ids_symmetry.py`
  - `scripts/check_canonical_paths.py`
  - `scripts/refresh_last_verified.py`
  - `scripts/fix_related_ids_symmetry.py`（auto-fix bidirectional）

### Fixed（Phase 6 多 agent 驗證後修正）
- **路徑錯誤**：`longphase-s` / `longphase-to-mod` 位置修正為 `/big8_disk/` 非 `/big7_disk/`（11 份文件）
- **行號錯誤**：
  - `RegionProcessor.cpp:806` → `:1066`（write_significance_summary）
  - `DistanceMatrix.cpp` 6 個距離度量行號全部更新（偏移 11-45 行）
- **閾值位置**：Cramér's V `Config.hpp:48` → `GlobalTest.hpp:32` (`min_cramers_v`)
- **Runtime override 標註**：PERMANOVA 999 struct default → 99 runtime（RegionProcessor.cpp:1573 override）
- **預設值修正**：`--threads` 預設 `1` → `16`（Config.hpp:43 struct default）
- **HPFineNGroups canonical filter**：補齊三條件 `NG=4 + AF<0.4 + NR≥80`（原漏 AF<0.4）
- **Research 路徑**：`research/loh_subclone_af_methylation/` → `loh_subclone_af/` + `loh_subclone_af_paired/`
- **HCC1395 BAM 實際路徑**：`/big8_disk/data/HCC1395/ONT_5khz_simplex_5mCG_5hmCG/HCC1395.bam`
- **7/7 → 5/7 精確化**：HPFineNGroups Part B 實際為 5/7 POS + 2/7 特殊
- **HPFineNGroups 機制更正**：從「subclone marker」改為「四 bucket occupancy / LOH-constrained phasing」
- **bidirectional related_ids**：58 個單向連結補回指（fix_related_ids_symmetry.py 批次處理）

### Validated
- `validate_frontmatter.py`：66/66 檔案通過
- `check_related_ids_symmetry.py`：66/66 文件雙向對稱
- `check_canonical_paths.py`：66/66 canonical_paths 存在
- 3 個多 agent 驗證（file:line 真偽 / 結論一致性 / 紅隊審查）全部完成

---

## [v0.3] — 2026-04-22（補齊用戶指出的缺漏 + ΔF1 SoT）

### Added
- **`06_workflows/07_cpp_change_pdd.md`** — C++ 修改 PDD 6 步驟完整協議（Hard rule；對應 `.claude/skills/cpp-change/SKILL.md` + CLAUDE.md「不可放寬的硬性規則」）
- **`00_governance/07_new_info_protocol.md`** — 新資訊補充與驗證協議（5 類別分類 × 3 層驗證 × SoT 規則 × multi-agent 驗證模板）
- **`03_pipelines/05_f1_baseline_canonical.md`** — ★ **F1 / ΔF1 SoT**（Single Source of Truth）
  - 含 CL-002 / CL-003 / CL-005 完整 provenance（樣本、method、運行條件、truth set、limitations、lock 日期、Chain Registry ID）
  - Reproducibility 規格（輸入 BAM/VCF 路徑、ISM 執行參數、下游 ML filter 配置）
  - Provenance 追溯鏈（→ `docs/CURRENT_FOCUS.md`、`10_Research_Chain_Registry.md#CL-002` 等 5 層權威）

### Changed（ΔF1 集中化）
- 7 份文件的 ΔF1 數字改為 SoT 連結，避免 22 處重複維護：
  - `03_pipelines/01_paired_full.md`
  - `03_pipelines/04_pipeline_comparison.md`
  - `03_pipelines/00_index.md`（加新檔入列）
  - `09_conclusions/01_positive_findings.md`
  - `10_research_status/01_current_focus_snapshot.md`
  - `08_truth_and_benchmark/02_f1_calculation.md`
  - `06_workflows/04_f1_benchmark.md`
  - `02_samples/01_hcc1395.md`
- `README.md` 快速導航表新增 3 個入口（PDD / F1 SoT / 新資訊協議）
- `00_governance/00_index.md`、`06_workflows/00_index.md`、`03_pipelines/00_index.md` 新增檔案列入
- 15 個 related_ids 雙向補齊（由 `fix_related_ids_symmetry.py` 自動化）

### Validated
- `validate_frontmatter.py`：69/69 檔案通過
- `check_related_ids_symmetry.py`：69/69 文件雙向對稱
- `check_canonical_paths.py`：69/69 canonical_paths 存在
- 文件總數：**72 md + 5 py**

---

## [v0.4] — 2026-04-22（補齊行為協議 + 腳本/報告索引）

### Added
- **`00_governance/08_hooks_and_automation.md`** — 10+ Hooks 配置（C++ commit Hard Gate、cpp_edit_guard、evidence_ledger_sync、md_link_check 等）
- **`00_governance/09_confirmation_protocol.md`** — ★ 確認時機協議速查：4 級暫停分類、兩種執行模式、20 個場景對照、Opus 4.7 subagent 觸發規則
- **`00_governance/10_think_before_code.md`** — ★ 實作前準則：假設陳述、2D 暫停矩陣（可逆性 override）、目標驅動驗證（Step→Verify）、高影響場景清單、首 turn 規格檢核
- **`06_workflows/08_analysis_scripts_index.md`** — 155 個 `scripts/analysis/*.py` 按前綴 6 大類索引 + 15 個高頻腳本詳細表 + 發現新腳本流程
- **`06_workflows/09_pptx_and_weekly_report.md`** — 3 個 skill（weekly-report / results-report / report）對照 + PPTX 規範摘要 + 7 Phase 週報流程索引 + 既有 manual 清單

### Changed
- `README.md` 快速導航新增 5 個入口
- `AGENT.md` 新增「行為協議快速入口」段落 + Opus 4.7 Prompt 策略段落 + MCP 知識庫對照段落
- `06_workflows/01_build_and_test.md` 新增 Git Workflow 簡版
- `00_governance/00_index.md`、`06_workflows/00_index.md` 同步
- 相關 related_ids 雙向補齊

### Validated
- `validate_frontmatter.py`：74/74 檔案通過
- `check_related_ids_symmetry.py`：74/74 雙向對稱
- `check_canonical_paths.py`：74/74 canonical_paths 存在
- 文件總數：**77 md + 5 py**

### SoT 原則遵循
- Confirmation protocol、PDD 6-steps、PPTX/週報皆採「本 KB 速查 + skill/manual 為權威」混合策略
- 避免複製 `.claude/skills/` 內容 → 未來 skill 更新時 KB 不需同步

---

## [v0.4.1] — 2026-04-22（多 agent 驗證後修正）

### Fixed（三個驗證 agent 發現的問題）

**🔴 事實錯誤**：
- `08_analysis_scripts_index.md` 腳本總數 155 → **151**（實測）；分類統計更新：`build_` 81→79、`run_` 17→7、`verify+validate+check+evaluate` 10→9、其他專題 22→~31
- `09_pptx_and_weekly_report.md` `build_pptx.py` 路徑錯誤 → 更正為 `scripts/analysis/build_weekly_report_pptx.py`（+ oral_v2 / LOH 變體）；註明 `docs/presentations/<run>/build_pptx.py` 為實例內一次性客製化檔
- 週報 **7 Phase → 8 Phase**（含 Phase 3.5 導演 Storyboard 審查）；PPTX 產生在 Phase 4 非 Phase 5

**🟠 SoT 違反修正**：
- `10_think_before_code.md` 拿掉 2D 矩陣副本（v0.4 版本已出現定義漂移：少「/ 影響結論」字樣），改為連結 `09_confirmation_protocol.md#-2d-暫停判斷矩陣` 為 SoT
- `09_confirmation_protocol.md` Hard Gate 清單分級註記：**CLAUDE.md 原生 (1-3)** vs **KB v0.2+ 擴充 (4-6)**
- 20 場景表擴充為 **25 場景**，標題註記「合併 skill 兩表 + CLAUDE.md 高影響 + v0.4 補強」；補入 CMakeLists/requirements/hooks/push 分類/.clang-format 場景
- `AGENT.md` 行為協議入口順序 10 → 09 → 07 改為 **09 → 10 → 07**（先政策、後方法）

**🟡 Literal 陷阱 / 盲測修補**：
- `README.md` 快速導航補入 `08_hooks_and_automation.md`（盲測 Q1「commit 被擋」原找不到）
- `README.md` 新增 Opus 4.7 AGENT.md 入口（盲測 Q7 原找不到）
- `09_confirmation_protocol.md` Hard Gate #5 `git push` 細分：受保護分支 / force push / feature branch 的差別處理
- `09` 長計算例外補灰色地帶範例（「看看 HCC1395 benchmark」無長度詞 → 仍暫停）
- `09_pptx_and_weekly_report.md` PPTX「禁斜體」加 rationale：CJK 字型可讀性
- `AGENT.md` 樣本名統一 `HCC1395_5kHz` → `HCC1395`（對齊 `02_samples/` 命名）

**🔵 Agent 發現但未處理（記錄於此供未來）**：
- `.claude/CLAUDE.md` line 209 誤植「AGENTS.md」（實為 `knowledge/AGENT.md`）— 屬 CLAUDE.md 問題，非 KB
- Agent B 誤搜 `/big8_disk/knowledge/`（實驗室全域 KB）而非本專案 KB — 驗證結果不可用，已由 Agent A 與 C 覆蓋

### Validated（再次）
- `validate_frontmatter.py`：74/74 通過
- `check_related_ids_symmetry.py`：74/74 雙向對稱
- `check_canonical_paths.py`：74/74 canonical_paths 存在
- 「155」出現次數：0（已全部替換）
- 「build_pptx.py（專案根類路徑）」出現次數：0（已更正）

---

## [v0.4.2] — 2026-04-22（業界標準比對後立即修正）

### Fixed（schema drift）
- `source_type` 擴充允許值 5→7，新增 `reference`（既有使用）與 `postmortem`（對齊 Google SRE blameless postmortem culture）
- `01_frontmatter_schema.md` 加「業界標準對照」章節：**Dublin Core / PROV-O / ADR (Nygard) / Diataxis** 四大標準 KB 欄位對照表
- `03_concluded_negative.md` 加 20 項關鍵字 → N 編號索引
- `06_merged_dataset_pitfalls.md` 連結 HCC1395 / blockers / master_dataset / per_region_files 雙向補齊

---

## [v0.5] — 2026-04-24（業界標準全面對齊 + 強化自動化）

### 🔧 Task 1: Schema rename source_type → content_nature
- **批量 rename 77 份 .md frontmatter** `source_type:` → `content_nature:`（對齊 Dublin Core `type` 延伸概念，避免與 Dublin Core `source` 語義衝突）
- `01_frontmatter_schema.md`：欄位表 / section heading / 業界對照表 / 反模式範例全面更新
- `validate_frontmatter.py`：schema constant `SOURCE_TYPE_VALUES` → `CONTENT_NATURE_VALUES`；加 backward compat warning

### 📋 Task 2: NEGATIVE 20 項全面 ADR 化（Nygard 5 欄）
- 新增 `09_conclusions/negative/` 子目錄 + 20 份 ADR 子檔（N01-N20）
- 每份 Nygard 5 欄：Status / Context / Decision-tested / Method / Result / Consequences + References
- 從 20 個 MEMORY `project_*.md` 萃取；48-67 行/檔；總計 1,077 行
- 舊 `03_concluded_negative.md` 改為**索引頁**：關鍵字快查 + 一句結論 + 連結子檔
- 對齊 [adr.github.io](https://adr.github.io/) + Google SRE blameless postmortem

### 🔔 Task 3: 3 個自動化 Hooks
- **`kb_schema_check.sh`**（PreToolUse/Bash）：`git commit` 含 `knowledge/` 變動時強制跑 3 驗證腳本；失敗 exit 2 阻擋
- **`kb_freshness_warn.sh`**（UserPromptSubmit）：`10_research_status/` 與 `hypothesis_queue_snapshot.md` `last_verified > 14 天` → 印提醒
- **`kb_sot_guard.sh`**（PostToolUse/Edit|Write）：寫入非 SoT 文件含 `+0.0112/0.9762/-0.0206/0.9650` → 提醒連結 SoT
- `.claude/settings.local.json` 加 3 hook 條目
- `08_hooks_and_automation.md` 同步擴充

### 📊 Task 4: evidence_ledger schema 擴充
- `03_evidence_ledger_format.md`：11 欄 → 13 欄，加 `operator` + `reviewer`（ELN accountability）
- `07_cpp_change_pdd.md` Step 6 範例 JSON 同步加兩欄
- 既有 19 cycles 豁免；新 cycles（2026-04-24+）必填
- 對齊 NIH reproducibility accountability

### Validated
- `validate_frontmatter.py`：**95/95** 通過
- `check_related_ids_symmetry.py`：**95/95** 雙向對稱
- `check_canonical_paths.py`：**95/95** canonical_paths 存在
- **總檔案：98 md + 5 KB py + 3 新 hook shell = 106 檔案**
- 3 hook 實測：`kb_freshness_warn` 靜默（all <14d 正確）；`kb_sot_guard` 偵測 + SoT 放行

### 業界對照升級

| 框架 | v0.4.2 | v0.5 |
|------|--------|------|
| ADR (Nygard) | 🟡（NEGATIVE 單薄）| ✅ **20 項全面 ADR 化** |
| Dublin Core | 🟡（加對照表）| ✅ **rename `source_type` → `content_nature` 避免衝突** |
| ELN | 🟡 | ✅ **evidence_ledger 加 operator/reviewer** |
| SoT | ✅ | ✅ + **pre-commit hook 強制** |
| Zettelkasten | ✅ | ✅ + **雙向由 pre-commit 保證** |
| Postmortem | ✅ schema | ✅ + pitfall 範例 |
| Diataxis | 🟡 | 🟡（tutorial 延後） |
| FAIR | 🟡 | 🟡（ELN 已補部分） |
| PROV-O | 🟡 | 🟡（JSON-LD 延後） |
| Antipatterns | 🟡 | 🟡（pitfall 模板延後） |

### 不做（本次 out of scope）
- Diataxis tutorial 入口
- FAIR creator/rights/language 欄位
- PROV-O JSON-LD 輸出
- backfill 既有 19 evidence_ledger cycles operator/reviewer

---

## 版本命名規則

- **v0.x** — Phase 1-4 逐步建置中
- **v1.0** — 全部 phase 完成且所有 scripts 驗證通過
- **v1.x** — 內容更新（不改架構）
- **v2.0** — 架構重大調整（例如：接入 MCP、新增頂層分組）
