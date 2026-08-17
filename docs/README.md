<!--
建立時間: 2026-01-11 10:00
更新時間: 2026-08-13（新增公開資訊全面校正、local-vs-live 狀態與完整稽核入口）
目標: 提供 docs/ 目錄的最新結構、命名規範與工作流程，並提供 AI 漸進查閱指引
處理範圍: docs/ 全目錄（archive/deep 歷史快照除外）
關聯檔案:
  - docs/CURRENT_FOCUS.md
  - docs/experiments/INDEX.md
  - docs/standards/20260228_文件命名與狀態管理規範_01.md
  - docs/standards/20260228_output軟連結與版本控管規範_01.md
  - docs/standards/20260303_文件盤點分類與歸檔流程規範_01.md
-->

# InterSubMod 文檔庫

> 🟢 **Canonical root（science authority 2026-08-01；handoff snapshot 2026-08-13）**：研究真值在 `/big7_disk/liaoyoyo2001/InterSubMod`；`/big8_disk/...` 是三月快照。08-13 snapshot 是最新導航／治理入口，但不改寫 08-01 frozen authority 日期。exact-PS×HP strict read-linkage L1 已 7/7 完成（W=85,621）；其後續 recurrence-allowed candidate analysis 有 98,955 final groups、71,955 read-AF ranked units 與 exact topology census，但 cellular clone count、CN/LOH-aware parent→child 與全樣本 biological tree 仍未建立。`research/20260710_layered_reconstruction_v2/current_layered_topology_v3_raw_all_v1.json` 是 historical pre-strict candidate census，不能當成最新 strict 結果。

## 快速導航

### 常用入口

- [當前目標](CURRENT_FOCUS.md)
- [**2026-08-13 完整研究資料與軟體交接 snapshot**](handoff/20260813_完整研究資料與軟體交接_01/00_INDEX.md) ← 新人／AI、資料治理、軟體 I/O、bip7/bip8 與 release gates 第一入口；尚非 release-ready
- [**Handoff bundles 索引**](handoff/README.md) ← 08-13 navigation + 08-01 frozen science authority
- [**Exact-PS 全 7 資料集最終 HTML observation**](/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/20260801_exact_ps_observation_report/all7_v1/report.standalone.html) ← 13/13 authority hashes、responsive/no-JS/A4 QA PASS；CN/LOH=`NOT_INTEGRATED`
- [**HTML builder／finalizer／QA 重生說明**](../research/20260801_exact_ps_observation_report/00_INDEX.md) ← Python 最後輸出階段與 fail-closed release 契約
- [**Exact-PS strict read-linkage 全資料集報告與完成層級稽核**](../research/20260723_production_exact_ps_strict_read_linkage/20260723_exactPS嚴格ReadLinkage全資料集報告_01/) ← L1 7/7；strict topology與clone/fusion 0/7
- [**Historical pre-strict LongPhase-S PASS candidate-tree census**](../research/20260710_layered_reconstruction_v2/20260714_LongPhaseS_PASS_sSNV共現與拓撲最新分析_01.md) ← 舊 grouping 的歷史比較，不是 current strict topology
- [**AI 啟動壓縮上下文與研究索引**](references/manual/20260424_AI啟動壓縮上下文與研究索引_01.md) ← 每次研究前快速掌握專案邏輯、重要數據、任務順序與待決策問題
- [**研究構想總索引**](concepts/2026/04/20260409_研究構想總索引_01.md) ← 研究大圖景、發展樹、理論基礎、論文大綱
- [研究歷史索引](experiments/INDEX.md)
- [研究方法與突破方向全域分析](references/20260324_InterSubMod研究方法與相關文獻突破方向全域分析_01.md)
- [五目標研究願景與 LOH 先導觀察策略](research/methylation_methodology/2026/03/20260326_InterSubMod五目標研究願景與LOH先導觀察策略_01.md)
- [LOH 盤點執行規格](plans/2026/03/20260326_LOH盤點執行規格_01.md)
- [方法學審查 closeout](methodology/20260324_方法學審查全域結論報告_01.md)
- [研究主線整合週報（含 phase 2 與 annotation）](reports/validated/2026/03/20260311_研究主線週報_20260305_20260311_phase2_annotation整合_01.md)
- [簡報專用入口](/big8_disk/liaoyoyo2001/InterSubMod/docs/presentations/README.md)
- [AI provenance 專用入口](/big8_disk/liaoyoyo2001/InterSubMod/docs/provenance/README.md)

### 研究構想入口（Layer 0.5）

- [研究構想總索引](concepts/2026/04/20260409_研究構想總索引_01.md) — 8 份構想文件的導覽入口
- [研究構想與理論基礎](concepts/2026/04/20260409_研究構想與理論基礎_01.md) — 生物學+技術+假設鏈
- [研究項目發展樹](concepts/2026/04/20260409_研究項目發展樹_01.md) — 6 分支 + 13 死路
- [待確認重要事項](concepts/2026/04/20260409_待確認重要事項_01.md) — P0 阻塞 + 隱含假設
- [論文寫作大綱](concepts/2026/04/20260409_論文寫作大綱_01.md) — 投稿策略+故事線

### 研究入口

- [研究主題入口](/big8_disk/liaoyoyo2001/InterSubMod/docs/research/README.md)
- [甲基方法學入口](/big8_disk/liaoyoyo2001/InterSubMod/docs/research/methylation_methodology/README.md)
- [五目標研究願景與 LOH 先導觀察策略](research/methylation_methodology/2026/03/20260326_InterSubMod五目標研究願景與LOH先導觀察策略_01.md)
- [LOH 盤點執行規格](plans/2026/03/20260326_LOH盤點執行規格_01.md)
- [突破方向全域分析](references/20260324_InterSubMod研究方法與相關文獻突破方向全域分析_01.md)
- [方法學審查 closeout](methodology/20260324_方法學審查全域結論報告_01.md)
- [5kHz 主實驗後主線藍圖](research/methylation_methodology/2026/03/20260307_5kHz主實驗與方法學驗證藍圖_01.md)
- [突破方向版執行計畫](plans/2026/03/20260307_純樣本甲基研究執行計畫_01.md)
- [2026-03-05 至 2026-03-10 研究主線週報](reports/validated/2026/03/20260310_研究主線週報_20260305_20260310_01.md)

### 手冊入口

- [參考資料入口](/big8_disk/liaoyoyo2001/InterSubMod/docs/references/README.md)
- [內部手冊入口](/big8_disk/liaoyoyo2001/InterSubMod/docs/references/manual/README.md)
- [研究推進與實驗觀察手冊](references/manual/20260307_研究推進與實驗觀察手冊_01.md)
- [研究報告 Agent 與 Skills 使用手冊](/big8_disk/liaoyoyo2001/InterSubMod/docs/references/manual/20260310_研究報告Agent與Skills使用手冊_01.md)
- [研究週報 PPTX 客製化設定與製作手冊](/big8_disk/liaoyoyo2001/InterSubMod/docs/references/manual/20260311_研究週報PPTX客製化設定與製作手冊_01.md)
- [個人 PPT 設計風格規範](/big8_disk/liaoyoyo2001/InterSubMod/docs/references/manual/20260311_個人PPT設計風格規範_01.md)

### 治理入口

- [docs 結構診斷與重整規劃](reports/validated/2026/03/20260311_docs架構診斷與重整規劃_01.md)
- [docs round 1 結構重整驗證報告](reports/validated/2026/03/20260311_docs_round1結構重整驗證報告_01.md)
- [docs round 2 結構重整驗證報告](reports/validated/2026/03/20260312_docs_round2結構重整驗證報告_01.md)
- [docs round 3 查詢入口與導航補強報告](reports/validated/2026/03/20260312_docs_round3查詢入口與導航補強報告_01.md)
- [外部參考入口](/big8_disk/liaoyoyo2001/InterSubMod/docs/references/external/README.md)
- `standards/`：規範與治理文件
- `archive/`：歸檔（含 `deep/` 歷史快照）

## AI Agent 建議起手式

```bash
scripts/analysis/check_ai_agent_readiness.sh
```

## AI 漸進查閱指引

本文件庫採用 4 層漸進披露架構，避免 AI 直接陷入大量細節文件：

| 層 | 文件 | 目的 |
|---|---|---|
| Layer 0 (地圖) | `docs/README.md`（本頁）| 全局導航，方向感 |
| Layer 0.5 (構想) | `docs/concepts/2026/04/` | 研究大圖景、發展樹、理論基礎、論文規劃 |
| Layer 1 (當下) | `docs/CURRENT_FOCUS.md` | 現在在做什麼、阻塞點 |
| Layer 2 (歷史) | `docs/experiments/INDEX.md` | 已試驗方向，成功/失敗總覽 |
| Layer 3 (細節) | 各 experiments/reports/solutions 文件 | 完整數據與推導 |

> **2026-03-24 後的使用方式**：`CURRENT_FOCUS.md` 與突破方向全域分析負責 live 主線；`experiments/INDEX.md` 只負責歷史索引與避免重踩失敗方向，不再直接充當近期任務清單。

### 建議查閱流程

**新任務起手式：**
1. 讀 `CURRENT_FOCUS.md` → 確認當前優先事項
2. 讀 `experiments/INDEX.md` → 避免重複已失敗的方向
3. 讀 `references/manual/20260424_AI啟動壓縮上下文與研究索引_01.md` → 快速取得整體邏輯、重要數據、任務順序、待確認決策矩陣
4. 讀 `references/20260324_InterSubMod研究方法與相關文獻突破方向全域分析_01.md` → 確認突破方向優先序
5. 讀 `research/methylation_methodology/2026/03/20260326_InterSubMod五目標研究願景與LOH先導觀察策略_01.md` → 確認最終研究願景、paired/TO 分線與 LOH 先導觀察定位
6. 讀 `plans/2026/03/20260326_LOH盤點執行規格_01.md` → 確認 LOH summary、decision ledger、figure 與 case panel 的固定輸出
7. 讀 `methodology/20260324_方法學審查全域結論報告_01.md` → 確認哪些方向已 closeout
8. 讀 `research/methylation_methodology/.../20260307_5kHz主實驗與方法學驗證藍圖_01.md` → 確認新的研究藍圖與樣本角色
9. 讀 `plans/2026/03/20260307_純樣本甲基研究執行計畫_01.md` → 確認當前執行階段與交付物
10. 讀 `references/manual/20260307_研究推進與實驗觀察手冊_01.md` → 對齊研究與紀錄流程
11. 執行 `scripts/analysis/check_ai_agent_readiness.sh` → 確認環境狀態

**特定問題查閱：**
- 找過去的解決方案 → `solutions/{topic}/`
- 找已驗證的實驗結果 → `experiments/validated/` 或 `experiments/finalized/`
- 找系統架構決策 → `architecture/` 或 `decisions/`
- 找研究歷史與方向優先級 → `experiments/INDEX.md`

## 目錄結構（現行）

| 目錄 | 用途 |
|---|---|
| `architecture/` | 長期架構設計文件 |
| `concepts/` | 構想與設計草稿 |
| `plans/YYYY/MM/` | 執行計畫與里程碑 |
| `reports/validated/YYYY/MM/` | 已驗證分析報告 |
| `reports/finalized/YYYY/MM/` | 最終結論與決策報告 |
| `presentations/validated/YYYY/MM/` | 簡報、版本檔、PDF 與 deck 資產 |
| `experiments/in_progress/YYYY/MM/` | 實驗草稿 |
| `experiments/validated/YYYY/MM/` | 可重現驗證結果 |
| `experiments/finalized/YYYY/MM/` | 實驗最終結論 |
| `solutions/{topic}/YYYY/MM/` | 問題解法與修復紀錄 |
| `research/README.md` | 研究主題總入口 |
| `research/{topic}/README.md` | 單一研究主題入口與邊界說明 |
| `research/{topic}/YYYY/MM/` | 專題研究文件 |
| `references/README.md` | 參考資料總入口 |
| `references/manual/README.md` | 內部手冊、agent/skill 規格與設定檔入口 |
| `references/manual/assets/` | profile、JSON config 等機器可讀設定檔 |
| `references/external/README.md` | 外部文獻與第三方指南入口 |
| `references/external/YYYY/MM/` | 文獻回顧、外部指南、第三方整合參考 |
| `provenance/ai_sessions/YYYY/MM/` | AI 對話 provenance 與執行報告 |
| `decisions/YYYY/MM/` | 重整與治理決策紀錄 |
| `archive/YYYY/MM/` | 一般歸檔 |
| `archive/deep/` | 歷史快照（保留原貌，不回溯改名） |

## 研究歷史與實驗索引

→ **[實驗總索引](experiments/INDEX.md)**：所有已嘗試方向的成功/失敗/建議後續

> ⚠️ **HISTORICAL／SUPERSEDED chronology**：下列是當時的工作狀態索引，不是 current
> claim registry，也不能繞過 2026-08-13 handoff 的 authority／superseded crosswalk。沒有在本索引
> 重驗 scope、denominator、producer 與 hash 的數字不得直接引用。

主要研究主題（依時間軸；歷史標籤）：
- 甲基化解析與 CIGAR 座標映射（2025-11）✅ 已完成
- 距離計算、聚類分析、Bernoulli 度量（2025-11 ~ 2025-12）✅ 已完成
- 統計顯著性分析（Fisher / PERMANOVA / Cramér's V）（2025-12 ~ 2026-01）— 歷史實作／測試里程碑；不代表 current full-scope science validation
- TP/FP 特徵富集分析與 F1 最佳化（2026-01）— historical F1=0.8481；本索引未重驗其 scope、denominator 或 producer，不得當 current baseline
- Subsample 混樣甲基化偏差分析（2026-02 ~ 2026-03）❌ NEGATIVE — tumor-normal 組織差異混淆
- Purity-Aware 策略驗證（2026-02 ~ 2026-03）❌ NEGATIVE — subsample 無法模擬純度效應
- 方法學審查 closeout 與突破方向 roadmap（2026-03-24）✅ 已收斂
- Phase 1A ML Read Classification（2026-03-25 ~ 2026-03-28）— historical paired-pure delta F1=+0.0112，屬 marginal observation、非最終結論
- 系統性觀察 O1-O15 + 因果鏈驗證（2026-03-31 ~ 2026-04-06）✅ ISM 定位轉向 characterization
- Self-Phasing 循環依賴與 PON-Only 修正（2026-04-02 ~ 2026-04-06）— historical「62% LOH 消失」觀察；本索引缺其分母／口徑，不得升格為 current confirmed rate
- Phase 2 A-D 程式碼實作（2026-04-12 ~ 2026-04-13）— code-structure milestone（Normal BAM／LOH／historical Subclone-named module），不等於 cellular-subclone architecture 或 science validation 完成
- Phase 2 全量驗證與分析（2026-04 ~）⏳ 進行中

## 研究啟動入口

### 現在該先看什麼

1. [當前目標](CURRENT_FOCUS.md)
2. [AI 啟動壓縮上下文與研究索引](references/manual/20260424_AI啟動壓縮上下文與研究索引_01.md)
3. [研究歷史索引](experiments/INDEX.md)
4. [研究方法與突破方向全域分析](references/20260324_InterSubMod研究方法與相關文獻突破方向全域分析_01.md)
5. [五目標研究願景與 LOH 先導觀察策略](research/methylation_methodology/2026/03/20260326_InterSubMod五目標研究願景與LOH先導觀察策略_01.md)
6. [LOH 盤點執行規格](plans/2026/03/20260326_LOH盤點執行規格_01.md)
7. [方法學審查 closeout](methodology/20260324_方法學審查全域結論報告_01.md)
8. [5kHz 主實驗後主線藍圖](research/methylation_methodology/2026/03/20260307_5kHz主實驗與方法學驗證藍圖_01.md)
9. [突破方向版執行計畫](plans/2026/03/20260307_純樣本甲基研究執行計畫_01.md)
10. [研究推進與實驗觀察手冊](references/manual/20260307_研究推進與實驗觀察手冊_01.md)
11. [研究主題入口](research/README.md)
12. [參考資料入口](references/README.md)

### 研究文件 Agent 與 Skills 入口

1. [主 Agent：InterSubMod 研究文件代理](/big8_disk/liaoyoyo2001/InterSubMod/.claude/agents/intersubmod-weekly-research-agent.md)
2. [Skill：研究脈絡整理](/bip7_disk/liaoyoyo2001/.codex/skills/intersubmod-context-synthesizer/SKILL.md)
3. [Skill：週報生成](/bip7_disk/liaoyoyo2001/.codex/skills/intersubmod-weekly-report-writer/SKILL.md)
4. [Skill：指令修正與偏好收斂](/bip7_disk/liaoyoyo2001/.codex/skills/intersubmod-report-prompt-refiner/SKILL.md)
5. [使用手冊](/big8_disk/liaoyoyo2001/InterSubMod/docs/references/manual/20260310_研究報告Agent與Skills使用手冊_01.md)
6. [週報專用指令手冊](/big8_disk/liaoyoyo2001/InterSubMod/docs/references/manual/20260310_研究週報撰寫指令與skill草案_01.md)
7. [本週研究主線週報](/big8_disk/liaoyoyo2001/InterSubMod/docs/reports/validated/2026/03/20260310_研究主線週報_20260305_20260310_01.md)
8. [最新整合週報（含 phase 2 與 annotation）](/big8_disk/liaoyoyo2001/InterSubMod/docs/reports/validated/2026/03/20260311_研究主線週報_20260305_20260311_phase2_annotation整合_01.md)
9. [個人 PPT 設計風格規範](/big8_disk/liaoyoyo2001/InterSubMod/docs/references/manual/20260311_個人PPT設計風格規範_01.md)
10. [個人 PPT profile 設定檔](/big8_disk/liaoyoyo2001/InterSubMod/docs/references/manual/assets/20260311_liao_research_ppt_profile_01.json)
11. [簡報專用入口](/big8_disk/liaoyoyo2001/InterSubMod/docs/presentations/README.md)

### 文件用途

1. `research/`：長期研究脈絡與高層策略
2. `plans/`：目前階段要做的任務與驗收條件
3. `references/manual/`：每次研究都要遵守的操作與紀錄流程、agent/skill 規格與設定檔
4. `references/external/`：外部文獻、方法與第三方工具整理
5. `provenance/ai_sessions/`：AI 對話過程與 provenance，不是正式結論主入口

## 命名規範

### 標準格式

```text
YYYYMMDD_主題_流水號.md
```

### 例外

1. 固定名稱：`README.md`、`CURRENT_FOCUS.md`、`INDEX.md`
2. 長期架構文件：`snake_case.md`
3. `archive/deep/` 歷史快照：維持原檔名

## 文件狀態規範

### reports/

1. `validated`：可重跑驗證，可供內部引用
2. `finalized`：最終對外口徑

### experiments/

1. `in_progress`：草稿/探索中
2. `validated`：驗證完成
3. `finalized`：整體結論定稿

## 建議工作流程

1. 新增文件：先決定狀態層（in_progress/validated/finalized）
2. 套用命名規範並填寫 metadata
3. 更新關聯報告索引
4. 需要歷史保留時移至 `archive/YYYY/MM/`

## 如何新增文件

1. 確定文件類型（實驗、報告、解決方案、計畫、AI 對話 provenance）
2. 選擇狀態層（in_progress / validated / finalized）
3. 套用命名：`YYYYMMDD_主題_流水號.md`
4. 填寫 metadata 區塊（見下方範本）
5. 若為實驗，同步更新 `experiments/INDEX.md`
6. 若解決了重要問題，在 `solutions/` 下補充紀錄

**Metadata 範本：**
```markdown
<!--
建立時間: YYYY-MM-DD HH:MM
目標: [本檔案的目標或用途]
處理範圍: [涵蓋的工作範圍]
關聯檔案:
  - [相關檔案路徑]
-->
```

## 2026-03-03 盤點輸出

1. 全專案檔案清冊：`reports/validated/2026/03/assets/20260303_repository_full_file_inventory_01.tsv`
2. 全專案目錄清冊：`reports/validated/2026/03/assets/20260303_repository_directory_inventory_01.txt`
3. docs 子樹清冊：`reports/validated/2026/03/assets/20260303_docs_file_inventory_01.tsv`
4. 本次歸檔待審區：`archive/2026/03/20260303_ai_sessions_raw_artifacts_pending_review_01/`

## 相容性說明

- 2026-02-28 後已啟用新分層與命名規範。
- 舊路徑可能已重整，請優先從 `reports/finalized` 與 `reports/validated` 查找。
