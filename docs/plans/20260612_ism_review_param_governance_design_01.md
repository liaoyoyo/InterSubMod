---
title: ISM 程式碼庫 review + 參數測試治理 — 整體設計（C→A→B 三階段）
date: 2026-06-12
status: draft（待用戶 review）
type: design-spec / brainstorming-output
branch: chore/ism-review-param-governance-202606
worktree: .claude/worktrees/ism-review-infra
base: develop (069cadb，已含 T1/T2 修正)
provenance:
  - 偵察 workflow wf_c6c0412d-517（5 agents / 70 萬 token / 2026-06-12，read-only）
  - 既有理解 SoT：InterSubMod/docs/method_comparison/20260609_ism_vs_external_methylation_tools/01_ism_method_spec_from_source.md
  - 既有稽核 SoT：InterSubMod/docs/methodology/20260610_code_methodology_detail_audit_01.md
  - git 治理 SoT：InterSubMod/docs/references/20260611_master_workflow_architecture_01.md §3-§4
---

# ISM 程式碼庫 review + 參數測試治理 — 整體設計

> **本檔角色**：這個 session 的母設計文檔（design spec）。記錄偵察結論、已確認方向、三階段 roadmap、各階段 approach 選項與推薦。各階段細化 plan 另開子文檔。

---

## L0 一句話

> 這個 session 要為 ISM 建一套**可持續運作的基礎設施**：在「理解層已 85% 完備」的事實上，**先打穩驗證地基（C）→ 建參數測試治理規則+機械守衛（A）→ 補理解缺口（B）**。最大增量在 A（參數治理目前幾乎全空），C 是其前置地基，B 把理解推到「改任何部份都能推論影響」。

---

## §1 偵察結論（現狀地圖 — fresh 2026-06-12）

> 來源：偵察 workflow wf_c6c0412d-517。**數字均來自 agent 對實際檔案的 grep**；標 ⚠ 者為 agent 自承未實跑確認的靜態計數。

### §1.1 理解層 — 🟢 已 85% 完備（最大 reframe）

ISM 的「程式碼理解 + 正確性稽核」已是**源碼規格級 + 全碼稽核級 + 外部對照級**三層成熟，最新到 2026-06-09~11：

| 既有資產 | 涵蓋 | 深度 |
|---|---|---|
| `01_ism_method_spec_from_source.md` | 6 核心模塊到 **file:line**（CORE1 距離矩陣 / CORE2 clustering / CORE3 PERMANOVA / CORE4 per-CpG ASM / CORE5 CramersV gate / CORE6 NormalBaseline cis-test）+ Python 層（Δβ/cis-test/copy-partition） | 規格級（最深）|
| `20260610_code_methodology_detail_audit_01.md` | 12-reviewer workflow 全 C++ 窮盡稽核 **105 findings**，CRITICAL/HIGH 經對抗驗證（7 CONFIRMED/3 REFUTED/6 UNCERTAIN）| 稽核級。**L0：無翻轉結論方向的錯誤** |
| method_comparison study 13 份 + KB knowledge/ | ISM vs 82 外部工具對照、三條 pipeline 串接、衍生特徵語意 | 對照+概念級 |

**理解層真正的 gap（B 階段要補）**：
- 🔴 **S3-S9 pipeline 逐步走查未完成**（S1-S9 只做到 S2）— §Methods 正文骨幹缺口
- ⚠ 幾個模塊無專屬深度文檔：HierarchicalClustering 內部、io 層座標慣例、LohBedAnnotator
- ⚠ Python 分析層（181 腳本）無專屬方法文檔，僅策略性子集被稽核
- ⚠ 多數結論限單樣本 HCC1395 ⭐3

### §1.2 測試/回歸 — 🟡 單元測試強、端到端回歸弱

- ⚠ **218 GoogleTest cases**（14 .cpp → `run_tests`，秒級，ctest）— agent 標：靜態 grep 計數，未實跑確認全綠
- 3 個 pipeline shell（chr19 秒級 → 單樣本 TP/FP 分鐘 → 雙 VCF 批次）+ 驗證階梯 0-4
- 🔴 **6 個 core 模塊無 gtest**：SomaticSnv / MatrixBuilder / **NormalBaseline** / **SubcloneAnalyzer** / **LohBedAnnotator** / RegionProcessor 主流程（粗體為研究主軸相關）
- 🔴 **無 golden output / regression baseline 機制** — 改動如何影響真實結果只能 pipeline 端到端跑 + 人工肉眼比對 CSV/heatmap

### §1.3 git/測試治理 — 🔴 參數測試場景幾乎全空（最大真實 gap）

現有 6 個 git hook（3 commit-time exit-2 Hard Gate + 2 advisory + 1 notify）+ pre_tier_upgrade gate + 主工作流 §3 四決策表，**為「聚焦 feature 開發、WIP≤2-3」設計**。針對「為掃描各種參數而開大量實驗分支」**幾乎沒有對應件**：

| gap | 現狀 | 影響 |
|---|---|---|
| 實驗分支命名/註冊 | 只有 type/topic（feat/* chore/*），無 param 維度 | 大量參數分支會爆炸成難辨識 |
| 參數↔commit↔結果 綁定 | `binary_versions.jsonl` 只記「C++ 改了哪些檔」，不記參數值/結果 | 無法追溯「每個改動對結果的影響」|
| 結果隔離 | worktree 隔離 CODE 不隔離 RUNTIME，結果目錄共享 | 同輸入不同參數→輸出互相覆蓋（最直接污染風險）|
| config/param 失效鏈 | invalidation ledger 只在 C++ src/include commit 觸發 | 改 config/CLI flag 完全不進 ledger → /check-staleness 偵測不到 |
| 多 worktree 編譯 | `pre_commit_compile_check` 用全域 `/tmp/ism_cpp_pending_compile.txt` marker | 多 worktree 並行編譯 arm/clear 互撞，compile gate 可能誤放行/誤擋 |
| 實驗分支批次收斂 | 只有 clean_gone 清 remote-deleted | 無「一輪 sweep 結束→贏家 merge、其餘記結論後刪」流程 |

### §1.4 網路最佳實踐 — 已收斂成清晰藍圖（A 階段直接對齊）

權威共識（DORA / Martin Fowler / DVC / PLOS Ten Simple Rules / cmake-git-version-tracking）→ 4 主軸 + 7 步工作流：

1. **分支策略**：trunk-based 心態的**短命實驗分支 + commit-per-run**；探索用 branch、定案用語義化 tag；≤3 active branch
2. **參數管理**：**config-as-code**（params.yaml）、params+metrics 同檔版控、改參數自動失效重跑下游（DVC 模型）、每個實驗=可還原 commit 快照
3. **可重現性鎖三樣**：binary 版本（CMake 內嵌 commit hash/dirty）、環境依賴（pin+lock）、random seed（預抽 N 個版控）+ **commit hash 寫進每筆輸出 provenance**（PLOS canonical）
4. **改參數驗證**：**golden output / run-vs-baseline 回歸對照**，浮點用 **tolerance 比對**（double ~1e-10 / float ~1e-5）而非 bitwise

> 關鍵：這些與本專案既有 §13 數據誠信、`binary_versions/invalidation` 機制**同源** → 對齊強化而非全盤外購（符合 harness restraint 原則）。

### §1.5 環境事實（fresh git 確認）

- T1/T2 修正（`6593f96`/`eed4300`）**已在 develop + 當前研究 branch**，僅未到 main（main 落後）→ **C 階段降級為非 blocker**，發布時才處理
- 23 本地 branch + 5 worktree + **4 個近 5min 活躍 session** → 本 session 已隔離到專屬 worktree
- build/ 存在 → ctest 可跑（C 階段先 rebuild 確認新鮮）

---

## §2 已確認方向（決策記錄）

| 決策點 | 用戶選擇 | 日期 |
|---|---|---|
| 切入順序 | **D. 按 C→A→B 系統化全做** | 2026-06-12 |
| 治理深度 | **規範 + 完整機械守衛 hook** | 2026-06-12 |
| 工作環境 | **A. 新 branch + 新 worktree**（已落地）| 2026-06-12 |
| Roadmap 骨架 | **照 C→A→B 走**（逐階段細化 + 逐階段 review）| 2026-06-12 |

---

## §3 整體 Roadmap

```
階段 C 前置整備地基   →   階段 A 參數治理（規範+機械守衛）  →   階段 B 補理解缺口
（驗證地基可靠）           （未來調參有章法·不衝突·可追溯）       （改任何部份可推論影響）
   依賴：無                    依賴：C（golden baseline 機制）       依賴：A（理解放進治理結構）
```

每階段 = 1 個 sub-project，各自 spec → plan → 實作 → review。**不會一次全做完才停。**

---

## §4 各階段設計 + approach 選項

> 以下每階段列 deliverable + 2-3 個 approach + 我的推薦。**待 review 時請對各階段的 approach 拍板**（或微調）。

### §4.C 階段 C — 前置整備地基

**目標**：讓「改動→驗證」的地基可靠（A 的 golden 對照需要它）。

| 任務 | 內容 | approach 選項 | 推薦 |
|---|---|---|---|
| **C1 ctest 實跑** | rebuild + 跑 218 測試確認全綠，記錄 baseline 數 | （無選擇，直接做）| — |
| **C2 補關鍵模塊測試** | 6 個無測試 core 模塊 | (a) 全補 6 個 / (b) **只補研究主軸 3 個**（NormalBaseline·SubcloneAnalyzer·LohBedAnnotator）+ 其餘列 backlog / (c) 只排序不補 | **(b)** — ROI 聚焦研究主軸，YAGNI |
| **C3 golden baseline 機制** | 已驗證輸出存檔 + 改動後比對 | (1) 純 shell + byte diff（浮點假陽性）/ (2) **Python 比對 + tolerance**（double 1e-10/float 1e-5，獨立於 build）/ (3) C++ golden test 進 ctest（最整合但工程大）| **(2)** — ISM 輸出是 CSV 數值，容差比對必須且可獨立跑 |

**完成判定**：ctest 全綠 + golden baseline（chr19 + 1-2 anchor loci）可重跑比對通過。

### §4.A 階段 A — 參數測試治理（規範 + 完整機械守衛）

**目標**：未來調參不衝突·可追溯·可重現。

**deliverables**：
- **A1 治理規範文檔** — 對齊 web 7 步：`exp/<param>` 分支命名 · commit-per-run · config-as-code（params.yaml）· 鎖三樣 · golden 對照 · 收斂歸檔
- **A2 永久記憶** — feedback memory：參數測試 git 工作流（你明確要求存永久記憶）
- **A3 機械守衛 hook**（完整層）：
  1. 實驗分支命名檢查（`exp/<param>=<value>` 約定）
  2. **param/config 改動納入 invalidation ledger**（補「只記 C++」的洞）
  3. 結果目錄按 param 自動隔離（防覆蓋）
  4. binary 內嵌 commit hash + dirty flag（cmake-git-version-tracking 對接 binary_versions.jsonl）
  5. 修多-worktree 編譯 marker 互撞（全域 `/tmp` marker → per-worktree）
  6. golden baseline 回歸 gate
- **A4 實驗分支註冊表** — param↔branch↔commit↔結果 對照矩陣

**核心 approach 選擇（機械守衛建在哪）**：

| approach | 做法 | 取捨 |
|---|---|---|
| **(1) 全 in-house 擴充**（推薦）| 擴充現有 invalidation/binary_versions + 新增 hook，借鑒 DVC 的 params.yaml 慣例但不裝 DVC | 對齊 harness restraint、無外部依賴、與既有機制同源 |
| (2) 引入 DVC | config-as-code + DAG 自動失效全外購 | 業界成熟但 overlap 既有機制 + 學習成本 + 違 restraint |
| (3) hybrid | 核心 in-house hook + 真的裝 DVC 管 param DAG | 中間，但 DVC 對 C++ pipeline 非原生 |

> 推薦 **(1)**，理由：`feedback_harness_restraint_over_adoption`（裝前先 grep 既有覆蓋 + 對抗驗證）；本專案已有 invalidation/binary_versions 同源機制，擴充比外購好。**A3 每個 hook 落地前會先 grep 既有覆蓋避免重造**，且機械守衛屬 harness 變動 → 收尾跑 `/harness-health` 確認無 neutering/drift。

**完成判定**：dogfood 一個真實參數實驗，走完整 7 步流程（開 exp 分支 → 改 params.yaml → commit → 跑 → 輸出嵌 provenance → golden 對照 → 歸檔），全程被機械守衛正確護住。

### §4.B 階段 B — 補 ISM 理解缺口

**目標**：改任何部份都能推論影響。

**deliverables**：
- **B1 S3-S9 pipeline 走查補完**（接續既有 S1-S2）
- **B2 「改 X 模塊→影響哪些下游結果/特徵」因果影響地圖**（核心交付）
- **B3 Python 分析層方法文檔**（Δβ/cis-test/copy-partition）
- **B4 總導覽 index** 串起所有既有 + 新增理解文檔

**B2 因果影響地圖 approach**：

| approach | 做法 | 取捨 |
|---|---|---|
| (1) 人工逐模塊追 | 逐個讀源碼追下游依賴 | 最準但慢 |
| **(2) workflow fan-out**（推薦）| 多 agent 各追一個模塊的下游依賴 + 合成，**每條影響鏈對證源碼驗證** | 快（本 session 已證 workflow 偵察有效）+ 驗證把關 |
| (3) 半自動靜態依賴 | 從 include 依賴 + 既有文檔生成 | 快但漏動態/語意依賴 |

> 推薦 **(2)**，但 fan-out 結果**每條「改 X→影響 Y」必須對證 file:line 源碼**才寫入（不可只靠 agent 推測，對齊 §13.7 完成宣稱 gate）。

**完成判定**：對隨機選 3 個模塊，能正確推論「改它影響哪些下游」並被源碼驗證。

---

## §5 Dogfood：本 session 自身 = 治理規則 A 的第一個活案例

本 session 的工作流本身就在演示 A 階段要建的規則：
- ✅ 開 `chore/` 專屬 branch（主題≠研究主軸）+ worktree 隔離（4 並行 session）
- 將遵循：commit-per-unit（每完成一驗證單元即 commit、只 stage 自己檔）
- 將遵循：完成後跑 `/harness-health`（因動 hook）、push/merge 回 develop/main 需用戶確認

→ A 的規範文檔會把「本 session 怎麼做」抽象成可複用規則。

---

## §6 待 review 的決策點（請你拍板）

1. **§4.C 的 C2 補測試範圍**：(a) 全補 6 / **(b) 只補主軸 3** / (c) 只排序 → 我推薦 (b)
2. **§4.C 的 C3 golden 機制**：(1) shell diff / **(2) Python tolerance** / (3) C++ ctest → 我推薦 (2)
3. **§4.A 的機械守衛架構**：**(1) in-house 擴充** / (2) DVC / (3) hybrid → 我推薦 (1)
4. **§4.B 的因果地圖**：(1) 人工 / **(2) workflow fan-out** / (3) 半自動 → 我推薦 (2)
5. **整份 spec 是否可作為後續逐階段細化的依據**

> 預設（你若全同意推薦）：C2=(b) / C3=(2) / A=(1) / B=(2)。下一步進階段 C 的 writing-plans。

---

## §7 後續文檔結構（將建）

```
docs/plans/20260612_ism_review_param_governance_design_01.md   ← 本檔（母設計）
docs/plans/20260612_phaseC_pretest_groundwork_plan_01.md       ← C 細化 plan（下一步）
docs/plans/20260612_phaseA_param_governance_plan_01.md         ← A 細化 plan
docs/plans/20260612_phaseB_understanding_gap_plan_01.md        ← B 細化 plan
docs/governance/param_testing_git_workflow.md                  ← A1 治理規範（A 階段產）
scripts/hooks/<新 hook>.sh                                      ← A3 機械守衛（A 階段產）
```
