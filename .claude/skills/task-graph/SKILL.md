---
name: task-graph
description: 研究任務檢視面板機制 — 隨時追蹤確認任務狀況與細節的單一入口。巢狀 WBS 任務樹（parent 含括 + depends_on 先後）→ 自動產生「焦點控制台」standalone HTML：頂部現在進度(焦點任務+已完成順序+git/各AI分支快照)、左樹導航、右焦點面板(輸入→任務→後續鄰域流程圖 + L0重點/L1摘要/L2細節分層 + 文件超連結 + 問題/驗證生命週期 + 逐項判讀存 localStorage)。唯一真值 state/tasks/graph*.json，TASKS*.md/board*.html 自動產生勿手改。人工分派地圖非自動執行器。USE WHEN：「任務圖」「task graph」「任務看板/控制台」「現在進度」「追蹤任務」「誰卡誰」「可以跑哪些」「整體進度」「主任務/子任務」「程式流程圖」「輸入必需/可選」「缺哪些資訊」「git 狀況/各 AI 分支」、要隨時確認 subclonal reconstruction 進度與細節時。SKIP WHEN：cycle 機械狀態（用 /cycle-state）、專案假說看板（用 /research-dashboard）、可重現性 script→figure DAG（用 /pipeline-manifest）、檔案組織稽核（用 /data-audit）、harness 自我稽核（用 /harness-health）、純 build/commit/docs。
allowed-tools: Read, Edit, Write, Glob, Grep, Bash
user-invocable: true
---

# Task Graph — 研究任務檢視面板機制

薄包裝 `scripts/tasks/task_graph.py`。**唯一機械真值** = `state/tasks/graph*.json`；`TASKS*.md`（git-diff 留底）與 `board*.html`（檢視面板）**由 script 自動產生，永不手改**（單一源 → 不會變第三個漂移 SoT）。

**定位**：人工分派地圖 + 隨時追蹤面板。標出「現在處理什麼 / 誰可跑 / 誰被卡 / 缺什麼 / 各 AI 在哪條分支」，**按鈕仍由使用者按**（全自動研究迴圈是專案 tombstone，本層不碰）。

## 兩個圖檔（schema 向後相容）

| 檔 | schema | 形態 | 看板 |
|----|--------|------|------|
| `state/tasks/graph.json` | v2（flat `component` 分組） | 30 節點全 subclonal 主軸 | `tasks_board.html` |
| `state/tasks/graph_ism.json` | v3（巢狀 `parent` 樹 + `focus`） | ISM 小樹聚焦範例 | `board_ism.html` |

## 指令（`--graph` 預設 graph.json；輸出依檔名衍生）

```bash
python3 scripts/tasks/task_graph.py --graph state/tasks/graph_ism.json all          # 全跑（validate+check+ready+render+render-html）
python3 scripts/tasks/task_graph.py --graph state/tasks/graph_ism.json validate     # schema/循環/邊/io/parent 一致性（exit 1 on error）
python3 scripts/tasks/task_graph.py --graph state/tasks/graph_ism.json check         # provenance / drift / 細節缺口 findings（advisory）
python3 scripts/tasks/task_graph.py --graph state/tasks/graph_ism.json ready         # 列可立即開跑（todo 且 deps 全 done，排除 milestone 容器）
python3 scripts/tasks/task_graph.py --graph state/tasks/graph_ism.json render        # → TASKS_ism.md
python3 scripts/tasks/task_graph.py --graph state/tasks/graph_ism.json render-html   # → board_ism.html（檢視面板）
```

### 生命週期 helper（v5；改 graph 後記得 render-html）

```bash
# F1 工作隔離：記任務在哪分支/worktree（板上「分支↔任務」對照 + 碰撞旗標；commit 仍由 git_branch_commit_guard 把關）
... claim --task T-x --branch feat/foo [--worktree ../wt-foo]
# F2 AI 交接包：產 state/tasks/handoff/<id>.{json,md}（return-contract + 已驗證輸出路徑 + refs），其他 AI grep 即讀
... handoff --task T-x
# F3 新任務建立+確認：generic 任務（真假說驗證改走 inject-hypothesis→cycle-init 並 links.cycle_id）
... add-task --id T-x --title ... --parent T-p --kind compute --owner claude [--depends-on a,b] [--io-in 'name:req:ref,...'] [--io-out 'a;b'] [--headline ...] [--dry-run]
# F4 驗證結果（實際 RUN 走 /run-evaluator 或 /verification-loop；此處只記結果+evidence）
... verify-result --task T-x --pass|--fail [--ref <evaluation.json/報告>]
# F5 驗證失敗→建子修任務（半自動，--dry-run 先確認；postmortem 走 /structured-tech-report）
... new-problem --task T-x --summary '<失敗現象>' [--cause ... --options 'a;b' --postmortem InterSubMod/docs/... --reopen '<可重試條件>'] [--dry-run]
```
（前綴皆 `python3 scripts/tasks/task_graph.py --graph state/tasks/graph_ism.json`）

## 面板（board HTML）有什麼

- **① 焦點控制台**
  - **頂部現在進度**：📍現在處理（focus 任務+狀態+重點）｜進度計數｜**已完成（依依賴順序，可點跳）**｜🔀 **Git 狀況 / 各 AI 分支**（render 時快照 branch + 最近 5 commit + `git worktree list`，附時間戳）。
  - **左欄（sticky）任務樹**：▸ 三角形原生展開/收合；點**任務名稱**才聚焦。葉節點末標輸入就緒 ✓就緒/⧖待產出/⊘需外部備；⚠=有問題。
  - **右欄焦點面板**（點左即更新）：**L0 重點一行** + **L1 摘要** +「**輸入 → 任務 → 後續**」鄰域流程圖（●必需 ○可選；只顯示對應任務+需要輸入+後續任務）+ 產出/交接(其他AI讀) + ▸ **L2 細節折疊**（觀察材料 / 問題 現象-根因-選項 / 驗證+on_fail / **文件超連結 file://** → 開對應 .md/.html/cycle / notes）+ **判讀鈕**（確認無誤/待查/有問題，存 localStorage 可匯出，其他 AI 可讀）。
- **② 其他檢視（折疊）**：完整依賴 DAG / 驗證結果 / Ready / 細節缺口 / 問題與驗證回饋彙整。

## 隨時更新與追蹤

- **刷新 git/進度** = 重跑 `render-html`（static HTML 無法真即時；git 是 render 當下快照 + 時間戳）。
- **改任務** → 編 `graph*.json`（加節點 / 改 status / 補 io/links/headline/summary/problem/verify）→ **再** `render && render-html`。**勿手改 TASKS/HTML**。

## 單一所有權（與既有機制分工）

| 本層擁有 | 既有機制擁有（本層只用 links/ref 引用，不重造） |
|---|---|
| 含括(parent) / 依賴(depends_on) / 分派(owner) / status / missing_info / io(CWL required/optional) / 分層敘述(headline/summary) | verdict + 證據 + 敘述 → cycle `state.json` / `evidence_ledger.jsonl` / memory；**問題**→postmortem、**驗證**→`/run-evaluator`/cycle eval、**失敗**→NEGATIVE ledger（節點 problem/verify/on_fail 只連過去） |

## 加節點規則（編 graph*.json）

- `id` 用 `T-<大寫英數與->`；v3 用 `parent`（含括樹，null=根）；`kind` ∈ data/compute/analysis/writing/milestone（milestone=容器，不進流程圖/ready）。
- **required 的 task-input（`io.inputs` ref=task 且 required:true）必須同時列進 `depends_on`**（validate 會擋）；optional 輸入 `required:false` 且**不**進 depends_on（否則卡 ready）。
- 分層：`headline`(L0 一行) / `summary`(L1 2-3 句)；缺則 fallback（headline=title、summary=notes）。
- 問題/驗證：`problem{summary,cause,options,ref,fix_child}`、`verify{how,ref}`、`on_fail`、`observe[]`（ref 連既有 postmortem/evaluator，勿在 graph 內重造演算法）。

## 完成宣稱前（§13.7 / §13.0）

跑 `validate` 綠 + `check` 看 findings；HTML 由 `build_workstation.py` 注入，缺必填值會 §13-A refuse（exit 3）。status/數字先有驗證過真值才寫 graph，勿填預期值。

## 與 build_workstation 的關係

`render-html` 把 graph 轉成 verify-workstation spec → 呼叫 `.claude/skills/verify-workstation/tools/build_workstation.py`（繼承逐項判讀 localStorage + JSON/CSV 匯出 + §13-A 反捏造）。**新增 `extra_js`**（additive，一次注入不進序列化）載入焦點切換 JS。流程圖/樹/鄰域為**純手刻 SVG + 原生 `<details>`**（零外部依賴、離線、grep-able）。

## 診斷

- validate FAIL：看訊息（bad id / parent 不存在 / 循環 / required input 未在 depends_on）→ 修 graph。
- HTML refuse exit 3：某節點缺必填 metric（如 `狀態`）→ 補 graph。
- 互動失效（focusFlow undefined）：多半是 extra_js 與 build_workstation 全域 `const` 撞名（如 LSKEY）→ 勿在 extra_js 重宣告既有全域。
- DRIFT finding：節點 links.cycle_id 指向 active.json 已 stale 的 cycle → 走 `/conclude-research` 退役或更新引用。
