---
title: 任務儀表板 + 資料/檔案結構整理 — 盤點、研究與 restraint-first 提案
date: 2026-06-02
type: governance / design-proposal
status: in_progress
evidence: 盤點層 = 8-agent workflow 實測（L1，附 source_path）+ 主 agent 親驗 4 條 load-bearing 事實（active.json / 2 skill / 76-task git / CURRENT_FOCUS 結構）；研究層 = WebSearch L2/L3
data_sources:
  - state/active.json
  - state/cycles/
  - research/autoresearch/evidence_ledger.jsonl
  - research/autoresearch/hypothesis_queue.json
  - docs/concepts/2026/05/20260529_project_tasks_01.json
  - docs/CURRENT_FOCUS.md
  - .claude/skills/cycle-state/SKILL.md
  - .claude/skills/research-dashboard/SKILL.md
workflow_run: wf_8ff0488b-3af (8 agents, 706s)
---

# 任務儀表板 + 資料結構整理 — 提案

> **一句話結論**：用戶要的「聚焦當前 + 背景任務」**最大價值不在做新儀表板**，而在補一個**根因落差**（G6/G1 主軸從未註冊成機械 cycle → `active.json` 空 → 任何讀它的 dashboard 都空）+ 一個近零成本的焦點區塊。底層**資料結構大致健全、不需重整**；真正的問題是 3 個「紀律/角色」缺口，全部低成本可補。新 HTML 儀表板降為「驗證既有工具不夠用後才做」的最後一步。
> **adversarial critique verdict**：`proposal_needs_trim`（原始 8-agent 提案 80% 工程花在次要的「碎片化統一」，吃掉了真實的「聚焦」需求）。本 doc 已套用 trim。

---

## 1. 用戶需求（4 則訊息整合）

| # | 原話要點 | 拆解 |
|---|---------|------|
| R1 | 加強儀表板「觀察與整理」，聚焦**當前任務 + 背景任務** | 焦點 vs 背景的清楚分離 |
| R2 | 先確認**如何整理統整與更新任務**、任務進度與關聯、現在資料狀態 | 盤點 + 整理/更新工作流 |
| R3 | HTML 顯示可大規模參考任務佈局功能，或**遊戲任務**整理顯示方式 | 顯示研究（含 game quest）|
| R4（偏好確認）| Hybrid 焦點+全景 / 全部統一範圍 / G6 主軸優先排序 / **結果譜系**（結果→子任務→結論→影響→相鄰任務可行性）| 已鎖定的偏好 |
| R5 | 同時驗證**檔案/數據/任務結構是否有更好的整理方式** | 資料架構評估（本 doc §4）|

---

## 2. 現況盤點（三層）

### 2.1 任務狀態資料層（state/ + research/autoresearch/）
盤點 agent 實測 + 主 agent 親驗：

| 檔 | 角色 | 現況 | provenance |
|----|------|------|-----------|
| `state/active.json` | active cycle 索引（≤5）| **0 active + 1 recently_concluded（cycle3 NEGATIVE）**，updated 2026-05-30 | 主 agent 親驗 ✓ |
| `state/cycles/{id}/state.json` | per-cycle 7-phase 狀態機 | 只有 **2 目錄**：cycle3(NEGATIVE) + Drill-2 fixture(合成) | 主 agent 親驗 ✓ |
| `evidence_ledger.jsonl` | 跨 cycle 結論 SoT（append-only，衝突時為準）| **75 entry**（主 agent 親數；近期含 2026-06-02 ASM-CN pilot + ASM magnitude PROVISIONAL）| 主 agent 親驗 ✓ |
| `hypothesis_queue.json` | 假說 backlog | **18 entry**（5 live = 2 adopted + 3 adopted_annotation；13 已 closed/rejected = 7 rejected + 6 closed）；H013-H018 於 5/31 closed | 主 agent 親驗 ✓ |

> 🔴 **§13.7 修正紀錄（2026-06-03，derive 腳本 Layer A 抓到）**：盤點 agent 原報「evidence_ledger 73 / hypothesis_queue **336 entry（18 live + 318 closed）**」**兩個數都錯** — ledger 實為 **75**；queue **336 是把 `wc -l`=335 行誤當 entry 數**，實際只 **18 個 entry**（H001-H018）。`scripts/build_focus_board.py`（讀真檔渲染）跑出 18/5/13 與盤點轉述衝突 → 親數證實腳本對。**教訓**：subagent 轉述數字未經親驗就寫進 doc = 同 §13 捏造風險；by-construction（derive 腳本讀真檔）是抓這類錯的機械防線。
| `state/schemas/` (8 檔) | JSON schema 驗證 SoT | 全 v1.0 鎖定，無 drift | 盤點 agent |

### 2.2 顯示層（既有，重疊！）
- **`/cycle-state` skill**（read-only，cycle 級）：讀 active.json + cycles/*，priority 排序 + phase/tier/gate/stale + 路由建議。← **但 active.json 空 → 渲染近乎空白**。主 agent 親讀 SKILL.md ✓
- **`/research-dashboard` skill**（read-only，專案級）：一頁式 = 方向 + 假說佇列 + 最近結果 + blocker + memory 統計。← **已覆蓋「現在在哪」大半需求**，只是 terminal 文字非視覺 HTML。主 agent 親讀 ✓
- **`harness_health.py` 6-燈儀表板**（HTML）：harness 自我稽核，非研究任務。
- **5/29 goal-landscape workboard**（`docs/reports/in_progress/.../20260529_..._dashboard_01.standalone.html`）：甘特 + 76-task DAG + localStorage 可編輯；**vanilla JS + inline、零 CDN/framework/build**（適合擴充）。
- **`fill_report.py`**（Layer A 反捏造 renderer）。

### 2.3 敘述/目標層（docs + memory）
- `CURRENT_FOCUS.md`（SessionStart 注入）：日期分段散文，**無固定「焦點/背景」區塊**（主 agent 親驗結構 ✓）；2 條 LIVE 主軸 = G6 LOH-phasing ⭐3 + G1 ZAR1L/BRCA2 ASM ⭐3。
- `MEMORY.md`：Active Research / Pending / Concluded 三分層（Pending 內有已 closed 項，與 queue 不一致）。
- `20260529_project_tasks_01.json`（**已 git commit** ✓）+ `.md` + `.html` 三件套：76 task / 6 goal / 11 blocker / 27 guardrail，凍結在 5/29（之後 4 個事件未反映）。

---

## 3. 核心發現 — 根因 = 註冊落差（非「缺儀表板」）

**親驗事實**：`state/cycles/` 無任何 phasing / ASM / LOH / G6 / G1 cycle 目錄。

→ G6（LOH-phasing 論文主軸）與 G1（ASM characterization）這兩條敘述層 **LIVE** 的主軸，**從未走 `/cycle-init` 註冊成機械 cycle**。所以：
- `active.json` 永遠是空的 → `/cycle-state` 永遠渲染空白 → 看起來「沒有任務在跑」，但實際研究是活的。
- 這是**因果倒置的根因**：不是「缺一個能釘住焦點的儀表板」，是**機械層根本沒有「焦點」可釘**（active.json 空）。
- **先補這個落差（1-2 個 `/cycle-init`）**，active.json 就反映真實 LIVE，既有 `/cycle-state` 立刻有東西可顯示 —— 可能讓 60% 的「新儀表板需求」直接消失。

---

## 4. 資料/檔案結構評估（R5 — 你問的「是否有更好結構」）

> **裁決：3 個機械 SoT 結構健全，不需重整。問題是 3 個紀律/角色缺口，全部低成本可補。**

### 4.1 用 5 lens 評（你上輪確認的標準）

| Lens | 現況評分 | 說明 |
|------|---------|------|
| ① 單一事實源/去冗餘 | 🟡 機械層好、敘述層冗餘 | 3 機械 SoT（ledger/active+cycles/queue）分工清楚、schema 驗證、衝突仲裁明確（ledger 為準）。**但**同一主軸又散在 CURRENT_FOCUS + MEMORY + 76-task json，口徑不一。 |
| ② 機器可讀 | 🟢 機械層佳 | 3 SoT 全 JSON + schema enum，可直接 parse。76-task json 也結構化。 |
| ③ 人可驗證 | 🟢 | tier 可回溯 ledger 行；evidence_ledger append-only + provenance。 |
| ④ 結果譜系可表達 | 🟢 **已存在，未被使用** | ledger 已有 `parent_cycles`（衍生關係）+ `corrects`（修正關係）+ `research_potential`（孵出子任務）欄位 → **結果譜系資料本來就在 ledger，不需新建 graph**（見 §7）。 |
| ⑤ restraint/遷移成本 | 🔴 主要風險點 | 任何「再加一份統一 tasks.json + HTML」會把描述同主軸的檔從 5 變 7，**加劇而非緩解碎片化**。 |

### 4.2 三個真實缺口（非「結構壞」）

1. **註冊落差**（§3）：G6/G1 不在 active.json。**修法：`/cycle-init`**（紀律問題，非結構問題）。
2. **76-task json 是「已 commit 的凍結快照」卻被當 live**：4 天未更（cycle3 conclude / ASM correction / queue hygiene / ASM-CN pilot 4 事件未反映）。**修法：明標 frozen-snapshot（handoff 用），不當 live 數據源；live 狀態回歸 3 機械 SoT。**
3. **無單一「焦點」surface**：焦點埋在 CURRENT_FOCUS 日期散文裡。**修法：CURRENT_FOCUS 頂部加固定區塊（見 §5 Step 1）。**

### 4.3 不要做（reframe 後）
- ❌ **不新增第 6/7 份 SoT**（tasks.json + HTML 各一）。
- ❌ **不把 derive ETL 塞進 `harness_health.py`**（會稀釋它「純 grep、無狀態、永遠 current」的價值；上游 schema 一變就 silent-wrong = §13 同構風險）。
- ❌ **不造新 status enum** 去合併 quest-state + Org-mode + phase（又一份 mapping = 新 drift 源）。
- ✅ 若真要 derive：必須 **ephemeral**（`.gitignore`、每次重生、標「derived, not SoT」），且只讀 3 機械 SoT，不讀已落後的 76-task json。

---

## 5. 顯示層提案（trim 後 — staged，restraint-first）

### Step 0（根因，最高優先）— 註冊 G6/G1 成機械 cycle
`/cycle-init` 把 G6 LOH-phasing（+ 視情況 G1 ASM）建成 `state/cycles/` cycle → active.json 反映真實 LIVE → **先跑既有 `/cycle-state` 看夠不夠**。**先驗證既有工具不夠用，才談新建。**

### Step 1（近零成本）— CURRENT_FOCUS 加固定焦點區塊
頂部加結構化區塊（SessionStart 已注入，零新檔、零 drift）：
```markdown
## ★ 當前焦點（pinned ≤2）
- G6 LOH-phasing ⭐3 | 下一步：formal Wilcoxon 全樣本 | prov:O
## ▸ 背景孵化（≤2 槽）
- ASM 全基因組 survey（chr5-22 跑中）
- HKU handoff commit（待決）
```
直接命中「聚焦現在 vs 背景」核心需求。

### Step 2（**僅當**你仍要瀏覽器視覺看板才做）— Option A 最小 HTML
擴充既有 5/29 vanilla workboard（**零新 framework / 零 graph 庫**）：頂部 Focus 大卡 + ≤2 背景槽 + Now/Next/Later + blocked-by 徽章 + tier/provenance badge + banner 顯示 gen_date↔ledger 末條 diff（反 staleness）。tasks.json 必 ephemeral。估工 ~半天。

### 已砍（critique 裁決）
- **Option B（skill-tree DAG）**：6-goal/數十節點對單人是 over-engineering；「哪些假說解鎖可跑」用 queue 過濾 + blocked-by 徽章就夠，不需 Mermaid/Cytoscape。
- **Option C（全沉浸 quest-log）**：3 天工 + 最易踩遊戲化 anti-pattern（celebration/集點/RNG/隱藏成就）違反 §13.7。**連提都收斂為一句「已評估為 over-build，不做」**。

---

## 6. 遊戲 quest 借鑑 — 折進 vs 丟棄

| 借鑑 | 處置 |
|------|------|
| Tracked-quest HUD pinning（釘當前焦點）| ✅ 折進 Step 1 焦點區塊 / Step 2 Focus 卡 |
| Main vs Side quest 分類 | ✅ 對應「焦點主軸 vs 背景孵化」|
| Objective checklist + 階段步驟 | ✅ 對應 cycle P0-P6 phase（已存在）|
| Thought-Cabinet 背景孵化 + WIP-limit ≤2 | ✅ 折進背景槽（機制性防分心）|
| tier = quest 難度加權（非活動量）| ✅ 沿用既有 ⭐ tier |
| 安靜 checkmark（無 celebration 特效）| ✅ 對齊 §13.7 未驗證滿足感戒律 |
| Skill-tree 三態依賴 graph | ❌ 丟（over-build，未證實需求）|
| celebration/集點/排行榜/RNG/隱藏成就/omniscient compass | ❌ 丟（anti-pattern，違反科學誠信）|

---

## 7. 結果譜系（你 Q2 的第五支柱）— 用既有 ledger 欄位，不新建

你要的「結果→子任務→結論→影響→相鄰任務可行性連鎖」**資料已經在 evidence_ledger**：
- `parent_cycles` = 上游衍生關係
- `corrects` = 修正關係（如 ASM magnitude correction 指向被修正條目）
- `research_potential` = 孵出的後續方向（= spawn 子任務）
- `decision` / `caveat` = 結論與影響
- 「相鄰任務可行性連鎖」例：G5 filter DEAD → H013-H018 closed（已在 queue `closed_reason` 記錄）

→ **顯示層只需 surface 這些既有欄位**，不需建新 graph 或新欄位。截斷句「相鄰任務可行性」即此 ripple，已有資料支撐。

---

## 8. Restraint flags（務必避免）
1. 勿引入需 build 的 framework（React/Vue/webpack）— 破壞既有 vanilla 離線單檔/可 grep/可版控。
2. 勿用 D3 手刻或 vis-network 力導向 — 數十節點用靜態足矣。
3. 勿把研究任務硬塞有日期甘特軸 — 研究多無確定 start/end，假精確 = fabricated 時程。
4. 勿新增常駐 hook/cron/daemon。
5. HTML localStorage 編輯**勿當 SoT** — 重蹈 `build_html.py` 覆寫 21 手改事故（[[feedback_self_phasing_PI_no_build_html_py]]）；若做 Step 2，pinned 偏好放獨立 focus.json，不放會被重生的 HTML。
6. 勿新增第 6/7 份持久 SoT。

---

## 9. Open decisions（待用戶拍板）
1. **路徑深度**：restraint 三段（Step 0→1→停看）/ 加做 Step 2 HTML / 直接做你早先選的 Hybrid 全景？
2. **根因**：現在用 `/cycle-init` 把 G6（+G1）註冊成機械 cycle？還是維持 docs+memory 追蹤、儀表板讀敘述層？
3. **76-task json 角色**：明標 frozen handoff snapshot（建議）/ 退役 / 維持 live 數據源（不建議）？

---

## 10. Lineage
- 盤點+研究：Dynamic Workflow `wf_8ff0488b-3af`（8 agents：3 inventory Explore + 3 web research + 1 synthesis + 1 critique），主迴圈 Opus 4.8 2026-06-02。
- critique verdict `proposal_needs_trim`，本 doc 已套用 trim（砍 Option B/C、根因升 Step 0、資料結構評估收斂為「不重整 + 補 3 缺口」）。
- 主 agent 親驗 4 條 load-bearing：active.json 空 / cycles/ 無 G6-G1 / 76-task json 已 commit / CURRENT_FOCUS 無焦點區塊 / 2 dashboard skill 存在且重疊。
- 對齊 [[feedback_harness_restraint_over_adoption]]（裝前 grep 既有覆蓋）+ §13 數據誠信 + §13.7 完成宣稱 gate。
