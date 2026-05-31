# Harness 全方位稽核 + 自我稽核儀表板設計（2026-05-31）

> **provenance**：Dynamic Workflow `harness-audit-2026`（run `wf_aa229291-ff8`，11 agents / 910K tokens / 138 tool calls / 8.5min）→ 主 agent L1 spot-verify。
> **方法**：6 路業界研究（restraint lens 內建）+ 3 路內部盤點（讀實檔）→ 獨立對抗 restraint critic（generator/evaluator 分離）→ 儀表板設計。
> **延續**：[[feedback_harness_restraint_over_adoption]]（上輪 21-agent，12 條建議存活 4 條）。本輪 19 條去重建議 → **0 個新外部工具**。

---

## §0 TL;DR — restraint 裁決

跨 6 研究軸 + 3 盤點軸去重 **19 條**：**ADD 6 / MODIFY 5 / ALREADY_COVERED 3 / SKIP_HYPE 2 / SKIP_RISK 3**。

**結論**：此 harness 8 個 orchestration pattern、6 個 verification pattern、observability SaaS、SDD/TDD/PRM 方法全部判 **ALREADY_COVERED 或 SKIP**（business-case 在 flat-rate / 單人 / 本地 repo 下塌掉，或正中 disk-full / rebuild-覆寫 / 繞過 exit-2 / always-load 四大已知失敗類）。**全部 6 條 ADD + 5 條 MODIFY 都是「修 latent bug / 補 trivial 缺口 / 修 doc drift」，不是裝工具。**

**最高 ROI 前 3（全 L1、全 trivial、全主動防已知失敗類）**：
1. 🔴 **wire `pre_tier_upgrade_check.sh`**（或明確降級為 advisory）— 關閉「⭐4/5 tier gate 被 4 個 SKILL.md 宣稱卻未執行」的研究誠信破口。
2. 🔴 **`headless-research` 加 `isolation: worktree`** — 封住唯一高 blast-radius agent（autonomous + Write/Edit 主樹 + 無限 Bash + 無 sandbox），直接防 rebuild-覆寫-21-commit 失敗類。
3. 🔴 **修 doc/agent drift**（Hard Gate 5→4、workflow 1→2、cycle3 state↔doc 衝突、3 個 agent drift bug、cache_telemetry pricing）。

---

## §1 Part A — 頂尖 AI 工作流模式（對應你 3 個核心問題）+ 你的覆蓋

> 全部 ALREADY_COVERED（L1 對磁碟資產驗證）。下表 = 「業界頂尖模式 → 你已用什麼覆蓋」。

### 軸 1：如何理解 AI 的產生（observability / 推理透明度）

| 業界模式 | 成熟度 | 你的覆蓋（L1） | 裁決 |
|---------|--------|---------------|------|
| LLM tracing SaaS（Langfuse / LangSmith / Braintrust / Helicone / Arize）| 高（企業 B2B）| `cache_telemetry.sh`(cache-hit%) + `subagent_completion_logger.sh`(token/cache) + 7 postmortem log | **SKIP_RISK**：self-host = always-on Postgres+ClickHouse（disk-full 類）；flat-rate 下 cost telemetry 邊際價值=0 |
| 推理過程透明化 / Step→Verify | 高 | `evidence_ledger` L1-L5 + `provenance-tier-audit` + AGENTS.md §6 Step→Verify | ✅ 已覆蓋 |
| OpenTelemetry GenAI semantic conventions | 中（spec 仍 experimental）| telemetry hook 已抓語意等價欄位（私有命名）| **SKIP_HYPE**：純命名對齊，無 export 需求前零現值 |

### 軸 2：如何確認 AI 工作正確、可驗證（這是你的 OVER-supplied 軸）

| 業界模式 | 成熟度 | 你的覆蓋（L1）| 裁決 |
|---------|--------|---------------|------|
| Generator/Evaluator 分離 | 高（Anthropic 3-agent harness）| **4 個 fresh-context read-only verifier**：`evaluator`(cycle 7-check) / `methodology-reviewer`(方法層) / `reproducibility-audit` / `security-reviewer` | ✅ 過度供給 |
| 對抗式 critic（refute-by-default）| 高 | evaluator + methodology-reviewer 對抗 prompt；本次 audit 即用獨立 critic | ✅ 已覆蓋 |
| Deterministic vs LLM-judge（Hybrid Norm）| 高 | `run-evaluator` 用 **deterministic 6-component 公式**（非 LLM 主觀打分）→ 正確閃過 LLM-judge self-preference bias | ✅ 已覆蓋（且選對 deterministic 半邊）|
| machine-checkable schema 驗證 | 高 | `kb_schema_check`(exit-2) gates knowledge/**；但 **state/cycles/** 的 8 個 schema 無 runner 驗證** | ⚠ 唯一殘留缺口 → 見 §5 #10 |
| Provenance / citation 驗證 | 高 | `evidence_ledger` + `citation-verification` + `provenance-tier-audit` + P-14 | ✅ 已覆蓋 |

> ⚠ critic 警告：再加通用 LLM-as-judge 會 import self-preference bias + 在 prose 欄位塌成 brittle keyword match。**勿加。**

### 軸 3：如何更好評估與計劃（你最深耕的軸）

| 業界模式 | 你的覆蓋（L1）| 裁決 |
|---------|---------------|------|
| Spec-Driven Development（Spec Kit / Kiro）| `research-loop` plan.json + `pre-decision-audit`（且多了 SDD 沒有的 Pre-registration）| ✅ 已覆蓋（加 SDD 會造成 dual-SoT 衝突，§9 已把入口 5→3）|
| TDD-for-agents | `cpp-change` PDD + `pre_commit_compile_check`(C++ TDD) | ✅ 已覆蓋（研究 script 上做 TDD = brittle degenerate assert，category-error）|
| 階層任務分解（ReAcTree）| 7-Phase + subagent return-contract + Dynamic Workflow | ✅ 已覆蓋 |
| Pre-Mortem / WRAP | `pre-decision-audit §2` Klein Pre-Mortem | ✅ 已覆蓋 |
| Cynefin + DACI/DECIDE | `confirmation-protocol §Cynefin` | ✅ 已覆蓋 |
| Self-reflection loop | 用 **external critic（run-evaluator + read-only verifier）** 取代 | ✅ 更優（單-LLM self-correct 業界證實有時退化）|

### 軸 0：orchestration 模式（8 個主流全覆蓋）
prompt-chaining(7-Phase) / routing(task_type_advisor + §8 表) / orchestrator-worker(research-orchestrator + parallel-* + workflow) / evaluator-optimizer(run-evaluator + verifier fleet) / reflection(scientific-rigor §8.3.1 + Opus 4.8 native) / plan-execute(research-loop + check-staleness + haiku tier-split) / Deep Agents(整個 harness 即是) / Dynamic Workflow(2 個 .js)。**全 ✅。**

---

## §2 Part B — agent / workflow 角色佈局調整

### B1. 18-agent 角色與權限矩陣（L1 萃取）

讀全部 18 個 agent .md 後的**真問題**（非重疊潔癖）：

**權限破口（建議修）**
- **`headless-research`**：`tools: Read,Write,Edit,Glob,Grep,Bash`（無限）+ `model:inherit` + **無 isolation 行** = 唯一 autonomous + 可寫主樹 + 無 sandbox 的 agent。安全完全靠 prompt 自律（「絕不 rm/C++/NO-GO」只是文字，Bash 機械層可 rm）。→ **加 `isolation: worktree`（一行）**。
- **`reviewer`**：`tools: Read, Write, Glob, Grep, Bash(python3:*)` — 帶 **Write**，但被分類為 fresh-context verifier；其餘 3 個 verifier（evaluator/methodology-reviewer/security-reviewer）都無 Write。半 generator 半 verifier = 業界警告的 self-confirmation anti-pattern。→ **去 Write，verdict 回主 agent 落地**。
- **`paper-miner`**：Bash 無限但只需 python3/curl → 收斂白名單（low priority，與其他一起改）。

**死 agent / drift bug（建議修）**
- **`literature-reviewer`**：50 處 `mcp__zotero__*` 引用，但你的 MCP 只有 knowledge/biorxiv/ensembl → **每次 zotero 呼叫都會 InputValidationError = 死 agent**。→ 標 deprecated 或改接 knowledge MCP。
- **`research-orchestrator`**：line 23 路由 `/research-ideation`（已更名 `/problem-framing-ideation`）。
- **`release.md`**：line 42 co-author `Claude Opus 4.5`（應 4.8）。
- **`researcher`**：description 承諾 PubMed/bioRxiv/Ensembl MCP，但 tools 只有 WebSearch/WebFetch（over-promise）。

**不要動（SKIP_RISK）**
- `parallel-benchmark` vs `parallel-analysis`：tools/isolation/model 相同、協議 90% 重疊，**但都能用**。合併 = 改 working contract + count drift + 失去 routing 清晰度，單人不值得。
- 4 個 claim-linter hook（evidence_level_lint / causal_claim_check / researcher_claim_evidence_check / terminology_guard）合併：working safety hook 重構有 regression 風險，flat-rate 下 latency 省下無意義。

**覆蓋缺口**：plan 軸無獨立 adversarial plan-critic（architect 自評自 plan），**但 skill 層 `/pre-decision-audit`（GO/PROBE/NO-GO + Pre-Mortem）已補足 → 不需新 agent。**

### B2. Dynamic Workflow 佈局
- 現有 2 個（`cross_sample_benchmark.js` + 本次 `harness_audit_2026.js`；§8 doc 仍寫 1 = drift）。
- 路由規則（§8）正確且是**最重要的業界級保護**：含 Hard Gate（C++ commit / 刪檔 / NO-GO / ledger 覆寫）的工作**絕不包進 workflow**（subagent 一律 acceptEdits + 繞過 exit-2 hook）。**勿放寬此規則。**

---

## §3 Part C — 當前任務盤點 + 權限劃分（L1）

> ⚠ 時序：CURRENT_FOCUS header 停在 **5/30**，evidence_ledger 有 **6 條 5/31 ASM entry 未回寫**。衝突時 ledger 為準（research/CLAUDE.md）。

### 當前 LIVE 主軸

| 主軸 | 狀態 | tier | 一句話 |
|------|------|------|--------|
| **① LOH-constrained phasing** | ✅ LIVE 論文主軸候選 | ⭐3 Grade B+→A | NG=2 Inner same-hap，n=7 全正向 Wilcoxon W=28 p=0.0078；剩 R-SELFREF 全 7-sample flag-on 負控（~25-50hr C++）升 full A |
| **② ZAR1L/BRCA2 ASM** | ✅ LIVE 但**已收斂為 characterization-only** | ⭐3 | BRCA2 HP-axis Δβ=−0.122；**5/31 ledger 進一步收斂**：ASM real but **non-directional + non-discriminative + coverage-modulated**（strong-ASM FP enrichment = regression-to-extreme artifact）。CURRENT_FOCUS 未反映 → 對外引用有 over-claim 風險 |
| **③ 甲基化當 FP filter** | ❌ DEAD | ⭐2 L4 | active.json `recently_concluded` = Phase2 Cycle3 P6 NEGATIVE（5/30）。direction EXHAUSTED |

### 3 處 task↔狀態不匹配（**全是文件同步滯後，非權限越界**）
1. **cycle3 衝突（最關鍵）**：CURRENT_FOCUS §D3 說 cycle3「P0→需 conclude（尚未）」，但 active.json 已列 `recently_concluded` P6_COMMIT NEGATIVE。→ active.json 已被 skill 更新但 CURRENT_FOCUS 文字未同步。
2. **ASM tier 語意過時**：⭐3 數字未變，但 5/31 ledger 已把含義從「方向 POSITIVE」收窄為「non-directional artifact regime」。
3. **stale queue**：H015/H016/H017/H018 status 仍 `queued`，但全建立於 filter LIVE 時期，filter 現已 DEAD → queue hygiene 缺口。

### 權限劃分裁決
所有 state/ledger/queue 寫入都有正確的 skill gate（conclude-research / cycle-init / inject-hypothesis）+ Hard Gate（C++ commit 編譯 / 長計算 / NO-GO）。**無授權越界**；問題集中在文件同步滯後 + 上述 agent 層權限破口（§2）。

### 每條 task → 建議路由
| Task | skill/agent/workflow | Gate |
|------|---------------------|------|
| 收 cycle3 + Phase2 archive | `/conclude-research`（唯一可寫 state.json）| NO-GO 判定 = Hard Gate 需確認 |
| LOH-phasing 升 A（R-SELFREF 7-sample C++ 重跑）| `cross_sample_benchmark.js` **或** `parallel-benchmark`；C++ 先 `/cpp-change`+`/methodology-audit` | 長計算需當輪明示；C++ commit 不可包進 workflow |
| ASM cross-sample（COLO829）| `/infra-ops`(disk)→`multi-sample-consistency` | disk + 長計算授權 |
| C-1 allele-not-HP BAM follow-up | `/feature-layered-observation`+`/auc-confound-guard` | 新假說先 `/inject-hypothesis` |
| 清 stale queue | `/pivot-direction` 降權 | 禁止直接 edit queue.json |
| 回寫 5/31 ASM + 解衝突 | `/provenance-tier-audit`（掃 orphan/over-claim）+ Edit CURRENT_FOCUS | read-only audit 無 Gate |

---

## §4 Part D — 自我稽核儀表板 + 固定觀察流程（核心交付）

> 判定：**可以整合**，且 **~80% 輸入已在磁碟（telemetry/state/ledger/INDEX），缺的只是「單一聚合 red/green 視圖」**。零新外部工具，擴充既有資產。

### L0 — 6 顆健康燈（每顆 = 一行可驗證指令，首跑即亮 4 紅 2 黃 = 證明有效）

| # | 燈 | 紅旗條件（單一指令）|
|---|----|---------------------|
| 1 | **Count drift** | `find .claude/skills -name SKILL.md\|wc -l` vs §3「45」；agents vs §8「18」；hooks vs §4「38」 |
| 2 | **Hard Gate 真實性** | 數 wired 且含 `exit 2` 的 vs §4 宣稱（**現紅**：宣稱 5，真 4；kb_sot_guard/verify_gate 是 advisory）|
| 3 | **Tier-gate 啟用** | `grep -c pre_tier_upgrade_check settings`（**現紅 =0**，但 4 SKILL.md 宣稱它執行）|
| 4 | **State↔doc 同步** | active.json `recently_concluded` vs CURRENT_FOCUS 待辦（**現紅**：cycle3 衝突）|
| 5 | **Ledger 新鮮度** | `tail -1 ledger` date vs CURRENT_FOCUS header（**現黃/紅**：5/31 vs 5/30 + ASM tier 未回寫）|
| 6 | **Queue hygiene** | `grep queued queue.json` 交叉 active.json concluded（**現紅**：H015-018）|

### 固定觀察流程（STEP 0-8，掛成 `/harness-audit --health` 或手動）
```
STEP 0 進入點：改了 skill/agent/hook · 改任務方向 · 月度 · >7天沒跑
STEP 1 COUNT     find SKILL.md|wc + ls agents/hooks/workflows  vs §3/§4/§8
STEP 2 GATE      settings hooks + 各 script grep "exit 2"       vs §4 宣稱
STEP 3 SETTINGS  model / effortLevel / permissions 條數          應 4-8[1m]/xhigh
STEP 4 STATE     active.json + cycle state.json                 vs CURRENT_FOCUS
STEP 5 LEDGER    tail ledger (date+tier)                        vs CURRENT_FOCUS header
STEP 6 QUEUE     queue status=queued                            vs concluded/DEAD
STEP 7 POSTMORTEM tail docs/postmortems/*.log (6 檔)            新 bypass/failure entry
STEP 8 RENDER    寫 state/health_snapshots/YYYYMMDD.json + diff 上次 → 重渲 HTML
紅旗總則：任一紅 → 卡片頂置 + L0 燈紅 + 給「建議 skill/動作」+ 檔案路徑
```

### 與既有資產整合（擴充 ≠ 新造）
| 既有資產 | 怎麼擴 |
|----------|--------|
| `harness-audit-2026` skill | 加 `--health` 子模式（純磁碟稽核，不做業界研究）|
| `skill_registry_sync.sh` | 現只查 skill/agent 計數 → 擴涵蓋 hook-event + Hard Gate 宣稱數 |
| `cycle-state` / `research-dashboard` / `provenance-tier-audit` skill | L0 燈 #4/#5/#6 直接吃其輸出 |
| `cache_telemetry.sh` | 修 Opus 4.7→flat-rate 註解，只留 cache-hit% |
| goal-landscape dashboard HTML | clone 骨架 + localStorage（天然「vs 上次 diff」）|

### 產出形式：**單一 standalone HTML（clone goal-landscape + localStorage diff）+ `/harness-audit --health` 產生器**；markdown companion 為 fallback；**不掛 PostToolUse hook**（17 hooks 已過載，health 該按需/月度跑）。

### ASCII 線框（L0）
```
╔════════════════════════════════════════════════════════════════════╗
║ InterSubMod HARNESS HEALTH    snapshot 2026-05-31   vs last (+4 red) ║
║ model claude-opus-4-8[1m] · effort xhigh                            ║
╠════════════════════════════════════════════════════════════════════╣
║  ●1 COUNT DRIFT      ●2 HARD-GATE TRUTH    ●3 TIER-GATE WIRED        ║
║  GREEN 45/18/38      RED claims 5, real 4  RED not wired (4 SKILL想) ║
║  ●4 STATE↔DOC        ●5 LEDGER FRESH       ●6 QUEUE HYGIENE          ║
║  RED cycle3 衝突     YELLOW 5/31>5/30      RED H015-018 stale        ║
║  ──────────────────────────────────────────────────────────────    ║
║  OVERALL ⚠ 4 RED · 2 YELLOW   → 點紅燈展開 L1 修復建議 ▼            ║
╚════════════════════════════════════════════════════════════════════╝
```

---

## §5 19 條 restraint 裁決全表

| # | verdict | area | tier | effort | 建議 | 已涵蓋/風險 |
|---|---------|------|------|--------|------|------------|
| 1 | **ADD** | R6 | L1 | trivial | wire `pre_tier_upgrade_check`（或降 advisory + 改 4 SKILL.md）| 未 wired，最高誠信破口 |
| 2 | **MODIFY** | R6 | L1 | trivial | 修 doc drift：§4 Hard Gate 5→4、§8 workflow 1→2、cycle3 reconcile | hook 確為 38（R6 的「39」是錯的，勿改 39）|
| 3 | **ADD** | agents | L1 | trivial | `headless-research` 加 isolation:worktree | 防 rebuild-覆寫失敗類 |
| 4 | **MODIFY** | agents | L1 | trivial | `reviewer` 去 Write | 對齊 3 個 read-only verifier |
| 5 | **MODIFY** | observ | L1 | trivial | `cache_telemetry` USD 標 flat-rate 無意義 | 只留 cache-hit% |
| 6 | **MODIFY** | agents | L1 | trivial | 修 research-orchestrator(/research-ideation)、release(Opus 4.5→4.8)、researcher(MCP over-promise)| 全 1 行 |
| 7 | **MODIFY** | agents | L1 | small | `literature-reviewer` 停用或改接 knowledge MCP | 死 agent（zotero MCP 缺）|
| 8 | **ADD** | agents | L1 | trivial | paper-miner Bash 收白名單 + 6 個無限-Bash agent pin TMPDIR | 防 /tmp disk-full |
| 9 | **ADD** | R6 | L1 | small | `/harness-audit --health` 聚合 skill（read-only grep 現成 log）| Part D 產生器 |
| 10 | **ADD** | eval | L1 | small | `skill_eval_check.sh`：state schema 對 golden fixture 驗證 | MEDIUM backlog，僅 machine-checkable 欄 |
| 11 | **ADD** | plan | L2 | trivial | plan.json 加 expected_files[]；§1 一行告知加「涉及檔案」| 官方 plan-mode best-practice |
| 12 | SKIP_HYPE | reflect | L3 | — | P3→P4 mid-check pitfall scorer | PRM 需 step-labeled 訓練資料，scenario-mismatch |
| 13 | SKIP_RISK | gov | L1 | — | 合併 4 claim-linter hook | working safety hook 重構有 regression 風險 |
| 14 | SKIP_RISK | agents | L1 | — | 合併 parallel-benchmark+analysis | working contract 改動不值得 |
| 15 | ALREADY | R1 | L1 | — | 8 orchestration pattern | 全映射既有資產 |
| 16 | ALREADY | R2 | L1 | — | 6 verification pattern | verifier 過度供給 |
| 17 | SKIP_RISK | R3 | L1 | — | observability SaaS（Langfuse 等）| 每個都撞 ≥1 已知失敗類 |
| 18 | SKIP_HYPE | R3 | L2 | — | OTel 欄位命名對齊 | 無 export 需求前零現值 |
| 19 | ALREADY | plan | L1 | — | SDD/TDD/PRM/Cynefin/premortem | 最深耕軸全覆蓋 |

---

## §6 落地優先序

- **🔴 P0（trivial，零風險，主動防失敗類）**：#1 tier-gate（需你選 wire vs advisory）· #2 doc drift · #3 headless isolation · #6 三 agent drift · #5 cache_telemetry 註解
- **🟡 P1（least-privilege + Part D 產生器）**：#4 reviewer 去 Write · #7 literature-reviewer · #9 `/harness-audit --health` + dashboard HTML
- **🟢 P2（backlog，非急）**：#8 paper-miner/TMPDIR · #10 skill_eval_check · #11 expected_files[]
- **任務側（用 skill 不手改）**：reconcile cycle3（`/conclude-research` 或 `/provenance-tier-audit` 核對）· 回寫 5/31 ASM tier · 清 stale queue（`/pivot-direction`）
- **不做**：#12 #13 #14 #17 #18（SKIP）· 0 個新外部工具

> 所有變動沿用 [[feedback_strategy_then_per_item_confirmation]]：逐項 diff preview + 原因，等 ack。
