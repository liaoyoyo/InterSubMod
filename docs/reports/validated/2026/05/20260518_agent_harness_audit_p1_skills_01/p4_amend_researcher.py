#!/usr/bin/env python3
"""Final amend P4 HTML with researcher agent results."""
from pathlib import Path
from datetime import datetime

OUTDIR = Path("/tmp/skill_audit_20260518")
existing = (OUTDIR / "p4_industry_deep_audit.html").read_text(encoding="utf-8")

RESEARCHER_SECTIONS = """
<h2>🔬 Researcher Agent 深度補強（355 sec runtime, 12 一手來源 + 4 WebSearch）</h2>

<h3>關鍵新發現</h3>
<ul>
<li><strong>Anthropic 實際有 2 篇 harness 文章</strong>（之前以為只有 1 篇）:
  <ul>
    <li>「Effective harnesses for long-running agents」<strong>2025-11-26</strong> — <strong>2-agent</strong>（Initializer + Coding）+ Feature list 200+ items + Puppeteer MCP</li>
    <li>「Harness design for long-running application development」<strong>2026-03</strong> — <strong>3-agent</strong>（Planner / Generator / Evaluator）+ Sprint Contracts + Playwright MCP</li>
  </ul>
</li>
<li><strong>量化證據</strong>: 相同 Claude Opus 不同 harness Terminal-Bench 2.0 差距 <strong>16pp</strong>（77% vs 93%）— 「harness 品質與 model 同等重要」</li>
<li><strong>Anthropic 對 long-task 結論</strong>: 「agent 跑超過 3-5 步沒用 prompt caching 就是燒錢」</li>
<li><strong>Prompt caching 經濟學</strong>: 寫入 1.25× standard、命中 0.10× standard（90% 折扣）；多 user 共用 stable template + 1hr TTL 可達 80-95% hit</li>
</ul>

<h3>3 大業界 source 補強重點</h3>

<details open><summary><strong>OpenAI Codex Harness（6 個 Codex 特有 pattern）</strong></summary>
<ol>
<li><strong>Repository as System of Record</strong>: knowledge 放結構化 docs/，配 ~100 行 AGENTS.md 當目錄。「agent 在 context 看不到的等於不存在」</li>
<li><strong>Agent Legibility Over Human Preferences</strong>: codebase 為 agent 可讀性優化，偏好 boring composable tech 而非黑盒函式庫</li>
<li><strong>Architectural Constraints as Code</strong>: dependency 流向用 enforced sequence (Types → Config → Repo → Service → Runtime → UI)，custom linter 把「taste invariants」（structured logging、命名、檔案大小）寫成 rule</li>
<li><strong>Continuous Garbage Collection</strong>: 「golden principles」週期性執行檢測 drift，自動 refactor PR 維持 coherence（非人為週五大掃除）</li>
<li><strong>Telemetry-Driven Loop</strong>: agent 透過 LogQL/PromQL 直接驗證 perf；Chrome DevTools Protocol + ephemeral observability stack 內建</li>
<li><strong>Progressive Disclosure</strong>: agent 從小入口起步，按需發現更深 context（mimicry human onboarding）</li>
</ol>
<p><strong>OpenAI vs Anthropic 觀點差異</strong>：OpenAI 強調 declarative intent + autonomous loop；Anthropic 強調 stress test 假設、model 變強就移除。OpenAI 無 multi-agent role separation（single autonomous loop + PR-as-handoff）；Anthropic 三 agent 分工。</p>
</details>

<details><summary><strong>Walking Labs 12 Lectures 完整 takeaway（quantified）</strong></summary>
<table>
<thead><tr><th>#</th><th>標題</th><th>Takeaway</th></tr></thead>
<tbody>
<tr><td>L01</td><td>Why Capable Agents Still Fail</td><td>同模型同 prompt，bare 20min/$9 出 broken vs harness 6hr/$200 出完整產品；五層失敗（task spec/context/exec env/verification/state）</td></tr>
<tr><td>L02</td><td>What a Harness Actually Is</td><td>五子系統：<strong>instructions / state / verification / scope / lifecycle</strong></td></tr>
<tr><td>L03</td><td>Repository as System of Record</td><td>「agent 看不到 = 不存在」；cold-start test 5 問題；ACID state（atomic commits / consistent state / isolation / durability）</td></tr>
<tr><td>L04</td><td>One Giant Instruction File Fails</td><td>Progressive disclosure 替代巨型 doc</td></tr>
<tr><td>L05</td><td>Keep Context Alive Across Sessions</td><td>PROGRESS.md + DECISIONS.md + git commits + init/shutdown routine；量化：rebuild time -78%、feature completion 58%→100%、hidden defect 43%→8%</td></tr>
<tr><td>L06</td><td>Initialization Needs Own Phase</td><td>Agent 動工前先驗證 env health</td></tr>
<tr><td>L07</td><td>Overreach and Under-Finish</td><td>One feature at a time + 明確完成定義</td></tr>
<tr><td>L08</td><td>Feature Lists as Primitives</td><td>Machine-readable scope boundary，agent 不能忽略</td></tr>
<tr><td>L09</td><td>Declare Victory Too Early</td><td>Verification gap：confidence ≠ correctness；需獨立評估器</td></tr>
<tr><td>L10</td><td>E2E Testing Changes Results</td><td>只有全管線跑通才算驗證</td></tr>
<tr><td>L11</td><td>Observability Inside Harness</td><td>沒有 agent action 可視性，debug 不可能</td></tr>
<tr><td>L12</td><td>Every Session Leaves Clean State</td><td>下個 session 成功取決於本 session 清理</td></tr>
</tbody></table>
<p><strong>核心公式</strong>: <code>Harness = Instructions + Tools + Environment + State + Feedback</code>，五部件缺一即崩。</p>
</details>

<h3>🌐 8 個結構同構的開源案例（InterSubMod 可直接 cross-reference）</h3>
<table>
<thead><tr><th>#</th><th>Repo</th><th>核心特徵</th><th>對 InterSubMod 的 parallel</th></tr></thead>
<tbody>
<tr><td>1</td><td><code>anthropics/cwc-long-running-agents</code></td><td>三層 quality loop（Default-FAIL contract + Fresh-Context Evaluator + Agent-Maintained Handoff）；<code>.claude/agents/evaluator.md</code> + 5 個 hooks</td><td><strong>幾乎一對一可參照</strong>：對應 InterSubMod 的 Hard Gate + evaluator + PROGRESS</td></tr>
<tr><td>2</td><td><code>HKUDS/OpenHarness</code></td><td>10 子系統 (engine/tools/skills/plugins/permissions/hooks/commands/MCP/memory/multi-agent)，43 tools + 54 commands</td><td><strong>結構完全同構</strong>：對應 .claude/skills/* + MEMORY.md 索引</td></tr>
<tr><td>3</td><td><code>code-yeongyu/oh-my-openagent</code></td><td>Sisyphus 主 orchestrator + Hephaestus/Prometheus/Oracle/Librarian/Explore 角色分工，54+ lifecycle hooks，hash-anchored edit 工具</td><td>Team Mode 1+8 並行 → 對應 parallel-benchmark subagent</td></tr>
<tr><td>4</td><td><code>BulloRosso/etienne</code></td><td>Triple-P adaptive memory (Picker/Packer/Ponderer loop)，prompt injection detection gateway</td><td>對 Hard Gate 矩陣 + injection 防禦有參照</td></tr>
<tr><td>5</td><td><code>microsoft/skills</code></td><td>跨 Claude Code / Cursor / Copilot 統一的 skill + MCP + AGENTS.md 規範</td><td>對應 InterSubMod skill governance 跨 agent harness</td></tr>
<tr><td>6</td><td><code>harness/harness-skills</code></td><td>SKILL.md + YAML frontmatter + MCP server 依賴宣告</td><td><strong>命名同 SKILL.md，完全同構</strong></td></tr>
<tr><td>7</td><td><code>lastmile-ai/mcp-agent</code></td><td>MCP 之上 workflow pattern + durable execution（persistent state）</td><td>對 evidence_ledger + cross-session continuity 有參照</td></tr>
<tr><td>8</td><td><code>ai-boost/awesome-harness-engineering</code></td><td>完整 ecosystem 索引（memory 三層、AgentMiddleware 六 hook、LLMLingua 20x 壓縮、codebase-memory-mcp 120× token 壓縮）</td><td>探索 P4+ 改進方向的目錄</td></tr>
</tbody></table>

<h3>🔧 業界主流 agent harness 對比矩陣</h3>
<table>
<thead><tr><th>工具</th><th>Memory</th><th>Sub-agent</th><th>Hooks</th><th>MCP</th><th>自治程度</th></tr></thead>
<tbody>
<tr><td><strong>Claude Code</strong></td><td>CLAUDE.md + auto memory</td><td>Agent Teams 並行</td><td>Pre/Post ToolUse + UserPromptSubmit</td><td>一等公民</td><td>中-高</td></tr>
<tr><td>Codex CLI</td><td>無持久</td><td>無原生</td><td>無原生 hook</td><td>支援</td><td>中（PR-loop）</td></tr>
<tr><td>Cursor 2.0</td><td>規則檔 + agent tabs</td><td>8 並行 (git worktree 隔離)</td><td>無</td><td>部分</td><td>IDE 內中度</td></tr>
<tr><td>Cline</td><td>Memory Bank 五檔</td><td>無</td><td>無</td><td>支援</td><td>Plan/Act 模式</td></tr>
<tr><td>Devin</td><td>雲端 sandbox</td><td>無暴露</td><td>黑盒</td><td>內建</td><td>最高（fully autonomous）</td></tr>
<tr><td>Aider</td><td>Repo map (PageRank)</td><td>無</td><td>無</td><td>部分</td><td>低</td></tr>
<tr><td>OpenCode/oh-my-openagent</td><td>AGENTS.md hierarchical</td><td>Team Mode 1+8 並行</td><td><strong>54-61 lifecycle hooks</strong></td><td>內建</td><td>高</td></tr>
<tr><td><strong>InterSubMod</strong></td><td>MEMORY.md (37 feedback) + auto memory + Concluded 區</td><td>14 active + Agent Teams</td><td><strong>22 hooks 跨 6 events</strong></td><td>knowledge + context7</td><td><strong>中-高</strong></td></tr>
</tbody></table>

<p><strong>對齊結論</strong>: InterSubMod 屬「agent orchestrator」類（同 Claude Code / Devin / Pi），結構參照對象為 <code>cwc-long-running-agents</code> + <code>OpenHarness</code>。</p>

<h3>📚 10 條跨來源綜合 takeaway（Researcher agent 整合）</h3>
<ol>
<li><strong>Harness ≠ prompt</strong>：所有強行為（permission deny、sandboxing、subagent isolation、Hard Gate）必須是 architectural 而非 prompt 層級。Prompt 可被 jailbreak，architecture 不能。</li>
<li><strong>Repository is single source of truth</strong>：「agent 看不到 = 不存在」。docs/、AGENTS.md、CLAUDE.md 不在 repo 的 knowledge 對 agent 等於零。InterSubMod 已落實，但需 cold-start test 驗證完整性。</li>
<li><strong>Three context techniques 並用</strong>：compaction + structured note-taking + sub-agent，缺一即 context anxiety。</li>
<li><strong>Fresh-context evaluator 是 long-task 救星</strong>：agent 自評永遠樂觀。需 separate agent，無 Write/Edit 權限，獨立 context，PASS/NEEDS_WORK 二元判定。<strong>InterSubMod 缺此角色</strong>。</li>
<li><strong>Default-FAIL contract + evidence-gated write</strong>：所有 feature 預設 passes: false，PreToolUse hook 阻擋未讀 evidence 的 write。對應 InterSubMod evidence_ledger 但需要 hook 強制。</li>
<li><strong>Prompt caching 是必要不是選項</strong>：3+ step agent 不用 prompt caching = 燒錢。InterSubMod 已 94.6% hit rate（業界 top）✅</li>
<li><strong>Progressive disclosure 優於 monolithic context</strong>：InterSubMod 41 skills 已對齊；rules/ path-scoped 已 P3 達成 ✅</li>
<li><strong>Harness 隨 model 升級應簡化</strong>：每個 harness 元件編碼一個 model 做不到的假設，這些假設值得 stress test。InterSubMod CLAUDE.md 中 4.6 遺留 scaffolding 已標記可移除，需執行。</li>
<li><strong>Multi-agent role separation 是 trade-off</strong>：cwc（2-agent）vs three-agent harness 顯示「角色越多 → context isolation 越強但 handoff artifact 要求越高」。InterSubMod 預設單回合 + 明示觸發 subagent 是正確的 Opus 4.7 策略。</li>
<li><strong>Observability 是 harness 內部組件</strong>：watch 命令即時監控 PROGRESS + git log + screenshots + evidence reads。InterSubMod 缺 evidence_ledger 即時可視化 dashboard。</li>
</ol>

<h3>🎯 10 條 InterSubMod 立即落地建議（Researcher vs Architect cross-validate）</h3>

<table>
<thead><tr><th>優先</th><th>項目</th><th>Researcher 建議</th><th>Architect 建議</th><th>Cross-validate?</th></tr></thead>
<tbody>
<tr><td>P1</td><td>Fresh-Context Evaluator subagent</td><td>✅ 新增 .claude/agents/evaluator.md (cwc pattern)</td><td>✅ E3 建 evaluator agent</td><td><strong>✅ 強共識</strong></td></tr>
<tr><td>P1</td><td>Default-FAIL evidence gate hook</td><td>✅ PreToolUse 攔截寫入 validated/，要求先 Read evidence_ledger</td><td>—</td><td>新發現</td></tr>
<tr><td>P1</td><td>Rules path-scoped frontmatter</td><td>✅ 4 個 rules 加 paths frontmatter，省 ~5KB</td><td>—（已 P3 d64c0e7 commit 完成）</td><td><strong>已完成 ✅</strong></td></tr>
<tr><td>P1</td><td>Prompt caching 落點檢核</td><td>✅ system + skills + 大 doc 加 cache_control 標籤</td><td>—（已驗證 94.6% hit rate）</td><td><strong>已最優 ✅</strong></td></tr>
<tr><td>P2</td><td>Watch dashboard (L11)</td><td>✅ scripts/watch_session.sh tail PROGRESS+log+evidence</td><td>—</td><td>新發現</td></tr>
<tr><td>P2</td><td>Cold-start test 5 問題</td><td>✅ AGENTS.md 加 5 問題；SessionStart hook 驗證可答</td><td>—</td><td>新發現</td></tr>
<tr><td>P2</td><td>CURRENT_FOCUS.md SessionStart hook</td><td>✅ Phase 5 規劃落實</td><td>✅ E4（commit ee648fb 已完成）</td><td><strong>已完成 ✅</strong></td></tr>
<tr><td>P2</td><td>5 入口降到 3</td><td>—</td><td>✅ E5 合併 research_direction + manual 入 CURRENT_FOCUS</td><td>Architect 獨特建議</td></tr>
<tr><td>P3</td><td>Opus 4.6 遺留 scaffolding 清理</td><td>✅ 掃 41 skills 移除「double-check」「每 N 步 status」</td><td>—</td><td>新發現</td></tr>
<tr><td>P3</td><td>Context reset vs compaction 試驗</td><td>✅ cycle 末用 fresh evaluator 而非 /compact 比較產出</td><td>—</td><td>新發現</td></tr>
<tr><td>P3</td><td>SKILL.md 依賴宣告（mcp_required / phase / dependencies）</td><td>✅ 對齊 harness/harness-skills 規範</td><td>✅ E8 hook 重疊重評 + 部分對齊</td><td>部分共識</td></tr>
<tr><td>P3</td><td>researcher_claim_evidence hook</td><td>—</td><td>✅ E7 偵測 L3 推測語言時提示 §2 等級</td><td>Architect 獨特建議</td></tr>
<tr><td>P3</td><td>memory_recall_logger</td><td>—</td><td>✅ E9 量化每 skill / memory 引用率</td><td>Architect 獨特建議</td></tr>
<tr><td>P4</td><td>Injection guard</td><td>✅ 結合 etienne pattern 設計</td><td>✅ E10 external_input_sanitizer.sh</td><td><strong>✅ 強共識</strong></td></tr>
</tbody></table>

<h3>🔥 Top 3 高 ROI 立即可做（兩 agent 強共識）</h3>
<ol>
<li><strong>🏗️ Fresh-Context Evaluator subagent</strong>（兩 agent 強推 + cwc-long-running-agents pattern + Walking Labs L09）
  <ul><li>新增 <code>.claude/agents/evaluator.md</code>，無 Write/Edit/Bash，僅 Read+Grep</li>
  <li>用於 cycle 末 verify P5/P6 產出（取代或補強自評）</li>
  <li>PASS / NEEDS_WORK 二元輸出</li>
  <li>工時 2 hr</li></ul>
</li>
<li><strong>🚨 Default-FAIL evidence gate hook</strong>（cwc pattern + Anthropic harness reliability）
  <ul><li>新 <code>.claude/hooks/verify-gate.sh</code>: PreToolUse 攔截寫入 reports/validated/**</li>
  <li>要求先 Read evidence_ledger.jsonl 才放行 Write/Edit</li>
  <li>用 .claude/.evidence-reads 計數檔監控</li>
  <li>工時 1 hr</li></ul>
</li>
<li><strong>📺 Watch session dashboard</strong>（Walking Labs L11 + observability）
  <ul><li>新 <code>scripts/watch_session.sh</code>: tail PROGRESS + git log + cache_telemetry + 最新 figures</li>
  <li>用戶啟動長 cycle 時開獨立 terminal 監控</li>
  <li>工時 45 min</li></ul>
</li>
</ol>

<h3>📊 InterSubMod 累計 audit 完成度</h3>
<p>跨 architect + researcher cross-validate，發現 InterSubMod 已執行 P0-P3 對齊度：</p>
<ul>
<li>✅ <strong>SessionStart hook</strong> (Walking Labs L06) — commit ee648fb</li>
<li>✅ <strong>Silent failure → fail-loud log</strong> (L11 + Anthropic) — commit d5db8dc</li>
<li>✅ <strong>compact_test.sh E2E</strong> (L09 + L10) — commit d64c0e7</li>
<li>✅ <strong>subagent_completion_logger</strong> (L11 observability) — commit d64c0e7</li>
<li>✅ <strong>Rules path-scoped</strong> (progressive disclosure) — commit d64c0e7</li>
<li>✅ <strong>Cache hit 94.6%</strong>（Anthropic 90% 目標超越）— P2 已 cache_telemetry 驗證</li>
<li>✅ <strong>Memory system</strong>（37 feedback + Concluded + spaced recall）— pre-existing best-in-class</li>
<li>✅ <strong>Hard Gate 結構性對抗 "declare victory too early"</strong> (L07/L09) — pre-existing</li>
<li>✅ <strong>Progressive simplification / scaffolding 清理規範</strong>（opus47-behavior.md 已列）— 規範存在但執行待</li>
</ul>

<p><strong>仍未做（兩 agent 共識）</strong>：</p>
<ul>
<li>🟡 Fresh-Context Evaluator subagent（最高 ROI）</li>
<li>🟡 Default-FAIL evidence gate hook</li>
<li>🟡 Watch dashboard</li>
<li>🟡 5 入口降到 3</li>
<li>🟡 Injection guard / researcher_claim hook / memory_recall_logger</li>
<li>🟡 SKILL.md 依賴宣告統一</li>
<li>🟡 Cold-start test 5 問題</li>
<li>🟡 Opus 4.6 遺留 scaffolding 實際清理</li>
</ul>
"""

amended = existing.replace("<hr>\n<p class=\"meta\">Generated:",
                            RESEARCHER_SECTIONS + "<hr>\n<p class=\"meta\">Generated:")
amended = amended.replace(
    "Researcher agent still running — will amend in subsequent commit if findings change verdict.",
    f"Researcher agent ✅ amended (355 sec, 12 一手來源 + 4 WebSearch). Audit complete."
)

(OUTDIR / "p4_industry_deep_audit.html").write_text(amended, encoding="utf-8")
print(f"Final HTML: {len(amended)} chars")
