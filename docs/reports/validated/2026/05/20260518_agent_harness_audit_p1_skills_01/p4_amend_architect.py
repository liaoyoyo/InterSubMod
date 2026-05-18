#!/usr/bin/env python3
"""Amend P4 industry deep audit HTML with architect agent results."""
from pathlib import Path
from datetime import datetime

OUTDIR = Path("/tmp/skill_audit_20260518")
existing = (OUTDIR / "p4_industry_deep_audit.html").read_text(encoding="utf-8")

# Architect agent A-E sections to append
ARCHITECT_SECTIONS = """
<h2>🤖 Architect Agent 補強分析（107 sec runtime）</h2>
<p class="meta">獨立 architect agent 分析結果，與本報告主體 cross-validate。Agent TL;DR: 「InterSubMod harness 在『研究記錄持久性』與『方法論結構化』遠超業界平均，但在『observability / context reset / 評估者分離』三個 Anthropic 強推維度仍有 gap」</p>

<h3>§A 對齊優勢（業界平均之上）— 7 項</h3>
<table>
<thead><tr><th>#</th><th>優勢</th><th>對應業界原則</th><th>InterSubMod 落地</th></tr></thead>
<tbody>
<tr><td>A1</td><td><strong>Repo as Record 徹底落實</strong></td><td>Walking Labs L3</td><td><code>research/&lt;topic&gt;/00_INDEX.md</code> Pre-reg + REFLECTION + evidence_ledger.jsonl + cycle state.json — 跨 41 sessions 持久</td></tr>
<tr><td>A2</td><td><strong>結構化交付 &gt; 對話</strong>（檔案優先）</td><td>Anthropic Harness Design</td><td>weekly-report → pptx/html/postmortem 全走 template；對話 ephemera 顯式分流</td></tr>
<tr><td>A3</td><td><strong>Hard Gate 明確分級</strong></td><td>Walking Labs L7（過早宣告完成）</td><td>3 條 Hard Gate（刪檔 / C++ commit / NO-GO）+ pre-commit hook 強制 exit 2 — 業界 anti-pattern「declare victory too early」有結構性對抗</td></tr>
<tr><td>A4</td><td><strong>Productive Failure 框架</strong></td><td>學術（Kapur 2008-2016）</td><td>§8.3.1 Reopen Threshold C1/C2/C3 — 業界 harness 罕見正式化</td></tr>
<tr><td>A5</td><td><strong>NEGATIVE 結論索引化</strong></td><td>Replication Crisis 教訓</td><td>20+ Concluded 區條目 + spaced recall 30d/90d — 防重複調查</td></tr>
<tr><td>A6</td><td><strong>Cache hit rate 94.6%</strong></td><td>Anthropic context economics</td><td>條件式 rules + skill description USE WHEN — 進入 progressive disclosure 設計</td></tr>
<tr><td>A7</td><td><strong>5 入口分層</strong></td><td>Cline Memory Bank</td><td>governance (AGENTS.md) / agent-specific (CLAUDE.md) / live state (CURRENT_FOCUS) / manual / queue — 比業界 single CLAUDE.md 更分離</td></tr>
</tbody></table>

<h3>§B 對齊 Gap（業界推薦但缺）— 7 項</h3>
<table>
<thead><tr><th>#</th><th>Gap</th><th>來源</th><th>證據</th></tr></thead>
<tbody>
<tr><td>B1</td><td><strong>Generator-Evaluator 分離不徹底</strong></td><td>Anthropic Harness Design</td><td>researcher / reviewer 命名與 plugin 重疊；scientist-reviewer 提案還在 A5-2 P4；無正式 evaluator agent 對 generator 輸出做獨立評分</td></tr>
<tr><td>B2</td><td><strong>Observability 不足</strong></td><td>Walking Labs L10</td><td>P5 audit A5-3 明確：「缺 cache_hit / cost / completion telemetry — 不知道哪些 agent 真有用」+ H4-5 缺 hook latency</td></tr>
<tr><td>B3</td><td><strong>Context Reset &gt; Compaction 未落實</strong></td><td>Anthropic harness-design 原文</td><td>仍仰賴 <code>/compact</code> + 5 入口；缺「完整 clean state」protocol（Walking Labs L11）</td></tr>
<tr><td>B4</td><td><strong>E2E Harness Testing 缺</strong></td><td>Walking Labs L9</td><td>無 <code>compact_test.sh</code>（已 P3 補）/ 無 <code>hook_concurrency_test.sh</code>（H4-4 P3）— 22 hooks 並行假設未驗證</td></tr>
<tr><td>B5</td><td><strong>Initialization 程序非標準化</strong></td><td>Walking Labs L6</td><td>SessionStart hook 已 P1 補 ✅；目前 AI 主動 Read CURRENT_FOCUS 自律</td></tr>
<tr><td>B6</td><td><strong>Silent failure 28% (已 P0 修)</strong></td><td>Anthropic harness reliability</td><td>P0 H4-3 已 commit d5db8dc — 16/21 hooks 改 fail-loud log ✅</td></tr>
<tr><td>B7</td><td><strong>Prompt Injection Guard 缺</strong></td><td>Anthropic security</td><td>P4 backlog <code>external_input_sanitizer.sh</code> 仍未實作</td></tr>
</tbody></table>

<h3>§C 過度工程 / Over-Engineering Risk — 6 項</h3>
<table>
<thead><tr><th>#</th><th>風險項</th><th>業界對比</th><th>緩解建議</th></tr></thead>
<tbody>
<tr><td>C1</td><td><strong>44 skills + 16 plugin = 60 觸發空間</strong></td><td>Voyager skill lib 通常 &lt; 30</td><td>已部分緩解（cache 94.6%）；但 cognitive load 對人類 maintainer 仍高。建議 quarterly 跑 memory_recall_logger 量化每 skill 實際觸發率，&lt;5%/月 → 候選 archive</td></tr>
<tr><td>C2</td><td><strong>§11 級聯圖 12 步</strong></td><td>Walking Labs L4「巨大 instruction file 失敗」反例</td><td>§0.5 已加最小子集 mitigation，但 12 步 cascade 仍是長尾觸發點；建議用「entrance test」自查</td></tr>
<tr><td>C3</td><td><strong>Meta-skill 層 (/scientific-rigor)</strong></td><td>README §8 自承「非業界共識，本專案創新」</td><td>ROI 邊際 — 整合 7 個底層 skill 但用戶可能直接呼叫子 skill；建議追蹤被觸發次數</td></tr>
<tr><td>C4</td><td><strong>PostToolUse Edit|Write 8 個 hooks</strong></td><td>Anthropic 「漸進式簡化」</td><td>H4-2 P3 已評估無重疊；但 architect 補：<code>terminology_guard</code> ∩ <code>md_link_check</code> ∩ <code>kb_sot_guard</code> 三者語意接近，再 review</td></tr>
<tr><td>C5</td><td><strong>17 段 structured-tech-report</strong></td><td>Diátaxis 4 種文件類型即可覆蓋</td><td>用 13 段 + 4 段 optional 分流，降強制段數</td></tr>
<tr><td>C6</td><td><strong>5 入口載入</strong></td><td>Anthropic single source of truth</td><td>CURRENT_FOCUS hook 注入後 §9 列表可刪 → 5 入口降到 3</td></tr>
</tbody></table>

<h3>§D 用戶風格適配建議 — 6 項</h3>
<p>用戶偏好：<em>業界嚴謹科學工程語言 / 全自動+並行 / 實證 L1 驗證 / 小規模 pilot &lt;2hr / 結構化決策表</em></p>
<table>
<thead><tr><th>#</th><th>建議強化</th><th>對應用戶偏好</th><th>落地位置</th></tr></thead>
<tbody>
<tr><td>D1</td><td><strong>「實證驗證 L1」hook 化</strong> — researcher claim 必升級實測才能定論</td><td>實證 L1 (memory: feedback_researcher_claim_needs_empirical_verification)</td><td>新增 hook <code>researcher_claim_evidence_check.sh</code>（PostToolUse，偵測「我推測 / 我認為 / probably」時提示用 §2 等級標記）</td></tr>
<tr><td>D2</td><td><strong>&lt;2hr pilot 強制 SMART</strong></td><td>小規模 pilot 偏好</td><td><code>/problem-framing-ideation</code> §1 已 P1 加 5W2H；建議加「pilot &lt; 2hr 否則 split」rule</td></tr>
<tr><td>D3</td><td><strong>結構化決策表 default 樣式</strong></td><td>偏好結構化</td><td>TL;DR + decision table + (影響, 信心) 標註 → 已在 AGENTS.md §15.1 落地，但需 hook 強化驗證每個 Tier 3+ 回應首句格式</td></tr>
<tr><td>D4</td><td><strong>全自動 + 並行明示觸發語自動偵測</strong></td><td>全自動偏好</td><td>建議加 <code>parallel_keyword_detector.sh</code> UserPromptSubmit hook 自動建議 spawn</td></tr>
<tr><td>D5</td><td><strong>業界嚴謹語言 = effect size + CI 必標</strong></td><td>業界嚴謹 (memory: feedback_outside_claim_must_query_kb)</td><td><code>/scientific-rigor §2.1 checklist</code> 已落地；建議補 PostToolUse 偵測「significant / improvement / better」未配 Cohen 等級時 stderr warn</td></tr>
<tr><td>D6</td><td><strong>KB 查詢義務自動化</strong></td><td>外部 claim 必先查 KB</td><td><code>kb_schema_guard.sh</code> 已存在；建議擴展為「外部工具名提及 → 強制 mcp__knowledge_search」hard gate</td></tr>
</tbody></table>

<h3>§E P4+ Fix Recommendations — 10 條（Architect agent 完整清單）</h3>
<table>
<thead><tr><th>優先</th><th>項目</th><th>影響 / 信心</th><th>工時</th><th>對齊框架</th></tr></thead>
<tbody>
<tr style="background:#dafbe1;"><td>已 P0</td><td><strong>E1</strong>: H4-3 silent failure → ✅ commit d5db8dc 已修</td><td>高/高</td><td>45m (DONE)</td><td>Walking Labs L10 + Anthropic fail-loud</td></tr>
<tr style="background:#dafbe1;"><td>已 P3</td><td><strong>E2</strong>: compact_test.sh ✅ commit d64c0e7 已建</td><td>高/中-高</td><td>1.5h (DONE)</td><td>Walking Labs L9 E2E testing</td></tr>
<tr style="background:#fff8c5;"><td>P1</td><td><strong>E3</strong>: 建 evaluator agent（與 researcher/reviewer 解耦），明確 generator-evaluator 分工</td><td>高/高</td><td>2h</td><td>Anthropic Harness Design 核心</td></tr>
<tr style="background:#dafbe1;"><td>已 ee648fb</td><td><strong>E4</strong>: SessionStart hook ✅ commit ee648fb 已建</td><td>中-高/高</td><td>1h (DONE)</td><td>Walking Labs L6</td></tr>
<tr style="background:#fff8c5;"><td>P1</td><td><strong>E5</strong>: 5 入口降到 3 — 合併 research_direction.md (9d stale) + manual.md 入 CURRENT_FOCUS</td><td>中/高</td><td>1h</td><td>Anthropic SoT</td></tr>
<tr style="background:#dafbe1;"><td>已 P3</td><td><strong>E6</strong>: subagent_completion_logger ✅ commit d64c0e7 已建</td><td>中/高</td><td>1.5h (DONE)</td><td>Walking Labs L10</td></tr>
<tr style="background:#fff8c5;"><td>P2</td><td><strong>E7</strong>: researcher_claim_evidence_check.sh — 偵測 L3 推測語言時提示 §2 等級標記</td><td>中-高/中</td><td>45m</td><td>用戶風格 D1 + Replication Crisis</td></tr>
<tr style="background:#fff8c5;"><td>P2</td><td><strong>E8</strong>: deprecate 候選 — terminology_guard ∩ md_link_check ∩ kb_sot_guard 三選一保留</td><td>中/中</td><td>45m</td><td>C4 過度工程緩解</td></tr>
<tr style="background:#fff8c5;"><td>P3</td><td><strong>E9</strong>: memory_recall_logger 量化 skill / memory 引用率，月度 &lt;5% 候選 archive</td><td>中/中</td><td>1h</td><td>C1 + Walking Labs context drift</td></tr>
<tr style="background:#fff8c5;"><td>P3</td><td><strong>E10</strong>: prompt injection guard external_input_sanitizer.sh</td><td>中/中</td><td>1h</td><td>B7 Anthropic security</td></tr>
</tbody></table>

<p><strong>Architect 補強的關鍵發現</strong>:</p>
<ul>
<li>✅ E1 / E2 / E4 / E6 已在 P0-P3 fix 完成（commits ee648fb / d5db8dc / d64c0e7）</li>
<li>🟡 剩 6 個未做：E3 evaluator agent / E5 5 入口降 3 / E7 researcher claim hook / E8 hook 重疊重評 / E9 recall logger / E10 injection guard</li>
<li>📊 對齊度差距主要在 <strong>觀察性 (observability)</strong> 與 <strong>generator-evaluator 分工</strong> — 與本報告主體一致</li>
</ul>

<h3>需進一步研究項（Architect 標記）</h3>
<ul>
<li><strong>Walking Labs L8 "Feature lists primitives"</strong>: 未取得完整 lecture 內容，無法判定 InterSubMod skill 設計是否符合 primitives 思維（vs feature lists）</li>
<li><strong>Anthropic 美學可分級化</strong>: 已知原則但 InterSubMod 對應落地點不明 — 可能 implicit 在 weekly-report 12 條 design principles，但未顯式宣告「美學分級」標準</li>
</ul>

"""

# Insert before final hr
amended = existing.replace("<hr>\n<p class=\"meta\">Generated:",
                            ARCHITECT_SECTIONS + "<hr>\n<p class=\"meta\">Generated:")

# Also amend with researcher agent placeholder (still running)
amended = amended.replace(
    "2 background agents (researcher + architect) still running for amendments",
    "Architect agent ✅ amended (107 sec). Researcher agent still running — will amend in subsequent commit if findings change verdict."
)

(OUTDIR / "p4_industry_deep_audit.html").write_text(amended, encoding="utf-8")
print(f"Amended HTML: {len(amended)} chars (was 19744)")
print(f"Net add: {len(amended) - 19744} chars")
