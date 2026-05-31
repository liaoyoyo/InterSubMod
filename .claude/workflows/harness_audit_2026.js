export const meta = {
  name: 'harness-audit-2026',
  description: '業界 AI 工作流研究 + 全 harness 盤點 + restraint 對抗驗證 gap 分析 + 自我稽核儀表板設計',
  phases: [
    { title: 'Research+Inventory' },
    { title: 'Restraint-Verify' },
    { title: 'Synthesize' },
  ],
}

// ── 內部 harness 摘要（餵給每個 research agent 做 overlap 自評）──
const HARNESS = `InterSubMod Claude Code harness（單人研究者，C++/Python 癌症長讀長定序研究）現況（2026-05-31 實測 ground truth）：
- 45 skills：7-Phase Waterfall (P0 cycle-init → P1 research-loop → P2 check-staleness → P3 feature-layered-observation → P4 multi-sample-consistency → P5 run-evaluator → P6 conclude-research) + 9 元方法論 (confirmation-protocol/known-pitfalls/cycle-state/research-context-loader/scientific-rigor/pre-decision-audit/problem-framing-ideation/provenance-tier-audit/fast-learning-coach) + 報告類 + narrative-frame catalog 50+ framework + pipeline-manifest(reproducibility DAG)
- 18 agents：architect developer optimizer tester researcher reviewer evaluator methodology-reviewer(方法層) parallel-benchmark parallel-analysis headless-research narrative-organizer(haiku) research-orchestrator(haiku) literature-reviewer paper-miner release reproducibility-audit(fresh-context read-only) security-reviewer(fresh-context read-only)
- 38 hooks 跨 7 events；5 個 exit-2 Hard Gate hooks：pre_commit_compile_check(C++ commit 必編譯) / kb_schema_check / pipeline_block_check(磁碟) / no_binary_commit / kb_sot_guard(F1 SoT 數字保護)；另 search_scope_guard(阻 grep -r . / find . 無 maxdepth) / verify_gate(Default-FAIL evidence gate) / task_type_advisor(6 類) / narrative_frame_advisor / precompact_autosave
- settings：model=claude-opus-4-8[1m] / effortLevel=xhigh / compactPrompt 已設 / sandbox=docker(debian:bookworm) / defaultMode 已移除(繼承 user auto) / 6 個 plugin 全 DISABLED / MCP=knowledge+biorxiv+ensembl(全 deferred-load 非 always-load) / remote env 已設
- 1 Dynamic Workflow (cross_sample_benchmark.js)
- 證據鏈：evidence_ledger.jsonl(L1-L5 provenance + ⭐1-5 tier) + hypothesis_queue.json + state machine(active.json/cycles/，由 skill 寫不可手改) + memory(file-based，frontmatter type=user/feedback/project/reference)
- 既有 Hard Gate 規則(§1)：刪檔永遠🔴 / C++ commit 必編譯 / 研究 NO-GO 判定需用戶確認 / evidence_ledger 覆寫🔴
- Cynefin domain + DACI/DECIDE 決策框架已落地(confirmation-protocol)；productive-failure reopen C1/C2/C3 門檻已落地(scientific-rigor §8.3.1)`

const RESTRAINT = `關鍵稽核 lens（務必套用）：此 harness 相對 2026 業界基線【已成熟】。上一輪 21-agent / 1.6M-token 研究中，12 條 high-fit 建議對抗驗證後只存活 4 條（全降 medium）。**預設 RESTRAINT（不裝）而非 adoption（裝）**。
對每個你找到的 pattern，必須明確指出：(a) 此 harness 是否【已涵蓋】（引用具體 skill/agent/hook 名）；(b) headline/business-case 是真實還是炒作（如 cost-watch 在 flat-rate 訂閱下 cost=0；deterministic-eval 在無 machine-checkable field 時塌成 brittle keyword match）；(c) 採納是否會【正中已知失敗類】（/tmp 800GB disk-full / rebuild 覆寫手動 commit / 繞過 exit-2 Hard Gate / always-load MCP context 膨脹）。
只有在【L1/L2 證據 + 具體缺口】時才建議 ADD。優先「已涵蓋」/「修 latent bug」/「補 trivial 缺口」勝過「裝新工具」。每條給 evidence tier（L1=讀過實檔/官方一手 / L2=官方文件 / L3=部落格/推論）。`

// ── 6 路業界研究（對應用戶 3 子問題：理解 AI 產生 / 驗證 AI 正確 / 評估與計劃）──
const research = [
  { key: 'R1-orchestration', q: `研究 2025-2026【快速竄升 / 被推崇的 AI agentic 工作流與編排模式】。WebSearch 查："agentic workflow patterns 2026"、"Anthropic building effective agents"、"deep agents subagent orchestration"、"plan and execute agent"、"orchestrator worker evaluator optimizer"、"context engineering 2026"。
聚焦：哪些 orchestration pattern（prompt chaining / routing / orchestrator-worker / evaluator-optimizer / reflection / plan-execute / deep agents / dynamic workflow）正在成為主流，誰在用（Anthropic/OpenAI/LangChain/Cursor/Devin/社群），成熟度。對每個 pattern 標來源 URL + evidence tier。` },
  { key: 'R2-verify', q: `研究【如何確認 AI 工作的細節是正確、合理、可驗證】的頂尖方法（這是用戶核心問題之一）。WebSearch 查："LLM as a judge 2026"、"generator evaluator separation agents"、"adversarial verification AI agents"、"self-consistency verification LLM"、"AI verifier model"、"deterministic eval vs LLM judge"、"machine checkable verification agent output"、"agent eval harness"。
聚焦：generator/evaluator 分離、LLM-as-judge、對抗式驗證（refute-by-default）、self-consistency、verifier model、machine-checkable assertion、provenance/citation 驗證。哪些對【單人研究 + C++/科學 pipeline + 證據鏈】場景真有用。標來源 + tier。` },
  { key: 'R3-observability', q: `研究【如何理解 AI 的產生 / AI 可觀測性與推理透明度】（用戶核心問題之一）。WebSearch 查："LLM observability 2026"、"agent tracing langfuse langsmith braintrust"、"OpenTelemetry GenAI semantic conventions"、"reasoning transparency LLM"、"agent eval dashboard"、"interpretability for practitioners 2026"。
聚焦：tracing/observability 平台（Langfuse/LangSmith/Braintrust/OTel GenAI）、推理過程透明化、token/cost telemetry、eval dashboard。對【本地單機 + flat-rate 訂閱 + 已有 cache_telemetry/subagent_completion_logger hook】場景哪些真有 ROI、哪些是企業 SaaS 過度。標來源 + tier。` },
  { key: 'R4-planning', q: `研究【如何更好評估、計劃、確認任務】的頂尖方法（用戶核心問題之一）。WebSearch 查："spec driven development AI 2026"、"TDD for AI agents"、"task decomposition agent planning"、"pre-mortem decision making AI"、"plan mode claude code best practices"、"agent self-evaluation reflection"。
聚焦：spec-driven development、TDD-for-agents、計劃分解、pre-mortem、decision framework（Cynefin/DACI 本 harness 已有）、self-eval/reflection loop。哪些補足本 harness 的「pre-decision-audit / implementation-notes / run-evaluator 三層樓」與 7-Phase。標來源 + tier。` },
  { key: 'R5-agent-roles', q: `研究【AI agent 各角色佈局、權限劃分、least-privilege、subagent 設計】的頂尖實踐。WebSearch 查："multi-agent role separation best practices"、"subagent design claude code"、"least privilege AI agent tools"、"context isolation subagent"、"agent permission scoping 2026"、"when to use subagents vs single agent"。
聚焦：角色分離原則（generator vs verifier 必分離、read-only reviewer）、權限/工具範圍最小化、context isolation、subagent fan-out 何時值得、model-tier 分派（便宜模型做機械型）。對照本 harness 18-agent 佈局給改進建議。標來源 + tier。` },
  { key: 'R6-self-audit', q: `研究【AI harness 自我稽核、健康儀表板、config drift 偵測、prompt/skill 版本治理、agent CI/regression gate】的頂尖實踐。WebSearch 查："agent harness health dashboard"、"prompt versioning governance"、"config drift detection AI agents"、"agent regression testing CI"、"skill registry management"、"eval driven development dashboard 2026"。
聚焦：如何建一套「每次更新/任務改動都能確認與驗證哪些是何、哪些要改進」的固定觀察流程與儀表板（用戶最終核心需求 Part D）。標來源 + tier。` },
]

// ── 3 路內部盤點（讀實檔，產生權威內部 map）──
const inventory = [
  { key: 'I1-agent-fleet', q: `讀遍 InterSubMod/.claude/agents/ 下全部 18 個 .md（architect developer optimizer tester researcher reviewer evaluator methodology-reviewer parallel-benchmark parallel-analysis headless-research narrative-organizer research-orchestrator literature-reviewer paper-miner release reproducibility-audit security-reviewer）。
對每個 agent 萃取：role 一句話 / tools 權限範圍 / model(haiku|inherit|sonnet|opus) / USE-WHEN / SKIP-WHEN。產出 markdown 表。
然後分析：(a) role 重疊（哪些 agent 做類似工作）；(b) 權限劃分問題（過度或不足授權；generator 與 verifier 是否正確分離）；(c) 覆蓋缺口（對「理解 AI 產生 / 驗證正確 / 計劃評估」三軸是否缺角色）；(d) model-tier 是否恰當（機械型該降 haiku）。套用 restraint。` },
  { key: 'I2-skill-hook-map', q: `讀 InterSubMod/.claude/CLAUDE.md（§0 task-type gate / §3 skills 分類 / §4 hooks / §5 rules / §8 agent+Dynamic Workflow 路由 / §11 narrative-frame / §12 搜尋紀律）。可抽樣讀 3-4 個 SKILL.md（如 .claude/skills/run-evaluator/SKILL.md, scientific-rigor/SKILL.md, pre-decision-audit/SKILL.md）確認結構。
產出：(1) 7-Phase Waterfall 鏈 + 進入點；(2) 元方法論/報告/視覺化 skill 分類；(3) 38 hooks 按 7 events 分布 + 5 個 Hard Gate；(4) 找出：orphan skill（無進入點）、功能重疊、缺漏的「觀察/驗證固定流程」。套用 restraint。` },
  { key: 'I3-task-state', q: `讀以下確認【當前任務與細節】（用戶 Part C）：
- InterSubMod/docs/CURRENT_FOCUS.md（最新段落 = 當前主軸）
- InterSubMod/research/autoresearch/hypothesis_queue.json
- InterSubMod/research/autoresearch/evidence_ledger.jsonl（用 tail 讀最後 ~15 行，勿全讀）
- InterSubMod/state/active.json
- InterSubMod/docs/experiments/INDEX.md（讀開頭即可）
產出：(1) 當前 LIVE 主軸與 tier（LOH-phasing ⭐3 / ZAR1L-BRCA2 ASM ⭐3 / FP-filter ❌DEAD）；(2) pending / blocker 清單；(3) 每條 task 適合哪個 agent/skill/workflow + 需要什麼權限；(4) 是否有 task 與權限劃分不匹配。客觀陳述，勿臆造數字。` },
]

phase('Research+Inventory')
const gathered = await parallel([
  ...inventory.map(t => () => agent(`${t.q}\n\n[內部現況參考]\n${HARNESS}`, { label: t.key, phase: 'Research+Inventory' })),
  ...research.map(t => () => agent(`${t.q}\n\n${RESTRAINT}\n\n[待對照的內部現況]\n${HARNESS}`, { label: t.key, phase: 'Research+Inventory' })),
])

const invMaps = gathered.slice(0, inventory.length).filter(Boolean)
const researchMd = gathered.slice(inventory.length).filter(Boolean)

// ── 對抗式 restraint 驗證（generator/evaluator 分離；獨立第二意見）──
phase('Restraint-Verify')
const VERDICT_SCHEMA = {
  type: 'object',
  additionalProperties: false,
  properties: {
    verdicts: {
      type: 'array',
      items: {
        type: 'object',
        additionalProperties: false,
        properties: {
          rec: { type: 'string', description: '建議一句話' },
          area: { type: 'string', description: '所屬研究軸 R1-R6' },
          already_covered_by: { type: 'string', description: '已涵蓋它的 skill/agent/hook 名；若無填 none' },
          headline_real: { type: 'boolean', description: 'business-case/headline 是否經得起對抗' },
          known_failure_risk: { type: 'string', description: '若採納會否正中已知失敗類；若無填 none' },
          verdict: { type: 'string', enum: ['ALREADY_COVERED', 'ADD', 'MODIFY', 'SKIP_HYPE', 'SKIP_RISK'] },
          evidence_tier: { type: 'string', enum: ['L1', 'L2', 'L3'] },
          effort: { type: 'string', enum: ['trivial', 'small', 'medium', 'large'] },
          rationale: { type: 'string' },
        },
        required: ['rec', 'area', 'already_covered_by', 'headline_real', 'verdict', 'evidence_tier', 'effort', 'rationale'],
      },
    },
    summary: { type: 'string', description: '一段總結：survive ADD 幾條 / ALREADY_COVERED 幾條 / SKIP 幾條 + 最高 ROI 的 3 條' },
  },
  required: ['verdicts', 'summary'],
}

const verdicts = await agent(
  `你是獨立的對抗式 restraint 審查員（generator/evaluator 分離原則 — 你沒看過上面 agent 的自評，要用懷疑態度重審）。
${RESTRAINT}

[內部 harness 權威盤點]
${invMaps.join('\n\n=====\n\n')}

[6 路業界研究與初步建議]
${researchMd.join('\n\n=====\n\n')}

任務：從上述全部研究中【抽出每一條具體建議】，去重，逐條給 verdict。預設懷疑：能歸為 ALREADY_COVERED 就不要給 ADD；headline 站不住就 SKIP_HYPE；會撞已知失敗類就 SKIP_RISK。只有真正有缺口 + L1/L2 證據才 ADD/MODIFY。summary 要點名「真正值得做的最高 ROI 前 3 條」。`,
  { schema: VERDICT_SCHEMA, label: 'restraint-critic', phase: 'Restraint-Verify' }
)

// ── 儀表板 + 固定觀察流程設計（用戶 Part D 核心交付）──
phase('Synthesize')
const dashboard = await agent(
  `設計一套【harness 健康自我稽核儀表板 + 固定觀察/驗證流程】，讓單人研究者每次「更新或任務改動」都能快速確認「哪些是何、哪些要改進」。這是用戶最終核心需求。

[內部現況盤點]
${invMaps.join('\n\n=====\n\n')}

[已通過 restraint 驗證的建議]
${JSON.stringify(verdicts.verdicts.filter(v => v.verdict === 'ADD' || v.verdict === 'MODIFY'), null, 1)}

要求（務必具體、可落地、低維護成本、單人友善 — 不要企業級 SaaS）：
1. **儀表板要追蹤哪些維度**（建議分層：L0 一眼健康燈號 / L1 drift 與待辦 / L2 明細）。涵蓋：skill/agent/hook 計數 vs 文件宣稱（drift）、settings 關鍵值、active cycle/tier、evidence_ledger 最新、待確認 Hard Gate、最近 postmortem。
2. **固定觀察流程**：一個可重複的 checklist / 命令序列（最好能掛成既有 hook 或一個新 skill /harness-audit），列出每步驟讀什麼、比對什麼、紅旗條件。
3. **與既有資產整合**：明確指出該擴充既有什麼（research-dashboard skill / skill_registry_sync.sh / goal-landscape dashboard HTML / cache_telemetry / state/cycle-state）而非新造輪子。
4. **產出形式建議**：單一 HTML（仿既有 goal-landscape dashboard，localStorage 可編輯）還是 markdown 還是 skill？給出推薦 + 理由。
5. **ASCII 線框**：畫出儀表板 L0 健康燈號區的版面示意。

輸出 markdown。`,
  { label: 'dashboard-design', phase: 'Synthesize' }
)

return {
  inventory_maps: invMaps,
  research_findings: researchMd,
  restraint_verdicts: verdicts,
  dashboard_design: dashboard,
}
