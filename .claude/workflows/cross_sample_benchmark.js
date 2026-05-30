// cross_sample_benchmark.js — InterSubMod 跨樣本 benchmark Dynamic Workflow 模板
//
// 用途：把 parallel-benchmark agent 處理的「N 個無資料依賴的同質樣本 benchmark」
//       改用 Dynamic Workflow 執行 — plan 從 context 移到 script、resumable、
//       一樣本掛掉不全重跑、中間結果不塞爆 context。
//
// 何時用（對齊 CLAUDE.md §8 Dynamic Workflow 路由規則）：
//   ✅ ≥3 樣本、同一 benchmark、無資料依賴、終態明確（一致性匯總表）
//   ❌ 含 Hard Gate（C++ commit/刪檔/NO-GO）→ 不可用 workflow（subagent acceptEdits 繞過 exit-2 hook）
//   ❌ 需中間結果決定下一步 → 用 parallel-benchmark agent 或主回合
//
// 用法：
//   Workflow({ scriptPath: ".../.claude/workflows/cross_sample_benchmark.js",
//              args: { samples: ["HCC1395_5kHz","COLO829","H2009"], benchmark: "<描述>" } })
//   或在對話中含 "workflow" keyword 觸發後請 Claude 載入本 script。
//
// ⚠ 跑前：(1) /model 確認模型；(2) 把 benchmark 需要的 Bash/python 指令先加進 allowlist
//        （否則 mid-run 權限提問打斷背景 run）；(3) 先用 1-2 樣本 scoped 校準 token。

export const meta = {
  name: 'cross-sample-benchmark',
  description: 'Run the same benchmark across N independent samples in parallel, then synthesize a cross-sample consistency table',
  phases: [
    { title: 'Benchmark', detail: 'one isolated agent per sample (worktree)' },
    { title: 'Synthesize', detail: 'cross-sample canonical ordering + consistency stats' },
  ],
}

const samples = (args && args.samples) || ['HCC1395_5kHz', 'COLO829', 'H2009']
const benchmarkDesc = (args && args.benchmark) || '<描述要跑的 benchmark：metric / baseline / dataset path>'

const RESULT = {
  type: 'object',
  additionalProperties: false,
  required: ['sample', 'metrics', 'verdict'],
  properties: {
    sample: { type: 'string' },
    metrics: { type: 'object', additionalProperties: true, description: 'F1 / AUC / TP / FP / delta vs baseline 等' },
    notes: { type: 'string' },
    verdict: { type: 'string', description: 'POSITIVE / NEGATIVE / MARGINAL + 一句理由' },
  },
}

phase('Benchmark')

// 每樣本一個隔離 agent（worktree 避免平行寫 figures/ 競爭）
const results = await parallel(
  samples.map((s) => () =>
    agent(
      `針對樣本 ${s} 執行以下 benchmark：${benchmarkDesc}\n`
      + `要求：(1) 讀對應 dataset；(2) 對 baseline 算 delta；(3) 標 effect size（Cohen's d / CI）；`
      + `(4) 回 metrics + verdict。不要宣告跨樣本結論（那是 synthesize 階段的事）。`
      + `對齊 /scientific-rigor §2：每個數字標來源檔案路徑。`,
      { label: `bench:${s}`, phase: 'Benchmark', schema: RESULT, isolation: 'worktree' }
    )
  )
)

const valid = results.filter(Boolean)

phase('Synthesize')

const summary = await agent(
  `以下是 ${valid.length}/${samples.length} 個樣本的 benchmark 結果（JSON）：\n${JSON.stringify(valid, null, 2)}\n\n`
  + `產出跨樣本一致性報告（對齊 /multi-sample-consistency）：\n`
  + `1. Canonical 排序表（樣本 × metric）\n`
  + `2. 方向一致性統計（幾/N 樣本同向 POSITIVE）\n`
  + `3. Confidence uplift 判定（≥5/7 同向才可升 ⭐3+）\n`
  + `4. 標出 outlier 樣本 + 可能原因\n`
  + `⚠ 不要在此 workflow 內升 tier 或寫 evidence_ledger（那是 Hard Gate / 主 agent 職責）。只產分析供主 agent 審。`,
  { label: 'synthesize', phase: 'Synthesize' }
)

return { summary, per_sample: valid, missing: samples.length - valid.length }
