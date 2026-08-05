export const meta = {
  name: 'loh-cnloh-external-evidence',
  description: '文獻+KB 證實：SEQC2 HCC1395 LOH/cnLOH 口徑、normal-germline 在腫瘤 cnLOH 區仍雜合的已知現象、cell-line 純度/subclonal LOH、甲基在 LOH/cnLOH 區的 haplotype 行為',
  phases: [
    { title: 'LitSearch', detail: '4 平行：KB-first 再 web' },
    { title: 'Synthesize', detail: '外部證據對本 pilot 假說的支持/反駁裁決' },
  ],
}
const KB = '/big8_disk/liaoyoyo2001/Knowledge'
const SCHEMA = {
  type: 'object', additionalProperties: false,
  properties: {
    summary: { type: 'string' },
    findings: { type: 'array', items: { type: 'object', additionalProperties: false,
      properties: {
        claim: { type: 'string' }, evidence: { type: 'string' }, source: { type: 'string' },
        relevance: { type: 'string' },
        tier: { type: 'string', enum: ['KB-verified','peer-reviewed','preprint','tool-doc','blog-weak'] },
      }, required: ['claim','evidence','source','relevance','tier'] } },
    key_takeaways: { type: 'array', items: { type: 'string' } },
    open_questions: { type: 'array', items: { type: 'string' } },
  }, required: ['summary','findings','key_takeaways','open_questions'],
}
const RULE = '紀律：先 grep 本地 KB(' + KB + ') 再 WebSearch/WebFetch(用 ToolSearch 載入)。每 claim 標 tier+來源。生醫主題優先 PubMed/期刊。不可編造 citation，找不到就說找不到。'
phase('LitSearch')
const dims = [
  { key: 'seqc2-loh-def', label: 'SEQC2 HCC1395 LOH/cnLOH 口徑',
    prompt: ['你是文獻調查員。主題：SEQC2 HCC1395 CNV/LOH benchmark (Masood et al. Genome Biology 2024) 的 LOH 定義與建立方法。',
      '要找：(1) SEQC2 LOH 是否=cnLOH(copy-neutral)？用什麼資料/技術界定(WGS callers + CytoScan/BeadChip/Bionano)？',
      '(2) LOH segment 邊界與信心如何定？是否有 subclonal/低信心 LOH？',
      '(3) SEQC2 用的是腫瘤 vs normal 配對如何界定 LOH？',
      '(4) HCC1395 已知的 LOH/cnLOH 範圍(覆蓋多少基因組)。',
      '先查 KB ' + KB + '/02_samples/hcc1395.md 與 papers/。', RULE, '回傳 SCHEMA。'].join('\n') },
  { key: 'normal-germline-loh', label: 'normal-germline 在腫瘤 LOH 區仍雜合',
    prompt: ['你是文獻調查員。主題：用「normal/germline 樣本」做 haplotype phasing，在「腫瘤有 cnLOH/LOH」的區域，normal germline 是否仍是平衡雜合(VAF~0.5)？這是本 pilot 的核心解釋。',
      '要找：(1) tumor-normal paired 分析中，cnLOH 是腫瘤體細胞事件、normal 不帶 → normal germline het 在該區仍雜合 — 這是否文獻共識？',
      '(2) 用 normal germline 給 tumor read 定相(haplotag)的標準做法(longphase/whatshap)在 LOH 區如何處理？',
      '(3) cnLOH 機制(mitotic recombination/gene conversion)為何 copy-neutral 但失雜合。', RULE, '回傳 SCHEMA。'].join('\n') },
  { key: 'cellline-purity-subclonal', label: 'HCC1395 cell-line 純度 + subclonal LOH',
    prompt: ['你是文獻調查員。主題：HCC1395 cell line 的純度與 subclonal/異質性，解釋「SEQC2 標 LOH 但 tumor read VAF~0.5(仍雜合)」。',
      '要找：(1) HCC1395 是純 cell line 嗎？有無 normal/間質污染？ONT BAM 純度。',
      '(2) cell line 仍可有 subclonal CNV/LOH(部分細胞才有) → bulk read VAF 介於 0.5 與 0/1。',
      '(3) SEQC2 LOH 用 WGS bulk consensus，subclonal LOH 可能被標為 LOH 但 read 層仍部分雜合。',
      '先查 KB ' + KB + '/02_samples/(hcc1395.md, subsample-purity.md)。', RULE, '回傳 SCHEMA。'].join('\n') },
  { key: 'methyl-loh-haplotype', label: '甲基在 LOH/cnLOH 區的 haplotype 行為',
    prompt: ['你是文獻調查員。主題：cnLOH/LOH 區的 allele-specific methylation 與 haplotype-specific 甲基行為。',
      '要找：(1) cnLOH 後兩 copy 來自同一親代 haplotype → 甲基是否該同質化(失去 ASM)？還是保留？',
      '(2) 文獻是否有「LOH/cnLOH 區仍觀察到 haplotype 間甲基差異」的報告？機制？',
      '(3) imprinting/X-inactivation 與 cnLOH 交互。', RULE, '回傳 SCHEMA。'].join('\n') },
]
const lit = await parallel(dims.map((d) => () =>
  agent(d.prompt, { label: 'lit:' + d.key, phase: 'LitSearch', schema: SCHEMA })
    .then((r) => ({ key: d.key, label: d.label, result: r }))))
const valid = lit.filter((e) => e && e.result)
phase('Synthesize')
const SYN = {
  type: 'object', additionalProperties: false,
  properties: {
    external_verdict: { type: 'string', description: '外部證據是否證實本 pilot 的核心解釋(層級錯配: normal germline 雜合 vs 腫瘤 cnLOH)' },
    hypothesis_support: { type: 'array', items: { type: 'object', additionalProperties: false,
      properties: { hypothesis: { type: 'string' }, external_verdict: { type: 'string', enum: ['supported','mixed','refuted','no-evidence'] }, evidence: { type: 'string' } },
      required: ['hypothesis','external_verdict','evidence'] } },
    new_possibilities: { type: 'array', items: { type: 'string' }, description: '文獻提出、本 pilot 尚未考慮的其他可能' },
    citations: { type: 'array', items: { type: 'string' } },
  }, required: ['external_verdict','hypothesis_support','new_possibilities','citations'],
}
const blob = valid.map((e) => '### [' + e.key + '] ' + e.label + '\n' + JSON.stringify(e.result, null, 2)).join('\n\n')
const syn = await agent(['你是文獻綜整員。本 pilot 發現「SEQC2 標純LOH 的區甲基卻能分開 HP1/HP2」，定性結論=層級錯配(phasing 用 normal germline 該區仍雜合 VAF~0.5; SEQC2 LOH=腫瘤 cnLOH)。tumor VAF 抽查分兩型: chr15 仍雜合(0.492無真cnLOH) / chr8 偏移(0.159真imbalance)。',
  '請用下方 4 份文獻調查輸出 SYN：',
  '1. external_verdict: 外部證據是否證實「層級錯配」核心解釋。',
  '2. hypothesis_support: 對 4 假說(H1 層級錯配/H2 SEQC2標記不準/H3 真imbalance/H4 兩拷貝已推翻)各給外部證據裁決。',
  '3. new_possibilities: 文獻提出、本 pilot 尚未想到的其他可能(直到找到多種可解釋可能)。',
  '4. citations: 可寫進報告的具體 citation。',
  '=== 4 份文獻 ===', blob].join('\n'), { label: 'synthesize-loh-lit', phase: 'Synthesize', schema: SYN })
return { lit: valid.map((e) => ({ key: e.key, summary: e.result.summary, takeaways: e.result.key_takeaways })), synthesis: syn }
