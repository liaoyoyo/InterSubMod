export const meta = {
  name: 'methyl-phasing-literature-prenarrative',
  description: '文獻預敘述 + 差驗證：甲基分群救相位 read 的可行性、ASM 量級、copy/alignment 多群 confound、nanopore 每-read 甲基工具、對抗反證 — 先查本地 KB 再 web，產出 pre-narrative + 設計假設合理性裁決',
  phases: [
    { title: 'LitSearch', detail: '5 平行文獻 agent (KB-first 再 web)' },
    { title: 'Synthesize', detail: 'pre-narrative + 設計假設逐條合理性裁決 + 須補的失效模式' },
  ],
}

const KB = '/big8_disk/liaoyoyo2001/Knowledge'

const LIT_SCHEMA = {
  type: 'object',
  additionalProperties: false,
  properties: {
    summary: { type: 'string' },
    findings: {
      type: 'array',
      items: {
        type: 'object',
        additionalProperties: false,
        properties: {
          claim: { type: 'string', description: '文獻支持/反對的具體陳述' },
          evidence: { type: 'string', description: '哪篇/哪個來源、數字、方法' },
          source: { type: 'string', description: 'KB 路徑 或 paper title+year+venue 或 URL' },
          relevance: { type: 'string', description: '對本任務(甲基救 unphase/HP3 相位)的意涵' },
          tier: { type: 'string', enum: ['KB-verified', 'peer-reviewed', 'preprint', 'tool-doc', 'blog-weak'] },
        },
        required: ['claim', 'evidence', 'source', 'relevance', 'tier'],
      },
    },
    key_takeaways: { type: 'array', items: { type: 'string' } },
    failure_modes: { type: 'array', items: { type: 'string' }, description: '文獻記載的失效/confound 模式（給設計加防線）' },
    open_questions: { type: 'array', items: { type: 'string' } },
  },
  required: ['summary', 'findings', 'key_takeaways', 'failure_modes', 'open_questions'],
}

const RULE = [
  '紀律：先 grep 本地 KB（' + KB + '，尤其 05_tools / 06_workflows / 03_file_formats）再用 WebSearch/WebFetch（用 ToolSearch 載入）。',
  '每個 claim 必標 tier（KB-verified > peer-reviewed > preprint > tool-doc > blog-weak）與來源。',
  '生醫主題優先找 PubMed / 期刊；工具行為找官方 doc。不可編造 citation；找不到就說找不到。',
  '聚焦「read-level / 單分子 甲基 → haplotype 分群」可行性，不要泛談 bulk 甲基。',
].join('\n')

phase('LitSearch')

const dims = [
  {
    key: 'methyl-read-phasing',
    label: '甲基分群 read→haplotype 可行性',
    prompt: [
      '你是文獻調查員。主題：能否用「單分子/每-read 甲基化模式」把缺乏 germline SNP 證據的 long read 分群到正確 haplotype？',
      '要找：',
      '1. 有沒有方法/論文用 per-read CpG 甲基把 read 分到 HP1 vs HP2（methylation-based read phasing / haplotype-aware methylation）。',
      '2. 報告的準確度（多少 read 可救、一致率、AUC/F1）。',
      '3. 用什麼距離/聚類（Hamming/Jaccard/Bernoulli/Euclidean、hierarchical/k-means/GMM）。',
      '4. nanopore 長讀長特有的做法（modkit / nanomethphase / methphaser / 類似）。',
      '特別查 nanomethphase、methphaser、Pyle、Akbari haplotype methylation 等關鍵字。',
      RULE,
      '回傳 LIT_SCHEMA。',
    ].join('\n'),
  },
  {
    key: 'asm-magnitude',
    label: 'ASM 量級夠不夠分群',
    prompt: [
      '你是文獻調查員。主題：allele-specific methylation (ASM) / haplotype 間甲基差異的典型「量級」夠不夠把 read 分群？',
      '要找：',
      '1. ASM 位點 haplotype 間 Δβ（甲基比例差）典型多大？imprinting 區 vs 一般 ASM vs somatic。',
      '2. 全基因組有多少比例位點呈顯著 ASM（這決定 unphase read 落在 ASM 區的機率）。',
      '3. somatic（腫瘤）ASM vs germline ASM 的差異與 confound（baseline allelic methylation）。',
      '4. 單一 read 上有幾個 CpG 才足以可靠分類（資訊量下限）。',
      '本專案已知 memory：HP-axis somatic-controlled OR≈1.79、ALLELE-axis 被 germline-het null confound(TP 11.1% < null 15.2%)、BRCA2 Δβ≈−0.122 — 請對照文獻看這量級在 read 分群是強還弱。',
      RULE,
      '回傳 LIT_SCHEMA。',
    ].join('\n'),
  },
  {
    key: 'copy-alignment-multimodal',
    label: 'copy/alignment 致同-tag 多群 confound',
    prompt: [
      '你是文獻調查員。主題：為什麼「同一 HP tag 的 read 在甲基空間會分成密集的多群」，以及這對甲基分群的 confound。',
      '要找：',
      '1. copy-number gain / 多拷貝 / segmental duplication 如何讓同位置 read 來自不同 paralog → 甲基模式分裂成多群。',
      '2. multi-mapping / mismapping / repeat 區 alignment artifact 如何造成假的甲基次群。',
      '3. cnLOH / 腫瘤異質性 (subclone) 如何讓同 haplotype 內甲基本就雙峰。',
      '4. 文獻怎麼偵測並排除這類多模態（modality test、coverage anomaly、mapping quality、paralog-aware）。',
      '這直接對應用戶新洞察：同 tag 但密集分兩/多群可能是 copy 或 alignment，須偵測避免汙染重指派。',
      RULE,
      '回傳 LIT_SCHEMA，failure_modes 要具體列偵測手段。',
    ].join('\n'),
  },
  {
    key: 'nanopore-methyl-tools',
    label: 'nanopore 每-read 甲基工具 + 視覺化',
    prompt: [
      '你是工具文獻調查員。主題：nanopore per-read 甲基「萃取 + 分群 + 視覺化（熱圖）」的標準工具鏈。',
      '要找（先查 KB ' + KB + '/05_tools）：',
      '1. MM/ML tag 解析：modkit（extract / pileup）、samtools、pysam mod API — 哪個產 per-read×CpG 矩陣最方便。',
      '2. per-read 甲基熱圖視覺化工具：methylartist、NanoMethViz、modbamtools、IGV methylation mode。',
      '3. IGV 對 modified base 的顯示（IGV.js / igv-reports 能否做 standalone 甲基 BAM view）。',
      '4. 標準閾值（如 modkit 5mC qual 門檻、本專案 longphase-s modcall qual>=204=mod/<=51=unmod）。',
      RULE,
      '回傳 LIT_SCHEMA，key_takeaways 給「最省工的 per-read 矩陣 + 熱圖 + IGV」工具建議。',
    ].join('\n'),
  },
  {
    key: 'adversarial-against',
    label: '對抗：反證甲基救相位會失敗',
    prompt: [
      '你是對抗性文獻調查員 — 任務是找「甲基分群救 unphase/HP3 read 會失敗或不可靠」的證據，挑戰本設計。',
      '要找：',
      '1. 甲基訊號在哪些情境對 read 分類無判別力（低 CpG 密度區、homogeneous methylation、CpG island vs open sea）。',
      '2. 甲基 noise / basecalling 錯誤率對單-read 分類的衝擊。',
      '3. 文獻中「甲基無法可靠 phase read」的反例或 caveat。',
      '4. 本專案已 concluded：pure-methylation TP/FP AUC<=0.58 天花板 — 文獻是否也有類似甲基判別力上限的記載。',
      '請盡量 refute「unphase read 可被甲基救」這個假設；找不到反證才回報 inconclusive。',
      RULE,
      '回傳 LIT_SCHEMA，failure_modes 是本 agent 的主要產出。',
    ].join('\n'),
  },
]

const lit = await parallel(
  dims.map((d) => () =>
    agent(d.prompt, { label: 'lit:' + d.key, phase: 'LitSearch', schema: LIT_SCHEMA })
      .then((r) => ({ key: d.key, label: d.label, result: r }))
  )
)
const valid = lit.filter((e) => e && e.result)

phase('Synthesize')

const SYNTH_SCHEMA = {
  type: 'object',
  additionalProperties: false,
  properties: {
    pre_narrative: { type: 'string', description: '現狀 + 文獻定位的預敘述（本任務在文獻地圖中的位置：誰做過、與我們的差異、我們補什麼）' },
    design_assumption_verdicts: {
      type: 'array',
      items: {
        type: 'object',
        additionalProperties: false,
        properties: {
          assumption: { type: 'string', description: '本設計的假設（如 unphase 可被甲基分群救 / HP-axis 比 ALLELE 乾淨 / 甲基量級夠分群）' },
          literature_verdict: { type: 'string', enum: ['supported', 'mixed', 'refuted', 'inconclusive'] },
          evidence: { type: 'string' },
          design_action: { type: 'string', description: '據此設計要加什麼防線/修正' },
        },
        required: ['assumption', 'literature_verdict', 'evidence', 'design_action'],
      },
    },
    failure_modes_to_add: { type: 'array', items: { type: 'string' } },
    recommended_tools: { type: 'array', items: { type: 'string' }, description: 'per-read 甲基矩陣/熱圖/IGV 工具建議' },
    citations_for_report: { type: 'array', items: { type: 'string' }, description: '可寫進報告的 citation（title year venue 或 KB path）' },
  },
  required: ['pre_narrative', 'design_assumption_verdicts', 'failure_modes_to_add', 'recommended_tools', 'citations_for_report'],
}

const blob = valid.map((e) => '### [' + e.key + '] ' + e.label + '\n' + JSON.stringify(e.result, null, 2)).join('\n\n')

const synth = await agent([
  '你是文獻綜整員（Opus 4.8）。背景：InterSubMod 新任務「用 per-read 甲基分群救 longphase-S 放棄的 unphase + 救不動的 H3 read 的相位歸屬」（metric=phasing rescue 非 variant F1；HP-axis 主；HCC1395 paired pilot）。下面是 5 份文獻調查。',
  '',
  '請輸出 SYNTH_SCHEMA：',
  '1. pre_narrative：本任務在文獻地圖的定位（已有 methylation read-phasing 工作嗎？methphaser/nanomethphase 等與我們差異？我們補什麼 gap？）。誠實標明哪些是別人做過的。',
  '2. design_assumption_verdicts：對本設計每條核心假設給文獻裁決 — (a) unphase read 可被甲基分群救、(b) HP-axis 比 ALLELE-axis 乾淨、(c) ASM 量級夠把 read 分群、(d) 同-tag 多群可偵測排除、(e) 甲基判別有 AUC 上限。每條給 supported/mixed/refuted/inconclusive + 證據 + 設計要加的防線。',
  '3. failure_modes_to_add：文獻記載、本設計該加的失效模式防線。',
  '4. recommended_tools：per-read 甲基矩陣 + 熱圖 + IGV 的最省工工具。',
  '5. citations_for_report：可寫進報告的具體 citation。',
  '',
  '=== 5 份文獻調查 ===',
  blob,
].join('\n'), { label: 'synthesize-literature', phase: 'Synthesize', schema: SYNTH_SCHEMA })

return { lit: valid.map((e) => ({ key: e.key, summary: e.result.summary, takeaways: e.result.key_takeaways, failure_modes: e.result.failure_modes })), synthesis: synth }
