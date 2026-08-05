export const meta = {
  name: 'longphase-s-inheritance-grounding',
  description: '唯讀深度分析 longphase-S somatic_haplotag 既有 H3->H1-1/H2-1 inheritance + 量化信心訊號 + RR/RA/AR/AA 兩點 read 結構 + 甲基注入點 + ISM 消費鏈，產出修正設計 delta + 方法抉擇 ledger',
  phases: [
    { title: 'DeepRead', detail: '5 個平行深讀 agent 各析 longphase-S 一個面向' },
    { title: 'Verify', detail: '對抗驗證 2 個關鍵 claim' },
    { title: 'Synthesize', detail: '修正設計 delta + 方法抉擇 ledger' },
  ],
}

const LPS = '/big8_disk/liaoyoyo2001/longphase-s/src'
const REPO = '/big7_disk/liaoyoyo2001/InterSubMod'

const DEEP_SCHEMA = {
  type: 'object',
  additionalProperties: false,
  properties: {
    summary: { type: 'string' },
    mechanisms: {
      type: 'array',
      items: {
        type: 'object',
        additionalProperties: false,
        properties: {
          name: { type: 'string', description: '機制/規則/訊號名稱' },
          how_it_works: { type: 'string', description: '具體怎麼運作（含閾值/公式/決策分支）' },
          file_refs: { type: 'array', items: { type: 'string' }, description: 'absolute path:line' },
          quantitative_values: { type: 'string', description: '涉及的閾值/常數；無則 N/A' },
          relevance_to_task: { type: 'string', description: '與甲基分群救 unphase/HP3 重指派的關係' },
        },
        required: ['name', 'how_it_works', 'file_refs', 'quantitative_values', 'relevance_to_task'],
      },
    },
    key_findings: { type: 'array', items: { type: 'string' } },
    open_questions: { type: 'array', items: { type: 'string' } },
  },
  required: ['summary', 'mechanisms', 'key_findings', 'open_questions'],
}

const DISCIPLINE = '搜尋紀律：用 Grep/Glob 工具廣搜，Read 只讀相關片段。所有 file_refs 必須 absolute path:line 且真的開檔讀過（不可臆測行號）。閾值/常數必須引用實際程式碼，不可編造。'

phase('DeepRead')

const dims = [
  {
    key: 'inheritance-algo',
    label: 'H3->H1-1/H2-1 inheritance 演算法',
    prompt: [
      '你是 longphase-S somatic haplotag inheritance 演算法深讀員。',
      '',
      '目標：完整逆向工程 longphase-S 如何決定一條 read 的 ReadHP（enum: unTag=0, H1=1, H2=2, H3=3, H4=4, H1_1=5, H1_2=6, H2_1=7, H2_2=8）。重點是 H3（somatic read 但未繼承）-> H1-1/H2-1（繼承到 germline haplotype family）的 inheritance 決策，以及一條 read 何時變成 unTag(unphase) / H3。',
      '',
      '必讀檔（用 Read 開來逐段看，不要只 grep）：',
      LPS + '/somatic_haplotag/SomaticVarCaller.cpp（2433 行，最大；inheritance + filter 核心）',
      LPS + '/somatic_haplotag/SomaticHaplotagProcess.cpp（654 行；主 workflow）',
      LPS + '/somatic_haplotag/SomaticHaplotag.cpp（142 行）',
      LPS + '/haplotag/HaplotagStrategy.cpp（germline tagging 決策邏輯）',
      LPS + '/haplotag/HaplotagType.h（enum + ReadHpResult/SomaticData struct）',
      '',
      '具體要回答（每點給 file:line + 確切閾值）：',
      '1. read->ReadHP 的完整決策樹：germline het 投票（H1 vs H2）如何決定？什麼條件->unTag？什麼->H3？percentageThreshold 之類的閾值確切是多少、在哪。',
      '2. inheritance：H3 read 如何繼承成 H1_1 / H2_1？依據什麼證據（somatic SNP 的 derive-by-HP？germline imbalance？）？欄位 deriveHP / existDeriveByH1andH2 / somaticReadDeriveByHP 怎麼算與使用？',
      '3. before inheritance vs after inheritance 的 read distribution（pipeline 報告提過 readDistri_beforeInheritance/afterInheritance.out）對應程式碼哪段？',
      '4. unTag/H3 read 此刻被丟棄還是仍寫進 BAM？（HP:i tag 輸出邏輯）',
      '5. 若要用甲基分群輔助 inheritance 決策，最自然的注入點是哪個函式/行？',
      '',
      DISCIPLINE,
      '回傳 DEEP_SCHEMA。mechanisms 把每條決策規則拆開列。',
    ].join('\n'),
  },
  {
    key: 'confidence-signals',
    label: '既有量化信心訊號清單',
    prompt: [
      '你是 longphase-S 量化信心訊號盤點員。',
      '',
      '目標：把 longphase-S 已計算的、可作 tag 信心數值化的 per-read / per-variant 量化訊號全部列出（用戶要求可量化判斷 tag 信心 + 抉擇方法要紀錄可量化依據）。',
      '',
      '必讀：',
      LPS + '/somatic_haplotag/SomaticHaplotagProcess.h + .cpp',
      LPS + '/somatic_haplotag/SomaticVarCaller.cpp',
      LPS + '/haplotag/HaplotagType.h（SomaticData / PosBase / ReadHpResult struct 欄位）',
      LPS + '/somatic_haplotag/TumorPurityEstimator.cpp（imbalance ratio 算法）',
      '',
      'HaplotagType.h 已知候選欄位（逐一去 .cpp 確認怎麼算、範圍、用途）：',
      'germlineHaplotypeImbalanceRatio, percentageOfGermlineHp (PosBase)；allelicImbalanceRatio, somaticHaplotypeImbalanceRatio, shannonEntropy (SomaticData)；pure_H1_1_readRatio / pure_H2_1_readRatio / pure_H3_readRatio / Mixed_HP_readRatio (SomaticData)；deriveHPsimilarVec, somaticSnpH3count, readHpCounter (ReadHpResult)；meanAltCountPerVarRead, zScore, minDistance；各 filteredByXXX flag。',
      '',
      '對每個訊號：(a) 怎麼算（公式/file:line）、(b) 值域、(c) 現用在哪（filter? purity? tag?）、(d) 能否當 read-level 或 variant-level 信心 proxy、(e) per-read 還是 per-variant 聚合。',
      '',
      DISCIPLINE,
      '回傳 DEEP_SCHEMA。key_findings 給最適合當 tag 信心數值的 top 3-5 個現成訊號。',
    ].join('\n'),
  },
  {
    key: 'rr-ra-ar-aa',
    label: 'RR/RA/AR/AA 兩點 read 結構',
    prompt: [
      '你是 longphase-S 多位點 read 結構深讀員。',
      '',
      '目標：搞清楚一條 read 同時跨過多個 somatic 位點時，longphase-S 怎麼記錄它在各位點是 Ref 還是 Alt（用戶要區分經過兩點 read 的 Ref-Alt vs Alt-Ref 即 RR/RA/AR/AA，並關聯 H1-1 多點 read 差異 + 甲基分群）。',
      '',
      '線索（HaplotagType.h 已見）：',
      'SomaticData.PosSomaticOffsetBase: array<vector<pair<int,char>>,2> 0:ref 1:alt；SomaticData.alleleCount array<int,2> 0:ref 1:alt；SomaticData.somaticReadHpCount；ReadHpResult.readHpCounter；enum Allele REF_ALLELE=0 ALT_ALLELE=1。',
      '',
      '必讀：',
      LPS + '/somatic_haplotag/SomaticVarCaller.cpp（PosSomaticOffsetBase/alleleCount 怎麼填）',
      LPS + '/somatic_haplotag/SomaticHaplotagProcess.cpp',
      LPS + '/haplotag/HaplotagType.h',
      '',
      '要回答：',
      '1. 一條 read 跨多 somatic 位點的有序 ref/alt 序列現在是否被記錄？在哪個資料結構？（PosSomaticOffsetBase 是 per-position 還是 per-read？）',
      '2. 若要做 read-level RR/RA/AR/AA 組合統計，現有結構夠不夠，缺什麼。',
      '3. 這與 H1_1（somatic 繼承到 H1 family）多點 read 的關係。',
      '4. 跨兩 somatic 點的 read 在此 codebase 是否常見（有無相關計數/log）。',
      '',
      DISCIPLINE,
      '回傳 DEEP_SCHEMA。',
    ].join('\n'),
  },
  {
    key: 'methyl-injection',
    label: '甲基注入點 + modcall 現況',
    prompt: [
      '你是 longphase-S 甲基化支援深讀員。',
      '',
      '目標：確認 longphase-S 的 haplotag/somatic_haplotag voting 流程目前是否讀取/使用 per-read 甲基（MM/ML / base modification）；modcall module 是否如 CLAUDE.md 說是 placeholder；若要把甲基分群接進 inheritance 決策，最小注入點在哪、工程量多大。',
      '',
      '必讀：',
      LPS + '/modcall/ModCall.cpp + ModCallParsingBam.cpp + ModCallProcess.cpp（甲基怎麼解析；是否用 hts_base_mod_state）',
      LPS + '/haplotag/Haplotag.cpp + HaplotagVcfParser.h（grep 顯示有 mod 字樣，確認是 MOD_FILE pass-through 還是真進 voting）',
      LPS + '/haplotag/HaplotagParsingBam.cpp（read 解析時有無碰 MM/ML）',
      LPS + '/somatic_haplotag/SomaticVarCaller.cpp（inheritance 有無用到甲基）',
      LPS + '/main.cpp（modcall command dispatch）',
      '',
      '要回答（給 file:line）：',
      '1. somatic_haplotag inheritance 決策有沒有用到任何甲基訊號？（預期：沒有，請確認）',
      '2. MOD_FILE / MOD 選項實際做什麼？把 modcall 結果 pass 進 phase 還是只 haplotag BAM 帶 tag？',
      '3. modcall module 用什麼 API 解甲基、輸出什麼格式？production-ready 還是 placeholder？',
      '4. 若要對 H3/unTag read 用甲基分群輔助 inheritance，最小注入點（哪函式）+ 需新接什麼 + 估計工程量（行數量級）。',
      '5. 既有 modcall 的甲基解析能否被 somatic_haplotag 重用，還是要重接。',
      '',
      DISCIPLINE,
      '回傳 DEEP_SCHEMA。',
    ].join('\n'),
  },
  {
    key: 'ism-consumption-data',
    label: 'ISM 消費鏈 + paired HCC1395 資料',
    prompt: [
      '你是 ISM 消費鏈與資料盤點員。repo root = ' + REPO + '；longphase-S 在 /big8_disk/liaoyoyo2001/longphase-s。',
      '',
      '目標：確認 paired 模式下 ISM 吃哪個 longphase-S tagged BAM、HCC1395 paired 資料在哪、既有 unphase/HP3 統計與 before/after inheritance 輸出檔在哪，為 Phase A 觀察鋪路。',
      '',
      '必讀/必查：',
      REPO + '/scripts/pipeline/config.sh（LONGPHASE_S_BIN 指向哪、HCC1395 sample config、NORMAL_BAM/TUMOR_BAM/germline path）',
      REPO + '/scripts/pipeline/steps/00_prepare_germline.sh',
      '/big8_disk/liaoyoyo2001/longphase_somatic_output_ssrs/（列出實際檔案 + HP3.out / readDistri_beforeInheritance.out / afterInheritance.out / somaticRead.out / totalRead.out 各檔格式；可快速統計則給 HP3/unTag 數量級）',
      REPO + ' 內 ISM 怎麼讀 tagged BAM 的 HP tag（ReadParser；grep HP tag 解析）',
      '任何既有 unphase/HP3 比例統計（docs/experiments / haplotag_qc.py 輸出）',
      '',
      '要回答：',
      '1. paired HCC1395 的 longphase-S tagged BAM 實際路徑 + ISM 消費它的進入點。',
      '2. longphase_somatic_output_ssrs/ 各 .out 檔的格式與可用統計（現成 before/after inheritance ground）。',
      '3. 既有 HP3 / unTag 數量或比例的任何數字。',
      '4. Phase A 觀察可直接用的現成 artifact 清單。',
      '',
      DISCIPLINE,
      '回傳 DEEP_SCHEMA。',
    ].join('\n'),
  },
]

const deep = await parallel(
  dims.map((d) => () =>
    agent(d.prompt, { label: 'deep:' + d.key, phase: 'DeepRead', schema: DEEP_SCHEMA })
      .then((r) => ({ key: d.key, label: d.label, result: r }))
  )
)
const validDeep = deep.filter((e) => e && e.result)

phase('Verify')

const VERIFY_SCHEMA = {
  type: 'object',
  additionalProperties: false,
  properties: {
    claim: { type: 'string' },
    verdict: { type: 'string', enum: ['confirmed', 'refuted', 'partial'] },
    evidence: { type: 'string', description: '含 file:line' },
  },
  required: ['claim', 'verdict', 'evidence'],
}

const inheritanceFinding = JSON.stringify(validDeep.find((d) => d.key === 'inheritance-algo') ? validDeep.find((d) => d.key === 'inheritance-algo').result : {}, null, 2)
const methylFinding = JSON.stringify(validDeep.find((d) => d.key === 'methyl-injection') ? validDeep.find((d) => d.key === 'methyl-injection').result : {}, null, 2)

const verifications = await parallel([
  () => agent(
    '對抗驗證以下 claim，試著找反證。claim：longphase-S 的 somatic_haplotag inheritance 決策完全不使用任何 per-read 甲基（MM/ML）訊號。\n\n去 ' + LPS + '/somatic_haplotag/ 與 ' + LPS + '/haplotag/ 親自開檔確認 inheritance/voting 路徑有無任何 mod/methyl/MM/ML/base_mod 觸及。\n\n甲基分析 agent 初步發現：\n' + methylFinding + '\n\n回傳 VERIFY_SCHEMA（evidence 含 file:line）。',
    { label: 'verify:methyl-not-in-voting', phase: 'Verify', schema: VERIFY_SCHEMA }),
  () => agent(
    '對抗驗證以下 claim，試著找反證或補正確值。claim（來自 inheritance agent）：longphase-S read->HP 與 H3->H1-1/H2-1 inheritance 的決策閾值如其所述。\n\ninheritance agent 發現（含宣稱閾值）：\n' + inheritanceFinding + '\n\n去 ' + LPS + '/somatic_haplotag/SomaticVarCaller.cpp + SomaticHaplotagProcess.cpp + ' + LPS + '/haplotag/HaplotagStrategy.cpp 親自核對：(1) percentageThreshold 之類閾值確切值與行號是否正確；(2) inheritance 條件分支是否被正確描述；(3) 有無遺漏關鍵決策分支。\n\n回傳 VERIFY_SCHEMA（confirmed=閾值與分支描述正確；partial=部分需修正並在 evidence 給正確值；evidence 含 file:line）。',
    { label: 'verify:inheritance-thresholds', phase: 'Verify', schema: VERIFY_SCHEMA }),
])
const validVerif = verifications.filter(Boolean)

phase('Synthesize')

const SYNTH_SCHEMA = {
  type: 'object',
  additionalProperties: false,
  properties: {
    corrected_understanding: { type: 'string' },
    design_delta: { type: 'string' },
    strategic_fork: {
      type: 'object',
      additionalProperties: false,
      properties: {
        question: { type: 'string' },
        options: { type: 'array', items: { type: 'string' } },
        recommendation: { type: 'string' },
        observable_basis: { type: 'string' },
      },
      required: ['question', 'options', 'recommendation', 'observable_basis'],
    },
    method_decision_ledger: {
      type: 'array',
      items: {
        type: 'object',
        additionalProperties: false,
        properties: {
          decision_id: { type: 'string' },
          fork: { type: 'string' },
          options: { type: 'array', items: { type: 'string' } },
          observable_data_needed: { type: 'string' },
          quantitative_criterion: { type: 'string' },
          current_recommendation: { type: 'string' },
          rationale: { type: 'string' },
        },
        required: ['decision_id', 'fork', 'options', 'observable_data_needed', 'quantitative_criterion', 'current_recommendation', 'rationale'],
      },
    },
    confidence_metrics_shortlist: {
      type: 'array',
      items: {
        type: 'object',
        additionalProperties: false,
        properties: {
          metric: { type: 'string' },
          source: { type: 'string' },
          meaning: { type: 'string' },
        },
        required: ['metric', 'source', 'meaning'],
      },
    },
    go_nogo_gate: { type: 'string' },
    spec_ready: { type: 'string', enum: ['yes', 'needs-user-input'] },
    remaining_user_question: { type: 'string' },
  },
  required: ['corrected_understanding', 'design_delta', 'strategic_fork', 'method_decision_ledger', 'confidence_metrics_shortlist', 'go_nogo_gate', 'spec_ready', 'remaining_user_question'],
}

const deepBlob = validDeep.map((e) => '### [' + e.key + '] ' + e.label + '\n' + JSON.stringify(e.result, null, 2)).join('\n\n')
const verifBlob = validVerif.map((v) => JSON.stringify(v, null, 2)).join('\n\n')

const synthPrompt = [
  '你是研究規劃綜整員（Opus 4.8，ultrathink）。背景：InterSubMod 新任務「用甲基分群拯救 unphase/HP3 read 的相位/tag 重指派」。先前一輪探索誤判改動目標為 longphase-to-mod；本輪已確認真正目標是 /big8 longphase-S（somatic_haplotag 模式，ISM 吃它產的 paired tagged BAM），且它已內建 H3->H1-1/H2-1 inheritance 邏輯與多個量化訊號。',
  '',
  '已鎖定的用戶決策（不要推翻）：',
  '- success metric = phasing rescue + 觀察對 LOH-constrained phasing 活軸的增益（絕不可回頭包裝成 variant TP/FP F1，否則落入已 concluded 死路）。',
  '- 先 ISM 原型（低風險可逆），聚焦 paired 模式 + 正確的 longphase-S 軟體（非 TO；TO 留到 paired 有結果後再評估）。',
  '- 首版能力：unphase/HP3 局部重指派（核心）+ per-read tag 信心量化 + RR/RA/AR/AA 兩點 read 區分；每步拆最小單元 + commit + 驗證 + 可回溯（TDD 風格）。',
  '- 驗證 scope：HCC1395 pilot（partial flag）。',
  '- 甲基證據軸：HP-axis 為主（somatic-controlled），ALLELE-axis 受 germline-het 基線 ASM confound 只作描述 + 強制 negative control。',
  '- 新要求：抉擇方法時要紀錄觀察數據 + 可量化的抉擇依據與決定供用戶檢核 -> 這是 method_decision_ledger 的用途。',
  '',
  '請輸出 SYNTH_SCHEMA：',
  '1. corrected_understanding：基於深讀+驗證，最準確描述 longphase-S 既有 inheritance 怎麼運作、甲基缺口在哪、ISM 端能拿到什麼。',
  '2. design_delta：相對先前段落1-6 設計，因 longphase-S 已有 inheritance + 量化訊號哪些要改/哪些不變。關鍵：先前以為要從零建重指派，現在更可能是 characterize longphase-S 既有 inheritance 的 before/after + 找它把 read 丟去 H3/unTag 的弱點 + 評估甲基能否補。',
  '3. strategic_fork：此發現開的最關鍵策略分岔（很可能是增強 longphase-S C++ inheritance vs ISM post-hoc 甲基重指派層），含 recommendation + 要看哪些觀察數據才能定。',
  '4. method_decision_ledger：>=8 個方法抉擇點，每個含 options / 要觀察什麼數據 / 量化判準 / 當前建議 / 理由。涵蓋：增強 longphase-S vs ISM post-hoc、HP-axis vs ALLELE-axis、信心指標選哪現成訊號、甲基 binary vs raw matrix、只救 H3 vs +unTag、重指派接受門檻、LOH/CN read 怎麼 gate、無甲基覆蓋 read 怎麼處理。',
  '5. confidence_metrics_shortlist：優先用現成 longphase-S 訊號（給 file:line）。',
  '6. go_nogo_gate：Phase A 最便宜前置驗證 + 量化通過門檻。',
  '7. spec_ready + remaining_user_question：判斷是否還需用戶拍板一個點（strategic_fork）才能寫 spec。用戶已說寫 spec，但若 strategic_fork 影響 spec 結構則標 needs-user-input。',
  '',
  '=== 5 份深讀 ===',
  deepBlob,
  '',
  '=== 2 份對抗驗證 ===',
  verifBlob,
].join('\n')

const synth = await agent(synthPrompt, { label: 'synthesize-design-delta', phase: 'Synthesize', schema: SYNTH_SCHEMA })

return {
  deep: validDeep.map((e) => ({ key: e.key, summary: e.result.summary, key_findings: e.result.key_findings })),
  verifications: validVerif,
  synthesis: synth,
}
