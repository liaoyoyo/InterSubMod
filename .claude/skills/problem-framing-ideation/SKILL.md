---
name: problem-framing-ideation
description: **P0 REGISTER 前置 — 問題框架化（problem framing）**：用 5W1H + gap analysis 把模糊想法收斂成 1-3 個可被 inject-hypothesis 接收的具體假說候選。**限投入 <2hr**（個人風格 anchor #5：避免無限迭代；超時則先選最好的進 inject-hypothesis）。觸發：「腦力激盪」「brainstorm」「研究構想」「gap analysis」「5W1H」「開始新方向」「框架化想法」「research-ideation」（前身）。**職責邊界**：本 skill 輸出 `候選假說清單 + 5W1H 表`；不註冊（→ inject-hypothesis）、不執行（→ research-loop Step 4+）、不建專案目錄（→ init-research）。
version: 0.2.0
---

# Research Ideation

Supports the complete workflow for the research project initiation phase, from literature review to research question definition, method selection, and research planning.

## Core Features

### 1. Idea Brainstorming (5W1H Framework)

Systematically brainstorm research ideas using the 5W1H framework:
- **What**: What problem or phenomenon to study
- **Why**: Why this problem is important
- **Who**: Target audience and stakeholders
- **When**: Time scope and context of the research
- **Where**: Application scenarios and domains
- **How**: Preliminary research methodology ideas

**Integration with superpowers:brainstorming**: Can invoke the superpowers:brainstorming skill for interactive brainstorming to help rapidly generate and evaluate research ideas.

### 2. Literature Review

Systematically search, analyze, and synthesize related literature:
- Build effective search keywords
- Search via WebSearch across academic databases (arXiv, Google Scholar, etc.)
- Screen and evaluate paper quality
- Identify research trends and gaps
- Generate structured literature reviews
- **Zotero Integration**: Papers are automatically added to Zotero via DOI, organized into topic-based collections, and open-access PDFs are auto-attached for full-text reading

### 3. Gap Analysis

Systematically identify and evaluate research gaps:
- **Literature gaps**: Identify topics or questions not yet sufficiently studied
- **Methodological gaps**: Discover limitations and improvement opportunities in existing methods
- **Application gaps**: Identify opportunities for theory-to-practice transfer
- **Interdisciplinary gaps**: Discover research opportunities at the intersection of different fields
- **Temporal gaps**: Identify new research needs arising from changes over time

**Analysis Dimensions:**
- Coverage of research topics
- Comparison of strengths and weaknesses of existing methods
- Completeness of experimental setups
- Availability of datasets and benchmarks
- Gap between theory and practice

### 4. Research Question Definition

Formulate specific research questions based on literature analysis:
- Identify research gaps and opportunities
- Apply SMART principles to formulate questions
- Evaluate importance, novelty, and feasibility
- Define research objectives and expected contributions

### 5. Method Selection

Select appropriate research methods:
- Analyze strengths and weaknesses of existing methods
- Evaluate method applicability
- Identify required technologies and resources
- Consider method feasibility

### 6. Research Planning

Develop detailed research plans:
- Plan research timeline
- Define milestones and deliverables
- Identify potential risks
- Assess resource requirements

## When to Use

### Scenarios for This Skill

Use the research-ideation skill in the following situations:

1. **Starting a new research project** - Have research interests but no clear research question yet
2. **Literature review** - Need to systematically understand a research field
3. **Research question formulation** - Need to transform vague ideas into specific research questions
4. **Method selection** - Need to choose appropriate research methods and technical approaches
5. **Research planning** - Need to plan research timeline and resources

### Typical Workflow

```
Research interest → Idea brainstorming (5W1H) → Literature review → Gap analysis → Define question → Select method → Create plan
```

**Output Files:**
- `literature-review.md` - Structured literature review
- `research-proposal.md` - Research proposal (including question, method, plan)
- `references.bib` - References in BibTeX format
- Zotero collection with organized papers and PDFs

## Integration with Other Systems

### Complete Research Workflow

```
research-ideation (Research initiation)
    ↓
Experiment execution (completed by user)
    ↓
results-analysis (Results analysis)
    ↓
ml-paper-writing (Paper writing)
```

### Data Flow

- **research-ideation output** → Guides experiment design and method selection
- **Experimental results** → results-analysis for statistical analysis
- **Analysis results** → Related Work and Methods sections of ml-paper-writing

### Zotero Integration

Through the Zotero MCP server, the research-ideation workflow automates literature management:

- **Paper Discovery**: WebSearch finds relevant papers across academic databases
- **Auto-Import**: Extract DOIs from search results, use `add_items_by_doi` to add papers with full metadata
- **Collection Organization**: `create_collection` creates topic-based collections with standard sub-collections (Core Papers, Methods, Applications, Baselines, To-Read)
- **PDF Attachment**: `find_and_attach_pdfs` automatically finds and attaches open-access PDFs via Unpaywall
- **Full-Text Reading**: `get_item_fulltext` reads indexed PDF content for analysis and note-taking
- **Library Search**: `search_library` and `get_collection_items` browse existing papers to avoid duplicates

### Key Configuration

- **Literature search scope**: Papers from the last 3 years by default, configurable
- **Output format**: Markdown format for easy editing and version control
- **Citation management**: Generates references in BibTeX format
- **Zotero collection naming**: `Research-{topic}-{YYYY}` format
- **PDF auto-attach**: Enabled by default for open-access papers via Unpaywall

## Additional Resources

### Reference Files

Detailed methodology guides, loaded on demand:

- **`references/5w1h-framework.md`** - 5W1H Framework Guide
  - What, Why, Who, When, Where, How — six dimensions
  - Systematic approach to brainstorming research ideas
  - Integration with superpowers:brainstorming
  - Usage examples and best practices

- **`references/literature-search-strategies.md`** - Literature Search Strategies
  - Keyword construction techniques
  - Academic database selection (arXiv, Google Scholar)
  - Search tips and screening criteria
  - Paper quality evaluation methods
  - DOI extraction and Zotero auto-import workflow

- **`references/zotero-integration-guide.md`** - Zotero MCP Integration Guide
  - Available Zotero MCP tools (browse, add, cite)
  - Collection organization strategy and naming conventions
  - Automated workflow: WebSearch → DOI → Zotero import → PDF attach
  - Full-text reading and structured note-taking
  - Common issues and troubleshooting

- **`references/gap-analysis-guide.md`** - Gap Analysis Guide
  - 5 types of Gap Analysis (literature, methodological, application, interdisciplinary, temporal)
  - 5 analysis dimensions
  - Systematic approach to identifying research opportunities
  - Usage examples and best practices

- **`references/research-question-formulation.md`** - Research Question Formulation
  - Applying SMART principles
  - Question type classification (exploratory, confirmatory, applied)
  - Evaluation criteria (importance, novelty, feasibility)
  - Defining research objectives and contributions

- **`references/method-selection-guide.md`** - Method Selection Guide
  - Common research method classification
  - Method applicability analysis
  - Strengths and weaknesses comparison
  - Resource requirement assessment

- **`references/research-planning.md`** - Research Planning
  - Timeline planning methods
  - Milestone definition techniques
  - Risk identification and mitigation
  - Resource allocation strategies

### Example Files

Complete working examples:

- **`examples/example-literature-review.md`** - Literature Review Example
  - Demonstrates structured literature review format
  - Includes research trend analysis and gap identification

- **`examples/example-research-proposal.md`** - Research Proposal Example
  - Demonstrates complete research proposal structure
  - Includes complete examples of question, method, and plan

---

## Phase & Chain Position

- **Phase**: **P0 REGISTER 前置**（pre-Phase 0；尚未進入 cycle state machine）
- **Chain**: forward-link chain #1 第 1 環
  ```
  problem-framing-ideation （5W1H 框架化）
      ↓
  inject-hypothesis （選定假說註冊到 hypothesis_queue.json）
      ↓
  /cycle-init （建 state/cycles/{id}/state.json）
      ↓
  research-loop / validation-protocol （P1 PLAN）
      ↓
  /check-staleness （P2 PRECHECK gate）
      ↓
  ... (P3 → P4 → P5 → P6)
  ```
- **上游觸發**: 用戶口述新想法 / pivot-direction 後重新定向 / review-evidence 發現未驗證 gap
- **下游 skill**: `inject-hypothesis`（必走）；專案級長期方向則接 `init-research`

## Dependencies

| 類別 | 項目 |
|---|---|
| **Uses** (本 skill 內部呼叫) | WebSearch / WebFetch（查文獻 gap）、Read（讀 docs/concepts/）、Grep（在 evidence_ledger 搜相似假說） |
| **Used by** (誰會觸發本 skill) | 用戶手動 / `pivot-direction`（定向後）/ `review-evidence`（發現 gap 後） |
| **Reads** | `docs/concepts/2026/04/20260409_研究構想總索引_01.md`、`docs/CURRENT_FOCUS.md`、`research/autoresearch/research_direction.md`、外部文獻 |
| **Writes** | `research/autoresearch/research_direction.md`（append candidate slot）。**不直接寫** hypothesis_queue.json（交由 inject-hypothesis） |

## Failure Mode & Diagnostics

| 失敗症狀 | 先看哪 | 排查步驟 |
|---|---|---|
| 用戶輸入太模糊無法收斂 | 對話歷史前 5 turn + `docs/CURRENT_FOCUS.md` | 用 `grill-me` skill 補充缺項；若仍無法 → 升級到 confirmation-protocol Hard Gate |
| 5W1H 表填不出來（缺 What/Why） | `docs/concepts/2026/04/20260409_研究構想總索引_01.md` | 查既有相似研究方向，借用其 Why 框架 |
| 候選假說 > 5 個（太發散） | 用戶優先序 + 個人風格 anchor #5（限 <2hr） | 限投入時間；超時則選 top 3 進 inject-hypothesis，其餘 archive 到 research_direction.md backlog |
| 與已 NEGATIVE 方向重複 | `docs/experiments/INDEX.md` 的 ❌ NO-GO 區、`MEMORY.md` Concluded 區 | 搜「相似 keywords」，找到歷史結論並引用；若重複 → 警告用戶並回到框架化 |
| 與當前 Active Cycle 衝突 | `state/active.json`（≤5 個 active）+ `docs/CURRENT_FOCUS.md` | 若 active 已滿 5 個 → 不註冊新假說，先收尾舊 cycle |

**何時升級到別的 skill / agent / 人工審查**：
- 連續 3 次無法收斂候選假說 → 跳到 `grill-me` skill（深度互審）
- 涉及 NO-GO 判定（推翻已驗證結論） → Hard Gate，必停問用戶
- 候選假說涉及 C++ 修改 → 預先呼叫 `methodology-audit` 評估方法學可行性

**個人風格適配**（依 `feedback_*` memory）：
- Anchor #5 「One-turn mechanism freeze」 → **限 <2hr 投入**；超時就把現有最好的 push 到 inject-hypothesis 而非繼續 ideate
- Anchor #3 「報告骨架」 → 5W1H 表輸出格式須能直接餵給 structured-tech-report 的 §Background+Mechanism 段
- Anchor #7 「pivot 容忍」 → 對與已撤回方向相似的提案，警告但不禁止（用戶有權重啟）

## DO NOT USE WHEN（v1.7 batch A）

- **已有具體假說** — 用 `/inject-hypothesis` 直接註冊到 hypothesis_queue.json
- **想啟動 cycle** — 用 `/cycle-init`（P0 REGISTER）
- **想跑分析** — 用 `/feature-layered-observation`（P3 PILOT）
- **已超 2hr brainstorm** — anchor #5 強制收斂；推現有最好的進 inject-hypothesis
- **想擴大研究 scope** — 用 `/init-research` 建多週專案

## Quality Checklist — 交付候選假說清單前自我檢查（v1.7 batch B）

- [ ] 候選假說 **1-3 個**（嚴格不超過 3，避免決策疲勞）
- [ ] 每個假說有完整 5W1H 表（What / Why / Who / When / Where / How）
- [ ] 與 MEMORY.md `## Concluded` 段已 NEGATIVE 假說對照（避免重蹈覆轍）
- [ ] 投入時間 < 2 hr（anchor #5）；超時不再 ideate
- [ ] 輸出 markdown 表格可直接餵 inject-hypothesis（含 priority / dataset / expected_effect 欄）
- [ ] gap analysis 含「為何這方向值得試」+ 「最壞情況 NEGATIVE 也學到什麼」
- [ ] 預期 effect size 範圍給上下界
- [ ] 不註冊（→ inject-hypothesis）、不執行（→ research-loop）、不建專案（→ init-research）

