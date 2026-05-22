# Framework Catalog — 50+ 業界敘述框架完整定義

> **單一 SoT**（source of truth）— `/narrative-frame` skill 的核心 catalog。
> 10 大類；每框架附：定義 / 結構 / 適合場景 / 不適合 / 業界源 / 範例 / template 路徑。

> **路徑硬性規則**：本檔內所有 .md ref 用 `InterSubMod/...` 前綴。

---

## 速查索引

| 類 | Frameworks |
|----|------|
| [§1 結論先行 / 商務報告](#1-結論先行--商務報告高層-exec-friendly) | SCQA / Pyramid / MECE / BLUF / Inverted Pyramid |
| [§2 經驗敘事 / 行為案例](#2-經驗敘事--行為案例履歷--面試--績效) | STAR / CAR / PAR / SOAR / STARR |
| [§3 短時口頭答辯](#3-短時口頭答辯30s-3min) | PREP / OREO / ABT / Elevator Pitch / Rule of Three |
| [§4 故事敘事](#4-故事敘事簡報--ted--影片) | Pixar Spine / 三幕劇 / 英雄之旅 / Freytag / Duarte Sparkline / Golden Circle / Monroe / Assertion-Evidence |
| [§5 問題分析因果](#5-問題分析因果) | 5W1H / 5 Whys / Fishbone / A3 / DMAIC / PDCA / OODA |
| [§6 教學解釋](#6-教學解釋概念) | Feynman / Diátaxis / Bloom / KWL / ELI5 / CPP |
| [§7 回饋對話 Coaching](#7-回饋對話--coaching) | SBI / DESC / GROW / Pendleton / Sandwich / Radical Candor |
| [§8 行銷説服](#8-行銷--説服--銷售) | AIDA / AIDCA / FAB / PAS / 4P / BAB |
| [§9 決策評估](#9-決策評估) | WRAP / Pre-Mortem / Cynefin / ADR / DACI / DECIDE / Eisenhower |
| [§10 簡報架構](#10-簡報架構--演講設計) | TED / 3-Act Pitch / Resonate / Presentation Zen / Lessig / Takahashi / PechaKucha / PEEL |

---

## §1 結論先行 / 商務報告（高層 exec-friendly）

### SCQA — Situation / Complication / Question / Answer

- **結構**: Situation（共識基線）→ Complication（張力 / 變化）→ Question（核心問題）→ Answer（解答 / 建議）
- **適合**: 諮詢報告開場、提案 hook、PI 報告 TL;DR、説服性簡報、abstract
- **不適合**: 純技術 RFC / spec、探索性開放對話、單純 status report
- **源**: Barbara Minto《The Pyramid Principle》(2020), McKinsey 1960s
- **Template**: `templates/scqa_skeleton.md`
- **範例**: `examples/pi_report_scqa_example.md`

### Pyramid Principle — Top-down conclusion-first

- **結構**: 結論（最頂）→ 3-7 個 supporting arguments → 證據（最底）
- **適合**: exec summary、決策 deck、PI 1-on-1、論文 abstract
- **不適合**: 探索性 / 開放式對話、敘事性 keynote
- **源**: Barbara Minto《The Pyramid Principle: Logic in Writing and Thinking》
- **核心原則**: 思考自下而上（findings → recommendations → 結論）；溝通自上而下（結論先講）

### MECE — Mutually Exclusive, Collectively Exhaustive

- **結構**: 分類維度互斥（無重疊）+ 窮盡（覆蓋全部）
- **適合**: issue tree、competitor 分類、root cause 分支、ML feature 設計
- **不適合**: 已知答案的線性敘述、純故事
- **源**: Barbara Minto / McKinsey

### BLUF — Bottom Line Up Front

- **結構**: 結論句（1 行）→ 細節（後）
- **適合**: 軍方/政府 brief、Slack 同步、email subject、issue comment、status update
- **不適合**: 説服性故事、敘事性報告
- **源**: US Army Field Manual 6-22

### Inverted Pyramid

- **結構**: 最重要 → 次要 → 細節（5W1H 在頂端 lead 段）
- **適合**: 新聞、blog post、status update、長 .md 開頭段
- **不適合**: 推理鋪陳論證、戲劇性敘事
- **源**: Journalism late-1800s
- **重要備註**: 讀者可在任一點離開仍能理解全文

---

## §2 經驗敘事 / 行為案例（履歷 / 面試 / 績效）

### STAR — Situation / Task / Action / Result

- **結構**: Situation（背景）→ Task（目標）→ Action（你做了什麼）→ Result（量化結果）
- **適合**: 面試行為題、PR self-review、case study、performance review
- **不適合**: 概念解釋、純技術 spec
- **源**: 1970s Behavioral Interviewing methodology

### CAR — Challenge / Action / Result

- **結構**: Challenge（問題）→ Action → Result
- **適合**: 履歷 bullet（壓縮版 STAR）、LinkedIn 成就條目
- **不適合**: 詳述任務分配

### PAR — Problem / Action / Result

- **結構**: Problem → Action → Result
- **適合**: 履歷成就條目、簡短回顧
- **不適合**: 含團隊脈絡的長敘事

### SOAR — Situation / Obstacle / Action / Result

- **結構**: Situation → Obstacle（強調障礙）→ Action → Result
- **適合**: 含「障礙」敘事的面試、創業故事
- **不適合**: 流暢無阻情境

### STARR — STAR + Reflection

- **結構**: STAR + Reflection（反思）
- **適合**: 教育界 portfolio、reflective practice、PhD viva
- **不適合**: 求快面試

---

## §3 短時口頭答辯（30s-3min）

### PREP — Point / Reason / Example / Point

- **結構**: Point（直答結論）→ Reason（理由）→ Example（證據）→ Point（重申）
- **適合**: Q&A 答辯、會議發言、教授追問、podcast soundbite
- **不適合**: 多重論點、長 keynote
- **Template**: `templates/prep_skeleton.md`

### OREO — Opinion / Reason / Example / Opinion

- **結構**: Opinion → Reason → Example → Opinion
- **適合**: 立場性發言、辯論、社論
- **不適合**: 客觀數據報告
- **差別**: vs PREP — OREO 強調主觀立場；PREP 強調事實結論

### ABT — And / But / Therefore

- **結構**: AND（設定共識）→ BUT（衝突 / 對比）→ THEREFORE（結論 / 後果）
- **適合**: 1 句話講完研究、科普 hook、verdict line、paper abstract
- **不適合**: 多重變數、複雜分析
- **源**: Randy Olson《Narrative Is Everything: The ABT Framework》(2019), 《Houston, We Have a Narrative》
- **核心**: 三幕劇精煉版；每個 scene、段、論文都有 ABT 結構
- **範例**: 「There are 3 people in a room AND it's late afternoon, BUT the woman appears to be questioning the man, THEREFORE the scene looks like an interrogation.」
- **Template**: `templates/abt_skeleton.md`

### Elevator Pitch（30-60s）

- **結構**: Hook（吸引）→ What（你是誰 / 做什麼）→ Why care（受眾為何在意）→ Ask（行動）
- **適合**: 電梯 / 雞尾酒會 / first contact、conference networking
- **不適合**: 深度技術討論

### Rule of Three

- **結構**: 3 個並列重點
- **適合**: 任何記憶優先場景（人腦對 3 點 retention 最高）、口頭結論、bullet list
- **不適合**: 需要窮盡分類（用 MECE）

---

## §4 故事敘事（簡報 / TED / 影片）

### Pixar Story Spine

- **結構**: "Once upon a time…" → "Every day…" → "One day…" → "Because of that…" × N → "Until finally…" → "And ever since…"
- **適合**: 科普影片、TED、impact story、品牌故事、變革敘事
- **不適合**: 純數據報告、quick brief
- **源**: Kenn Adams (improv community 1990s); Pixar 4th rule of storytelling
- **基礎**: 建立在 ABT 之上 — adds twists 在 BUT 與 THEREFORE 之間
- **Template**: `templates/pixar_spine_skeleton.md`

### 三幕劇（Three-Act Structure）

- **結構**: Act 1 Setup（25%）→ Act 2 Confrontation（50%）→ Act 3 Resolution（25%）
- **適合**: 90 分鐘演講、紀錄片、keynote、長 case
- **不適合**: 短講、技術文件
- **源**: Syd Field《Screenplay: The Foundations of Screenwriting》

### 英雄之旅（Hero's Journey）— 12 步

- **結構**: Ordinary World → Call to Adventure → Refusal → Meeting Mentor → Crossing Threshold → Tests/Allies/Enemies → Approach → Ordeal → Reward → Road Back → Resurrection → Return with Elixir
- **適合**: 創業故事、轉型敘事、品牌神話、人物 portrait
- **不適合**: 學術報告、客觀分析
- **源**: Joseph Campbell《The Hero with a Thousand Faces》(1949); Christopher Vogler《The Writer's Journey》

### Freytag's Pyramid

- **結構**: Exposition → Rising Action → Climax → Falling Action → Resolution
- **適合**: 戲劇、小説、長篇 case、紀錄影片
- **不適合**: 短講、商務 brief
- **源**: Gustav Freytag《Die Technik des Dramas》(1863)

### Duarte Sparkline / Persuasive Story Pattern

- **結構**: What is（現狀）↔ What could be（願景）反覆對比 + Call to Adventure + New Bliss
- **適合**: 説服性 keynote、變革倡議、TED、product launch
- **不適合**: 數據 report、技術深探
- **源**: Nancy Duarte《Resonate》(2010)
- **核心**: 經典演講（MLK "I Have a Dream", Steve Jobs iPhone launch）的共同 pattern

### Golden Circle — Why / How / What

- **結構**: Why（核心信念）→ How（做法）→ What（產品 / 結果）
- **適合**: 願景演講、創業 pitch、品牌、leadership
- **不適合**: 純技術文件、客觀比較
- **源**: Simon Sinek《Start with Why》(2009) + TED talk (2009)
- **核心**: 偉大領導者 / 品牌（Apple / Wright Bros / MLK）都從 Why 開始
- **Template**: `templates/golden_circle_skeleton.md`

### Monroe's Motivated Sequence

- **結構**: Attention → Need → Satisfaction → Visualization → Action
- **適合**: 募款、政策倡議、行銷、call-to-action 演講
- **不適合**: 純資訊報告
- **源**: Alan Monroe, Purdue University (1935)

### Assertion-Evidence

- **結構**: Slide title = 一句完整斷言（聲明結論）；body = 視覺證據支持
- **適合**: 科學簡報、學術 conference、PI 報告、academic defense
- **不適合**: 探索性開放討論
- **源**: Michael Alley《The Craft of Scientific Presentations》, Penn State, NSF-backed
- **效果**: Kuwait Univ 2025 實驗 108 學生確認 lower cognitive load
- **反例**: ❌「Results」 → ✅「V6 修正 100% chr19 priority bug victims (752/752 reads)」

---

## §5 問題分析因果

### 5W1H — Who / What / When / Where / Why / How

- **結構**: 6 個維度問題
- **適合**: 事件 brief、issue 描述、bug report、新聞 lead、journalism
- **不適合**: 已知 5W1H 後的深度分析
- **源**: 古希臘修辭學 → 現代新聞學 standard

### 5 Whys

- **結構**: 連問 5 次「為什麼」逐層追根因
- **適合**: 根因分析、postmortem、Toyota gemba（現場 walk）、debug
- **不適合**: 多因素並存（用 Fishbone）
- **源**: Toyota Production System / 大野耐一

### Fishbone / Ishikawa Diagram

- **結構**: 6M 分類 — Man / Machine / Method / Material / Measurement / Milieu（環境）
- **適合**: 複雜根因（多因素並存）、工程 RCA、product defect
- **不適合**: 單一線性原因（用 5 Whys）
- **源**: 石川馨 1968

### A3 Report

- **結構**: Background → Current State → Goal → Analysis（5 Whys + Fishbone）→ Proposal → Plan → Follow-up（全寫在一張 A3 紙）
- **適合**: Toyota Kaizen、工程改善敘事、bug fix postmortem
- **不適合**: 純概念解釋、純結果展示
- **源**: Toyota internal practice
- **InterSubMod 對應**: 既有 `/structured-tech-report` 13 段 = A3 + ADR + Postmortem hybrid

### DMAIC — Define / Measure / Analyze / Improve / Control

- **結構**: 5 階段 Six Sigma 流程
- **適合**: 流程改善（process improvement）、quality control、operational excellence
- **不適合**: 探索性研究、創新性問題

### PDCA — Plan / Do / Check / Act

- **結構**: 4 階段持續改善 cycle
- **適合**: 持續改善、實驗設計、iterative development
- **不適合**: 一次性決策
- **源**: W. Edwards Deming（基於 Walter Shewhart）

### OODA Loop — Observe / Orient / Decide / Act

- **結構**: 4 階段快速決策 loop
- **適合**: 高速決策、紅藍對抗、incident response、軍事戰術
- **不適合**: 慢思考、長期戰略
- **源**: John Boyd, US Air Force papers
- **核心**: 比對手更快跑完 OODA 即勝出

---

## §6 教學解釋（概念）

### Feynman Technique

- **結構**: (1) 選概念 → (2) 用 5 歲小孩語言重述 → (3) 找漏洞 → (4) 簡化 + 比喻 → (5) Repeat
- **適合**: 概念深度測試、學習自我評估、跨領域教學
- **不適合**: 高層 exec summary（過度簡化）
- **源**: Richard Feynman / James Gleick《Genius》
- **Template**: `templates/feynman_skeleton.md`

### Diátaxis — Tutorial / How-to / Reference / Explanation（2×2）

- **結構**: 2×2 framework
  - Tutorial（learning-oriented, hands-on）
  - How-to（problem-oriented, recipe）
  - Reference（information-oriented, lookup）
  - Explanation（understanding-oriented, discursive）
- **適合**: 技術文件、API doc、SDK、developer portal
- **不適合**: 一般敘事 / 故事
- **源**: Daniele Procida [diataxis.fr](https://diataxis.fr/)

### Bloom's Taxonomy

- **結構**: Remember → Understand → Apply → Analyze → Evaluate → Create
- **適合**: 教案、學習路徑、課程設計、test bank
- **不適合**: 1 次性溝通
- **源**: Benjamin Bloom 1956（2001 revised by Anderson & Krathwohl）

### KWL — Know / Want / Learned

- **結構**: K（已知）→ W（想知道）→ L（學到）
- **適合**: 工作坊、reading group、研究調查 entry point、自學日記
- **不適合**: 説服性簡報

### ELI5 — Explain Like I'm 5

- **結構**: 用 5 歲小孩能懂的語言 + 比喻
- **適合**: reddit 風科普、跨領域同步、技術 → 業務溝通
- **不適合**: 學術精確場合（會失嚴謹）
- **源**: Reddit /r/explainlikeimfive

### CPP — Concept / Procedure / Principle

- **結構**: 是什麼（Concept）→ 怎麼做（Procedure）→ 為什麼（Principle）
- **適合**: 教學講義、tutorial、手冊
- **不適合**: 説服性 keynote
- **源**: Instructional design theory

---

## §7 回饋對話 / Coaching

### SBI — Situation / Behavior / Impact

- **結構**: Situation（當時情境）→ Behavior（觀察到的行為）→ Impact（造成的影響）
- **適合**: 一對一回饋、performance review、難對話
- **不適合**: 公開表揚（過於結構化）
- **源**: Center for Creative Leadership (CCL)
- **Template**: `templates/sbi_skeleton.md`

### DESC — Describe / Express / Specify / Consequences

- **結構**: Describe（描述）→ Express（表達感受）→ Specify（具體要求）→ Consequences（後果）
- **適合**: 衝突管理、難對話、assertive communication
- **不適合**: 一般正面回饋
- **源**: Bower & Bower《Asserting Yourself》

### GROW — Goal / Reality / Options / Will

- **結構**: Goal（目標）→ Reality（現況）→ Options（選項）→ Will（行動意願）
- **適合**: 教練對話、mentoring、self-coaching
- **不適合**: 即時回饋、危機處理
- **源**: John Whitmore《Coaching for Performance》

### Pendleton's Model

- **結構**: 學員講好的 → 教練講好的 → 學員講可改 → 教練講可改
- **適合**: 醫學教育、臨床帶教、formative assessment
- **不適合**: 1 對 1 績效 review（先正後負模式已被 Radical Candor 批評為 inauthentic）

### Sandwich Feedback

- **結構**: 正面 → 改進 → 正面
- **適合**: 入門級 1-on-1（不熟對方時）
- **不適合**: 高信任關係（被批 inauthentic）

### Radical Candor — 2×2 (Care Personally × Challenge Directly)

- **結構**: 2 軸 4 象限
  - Care + Challenge = **Radical Candor**（理想）
  - Care + No Challenge = Ruinous Empathy
  - No Care + Challenge = Obnoxious Aggression
  - No Care + No Challenge = Manipulative Insincerity
- **適合**: 直接回饋、扁平團隊、Silicon Valley culture
- **不適合**: 階層嚴格組織、低信任場合
- **源**: Kim Scott《Radical Candor》(2017)

---

## §8 行銷 / 説服 / 銷售

### AIDA — Attention / Interest / Desire / Action

- **結構**: Attention（吸引注意）→ Interest（建立興趣）→ Desire（激發渴望）→ Action（呼籲行動）
- **適合**: 廣告 copy、landing page、email、direct response
- **不適合**: 學術論文、客觀報告
- **源**: E. St. Elmo Lewis (1898)
- **Template**: `templates/aida_skeleton.md`

### AIDCA — AIDA + Confidence

- **結構**: A → I → D → C（建立信心 / 信任）→ A
- **適合**: 高單價產品、B2B sales、複雜服務
- **不適合**: 衝動購買 / 快消品

### FAB — Feature / Advantage / Benefit

- **結構**: Feature（產品特性）→ Advantage（相對優勢）→ Benefit（給用戶帶來什麼）
- **適合**: sales pitch、產品介紹、技術 → 業務翻譯
- **不適合**: 純功能列表、規格書

### PAS — Problem / Agitate / Solution

- **結構**: Problem（痛點）→ Agitate（放大痛感）→ Solution（解決方案）
- **適合**: copywriting、痛點行銷、insurance / financial
- **不適合**: B2B technical sales（過於 emotional）

### 4P — Promise / Picture / Proof / Push

- **結構**: Promise（承諾）→ Picture（描繪場景）→ Proof（證據）→ Push（推動行動）
- **適合**: 直銷文案、long-form sales letter
- **不適合**: 短 banner ad
- **源**: Henry Hoke / Direct Mail copywriting

### BAB — Before / After / Bridge

- **結構**: Before（現況痛）→ After（理想世界）→ Bridge（你的方案是橋樑）
- **適合**: 轉型故事、產品改變敘事、變革倡議
- **不適合**: 客觀數據比較

---

## §9 決策評估

### WRAP — Widen options / Reality-test / Attain distance / Prepare to be wrong

- **結構**: Widen（拓寬選項，避開 narrow framing）→ Reality-test（測試假設）→ Attain distance（情緒距離）→ Prepare to be wrong（預備出錯）
- **適合**: 重大個人 / 事業決策、四象限思考
- **不適合**: 高速 OODA 場合
- **源**: Chip & Dan Heath《Decisive》(2013)
- **Template**: `templates/wrap_skeleton.md`

### Pre-Mortem

- **結構**: 假設 6 月後失敗 → 反推 risky assumption → 防範 action
- **適合**: 計畫風險檢核、project kickoff、大決策前
- **不適合**: 純樂觀 brainstorm
- **源**: Gary Klein, HBR 2007《Performing a Project Premortem》

### Cynefin — 5 域

- **結構**: Clear / Complicated / **Complex** / Chaotic / Disorder
  - **Clear**: 已知已知；best-practice
  - **Complicated**: 未知已知；expert + good-practice
  - **Complex**: 未知未知；probe-first
  - **Chaotic**: 無因果；act-sense-respond
  - **Disorder**: 無法分類
- **適合**: 因果性質判斷、危機分類、組織決策
- **不適合**: 簡單線性問題
- **源**: Dave Snowden, IBM Research (1999) → [thecynefin.co](https://thecynefin.co/)

### ADR — Architecture Decision Record

- **結構**: Context → Decision → Status → Consequences
- **適合**: 工程架構決策紀錄、技術選型、長期可追溯
- **不適合**: 短期 ad-hoc 決策
- **源**: Michael Nygard 2011

### DACI — Driver / Approver / Contributor / Informed

- **結構**: 4 角色 RACI 變體
- **適合**: 角色釐清、責任分配、cross-team project
- **不適合**: 1-2 人決策
- **源**: Intuit internal practice

### DECIDE — Define / Establish criteria / Consider alternatives / Identify best / Develop plan / Evaluate

- **結構**: 6 階段結構化決策
- **適合**: 醫療緊急、急救、生死決策
- **不適合**: 創意 / brainstorm

### Eisenhower Matrix — Urgent × Important（2×2）

- **結構**: 4 象限
  - Urgent + Important = Do now
  - Important + Not Urgent = Schedule
  - Urgent + Not Important = Delegate
  - Neither = Drop
- **適合**: 時間管理、優先序、todo list
- **不適合**: 戰略思考
- **源**: Dwight Eisenhower / Stephen Covey《7 Habits》

---

## §10 簡報架構 / 演講設計

### TED 18-min format

- **結構**: 1 idea worth spreading + narrative arc + emotional hook
- **適合**: TED / TEDx / 一般 keynote
- **不適合**: 技術 deep dive、產品 demo
- **時長**: 18 分鐘是 sweet spot（cognitive science backed）

### 3-Act Pitch（Duarte）

- **結構**: Beginning（What is）→ Middle（Call to Adventure + 反覆 What is ↔ What could be）→ End（New bliss）
- **適合**: 變革推銷、轉型 deck、product launch
- **不適合**: 數據 report
- **源**: Nancy Duarte《Resonate》

### Resonate / Slide:ology

- **結構**: Audience hero × Mentor presenter
- **適合**: 全本簡報設計、品牌大會
- **不適合**: 純技術文件
- **源**: Nancy Duarte《slide:ology》(2008) +《Resonate》(2010)

### Presentation Zen

- **結構**: 簡 → 清 → 留白 → Story（敘事）
- **適合**: 視覺優先簡報、Apple-style keynote
- **不適合**: 數據密集 dashboard
- **源**: Garr Reynolds《Presentation Zen》(2008)

### Lessig Method

- **結構**: 1 word / 1 image per slide × 100+ slides，快節奏切換
- **適合**: 快節奏 keynote、活力演講、議題倡議
- **不適合**: 教學深度
- **源**: Lawrence Lessig, Stanford

### Takahashi Method

- **結構**: 大字 / 1-2 詞 / slide
- **適合**: 日式極簡、focus 演講
- **不適合**: 紀錄性簡報
- **源**: Masayoshi Takahashi (Japanese Ruby community)

### PechaKucha — 20×20

- **結構**: 20 張 × 20 秒 = 6:40 固定
- **適合**: 設計 / 創意社群活動、networking event
- **不適合**: 深度技術分享

### PEEL — Point / Evidence / Explain / Link

- **結構**: Point（論點）→ Evidence（證據）→ Explain（解釋連結）→ Link（連到下段）
- **適合**: 學術 paragraph、essay、文章段落
- **不適合**: 口頭簡報

---

## §11 跨類整合套餐（hybrid framework — 既有 InterSubMod skill 對應）

這 7 個 hybrid framework 是既有 7 個報告類 skill 的核心 — 在 thin wrapper 化後仍是預設套餐。

### A3+ADR+Postmortem-hybrid（13 段）

- **結構**: TL;DR / 報告目的 / 系統背景 / 原本流程 / 問題描述 / 根本原因 (5 Whys) / 修改方向 (ADR Decision) / 修改內容 7.1 非工程 / 7.2 工程 / 新舊比較 / 驗證 Step→Verify / 影響範圍 / 風險限制 / 後續工作 / 結論
- **基底**: Toyota A3 + Michael Nygard ADR + Google SRE Postmortem + Daniele Procida Diátaxis
- **適合**: 單一工程 / pipeline / 方法學改動深度敘事
- **InterSubMod 對應**: `/structured-tech-report` thin wrapper 預設

### Multi-Thread-Narrative（17 段 4 主線）

- **結構**: W1 raw → W2 主線（進展 / 問題 / 求協助 / 探索 4 選 1）→ W3 內容 4 層 [F/O/I/U] → W4 4 桶分流 → W5 紅旗 → W6 教授問答 → W7 母稿 Layer 0-4
- **基底**: SCQA + STAR per item + Bloom's 4 層分類 + 教授視角預測
- **適合**: 每週多主題彙報、lab meeting、PI weekly
- **InterSubMod 對應**: `/weekly-report` thin wrapper 預設

### Audience-Scenario-Pitch（5 階段 P1-P5）

- **結構**: P1 audience + thesis → P2 outline → P3 section → P4 slide build → P5 speaker script
- **6 報告模板**:
  - improvement_report: Before → Problem → Cause → Solution → Verify → Impact
  - comparison_report: Setup → A → B → Side-by-side → Verdict
  - executive_summary: TL;DR → Top 3 → Risks → Asks
  - data_showcase: Hypothesis → Data → Stats → Caveats
  - concept_walkthrough: Why → Define → Mechanism → Examples → Boundary
  - academic_defense: Background → Question → Method → Result → Discussion → Future
- **適合**: PI / lab meeting / conference 簡報
- **InterSubMod 對應**: `/pptx-build` thin wrapper 預設

### Data-Showcase

- **結構**: Hypothesis → Data → Stats（含 Cohen ribbon + CI）→ Caveats
- **適合**: 單實驗結果分析、metric 解讀
- **InterSubMod 對應**: `/results-report` thin wrapper 預設

### Verdict-Pyramid

- **結構**: Verdict（最頂結論）→ Supporting evidence（Pyramid 結構）→ WRAP falsifier observable
- **適合**: 假説最終判定（POSITIVE / NEGATIVE / CONDITIONAL / NO-GO）、研究 cycle 收尾
- **InterSubMod 對應**: `/conclude-research` thin wrapper 預設

### AI-Session-Companion

- **結構**: Timeline → Key decisions → Provenance（commit hash + cycle_id + 生成時間）
- **適合**: AI 對話過程本身的決策紀錄
- **InterSubMod 對應**: `/report` thin wrapper 預設

### Hybrid-PI-Report (Hybrid HTML 4 層 L0-L3)

- **結構**: L0 Headline（一句 verdict + 4-6 焦點數字）→ L1 Top Findings（3-5 cards）→ L2 Evidence Cards（折疊 + tier badge）→ L3 Raw Data（DataTables / Plotly）
- **適合**: PI 終版 HTML、複雜 evidence 整合
- **InterSubMod 對應**: `/html-report-build` standalone mode

---

## §12 跨類補充（非核心但有用）

| Framework | 結構 | 適合 |
|-----------|------|------|
| **PESTLE** | Political / Economic / Social / Technological / Legal / Environmental | 外部環境分析 |
| **SWOT** | Strengths / Weaknesses / Opportunities / Threats | 內部 + 外部評估 |
| **First Principles** | 拆到最基本不可分原則 → 重建 | Elon Musk 風推理 |
| **Six Thinking Hats** | 6 帽子（白 / 紅 / 黑 / 黃 / 綠 / 藍）依序戴 | brainstorm 防偏見 |
| **Hourglass Writing** | Wide → Narrow → Wide | essay / paper |
| **CRABS** | Coherence / Relevance / Accuracy / Brevity / Structure | writing checklist |
| **5C Analysis** | Customer / Company / Competitor / Collaborator / Context | marketing situational |
| **OKR** | Objectives + Key Results | goal-setting |
| **STP** | Segmentation / Targeting / Positioning | marketing strategy |
| **Job-To-Be-Done** | When ___, I want ___, so I can ___ | product design |

---

## §13 場景 → 框架快速查表

完整對照見 `scenario_to_framework.md`；本表為入口 sample。

| 想表達 | 首選 | 備選 |
|--------|------|------|
| 報告新發現給 PI（≤5min） | SCQA + ABT | Assertion-Evidence + PREP |
| 回答教授追問（30s） | PREP | OREO + Rule of Three |
| 本週進度 | Multi-Thread-Narrative | SCQA + STAR per item |
| bug postmortem | A3+ADR+Postmortem-hybrid | Fishbone + 5 Whys |
| 架構決策 | ADR | Pre-Mortem + WRAP |
| 研究 keynote | Duarte Sparkline + Pixar Spine | Golden Circle + Assertion-Evidence |
| 論文 abstract | SCQA + Pyramid | BLUF + Inverted Pyramid |
| science 科普影片 | Pixar Spine + ABT | Hero's Journey |
| 產品 demo / pitch | Golden Circle + FAB | AIDA + BAB |
| slack 同步 | BLUF | TL;DR + 5W1H |
| 行為面試 | STAR | CAR / SOAR |
| 1-on-1 回饋 | SBI | Radical Candor |
| 救急決策（高壓） | OODA | DECIDE |
| 不確定方向選擇 | WRAP + Pre-Mortem + Cynefin | Bland 2×2 (Assumption Map) |
| 教 1 個概念 | Feynman + ELI5 | Concept-Procedure-Principle |
| doc / SDK 文件 | Diátaxis | Inverted Pyramid |
| 改善企劃 | DMAIC / PDCA | A3 |
| issue tree / case interview | MECE + Pyramid | 5W1H |

---

## §14 業界源頭精要

完整 URL + ISBN + 一句引用見 `framework_business_sources.md`。

| 核心源 | 代表書 / URL | 框架貢獻 |
|-------|------------|--------|
| **Barbara Minto** | 《The Pyramid Principle》(2020) | SCQA + Pyramid + MECE |
| **Randy Olson** | 《Narrative Is Everything: The ABT Framework》(2019) | ABT |
| **Joseph Campbell** | 《The Hero with a Thousand Faces》(1949) | 英雄之旅 |
| **Syd Field** | 《Screenplay》(1979) | 三幕劇 |
| **Nancy Duarte** | 《Resonate》(2010) +《slide:ology》(2008) | Sparkline + 3-Act |
| **Simon Sinek** | 《Start with Why》(2009) + TED 2009 | Golden Circle |
| **Michael Alley** | 《The Craft of Scientific Presentations》 | Assertion-Evidence |
| **Toyota / 大野耐一** | TPS internal | A3 + 5 Whys |
| **石川馨** | 1968 工程文獻 | Fishbone |
| **John Boyd** | USAF papers | OODA Loop |
| **Daniele Procida** | [diataxis.fr](https://diataxis.fr/) | Diátaxis |
| **Benjamin Bloom** | 1956《Taxonomy of Educational Objectives》 | Bloom's |
| **CCL** | CCL papers | SBI |
| **Heath brothers** | 《Decisive》(2013) | WRAP |
| **Gary Klein** | HBR Sep 2007 | Pre-Mortem |
| **Dave Snowden** | [cynefin.io](https://thecynefin.co/) | Cynefin |
| **Kenn Adams** | improv community | Pixar Story Spine |
| **Alan Monroe** | Purdue 1935 | Monroe's Motivated Sequence |
| **Michael Nygard** | 2011 ThoughtWorks | ADR |
| **Kim Scott** | 《Radical Candor》(2017) | Radical Candor 2×2 |
| **Lewis** | 1898 marketing | AIDA |
| **Garr Reynolds** | 《Presentation Zen》(2008) | Zen |
| **Edward Tufte** | 《The Visual Display》(1983/2001) | Data-ink ratio |
| **Robin Williams** | 《The Non-Designer's Design Book》(1993) | CRAP |

---

## §15 Catalog 維護準則

1. **添加新 framework**: 必填「結構 / 適合 / 不適合 / 源」4 項；缺項拒收
2. **更新業界源 URL**: 同步更新 `framework_business_sources.md`
3. **drift check**: 每月 audit — `scenario_to_framework.md` 推薦項全部能 grep 到本檔
4. **template 追加**: 新 framework 收進來必同時新增 `templates/<name>_skeleton.md`
5. **去重**: 已 superseded 框架（如 sandwich feedback → radical candor）標 "deprecated; see X"

---

## §16 不在 catalog 的（intentionally excluded）

| Framework | 為什麼不收 |
|-----------|----------|
| TDD test cycle | 工程方法論，非敘述框架 |
| MoSCoW | 優先序方法，非敘述 |
| INVEST | user story 規範，非敘述 |
| RACI（含 R 而非 D） | DACI 已涵蓋 |
| Lean Canvas | 商業模型畫布，非敘述 |
| Story Mapping | 產品規劃，非敘述 |
| 5E（教育）| 部分被 Bloom 涵蓋 |

如需擴充上述為敘述用途，重寫成「結構 / 適合 / 不適合 / 源」4 項後可加入 §12 跨類補充。
