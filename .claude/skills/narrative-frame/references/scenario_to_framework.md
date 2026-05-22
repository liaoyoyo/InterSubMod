# Scenario → Framework 對照表

> `/narrative-frame` N1 場景識別後從本表查推薦框架。**5W 維度 + 自然語問句**雙索引，用戶不必背 framework 名。

---

## §1 5W → Framework 快速對照矩陣

### Who（受眾）→ formality / 術語層

| Who | 主框架傾向 | 避免 |
|-----|---------|-----|
| **PI / 教授** | SCQA / Pyramid / Assertion-Evidence | ELI5（過簡）、Sandwich（不真誠）|
| **同儕** | ABT / Pixar Spine / Multi-Thread-Narrative | 純技術 ADR（除非工程同儕） |
| **大眾 / 科普** | Pixar Spine / Hero's Journey / ELI5 / ABT | A3（太工程）、Diátaxis（過於分類）|
| **自己 / 草稿** | KWL / Feynman / 5W1H | TED 18min（過度設計） |
| **學生 / 教學** | Feynman / Bloom / CPP / Diátaxis | BLUF（無 buildup）|
| **投資人 / 評審** | Golden Circle / SCQA / FAB | Feynman（過簡）|
| **下屬 / mentee** | SBI / GROW / Radical Candor | DACI（過度結構化） |
| **客戶 / sales** | AIDA / FAB / PAS / BAB | A3（過工程） |

### Why（目的）→ 框架類別

| Why | 主框架類 | 例 |
|-----|--------|----|
| **説服** | §4 故事敘事 + §8 行銷 | Duarte Sparkline / AIDA / Golden Circle |
| **解釋** | §6 教學解釋 | Feynman / ELI5 / Diátaxis / CPP |
| **報告進度** | §1 結論先行 + §2 經驗敘事 | SCQA / STAR / Multi-Thread-Narrative |
| **探索 / 腦力激盪** | §5 問題分析 + §9 決策 | 5W1H / 5 Whys / WRAP / Cynefin |
| **紀錄** | §1 結論先行 | BLUF / Inverted Pyramid / ADR / A3 |
| **答辯 / Q&A** | §3 短時口頭 | PREP / OREO / ABT |
| **教學** | §6 教學解釋 | Feynman / Bloom / KWL |
| **比較 / decision** | §9 決策評估 | WRAP / Pre-Mortem / DECIDE |
| **回饋** | §7 對話 / Coaching | SBI / DESC / GROW |
| **行銷 / pitch** | §8 行銷 + §4 故事 | AIDA / FAB / Golden Circle / Pixar |

### What（內容類型）→ specific framework

| What | 推薦 |
|------|-----|
| **新發現** | SCQA + STAR per item / ABT |
| **進度更新** | Multi-Thread-Narrative (W1-W7) / Inverted Pyramid |
| **概念解釋** | Feynman + ELI5 / CPP / Diátaxis |
| **A vs B 比較** | Comparison Report (6 模板之一) / FAB / MECE |
| **故事 / case** | Pixar Spine / Hero's Journey / Freytag |
| **決策** | ADR / WRAP / DECIDE / Pre-Mortem |
| **數據 / 結果** | Data-Showcase / Assertion-Evidence / BLUF |
| **流程** | DMAIC / PDCA / A3 |
| **bug / 失敗** | A3 + 5 Whys + Postmortem / Fishbone |

### When（時長）→ 細節粒度

| When | 推薦 framework 類 | 段數上限 |
|------|--------------|--------|
| **30s（pitch / Slack）** | ABT / Elevator / BLUF / PREP | 1 段（≤3 句） |
| **2-5min（lab quick）** | SCQA / STAR / PREP × 3 | 3-5 段 |
| **18min（TED / lab meeting）** | Pixar Spine / Duarte Sparkline / Multi-Thread-Narrative | 5-7 段 |
| **45-60min（lecture / academic defense）** | Hero's Journey / 三幕劇 / academic_defense template | 7-13 段 |
| **紙本（paper / report）** | 13 段 hybrid / Pyramid 深巢 | 13+ 段 |

### How（輸出形式）→ rendering style

| How | 推薦 framework | 注意 |
|-----|-------------|-----|
| **口頭** | PREP / ABT / Pixar Spine / TED arc | 短句 + 三點法則 |
| **slide deck** | Assertion-Evidence + Audience-Scenario-Pitch | 標題 = 結論句 |
| **.md 報告** | 13 段 hybrid / Multi-Thread / SCQA | 含 frontmatter |
| **HTML standalone** | Hybrid-PI-Report 4 層 (L0-L3) | sticky TOC + 折疊 |
| **paper / academic** | SCQA + Pyramid / academic_defense | Pre-reg 對照表 |
| **email** | BLUF / PREP | 1 段結論 + 細節 |
| **Slack** | BLUF + Inverted Pyramid | 1 行 + 細節 |
| **對話 inline** | Tier 2 首行宣告 framework | 「用 PREP：」短句 |

---

## §2 自然語問句 → Framework（用戶不必背名）

「**我想 / 我要 / 該怎麼...**」→ 對應框架

| 自然語句 | 對應框架 |
|---------|--------|
| 「整理本週進度」 | Multi-Thread-Narrative (= /weekly-report) |
| 「報告新發現」 | SCQA + ABT |
| 「説服老闆批准 X」 | Golden Circle + Pre-Mortem + WRAP |
| 「説服 PI 換方向」 | Duarte Sparkline + WRAP |
| 「教新人 Y 概念」 | Feynman + ELI5 + CPP |
| 「整理 commit fix」 | A3 + 5 Whys + Postmortem (= /structured-tech-report) |
| 「答辯這個 Q」 | PREP |
| 「pitch 30 秒」 | ABT / Elevator |
| 「做 PI 簡報」 | Audience-Scenario-Pitch (= /pptx-build) |
| 「整合 3 份報告」 | Pyramid + MECE + narrative-organizer agent |
| 「對比 A vs B」 | Comparison Report / FAB / MECE |
| 「找 bug 根因」 | 5 Whys + Fishbone + A3 |
| 「決定要不要做 X」 | WRAP + Pre-Mortem + Cynefin |
| 「給同事回饋 Z」 | SBI / DESC / Radical Candor |
| 「寫實驗結果」 | Data-Showcase (= /results-report) |
| 「研究收尾 verdict」 | Verdict-Pyramid (= /conclude-research) |
| 「寫 paper abstract」 | SCQA + Pyramid + BLUF |
| 「TED-style keynote」 | Pixar Spine + Duarte Sparkline |
| 「Slack 同步」 | BLUF |
| 「email 給教授」 | BLUF + PREP |
| 「postmortem」 | A3 + Inverted Pyramid + Blameless 5-段 |
| 「做架構決策紀錄」 | ADR |
| 「分配責任 / who's accountable」 | DACI / RACI |
| 「優先序」 | Eisenhower 2×2 / Bland 2×2 |
| 「教這個工具 / library」 | Diátaxis (4 quadrant) |
| 「腦力激盪 / explore options」 | Cynefin → Complex 域 PROBE / WRAP Widen |
| 「説故事 / impact story」 | Pixar Spine / Hero's Journey |
| 「給管理層摘要」 | Executive Summary template / Pyramid |
| 「比較選項」 | MECE + Pyramid / DECIDE |
| 「找痛點 + 賣方案」 | PAS + AIDA |
| 「救急決策 / 高壓 / incident response」 | OODA Loop (Boyd) |
| 「持續改善 / quality cycle」 | PDCA / DMAIC |

---

## §3 場景特例（直接路由 — 不需 5W 全跑）

對應 SKILL.md §4 場景特例。**用戶説這些觸發詞 → 直接路由 thin wrapper skill**：

| 觸發詞 | 路由到 | 預設 framework |
|--------|-----|--------------|
| 「整理本週 / 週報 / lab meeting」 | `/weekly-report` thin wrapper | Multi-Thread-Narrative |
| 「整理 commit / fix / 修補 / bug postmortem」 | `/structured-tech-report` thin wrapper | A3+ADR+Postmortem-hybrid |
| 「整理實驗結果 / metric 分析」 | `/results-report` thin wrapper | Data-Showcase |
| 「整理 PI 簡報 / pitch / 教授級簡報」 | `/pptx-build` thin wrapper | Audience-Scenario-Pitch |
| 「整理 AI session」 | `/report` thin wrapper | AI-Session-Companion |
| 「研究收尾 / 假説 verdict / NO-GO」 | `/conclude-research` thin wrapper | Verdict-Pyramid |
| 「整理 PI 終版 HTML」 | `/html-report-build` standalone | Hybrid-PI-Report (L0-L3) |

---

## §4 Cynefin Domain 強制路由

對應 SKILL.md §5。**框架選擇前必先過 Cynefin domain gate**：

| Cynefin Domain | 允許 framework 類 | 禁止 framework 類 |
|---------------|----------------|----------------|
| **Clear** | 全可用 — best-practice 直套（SCQA / STAR / PREP / FAB） | 無 |
| **Complicated** | 多數可用 — 建議 §5 problem analysis（5 Whys / Fishbone / A3）+ §6 教學 | 純故事敘事 §4（過度簡化）|
| **Complex** | **強制 probe-first** — §9 決策（WRAP / Cynefin 自己 / Pre-Mortem）+ §6 教學（Feynman 反問） | **所有 deterministic best-practice framework**（SCQA / Pyramid / STAR — 假設因果未確立時誤導）|
| **Chaotic** | 先穩定（act-sense-respond） | 任何敘事 framework |

**判斷 Cynefin domain**：問「相同行動是否曾重複產生**可預測**結果？」
- Yes → Clear / Complicated
- Maybe → Complex
- No → Chaotic

借用 `/pre-decision-audit §0` 邏輯。

---

## §5 InterSubMod 專案特定場景

研究專案內常見場景 + 對應 framework：

| 研究階段 | 場景 | 推薦 framework |
|---------|------|--------------|
| **P0 cycle init** | 新 cycle 規格説明 | SCQA + Pre-reg 對照表 |
| **P1 PLAN** | 假説 ORIENT + 驗證設計 | Multi-Thread / WRAP |
| **P3 PILOT** | 單樣本特徵觀察 | Data-Showcase + Assertion-Evidence |
| **P4 multi-sample** | 7 樣本一致性 | Multi-Thread + cross-sample comparison table |
| **P5 EVALUATE** | tier 升級 retraction risk | Verdict-Pyramid + WRAP falsifier |
| **P6 CONCLUDE** | 最終 verdict 宣告 | Verdict-Pyramid + Postmortem (NEGATIVE) |
| **PI 1-on-1 (5min)** | 上週進度 + 本週計畫 | Multi-Thread + SCQA |
| **Lab meeting (15min)** | 詳細 case 講解 | Pixar Spine + Assertion-Evidence |
| **Conference talk (18min)** | 突破性結論 | Duarte Sparkline + Golden Circle |
| **論文 abstract** | 200-word summary | SCQA + Pyramid + BLUF |
| **論文 results section** | 結果展示 | Data-Showcase + Assertion-Evidence |
| **論文 discussion** | 機制 + limitation | 5 Whys + WRAP falsifier + Hourglass |
| **bug postmortem** | 失敗事件根因 | A3 + 5 Whys + Fishbone + Blameless 5-段 |
| **C++ refactor 説明** | 工程改動敘事 | A3+ADR+Postmortem-hybrid (= 13 段) |
| **新 skill / hook spec** | 設計決策紀錄 | ADR + implementation-notes 4 sections |
| **跨樣本 framework drift** | 5 個樣本不一致 root cause | Cynefin (Complex) → WRAP + Pre-Mortem |
| **decision audit (前)** | 進入 cycle 前審計 | `/pre-decision-audit` 7 sections |

---

## §6 框架組合（hybrid usage）

主框架 + 內嵌 sub-framework — 大型報告常用組合：

| 主框架 | 內嵌 sub | 適用 |
|--------|--------|------|
| **SCQA** + 每段 **STAR** | section-by-section behavioral evidence | PI 報告新發現 |
| **Pixar Spine** + 內嵌 **ABT** in each "because of that" | 故事 + 因果 | 科普 keynote |
| **Multi-Thread-Narrative** + per-thread **PREP** | 週報快速回答型 | weekly + Q&A 預備 |
| **13 段 hybrid** + §6 用 **5 Whys** + §7 用 **ADR** | 完整工程敘事 | structured-tech-report |
| **TED 18min** 外殼 + 內裏 **Duarte Sparkline + Pixar Spine** | 説服性 keynote | conference talk |
| **Assertion-Evidence** slide title + 內 **FAB** body | 産品 / 技術簡報 | pitch deck |
| **WRAP** 主流程 + 每選項 **Pre-Mortem** | 決策評估 | 大方向 pivot |
| **Diátaxis** + 每 quadrant **CPP** | tutorial + 概念解釋 | SDK doc |
| **Inverted Pyramid** + 結尾 **BLUF** | 新聞 / blog | status update |
| **A3** + §Analysis 用 **Fishbone + 5 Whys** + §Proposal 用 **ADR** | postmortem | bug fix |

---

## §7 反 pattern（不推薦組合）

| 不推薦組合 | 原因 |
|---------|------|
| TED 18min keynote + 全程 BLUF | 缺敘事弧；觀眾累 |
| Slack 同步 + Pixar Spine | 過度設計 30s 訊息 |
| Bug postmortem + AIDA | 行銷框架不適合根因 |
| 學術 paper + Pure Hero's Journey | 個人化敘事 vs 客觀分析衝突 |
| Q&A 答辯 + 三幕劇 | 答 30s 不需 setup-confrontation-resolution |
| 概念教學 + BLUF | 缺 buildup 學習者不懂 why care |
| 行為面試 + Pyramid | 行為題期望具體 case 而非結論先行 |
| 配對 SBI + Sandwich | 兩者衝突（前者直接 / 後者迂迴） |
| Cynefin Complex + SCQA | Complex 域因果未確立，不可套 deterministic framework |

---

## §8 N2 推薦邏輯（給 AI 內部用）

當 AI 在 N2 推薦框架時，按以下順序：

1. **檢查觸發詞 / 自然語句** — 命中 §2 / §3 → 直接路由
2. **解析 5W 維度** — Who / Why / What / When / How
3. **跑 Cynefin gate** — Complex 強制 PROBE，否則進 4
4. **多維度交叉查詢** — §1 + §6 hybrid 組合
5. **產出格式**：

```
[narrative-frame N2]
主框架: <name>
  - 結構: <一句講完>
  - 為什麼選: <對應 5W 哪幾維>
  - 業界源: <author + 書/年>

備框架 1: <name>
  - 差異: <何時切換更好>

備框架 2: <name>
  - 差異: <同上>

不選擇:
- <name>: <為什麼不適合（針對 5W）>
- <name>: <為什麼不適合>
```

6. **暫停 C2** — 等用戶 ack / override

---

## §9 維護準則

- 新增 framework 進 catalog → 必同步本檔加 row（至少 §2 自然語問句 1 條）
- 場景對應改變（如 weekly-report W2 變動）→ 更新 §3 路由與 §5 InterSubMod 場景
- Cynefin domain 路由（§4）跟 `/pre-decision-audit §0` 同步，single SoT
- §7 反 pattern 持續累積（每月 audit 一次）
