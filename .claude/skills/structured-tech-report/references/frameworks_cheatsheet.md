# 框架速查表（structured-tech-report 引用）

13 段範本背後的方法學來源。寫到對應段時可回查本檔挑形式。

---

## A. 問題分析類

### Toyota A3 Thinking

- **出處**：Sobek & Smalley, *Understanding A3 Thinking* (2008); Lean Enterprise Institute
- **結構**：Background → Current Condition → Goal → Cause Analysis → Countermeasures → Plan → Follow-up
- **適用**：工程缺陷修復、流程改進；本 skill 13 段最貼合的母版
- **連結**：https://sloanreview.mit.edu/article/toyotas-secret-the-a3-report/

### 5 Whys

- **出處**：Sakichi Toyoda → Taiichi Ohno (Toyota Production System)
- **結構**：連問五次「為什麼」直到根因
- **適用**：單因果鏈深挖；§5 預設方法
- **與 fishbone 搭配**：先 fishbone 列分支 → 再對主分支跑 5 Whys

### Ishikawa Fishbone (RCA)

- **出處**：Kaoru Ishikawa, 1960s Kawasaki
- **結構**：6M（Man / Machine / Material / Method / Measurement / Environment）多因素分支
- **適用**：多因素問題；§5 替代 5 Whys 的選項

### Google SRE Postmortem

- **出處**：sre.google/sre-book
- **結構**：Summary → Impact → Root Causes → Trigger → Resolution → Detection → Action Items → Lessons (well/badly/lucky) → Timeline
- **適用**：事故覆盤（blameless）；§4 + §5 + §12 可參考
- **連結**：https://sre.google/sre-book/postmortem-culture/

### Atlassian Postmortem

- **出處**：atlassian.com/incident-management
- **結構**：Summary → Leadup → Fault → Impact → Detection → Response → Recovery → 5 Whys RCA → Timeline
- **適用**：ITSM 場景；本 skill 不主推但可借段落順序

---

## B. 決策／架構類

### ADR (Architecture Decision Record)

- **出處**：Michael Nygard, 2011 (Cognitect)
- **結構**：Title → Status → Context → Decision → Consequences
- **適用**：架構決策留痕；§6 + §10 + §11 對應
- **連結**：https://www.cognitect.com/blog/2011/11/15/documenting-architecture-decisions
- **本專案約定**：archi 決策應另落 `InterSubMod/docs/decisions/YYYY/MM/`，本 skill 在 §6 摘要引用

### MADR / Y-statement

- **出處**：adr.github.io
- **MADR 結構**：Context → Decision Drivers → Considered Options → Decision Outcome → Pros/Cons of Each
- **Y-statement**：「In context of X, facing Y, we chose Z to achieve W, accepting V」
- **適用**：§6 候選方案表的單句摘要

### Rust RFC / Python PEP

- **結構**：Summary → Motivation → Guide-level → Reference-level → Drawbacks → Alternatives → Unresolved
- **適用**：跨版本提案；§6 + §11 可借「Drawbacks」與「Alternatives」段名

---

## C. 變更敘事類

### Before/After/Why

- **出處**：GitHub / Conventional Commits 社群慣例（無單一權威）
- **結構**：前態 → 後態 → 變更動機
- **適用**：§3（Before）與 §8（新舊比較）對應；commit message body 通用

### STAR

- **出處**：DDI (Development Dimensions International), 1974
- **結構**：Situation → Task → Action → Result
- **適用**：個人貢獻敘事；**不適合**根因類技術報告（缺 cause/verify）；本 skill 13 段不採用 STAR
- **連結**：https://en.wikipedia.org/wiki/Situation,_task,_action,_result

---

## D. 受眾分層類

### Diátaxis

- **出處**：Daniele Procida, https://diataxis.fr/
- **結構**：Tutorial / How-to / Reference / Explanation 四象限
- **適用**：文件**整體**架構，非單篇報告
- **本 skill 內對應**：§7.1 ≈ Explanation（給外部讀者）；§7.2 ≈ Reference（給工程師查找）

### Amazon 6-pager

- **出處**：Bezos 2004 內部 memo
- **結構**：Intro → Goals → Tenets → State of Business → Lessons Learned → Strategic Priorities
- **適用**：高層決策、靜默會議閱讀；本 skill 借「TL;DR」（§0）與「Strategic」（§13）

---

## E. 規範類（介面／規格層）

### IETF RFC 2119

- **結構**：MUST/SHOULD/MAY 規範語；前言 → Terminology → Spec → Security → IANA
- **適用**：協定／介面規格層；本 skill 不直接套用，但 `data_specs/` 類規範文件可借

---

## F. 框架選用快速決策表

| 報告類別 | §5 用法 | §6 用法 | 補強段 |
|---------|--------|---------|--------|
| bug fix | 5 Whys | A/B/C 候選表 | §9 Step→Verify、§12 SRE Action Items |
| pipeline 變更 | 5 Whys + Ishikawa | A/B/C 候選表 + Y-statement | §8 Before/After diff、§9 7-sample 一致性 |
| 方法學改進 | Ishikawa（多因素） | MADR Decision Drivers | §9 三關 confound guard、§11 限制邊界 |
| 架構決策 | A3 Current Condition | ADR Context+Decision+Consequences | §10 + §11 合併為「Consequences」 |
| 事故覆盤 | SRE Postmortem | — | §4 + §5 + §12 + Lessons |

---

## G. 推薦延伸閱讀

1. [Toyota's Secret: The A3 Report (MIT SMR)](https://sloanreview.mit.edu/article/toyotas-secret-the-a3-report/) — 與 13 段最對齊
2. [Nygard ADR original (Cognitect)](https://www.cognitect.com/blog/2011/11/15/documenting-architecture-decisions) — ADR 起源
3. [Google SRE Book — Postmortem Culture](https://sre.google/sre-book/postmortem-culture/) — 事後覆盤黃金標準
4. [Diátaxis — Start Here](https://diataxis.fr/start-here/) — 受眾分層
5. [Atlassian Postmortem Templates](https://www.atlassian.com/incident-management/postmortem/templates) — 表格化骨架
