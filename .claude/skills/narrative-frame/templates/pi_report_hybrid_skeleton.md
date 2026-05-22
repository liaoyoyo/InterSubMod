# Hybrid-PI-Report Skeleton (= 既有 html-report-build standalone)

> 4 層 L0-L3 + DataTables + Plotly + Alpine + MVP.css

```markdown
---
framework: Hybrid-PI-Report (4 層 L0-L3)
when: <複雜 PI 報告 / 大量 evidence + 終版 / PI 反覆讀>
對應既有 skill: /html-report-build (standalone mode)
output: HTML standalone (.standalone.html)
---

# 4 層架構

## L0: Headline（最頂 — 5 秒讀完）

- 一句 verdict（≤30 字）
- 4-6 焦點數字（stat-grid）
- Sticky TOC（左欄）

## L1: Top Findings（3-5 cards open by default）

每 card:
- Assertion-Evidence title (結論句)
- 1 大圖 (PNG or inline SVG)
- 3-5 bullet supporting
- 1-2 line caveat

## L2: Evidence Cards（折疊展開）

每 card:
- Tier badge ⭐⭐⭐ (L1-L5 from /scientific-rigor §2)
- Cohen ribbon + CI 欄
- DAG mermaid diagram (if causal claim)
- Pre-reg 對照表 (if validated)

## L3: Raw Data（最底 — 給 power user）

- DataTables interactive table（sortable / searchable / 1000+ rows OK）
- Plotly interactive plot (zoom / hover tooltip)
- Alpine.js filters
- Provenance footer (commit hash + cycle_id + 生成時間)

---

# 技術 stack

| Layer | Lib | CDN |
|-------|-----|-----|
| Markup | MVP.css | mvp.css |
| Table | DataTables 1.13 | jsDelivr |
| Plot | Plotly 2.27 | plot.ly |
| Reactive | Alpine 3.13 | alpinejs.dev |
| CSS scope | Inline `<style>` design_tokens.css | local |

---

# 6-Taboo Audit (output 前必跑)

1. 多重 gradient
2. Glass morphism
3. 多重 indigo
4. Emoji 開頭 h1/h2/h3
5. text-shadow
6. box-shadow 0 0 glow

詳見 InterSubMod/.claude/skills/html-report-build/SKILL.md §6-Taboo Audit

---

# 嚴謹度繼承（/scientific-rigor）

standalone PI 終版必繼承:
- §2: claim card 必標 ⭐⭐⭐⭐⭐ tier badge
- §3: metric table 必含 Cohen ribbon + CI
- §4: standalone 必嵌 SVG DAG（mermaid CDN 或 inline SVG）
- §7: validated 報告必含 Pre-registration 對照表
- §8.4: footer 必有 commit hash + cycle_id + 生成時間

---

Framework: Hybrid-PI-Report (InterSubMod html-report-build standalone)
對應既有 skill: /html-report-build (standalone mode)
```
