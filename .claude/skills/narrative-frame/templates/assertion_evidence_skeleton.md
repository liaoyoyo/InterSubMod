# Assertion-Evidence Skeleton

> Michael Alley《The Craft of Scientific Presentations》, Penn State

```markdown
---
framework: Assertion-Evidence (slide-level)
when: <科學簡報 / 學術 conf / PI 報告 / academic defense>
---

# Slide Design Rule

Each slide:
- **Title** = a complete sentence asserting the conclusion (NOT just topic label)
- **Body** = visual evidence supporting the assertion

❌ Bad: `Results`
✅ Good: `V6 修 100% chr19 priority bug victims (752/752 reads)`

❌ Bad: `Methods`
✅ Good: `Cross-sample LOSO with 5-fold CV shows HCC1395 +0.022 F1`

❌ Bad: `Discussion`
✅ Good: `caller_af direction-inconsistency explains 4/5 sample 卡住`

---

## Slide Skeleton

```html
<article class="slide-canvas">
  <h2>{Assertion — 一句完整斷言}</h2>
  <figure>
    <img src="evidence.png" alt="...">
    <figcaption>{Supporting data + source}</figcaption>
  </figure>
  <p>{Optional 1-2 sentence elaboration}</p>
</article>
```

---

業界驗證: Kuwait Univ 2025 實驗 108 學生 — Assertion-Evidence slides 顯著降低 cognitive load (vs traditional topic-bullet slides)

---

Framework: Assertion-Evidence (Michael Alley《Craft of Scientific Presentations》Penn State)
業界引用: 「Slide title = a complete sentence asserting the conclusion; body = visual evidence」
```
