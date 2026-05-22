# Pendleton's Model Skeleton

> David Pendleton《The Consultation》(1984) — 教學回饋

```markdown
---
framework: Pendleton's Model
when: <醫學教育 / 臨床帶教 / formative assessment>
---

# 回饋對話結構

## 1. 學員自評: 哪裡好？

「<學員>，你覺得這次 X 做得好的地方是什麼？」

## 2. 教練補充: 還有哪些好的？

「我觀察到你還做了 Y 也很好，因為 ...」

## 3. 學員自評: 哪裡可改進？

「你覺得哪裡可以做得更好？」

## 4. 教練補充: 還有哪些可改進？

「另外我建議下次可以試 Z，因為 ...」

## 5. 結論 + 行動點

「整體來説今天的重點是 W；下次重點放 ...」

---

關鍵原則:
- **學員先講，教練後講** — 保護心理安全
- **「好的」先講，「可改進」後講** — 不是 sandwich（更直接）
- **避免 evaluative language**「你做得不夠」→「下次可試 X」

---

範例（C++ code review）:

1. 學員: 「我覺得 V6 的 hp tag 邏輯重構得比 V5 清楚」
2. 教練: 「對；另外你加的 per-region monitoring 很好，避免下次又盲」
3. 學員: 「我覺得 commit message 寫得不夠細」
4. 教練: 「對；另外 unit test 可以再加 chr8 edge case 覆蓋」
5. 結論: 「整體重點是 metric 設計 + 完整 test；下次 C++ change PR 前先跑 per-region」

---

Framework: Pendleton's Model (Pendleton《The Consultation》1984)
```
