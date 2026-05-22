# MECE Skeleton

> Barbara Minto / McKinsey — Mutually Exclusive, Collectively Exhaustive

```markdown
---
framework: MECE
when: <issue tree / 分類分析 / case interview / ML feature 設計>
---

# <Title — 分析主題>

## 維度（必互斥 + 共同窮盡）

維度選擇: <選一個維度作分類軸 — e.g., 時間 / 對象 / 流程 / 地理>

| Category | Definition | Examples | 不屬於此類 |
|----------|------------|----------|-----------|
| **A** | <定義 1> | <e.g., ...> | <反例> |
| **B** | <定義 2> | <e.g., ...> | <反例> |
| **C** | <定義 3> | <e.g., ...> | <反例> |

## MECE 自審

- ✓ **互斥**: 任一物件只屬於 1 個 category？驗證: 列 5 個樣本，每個只 fit 1 個
- ✓ **窮盡**: A+B+C 覆蓋所有可能？驗證: 想 5 個邊緣案例，全部能歸類
- ⚠ 若有「其他」category → 通常表示分類軸未對；重選維度

---

Framework: MECE (Barbara Minto《Pyramid Principle》1978, McKinsey)
```
