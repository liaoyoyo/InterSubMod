# ELI5 Skeleton

> Reddit /r/explainlikeimfive (2011)

```markdown
---
framework: ELI5
when: <reddit 風科普 / 跨領域同步 / 技術 → 業務溝通>
---

# <Concept>

## 1 句話講給 5 歲小孩

<避免術語；用具體日常物比喻>

## 為什麼這比喻 work

<解釋比喻對應到哪個 mechanism>

## 比喻的局限（誠實揭露）

<這個比喻在哪裡失效？什麼時候不能用？>

---

範例:

**Concept**: Phasing

**ELI5**: 「想像你打牌時撲克牌兩副混在一起 — phasing 就是把它們分回兩堆。」

**為什麼 work**: 兩條染色體 = 兩副牌；每個位點 ATCG = 牌面；phasing algorithm 看 read 共享哪些位點來判斷屬於哪副。

**局限**: 真實 phasing 有 noise（牌面有時看錯）；ONT long read 像「能一次看 50 張牌」（強），short read 像「只看 1-2 張」（弱）。所以這個比喻適合解釋 long-read phasing；short-read 需另一個比喻。

---

Framework: ELI5 (Reddit /r/explainlikeimfive 2011)
```
