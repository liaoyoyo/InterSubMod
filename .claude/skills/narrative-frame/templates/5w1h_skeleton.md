# 5W1H Skeleton

> 古希臘修辭學 → 新聞學 standard

```markdown
---
framework: 5W1H
when: <事件 brief / issue 描述 / bug report / 新聞 lead>
---

# <Title>

**Who**: <主角 / 涉事者>
**What**: <發生了什麼 / 動作>
**When**: <時間 / 日期>
**Where**: <地點 / 系統 / 模組>
**Why**: <原因 / 動機>
**How**: <方式 / 機制>

---

範例（bug report）:

**Who**: V3F binary (longphase-to-mod v1.7-I)
**What**: chr19 752 reads 標 hp=33 (非 1/2)
**When**: 2026-Q1 首次觀察；2026-05-09 V6 修
**Where**: longphase HaplotagProcess.cpp judgeHaplotype:533
**Why**: priority bug — countMap[HP2]=3 邏輯 + assignment 1.77× 上游偏置
**How**: V6 同時修 assignment + priority 兩層，跑 paired pileup 驗證

---

Framework: 5W1H (古希臘修辭學 / journalism standard)
```
