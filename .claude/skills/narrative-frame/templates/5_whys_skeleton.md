# 5 Whys Skeleton

> 大野耐一 / Toyota Production System

```markdown
---
framework: 5 Whys
when: <根因分析 / postmortem / Toyota gemba / debug>
---

# <Problem Statement>

<明確 problem — 1 句>

## 5 Whys

- **Why 1**: 為什麼 <problem>？
  → Because <answer 1 — 表象原因>

- **Why 2**: 為什麼 <answer 1>？
  → Because <answer 2 — 深一層>

- **Why 3**: 為什麼 <answer 2>？
  → Because <answer 3>

- **Why 4**: 為什麼 <answer 3>？
  → Because <answer 4>

- **Why 5**: 為什麼 <answer 4>？
  → **Root cause**: <answer 5>

## Action

<針對 root cause（不是表象）的修正>

---

範例:

Problem: V1.7-I 修正後仍見 17.3× HP1 偏置

- Why 1: 為什麼仍偏置？→ 既有 model 假設 priority bug 是唯一 amplifier
- Why 2: 為什麼 priority bug 單獨不夠？→ 上游 assignment 已 1.77× 偏 HP1
- Why 3: 為什麼 assignment 偏 HP1？→ longphase getVote() countMap 在 per-read 內 reset 但 cumulative across reads 偏 HP1
- Why 4: 為什麼 cumulative 偏 HP1？→ germline het positions 在 chr19 q-arm 密度遞減 (LOH 區)
- Why 5: 為什麼沒早發現？→ **既有 metric 只看 chr-level ratio 不看 per-region 分佈**

Root cause: metric 設計盲點 — 修正: 加 per-region monitoring

---

⚠ Caveats:
- 不是每個問題都剛好 5 個 why；3-7 個都合理
- 答案要可驗證（不是猜測）
- 避免問人因（「為什麼 X 沒注意」）— 改問系統因（「為什麼 monitoring 沒抓到」）

---

Framework: 5 Whys (大野耐一 / Toyota TPS 1950s-60s)
```
