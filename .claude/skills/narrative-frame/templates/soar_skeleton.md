# SOAR Skeleton

> Situation / Obstacle / Action / Result — 含「障礙」敘事

```markdown
---
framework: SOAR
when: <面試含障礙敘事 / 創業故事 / 突破型 case>
---

【S】Situation: <情境背景>

【O】Obstacle: <**障礙** — 為什麼這事不簡單；意外 / 阻力 / constraint>

【A】Action: <如何克服 — 強調 problem-solving>

【R】Result: <量化結果 + 學到的教訓>

---

範例:

【S】2026-Q1 設計 V3F 修正 chr19 priority bug
【O】V3F 修 chr19 但 V5 在 germline-absent 區仍 4.19:1 偏 HP1 — 發現 Layer 1.5 設計缺陷繼承 priority bug
【A】2026-05 設計 V6 同時處理 assignment + priority 兩層；跑 paired pileup + LOSO 驗證
【R】100% chr19 修復 + HCC1395 F1 ↑0.022（vs V3F）；學到 「metric 設計盲點」是 12 月 audit 失效根因

---

Framework: SOAR (Behavioral Interviewing 變體)
```
