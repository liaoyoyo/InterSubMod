# PDCA Skeleton

> W. Edwards Deming (基於 Walter Shewhart) — 持續改善 cycle

```markdown
---
framework: PDCA
when: <持續改善 / 實驗設計 / iterative development>
---

# Cycle <N>: <Title>

## 1. Plan（計畫）

- **Hypothesis**: <假設>
- **Expected outcome**: <預期 metric>
- **Procedure**: <步驟>
- **Success criteria**: <pass / fail threshold>

## 2. Do（執行）

- <實際執行紀錄 — 含 deviation if any>
- <data collected>

## 3. Check（檢核）

- **Actual outcome**: <實際 metric>
- **vs Expected**: PASS / FAIL / PARTIAL
- **Surprises**: <意外發現>

## 4. Act（行動 — 下一輪輸入）

- **Standardize** (如 PASS): <納入 SOP>
- **Adjust** (如 PARTIAL): <修正後重跑>
- **Abandon** (如 FAIL): <NO-GO + 記教訓>

→ Cycle <N+1> 從哪裡接續

---

範例（Cycle 4 H_NEW_4 LOSO 驗證）:

1. **Plan**: drop caller_af 後 cross-sample LR 應 ≥2/5 ΔF1 > +0.002 (prior 25%)
2. **Do**: 跑 5-sample LOSO with 9 features
3. **Check**: 0/5 above threshold；但 HCC1395 unexpected +0.00699（VIOLATED prior）
4. **Act**: PARTIAL — Cycle 5 pivot 至 chr8 zone gate (Path B) 驗證 +0.00699 generalize

---

Framework: PDCA (Deming《Out of the Crisis》1982; Shewhart 1939)
```
