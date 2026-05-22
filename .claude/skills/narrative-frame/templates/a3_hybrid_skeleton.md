# A3+ADR+Postmortem-hybrid Skeleton (= 既有 structured-tech-report 13 段)

> Toyota A3 + ADR + Google SRE Postmortem + Diátaxis 多源融合

```markdown
---
framework: A3+ADR+Postmortem-hybrid (13 段)
when: <單一工程改動 / pipeline 變更 / bug fix>
對應既有 skill: /structured-tech-report (thin wrapper)
---

# <Title — Assertion-Evidence 一句結論>

## 0. TL;DR

<3-5 行 BLUF — 沒空讀全文的人>

## 1. 報告目的

<為什麼寫這份報告？受眾 + 目標>

## 2. 系統背景

<被改動的模組原本負責什麼>

## 3. 原本流程（Before）

<改前的工作方式>

## 4. 問題描述

<什麼情境出什麼錯；量化>

## 5. 根本原因

5 Whys / Fishbone（見 5_whys_skeleton.md / fishbone_skeleton.md）

## 6. 修改方向（ADR Decision）

候選方案 + Decision Drivers + Chosen + Rationale（見 adr_skeleton.md）

## 7. 修改內容

### 7.1 非工程版（白話）
<沒提 file path / function name / commit hash>

### 7.2 工程版（精確）
<file:line + commit hash + API change>

## 8. 新舊比較

| Aspect | Before | After |
|--------|--------|-------|
| ... | ... | ... |

## 9. 驗證方式（Step→Verify）

| Step | Action | Verify |
|------|--------|--------|
| 1 | ... | <external observable> |

## 10. 影響範圍

<受影響使用者 / 資料 / 下游>

## 11. 風險與限制

<未解 confound / 樣本限制 / 假設未驗>

## 12. 後續工作

<Action items + owner + deadline>

## 13. 結論

<3-5 句總結 — 對照 §1 報告目的>

---

Framework hybrid: Toyota A3 + Nygard ADR + Google SRE Postmortem + Procida Diátaxis
對應既有 skill: /structured-tech-report (thin wrapper 預設套此 framework)
```
