<!--
建立時間: 2026-05-18
parent: P3 audit M10 fix (Eisenhower 四象限)
purpose: 用戶端 task triage template — 緊急 × 重要 4 象限
status: optional template, copy when needed
-->

# Task Priority Matrix Template — Eisenhower 四象限

> **使用方式**: 對於 backlog / pending items / 多任務優先級決策時，複製此 4 象限後填空。
> **不適用**: 單一任務（直接做即可）/ 線性 plan workflow（用 plan mode）/ Hard Gate 任務（依 CLAUDE.md §1）。
> **預期效果**: 5 分鐘內把模糊 backlog 分類完成，明確下個動作。

---

## Matrix Template

```
                  緊急 (Urgent)              不緊急 (Not Urgent)
              ┌────────────────────────┬────────────────────────┐
              │ Q1: DO NOW (危機)       │ Q2: SCHEDULE (戰略)    │
   重要        │  - <task>              │  - <task>              │
 (Important)  │  - <task>              │  - <task>              │
              │  ⇒ 立刻處理            │  ⇒ 排期 / 主動投入     │
              ├────────────────────────┼────────────────────────┤
              │ Q3: DELEGATE (打擾)     │ Q4: ELIMINATE (浪費)   │
  不重要       │  - <task>              │  - <task>              │
(Not Important│  - <task>              │  - <task>              │
              │  ⇒ 委派 / 簡化         │  ⇒ 刪除 / 延後不做     │
              └────────────────────────┴────────────────────────┘
```

---

## 4 象限定義（業界共識）

| 象限 | 名稱 | 判別條件 | 處置 |
|------|------|---------|------|
| **Q1** | DO NOW / 危機處理 | 緊急 + 重要 — deadline 24-48 hr / 有外部依賴 / 阻擋他人 | 立刻處理，最高 priority |
| **Q2** | SCHEDULE / 戰略投入 | 不緊急 + 重要 — 長期影響 / 預防性工作 / 能力建構 | **真正生產力來源**，排期主動投入（每週 fixed time slot） |
| **Q3** | DELEGATE / 緊急打擾 | 緊急 + 不重要 — 他人請求 / 短期回應壓力 | 委派 / 自動化 / 設標準回應 |
| **Q4** | ELIMINATE / 浪費 | 不緊急 + 不重要 — 干擾 / 低 ROI 活動 | 刪除 / 延後 / 拒絕 |

---

## InterSubMod 場景對應

| 任務類型 | 典型象限 | 範例 |
|---------|---------|------|
| Hard Gate (C++ commit / NO-GO / 刪檔) | **Q1** | V6 production tag W3 deadline |
| Fix backlog (P3/P4 from audit) | **Q2** | router_pre_bash 整合 / memory archive |
| Plugin 通知 / 用戶反問 | **Q3** | 簡單 confirm dialog / 路徑詢問 |
| 過時 stale conclusion 重訪 | **Q4** | 已 NEGATIVE 30d+ 想再試（需 §8.3.1 reopen 三條件才升 Q2） |
| 研究方向 pivot | Q1 或 Q2 | 依 deadline 緊急度 |

---

## 與 InterSubMod skills 的關係

- **`/scientific-rigor §7.1 Pre-registration`**: confirmatory cycle 是 Q2 戰略；exploratory pilot 可能是 Q3（不該佔太多時間）
- **`/scientific-rigor §8.3.1 Reopen Threshold`**: 想 reopen NEGATIVE = 預設 Q4，需 C1/C2/C3 一條才能升 Q2
- **`feedback_small_scale_validation_first`**: pilot < 2hr 對應 Q3 快速委派風格（不投入 Q2 戰略級資源）
- **`feedback_execution_mode_hierarchy`**: 全自動 fan-out 適合 Q2 已驗證階段；Q1 新概念需確認

---

## 判別流程（給 AI 用）

當用戶丟出 backlog 或多任務時：

1. **詢問 deadline**（如無 → 默認 Q2 或 Q4）
2. **評估影響範圍**（會否阻擋他人 / 影響結論 → Q1 / Q2）
3. **填入 4 象限表**
4. **Q1 優先；Q2 排期；Q3 委派；Q4 拒絕**
5. **每週審視 Q4 → 若仍存在 30d+ 強制刪除或重評估**

---

## 範例（2026-05-18 InterSubMod backlog 應用）

```
                  緊急 (W3 5/22 內)          不緊急 (W4+)
              ┌────────────────────────┬────────────────────────┐
              │ Q1:                     │ Q2:                    │
   重要        │ - V6 production tag    │ - P3 fix items (router │
              │ - PI errata email       │   / worktree / paths)  │
              │ - COLO829 truth set     │ - P4 fix (injection   │
              │   權限決策              │   guard / memory      │
              │                         │   archive)            │
              ├────────────────────────┼────────────────────────┤
              │ Q3:                     │ Q4:                    │
  不重要       │ - skill_change_audit   │ - 6S 現場管理整合      │
              │   log monthly review    │   (重複既有)           │
              │                         │ - 7 habits framework   │
              │                         │   對齊 (scope 不符)   │
              └────────────────────────┴────────────────────────┘
```

---

## 元層說明

本模板源自 2026-05-18 P3 Methodology Audit M10 fix（priority P3）。原始分析見 `InterSubMod/docs/reports/validated/2026/05/20260518_agent_harness_audit_p1_skills_01/p3_methodology_alignment.standalone.html`。
業界參考：Stephen Covey *7 Habits* Quadrant II + David Allen *Getting Things Done* (GTD) capture-clarify-organize-reflect-engage。
