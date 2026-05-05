# Personal Style Log — 用戶個人 PPT 風格累積

> 每次用戶在 C3/C4/C5 提出**通用必要**修正（依 `prompts/feedback_classification.md` 分類）累積至此。
> AI 在每張 slide build 前**自動讀取所有 `status: active` 規則**，加進 §20.E 6 問 audit / visual_review 10-check / multi_agent_review Wave 1 specific check。
>
> 啟動日：2026-05-05（pptx-build v2.4 機制建立）

---

## 規則紀錄格式

```markdown
### {YYYY-MM-DD} — {規則簡述}（≤ 30 字）

- **觸發來源**：{slide 出處 / master_draft / 場景}
- **規則細節**：{該做 / 不該做}
- **反例**：{違反此規則的具體 case}
- **與既有規則關係**：{補強 §X / 覆蓋 §Y / 新增}
- **檢核方式**：{哪個 Agent / visual_review 加哪一條}
- **狀態**：active / [PROVISIONAL count=N] / archived
```

---

## Active 規則

（尚無累積 — 此檔於 2026-05-05 建立，等待用戶在後續 PPT 製作中累積）

---

## [PROVISIONAL] 觀察中規則

（尚無 — [PROVISIONAL] 累積 ≥3 次同類修正後升級到 active）

---

## Archived 規則（歷史 / 已棄用）

（尚無）

---

## 累積統計

| 指標 | 數值 |
|------|------|
| Active 規則數 | 0 |
| [PROVISIONAL] 規則數 | 0 |
| Archived 規則數 | 0 |
| 最近更新 | 2026-05-05 |

## 與既有規則關係（衝突 / 覆蓋 / 補強）

當用戶通用規則與 playbook 既有規則（§5/§6/§14/§16/§20 等）衝突時：

| 既有規則 | 用戶覆蓋 | 處理 |
|---------|---------|------|
| §6 中文 ≤ 60 字 | 用戶要 ≤ 50 字 | 標 [OVERRIDE]，套用更嚴標準 |
| §14 紅色僅 caveat | 用戶要紅色強調 | 標 [OVERRIDE]，但提示違反 colorblind safe |
| ... | ... | ... |

[OVERRIDE] 標記讓 AI 知道是用戶有意識的覆蓋，不是 bug。

## 後續檢核整合（自動）

每張 slide build 前 AI 執行：

```python
# pseudo-code
def pre_build_check(slide_n):
    personal_rules = load_active_rules('style_library/personal_style_log.md')

    # 加進 §20.E 6 問補充
    extended_audit = base_6_questions + [
        f"是否符合用戶規則 '{rule['name']}'？" for rule in personal_rules
    ]

    # 加進 visual_review 10-check 擴充
    extended_checks = base_10_checks + [
        rule['check_method'] for rule in personal_rules
    ]

    # 加進 multi_agent_review Wave 1
    for agent in ['T', 'C', 'L', 'B']:
        agent_specific = filter_by_agent(personal_rules, agent)
        extend_agent_prompt(agent, agent_specific)

    return run_all_checks(...)
```

## 規則生命週期

```
新規則
   ↓ feedback_classification 分類 = 通用
Active (status: active)
   ↓ 用戶長期未提及 + AI 無觸發
經過 N 次 PPT 都未引用
   ↓ memory-consolidation 建議 archive
Archived (歷史保留，不主動套用)
   ↓ 用戶手動恢復
Active 再次套用
```

## 範例：未來規則寫法（示意，2026-05-05 尚無真實案例）

```markdown
### 2026-XX-XX — 中文標題禁用問號

- **觸發來源**：v3 24-slide audit 中 S11 標題「為何 self-phasing 會卡 ISM？」
- **規則細節**：assertion-evidence 標題用陳述句，不用疑問句（「Self-phasing 卡住 ISM 5 目標」優於問號）
- **反例**：「為何 X 會 Y？」「How does X work?」
- **與既有規則關係**：補強 §5 第 1 條 heading=thesis sentence
- **檢核方式**：Agent-T 字體 audit 加「標題無 ?」項
- **狀態**：active
```

## 與 memory 系統互動

InterSubMod 已有 user-level Memory（`/bip7_disk/.../memory/`），用 `feedback_*` 類型已紀錄 PPT 偏好（如 `feedback_pptx_director_storyboard`、`feedback_pptx_visual_first_philosophy`）。

**分工**：
- **Memory `/feedback_*.md`** — 高層哲學偏好（如「視覺優先」「導演審查」）
- **本檔 personal_style_log.md** — 細粒度可檢查規則（字數 / 字級 / 顏色 / 對齊）

兩者互補：Memory 是「為什麼」，本檔是「具體怎麼做」。
