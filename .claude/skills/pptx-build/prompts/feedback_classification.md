# feedback_classification.md — 用戶修正分類 + 個人風格累積（v2.4 新增）

> 用戶提出修正/改善建議時，AI 必須**強制分類**為「通用」vs「本次特定」，避免 over-fit 單張 slide 污染通用規則。通用修正累積到 personal_style_log，作為後續檢核標準。

## 觸發時機

- C3 section confirm 後用戶說「這張不對」
- C4 slide confirm 後用戶提修正
- C5 speaker review 後用戶調整
- Wave 1+2 多 Agent 結果用戶不認同某條
- 任何用戶說「這個排版要改」「文字太密」「顏色不對」「我不喜歡 X」時

## 強制分類流程

### Step 1：AI 收到用戶修正要求 → 立即 AskUserQuestion 分類

```yaml
question: "您剛提的修正是哪一類？"
header: "修正分類"
multiSelect: false
options:
  - label: "通用必要（所有 PPT 都該如此）"
    description: "→ 累積到 InterSubMod/.claude/skills/pptx-build/style_library/personal_style_log.md，下次自動套用"
  - label: "本次特定（只是這張 slide / 這份 PPT 的特殊情境）"
    description: "→ 只本次修正，不污染通用規則"
  - label: "不確定，請 AI 協助判斷"
    description: "→ AI 用 5 維度評估後給建議分類"
  - label: "兩者皆是（先本次修正，但若 N 次再出現則升級通用）"
    description: "→ 標 [PROVISIONAL]，等 ≥3 次出現自動升級到 personal_style_log"
```

### Step 2：依分類處置

#### 通用修正
- 寫入 `InterSubMod/.claude/skills/pptx-build/style_library/personal_style_log.md`
- 加日期 + slide 出處（哪張 slide 觸發）+ 規則描述 + 反例
- 後續所有 PPT 自動套用（在 §20 紅旗清單中加新項）
- 若衝突既有規則，標 [OVERRIDE] 並說明

#### 本次特定
- 只在本次 PPT 修正
- 不寫入 personal_style_log
- 留 inline 註解（如：「此 slide caveat 字體加大為 14pt — 因為 audit 主題特殊」）

#### AI 協助判斷（5 維度）

| 維度 | 通用訊號 | 本次特定訊號 |
|------|---------|------------|
| 觸發頻率 | 過去多份 PPT 都有此問題 | 只此張獨特 |
| 普適性 | 適用所有 audience tier | 限特定教授/同儕 |
| 與既有規則關係 | 補強既有規則 | 違反既有規則但合理 |
| 修正成本 | 寫 1 條規則永久解決 | 每次手動處理 OK |
| 用戶語氣 | 「以後都這樣」「我喜歡」 | 「這次特殊」「這張比較複雜」 |

≥ 3 通用訊號 → 建議「通用」；否則「本次特定」。

#### [PROVISIONAL] 升級規則

每次標 [PROVISIONAL] 累積計數，當：
- 同一規則在 ≥ 3 次不同 PPT 出現
- → AI 主動建議用戶升級為通用，寫入 personal_style_log

## Step 3：紀錄格式

寫入 personal_style_log 時用標準格式：

```markdown
### {YYYY-MM-DD} — {規則簡述}（≤ 30 字）

- **觸發來源**：{slide 出處 / 哪份 master_draft / 何種場景}
- **規則細節**：{完整描述，含「該做」「不該做」}
- **反例**：{違反此規則的具體 case}
- **與既有規則關係**：{補強 §X / 覆蓋 §Y / 新增}
- **檢核方式**：{哪個 Agent 該檢查 / visual_review 加哪一條}
- **狀態**：active / [PROVISIONAL count=N] / archived
```

## Step 4：後續檢核

每張 slide build 前自動讀取 `personal_style_log.md` 中所有 `status: active` 規則：
- 加進 §20.E 6 問 audit 的「補充問」
- 加進 visual_review 10-check 的擴充項
- 加進 multi_agent_review Wave 1 各 Agent 的 specific check

逐步累積 → 個人化 checklist 自動套用。

## 失敗模式

| 情境 | 處置 |
|------|------|
| 用戶不分類，直接「都改」 | AI 強制 AskUserQuestion 分類後才修；不允許跳過 |
| 用戶選通用但規則太特殊（如「這份用紅色強調」）| AI 反問：「下份 PPT 也用紅色嗎？」確認普適性 |
| personal_style_log 規則衝突（如「強調用紅」vs 既有 §14 紅僅 caveat）| 標 [OVERRIDE]，註明用戶覆蓋既有規則 |
| 規則累積過多（>30 條）| memory-consolidation 協助合併 / 過時 archive |

## 範例對話

```
用戶：「這張 slide 文字太密，請減少」
AI: [強制分類]「請問這是：(A) 通用要求 (B) 本次特定 (C) 不確定」

情境 A — 通用：
  AI: 已寫入 personal_style_log：「中文 ≤ 50 字 per slide（嚴於既有 60 字上限）」
       下次所有 PPT 自動套用此緊縮上限。

情境 B — 本次特定：
  AI: 本次修正此張至 50 字以內，不寫 log。

情境 C — AI 判斷：
  AI: 5 維度評估：觸發頻率 2/3 過去 PPT 也有 / 普適性 高 / ... → 建議「通用」
       是否確認？(yes/no)
```

## 與其他 prompts/ 關聯

- 觸發於：每 C3/C4/C5 用戶 free-text 修正後
- 上游：`slide_confirm.md` / `visual_review.md` / `multi_agent_review.md` 結果
- 下游：寫入 `style_library/personal_style_log.md`
- 後續檢核：每張 slide build 前讀取 personal_style_log
