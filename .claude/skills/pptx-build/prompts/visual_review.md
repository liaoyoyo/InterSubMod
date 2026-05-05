# visual_review.md — Claude Vision 10-Checkpoint

## 使用時機
每張 slide build 後立即跑（C4 前）。

## 10-Check 表格

| # | 項目 | PASS / FAIL / PARTIAL | 理由 |
|:-:|------|---------------------|------|
| 1 | Heading = thesis sentence | | |
| 2 | 1 idea/slide | | |
| 3 | ≤ 6 elements | | |
| 4 | 文字密度（中 ≤ 60 / 英 ≤ 30） | | |
| 5 | Distracted takeaway | | |
| 6 | 視覺對比（顏色/字型/字級） | | |
| 7 | 圖表 vs 文字 ≥ 60% 視覺 | | |
| 8 | 引文 / 數據來源 | | |
| 9 | 1 minute timing check | | |
| 10 | 技術災難備援 | | |

## Pass 標準
- ≥ 8 PASS = 通過
- 1-2 FAIL = 修正後再 review
- ≥ 3 FAIL = 重新設計 slide

## 從 PNG 讀取
```python
from ppt_toolkit.claude_vision_review import review_slide_png
result = review_slide_png(png_path, focal_point="...", section_thesis="...")
print(result['pass_rate'], result['issues'])
```

## P5 整體 timing alarm
P5 結束時加總所有 slide 的 1 minute check：
- 24 張 × 1 min = 24 min vs 預設 30 min → OK
- 字數 24,412 → 60 min vs 30 min → 警告，建議拆 Tier 3

---

## v2.3 補強：Vision 圖像檢視具體流程（單 Agent 10-check）

### 6 步驟標準流程

1. **Read PNG**：用 Read tool 讀 wireframe PNG（Claude 多模態自動視覺化）
2. **逐項對照 10 條**：每條問「我從圖片看到什麼？符合標準嗎？」
3. **PASS/FAIL/PARTIAL**：每條給判定 + 一句具體理由（不能空泛「不錯」）
4. **issue 清單**：FAIL 列具體位置（如「右下角 caveat 字體 < 9pt」）
5. **修正建議**：給具體 ppt_toolkit API call（如「改用 caveat_red_strip.render(font_pt=10)」）
6. **整合到 slide_confirm.md C4**：Vision 結果作為 C4 AskUserQuestion 顯示內容

### 整份 PPTX 後升級到 multi_agent_review.md（並行 6 Agent 深度檢核）

當每張 slide 通過 visual_review 10-check 後，整份 PPTX 完成時：
- **Wave 1 並行 4 Agent (T/C/L/B)** — 結構+視覺（必用 Vision 看 PNG）
- **Wave 2 並行 2 Agent (S/D)** — 內容+整合（看文字檔）

詳見 `prompts/multi_agent_review.md`。

### 三層保險機制

| 層 | 階段 | 工具 | 範圍 |
|:-:|------|------|------|
| 1 | build 前 | §20.E 6 問 self-audit | 每張 slide |
| 2 | build 後 | visual_review 10-check（單 Agent + Vision）| 每張 slide |
| 3 | 整份完成 | multi_agent_review Wave 1+2（6 Agent 並行）| 每張 + 整份 |

任一層 FAIL → 修正後重跑該層。
