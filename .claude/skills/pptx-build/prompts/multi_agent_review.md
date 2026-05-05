# multi_agent_review.md — Phase 5 6 Agent 多代理驗證（並行）

> 補強 playbook §24「Phase 5 多 Agent 驗證」具體 context + 步驟 + Vision 整合。
> 用 `Agent` tool spawn `feature-dev:code-reviewer` subagent_type，每個 Agent 獨立 context 並行執行。

## 觸發時機

- **每張 slide build 完成 + wireframe 截圖後** → Wave 1 並行 4 Agent
- **整份 PPTX 完成 + 所有 speaker note 寫好後** → Wave 2 並行 2 Agent
- **C5 speaker script timing 後** → 必跑兩 Wave 才能交付

## Wave 1：結構 + 視覺（並行 4 Agent）

每張 slide 同時 spawn 4 個 subagent，各自獨立讀 PNG + style_library YAML，產出 PASS/FAIL/PARTIAL 表。

### Agent-T (Typography 字體)

```yaml
subagent_type: feature-dev:code-reviewer
context_inputs:
  - {slide_path}.png  (Vision: 必看)
  - InterSubMod/.claude/skills/pptx-build/style_library/objects/thesis_title_bar.yaml
  - InterSubMod/.claude/skills/pptx-build/playbook.md §6 字體規範
prompt: |
  你是 Typography Agent。檢查 slide PNG 的字體：
  1. 標題 ≥ 32pt Arial Bold + CJK fallback Droid Sans Fallback？
  2. 內文 ≥ 14pt Arial + CJK fallback？
  3. 字級下限 9pt（footnote/citation 例外）？
  4. Latin 字符在 CJK 字型下無方塊（每個字必檢查）？
  5. 標題粗體、引文 italic 規則一致？
output_format: PASS/FAIL/PARTIAL × 5 項 + 具體位置 + 修正建議
```

### Agent-C (Color 色彩)

```yaml
context_inputs:
  - {slide_path}.png  (Vision: 必看)
  - InterSubMod/.claude/skills/pptx-build/style_library/colors/palette.yaml
prompt: |
  你是 Color Agent。檢查 slide PNG 的色彩：
  1. 用色是否在 6 色 token（red/green/blue/amber/grey/black）內？
  2. 紅色用於 caveat / 綠色用於 insight / 藍色用於 method 是否正確？
  3. WCAG AA 對比 ≥ 4.5:1（body text on white）？
  4. 紅綠相鄰嗎（colorblind 風險）？
  5. ≤ 3 colors per slide（不含黑白）？
output_format: PASS/FAIL/PARTIAL × 5 項
```

### Agent-L (Layout 佈局)

```yaml
context_inputs:
  - {slide_path}.png  (Vision: 必看)
  - InterSubMod/.claude/skills/pptx-build/style_library/layouts/{used_layout}.yaml
  - InterSubMod/.claude/skills/pptx-build/playbook.md §16 對齊規則 + §20.B focal_point_zone
prompt: |
  你是 Layout Agent。檢查 slide PNG 的版面：
  1. 元素數量 ≤ 6（Tier 1 element count）？
  2. 對齊規則符合 §16（標題距上 0.4in / 視覺中心 / 多欄寬度）？
  3. focal_point_zone 與 layout YAML 標註一致？
  4. 留白 20-30% / 視覺 55-65% / 文字 ≤ 15%？
  5. caveat strip 對齊被警示元素左邊緣？
output_format: PASS/FAIL/PARTIAL × 5 項 + 對齊偏移量化（in 為單位）
```

### Agent-B (Bilingual 雙語)

```yaml
context_inputs:
  - {slide_path}.png  (Vision: 必看)
  - InterSubMod/.claude/skills/pptx-build/playbook.md §6 雙語規範
prompt: |
  你是 Bilingual Agent。檢查 slide PNG 的中英文格式：
  1. 中文主、英文副（英文 60% 字級 + 縮排 0.25"）？
  2. 中文 ≤ 60 字 / 英文 ≤ 30 word per slide（不含標題）？
  3. 標題雙語格式一致？
  4. emoji（⭐🏅）有無（CJK 字型可能渲染失敗）？
  5. 引文格式中英一致？
output_format: PASS/FAIL/PARTIAL × 5 項
```

## Wave 2：內容 + 整合（並行 2 Agent）

整份 PPTX 完成後同時 spawn 2 subagent，看全份檔案。

### Agent-S (Speaker Notes)

```yaml
context_inputs:
  - {pptx_dir}/03_slide_layout_and_script.md
  - {master_draft_dir}/master_draft.md (frontmatter 取 target_duration_min)
  - InterSubMod/.claude/skills/pptx-build/playbook.md §6 timing
prompt: |
  你是 Speaker Notes Agent。檢查整份 speaker note：
  1. 每張 75-90 sec ≈ 中 300-360 字（量化：每張字數）？
  2. 加總 vs target_duration_min（25 min × 400 字/min = 10,000 字）？
  3. [ORAL-OPTIONAL] 標記在 Tier 3 段落出現？
  4. 雜訊用語（「順便提一下」「附帶說明」「另外」）已削減？
  5. 每張 speaker note 有 evidence path 引用？
output_format: 每張 sec 估計表 + 總時長 vs target + 雜訊紅旗清單
```

### Agent-D (Data 數據準確性)

```yaml
context_inputs:
  - {pptx_dir}/01_full_narrative.md + 02_slide_outline.md + 03_slide_layout_script.md
  - master_draft.md（取所有 source_artifacts）
  - 所有母稿引用的 evidence path（必讀至少 3 個）
prompt: |
  你是 Data Accuracy Agent。檢查 PPTX 中所有數據：
  1. 每個數字（如 17.3:1 / 0.7166 / 938f0df）是否與 source 一致？
  2. caveat 完整（如「PI 報告 Pass 1 only」是否每處引用都有）？
  3. evidence path 真實存在（grep 確認）？
  4. 無過度宣稱（[F] 的數據不超出 source 範圍）？
  5. 無流水帳（每張 slide 有具體因果連接）？
output_format: 數據錯誤清單（slide # / 錯誤數字 / 正確值 / source）+ caveat 缺失清單
```

## Pass 標準

每張 slide：
- Wave 1 ≥ 16/20 PASS（4 Agent × 5 項）= 通過
- Wave 1 ≥ 12 PASS = 修正後 re-review
- Wave 1 < 12 PASS = 重新設計 slide

整份 PPTX：
- Wave 2 數據錯誤 = 0 + caveat 缺失 = 0 + 雜訊紅旗 = 0 → 交付
- 任一 > 0 → 修正後再 review

## Vision 圖像檢視整合（Agent-T/C/L/B 必跑）

每個 Wave 1 Agent 都用 Claude Vision 讀 PNG：

```python
# 流程（playbook §7 已定義 ppt_toolkit.claude_vision_review）
from ppt_toolkit.claude_vision_review import review_slide_png

# 並行 4 Agent，各自獨立 Vision 看 PNG
def wave1_review(png_path, layout_yaml, palette_yaml):
    results = {}
    # spawn 4 subagent in parallel via Agent tool
    for agent in ['T', 'C', 'L', 'B']:
        # 每個 Agent 獨立 Vision call + 各自 PASS/FAIL 表
        results[agent] = spawn_agent(
            type='feature-dev:code-reviewer',
            png=png_path,
            specific_context=f"prompts/multi_agent_review.md#agent-{agent}"
        )
    return integrate(results)
```

關鍵：**每個 Agent 獨立 Vision call**，避免 context 污染。

## 並行調度（用 Agent tool 同時 spawn）

```
Wave 1 (per slide):
    parallel [Agent-T, Agent-C, Agent-L, Agent-B]
    → 4 subagent 同時跑 Vision + 結構檢查
    → 各自獨立 PASS/FAIL 表
    → 主 conversation 整合（取交集 = critical issues）

Wave 2 (整份 PPTX):
    parallel [Agent-S, Agent-D]
    → 2 subagent 同時跑全文掃描
    → 各自獨立輸出
```

## 修正循環

| 階段 | 觸發 | 動作 |
|------|------|------|
| 第 1 round Wave 1 | 每張 slide build 後 | 並行 4 Agent → 整合 issue 清單 |
| 修正 | issue > 0 | 主 conversation 修正 PPTX 程式碼或內容 |
| 第 2 round Wave 1 | 修正後 | 重跑 4 Agent → 確認 PASS |
| 第 1 round Wave 2 | 整份完成 + Wave 1 全 PASS | 並行 2 Agent → 整合 issue |
| 修正 | issue > 0 | 修正後重跑 |
| 第 2 round Wave 2 | 修正後 | 確認 PASS → 交付 |

## 與既有 prompts/ 的關聯

- `visual_review.md` 是 **單 Agent 10-checkpoint**（每張快速）
- `multi_agent_review.md` 是 **多 Agent 並行深度檢核**（每張 Wave 1 + 整份 Wave 2）
- 流程：visual_review (10-check) → 通過後 multi_agent_review (Wave 1) → 整份 multi_agent_review (Wave 2)

每張 slide 三層保險：
1. §20.E 6 問 audit（self-audit，build 前）
2. visual_review 10-check（單 Agent，build 後）
3. multi_agent_review Wave 1（4 Agent 並行，build 後）

## 失敗模式處理

| 情境 | 處置 |
|------|------|
| Agent 之間結論衝突（如 T 說 PASS、L 說 FAIL）| 以更嚴格者為準（FAIL 優先），主 conversation 確認 |
| Vision 看不清 PNG（解析度不足）| 重產 PNG with `--dpi 300` + 重跑 Agent |
| Agent 跑超時（>2 min/Agent）| 拆 prompt 為更小子任務 |
| 6 Agent 全 PASS 但用戶肉眼看仍有問題 | 補新 Agent rule（迭代改進）|
