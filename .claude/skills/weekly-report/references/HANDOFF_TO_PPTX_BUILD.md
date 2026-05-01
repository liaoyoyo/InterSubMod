# 母稿 → pptx-build skill handoff 規範

本檔定義 weekly-report W7 → C4 之後，將 master_draft.md 銜接給 pptx-build sub-skill 的協議。

---

## 1. 母稿輸出位置

```
InterSubMod/docs/reports/validated/YYYY/MM/YYYYMMDD_週報_主題/
├── master_draft.md                     # 主檔（Layer 0-4 + 17 段）
├── next_week_plan.md                   # 選 D 時加產
└── (later: pptx-build artifacts)       # 選 A 時 pptx-build 產出
```

**版本控制**：
- 第 1 版：`master_draft.md`
- 用戶要求修改時：`master_draft_v2.md`、`master_draft_v3.md`...
- 最終版（C4 確認）：移除版本號或保留最新；舊版移到 `archive/` 子資料夾

---

## 2. 母稿 frontmatter schema（必填）

```yaml
---
title: 週報 YYYYMMDD - <主題>
type: weekly_master_draft
date: YYYY-MM-DD
status: ready_for_handoff   # ready_for_handoff | draft | finalized
report_type: progress | problem | consult | exploration   # W2 → C1 主線類型
main_statement: "≤ 30 字"   # W2 → C1 主線一句話
audience: advisor             # advisor | PI | committee | peer
target_duration_min: 20       # 預期報告時長（教授週報通常 20-30 min）
source_artifacts:
  - InterSubMod/docs/experiments/in_progress/...
  - InterSubMod/research/autoresearch/evidence_ledger.jsonl#L1234
material_classification:
  facts: 5         # [F] count
  observations: 3  # [O] count
  inferences: 2    # [I] count
  unconfirmed: 1   # [U] count
priority_buckets:
  ppt: 6           # PPT 桶數量
  speaker_note: 8  # 講稿桶
  appendix: 3      # 備註桶
  shelf: 2         # 暫存桶
suggested_pptx_template: improvement_report   # W7 §9 建議 pptx-build 模板
estimated_pptx_slides: 18    # W7 §10 建議投影片張數
professor_qa_count: 6        # W6 → C3 教授追問數
handoff_choice: pending      # pending | A | B | C | D（C4 後填）
---
```

---

## 3. C4 後的 4 選 handoff prompt 模板

```
AskUserQuestion:
question: "母稿已完成（17 段，N 字，Layer 0-4 結構齊全）。下一步？"
header: "Handoff"
multiSelect: false
options:
  - label: "A. 立即觸發 pptx-build"
    description: "自動進入 pptx-build skill，讀取本母稿 frontmatter，跳過 main thesis 鎖定。預估 1-2 hr 完成 PPTX。"
  - label: "B. 母稿留檔，下次手動觸發"
    description: "母稿已存於 InterSubMod/docs/reports/validated/.../master_draft.md。下次用 /pptx-build --from-draft <path> 啟動。"
  - label: "C. 母稿即終點"
    description: "本次任務不產 PPT（如 internal status check）。週報任務結束。"
  - label: "D. 母稿留檔 + 加寫下週計畫"
    description: "母稿保留，並加產 next_week_plan.md（從母稿 §16 下一步行動 + §17 教授追問 中萃取）。"
```

**互動模式**：用戶必選 1。
**全自動模式**：預設選 B（留檔），不自動觸發 pptx-build（避免連續長時間工作）。

---

## 4. pptx-build 從母稿讀取規格

當用戶選 A，或手動執行 `/pptx-build --from-draft <path>` 時：

### 4.1 自動讀取的欄位

| 母稿欄位 | pptx-build 用途 | 對應 P 階段跳過項 |
|---------|---------------|----------------|
| `report_type` | 觸發對應 templates/{type}.md | P1 報告類型識別 |
| `main_statement` | 鎖定 main thesis | P1 main thesis 鎖定 |
| `audience` | 設定 audience tier | P1 Audience & Goal 部分 |
| `target_duration_min` | 估算 slide 數 + speaker note 字數 | P1 / P5 計算 |
| `suggested_pptx_template` | 6 模板（improvement/comparison/...）| P1 模板選擇 |
| Layer 0.1 + §1 + §2 | thesis_title_bar + cover slide | P2 outline 起點 |
| Layer 2 Thread × N | 1-2 張 slide per Thread | P3 section 結構 |
| Layer 2 §3 [F] / §4 [O][I] / §5 [U] | slide 內容 + 4 層分類標籤帶過 | P4 slide build |
| Layer 3 整合 | 1 張 integration slide | P3 section |
| Layer 4 §16 + §17 | future tree + Q&A backup slides | P3 section + appendix |

### 4.2 pptx-build 跳過 P 階段項目

從母稿讀取時：
- ❌ 跳過 P1 Audience & Goal 中的 main thesis 鎖定（已由 W2 完成）
- ❌ 跳過 P1 報告類型識別（已由 W2 完成）
- ❌ 跳過 §20 Tier 評分（4 層分類已由 W3 完成，[F]→S/[O]→A/[I]→B/[U]→C 自動 mapping）

### 4.3 pptx-build 仍保留的 P 階段項目

從母稿讀取時，pptx-build 仍須走：
- ✅ P2 Outline checkpoint（拆 5-7 段，每段 thesis）
- ✅ P3 Section checkpoint（每 section 5-7 slide 標題 + focal point）
- ✅ P4 Slide checkpoint（逐張 build → wireframe → Vision 10-check）
- ✅ P5 Speaker script（從母稿 §14 講稿例子延伸）
- ✅ §22 v3 audit + 雜訊紅旗檢查

---

## 5. D 選項：next_week_plan.md 自動產生規格

選 D 時，weekly-report 自動萃取以下內容寫入 `next_week_plan.md`：

```markdown
---
title: 下週計畫 YYYYMMDD（基於 master_draft.md）
type: next_week_plan
parent_master_draft: <path>
date: YYYY-MM-DD
---

# 下週優先行動

## 1. 本週母稿 §16 萃取

[從 master_draft.md §16 下一步行動清單 全部複製]

## 2. 教授追問預備（§17 萃取）

[從 master_draft.md §17 教授可能提問 + 回答準備]

## 3. 待補資料追蹤（§11 + §12 + §13）

[從 master_draft.md §11/§12/§13 萃取]

## 4. 暫存項目延後評估（§15）

[從 §15 列項，標 [SHELVED] 等下週判斷是否升級]

## 5. blockers（從 §5 待確認項目）

[從 §5 [U] 標籤項目轉成 blocker 追蹤]

# 下週時間預估

| 項目 | 預估 hr |
|------|---------|
| ... | ... |
```

---

## 6. 銜接後的 pptx-build outputs

最終資料夾結構（用戶選 A 後）：

```
InterSubMod/docs/reports/validated/YYYY/MM/YYYYMMDD_週報_主題/
├── master_draft.md
├── pptx_build/
│   ├── 00_storyboard.md            # pptx-build P3.5 Storyboard
│   ├── 01_full_narrative.md        # pptx-build P2 outline
│   ├── 02_slide_outline.md         # pptx-build P3 section list
│   ├── 03_slide_layout_script.md   # pptx-build P4 + P5
│   ├── build_pptx.py               # 生成腳本
│   ├── output.pptx                 # 最終 PPTX
│   └── wireframes/                 # 截圖驗證
└── audit/
    ├── visual_review_results.csv   # 10-check 結果
    └── focal_point_audit.csv       # §20 6 問 audit
```

---

## 7. Failure modes

| 情境 | weekly-report 處置 |
|------|-------------------|
| 用戶選 A 但 pptx-build skill 未載入 | 自動觸發失敗 → 退回 B（留檔）+ 提示「請手動 /pptx-build --from-draft」|
| 母稿 frontmatter 不完整 | C4 確認時擋下，要求補欄 |
| `report_type` 不在 6 模板列表 | 預設用 improvement_report + 提示用戶調整 |
| 用戶選 A 但 frontmatter 中 `status != ready_for_handoff` | 警告 + 要求先完成 C4 確認 |

---

## 8. 與 pptx-build skill 的契約

pptx-build skill 必須支援：

1. CLI 觸發：`/pptx-build --from-draft <path>`
2. 自動偵測：用戶說「我要做簡報」+ 提供母稿 path 時自動 from-draft
3. frontmatter 解析：讀取本檔 §2 schema 中所有欄位
4. 跳過項目處理：P1 main thesis 鎖定 + P1 報告類型識別 + §20 Tier 評分
5. P2-P5 仍走完整流程
6. 輸出位置：`{master_draft_dir}/pptx_build/`（與母稿同目錄下的子資料夾）

詳細 pptx-build P 階段流程 → `InterSubMod/.claude/skills/pptx-build/playbook.md`（Phase 2 落地後可用）。
