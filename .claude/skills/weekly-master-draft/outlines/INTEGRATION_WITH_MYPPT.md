---
title: weekly-master-draft 與 myPPT / weekly-report 三 skill 整合設計
date: 2026-05-01
status: outline-for-review
---

# 三 skill 整合設計

## 1. 既有 skill 盤點

```
[既有 1] weekly-report skill
        位置：（既有，已啟用）
        職責：InterSubMod 週進度自動收集 → markdown 週報 → PPTX → 索引更新
        觸發：週報、weekly report、整理本週、報告
        缺口：不做 4 層分類、不做教授問答預測、不做主線類型識別、PPT 直接生成（無母稿確認步驟）

[既有 2] structured-tech-report skill
        職責：單一工程改動 deep dive 13 段報告
        觸發：寫技術報告、整理修改、self-phasing 報告
        cadence：單次事件，非週期性

[新建 1] myPPT skill（DRAFT，38 檔已落地，16 檔未落地）
        職責：簡報 deck 生成主流程
        觸發：簡報、PPT、教授報告

[新建 2] weekly-master-draft skill（本檔規劃）
        職責：週報母稿（先於 PPT）
        觸發：週報、向教授報告
```

## 2. 三 skill 衝突診斷

**核心衝突**：weekly-report 與 weekly-master-draft 觸發 keyword 高度重疊，且職責有部分重疊（皆涉及週報）。

| 維度 | weekly-report（既有）| weekly-master-draft（新）|
|------|---------------------|------------------------|
| 出發點 | 自動收集 + 草稿 | 互動確認 + 邏輯檢查 |
| 主軸識別 | ❌ 無 | ✅ 4 種主線類型 |
| 內容分類 | ❌ 無 | ✅ [F]/[O]/[I]/[U] 4 層 |
| 教授問答預測 | ❌ 無 | ✅ 5-7 個追問 |
| 紅旗檢查 | ❌ 無 | ✅ 過度宣稱 / 流水帳 |
| 直接產 PPT | ✅ 是 | ❌ 不（僅產母稿） |
| 17 段母稿格式 | ❌ 無 | ✅ 標準格式 |

**結論**：weekly-master-draft 是 weekly-report 的「品質升級層」，不是替代品。需明確分工。

## 3. 三選項：請用戶決策

### 選項 X：**整併方案**（推薦，但需評估 weekly-report 既有用法）
- 把 weekly-master-draft 的 7 階段 + 4 層分類 + 教授問答預測直接升級進 weekly-report
- weekly-report 維持現名（避免外部依賴中斷），內部增加 W1-W7 + C0-C6
- 完成後觸發「要繼續產 PPT 嗎？」 → 銜接 myPPT
- **優點**：不增加 skill 數量、不分裂觸發 keyword
- **缺點**：weekly-report 現有實作可能要大改；既有用戶會感受 prompt 流程變慢

### 選項 Y：**並行方案**（本大綱目前設計）
- weekly-master-draft 與 weekly-report 並存，觸發優先序：
  - 「週報」+「向教授」→ weekly-master-draft（優先）
  - 純「週報」/「整理本週」→ weekly-report（既有快速版）
- weekly-master-draft 內部 W1 raw data 收集可呼叫 weekly-report 既有 collector
- **優點**：用戶可選擇「快速」vs「深度」兩種模式
- **缺點**：兩 skill keyword 衝突，需精準描述觸發優先序

### 選項 Z：**廢棄方案**
- 廢棄 weekly-report，僅留 weekly-master-draft
- 風險：既有用法被打斷
- **不推薦**

## 4. 銜接協議：weekly-master-draft → myPPT

### 4.1 母稿輸出 schema

```yaml
# 母稿 .md frontmatter
---
title: 週報母稿 YYYYMMDD
type: weekly_master_draft
date: YYYY-MM-DD
status: ready_for_handoff
report_type: progress|problem|consult|exploration  # W2 主線類型
main_statement: "≤ 30 字"  # W2 鎖定
audience: advisor  # 預設教授
target_duration_min: 20  # 預設 20 min
source_artifacts:
  - InterSubMod/docs/experiments/in_progress/...
  - research/autoresearch/evidence_ledger.jsonl#L1234
material_classification:
  facts: 5      # [F] count
  observations: 3  # [O] count
  inferences: 2  # [I] count
  unconfirmed: 1  # [U] count
priority_buckets:
  ppt: 6
  speaker_note: 8
  appendix: 3
  shelf: 2
---

## §1 - §17 母稿正文
...
```

### 4.2 myPPT 接母稿時的行為

myPPT 偵測到 `--from-draft <path>` 或母稿 frontmatter 含 `type: weekly_master_draft`：

1. **跳過項目**：
   - §1 Audience & Goal 的 main thesis 鎖定（已由 W2 完成）
   - §1 主線類型識別（已由 W2 完成）
   - §20 Tier 評分（已由 W3-W4 完成，直接 mapping：[F]→S, [O]→A, [I]→B, [U]→C）

2. **保留項目**：
   - §2 Outline（仍要拆 5-7 段）
   - §3 Section（仍要逐 section 確認）
   - §4 Slide build（仍要逐張視覺化）
   - §5 Visual review 10-checkpoint
   - §6 Speaker script（從母稿 §14「可用於講稿的例子」延伸）

3. **強化項目**：
   - 母稿 §17「教授可能提問」自動進 backup slide / Q&A slide
   - 母稿 §11「需要補充的資料」變成 P3 Section 階段的「待補空格」提示

### 4.3 handoff 時點

母稿 C6 確認後，AskUserQuestion：

```
母稿已完成（17 段，N 字）。下一步？
- (A) 繼續銜接 myPPT 產 PPT（自動帶入母稿）
- (B) 母稿留檔，下次再用 /myPPT --from-draft <path> 啟動
- (C) 母稿即終點，不產 PPT（如僅 internal status check）
```

## 5. 與其他 skill 的關係（不需整合）

| Skill | 關係 |
|-------|------|
| structured-tech-report | 平行 — 處理單一工程改動，不在週期性 cadence |
| review-evidence | 上游工具 — W1 raw data 收集可呼叫 |
| provenance-tier-audit | 上游工具 — W1 可呼叫做 evidence ledger sanity check |
| confirmation-protocol | 規範來源 — C0-C6 對應 Hard Gate / Gate / Review 級別 |
| doc-standards | 規範來源 — 母稿 .md 落點符合命名規範 |

## 6. 部署順序（4 階段）

```
Phase 1：先審所有 4 份大綱（本檔 + SKILL_OUTLINE + PLAYBOOK_OUTLINE + SUPPORTING_FILES_OUTLINE）
        ↓
Phase 2：用戶決定整合方案（選項 X / Y / Z）
        ↓
Phase 3：依選項落地實際檔案：
        - 選項 X：升級既有 weekly-report skill（內部新增 W1-W7）
        - 選項 Y：新建 weekly-master-draft skill 並調整 weekly-report description
        - 選項 Z：替換 weekly-report
        ↓
Phase 4：myPPT skill 加 --from-draft 接口（在 myPPT prompts/outline_confirm.md 增加母稿讀取分支）
```

## 7. 大綱審查重點

1. **整合方案決策**：選 X / Y / Z 哪一個？我傾向 X（避免 keyword 衝突），但需先看 weekly-report 既有實作能否支援 W1-W7 流程
2. **母稿輸出路徑**：`InterSubMod/docs/weekly_reports/YYYYMMDD/` vs 既有 weekly-report 路徑（請告知既有路徑）
3. **handoff prompt 預設**：(A)(B)(C) 三選項是否足夠？
4. **Tier mapping 規則**：[F]→S / [O]→A / [I]→B / [U]→C 是否合理？還是另設 tier？
5. **myPPT --from-draft 介面**：是否同意加這個 flag 或改用 frontmatter 自動偵測？
6. **selectivity**：weekly-master-draft 是否要支援「非週期性」的 advisor briefing（如 ad-hoc PI 報告）？還是嚴格限定週期性？
