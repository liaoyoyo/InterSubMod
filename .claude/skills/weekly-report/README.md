# weekly-report skill (v2)

> InterSubMod 每週研究週報整理工作流。
> 由 v1「直接產 PPTX」升級為 v2「先產母稿，再 handoff 給 pptx-build」。

## 解決什麼

- 每週把研究進展、學習、問題、資料佐證、後續計畫整理成 **可向教授報告的 17 段母稿**
- 透過 4 主線類型識別、4 層分類 [F]/[O]/[I]/[U]、5-7 個教授追問預測，確保 **邏輯嚴謹 + 不過度宣稱 + 不流水帳**
- 母稿產出後 4 選 handoff（立即接 pptx-build / 留檔 / 終點 / 加寫下週計畫），給用戶中斷空間

## 何時觸發

| 觸發語境 | 範例 |
|---------|------|
| 週報整理 | 「本週週報」「整理本週」「weekly report」 |
| 向教授彙報 | 「向教授報告」「PI 週彙報」「lab meeting 用」 |
| 研究進度回顧 | 「研究進度報告」「本週進展」 |

## 5 個 checkpoint 流程

```
W1 Raw Data 收集（git+Memory+KB+上週 Layer 3-4） → C0 確認
W2 4 主線類型識別（進展/問題/求協助/探索）       → C1 確認 ★ fast-track 必停
W3 內容 4 層分類 [F]/[O]/[I]/[U]
W4 重點排序 + 4 桶分流（PPT/講稿/備註/暫存）      → C2 確認
W5 邏輯紅旗檢查（過度宣稱/流水帳/教授視角缺）
W6 教授問答預測 5-7 個 + 預備回答                  → C3 確認
W7 17 段母稿產出（Layer 0-4 結構 + 17 段內部標籤） → C4 確認 ★ fast-track 必停
[handoff 4 選]
```

## 與全域 docs 連結

| 主題 | 文件 |
|------|------|
| 5 階段詳細方法論 | `playbook.md` §1-§16 |
| Layer 0-4 + 17 段 mapping | `references/LAYER_STRUCTURE.md` §B |
| [F]/[O]/[I]/[U] 規則 + 與 Tier 1/2/3 並用範例 | `references/LAYER_STRUCTURE.md` §C |
| 4 桶分流 5 維度評分 | `references/LAYER_STRUCTURE.md` §D |
| W1 raw data 自動收集細則 | `references/COLLECTION_PROTOCOL.md` |
| 母稿 frontmatter schema + handoff 4 選模板 | `references/HANDOFF_TO_PPTX_BUILD.md` |
| 4 主線類型 narrative skeleton | `templates/{progress|problem|advisor|new_direction}_focus.md` |
| 5 個 checkpoint 互動 prompt | `prompts/*.md` |
| 完整母稿範例 | `examples/master_draft_example.md` |

## 與其他 skill 關聯

- **下游接棒** → `pptx-build` skill（C4 後 handoff A 觸發；自動讀母稿跳過 main thesis 鎖定）
- **上游觸發**（從 myPPT 總入口 場景識別後委派）
- **上游工具** → `review-evidence` / `provenance-tier-audit`（W1 raw data 收集）
- **規範來源** → `confirmation-protocol`（C0-C4 對應 Hard Gate 級別）/ `doc-standards`（母稿命名）
- **平行不重疊** → `structured-tech-report`（單一工程改動 deep dive，非週期性）

## 母稿輸出位置

```
InterSubMod/docs/reports/validated/YYYY/MM/YYYYMMDD_週報_<主題>/
├── master_draft.md           # W7 → C4 主檔
├── next_week_plan.md         # 選 D 時加產
└── pptx_build/               # 選 A 觸發後 pptx-build 產出
```

## v1 → v2 主要差異

| 維度 | v1（已備份 SKILL.md.v1.bak）| v2 |
|------|------|------|
| 階段命名 | Phase 0-6 | W1-W7 |
| Checkpoint 數 | 7 個 | 5 個（合併）|
| 主線類型識別 | ❌ | ✅ 4 選 1 |
| 4 層分類 [F]/[O]/[I]/[U] | ❌ | ✅ 與 Tier 1/2/3 並用 |
| 教授問答預測 | ❌ | ✅ 5-7 個 |
| 過度宣稱紅旗 | ❌ | ✅ |
| 直接產 PPTX | ✅ | ❌（移到 pptx-build）|
| 17 段母稿格式 | ❌ | ✅ 內部標籤 |
| handoff 4 選 | ❌ | ✅ |

詳細 v2 升級理由與決策矩陣 → `InterSubMod/.claude/skills/weekly-master-draft/outlines/UPGRADE_PLAN_FOR_WEEKLY_REPORT.md`
