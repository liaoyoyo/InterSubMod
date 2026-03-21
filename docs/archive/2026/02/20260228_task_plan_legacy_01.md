<!--
建立時間: unknown (legacy)
metadata補齊時間: 2026-03-01 03:02
狀態: legacy
資料來源: 歷史文件補齊（內容未改寫）
驗證命令: python3 scripts/analysis/validate_docs_structure.py
關聯文件: docs/standards/20260228_文件命名與狀態管理規範_01.md
-->

# Task Plan: Git 變更整理與分批提交策略

<!--
建立時間: 2026-01-23
目標: 分析當前 Git 變更，根據開發目標分類並建立合理的分支和 commit 策略
處理範圍: 當前 develop 分支上的所有已暫存和未追蹤變更
-->

## Goal
將當前混雜的 Git 變更分類整理，根據功能目標分批建立 feature 分支，撰寫清晰的 commit message，並依序合併回 develop 分支。

## Context
- **當前分支**: `develop`
- **主要分支**: `main`
- **已暫存變更**: 13 個檔案，+3784/-53 行
- **未暫存變更**: 5 個檔案 (主要是刪除檔案)
- **未追蹤檔案**: 4 個新增的分析報告檔案

## Phases

### Phase 1: 分析變更分類 [complete] ✓
**Goal**: 將所有變更按照功能目標分類

**Actions**:
- [x] 檢視 git status 和 diff 統計
- [x] 解決 AD (add-delete) 檔案衝突
- [x] 讀取核心程式碼變更內容 (LabelTest, RegionProcessor, Stats)
- [x] 讀取文檔變更內容 (研究報告)
- [x] 確認刪除檔案的原因 (舊工具已過時)
- [x] 將變更分為 5 個邏輯群組

**Output**: `findings.md` 中的詳細分類表 + 執行順序

**Completed**: 2026-01-23

---

### Phase 2: 制定分支策略 [complete] ✓
**Goal**: 為每個變更群組設計分支名稱和提交順序

**Actions**:
- [x] 根據變更性質決定分支命名 (feature/docs/cleanup/chore)
- [x] 決定分支的依賴關係和合併順序
- [x] 設計每個 commit 的範圍和 message

**分支執行計劃** (詳見 findings.md):
1. `feature/multi-stage-hp-verification-2026-01-23` - 核心功能 (6 檔案)
2. `docs/methylation-significance-research-2026-01-23` - 研究文檔 (4 檔案)
3. `docs/multi-stage-verification-reports-2026-01-23` - 實作報告 (4 檔案)
4. `cleanup/remove-obsolete-f1-tool-2026-01-23` - 清理 (1 檔案)
5. `chore/update-claude-config-2026-01-23` - 設定 (1 檔案)

**Dependencies**: Phase 1

**Completed**: 2026-01-23

---

### Phase 3: 執行第一批變更 (核心功能) [pending]
**Goal**: 處理核心程式碼的功能增強

**Actions**:
- [ ] 建立 feature 分支
- [ ] 暫存相關檔案
- [ ] 撰寫 commit message (包含 Co-Author)
- [ ] 執行編譯測試
- [ ] 合併回 develop

**Dependencies**: Phase 2

---

### Phase 4: 執行第二批變更 (文檔) [pending]
**Goal**: 處理文檔和分析報告的新增

**Actions**:
- [ ] 建立 docs 分支
- [ ] 暫存文檔檔案
- [ ] 提交變更
- [ ] 合併回 develop

**Dependencies**: Phase 3

---

### Phase 5: 執行第三批變更 (清理) [pending]
**Goal**: 處理檔案刪除和設定更新

**Actions**:
- [ ] 建立 cleanup 分支
- [ ] 暫存刪除操作
- [ ] 提交變更
- [ ] 合併回 develop

**Dependencies**: Phase 4

---

### Phase 6: 最終驗證 [pending]
**Goal**: 確認所有變更已正確提交

**Actions**:
- [ ] 檢查 git status (應該是乾淨的)
- [ ] 檢視 git log 確認 commit 歷史合理
- [ ] 執行完整測試
- [ ] 更新 `progress.md` 總結

**Dependencies**: Phase 5

---

## Decision Log

| Decision | Rationale | Date |
|----------|-----------|------|
| 使用多個 feature 分支而非單一 commit | 變更涉及 3 個不同領域 (核心功能/文檔/清理)，分開處理更清晰 | 2026-01-23 |

## Errors Encountered

| Error | Attempt | Resolution |
|-------|---------|------------|
| - | - | - |

## Files Modified

| File | Purpose | Phase |
|------|---------|-------|
| task_plan.md | 專案規劃 | 1 |
| findings.md | 變更分析結果 | 1 |
| progress.md | 執行日誌 | 1-6 |

## Notes

- 遵循專案的 commit message 規範 (包含 Co-Authored-By)
- 每個分支在合併前都要確認編譯和測試通過
- 分支命名遵循 `feature/` 或 `docs/` 前綴
