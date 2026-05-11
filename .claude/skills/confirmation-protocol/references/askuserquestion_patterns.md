# Hard Gate AskUserQuestion 範本

5 個 Hard Gate 場景的 AskUserQuestion 結構化問題模板。模型呼叫 AskUserQuestion 工具時應參照此處對應 §N。

## §1. 刪除/搬移檔案

**何時觸發**：用戶或模型決定要 `rm` / `mv` 任何檔案（含 output/、docs/、research/ 內任何路徑）。

**AskUserQuestion 結構**：
```
question: 「確認 [刪除|搬移] 以下 N 個檔案？\n\n<列出前 5 個 + 總數>」
header: 「刪檔確認」
options:
  - label: "✅ 確認執行"
    description: "依列表 [刪除|搬移] 所有檔案，無 backup"
  - label: "📦 改為 archive"
    description: "搬移到 .archive/ 子目錄而非 rm，保留可救回"
  - label: "❌ 取消"
    description: "中止本次操作，保留所有檔案"
multiSelect: false
```

**後續邏輯**：
- ✅ → 執行 + append .claude/state/destructive_ops.jsonl audit log
- 📦 → mv 到 archive 目錄 + audit log
- ❌ → 中止；不寫 audit log

---

## §2. C++ commit

**何時觸發**：用戶要 commit 含 `.cpp` / `.hpp` / `.h` 改動，或 PreToolUse hook 觸發 git commit 攔截後的後續決策。

**AskUserQuestion 結構**：
```
question: "C++ 修改 commit 前確認：已編譯通過 + 已跑 tests？\n\n改動檔案：<list .cpp/.hpp 檔案>\n編譯狀態：<from .needs_compile flag>"
header: "C++ commit"
options:
  - label: "✅ Commit 現在"
    description: "編譯+測試已通過，直接 commit"
  - label: "🔄 先跑 batch test"
    description: "執行 ./scripts/run_batch_vcf_analysis.sh 後再決定"
  - label: "❌ Rollback 改動"
    description: "git checkout 還原 C++ 檔案，不 commit"
multiSelect: false
```

**後續邏輯**：
- ✅ → git commit 執行
- 🔄 → 觸發 batch test，等結果再回到本決策
- ❌ → git checkout 對應檔案

---

## §3. 研究方向 NO-GO 判定

**何時觸發**：模型/用戶要把某研究方向標記為 NO-GO（寫入 evidence_ledger.jsonl verdict=NO-GO）。

**AskUserQuestion 結構**：
```
question: "確認研究方向「<direction_name>」判定為 NO-GO？\n\n依據：<簡述 1-2 行 evidence>\n影響：撤回所有依賴此方向的 hypothesis（共 N 個）"
header: "NO-GO 判定"
options:
  - label: "✅ NO-GO 確定"
    description: "寫入 evidence_ledger verdict=NO-GO，撤回依賴假說"
  - label: "🔄 再試 1 cycle"
    description: "保留方向 active，再跑 1 個 cycle 收集數據"
  - label: "❌ 撤回判定"
    description: "本次不做 NO-GO 決定，回到 active 狀態"
multiSelect: false
```

**後續邏輯**：
- ✅ → append evidence_ledger + 更新 hypothesis_queue dependencies
- 🔄 → 不寫 ledger，建立新 cycle
- ❌ → 維持現狀

---

## §4. 覆寫 evidence_ledger / MEMORY.md 既有記錄

**何時觸發**：模型要修改 evidence_ledger.jsonl 既有 entry 的核心欄位（verdict / tier / hypothesis_id），或修改 MEMORY.md 既有 entry 的核心結論。注意：append-only 操作（新增 entry）不在此列。

**AskUserQuestion 結構**：
```
question: "確認 [覆寫|修改] 既有記錄？\n\n檔案：<file>\n原內容：<key fields>\n新內容：<diff>"
header: "覆寫紀錄"
options:
  - label: "✅ Append + retract 舊"
    description: "新 entry 標 supersedes，舊 entry 標 retracted（保留歷史）"
  - label: "🔄 只 amend 新欄位"
    description: "增加新欄位但不動既有欄位（最保守）"
  - label: "❌ 取消"
    description: "不修改"
multiSelect: false
```

**後續邏輯**：
- ✅ → append + 在舊 entry 加 retracted_by 欄位
- 🔄 → JSON merge 新欄位到既有 entry
- ❌ → 中止

---

## §5. 遠端 push (git push / PR)

**何時觸發**：用戶或模型要 `git push` / `gh pr create`。

**AskUserQuestion 結構**：
```
question: "確認 push 到遠端？\n\n分支：<branch>\n上游：<upstream>\nN commits ahead，包含：<前 3 commit subjects>"
header: "Remote push"
options:
  - label: "✅ Push"
    description: "直接 push，commits 立即可見於遠端"
  - label: "🔄 改 PR draft"
    description: "建立 draft PR（不 ready-for-review），給時間 review 再升級"
  - label: "❌ 取消"
    description: "保持 local，不 push"
multiSelect: false
```

**後續邏輯**：
- ✅ → git push
- 🔄 → gh pr create --draft
- ❌ → 中止

---

## 通用約定

1. **option label 字數限制**：≤ 1-5 詞 + 1 emoji（per AskUserQuestion tool spec）
2. **首選項標 (Recommended)**：若有明確推薦，第一項加 "(Recommended)" 後綴
3. **「Other」自動加入**：AskUserQuestion 系統會自動提供 Other 讓用戶自由輸入，不需在 options 列
4. **multiSelect**：Hard Gate 永遠 `false`（單選），避免「同時選擇衝突動作」
5. **不能用 AskUserQuestion 來問「Plan 通過嗎？」**：plan 用 ExitPlanMode 工具，不是 AskUserQuestion
