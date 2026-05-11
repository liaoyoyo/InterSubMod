---
id: ism-kb-00-governance-confirmation-protocol
name: "確認時機協議（速查）"
description: "KB 速查：Hard Gate / Gate / Review / FYI 4 級暫停分類、互動/全自動模式切換、20 個典型場景對照、Opus 4.7 subagent 觸發規則。權威定義在 .claude/skills/confirmation-protocol/SKILL.md。"
status: active
last_verified: 2026-04-22
content_nature: frozen-decision
doc_type: reference
verified_scope: "索引對應 .claude/skills/confirmation-protocol/SKILL.md（219 行權威版）"
related_ids:
  - ism-kb-00-governance-index
  - ism-kb-00-governance-query-protocol
  - ism-kb-00-governance-new-info-protocol
  - ism-kb-06-workflows-cpp-change-pdd
  - ism-kb-00-governance-think-before-code
tags: [governance, confirmation, hard-gate, mode, autonomy, opus-4.7]
canonical_paths: [00_governance/09_confirmation_protocol.md]
alias_paths: []
---

# 確認時機協議（速查）

- 一句結論：AI 執行前回答「不可逆？→ 🔴 必停；否則查 2D 矩陣（影響×信心）」；本頁速查，權威在 skill
- 適用對象：AI agent 決策時、用戶切換執行模式時
- 可直接執行命令（驗證日期：2026-04-22）：
  ```bash
  cat /big7_disk/liaoyoyo2001/InterSubMod/.claude/skills/confirmation-protocol/SKILL.md
  ```

---

## 🎚️ 兩種執行模式

| 模式 | 觸發詞 | 行為 |
|------|--------|------|
| **互動模式**（預設） | `互動模式`、`interactive` | Gate 和 Review 暫停等確認 |
| **全自動模式** | `全自動`、`auto`、`自動跑到底` | 只 Hard Gate 停；Review/FYI 自動通過，結束出報告 |

**範圍**：模式切換只在當前對話有效；新對話預設互動模式。

---

## 🚦 4 級暫停分類

| 級別 | 互動模式 | 全自動模式 | 用途 |
|------|---------|-----------|------|
| 🔴 **Hard Gate** | 暫停 | **仍暫停** | 不可逆操作（任何模式都必停） |
| 🟠 **Gate** | 暫停 | 自動過（預設選項） | 重要決策點 |
| 🟡 **Review** | 暫停 | 自動過 | 中間產出展示 |
| 🟢 **FYI** | 顯示 | 靜默 | 低風險通知 |

---

## 🔴 Hard Gate 清單（不可逆，永不跳過）

**CLAUDE.md 原生硬性規則（1-3）**：
1. 刪除 / 搬移任何檔案（含 `output/` 自動產出）
2. **C++ 原始碼 commit**（對應 pre-commit hook `exit 2`）
3. **研究方向 NO-GO 判定**（寫入 evidence_ledger `verdict=NO-GO`）

**KB v0.2+ 擴充（4-6）**：
4. **覆寫** `evidence_ledger.jsonl` 或 `MEMORY.md`（append 不在此列）
5. **Git push 到受保護分支**（`main` / `develop` / `main`/`master`）或 **force push**（`push --force`）；**建立 PR**；push 到自己的 feature branch 為 🟡/🟢
6. 啟動 >10 min 計算且用戶**未明示**

**長計算有條件例外**：用戶當輪明示含長度指示（「跑全量」「平行 7 樣本」「>10min」「過夜」）→ 🟢 告知直接執行。若只說「看看 HCC1395 benchmark」（無長度詞）→ 仍 🔴 暫停確認時間。

---

## 📊 2D 暫停判斷矩陣

```
影響度 ＼ 信心度    高信心        中信心        低信心
───────────────────────────────────────────────────
低（<10min 可逆）   🟢 告知       🟡 列假設     🟠 節點暫停
中（10min-1h）      🟡 列假設     🟠 節點暫停   🔴 立即暫停
高（>1h / 影響結論）🟠 節點暫停   🔴 立即暫停   🔴 立即暫停
```

**可逆性 override**：任一不可逆（見 Hard Gate 清單）→ 永遠 🔴

### 三維度定義
- **影響**：低（可逆/局部/<10min）/ 中（多檔案/重跑/10min-1h）/ 高（跨模組/影響結論/>1h）
- **信心**：高（先例+單一解讀）/ 中（2-3 解讀偏向明確）/ 低（多解讀/無先例/規格缺項）

---

## 💬 一行告知格式

`[決策]（影響: X, 信心: Y, 理由: 一句）`

**範例**：
- 🟢 `以 within-group OLS 殘差化 AF（影響: 低, 信心: 高, 理由: O12 已驗證 pooled OLS collider bias）`
- 🟡 `選 NHD 距離（影響: 低, 信心: 中, 理由: 與 O9 baseline 對齊，可直接比較）`

---

## 📋 25 個場景對照表

> **來源合併說明**：本表合併自 `.claude/skills/confirmation-protocol/SKILL.md` §「跨場景快速查表」(12 項) + §「各場景暫停類型」(12 項) + `CLAUDE.md` §「高影響場景清單」(5 類) + KB v0.4 Review 補強（CMakeLists/requirements/hooks 編輯）。

| 場景 | 級別 |
|------|------|
| 加圖表標註 / 改 matplotlib 參數 | 🟢 |
| 選擇統計方法且有先例 | 🟢 |
| 改單一 Python 分析腳本邏輯 | 🟡 |
| 新增 Python 特徵欄位 | 🟡 |
| 微調 C++ 閾值（已有 baseline） | 🟠 |
| 新增距離度量算法 | 🟠 |
| 決定是否 escalate 至 full benchmark | 🟠 |
| 研究重點排序 | 🔴（高影響） |
| 假說選擇（新主題） | 🔴（高影響） |
| 統計方法切換（影響結論） | 🔴 |
| 「改進」/「優化」模糊指令 | 🔴（先問定義） |
| 多檔案重構 | 🔴（中-高） |
| 判定研究方向 NO-GO | 🔴（不可逆） |
| 刪除 output/ 下任何檔案 | 🔴（不可逆） |
| C++ 原始碼 commit | 🔴（不可逆） |
| 覆寫 evidence_ledger.jsonl | 🔴（不可逆） |
| 啟動 >10min 計算（用戶明示） | 🟢 |
| 啟動 >10min 計算（模型自判） | 🔴 |
| git commit（非 C++） | 🟢 FYI |
| git push 到 feature branch | 🟡 告知 |
| git push 到 main/develop 或 --force | 🔴 不可逆 |
| 修改 `CMakeLists.txt` | 🟠（影響 build 結構） |
| 修改 `requirements.txt` / `conda env` | 🟠（跨機器相依） |
| 編輯 `scripts/hooks/*.sh`（hook 本身） | 🔴（影響所有後續 tool use）|
| 修改 `.clang-format` | 🟠（影響 diff 大小） |
| 記憶更新（append to MEMORY.md） | 🟢 FYI |

**完整表 + 「與既有權限邊界的關係」** → `.claude/skills/confirmation-protocol/SKILL.md` §「各場景暫停類型」

---

## 🤖 Opus 4.7 Subagent 觸發規則

**預設行為**（不需明指）：
- 單檔案修改 → 模型直接做，**不 spawn subagent**
- ≤3 查詢探索 → 直接 Grep/Glob，不呼叫 Explore
- 單樣本測試 → 直接執行，不呼叫 parallel-benchmark

**必須明寫觸發語才啟用**：

| 情境 | 觸發語 | Agent |
|------|--------|-------|
| 跨樣本平行 benchmark | 「平行驗證 HCC1395 / COLO829 / H2009」 | `parallel-benchmark` |
| 深度程式碼探索（>3 查詢） | 「深度探索 XXX 實作鏈」 | `Explore` / `general-purpose` |
| PR / 大型修改審查 | 「code review 這次 commit」 | `pr-review-toolkit:code-reviewer` |
| 獨立第二意見 | 「獨立審查 / 第二意見」 | `feature-dev:code-reviewer` |
| 方法學深度分析 | 「啟動 methodology-audit」 | `methodology-reviewer` |

**例外（結構性必要，不跳過）**：`/cpp-change` Step 4.5 的 code-reviewer 呼叫

---

## 🎯 決策流程（每次決策依序檢查）

```
1. 這是不可逆操作嗎？（見 Hard Gate 清單）
   ├─ 是 → 🔴 立即暫停
   └─ 否 → 進入 2

2. 評估影響度（低 / 中 / 高）
3. 評估信心度（高 / 中 / 低）
4. 查 2D 矩陣 → 得出級別

5. 依級別行動：
   🟢 → 一行告知 + 繼續
   🟡 → 列 2-3 假設 + 繼續第一步
   🟠 → 完成當前節點 + 暫停報告
   🔴 → 立即停下報告 + 等確認
```

---

## 🔗 AI 自主權限邊界

**可自主**（不需暫停；需記錄理由）：
- 選擇分析腳本 / 統計方法
- 圖表格式、排版
- 非 C++ git commit
- 檔案命名、目錄組織
- 中間數據格式

**需暫停（Gate 以上）**：
- 研究方向、假說確認
- C++ 程式碼修改
- 結論判定（positive / NEGATIVE / NO-GO）
- 刪除 / 搬移檔案
- >10 min 計算（未明示時）

---

## 權威來源

- **完整定義**：[`.claude/skills/confirmation-protocol/SKILL.md`](../../.claude/skills/confirmation-protocol/SKILL.md)（219 行）
- **CLAUDE.md**：「確認時機協議」章節（硬性要求聲明）
- **相關 hooks**：[08_hooks_and_automation.md](08_hooks_and_automation.md)（C++ commit Hard Gate hook 實作）

---

## 相關

- Query protocol：[05_query_protocol.md](05_query_protocol.md)（查詢情境對照）
- Think Before Code：[10_think_before_code.md](10_think_before_code.md)（假設陳述 + 暫停矩陣深化）
- PDD 6 steps：[../06_workflows/07_cpp_change_pdd.md](../06_workflows/07_cpp_change_pdd.md)（Hard Gate 典型場景）
- Hooks 自動化：[08_hooks_and_automation.md](08_hooks_and_automation.md)
