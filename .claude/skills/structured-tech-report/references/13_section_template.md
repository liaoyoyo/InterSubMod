# 13 段式技術報告範本（structured-tech-report）

> 複製本檔內容到目標 .md，逐段填入。資料不足處標 `> ⚠ 待確認：<具體問題>`，不要因此停工。
> 路徑列給用戶時必以 `InterSubMod/...` 前綴。

---

```html
<!--
build_date: YYYY-MM-DD HH:MM
agent: <填入 session id 或 'manual'>
status: draft | validated | finalized
report_class: bug-fix | pipeline-change | methodology-improvement | architecture-decision
audience: PI | engineer | reviewer | ai-session
inputs:
  - InterSubMod/path/to/source1.md
  - InterSubMod/src/core/File.cpp:L120-L188
outputs:
  - InterSubMod/docs/reports/.../<this_report>.md
  - figures/<section>/figXX_*.png
related_memory:
  - memory/project_<topic>.md
verdict: pending | POSITIVE | NEGATIVE | CONDITIONAL
last_verified: YYYY-MM-DD
-->
```

# {報告主題}

> **TL;DR**（3-5 行）｜寫給「沒空看全文，只想知道發生什麼／改了什麼／結果如何」的人。
> 範例格式：
> - **問題**：……（量化）
> - **修改**：……（一句話）
> - **結果**：…F1=0.71xx，AMB% xx→xx%……
> - **狀態**：validated / pending / conditional

---

## 1. 報告目的

**為什麼寫這份報告**（限 3-5 句）：
- 觸發背景（哪個 PR／哪份觀察／哪個 issue 啟動）
- 預期讀者收穫（讀完能回答什麼問題）

## 2. 系統背景

**此模組／流程原本負責什麼**：
- 一句話定位
- 與上下游的關係（誰餵它？它餵誰？）
- 關鍵詞表（首次出現的縮寫展開：例 `ISM = Inter-Subclonal Methylation`）

## 3. 原本流程（Before）

用條列或流程箭頭：

```
Step 1: <input> ─→ <tool A> ─→ <intermediate>
Step 2: <intermediate> ─→ <tool B> ─→ <output>
```

> 必要時加 pseudo-code 或 sequence diagram；嚴禁丟整段原始碼，重點是流程而非實作。

## 4. 問題描述

**情境 → 現象 → 影響**（必量化）：

| 維度 | 內容 |
|-----|-----|
| 觸發情境 | 在什麼樣本／什麼參數／什麼版本下出現 |
| 錯誤現象 | 具體可觀察（log、圖表、metric） |
| 影響範圍 | 哪些 dataset／feature／使用者受影響 |
| 量化指標 | 例：62% LOH 消失、F1 下降 0.05、AMB% 17.5% |

**圖表佐證**（若有）：
```markdown
![Fig 4a — <標題>](figures/04_problem/fig4a_<topic>.png)
```

## 5. 根本原因

依「報告類別」二選一：

### 5.A 5 Whys（適合 bug fix／pipeline 變更）

1. 為何 X？→ 因為 Y
2. 為何 Y？→ 因為 Z
3. 為何 Z？→ 因為 W
4. 為何 W？→ 因為 V
5. 為何 V？→ 因為 U（**根因**）

### 5.B Ishikawa Fishbone（適合多因素）

```
                Methods            Measurement
                   |                    |
                   ↓                    ↓
        ─────────────────────────────────────→ 問題
                   ↑                    ↑
                   |                    |
                Material              Machine
```

### 5.C Decision Drivers（適合架構決策）

| Driver | 描述 | 權重 |
|--------|------|-----|
| ……   | ……  | …  |

## 6. 修改方向

列**所有**候選方案（不只採用的）：

| 候選 | 內容 | Pros | Cons | 採用？ |
|-----|------|------|------|-------|
| A | …… | …… | …… | ❌（理由：……） |
| B | …… | …… | …… | ✅ |
| C | …… | …… | …… | ❌（理由：……） |

**設計目標**（success criteria）：
- ……
- ……

## 7. 修改內容（雙語）

### 7.1 非工程版（給 PI／外部讀者）

用白話 3-5 句說明改了什麼。**禁用** file path、function name、commit hash、API 名。

> 範例：「在貼 HP 標籤時，碰到體細胞變異附近的 read 改投信心票；票數不足就標『不確定』而非硬指派。」

### 7.2 工程版（給工程同事）

具體技術調整：

| 維度 | 變更 |
|-----|-----|
| Input validation | …… |
| Parser logic | `InterSubMod/src/core/File.cpp:L120-L188`（commit `<hash>`） |
| Database schema | …… |
| API signature | `before(...) → after(...)` |
| Error handling | …… |
| Logging | …… |
| Tests | `InterSubMod/tests/test_<name>.cpp`（new / modified） |

> 涉 .cpp 改動：commit hash **必填**；如尚未 commit 標 `> ⚠ 待補 commit hash`。

## 8. 新舊流程比較

| 維度 | 原本流程 | 新流程 | 改變原因 |
|-----|---------|--------|---------|
| Input | …… | …… | …… |
| Step 1 | …… | …… | …… |
| Output | …… | …… | …… |
| 性能 | xx ms | yy ms | …… |
| 結果指標 | F1=xx | F1=yy | …… |

## 9. 驗證方式

採 **Step → Verify** 格式（每步驗證必須**外部可觀察**）：

```
1. <步驟> → 驗證: <具體可觀察的檢查方式：命令／檔案／數值範圍>
2. <步驟> → 驗證: <…>
3. <步驟> → 驗證: <…>
```

**測試矩陣**：

| 案例 | 輸入 | 預期 | 實際 | Pass? |
|-----|------|------|------|-------|
| 正常 | …… | …… | …… | ✅ |
| 邊界 | …… | …… | …… | ✅ |
| 異常 | …… | …… | …… | ✅ |
| 回歸 | 修改前 baseline | F1 ≥ baseline-0.01 | …… | ✅ |
| 效能 | …… | …… | …… | ✅ |

> 涉統計驗證 → 必呼叫 `/auc-confound-guard` 三關（within-group OLS + AF-bin × + permutation）。

## 10. 影響範圍

| 受影響對象 | 影響 |
|-----------|-----|
| 使用者 | …… |
| 資料格式 | …… |
| API | …… |
| 資料庫 | …… |
| 部署 | …… |
| 舊資料 | 是否需重跑？哪些版本範圍？ |
| 下游分析 | 哪些 ISM／VCF feature 計算需重跑？ |
| 維護流程 | …… |

## 11. 風險與限制

- 仍待驗證的情境：……
- 已知副作用：……
- 適用邊界：……
- 與既有功能的相容性風險：……

## 12. 後續工作

對應 SRE Action Items 格式：

| ID | 動作 | 負責人 | 期限 | 連結 |
|----|------|-------|------|------|
| F1 | 補測試 | …… | …… | …… |
| F2 | 更新文件 | …… | …… | …… |
| F3 | 監控 log | …… | …… | …… |
| F4 | 重跑舊資料 | …… | …… | …… |
| F5 | rollback plan 撰寫 | …… | …… | …… |

## 13. 結論

3-5 句，回答：
- 問題是什麼
- 為什麼改
- 改了什麼
- 目前結果如何
- 後續要注意什麼

---

## 附錄

### A. 引用清單

- 文件：`InterSubMod/docs/...`
- 程式碼：`InterSubMod/src/core/...:L<line>`
- 數據：`InterSubMod/output/.../*.tsv`
- 圖表：本報告 `figures/` 子目錄
- 外部論文：[作者 et al., year](URL)
- MEMORY：`memory/project_<topic>.md`

### B. 變更歷史

| 日期 | 變更 | 作者 |
|-----|------|-----|
| YYYY-MM-DD | 初稿 | …… |
| YYYY-MM-DD | 補 §9 7 樣本驗證 | …… |
