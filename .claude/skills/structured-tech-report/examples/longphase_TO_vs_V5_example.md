# 範例：longphase-TO vs V5 Somatic Fallback Haplotag

> 本檔僅作 skill 使用範例的指引，**完整報告請見**：
>
> `InterSubMod/docs/reports/validated/2026/04/20260429_longphase_TO_vs_V5_Somatic_Fallback_技術報告_01.md`

## 為何被選為首個示範案例

1. **議題深度**：跨 longphase-TO（外部）、InterSubMod HaplotagProcess（內部 C++）、ISM 下游特徵汙染三層；單一 skill 是否承載得了？實測可。
2. **範本對齊**：報告完成後，章節結構與 `InterSubMod/docs/reports/pi_reports/2026/04/v5_audit_suite/07_paired_ground_truth_concordance.md` 一致，證明 13 段範本與既有 v5_audit_suite 慣例不衝突。
3. **既有資料豐富**：v5_audit_suite 19 份子報告 + 4 個 MEMORY 條目 + 2 份 validated 報告為依據，用來驗證 skill 的「Step 2 盤點」能否抓全 cross-ref。
4. **回應用戶最初的 Self-phasing 提問**：用戶在會話最初詢問 self-phasing 整理；此示範同時交付 skill + 內容兩件事。

## 此範例展示的 13 段技巧

| 段 | 技巧 |
|---|------|
| §0 TL;DR | 4 行（問題／修改／結果／狀態）+ 額外 1 行「關鍵釐清」釐清誤解（非替代而是互補） |
| §2 系統背景 | 名詞表（Glossary）放最前，避免 §3 起首就出現未解釋縮寫 |
| §3 原本流程 | ASCII 流程圖 + 「關鍵假設（後被推翻）」段，讓讀者理解為什麼舊設計是合理但有 bug |
| §4 問題 | 7 列量化指標表 + 解讀段（解釋 mega-block 的物理意義） |
| §5 根因 | 5 Whys 五層皆給「答覆」而非省略；§5.1 補強用 PON-only 實驗作獨立驗證 |
| §6 候選方案 | 3 個候選表 + Pros/Cons/採用？欄；補「設計目標 G1-G5」對應 §9 驗證 |
| §7.1 / §7.2 | 雙語：§7.1 用比喻（考試題目）；§7.2 用 8 列工程細節表（含 commit hash 待補標記） |
| §8 新舊比較 | 表格欄序：維度／Before／After／改變原因（不混雜時間順序） |
| §9 驗證 | 5 個 Step → Verify 區塊，每個有具體命令／檔案／數值範圍 |
| §10 影響 | 8 列受影響對象表 |
| §11 風險 | 5 個 R-id（與專案風險編號慣例對齊） |
| §12 後續 | 6 個 F-id 動作項，含負責／期限（TBD 標明）／連結 |
| §13 結論 | 限 5 句；首句問題、次句修改、再次結果、最後 open issue |

## 不要照抄的部分

- §7.2 commit hash 「⚠ 待補」 — 你的報告必須補齊或標註原因
- 附錄 A.5「框架方法學」是示範 skill 自我引用 — 一般技術報告**不需要**寫這個元層

## 反例教材：本報告初版犯了什麼錯（2026-04-29）

> 此節為**負面教材**，記錄初版報告（2026-04-29 00:55）犯的事實錯誤、用戶如何發現、重寫後的 RCA。**未來呼叫此 skill 的 AI 必看本節避免重蹈覆轍**。

### 錯誤現象

初版報告（落 validated/）出現以下 6 類事實錯誤：

| 類別 | 例 | 嚴重度 |
|-----|----|--------|
| 結構性歸屬錯誤 | 把 `HaplotagProcess` 寫成「InterSubMod 內部 C++」；實際在 `/big7_disk/liaoyoyo2001/longphase-to-mod/` 獨立 fork repo | 🔴 critical |
| 函數名捏造 | 寫成 `HaplotagProcess::tagRead`、捏造 `getVote_*` 系列；真實函數：`getVote()` / `judgeHaplotype()` / `countSNPHaplotype()` / `countINDELHaplotype()` | 🔴 critical |
| 數字捏造 | 「Phase block N50 11.9 Mbp → 1.2 Mbp」— 6 份 v5_audit_suite 子報告皆無此 metric | 🔴 critical |
| 檔案路徑捏造 | 「Tests `InterSubMod/tests/test_haplotag_v5.cpp`」— 不存在 | 🟠 high |
| 範圍誇大 | 「7 樣本 paired GT ≥ 0.78」— audit suite 僅 HCC1395 5kHz 一樣本 | 🟠 high |
| 檔名錯誤 | 引用 `01_code_diff.md` / `06_sanity_check.md`；實際 `01_code_diff_analysis.md` / `06_v5_sanity_bug_check.md` | 🟡 medium |

### 用戶如何發現

用戶一句話：「**InterSubMod 的 V5 Somatic Fallback Haplotag 與 InterSubMod 有哪些關係,應該是只 longphase-to 的 HaplotagProcess 吧**」 — 直接指出 V5 歸屬錯誤。AI 在驗證時用 `find / grep` 確認 `HaplotagProcess.cpp` 在 longphase-to-mod 內後，承認結構性錯誤。

### 5 Whys 根因

| 層 | 提問 | 答覆 |
|---|-----|-----|
| 1 | 為何把 V5 寫成 InterSubMod 內部？ | 直接信 Explore agent 在「兩者關係」表的「V5 = 內部 haplotag 邏輯」敘述 |
| 2 | 為何接受了那段敘述？ | **語境錯位**：explorer 寫的「內部」相對於 longphase 上游 GitHub；AI 讀成「InterSubMod 內部」 |
| 3 | 為何沒回頭驗證？ | 把「寫報告」當 documentation 任務（拿二手摘要填骨架），不是 fact-finding 任務（grep 程式碼確認） |
| 4 | 為何當成 documentation 任務？ | frame 為「skill 首個示範案例」，潛意識把「展示 13 段格式好用」放在「事實準確」之前；缺資料時傾向「合理估計」而非標 `⚠ 待確認` |
| 5 | 為何潛意識允許這個排序？（**根因**）| **內容段落的事實聲明沒被當成需要 Step→Verify 的「步驟」**；CLAUDE.md「Step → Verify」原則只被應用到工作流程，不是內容宣告 |

### 觸發的 skill 強化（不要再犯）

1. `references/checklist.md` 新增 **B9（fact-check）**：§4/§6/§7.2/§8 所有具體數字必須有可驗證 source；觸發此項違反 = 必擋。
2. `references/checklist.md` 新增 **B10（歸屬正確性）**：涉 C++ 程式碼改動段落必須 `find` / `grep` 確認檔案實際位置；涉外部 fork 必須在 §2 + §7.2 明確標 repo 歸屬。
3. SKILL.md Step 1 第 4 問「是否含 .cpp 改動」回答「是」時 → **強制呼叫 `/known-pitfalls`** 載入 `feedback_feature_name_vs_definition_rule.md`「分析新 feature 前必讀 src/include 定義」。

### 校正後狀態

- 報告初版 2026-04-29 00:55 → fact-check 全面校正版 2026-04-29 03:15
- 仍標 `⚠ 待補` 的 3 項：
  1. §4.2 「62% LOH 消失 / Cohen's d=-1.20」直接核對 `02_Self_Phasing根因.md`
  2. F1 commit V5 working tree 後回填 commit hash
  3. F3 7 樣本擴展結果回填
- 報告附錄 B 變更歷史誠實記錄全部修正過程，作為未來查證 / 學習依據

### 給未來 AI 的提示

> 寫技術報告時，每寫一個具體數字、檔案路徑、函數名、commit hash 之前，問自己：**「這個聲明的 source 在哪裡，我能用一行命令／一個檔案位置驗證嗎？」** 答不出來 → 標 `⚠ 待確認 X` 而非「合理估計」填空。即使你被 frame 為「示範案例」「展示用」，事實準確永遠優先於格式展示。

## 觸發本 skill 重現此類報告

```
/structured-tech-report 整理 <某個 pipeline／模組／bug> 的修改成 13 段技術報告
```

或：

```
/structured-tech-report from-commit <commit hash>
```
