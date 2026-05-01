---
title: weekly-master-draft playbook.md 大綱（待用戶審查）
date: 2026-05-01
status: outline-for-review
target_lines: 1000-1200
---

# playbook.md 大綱（§1-§16 主方法論）

## 章節結構（預估行數）

### §1 系統理念（≈ 60 行）
- 為何先做母稿再做 PPT（核心原則）
- Human-in-the-loop 三角色：AI 整理者 / AI 追問者 / 用戶判斷者
- 不應該：AI 不替用戶決定研究立場 / 不把推論寫成事實 / 不過度美化
- 應該：AI 扮演「週報整理助理」+「簡報規劃教練」
- 母稿 vs PPT vs 講稿：三者關係圖

### §2 W1 Raw Data 收集（≈ 80 行）
- 預期輸入來源：evidence_ledger.jsonl / docs/experiments/in_progress/ / git log / commit message / CURRENT_FOCUS.md
- AI 自動掃描清單（grep keyword 列表）
- 用戶補充輸入：本週實際做了什麼（口述 / 條列）
- C0 確認 prompt：列出找到的 raw data，請用戶補充缺漏 / 排除誤抓
- 輸出格式：時間軸 + 每日活動 + 產出 artifacts 清單

### §3 W2 主線類型識別（≈ 90 行）
- 4 種主線類型詳細定義（進展 / 問題 / 求協助 / 探索）
- 識別依據：raw data 屬性比例（產出多 → 進展；blocker 多 → 問題；不確定多 → 求協助；pilot 多 → 探索）
- 混合場景處理：以教授最關心的為主，其他降 sub-thread
- 一句話 main statement（≤ 30 字）強制鎖定
- C1 確認 prompt：4 選 1 + main statement 撰寫
- 範例：4 種主線各舉 1 個虛擬週報範例

### §4 W3 內容 4 層分類（≈ 120 行）
**核心方法論段落，最詳細的章節之一**

4 層定義（與 myPPT Tier S/A/B/C/D 不同維度）：
- **[F] Fact 確定事實**：有具體 source（檔案 path / commit hash / output csv）
- **[O] Observation 初步觀察**：有結果但 N 不足或未驗證
- **[I] Inference 合理推論**：根據資料推測，需保留語氣
- **[U] Unconfirmed 待確認**：有疑問或不確定

每筆素材標籤決策樹（mermaid）：
1. 有具體檔案 source？ → 是→F 否→O
2. F：N 樣本 ≥ confidence threshold？ → 否→降 O
3. O/I：是否單一樣本？ → 是→I 否→看其他證據
4. U：標出「需 X 才能確認」

範例：v3 self-phasing 的 7 筆素材如何分類

C2 確認 prompt：每筆素材列當前 AI 標籤 + 用戶可改

### §5 W4 重點排序與 4 桶分流（≈ 110 行）
- 排序評分標準（5 維度）：
  1. 研究重要性（1-5）
  2. 資料證據強度（F=5/O=3/I=2/U=1）
  3. 教授可能關心程度（1-5）
  4. 是否影響下週計畫（1-5）
  5. 是否適合簡報呈現（1-5）
- 加總分數 → 4 桶分流：
  - 18-25：PPT（slide 上）
  - 13-17：講稿（speaker note）
  - 8-12：備註（appendix / 口頭預備）
  - <8：暫存（不放本次報告）
- 每桶上限規則（避免 PPT 過載）：PPT ≤ 8 筆 / 講稿 ≤ 15 筆
- C3 確認 prompt：列每筆素材 + 預設桶 + 用戶可調

### §6 W5 邏輯檢查（紅旗）（≈ 80 行）
- 過度宣稱紅旗 6 條 + 範例
- 流水帳紅旗 4 條 + 範例
- 教授視角缺紅旗 5 條 + 範例
- AI 自動掃描方式（keyword 比對 + 結構檢查）
- C4 確認 prompt：列出觸發紅旗的句子 + 改寫建議

### §7 W6 教授問答預測（≈ 80 行）
**核心方法論段落之二**

預測來源：
- 本週主線類型 → 對應問題模板
- 4 層分類中的 [I] / [U] → 教授會追問
- 「過度宣稱」紅旗未完全修正 → 教授會質疑
- 與過去報告的不一致 → 教授會記得

5-7 個追問模板（依主線類型）：
- 進展型：「為什麼是這個方法？」「跟既有 baseline 比？」「下一步邊界？」
- 問題型：「根因確認了嗎？」「短期 workaround？」「影響範圍？」
- 求協助型：「你傾向哪個？」「不選的代價？」「怎麼決定？」
- 探索型：「pilot 結果可信嗎？」「scale up 風險？」「資源？」

每個追問必須有「預備回答 1 段」（含 evidence 引用）

C5 確認 prompt：列預測追問 + 用戶補充教授個人偏好

### §8 W7 17 段母稿產出（≈ 100 行）
**最終輸出格式定義**

17 段標準格式（基於您的敘述）：
```
1. 本週報告主線（一句話）
2. 本週一句話重點
3. 已確認內容 [F]
4. 初步觀察與合理推論 [O][I]
5. 待確認內容 [U]
6. 不建議放入 PPT 的內容（暫存）
7. 報告重點優先順序
8. 建議報告順序
9. 建議 PPT 模板類型（指向 myPPT 6 模板）
10. 建議投影片架構
11. 需要補充的資料
12. 需要製作的圖表
13. 需要補充的定義或解釋
14. 可用於講稿的例子
15. 暫存紀錄
16. 下一步行動清單
17. 教授可能提問與回答準備
```

每段填寫規則 + 範例 + 字數建議

C6 確認 prompt：列母稿全文 + 用戶逐段批准 / 修改

### §9 §10 與 myPPT 的銜接設計（≈ 60 行）
- handoff 協議（output schema）
- 母稿 frontmatter 規範
- myPPT 讀取母稿後跳過的 §（§1 部分項目）
- myPPT 讀取母稿後保留的 §（§2 outline / §3 section / §4 slide / §5 visual / §6 speaker）
- AskUserQuestion 銜接模板

### §10 範例：1 份完整母稿（≈ 80 行）
- 虛構場景：某週 InterSubMod 進展（包含 1 個 [F]、2 個 [O]、1 個 [I]、1 個 [U]）
- 完整跑 W1-W7
- 最終 17 段母稿全文 + 銜接 myPPT 後的 outline

### §11 反例（≈ 60 行）
列 5-6 個典型錯誤：
- 把 [I] 寫成 [F]（過度宣稱）
- 流水帳堆 8 件事
- 母稿無教授問答段
- 主線類型混雜（強進展+強問題）
- raw data 收集遺漏（git log 沒掃）

### §12 互動 protocol 模板（≈ 60 行）
- W1-W7 對應 prompts/ 檔案位置
- 每輪 ≤ 5 個問題的設計理由
- AskUserQuestion 4 options 模板（針對 7 個 checkpoint）

### §13 fast-track 模式（≈ 40 行）
- 全自動觸發條件
- 跳過項目清單
- 仍強制項目清單

### §14 與既有 skill 互引用（≈ 40 行）
- weekly-report skill 的關係（W1 升級為 weekly-report）
- structured-tech-report skill 的關係（單一工程改動 deep dive 不適用本 skill）
- evidence-review skill 的關係（W1 自動 grep evidence_ledger.jsonl）

### §15 輸出檔案管理（≈ 30 行）
- 路徑規範：`InterSubMod/docs/weekly_reports/YYYYMMDD/master_draft_v{N}.md`
- frontmatter 規範
- 版本控制：v1 / v2 / v_final 命名
- 與 myPPT 的 PPT 輸出路徑分開

### §16 success criteria（≈ 30 行）
- 母稿是否能讓教授獨立讀懂
- 5-7 個教授追問是否覆蓋 80% 實際追問
- [F]/[O]/[I]/[U] 標籤是否準確（用戶 audit 通過率 ≥ 90%）
- handoff 到 myPPT 是否平順（用戶不需重新填 main thesis）

## 章節合併候選（請評估是否要合）

- §9 銜接設計 + §15 輸出檔案管理（共 90 行 → 可合 70 行）
- §11 反例 + §6 紅旗（共 140 行 → 可合 120 行）
- §13 fast-track 移到 SKILL.md（playbook 不重複）

## 大綱審查重點

1. 16 章節是否合理？是否需要刪減 / 合併？
2. §4 內容 4 層分類（F/O/I/U）vs myPPT 的 Tier S/A/B/C/D 是否會混淆？要不要統一？
3. §8 17 段母稿格式是否完整？是否要刪減為 12-13 段？
4. §7 教授問答預測 5-7 個是否合理？要更多 / 更少？
5. §10 範例要用真實 InterSubMod 案例還是虛構？
6. §14 與 weekly-report skill 的關係，您傾向 (a) 合併、(b) weekly-report 變 W1、(c) 並行不互動？
7. 預估 ~1000-1200 行是否合理？要更精簡（800 行）還是更詳盡（1500 行）？
