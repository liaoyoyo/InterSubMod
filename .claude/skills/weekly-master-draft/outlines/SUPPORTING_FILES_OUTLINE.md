---
title: weekly-master-draft 其他支援檔案大綱（待用戶審查）
date: 2026-05-01
status: outline-for-review
---

# 其他支援檔案大綱（templates + prompts + examples + README）

## A. README.md（≈ 70 行）

skill 入口；4 區塊：
1. 一段 200 字概覽（解決什麼、為何先做母稿才做 PPT、預期 30-60 min 流程）
2. 何時觸發（5-7 條 keyword + 場景對應）
3. 7 階段 W1-W7 使用流程速覽
4. 與 myPPT / weekly-report / structured-tech-report 的關聯與分工

## B. templates/ 4 檔（每檔 ≈ 70-90 行）

母稿主線類型對應的 narrative skeleton。每檔結構統一：
1. 用途定義（1 段）
2. 觸發 keyword（4-5 條）
3. 17 段填寫指引（針對該主線類型，哪幾段該詳寫、哪幾段該略寫）
4. 範例：1 個虛構週報範例母稿（300-500 字）
5. 與其他主線類型的差別

### templates/progress_focus.md（進展型）
- 重點段：3, 4, 7, 8, 16
- 略寫段：5（待確認少）、9（PPT 模板會 default executive_summary）
- 範例：「本週驗證 X 假設，3 樣本均為 +0.05 ΔF1」

### templates/problem_focus.md（問題型）
- 重點段：5, 11, 17
- 略寫段：3（事實少）、16（下週計畫待教授判斷）
- 範例：「本週發現 self-phasing 17.3:1 異常，根因待確認」

### templates/advisor_consult.md（求協助型）
- 重點段：5, 17（教授問答最重要）
- 必補段：每個方向的利弊分析
- 範例：「本週 3 個 candidate direction，請教授判斷」

### templates/new_direction_explore.md（探索型）
- 重點段：4（pilot 結果）、5（待確認多）、11（補資料）
- 略寫段：3（事實少）
- 範例：「本週 pilot HPFineNGroups，AUC 0.62 但需 7 樣本驗證」

## C. prompts/ 7 檔（每檔 ≈ 50-70 行）

每檔結構：
1. 使用時機（1 段）
2. 觸發前置條件（哪個 W 階段 / C checkpoint）
3. AskUserQuestion 參數模板（≤ 5 個 question / options）
4. 預期輸出格式
5. 用戶 cancel / 修改 / 再迭代的處置

### prompts/raw_data_collect.md（W1 → C0）
- 時機：skill 啟動後第一步
- 自動掃描來源：`grep keyword` 模板（git log / evidence_ledger / experiments 目錄）
- 4 options：所掃結果完整 / 補漏 / 移除誤抓 / 全部重來

### prompts/main_thread_identify.md（W2 → C1）
- 時機：raw data 確認後
- 4 options（progress / problem / consult / explore）
- 強制要求：main statement ≤ 30 字

### prompts/content_classify_4tier.md（W3 → C2）
- 時機：每筆素材分類前
- AI 預設標 [F]/[O]/[I]/[U]，列表呈現給用戶
- 5 options：標籤正確 / 改成 X / 補 source / 拆分多項 / 棄用

### prompts/priority_sort_4bucket.md（W4 → C3）
- 時機：每筆素材分桶前
- AI 計算 5 維度評分 → 預設 PPT/講稿/備註/暫存
- 4 options：分桶正確 / 升桶 / 降桶 / 暫存

### prompts/logic_check_redflag.md（W5 → C4）
- 時機：分桶完成後
- AI 掃描全文紅旗（過度宣稱 / 流水帳 / 教授視角缺）
- 列出觸發句 + 改寫建議
- 4 options：採納改寫 / 自己改 / 此 flag 為誤判 / 跳過

### prompts/professor_qa_predict.md（W6 → C5）
- 時機：紅旗檢查完成後
- AI 列 5-7 個追問 + 預備回答
- 4 options：完整 / 補漏 / 改寫某題 / 加教授個人偏好

### prompts/master_draft_finalize.md（W7 → C6）
- 時機：所有 checkpoint 完成
- AI 組裝 17 段母稿 + 用戶逐段批准
- 4 options：批准全文 / 修改 §N / 重跑某 W / 棄用整份

## D. examples/ 1 檔（≈ 200 行）

### examples/master_draft_example.md
- 完整虛構案例：2026-05-XX 週報
- 場景：HPFineNGroups TO 模式 phasing 信號（探索型 + 問題型混合）
- W1-W7 全程展示
- 17 段母稿全文
- handoff 到 myPPT 後預期 outline

## 大綱審查重點

1. README.md 是否要加「Quick start example」（30 sec 啟動指引）？
2. 4 個 templates 是否覆蓋您所有週報場景？是否需新增 / 合併？
3. 7 個 prompts 是否過多？(C0-C6 七個 checkpoint 對應 7 個 prompt)能否合併（如 C2+C3 合一個 prompt）？
4. examples/ 用 1 個還是 2 個範例（一個進展型 + 一個求協助型）？
5. 範例週報用真實 InterSubMod 案例還是純虛構？
