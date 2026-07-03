<!--
建立時間: 2026-06-09
狀態: hub (G6 論文分析聚焦確認資料夾 — 導航 index)
報告類型: paper_focus_hub_index
受眾: 廖子游（聚焦確認）· PI · 其他 AI session
provenance_note: 本資料夾不產新數字；所有 metric 沿用 session 已驗證集合（workflow wf_a8ccbb34-3f7 原檔 grep-back，2026-06-09）+ knowledge/11_external_literature/10。每數字標 🟢P/🟡S/🔵framing。
-->
<!-- provenance-verified: 三色分級來自 wf_a8ccbb34-3f7（8 組 auditor 原檔單檔 grep）+ wf_9e169112-573 收斂稽核；本 hub 為導航，非新分析。內部 L1-L2 / 外部 L3。 -->

# G6 論文分析聚焦確認 — 資料夾索引（hub）

> **這個資料夾是什麼**：把「論文要做什麼、哪些可說明/可解釋有證據、哪些方向要確認、用什麼驗證機制」收斂成一個**可持續編輯的工作中樞**。和 `docs/concepts/`（給其他 AI 的單頁地圖）、`docs/plans/`（執行 DAG）互補 —— 這裡是**聚焦確認 + 細節任務樹**的家。
> **怎麼用**：先讀 `00`（本檔）→ `01 背景`定錨題目 → `02 文件庫`查權威來源 → `03 方向卡`看每個方向的 verdict + 燒烤 → `04 對應表`查每個結論的驗證機制 → `05 任務樹`領要做的事。

---

## ⭐ 如何讀（認知負擔最小化 — 4 層揭露機制）

> 你的回饋：報告缺「重點邏輯框架 + 重要性分層揭露」。本資料夾**每份文件一律照下面 4 層寫，讀到夠就停**：

| 層 | 是什麼 | 讀它的時機 |
|----|--------|-----------|
| **L0 一眼結論** | 開頭 1 句粗體（blockquote）| 只想知道結果 → 讀完就走 |
| **L1 重點邏輯** | 3–5 條 bullet（為什麼／推理骨架）| 30 秒抓推理鏈 |
| **L2 細節展開** | 表格／方向卡 | 要做／要確認時才讀 |
| **L3 原始溯源** | 🟢🟡🔵 + 檔:行 | 引用／驗證時才查 |

**行內重要性標記**（掃描用）：🔴 決定性（卡論文／要你決策）｜⭐ 重要有價值｜◽ 背景／已穩定｜🟢 有證據 / 🟡 待對賬 / 🔵 framing。

> **這頁本身就是 L0+L1**：只讀 `00_INDEX` 就能掌握 80%（現況一句 + 燒烤 3 點 + 可行性一覽 + 故事線）；要細節再進 `01`–`06`。

---

## 檔案地圖（子資料夾分類）

| 路徑 | 分類 | 內容 |
|------|------|------|
| `00_INDEX.md` | hub | 本檔：導航 + 4 層揭露機制 + 可行性一覽 + 故事線 |
| `01_focus_notes/` | **.md 敘述（聚焦確認）** | `01` 背景/關鍵字/共驗證 · `02` 文件庫/🟡待對賬（✅2026-06-09 T-PROV 結案）· `03` 方向卡 · `04` 知識↔驗證機制↔研究機制 · `05` 任務樹 · `06` 目標/故事線 · `07` 研究觀察任務清單 · **`08` 任務執行回報（catalog + provenance）⭐2026-06-09** |
| `02_paper_framework/` | **.md 敘述（論文框架）** | 論文 IMRaD 逐節寫作重點與確認 · catalog schema 提案 · **位點甲基分群catalog_結果_R6（332,705 loci × 7 TAG，骨架已建✅）⭐2026-06-09** · **論文架構_正式學術版（Slide2Thesis Nature-journal 格式，Pandoc-ready）** |
| `03_references/` | **論文參考外部資料** | 外部文獻索引（→ `knowledge/11` ~38 篇）+ 引用清單 + 3 不對齊 + PDF 放置處 |
| `04_figures/` | **圖片** | 論文需要的圖 manifest（哪張圖／料／tier／狀態）+ 生成的 PNG 放這 |
| `05_html_staging/` | **可檢視確認的暫存 HTML** | 互動確認 HTML（勾選／回答），如 `論文框架確認.standalone.html` |

---

## 三色圖例（貫穿全資料夾 — 回答「可說明 vs 可解釋有證據」）

| 色 | 意義 | 可信度 |
|----|------|--------|
| 🟢 **可解釋有證據（P-level）** | 原始 json/jsonl/md 單檔 grep 對賬 exact | 投稿可直接引用 |
| 🟡 **原檔待對賬（S-level）** | 結論寫在 synthesis 文件，但指定源檔本輪 grep 未命中（檔不存在/值在另一檔/derived） | 投稿前必定位正確源檔再 grep |
| 🔵 **可說明（framing 立場）** | 論文敘述/定位/角色分配；由整體證據結構支撐，本身非單一硬數字 | 是 framing 不是證據 |

---

## 現況 BLUF（一句）

一篇能過 review 的論文今天就存在 = read-level LOH/haplotype + 甲基的 **characterization + tooling 論文（不是 variant filter）**；主體 = **三~四道防彈 NEGATIVE** + phasing 脊柱（Grade B+ 非 A）+ ASM copy-confounded 支撐。**HD-1 已收斂為非承重（骨幹 germline-anchored 非循環 by design + 循環 phasing-spine 降 §2.6 支撐 observation；R-SELFREF 降 optional，2026-07-03）**。

**🔴 本資料夾最需要你確認/被燒烤的三件事**（詳見 `03`）：
1. **「共驗證」要定義成「正交佐證」不是「互相篩選」** —— 後者已死四道。
2. **AF↔NGroups 不是甲基 subclone marker** —— HD-4（06-08）已判定是 phasing（NGroups=HP-tag count）。你方向卡 D4 的原假設這週已被你自己的分析推翻。
3. **「ASM 3.95%>1.07% 能找到明顯 TP」要分清 enrichment（characterization）≠ retrieval（filter）** —— 4% 是低 sensitivity，FP 也有 1.07%，COLO829 TP≈FP。

---

## ⭐ 可行性一覽（L1 — 讀這張就知道 11 方向全貌）

> 裁決：✅ CONFIRM 可行/已成立｜⭐ GO 高價值可做｜🟡 PROBE 要先確認/有條件｜🔴 NO-GO 否決（附原因）。細節+證據+外部佐證 → `03 方向卡`。

| # | 方向（你提的）| 裁決 | 一句原因（證據級）|
|---|------|------|---------|
| 你-1 | 甲基輔助 phasing／救 unphase (umtag) | 🟡 PROBE | 白地真實+外部已驗，缺 yardstick／單樣本／未註冊 → future-work（🟡）|
| 你-2 | 甲基對 read 分群 | 🟡 PROBE | ISM 已做；但揭示 germline allelic 非 somatic subclone（🟢）|
| 你-3 | cis-candidate 大規模驗證 | 🟡 PROBE | 真開放方法題；chr17 唯一乾淨、6 個 untestable（🟡）|
| 你-4 | AF↔NGroups 當 subclone 甲基 marker | 🔴 **NO-GO** | **HD-4 已推翻＝phasing 非甲基**（NGroups=HP-tag count，🟢）|
| 你-5 | strong-ASM FP 富集 特殊刪除 | 🔴 NO-GO(filter) | regression-to-extreme＋反判別 OR=0.194＋與覆蓋/LOH filter 冗餘（🟡）|
| 你-6 | 甲基做 subclone 分析 | 🔴 NO-GO(原法) | B2 clustering NEG：ARI 0.135<null 0.177＝germline-allelic（🟢）|
| ISM-1 | ISM 互補性（CramersV vs Δβ）| ✅ CONFIRM | 互補成立（framing 🔵 + 例 🟡）|
| ISM-2 | ISM 過嚴/過鬆 | 🟡 要確認 | CramersV gate 偏嚴 → 487 latent 真分群被 gate 掉（要 audit）|
| ISM-3 | 合併 Δβ 找位點 | ⭐ PROBE/code | 可當**互補發現通道**（非取代 CramersV）|
| ISM-4 | 位點甲基分群標籤 catalog | ⭐ **GO** | 高價值 characterization 交付物（= 你的核心需求）|
| 困難-1 | paired TP/FP 篩選 | 🔴 NO-GO(filter) | FP 太少 power 不足＋filter 死＋normal 甲基加重 confound（🟢）|

## 故事線一句（L0 — 如何敘述）

> **在 LOH 約束下，somatic／haplotype／甲基三訊號 read-level 交織、互相「正交佐證」同一個 haplotype 重塑事件**（不是互相篩選）：somatic SNV 在保留的單一 haplotype 上形成子家族（phasing 佐證 somatic↔haplotype），其中少數位點伴隨真實 ASM（甲基 characterize 該事件）；交付 = read-level characterization + 位點甲基分群 catalog + 一份誠實的「甲基能佐證、不能篩選」機制目錄。完整 → `06`。

---

## 權威真值順序（衝突時的裁決）

1. `InterSubMod/knowledge/11_external_literature/10_paper_readiness_convergence.md`（論文就緒收斂，最完整）
2. `InterSubMod/research/autoresearch/evidence_ledger.jsonl`（append-only 歷史真理；**衝突時以 ledger 為準**）
3. `InterSubMod/state/active.json`（active cycle 快照）
4. 視覺版：`InterSubMod/docs/concepts/2026/06/20260609_G6_研究構想三層架構_漏斗研究卡證據鏈_01.standalone.html`

⚠ tsg 專案（`research/tsg_promoter_asm_reviewer/genome_survey_v2/`）06-07 仍在寫 → T2/T3/T5 數字可能再動，投稿前以該專案定稿為準。

---

## 配套文件（不在本資料夾，但相關）

- 三層架構視覺版（漏斗/研究卡/證據鏈/roadmap）：`InterSubMod/docs/concepts/2026/06/20260609_G6_研究構想三層架構_漏斗研究卡證據鏈_01.standalone.html`
- 執行 DAG：`InterSubMod/docs/plans/20260608_G6_paper_執行計劃_DAG_重點任務與依賴_01.standalone.html`
- 給其他 AI 的現況地圖：`InterSubMod/docs/concepts/2026/06/20260608_研究現況地圖_整體目標與流程_給其他AI_01.md`
- NEGATIVE 底座草稿：`InterSubMod/docs/reports/in_progress/2026/06/20260608_G6_methods_negative_backbone_draft_01.md`
