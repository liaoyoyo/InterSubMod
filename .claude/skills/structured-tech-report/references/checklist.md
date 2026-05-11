# 完成前自我檢核 25 項（structured-tech-report）

每項標 ✅ / ⚠（需補但不擋） / ❌（必擋）。

**門檻**：≥ 22 ✅ 才能進 Step 5（登錄）；< 18 ✅ 退回 Step 3（重填骨架）。

---

## A. 結構完整性（5 項）

- [ ] **A1** 13 段全到（§0 TL;DR ~ §13 結論），順序未亂
- [ ] **A2** §7 拆 §7.1（非工程版）+ §7.2（工程版）兩小節
- [ ] **A3** 開頭 HTML comment metadata 含 `build_date / agent / status / inputs / outputs / verdict`
- [ ] **A4** 結尾「附錄 A 引用清單」與「附錄 B 變更歷史」就位
- [ ] **A5** 章節標題使用 `## 1. 報告目的`、`## 2. 系統背景` … 格式（編號連續）

## B. 內容品質（8 項）

- [ ] **B1** §0 TL;DR 限 3-5 行，含問題／修改／結果／狀態四要素
- [ ] **B2** §4 問題描述**有量化指標**（具體數字／百分比／指標名）
- [ ] **B3** §5 根因採用 5 Whys 或 Ishikawa 或 ADR Decision Drivers，**未停在「修了就好」**
- [ ] **B4** §6 至少列 2 個候選方案（不只採用的），每個含 Pros/Cons 與不採用理由
- [ ] **B5** §7.1 非工程版**未出現** file path、function name、commit hash、API 名
- [ ] **B6** §7.2 工程版含 file:line（涉 .cpp 改動 → commit hash 必填或標 `⚠ 待補 commit hash`）
- [ ] **B7** §8 新舊比較**有表格**，欄位含「原本／新／改變原因」
- [ ] **B8** §13 結論限 3-5 句，回答「問題／為什麼／改了什麼／結果／後續」
- [ ] **B9 ⭐**（fact-check）§4／§6／§7.2／§8 所有具體數字（百分比、metric 值、行數、commit hash、檔案路徑）必須有可驗證的 source 標註於該段落或附錄 A。**禁止**：(a) 引述 explorer / 二手摘要的數字而未開讀原始 .md／.cpp；(b) 用「合理估計」填空而不標 `⚠ 待確認`；(c) 數字無 source 寫死在表格內。**正確做法**：每個具體數字必須能對應到「附錄 A 引用清單某檔某段／某 file:line」。觸發此項違反 = 必擋。
- [ ] **B10**（歸屬正確性）涉及 C++ 程式碼改動的段落，必須 grep `find <repo> -name "<file>"` 確認檔案實際位置；不可假設「`HaplotagProcess` 在 InterSubMod」「`getVote` 在 src/core/」這類字面推論。涉外部 fork（如 longphase-to-mod）時，**必須在 §2 系統背景與 §7.2 工程版明確標明 repo 歸屬**。

## C. 驗證可觀察性（4 項）

- [ ] **C1** §9 採 `Step → Verify` 格式，每步驟有**外部可觀察**驗證（檔案／命令／數值範圍）
- [ ] **C2** §9 **無**弱驗證標準（「讓它能動」「看起來合理」「double-check」）
- [ ] **C3** 涉統計顯著性 → 已呼叫 `/auc-confound-guard` 三關，並引用結果
- [ ] **C4** 涉 7 樣本驗證 → 已呼叫 `/multi-sample-consistency`，並引用 canonical 排序表

## D. 引用與連結（5 項）

- [ ] **D1** 所有 .md 路徑前綴 `InterSubMod/...`（對齊專案 hook）
- [ ] **D2** 所有引用的 .md／.cpp／.tsv 路徑經 `ls` / `Read` 驗證存在
- [ ] **D3** 涉 MEMORY 引用時，引用 ID 匹配 `MEMORY.md` 索引（非孤兒 link）
- [ ] **D4** 圖表落點 `figures/{section}/figXX_*.png` 命名規則一致；相對路徑 ≤ 2 層
- [ ] **D5** 至少 1 條反向連結（INDEX 或上層報告 link 回本報告）

## E. 落點與登錄（3 項）

- [ ] **E1** 檔案落點符合「報告類別 → 目錄」對應（validated／in_progress／decisions）
- [ ] **E2** 對應 INDEX.md 已新增 1 行（含日期／路徑／一句結論）
- [ ] **E3** 若涉 14 結論之一 → `research_landscape/00_INDEX.md` 同步更新

---

## 備註欄

| Item | ✅/⚠/❌ | 補充說明 |
|------|--------|---------|
| A1 | | |
| A2 | | |
| A3 | | |
| ... | | |
| E3 | | |

**總分**：__ / 25 ✅
