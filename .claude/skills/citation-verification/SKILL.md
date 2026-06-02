---
name: citation-verification
description: 學術引用驗證：每個 citation 必須 WebSearch + Google Scholar 驗證後才寫入 .bib。USE WHEN：撰寫論文、檢查既有引用、產出 bibliography。 SKIP WHEN 內部 docs（無學術引用需求）、非 .bib 引用格式、純註解編輯、純 .py 程式碼開發、既有 .bib 已驗證過且無新增條目。
allowed-tools: Read, Write, WebSearch, WebFetch, Grep
user-invocable: true
---

# Citation Verification

**核心原則**：AI 從記憶生成的引用錯誤率約 40%。**每筆引用必須當下驗證**，不要寫完才檢查。

## 必走步驟

```
1. WebSearch 找論文（標題 + 第一作者 + 年份）
2. Google Scholar 查 site:scholar.google.com 確認存在
3. 對照 4 項：標題 / 作者（至少第一）/ 年份（±1 容許 preprint）/ 期刊或會議
4. 從 Google Scholar 「Cite」複製 BibTeX
5. 若引用具體 claim → WebFetch PDF 用關鍵字搜尋確認該 claim 真的出現
6. 通過 → 寫入 .bib；失敗 → 標 `[CITATION NEEDED]` 並告知用戶
```

## 驗證查詢範本

```
WebSearch: "<title>" "<first author>" <year>
WebSearch: site:scholar.google.com "<title>" <first author>
```

## 失敗處理（不要編造）

| 情境 | 動作 |
|------|------|
| Google Scholar 找不到 | 試 arXiv / DOI；皆無 → `[CITATION NEEDED]` 並回報用戶 |
| 標題或作者拼字不確定 | 用部分字串 + 年份再查；勿憑記憶補全 |
| Claim 在論文中找不到 | 移除該 claim 或改述為 weaker form；勿假設論文支持 |
| Citation count 異常低（<5） | 警示用戶可能是 obscure 或被撤稿 |

## 每個文獻 fact 附 source（2026-06-02 借鑑 PaperQA2 citation-by-construction）

延伸 CLAUDE.md §13 反捏造到**文獻 claim**：報告/論文裡每個引用某文獻的**具體陳述**（數字 / 結論 / 方法），必標**可追溯 source**（`(Author Year, §N/p.N)` 或 DOI + 該 claim 出現的段落）。原文 grep 不到支持 = 等同捏造（PaperQA2 即使 ~9% err 仍以「每個 fact 附 source」為設計核心）。與 `number_provenance`（內部產出數字）**互補分工**：B8 管外部文獻 claim、number_provenance 管內部數字。

## 反模式（紅線）

- ❌ 從記憶生成 BibTeX → 40% 錯誤率
- ❌ 跳過 Google Scholar 直接用 WebSearch 結果 → 缺少存在性驗證
- ❌ 假設「這篇有名我記得」→ 即使 Transformer 也要查
- ❌ Claim 引用未實際對照原文 → reviewer 會抓

## 嚴謹度繼承（/scientific-rigor）

citation-verification 是 /scientific-rigor 的「引用層」防線。本 skill 直接執行 §2 證據分級 + §4 DAG 的引用驗證部分：

- **§2 Evidence Tier 對應「引用分級」**:
  - L1 引用: KB primary source + verified URL + 本地 mirror
  - L2 引用: 一手論文 + DOI 可訪問
  - L3 引用: 二手 review + 業界 blog
  - L4 引用: 推測無源（標 `[需驗證]`）
  - L5 引用: 已知矛盾來源
- **§4 DAG 引用 confound 識別**:
  - 引用 paper claim 必對應其 DAG 假設（如 confounder 對應）
  - 若 paper 與本專案 DAG 不對齊 → flag「外部引用 vs 本地驗證 gap」
- **§8.4 Provenance 強制**: 每個引用必含 access_date + verified_at + L1-L5 tier
- **§7 Pre-reg**: 跨領域引用必對應事先註冊的「相關假說」欄

## 與其他 skill 的關係

獨立 skill；本專案 Phase 2 撰寫論文時使用。撰寫流程中，每加一個 `\cite{...}` 即觸發本 skill 一次驗證。
