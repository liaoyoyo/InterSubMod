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

## 反模式（紅線）

- ❌ 從記憶生成 BibTeX → 40% 錯誤率
- ❌ 跳過 Google Scholar 直接用 WebSearch 結果 → 缺少存在性驗證
- ❌ 假設「這篇有名我記得」→ 即使 Transformer 也要查
- ❌ Claim 引用未實際對照原文 → reviewer 會抓

## 與其他 skill 的關係

獨立 skill；本專案 Phase 2 撰寫論文時使用。撰寫流程中，每加一個 `\cite{...}` 即觸發本 skill 一次驗證。
