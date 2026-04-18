# 20260304_連接knowledge.pptx

- Slides: 16
- Date: 2026-03-04

## Slide 1: 知識庫連接指南

知識庫連接指南

Knowledge Base Integration Guide

CCU Bioinformatics Lab · 2026-03

38
篇文件

9
大分類

4.7TB
資料量

8
分析工具

6
ONT 樣本

---

## Slide 2: 目錄 / Agenda

目錄 / Agenda

01

什麼是知識庫？
What is the Knowledge Base?

02

手動瀏覽方式
Manual Navigation

03

目錄結構與文件格式
Directory Structure

04

新成員建議閱讀順序
Onboarding Reading Path

05

AI 查詢 + MCP 設定
AI Query + MCP Setup

06

MCP 四個核心工具
Four MCP Core Tools

07

MCP 設定與驗證
MCP Configuration & Verify

08

近期更新與改善成果
Recent Updates & Results

---

## Slide 3: 什麼是知識庫？  What is the Knowledge Base?

什麼是知識庫？  What is the Knowledge Base?

一套結構化的 Markdown 文件系統，記錄 ONT 癌症長讀長測序研究所需的所有資訊。
A structured Markdown documentation system for ONT cancer long-read sequencing research.

規模 / Scale

38
篇文件

9
大分類目錄

4.7TB
資料量

8
分析工具

6
ONT 癌症樣本

雙重用途 / Dual Purpose

人工可翻閱
任何文字編輯器 / GitHub 皆可開啟
Human-readable via any text editor or GitHub

AI 可查詢
透過 MCP 協定直接查詢 38 篇文件
AI-queryable via MCP protocol — direct access

---

## Slide 4: 目錄結構 / Directory Structure

目錄結構 / Directory Structure

/big8_disk/liaoyoyo2001/Knowledge/

├── 01_data_overview/

資料版圖與統計  Data overview & inventory

├── 02_samples/

各癌症樣本詳細說明  Cancer sample details

├── 03_file_formats/

VCF / BAM 格式規格  File format specifications

├── 04_databases/

PON、SEQC2、參考基因組  Databases & references

├── 05_tools/

工具說明（LongPhase 家族等）  Tool documentation

├── 06_workflows/

完整分析流程  Analysis workflows

├── 07_script_docs/

腳本文件說明  Script documentation

├── 08_references/

路徑索引 / 論文索引  Path & paper index

├── 09_standards/

規範與決策記錄  Standards & decision log

---

## Slide 5: 每份文件的結構 / Document Structure

每份文件的結構 / Document Structure

每份文件開頭固定有三項，讓讀者快速判斷是否需要深入閱讀

1

一句結論  One-liner Conclusion
本文件是什麼、記錄了哪些內容

2

適用對象  Target Audience
誰應該讀這份文件

3

可直接執行命令  Ready-to-run Command
含驗證日期，可直接複製貼上執行

---

## Slide 6: 新成員建議閱讀順序 / Onboarding Reading Path

新成員建議閱讀順序 / Onboarding Reading Path

①  整個知識庫的地圖
README.md
The map of everything

②  6 個 ONT 樣本一覽
cancer_samples.md
6 ONT cancer samples at a glance

③  主分析流程
somatic_variant_calling.md
Main analysis workflow

④  查找特定工具參數
05_tools/（依需求）
Look up specific tool parameters

💡 新人 30 分鐘獨立上手，不占用資深成員時間 — Onboard in 30 minutes independently

---

## Slide 7: AI 查詢運作原理 / How It Works

AI 查詢運作原理 / How It Works

研究員
Researcher

Claude Code
AI Interface

MCP Server
knowledge_server.py

38 篇文件
Knowledge Base

← AI 回答 + 來源文件路徑（可追溯驗證）

即時查詢
不需等學長姐
Instant answers

中英文皆可
甲基化 = methylation
Bilingual search

附來源路徑
每個回答可追溯
Source-traceable

---

## Slide 8: MCP 四個核心工具 / Four MCP Core Tools

MCP 四個核心工具 / Four MCP Core Tools

搜尋

knowledge_search
關鍵字搜尋（中英文皆可）

knowledge_search("甲基化 somatic")

文件

knowledge_get_doc
取得文件全文

knowledge_get_doc("kb-05-tools-longphase-s")

路徑

knowledge_resolve_path
解析別名 → 正規路徑

knowledge_resolve_path("體細胞定相")

瀏覽

knowledge_list  ★NEW
依 tag / status 瀏覽文件

knowledge_list(tag="tool")

---

## Slide 9: 搜尋語法 / Search Syntax

搜尋語法 / Search Syntax

中文搜尋

甲基化

自動 bigram 展開：甲基化 → 甲基 + 基化

英文搜尋

phasing

含詞形變化 stem matching：phasing → phase

多詞搜尋

LongPhase somatic

全詞命中優先排序（×1.5 bonus）

排除詞

-germline

NOT operator：排除含此詞的文件

---

## Slide 10: MCP 設定方式 / Configuration

MCP 設定方式 / Configuration

方式 A：請 AI 幫你設定（最簡單）

在 Claude Code 中輸入：

幫我在這個專案設定 Knowledge MCP，
server 路徑是 /big8_disk/liaoyoyo2001/
Knowledge/scripts/mcp/knowledge_server.py

✅ AI 自動修改 .mcp.json，重載後即可使用

方式 B：手動設定（全域推薦）

codex mcp add knowledge -- python3
.../knowledge_server.py

✅ 所有新專案自動共用，無需重複設定

✅ 驗證安裝 / Verify Installation

codex mcp list
python3 .../knowledge_server.py --test
# 預期輸出：43/43 passed

43/43 ✓

---

## Slide 11: 專案層設定 .mcp.json / Project Config

專案層設定 .mcp.json / Project Config

放在專案根目錄，僅對該專案生效（可與全域設定並存）

{
"mcpServers": {
"knowledge": {
"type": "stdio",
"command": "python3",
"args": [
".../knowledge_server.py"
]
}
}
}

適用場合
專案需要獨立的 MCP 設定，或版本控制時需納入設定

與全域設定的關係
可同時存在；.mcp.json 設定只對當前專案有效

檔案位置
放在專案根目錄（與 .git 同層）

---

## Slide 12: 近期 Git 更新時序 / Development Timeline

近期 Git 更新時序 / Development Timeline

2026-02-10

6724b4a  初始建立：文件架構、CLAUDE.md、審計報告

2026-02-12

925c7d2  目錄重構：shell→scripts；補 README、索引

2026-02-27

e3df65f  規範層：新增 09_standards，決策記錄機制

2026-02-28

26c95e3  Phase 1：MCP Server 建立（FastMCP stdio，21 項測試）

2026-03-01

64698e1  Phase 2：MCP bug 修復（A1–A6）+ 12 份文件 description 改寫

2026-03-01

3337993  Phase 3：搜尋引擎升級（bigram+stemming+NOT）；測試 21→43

2026-03-01

683aef6  21 份文件補齊中英文 alias_paths，tags 語意化

---

## Slide 13: 改善成果 / Improvement Results

改善成果 / Improvement Results

指標 / Metric

改善前 / Before

改善後 / After

索引文件數

35 篇

38 篇（含 09_standards）

MCP 自動測試

21 項

43 項全通過 ✓

搜尋邏輯

OR（偽匹配）

OR + AND bonus（精準排序）

中文搜尋

整詞完整匹配

Bigram 展開（甲基化→甲基+基化）

英文搜尋

完整詞匹配

Stemming（phasing→phase）

Alias paths

多數為空

21 份文件補齊中英文別名

Tags 品質

通用（workflow）

語意化（phasing, haplotagging…）

---

## Slide 14: 對工作流程的影響 / Workflow Impact

對工作流程的影響 / Workflow Impact

⛔ 設定前 / Before

✅ 設定後 / After

「這個 BAM 路徑？」→ 問學長姐

「LongPhase-S 的參數？」→ 翻舊腳本

新人需要長時間 shadow 才能上手

回答不附來源，難以追溯驗證

Claude Code 直接問 → 秒得有來源的答案

中英文查詢皆可（甲基化 = methylation）

新人 30 分鐘獨立上手，不占資深成員時間

AI 回答附來源文件路徑，可追溯驗證

---

## Slide 15: 總結 / Summary

總結 / Summary

📚
38 篇結構化文件
38 docs in 9 categories

🤖
MCP 協定 AI 直接查詢
AI-queryable via MCP protocol

🔍
智慧搜尋引擎
Bigram + Stemming + NOT

⚡
30 分鐘快速上手
New members onboard in 30 min

Q & A  ·  問題討論

CCU Bioinformatics Lab · Knowledge Base v2 · 2026-03

---

## Slide 16: 附錄：快速參考卡片 / Quick Reference

附錄：快速參考卡片 / Quick Reference

📂 路徑：  /big8_disk/liaoyoyo2001/Knowledge/

中文搜尋

knowledge_search("甲基化")

英文+排除

knowledge_search("phasing -germline")

多詞搜尋

knowledge_search("LongPhase somatic")

取得文件

knowledge_get_doc("kb-05-tools-longphase-s")

瀏覽工具

knowledge_list(tag="tool")

解析別名

knowledge_resolve_path("體細胞定相")

---
