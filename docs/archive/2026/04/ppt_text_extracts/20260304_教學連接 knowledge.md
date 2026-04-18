# 20260304_教學連接 knowledge.pptx

- Slides: 14
- Date: 2026-03-04

## Slide 1: 知識庫連接指南

知識庫連接指南

Knowledge Base Integration Guide

►  如何手動查閱知識庫
How to browse the knowledge base manually
►  如何透過 AI 自動查詢
How to query with AI assistance
►  如何設定 MCP 連線
How to configure the MCP connection

CCU Bioinformatics Lab  |  2026-03

> **Notes:** 【投影片目的】
封面。說明本簡報三大主題。

【備忘要點】
1. 本簡報適合對象：新進實驗室成員、或現有成員想導入 AI 查詢功能者。
2. 知識庫（Knowledge Base）= /big8_disk/liaoyoyo2001/knowledge/
3. MCP（Model Context Protocol）是 Anthropic 定義的標準協定，讓 Claude Code 可連接外部資料源。
4. 本簡報時間估計：15–20 分鐘。

【開場建議】
「今天分三個部分：第一，知識庫是什麼；第二，怎麼用；第三，用了之後有什麼差別。」

---

## Slide 2: 本場重點

本場重點

Agenda

2

1

一、什麼是知識庫？

What is the Knowledge Base?

2

二、如何手動瀏覽

How to Browse Manually

3

三、AI 查詢 + MCP 設定

AI Query + MCP Configuration

4

四、工作成效與近期更新

Impact & Recent Updates

CCU Bioinformatics Lab  |  2026-03

> **Notes:** 【投影片目的】
概覽四大章節，讓聽眾建立整體框架。

【各章節說明】
01 什麼是知識庫 — 規模數字（38 篇文件、4.7TB、8 個工具）+ 雙重用途。
02 如何手動瀏覽 — 目錄結構、新成員建議閱讀順序。
03 AI 查詢 + MCP — 四個工具、搜尋語法、兩種設定方法。
04 成效 — Before/After 對比、Git 時序、搜尋引擎升級數字。

【時間分配建議】
01：3 分鐘 / 02：3 分鐘 / 03：6 分鐘 / 04：3 分鐘 / 問答：5 分鐘

---

## Slide 3: 一、什麼是知識庫？

一、什麼是知識庫？

What is the Knowledge Base?

3

38

知識文件

Documents

9

分類目錄

Categories

4.7 TB

資料量

Data Covered

8

分析工具

Tools

►

結構化 Markdown 文件系統

Structured Markdown documentation system

►

人工可翻閱（文字編輯器 / GitHub）

Human-readable via text editor or GitHub

►

AI 可透過 MCP 協定直接查詢

AI-queryable via MCP — no manual search needed

►

涵蓋 6 個 ONT 癌症樣本與完整流程

6 ONT cancer samples + full analysis workflows

CCU Bioinformatics Lab  |  2026-03

> **Notes:** 【投影片目的】
用數字讓聽眾快速掌握知識庫規模，並說明雙重用途。

【數字說明】
- 38 篇文件：涵蓋 01–09 九大目錄
- 9 個分類：data_overview / samples / file_formats / databases / tools / workflows / script_docs / references / standards
- 4.7 TB：6 個 ONT 癌症樣本（HCC1395、COLO829、H1437、H2009、HCC1937、HCC1954）
- 8 個工具：LongPhase、LongPhase-S、LongPhase-TO、ClairS、DeepSomatic、DeepVariant、WhatsHap、InterSubMod

【雙重用途】
- 人工：任何文字編輯器或 GitHub 均可閱讀
- AI：透過 MCP Server（knowledge_server.py）讓 Claude Code 直接查詢 38 篇文件
  不需要人工描述背景，AI 自動找到最相關的文件並附來源路徑

---

## Slide 4: 二、知識庫目錄結構

二、知識庫目錄結構

Knowledge Base Directory Structure

4

01

資料總覽與統計

Data overview & inventory

02

各癌症樣本詳細說明

Cancer sample details

03

VCF / BAM 格式規格

File format specifications

04

PON、SEQC2、參考基因組

Databases & reference genomes

05

工具說明（LongPhase 家族）

Tool documentation

06

完整分析流程

Analysis workflows

07

腳本文件說明

Script documentation

08

路徑索引 / 論文索引

Path & paper index

09

規範與決策記錄

Standards & decision log

■ 每份文件開頭：一句結論 + 適用對象 + 可直接執行命令（含驗證日期）  |  Every doc: one-liner conclusion + target audience + ready-to-run command

CCU Bioinformatics Lab  |  2026-03

> **Notes:** 【投影片目的】
讓聽眾了解知識庫的九大目錄架構。

【各目錄說明】
01 data_overview — 資料版圖、HG002 驗證集、儲存統計
02 samples — HCC1395、COLO829、H1437、H2009、HCC1937、HCC1954 六個樣本各有獨立說明
03 file_formats — VCF 格式（somatic / germline 差異）、BAM / CRAM、Modcall 格式
04 databases — Panel of Normals (PON)、SEQC2 truth set、GRCh38 參考基因組路徑
05 tools — LongPhase（germline）、LongPhase-S（tumor-normal）、LongPhase-TO（tumor-only）、ClairS、DeepSomatic
06 workflows — somatic calling 完整流程、phasing 流程、甲基化分析、benchmark 方式
07 script_docs — auto_run、run_benchmark 等腳本的說明文件（非腳本本身）
08 references — 所有伺服器路徑索引、論文索引、腳本索引
09 standards — Anthropic 規範、命名標準、MCP onboarding 指南、決策記錄

【文件結構規則】
每份文件開頭固定三欄：①一句結論 ②適用對象 ③可直接執行命令（含驗證日期）

---

## Slide 5: 新成員入門路徑（建議閱讀順序）

新成員入門路徑（建議閱讀順序）

Recommended Onboarding Reading Path

5

1

README.md

整個知識庫的地圖，30 分鐘掌握全局

The map of everything — overview in 30 minutes

2

02_samples/cancer_samples.md

六個 ONT 癌症樣本一覽（路徑、深度、工具版本）

Overview of 6 ONT cancer samples: paths, depth, tool versions

3

06_workflows/somatic_variant_calling.md

體細胞變異分析完整流程（ClairS → LongPhase-S → benchmark）

Full somatic calling pipeline: ClairS → LongPhase-S → benchmark

4

05_tools/（依需求查閱）

依任務查找特定工具的參數、陷阱與執行範例

Look up tool parameters, pitfalls, and examples as needed

CCU Bioinformatics Lab  |  2026-03

> **Notes:** 【投影片目的】
提供新成員第一週的具體閱讀行動計畫。

【各步驟補充】
Step 1 README.md
  - 包含：資料版圖、目錄速查、CLAUDE.md 連結
  - 閱讀後應能回答：「資料放哪裡？」「工具在哪？」
Step 2 cancer_samples.md
  - 6 個樣本各有：BAM 路徑、深度、ClairS 版本、可用工具版本
  - 重要：ClairS v0.4.0 只有 HCC1395 和 COLO829；v0.4.1 才有全部 7 個樣本
Step 3 somatic_variant_calling.md
  - 完整流程：原始 BAM → ClairS → phasing（LongPhase-S）→ benchmark
  - 包含常見陷阱：--ont 參數、bgzip vs gzip 差異
Step 4 05_tools/
  - 查完整體流程後，再依任務需求查個別工具
  - 每個工具文件都有：一句結論 + 可直接複製的執行命令

【常見問題】
Q: 需要全部讀完嗎？
A: Step 1-3 約 1-2 小時，之後按需查閱即可。

---

## Slide 6: 三、AI 如何查詢知識庫？

三、AI 如何查詢知識庫？

How AI Queries the Knowledge Base

6

研究員

Researcher

在 Claude Code 輸入問題

Types question
in Claude Code

MCP Server

knowledge_server.py

自動搜尋 38 篇知識文件

Auto-searches
38 documents

知識文件庫

Knowledge Base

Markdown 文件 01–09 目錄

Structured Markdown
directories 01-09

→

→

←  AI 回答附來源文件路徑（可追溯驗證）  |  AI answers include source document path

■ 實際案例  Real Example

問：「LongPhase-S somatic_haplotag 支援 --ont 嗎？」

→ AI：不支援（陷阱！）—— 來源：05_tools/longphase_s.md

Q: "Does LongPhase-S somatic_haplotag support --ont?"  →  AI: No — Source: 05_tools/longphase_s.md

CCU Bioinformatics Lab  |  2026-03

> **Notes:** 【投影片目的】
說明 AI 查詢的運作原理，讓聽眾了解 MCP 的角色。

【MCP 工作原理】
1. 研究員在 Claude Code 輸入問題（中英文均可）
2. Claude Code 呼叫 knowledge_server.py 提供的 MCP 工具
3. knowledge_server.py 對 38 篇文件執行關鍵字搜尋（OR+AND bonus 演算法）
4. 回傳最相關文件的摘要 + 路徑給 Claude Code
5. Claude Code 整合知識庫內容後回答，並附上來源文件路徑

【為什麼需要附來源路徑？】
- 防止 AI 幻覺（hallucination）
- 讓使用者可以直接查閱原始文件驗證
- 建立「AI 答案有依據」的信任

【MCP Server 技術細節】（供技術聽眾備用）
- 實作：FastMCP + stdio 傳輸
- 路徑：/big8_disk/liaoyoyo2001/knowledge/scripts/mcp/knowledge_server.py
- 測試：python3 knowledge_server.py --test → 預期 43/43 通過

---

## Slide 7: MCP 四個查詢工具

MCP 四個查詢工具

Four MCP Query Tools

7

knowledge_search

關鍵字搜尋（中英文皆可）

Keyword search — Chinese or English

knowledge_search("LongPhase somatic")

knowledge_get_doc

取得文件完整內容

Retrieve full document content

knowledge_get_doc("kb-05-tools-longphase-s")

knowledge_resolve_path

別名 → 正規路徑解析

Resolve alias to canonical path

knowledge_resolve_path("體細胞定相")

knowledge_list  ★

依 tag / status 瀏覽

Browse docs by tag or status

knowledge_list(tag="tool")

CCU Bioinformatics Lab  |  2026-03

> **Notes:** 【投影片目的】
介紹四個 MCP 工具的用途與範例。

【各工具細節】
1. knowledge_search(query, tags?, status?, max_results?)
   - 支援 OR+AND bonus 搜尋演算法
   - 中文自動 Bigram 展開：甲基化 → {甲基, 基化}
   - 英文自動詞形還原：phasing → {phase, phas}
   - NOT 運算子：-germline 排除含此詞的文件
   - tags 參數：可過濾特定類別，e.g. tags=["tool"]
   - max_results 預設 5

2. knowledge_get_doc(doc_id_or_path, section?)
   - doc_id 格式：kb-05-tools-longphase-s（見 knowledge://catalog）
   - section 可指定章節，避免回傳過長

3. knowledge_resolve_path(path_or_alias)
   - 支援中英文別名：「體細胞定相」→ 05_tools/longphase_s.md
   - 21 份文件已補齊 alias_paths（中英文）

4. knowledge_list(tag?, status?)
   - tag 可選：tool / workflow / reference / standard / sample / format / database
   - status 可選：active / legacy / wip

---

## Slide 8: 搜尋語法與進階功能

搜尋語法與進階功能

Search Syntax & Advanced Features

8

甲基化

中文搜尋（Bigram 展開）

甲基化 → {甲基, 基化}  自動展開多字詞

Auto bigram: 甲基化 → {甲基, 基化}

phasing

英文搜尋（詞形還原）

phasing → {phase, phas}  無需安裝外部套件

Suffix stemming: phasing → {phase, phas}

LongPhase somatic

多關鍵字：全詞命中優先

全部詞命中 → 分數 ×1.5，精準排序消除偽匹配

OR + AND bonus: all-match ranked ×1.5 over partial

-germline

NOT 運算子：排除特定文件

在詞前加 - 可排除含該詞的所有文件

Prefix - to exclude documents containing the term

CCU Bioinformatics Lab  |  2026-03

> **Notes:** 【投影片目的】
教導聽眾如何有效使用知識庫的搜尋語法。

【搜尋引擎技術細節（Phase 3 升級）】
- 演算法：OR with AND bonus
  每個關鍵字獨立計分（確保召回率 recall）
  若所有關鍵字都命中同一文件，分數乘 1.5（提升精確率 precision）

- 中文 Bigram 展開：
  「甲基化分析」會被展開為 {甲基, 基化, 化分, 分析}
  使得部分詞匹配也能找到相關文件

- 英文 Stemming（規則式，無需 NLTK/spaCy）：
  剝除常見詞尾 -ing, -ed, -er, -tion 等
  phasing → phase → phas（多層剝除）

- NOT 運算子：
  -germline LongPhase  → 找 LongPhase 相關但排除 germline 的文件

【實際使用建議】
- 找工具文件：knowledge_search("工具名稱") + knowledge_list(tag="tool")
- 找流程文件：knowledge_search("somatic calling") + knowledge_list(tag="workflow")
- 找特定路徑：knowledge_resolve_path("樣本名稱" 或 "工具名稱")

---

## Slide 9: 設定方法 A：讓 AI 幫你設定

設定方法 A：讓 AI 幫你設定

Setup Method A — Let AI Configure It

9

最簡單的方式：直接在 Claude Code 中說以下指令

Easiest method: just tell Claude Code to set it up for you

幫我在這個專案設定 Knowledge MCP，
server 路徑是
/big8_disk/liaoyoyo2001/knowledge/scripts/mcp/knowledge_server.py

✓ AI 自動修改 .mcp.json

AI auto-edits the .mcp.json file for this project

✓ 重新載入後即可使用 MCP 工具

MCP tools available after reloading Claude Code

CCU Bioinformatics Lab  |  2026-03

> **Notes:** 【投影片目的】
示範最簡單的 MCP 設定方式。

【完整操作步驟】
1. 開啟 Claude Code 並進入目標專案目錄
2. 在對話框輸入（複製貼上上方指令）
3. Claude 會：
   a. 檢查專案根目錄是否有 .mcp.json
   b. 若無：自動建立 .mcp.json 並寫入設定
   c. 若有：讀取現有設定並新增 knowledge server
4. 重新載入 Claude Code（Ctrl+R 或重新開啟）
5. 測試：輸入 knowledge_search("test") 確認可用

【.mcp.json 內容（AI 會自動寫入）】
{
  "mcpServers": {
    "knowledge": {
      "type": "stdio",
      "command": "python3",
      "args": ["/big8_disk/liaoyoyo2001/knowledge/scripts/mcp/knowledge_server.py"]
    }
  }
}

【注意事項】
- knowledge_server.py 需要 Python 3.8+ 且 FastMCP 套件已安裝
- 伺服器路徑不能有空格（目前路徑 /big8_disk/... 無空格，正常）

---

## Slide 10: 設定方法 B：手動設定 MCP

設定方法 B：手動設定 MCP

Setup Method B — Manual Configuration

10

方法 1：全域註冊（推薦）

Global registration — shared across all projects

codex mcp add knowledge -- python3 \
  /big8_disk/liaoyoyo2001/knowledge/\
  scripts/mcp/knowledge_server.py

驗證是否成功  Verify

codex mcp list
python3 .../knowledge_server.py --test
# 預期：43/43 passed

■ 所有新專案自動共用  |  One-time setup for all future projects

方法 2：專案層 .mcp.json

Project-level — version-controlled with git

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

■ 適合團隊環境，可提交到 git 版本控制

Suitable for teams; commit .mcp.json to git repo

CCU Bioinformatics Lab  |  2026-03

> **Notes:** 【投影片目的】
說明兩種手動設定 MCP 的方式。

【方法 1：全域（codex mcp add）】
- 適合：個人使用，不想每個專案重複設定
- 儲存位置：~/.config/codex/settings.json（或類似路徑）
- 完整指令（複製用）：
  codex mcp add knowledge -- python3 /big8_disk/liaoyoyo2001/knowledge/scripts/mcp/knowledge_server.py
- 驗證：codex mcp get knowledge 應顯示設定內容

【方法 2：專案層（.mcp.json）】
- 適合：團隊共用，確保所有人使用相同設定
- 建立方式：在專案根目錄建立 .mcp.json（見投影片程式碼）
- 注意：args 中的路徑需使用完整絕對路徑

【兩種方法的選擇建議】
- 個人研究機器 → 方法 1（全域，最方便）
- 多人共用的分析環境 → 方法 2（.mcp.json commit 到 git）
- 兩種可同時存在，.mcp.json 優先於全域設定

【常見問題】
Q: 設定後仍顯示 "tool not found"？
A: 需要重新載入 Claude Code（Ctrl+R 或重啟），MCP server 才會生效。

Q: --test 顯示測試失敗？
A: 確認 knowledge_server.py 的 Python 環境有安裝 fastmcp：pip install fastmcp

---

## Slide 11: 四、對工作流程的實際改善

四、對工作流程的實際改善

Impact on Research Workflow

11

設定前  Before

✗  BAM 路徑在哪？→ 問學長姐

BAM path? → Ask someone

✗  工具參數？→ 翻舊腳本和 GitHub

Tool params? → Dig old scripts / GitHub

✗  哪個版本支援此功能？→ 逐一測試

Which version? → Test one by one

✗  分析前找資料的時間 > 分析本身

Time searching > time analyzing

設定後  After

✓  Claude Code 直接問 → 秒得來源依據的答案

Ask Claude Code → instant answer with source path

✓  AI 附來源文件路徑，可點擊驗證

AI cites source doc — every answer is traceable

✓  中英文查詢皆可（甲基化 = methylation）

Chinese or English queries both work equally

✓  新人 30 分鐘獨立上手，不佔用學長姐時間

New members onboard in 30 min independently

CCU Bioinformatics Lab  |  2026-03

> **Notes:** 【投影片目的】
用 Before/After 對比說明知識庫 + MCP 帶來的實際效益。

【實際對比案例（可現場示範）】
情境：研究員想知道 HCC1395 有哪些甲基化資料

設定前：
1. 打開 /big8_disk/data/ 逐目錄找
2. 問學長姐「ONT_5kHz 有沒有 5mCG？」
3. 查舊腳本看有沒有執行過相關分析
→ 約 5–10 分鐘，且不確定答案是否最新

設定後：
1. 在 Claude Code 輸入：「HCC1395 有哪些甲基化資料？」
2. AI 查詢 knowledge_search("HCC1395 methylation 5mCG")
3. 回傳：ONT_5kHz 和 ONT_Dorado 有 5mCG tag；原始 ONT 無
→ 約 10–15 秒，且附有來源文件路徑供驗證

【量化指標】（Phase 1–3 改善後）
- MCP 測試通過率：43/43（100%）
- 支援搜尋語言：中文 + 英文（Bigram + Stemming）
- 文件涵蓋率：38/38 篇均已索引
- 中英文別名：21 份文件補齊 alias_paths

---

## Slide 12: 近期 Git 更新時序

近期 Git 更新時序

Recent Git Update Timeline

12

2026-02-10

初始建立

Initial Setup

文件架構 · CLAUDE.md AI 指引 · 資料審計

2026-02-27

規範層建立

Standards Layer

09_standards 目錄 · 命名規範 · 決策記錄 log

2026-02-28

Phase 1 — MCP 上線

MCP Gateway Live

FastMCP stdio · 38 篇文件索引 · 21 項測試通過

2026-03-01

Phase 2 — Bug 修復

Bug Fixes A1–A6

09_standards 納入索引 · 12 份 description 改寫

2026-03-01

Phase 3 — 搜尋升級

Search Engine v3

OR+AND bonus · CJK Bigram · Stemming · NOT · 43/43 通過

CCU Bioinformatics Lab  |  2026-03

> **Notes:** 【投影片目的】
說明知識庫從建立到 MCP 上線的演進歷程。

【各里程碑詳細說明】
2026-02-10 初始建立（Commit: 6724b4a, d962755）
  - 建立 01-09 目錄結構
  - CLAUDE.md：AI Agent 使用指引
  - 資料審計報告：確認 4.7TB 資料完整性

2026-02-27 規範層（Commit: e3df65f）
  - 新增 09_standards/ 目錄
  - 07_scripts 重命名為 07_script_docs（區分文件 vs 腳本）
  - 建立 decision-log：記錄架構決策的「單一來源」

2026-02-28 Phase 1 MCP（Commit: 26c95e3, 81e70bb）
  - 實作 knowledge_server.py（FastMCP + stdio）
  - 三個工具：knowledge_search / knowledge_get_doc / knowledge_resolve_path
  - 21 項自動測試全通過

2026-03-01 Phase 2 Bug 修復（Commit: 64698e1）
  - A1: max_results=-1 回傳全部文件的 bug
  - A2-A6: 其他邊界案例修復
  - 12 份文件 description 改寫（模板 → 含觸發詞描述）

2026-03-01 Phase 3 搜尋升級（Commit: 3337993）
  - 搜尋演算法：OR with AND bonus（×1.5 全詞命中）
  - 中文 CJK Bigram 展開
  - 英文規則式 Stemming（無需外部套件）
  - NOT 運算子 (-term)
  - 新增 knowledge_list 工具
  - 測試從 21 項擴充到 43 項，全數通過

---

## Slide 13: 改善成果一覽

改善成果一覽

Summary of Improvements

13

■ 搜尋引擎升級（Phase 3）

Search Engine Upgrade

OR + AND Bonus

全詞命中排名 ×1.5，消除偽匹配

Full match ×1.5; eliminates false positives

CJK Bigram

中文多字詞自動展開搜尋

Auto-expand Chinese multi-char terms

英文 Stemming

詞形還原，不需外部套件

Rule-based morphology, no extra packages

NOT 運算子

-term 排除特定文件

Exclude docs with -term prefix

■ 量化改善成果

Quantified Improvement Results

指標 / Metric

改善前 / Before

改善後 / After

索引文件數
Docs indexed

35 篇

38 篇

MCP 測試項目
MCP test cases

21 項

43 項 ✓

中英文別名
Alias paths

空 /
none

21 份

Description 品質

模板
template

語意化
semantic

搜尋精準度
Search precision

OR only

+AND bonus

CCU Bioinformatics Lab  |  2026-03

> **Notes:** 【投影片目的】
用數字和清單呈現三個 Phase 的具體改善成果。

【量化數字說明】
- 35 → 38 篇文件：加入 09_standards 三份文件（命名規範、MCP onboarding、决策记录）
- 21 → 43 項測試：Phase 3 新增 22 項測試，覆蓋 Bigram、Stemming、NOT、knowledge_list 等功能
- 21 份 alias_paths：讓中英文別名均能找到正確文件（例：「體細胞定相」= somatic-haplotagging）
- Description 語意化：12 份原本是模板說明的文件，改寫為含關鍵觸發詞的描述
  （例：pon_databases 加入 gnomAD、dbSNP 等關鍵詞，讓搜尋時能找到）
- 搜尋精準度：OR only → OR + AND bonus，實測多關鍵字搜尋的誤匹配率大幅下降

【Phase 3 測試覆蓋範圍（43 項）】
- 基礎搜尋（英/中/混合）
- AND bonus 排序驗證
- CJK Bigram 展開驗證
- 英文 Stemming 驗證
- NOT 運算子邊界案例
- knowledge_list tag/status 篩選
- knowledge_resolve_path 中英文別名解析
- knowledge_get_doc 章節擷取
- max_results 邊界（0、1、負數）

---

## Slide 14: 立即開始使用

立即開始使用

Get Started Now

14

手動瀏覽

Manual Browsing

/big8_disk/liaoyoyo2001/knowledge/

從 README.md 開始，即刻掌握全局

Start with README.md — overview in 30 min

AI 查詢

AI-Assisted Query

knowledge_search("你的問題")

中英文均可，回答附來源路徑

Chinese or English — answers include source

提交問題

Report or Improve

issues/ 目錄（Markdown）

歡迎回報錯誤或提交改善 PR

Bug reports and PRs are always welcome

CCU Bioinformatics Lab  |  Knowledge Base v2  |  2026-03

> **Notes:** 【投影片目的】
提供行動呼籲（CTA），讓聽眾知道下一步。

【三個行動選項說明】
1. 手動瀏覽
   - 路徑：/big8_disk/liaoyoyo2001/knowledge/
   - 建議第一個打開：README.md
   - 無需任何設定，立即可用

2. AI 查詢
   - 前提：完成 MCP 設定（方法 A 或 B）
   - 直接在 Claude Code 輸入問題
   - 範例：「HCC1395 有 5mCG 甲基化資料嗎？」

3. 提交問題 / 改善
   - 發現文件有誤：在 issues/ 目錄建立 Markdown 問題記錄
   - 想更新內容：直接編輯 .md 文件並提交 PR
   - GitHub: CCU-Bioinformatics-Lab/knowledge

【結語建議】
「知識庫是活文件。任何人發現錯誤或想補充，歡迎直接修改並提 PR。
 AI 查到錯誤答案時，最好的修復方式就是去更新源頭的知識文件。」

【問答常見問題清單】
Q1: 設定 MCP 後多久生效？→ 重新載入 Claude Code 後立即生效
Q2: 知識庫多久更新一次？→ 無固定週期，隨研究進展更新
Q3: 可以加入自己的文件嗎？→ 可以，按照 09_standards 的格式規範撰寫即可
Q4: AI 回答有錯怎麼辦？→ 找到對應的 .md 文件修正，來源準確了 AI 就準確
Q5: 在 Windows / Mac 開啟 Claude Code 可以用嗎？→ 可以，MCP 設定跨平台

---
