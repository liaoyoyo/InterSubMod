---
name: data-audit
description: 研究輸出的組織完整性檢核 — 檢查圖片連結、索引覆蓋、命名合規、元數據、gitignore、散落檔案。USE WHEN：研究完成後整理 docs/experiments/ 或 docs/reports/ 輸出、「檢查檔案組織」「audit」。涉及 .md/.png/.csv 檔案。 SKIP WHEN 單檔案編輯無組織變動、程式碼級驗證（用 verification-loop）、純 build、純 .py 開發、不涉及 docs/experiments 或 docs/reports 結構變更。
allowed-tools: ["Bash", "Glob", "Grep", "Read"]
user-invocable: true
tags: ["audit", "organization", "research"]
---

# 資料整理檢核 (Data Audit)

對 InterSubMod 專案的研究輸出進行 6 項組織完整性檢查。依序執行以下檢查，輸出結構化報告。

## 規範參考

檢查標準定義於 `AGENTS.md` 的「研究輸出組織規範」和「主要查詢路徑與重點資訊」區段。執行前先讀取 `AGENTS.md` 確認最新規範。

## 檢查項目

### 1. 圖片連結完整性

掃描所有 `.md` 檔案中的圖片引用（`![...](path)` 格式），驗證：
- 目標檔案是否存在
- 是否使用相對路徑（禁止絕對路徑）
- 相對路徑深度是否 <= 2 層（禁止 `../../../`）

```bash
# 找出所有圖片引用
grep -rn '!\[.*\](.*\.png\|.*\.svg\|.*\.jpg)' docs/ research/ --include='*.md'
# 對每個引用，從 .md 所在目錄計算實際路徑，檢查檔案是否存在
```

**輸出格式**：
```
[圖片連結] 總引用: N, 存在: N, 斷裂: N, 絕對路徑: N, 過深: N
  斷裂: file.md:line -> missing_image.png
  絕對路徑: file.md:line -> /abs/path/image.png
  過深: file.md:line -> ../../../too/deep.png
```

### 2. 索引覆蓋率

掃描 `research/` 和 `docs/reports/` 下的所有 `.md` 檔案，確認每個都被至少一個索引檔引用：
- `research/{study}/README.md`
- `docs/reports/{topic}/00_INDEX.md`
- `docs/experiments/INDEX.md`
- `docs/README.md`

```bash
# 列出所有索引檔
find research/ docs/reports/ -name 'README.md' -o -name '00_INDEX.md' -o -name 'INDEX.md'
# 列出所有 .md 檔案
find research/ docs/reports/ -name '*.md' ! -name 'README.md' ! -name '00_INDEX.md' ! -name 'INDEX.md'
# 對每個非索引 .md，檢查是否被某個索引引用
```

**輸出格式**：
```
[索引覆蓋] 總檔案: N, 已索引: N, 未索引: N (覆蓋率 XX%)
  未索引: path/to/orphan_file.md
```

### 3. 命名合規

檢查檔案命名是否符合規範：
- `.md` 檔案：`{YYYYMMDD}_{中文說明}_{NN}.md`（INDEX/README 例外）
- `.png` 檔案：`{NN}_{英文描述}.png`
- 圖片子目錄：必須叫 `figures/`（不允許 `images/`、`assets/`、`plots/`）

```bash
# 找出不符合命名規範的 .md（排除 README/INDEX）
find docs/ research/ -name '*.md' ! -name 'README.md' ! -name 'INDEX.md' ! -name '00_INDEX.md' | grep -v -E '/[0-9]{8}_.*_[0-9]{2}\.md$'
# 找出非標準圖片子目錄
find docs/ research/ -type d \( -name 'images' -o -name 'assets' -o -name 'plots' \)
```

**輸出格式**：
```
[命名合規] .md 合規: N/N (XX%), .png 合規: N/N (XX%)
  不合規 .md: path/to/bad_name.md
  非標準目錄: path/to/images/
```

### 4. 元數據存在

檢查每個 `.md` 檔案開頭是否有 HTML 註解格式的元數據（`<!-- 建立時間:`）。

```bash
# 對每個 .md 檔案，檢查前 10 行是否包含 '建立時間'
find docs/ research/ -name '*.md' -exec sh -c 'head -10 "$1" | grep -q "建立時間" || echo "$1"' _ {} \;
```

**輸出格式**：
```
[元數據] 有元數據: N/N (XX%)
  缺少: path/to/no_metadata.md
```

### 5. gitignore 覆蓋

確認 `.gitignore` 包含以下規則：
- `research/*/figures/`
- `research/*/data/`
- `/*.pdf`

```bash
grep -n 'research/\*/figures/' .gitignore
grep -n 'research/\*/data/' .gitignore
grep -n '/\*.pdf' .gitignore
```

**輸出格式**：
```
[gitignore] 規則覆蓋: N/3
  缺少: research/*/figures/
```

### 6. 散落檔案偵測

偵測不應存在於當前位置的檔案/目錄：
- 根目錄的 `.pdf` 檔案
- 非標準頂層目錄（`images/`、`Testing/`）
- `docs/` 或 `research/` 外的散落 `.png` 檔案

```bash
ls *.pdf 2>/dev/null
ls -d images/ Testing/ 2>/dev/null
find . -maxdepth 1 -name '*.png' 2>/dev/null
```

**輸出格式**：
```
[散落檔案] 問題數: N
  根目錄 PDF: deck_v2.pdf
  非標準目錄: images/, Testing/
```

## 彙總報告格式

執行完 6 項檢查後，輸出彙總：

```
========== 資料整理檢核報告 ==========
日期: YYYY-MM-DD

1. [圖片連結]  通過/警告  (斷裂: N, 絕對路徑: N, 過深: N)
2. [索引覆蓋]  通過/警告  (覆蓋率: XX%)
3. [命名合規]  通過/警告  (.md: XX%, .png: XX%)
4. [元數據]    通過/警告  (覆蓋率: XX%)
5. [gitignore] 通過/警告  (規則: N/3)
6. [散落檔案]  通過/警告  (問題數: N)

總評: X/6 通過
=====================================
```

通過標準：
- 圖片連結：0 個斷裂 + 0 個絕對路徑
- 索引覆蓋：> 90%
- 命名合規：> 80%
- 元數據：> 80%
- gitignore：3/3 規則
- 散落檔案：0 個問題
