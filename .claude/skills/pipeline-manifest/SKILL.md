---
name: pipeline-manifest
description: 可重現性 provenance 清單 — 掃描分析腳本，建立 script → inputs → outputs → figures/tables → reports 的對應 DAG，偵測 orphan（無腳本產生的圖、無人消費的輸出、報告引用但不存在的檔案）。USE WHEN：「pipeline manifest」「provenance map」「腳本對應圖表」「哪個腳本產生這張圖」「reproducibility map」「投稿前可重現性檢查」「supplementary methods 表」、cycle 收尾前盤點輸出鏈、論文 submission 前產 reproducibility 附錄。SKIP WHEN 單一腳本除錯（直接 Read）、純檔案組織檢查（用 /data-audit）、程式碼正確性驗證（用 /verification-loop）、純 build / commit / docs 寫作、無 scripts→figures 鏈的任務。
allowed-tools: Read, Write, Glob, Grep, Bash(ls:*), Bash(find:*), Bash(grep:*), Bash(wc:*), Bash(head:*)
user-invocable: true
---

# Pipeline Manifest — Reproducibility Provenance Map

**核心原則**：每一張進報告的圖 / 每一個 TSV 結論數字，都必須能反查「哪個腳本、吃哪個輸入、在哪個 commit 產生」。reviewer 與未來的你都會問這個。本 skill 把分散在 `scripts/` / `research/` / `docs/` 的隱性對應關係**顯性化成一張表**，並抓出斷鏈。

> **填補空缺**：對齊社群 `claude-research/pipeline-manifest`（map scripts → inputs/outputs → paper figures/tables）。本專案原本只有 `/data-audit`（查檔案組織完整性）但**不建立 script→figure provenance DAG** — 本 skill 補上這層。

---

## Phase & Chain Position

- **Phase**：P5-P6（cycle 收尾 / 報告產出前）與**論文 submission 前**為主要觸發點；亦可任何時間 ad-hoc 盤點。
- **Chain 上游**：`/feature-layered-observation`（P3）/ `/results-analysis` 產出 figures + TSV 後。
- **Chain 下游**：`/conclude-research`（P6 收尾引用本 manifest 當 reproducibility 附錄）/ `/results-report`（decision section 引用 provenance）/ 論文 supplementary methods。
- **與 `/data-audit` 分工**：data-audit = 「檔案存在嗎、命名合規嗎、被索引嗎」（組織層）；pipeline-manifest = 「這個輸出**從哪來、誰消費它**」（因果鏈層）。兩者互補不重疊。

## Dependencies

- **Uses**：Glob / Grep / Read（掃描）；Bash(find/grep)（建索引）。
- **Used-by**：`/conclude-research`（附錄）、`/results-report`（provenance）、論文投稿 reproducibility 檢查。
- **Reads**：`scripts/**/*.py`、`scripts/**/*.sh`、`research/**/figures/*.png`、`research/**/*.tsv`、`docs/experiments/**/*.md`、`research/autoresearch/evidence_ledger.jsonl`。
- **Writes**：`<target_dir>/PIPELINE_MANIFEST.md` + `<target_dir>/pipeline_manifest.tsv`（機器可讀）。

## Failure Mode & Diagnostics

| 失效模式 | 症狀 | 診斷 |
|---|---|---|
| 腳本用變數路徑（`f"{outdir}/x.png"`） | grep 抓不到字面路徑 | 標 `[DYNAMIC]`，人工確認 outdir 來源；不臆測 |
| 圖由 notebook / 互動 session 產生 | 圖存在但無對應 .py | 標 `ORPHAN_FIGURE`，要求補腳本或註記「manual」 |
| 報告引用已刪除的圖 | .md 有 `![](...)` 但檔案不存在 | 標 `BROKEN_REF`，列出引用行號 |
| 跨 disk 絕對路徑輸入（/big8 NFS） | 輸入在 read-only NFS | 標 `EXTERNAL_INPUT`，記錄但不視為斷鏈 |

---

## 執行步驟（Step 1-6）

### Step 1 — 界定範圍
確認 `<target_dir>`（預設整個 repo，或某 cycle 的 `research/<project>/cycleN/`）。問清楚是「全 repo 盤點」還是「單 cycle/單報告的 provenance」。

### Step 2 — 建立三層索引
```
A. 腳本層：Glob scripts/**/*.{py,sh} + research/**/*.py
   → 對每個腳本 Grep 抓 input/output 字面路徑（open()/read_csv/savefig/to_csv/--output/-o）
B. 產物層：Glob research/**/figures/*.png + **/*.tsv + **/*.csv
C. 消費層：Grep docs/**/*.md + research/**/*.md 的 ![](...) 圖片引用 + 表格 source 註記
```

### Step 3 — 連邊（建 DAG）
對每個產物，反查 Step 2A 哪個腳本的 output 命中它；對每張被報告引用的圖，連到產生它的腳本。輸出三元組 `script → output → report`。

### Step 4 — 抓 orphan / 斷鏈
- `ORPHAN_FIGURE`：圖存在，無腳本產生 → 標記
- `UNUSED_OUTPUT`：腳本產出，無報告/下游消費 → 標記（可能是中間檔，註記即可）
- `BROKEN_REF`：報告引用，檔案不存在 → 🔴 必修
- `DYNAMIC` / `EXTERNAL_INPUT`：見上表

### Step 5 — 標 provenance tier + commit
每條鏈標：腳本最後修改 commit（`git log -1 --format=%h -- <script>`，若 git 可用）、輸入資料 `dataset_id`（對照 evidence_ledger）、產生日期。對齊 `/scientific-rigor §8.4`（provenance footer 必含 commit hash）。

### Step 6 — 產出
寫兩份：
- `PIPELINE_MANIFEST.md`：人讀，含 DAG 表 + orphan 清單 + 🔴 BROKEN_REF 待修清單 + 5-行 summary。
- `pipeline_manifest.tsv`：機器讀，欄位 `report\tfigure_or_table\tproducing_script\tinput_data\tscript_commit\tstatus`。

---

## 輸出範本（PIPELINE_MANIFEST.md 骨架）

```markdown
# Pipeline Manifest — <target> (<date>)

## Summary（≤5 行）
- 涵蓋 N 個腳本 / M 張圖 / K 個 TSV / R 份報告
- 完整鏈：X 條 | ORPHAN_FIGURE：Y | BROKEN_REF：Z 🔴 | UNUSED_OUTPUT：W

## Provenance DAG
| Report | Figure/Table | Producing Script | Input Data | Script Commit | Status |
|--------|-------------|------------------|------------|---------------|--------|
| ...    | fig_3a.png  | scripts/analysis/o12.py | HCC1395_master.tsv | a3f9c1 | ✅ COMPLETE |

## 🔴 BROKEN_REF（必修）
- docs/.../report.md:142 引用 `fig_missing.png` — 檔案不存在

## ORPHAN_FIGURE（補腳本或註記 manual）
- research/.../figures/adhoc_plot.png — 無對應腳本
```

---

## 反模式（紅線）

- ❌ 臆測動態路徑指向哪個檔 → 標 `[DYNAMIC]` 交人工，不猜
- ❌ 把 read-only NFS 外部輸入當斷鏈 → 那是合法 `EXTERNAL_INPUT`
- ❌ 只列完整鏈不列 orphan → orphan 才是 reproducibility 風險所在
- ❌ 省略 commit hash → 違反 §8.4，reviewer 無法重現

## 嚴謹度繼承（/scientific-rigor）

- **§8.4 Provenance**：每條鏈含 commit hash + dataset_id + date
- **§2 Evidence Tier**：BROKEN_REF 的結論數字降為 L4（無法重現）直到修復
- **§7.2 可重現性 7 項**：本 manifest 直接產出其中「腳本→輸出對應」與「輸入資料溯源」兩項
