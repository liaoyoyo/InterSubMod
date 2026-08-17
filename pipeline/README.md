# SUPERSEDED／HISTORICAL — InterSubMod × LongLineage 早期整合草案

> **不可作為新人執行入口。** 本檔保留早期 feature-preview 介面與 tag 語意作 provenance；
> 現行入口是 [`docs/handoff/20260813_完整研究資料與軟體交接_01/00_INDEX.md`](../docs/handoff/20260813_完整研究資料與軟體交接_01/00_INDEX.md)
> 與 [`QUICKSTART.md`](../QUICKSTART.md)。LongLineage private baseline/main snapshot `5daf50f`
> 沒有 tagged-BAM writer；private public-preview candidate `b9aaa12` 雖包含
> `longlineage-tag-bam`，仍為 **NOT_READY／non-production**，P3/P4/P5/P7/P8 均 `BLOCKED`。
> Feature presence 不等於 supported workflow，production `run` 仍刻意 fail closed。

## 這份文件回答什麼

InterSubMod 分析「甲基化 × read 標籤」的 association。核心流程不要求 lineage tags；
本檔只解釋 private preview 曾提出的 optional tag contract，不證明該整合流程可公開重跑。

---

## 1. InterSubMod 需要什麼樣的輸入 BAM

| 必要 | 說明 |
|---|---|
| `MM` / `ML` tag | 甲基化呼叫（dorado 的 `C+m?`）。**沒有它 ISM 無法做任何甲基分析** |
| `HP:Z` | LongPhase-S 可能出現 `. / 1 / 2 / 3 / 4 / 1-1 / 2-1 / 1-2 / 2-2`；frozen baseline `ddd8909a` 的 ISM 有效 HP 分組只辨識 `1/2/1-1/2-1/3` 與 `0/unphased`，`4/1-2/2-2` 落 `HP_OTHER`／被 HP 檢定排除。實跑前必做 HP value census。 |
| 座標排序 + `.bai` | 一般 BAM 要求 |

| 歷史 preview 可選欄位（尚非 supported public workflow） | 說明 |
|---|---|
| `lc:Z` | lineage component（unit_id） |
| `lu:Z` | lineage block（block_id） |
| `lv:Z` | **階層路徑**，如 `HP2-1-1`；`+` 後綴表示「該節點或其後代」 |
| `lp:Z` | 該 read 在 block 上的觀察 pattern（`R`/`A`/`X`） |
| `lo:Z` | 突變順序，如 `2034084>2057742` |
| `ls:A` | 標籤可信度：`U` 唯一／`M` 拓撲並列／`P` 部分覆蓋／`A` 拓撲未定 |

**沒有 `lc/lu/lv` 也能跑** —— ISM 會退回只用 `HP` 與 `ALT` 兩個軸。

---

## 2. Tagged-BAM 的 revision-scoped 狀態

以下是 private candidate `b9aaa12` 的**預覽 artifact 名稱**，不是 production contract：
```
<輸出目錄>/paths/*.unit_lineage_paths.tsv.gz          階層路徑 + 突變順序
<輸出目錄>/assign/*.read_lineage_assignments.tsv.gz   read → block 指派
<輸出目錄>/bam/*.lineage_tagged.bam                   preview-only tagged BAM
```

不得依此段自行執行 private repo 或宣稱輸出已驗證；請先查交接包的 capability matrix、
blocker set 與 site preflight。若 P3/P4/P5/P7/P8 任一仍 blocked，這條 chain 就不能升格為 supported。

---

## 3. InterSubMod 的 supported 執行路徑

Supported core smoke／site-profile 命令只維護於 `QUICKSTART.md` 與 handoff package；
目前只把 HP/ALT 與 read-level methylation association 列為 supported。`lv/lc/lu` 整合屬
private preview，不能從本歷史檔複製命令執行。

### 兩個 lineage 相關參數

| 參數 | 值 | 說明 |
|---|---|---|
| `--group-by-tag` | `HP` / `ALT` / `lc` / `lu` / `lv`（feature syntax可重複／逗號分隔） | **歷史文字曾誤稱多軸並行。** Inspected feature loop 目前只評估第一個 lineage axis；後續 lineage-axis values 不會各自獨立產生檢定。HP／ALT 是另一路 grouping，不能用來推論所有 lineage axes 已 multi-axis。此介面不在 `ddd8909a` supported baseline。 |
| `--require-tag-status` | `U`（預設）/ `UM` / `any` | 使用 lineage 軸時的最低可信度 |

⚠ **樣本量警告**：HCC1395 全基因組的 `ls` 分佈為
`P 88.3%｜U 7.7%｜M 3.9%｜A 0.1%`。
用 `--require-tag-status U` 時可用 read 大幅減少，許多 region 會因樣本不足而無法檢定。
這是 ONT read 長度 vs block 跨度的資料本質，不是設定錯誤。

### 輸出

| 檔案 | 內容 |
|---|---|
| `significance_summary.csv` | per-region 統計；feature lineage-axis 僅第一軸被評估，不可稱完整多軸輸出 |
| `<region>/reads/reads.tsv` | per-read，**後 7 欄**為 lineage 資訊（無標籤時寫 `.`） |
| `<region>/methylation/methylation.csv` | read × CpG 甲基矩陣 |

---

## 4. 整合觀察報告（HTML；歷史構想）

早期草案曾規劃由 LongLineage presentation 整合；它未被本 handoff 驗證為 public supported
workflow。公開 HTML 是 validated資料的呈現層，不重算 science；可用頁面與 receipt 請從 handoff index 導航。

---

## 5. 三專案的邊界

| 專案 | 職責 | 甲基資料 |
|---|---|---|
| **LongLineage private preview** | read linkage／candidate reconstruction；tagged-BAM 只存在 candidate revision | **不需要** |
| **InterSubMod** | 依 BAM 標籤做甲基關聯分析 | **必要** |
| **LongLineage/presentation** | 整合兩者輸出產 HTML | 可缺（降級） |

⇒ 兩條 provenance chain 必須分開驗證；不得因 private candidate 有某個檔名，就宣稱完整鏈已 supported 或 production-ready。
