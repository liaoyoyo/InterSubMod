# InterSubMod 的前置步驟與執行流程

> 本目錄**只有說明，沒有程式**。
> 產生輸入 BAM 的實作在 **LongLineage**，不屬於 InterSubMod。

## 這份文件回答什麼

InterSubMod 分析「甲基化 × read 標籤」的關聯。要跑它，你的 BAM 必須帶上
lineage 標籤。這份文件說明那些標籤是什麼、從哪來、怎麼準備。

---

## 1. InterSubMod 需要什麼樣的輸入 BAM

| 必要 | 說明 |
|---|---|
| `MM` / `ML` tag | 甲基化呼叫（dorado 的 `C+m?`）。**沒有它 ISM 無法做任何甲基分析** |
| `HP:Z` | 單倍型標籤（九態：`1`/`2`/`1-1`/`2-1`/`3`） |
| 座標排序 + `.bai` | 一般 BAM 要求 |

| 可選（啟用 lineage 軸才需要） | 說明 |
|---|---|
| `lc:Z` | lineage component（unit_id） |
| `lu:Z` | lineage block（block_id） |
| `lv:Z` | **階層路徑**，如 `HP2-1-1`；`+` 後綴表示「該節點或其後代」 |
| `lp:Z` | 該 read 在 block 上的觀察 pattern（`R`/`A`/`X`） |
| `lo:Z` | 突變順序，如 `2034084>2057742` |
| `ls:A` | 標籤可信度：`U` 唯一／`M` 拓撲並列／`P` 部分覆蓋／`A` 拓撲未定 |

**沒有 `lc/lu/lv` 也能跑** —— ISM 會退回只用 `HP` 與 `ALT` 兩個軸。

---

## 2. 帶標籤的 BAM 從哪來

由 **LongLineage** 產生（`/big7_disk/liaoyoyo2001/LongLineage`）：

```bash
cd /big7_disk/liaoyoyo2001/LongLineage
/usr/bin/cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
/usr/bin/cmake --build build -j16

LONGLINEAGE_BUILD=$PWD/build bash scripts/run_sample.sh \
    --sample HCC1395 \
    --in-bam  <已 haplotag 的 BAM 或 raw BAM> \
    --out-root <輸出目錄> \
    --threads 16
```

產出：
```
<輸出目錄>/paths/*.unit_lineage_paths.tsv.gz          階層路徑 + 突變順序
<輸出目錄>/assign/*.read_lineage_assignments.tsv.gz   read → block 指派
<輸出目錄>/bam/*.lineage_tagged.bam                   ← 這就是 ISM 的輸入
```

⚠ 詳細規格見 `LongLineage/docs/lineage/`（9 份文件，含 `HIERARCHICAL_TAG_SPEC.md`
與 `KNOWN_ISSUES.md`）。

---

## 3. 執行 InterSubMod

```bash
cd /big7_disk/liaoyoyo2001/InterSubMod
./build/bin/inter_sub_mod \
    -t <輸出目錄>/bam/HCC1395.chr1.lineage_tagged.bam \
    -r <reference.fa> \
    -v <PASS.vcf.gz> \
    -o <ISM 輸出目錄> \
    --group-by-tag HP,ALT,lv \
    --require-tag-status U \
    -j 16 -w 5000
```

### 兩個 lineage 相關參數

| 參數 | 值 | 說明 |
|---|---|---|
| `--group-by-tag` | `HP` / `ALT` / `lc` / `lu` / `lv`（可多選，逗號分隔） | 要檢定哪些分組軸。多軸並行，讓每個位點能回答「甲基差異與**哪個**軸有關」 |
| `--require-tag-status` | `U`（預設）/ `UM` / `any` | 使用 lineage 軸時的最低可信度 |

⚠ **樣本量警告**：HCC1395 全基因組的 `ls` 分佈為
`P 88.3%｜U 7.7%｜M 3.9%｜A 0.1%`。
用 `--require-tag-status U` 時可用 read 大幅減少，許多 region 會因樣本不足而無法檢定。
這是 ONT read 長度 vs block 跨度的資料本質，不是設定錯誤。

### 輸出

| 檔案 | 內容 |
|---|---|
| `significance_summary.csv` | per-region 多軸統計 |
| `<region>/reads/reads.tsv` | per-read，**後 7 欄**為 lineage 資訊（無標籤時寫 `.`） |
| `<region>/methylation/methylation.csv` | read × CpG 甲基矩陣 |

---

## 4. 整合觀察報告（HTML）

由 LongLineage 的呈現層產生，**不在 InterSubMod**：

```bash
python3 /big7_disk/liaoyoyo2001/LongLineage/presentation/build_lineage_report.py \
    --out-root <LongLineage 輸出目錄> --sample HCC1395 \
    --ism-summary <ISM 輸出目錄>/.../significance_summary.csv \
    --output HCC1395.report.html
```

缺任何輸入時，對應面板會標示「不可用 + 原因」，其餘照常產出。

---

## 5. 三專案的邊界

| 專案 | 職責 | 甲基資料 |
|---|---|---|
| **LongLineage** | 純遺傳 read linkage → 拓撲 → tagged BAM | **不需要** |
| **InterSubMod** | 依 BAM 標籤做甲基關聯分析 | **必要** |
| **LongLineage/presentation** | 整合兩者輸出產 HTML | 可缺（降級） |

⇒ 甲基不一定有，所以 LongLineage 是基礎層；InterSubMod 是可選的加值分析。
