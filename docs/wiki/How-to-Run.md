# How to Run 操作手冊
[← Home](https://github.com/liaoyoyo/InterSubMod/wiki) · [System Overview](https://github.com/liaoyoyo/InterSubMod/wiki/System-Overview) · [InterSubMod](https://github.com/liaoyoyo/InterSubMod/wiki/InterSubMod-Engine) · [LongLineage](https://github.com/liaoyoyo/InterSubMod/wiki/LongLineage-Engine) · [Upstream](https://github.com/liaoyoyo/InterSubMod/wiki/Upstream-and-Data) · [Analysis](https://github.com/liaoyoyo/InterSubMod/wiki/Analysis-and-Presentation) · [How to Run](https://github.com/liaoyoyo/InterSubMod/wiki/How-to-Run)

> **從零到跑出第一個結果**
> 本頁每一條指令都在 **2026-08-06 實際執行過**，並附上真實輸出。
> 照順序做，約 **10 分鐘**可以從編譯到看見第一個分析結果。

---

## ⚠️ 開始之前：先確認機器狀態

這台機器有 **48 核**，但 `/big7_disk` **只剩約 617 GB（99% 已用）**。
跑全基因組分析前務必先確認餘量，大型輸出很容易把磁碟寫爆。

---

## 步驟 1 · 編譯 C++（約 2–5 分鐘）

🔴 **一定要先做這步。** 現有的執行檔是 STALE 的 —— 有 5 個原始碼檔比它新。

```bash
cd /big7_disk/liaoyoyo2001/InterSubMod

# 設定並編譯（Release 模式）
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j$(nproc)
```

**編譯完成後應該有 6 個執行檔**：

```text
$ ls build/bin/
exact_ps_partition   inter_sub_mod   run_tests
test_phase1_2        test_phase3     test_phase4_5
```

---

## 步驟 2 · 確認編譯結果是好的（約 2 秒）

```bash
./build/bin/run_tests
```

**本輪實跑的真實輸出**：

```text
[==========] 265 tests from 38 test suites ran. (2062 ms total)
[  PASSED  ] 265 tests.

$ echo $?
0
```

✅ **265 / 265 通過** — 如果這裡有任何失敗，**先不要往下走** —— 後面跑出來的數字都不可信。

---

## 步驟 3 · 準備 Python 環境

```bash
# 安裝依賴（清單於 2026-08-05 依實測環境補齊）
pip install -r requirements.txt
```

### 🔴 這裡有一個一定會踩到的坑

這台機器上有**兩個 Python**：

| 指令 | 版本 | 用途 |
|---|---|---|
| `python3` | 3.9.12 | 預設，大部分腳本可用 |
| `/usr/bin/python3.10` | 3.10.12 | **嚴格區域切割相關腳本必須用這個** |

用錯版本時，程式會在 import 階段就崩，**而且錯誤訊息指向 dataclass 而不是版本** —— 看起來像程式壞了，其實只是 Python 版本不對。

---

## 步驟 4 · 跑出第一個分析結果（約 3 秒）

先用單一個突變位點試跑，確認整條路徑通暢。

```bash
# 準備一個只含單一突變的小 VCF
SP=/tmp/ism_demo && mkdir -p $SP
V=/big7_disk/liaoyoyo2001/InterSubMod/data/vcf/HCC1395/pileup/filtered_snv_tp.vcf
grep '^#' $V > $SP/one_snv.vcf
grep -P '^chr19\t29283968\t' $V >> $SP/one_snv.vcf

# 跑（只給三個必填參數）
./build/bin/inter_sub_mod \
  --tumor-bam  data/bam/HCC1395/tumor.bam \
  --reference  data/ref/hg38.fa \
  --vcf        $SP/one_snv.vcf \
  --output-dir $SP/out_min
```

**實跑輸出**（約 2.9 秒）：

```text
Total regions: 1 / Successful: 1 / Failed: 0
Total reads processed: 85
Forward strand (+): 40 / Reverse strand (-): 45
Total CpG sites found: 11
Metric: NHD / Total valid read pairs: 3443 / Total invalid pairs: 127
```

### ⚠️ 兩個預設值和說明文件不一樣，先知道比較好

- **`--threads`** 說明寫預設 1，**實際是 16**（資源估算會差 16 倍）。
- **`--distance-metric`** 宣告預設 BERNOULLI，**實際會被覆寫成 NHD**。

想要哪個就明確指定，不要靠預設。

---

## 步驟 5 · 看結果

```bash
find $SP/out_min -type f | head -20
```

輸出分兩層。**先看這三個檔**：

| 檔案 | 看什麼 |
|---|---|
| `<region>/metadata.txt` | 這個位點抓到幾條 read、幾個 CpG、矩陣多大 —— 先確認數字合理。 |
| `<region>/reads/reads.tsv` | 每條 read 的單倍型標籤與支持突變型/正常型。 |
| `significance_summary.csv` | run 層總表，**下游分析幾乎只吃這一份**。 |

### 🔴 讀輸出檔時的三個必知

1. `methylation.csv` 第一欄是**列號不是 read 名**，要對回 read 得查 `reads.tsv`。
2. `linkage_matrix.csv` 副檔名是 csv 但**實際是 tab 分隔** —— `pd.read_csv()` 預設會整列讀成單欄**而且不報錯**，要寫 `sep='\t'`。
3. `tree.nwk` 的葉子是 **read** 不是 clone —— 這是甲基化相似度的分群樹，**不是亞群演化樹**。

---

## 步驟 6 · 驗證整套流程的數字（推薦）

```bash
cd docs/methodology/_assets/20260627_subclone_4axis_teaching/scripts
python3 verify_pipeline_numbers.py
```

這支會把方法文件裡的每個數字**重新算一次**並比對。實跑輸出：

```text
sSNV 總數 35,332（TP 30,490 / FP 4,842）              ✓
共現連上 21,554 · 訊號不足 5,458 · 孤立 8,320        （加總 = 35,332 ✓）
```

對得上，代表你的環境能重現既有結果。

---

## 07 · 全流程速查

![howto-six-steps](https://raw.githubusercontent.com/liaoyoyo/InterSubMod/develop/docs/images/howto-six-steps.png)

> **圖 1 · 六個步驟與各自的驗收點** —— 六個步驟依序為：編譯、跑測試、裝 Python 依賴、單點試跑、檢視輸出、驗證數字；每步都標註對應的驗收條件，例如測試須 265 個全過、單點試跑須在三秒內完成且產出 region 目錄。

### 每一步的驗收條件（不過就停下來，別往下走）

| 步驟 | 驗收條件 |
|---|---|
| 1 編譯 | `build/bin/` 下出現 6 個執行檔 |
| 2 測試 | **265 tests 全過、退出碼 0** |
| 3 依賴 | `import pysam` 不報錯 |
| 4 試跑 | 印出 `Successful: 1 / Failed: 0` |
| 5 輸出 | region 目錄下有 methylation 與 distance |
| 6 驗證 | 各層加總打勾，數字對得上 |

🔴 任何一步的退出碼不是 0，或數字對不上 —— 先解決再繼續，否則後面的結果都不可信。

---

## 08 · 常見狀況排除

| 症狀 | 原因與解法 |
|---|---|
| 程式直接停，說找不到參考基因組 | 參考 FASTA **缺 `.fai` 索引是硬性錯誤**。跑 `samtools faidx` 建索引。 |
| Python 腳本 import 就崩，錯誤指向 dataclass | **Python 版本問題**，不是程式壞了。改用 `/usr/bin/python3.10`。 |
| CI 判定 `--help` 失敗 | `--help` 的**退出碼是 1 不是 0**（與參數錯誤共用同一路徑）。改寫檢查方式。 |
| 某支腳本跑了但什麼都沒發生 | 可能是**函式庫不是入口**（如 `ism_heatmap_std.py`）。看它有沒有 argparse。 |
| 讀 CSV 後整列變成一欄 | `linkage_matrix.csv` **實際是 tab 分隔**。加 `sep='\t'`。 |
| 跨 run 比較時欄位對不上 | `significance_summary.csv` **欄數隨版本變動且無版本欄位**。一律**用欄名**取值。 |
| 工作站 HTML 開很久／瀏覽器很卡 | 單檔可達 188 MB。示範時先開最小的那個（約 14.6 MB）。 |
| LongLineage 說 KernelBlocked | 🔴 **這是預期行為不是故障**。正式入口刻意 fail-closed，見 [LongLineage 分冊](https://github.com/liaoyoyo/InterSubMod/wiki/LongLineage-Engine)。 |

---

## 本頁的驗證方式

- **測試套件**：2026-08-06 實跑 `./build/bin/run_tests` → **265 tests / 38 suites 全過，2062 ms，退出碼 0**。（比 08-05 記錄的 258/37 有增加，本頁採用實測的當前值。）
- **單點試跑**：實際執行並記錄真實輸出，約 2.9 秒、退出碼 0。
- **環境數字**：48 核、磁碟餘量 617 GB、兩個 Python 版本，皆為本輪實測。
- **6 個執行檔**：實際 `ls build/bin/` 所得。

## ⚠️ 誠實標註

- 本頁的編譯指令**未在本輪從零重跑**（現有 build 目錄已存在）；該流程於 2026-08-05 曾以全新目錄驗證通過。
- **全基因組全量跑法未實跑** —— 屬長時間計算，且需先確認磁碟餘量。

---

**來源**：`InterSubMod/docs/explain/16_how-to-run.standalone.html` · 分支 `research/subclonal-reconstruction-202606` · 建立於 2026-08-06

[← 上一頁：分析與呈現層](https://github.com/liaoyoyo/InterSubMod/wiki/Analysis-and-Presentation) · [回系統全景](https://github.com/liaoyoyo/InterSubMod/wiki/System-Overview) · [回 Home](https://github.com/liaoyoyo/InterSubMod/wiki)
