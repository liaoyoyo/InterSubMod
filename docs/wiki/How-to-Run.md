# How to Run 操作手冊
[← Home](https://github.com/liaoyoyo/InterSubMod/wiki) · [System Overview](https://github.com/liaoyoyo/InterSubMod/wiki/System-Overview) · [InterSubMod](https://github.com/liaoyoyo/InterSubMod/wiki/InterSubMod-Engine) · [LongLineage](https://github.com/liaoyoyo/InterSubMod/wiki/LongLineage-Engine) · [Upstream](https://github.com/liaoyoyo/InterSubMod/wiki/Upstream-and-Data) · [Analysis](https://github.com/liaoyoyo/InterSubMod/wiki/Analysis-and-Presentation) · [How to Run](https://github.com/liaoyoyo/InterSubMod/wiki/How-to-Run)

> **從零到跑出第一個結果**
> 編譯與測試可由 public clone 執行；真正的 `inter_sub_mod` 分析另需使用者合法取得的 BAM、FASTA 與 VCF。
> 下方輸出是 historical internal HCC1395 single-locus smoke receipt；本頁沒有完整釘定
> hardware、input locus、commit 與 date，因此不報秒數，也不把它當成一般 runtime。

> **公開可重現性邊界**：Git 物件不含 `data/bam/HCC1395/tumor.bam`、`data/ref/hg38.fa` 或示例 `one_snv.vcf`。在發布有授權、附 checksum 的小型 fixture 之前，步驟 4–5 必須標為 **internal-data example**；讀者不能只靠 clone 直接重現該數值輸出。

---

## ⚠️ 開始之前：先確認機器狀態

2026-08-06 內部實測機有 48 核，當時 `/big7_disk` 只剩約 617 GB；這不是讀者環境的固定值。
跑分析前請在自己的輸出磁碟執行 `df -h` 並確認餘量，大型輸出很容易把磁碟寫爆。

---

## 步驟 1 · 編譯 C++（約 2–5 分鐘）

🔴 **一定要先做這步。** 不要依賴工作目錄裡既有 binary 的 mtime 狀態；每次 build/run
都應記錄 source commit、dirty diff、binary SHA-256、compiler 與 build command。

```bash
git clone https://github.com/liaoyoyo/InterSubMod.git
cd InterSubMod

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

**版本限定的稽核輸出**（tracked core `73afaeac`，2026-08-12；未來 commit 請以實際輸出為準）：

```text
[==========] Running 270 tests from 39 test suites.
[  PASSED  ] 270 tests.

$ echo $?
0
```

✅ 驗收的是**目前 checkout 的實際退出碼為 0，且輸出沒有 failed tests**；不要把快照的
270/39 硬編成未來版本的固定契約。如果有任何失敗，先不要往下走。

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

## 步驟 4 · 跑出第一個分析結果（internal-data example）

先用單一突變位點試跑，確認整條路徑通暢。下列三個輸入必須由使用者自行提供；public repo
目前沒有可直接代入的 fixture。Historical internal smoke 只能證明該次 input exit 0，不能提供
一般 runtime；要報時程必須另附 hardware、input locus、commit、date 與 repetitions。

```bash
# 填入自己合法取得的絕對路徑
ISM_TUMOR_BAM=/absolute/path/to/tumor.bam
ISM_REFERENCE=/absolute/path/to/reference.fa
ISM_INPUT_VCF=/absolute/path/to/somatic.vcf.gz
ISM_DEMO_DIR=/tmp/ism_demo
mkdir -p "$ISM_DEMO_DIR"

# 先驗證外部輸入；reference 需有 .fai，BAM 需可被 samtools 讀取
test -r "$ISM_TUMOR_BAM"
test -r "$ISM_REFERENCE" && test -r "${ISM_REFERENCE}.fai"
test -r "$ISM_INPUT_VCF"
samtools quickcheck "$ISM_TUMOR_BAM"

# 取輸入 VCF 的第一筆 biallelic SNV 作小型 smoke test
bcftools view -h "$ISM_INPUT_VCF" > "$ISM_DEMO_DIR/one_snv.vcf"
bcftools view -m2 -M2 -v snps -H "$ISM_INPUT_VCF" | head -n 1 >> "$ISM_DEMO_DIR/one_snv.vcf"
test "$(bcftools view -H "$ISM_DEMO_DIR/one_snv.vcf" | wc -l)" -eq 1

# 跑（只給三個必填參數）
./build/bin/inter_sub_mod \
  --tumor-bam  "$ISM_TUMOR_BAM" \
  --reference  "$ISM_REFERENCE" \
  --vcf        "$ISM_DEMO_DIR/one_snv.vcf" \
  --output-dir "$ISM_DEMO_DIR/out_min"
```

**Historical internal HCC1395 single-locus receipt 的實跑輸出**（只保留功能輸出，不報一般 runtime）：

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
find "$ISM_DEMO_DIR/out_min" -type f | head -20
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

## 步驟 6 · 驗證 historical 35,332-site pipeline 數字（推薦）

```bash
cd docs/methodology/_assets/20260627_subclone_4axis_teaching/scripts
python3 verify_pipeline_numbers.py
```

這支只重算 2026-06-27 教材所用的 **historical 35,332-site pipeline** 指標並比對；它不驗證 2026-07-24 exact-PS、LongLineage、儲存量、code count 或目前 test count。實跑輸出：

```text
sSNV 總數 35,332（TP 30,490 / FP 4,842）              ✓
共現連上 21,554 · 訊號不足 5,458 · 孤立 8,320        （加總 = 35,332 ✓）
```

對得上，只代表 tracked historical data 能重現上述 35,332-site 指標。exact-PS 應另查 `docs/handoff/20260801_exactPS_readAF_CNV_AI交接_01/denominator_registry.tsv`、`authority_manifest.json` 與 solver receipts；LongLineage 則須在指定 commit 分別執行其 gates／binary inventory。

---

## 07 · 全流程速查

![howto-six-steps](https://raw.githubusercontent.com/liaoyoyo/InterSubMod/develop/docs/images/howto-six-steps.png)

> **圖 1 · 六個步驟與各自的驗收點** —— 六個步驟依序為：編譯、跑測試、裝 Python
> 依賴、單點試跑、檢視輸出、驗證數字。圖內若出現舊測試數或固定秒數，視為 historical
> illustration；current acceptance 以本頁文字、實際 exit code 與當次 command output 為準。

### 每一步的驗收條件（不過就停下來，別往下走）

| 步驟 | 驗收條件 |
|---|---|
| 1 編譯 | `build/bin/` 下出現 6 個執行檔 |
| 2 測試 | 當前 checkout 的 test command **退出碼 0，且 0 failed tests**；不硬編固定 test count |
| 3 依賴 | `import pysam` 不報錯 |
| 4 試跑 | 印出 `Successful: 1 / Failed: 0` |
| 5 輸出 | region 目錄下有 methylation 與 distance |
| 6 驗證 | historical 35,332-site 指標加總打勾；不代表其他 claim family 已驗證 |

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
| 跨 run 比較時欄位對不上 | `significance_summary.csv` 欄數隨版本變動。commit `73afaea` 為 199 欄，含 `VerificationSchemaVersion=2` 與 `RegionStratificationSchemaVersion=1` 兩個 component-level 欄，但仍無 single whole-file layout version；一律**用欄名並檢查 schema 欄**。 |
| 工作站 HTML 開很久／瀏覽器很卡 | 單檔可達 188 MB。示範時先開最小的那個（約 14.6 MB）。 |
| LongLineage 說 KernelBlocked | 🔴 **這是預期行為不是故障**。正式入口刻意 fail-closed，見 [LongLineage 分冊](https://github.com/liaoyoyo/InterSubMod/wiki/LongLineage-Engine)。 |

---

## 本頁的驗證方式

- **測試套件**：最近一次公開稽核快照為 tracked core `73afaeac`、2026-08-12：
  `./build/bin/run_tests` 得 270 tests / 39 suites、退出碼 0；CTest 270/270。
  這不是 future commits 的固定期望值，每次應保留當次 command output。
- **單點試跑**：historical internal HCC1395 single-locus receipt 的退出碼為 0；因 provenance
  不足以支援一般效能 claim，本頁不報秒數。
- **環境數字**：48 核、磁碟餘量 617 GB、兩個 Python 版本，皆為本輪實測。
- **6 個執行檔**：實際 `ls build/bin/` 所得。

## ⚠️ 誠實標註

- 本頁的編譯指令**未在本輪從零重跑**（現有 build 目錄已存在）；該流程於 2026-08-05 曾以全新目錄驗證通過。
- **全基因組全量跑法未實跑** —— 屬長時間計算，且需先確認磁碟餘量。

---

**來源**：`InterSubMod/docs/explain/16_how-to-run.standalone.html` · 分支 `research/subclonal-reconstruction-202606` · 建立於 2026-08-06

[← 上一頁：分析與呈現層](https://github.com/liaoyoyo/InterSubMod/wiki/Analysis-and-Presentation) · [回系統全景](https://github.com/liaoyoyo/InterSubMod/wiki/System-Overview) · [回 Home](https://github.com/liaoyoyo/InterSubMod/wiki)
