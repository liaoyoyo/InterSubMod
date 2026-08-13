# How to Run 操作手冊
[← Home](https://github.com/liaoyoyo/InterSubMod/wiki) · [System Overview](https://github.com/liaoyoyo/InterSubMod/wiki/System-Overview) · [InterSubMod](https://github.com/liaoyoyo/InterSubMod/wiki/InterSubMod-Engine) · [LongLineage](https://github.com/liaoyoyo/InterSubMod/wiki/LongLineage-Engine) · [Upstream](https://github.com/liaoyoyo/InterSubMod/wiki/Upstream-and-Data) · [Analysis](https://github.com/liaoyoyo/InterSubMod/wiki/Analysis-and-Presentation) · [How to Run](https://github.com/liaoyoyo/InterSubMod/wiki/How-to-Run)

> **從零到跑出第一個結果**
> 編譯、測試與 tracked tiny synthetic DEMO 可由 public clone 執行；真正的 `inter_sub_mod` 分析另需使用者合法取得的 BAM、FASTA 與 VCF。
> 下方 HCC1395 輸出是 **2026-08-06 內部資料收據**，不是 public-clone 的 runtime 承諾。

> **公開可重現性邊界**：Git 物件已追蹤 `tests/fixtures/tiny_public/` 的 synthetic FASTA／VCF／SAM source 與 schema，runner 會在暫存目錄建立 BAM/index。Git 不含 HCC1395 tumor BAM、hg38 FASTA 或 real-data VCF；tiny PASS 只證明軟體接線，不是生物驗證。

---

## ⚠️ 開始之前：先確認機器狀態

執行時間與容量依硬體、mount 與輸入而異。跑分析前請執行 site doctor，並在自己的
`$DATA_ROOT`／輸出 mount 檢查可用容量、索引與 tool hash；本頁不承諾固定 runtime。

---

## 步驟 1 · 編譯 C++（runtime 依環境而異）

🔴 **一定要先做這步。** Git 不發布受版控的 build output。請從 release manifest 複製 immutable commit SHA，detached checkout 後在 repo 外的新目錄 clean build；不要沿用 clone 內的舊 `build/`。

```bash
git clone https://github.com/liaoyoyo/InterSubMod.git
cd InterSubMod
HANDOFF_COMMIT="<IMMUTABLE_HANDOFF_COMMIT_SHA>"
git checkout --detach "$HANDOFF_COMMIT"
test "$(git rev-parse HEAD)" = "$HANDOFF_COMMIT"
test -z "$(git status --porcelain)"

# repo 外的全新 Release build；build output 不進版本控制
REPO_ROOT="$(pwd -P)"
BUILD_ROOT="$(mktemp -d "${TMPDIR:-/tmp}/ism-build.XXXXXXXX")"
cmake -S "$REPO_ROOT" -B "$BUILD_ROOT" -DCMAKE_BUILD_TYPE=Release
cmake --build "$BUILD_ROOT" -j$(nproc)
test -z "$(git -C "$REPO_ROOT" status --porcelain)"
```

**執行檔數量與名稱是 commit-specific**，由該次 clean-build receipt 動態保存；本流程只固定檢查兩個必要入口：

```bash
find "$BUILD_ROOT/bin" -maxdepth 1 -type f -perm -u+x -printf '%f\n' | sort
test -x "$BUILD_ROOT/bin/inter_sub_mod"
test -x "$BUILD_ROOT/bin/run_tests"
```

---

## 步驟 2 · 確認編譯結果是好的（不宣稱固定 runtime）

```bash
"$BUILD_ROOT/bin/run_tests"
ctest --test-dir "$BUILD_ROOT" --output-on-failure
```

**本輪實跑的真實輸出**：

```text
[==========] <N> tests from <S> test suites ran.
[  PASSED  ] <N> tests.

$ echo $?
0
```

✅ **以當次輸出動態取得 N／S，且 N/N 通過、0 failure** — 如果這裡有任何失敗，**先不要往下走** —— 後面跑出來的數字都不可信。

---

## 步驟 3 · 準備 Python 環境

```bash
# 建立與 hosted acceptance 相同的 hash-locked Python 3.10 環境
PYTHON_ENV="$(mktemp -d "${TMPDIR:-/tmp}/ism-python.XXXXXXXX")/venv"
python3.10 -m venv "$PYTHON_ENV"
"$PYTHON_ENV/bin/python" -m pip install --require-hashes \
  --requirement "$REPO_ROOT/requirements-ci.lock"
export PATH="$PYTHON_ENV/bin:$PATH"
```

### 🔴 Python 版本是 hard requirement

Public acceptance workflow 與 hash-locked CI 使用 **Python 3.10**：

| 指令 | 版本 | 用途 |
|---|---|---|
| `python3.10` | 3.10 | site tooling、portable workflow、fixture validation 與嚴格區域切割 |

先跑 `python3.10 --version`；缺少 3.10 應 fail closed，不把 import/dataclass 錯誤誤判成研究程式缺陷。
一般研究繪圖可另依 `requirements.txt` 建環境，但那不是 handoff acceptance receipt 的環境。

---

## 步驟 4 · Public tiny synthetic E2E（DEMO）

public repo 追蹤 synthetic FASTA／VCF／SAM source；runner在新的暫存目錄建立BAM/index，
完成 build→run→schema validation：

```bash
"$REPO_ROOT/scripts/handoff/run_tiny_public_e2e.sh" --repo-root "$REPO_ROOT" --jobs 4
```

驗收輸出含 `TINY_E2E_RESULT`、`"all_pass": true` 與
`"tree_semantics": "read_dendrogram_from_methylation_distance_not_cellular_lineage"`。這是醒目 **DEMO**：
只驗軟體接線，不是生物資料、benchmark或science validation，不寫入science ledger。

### 進階：自備／內部真實資料（internal-data example）

下列三個輸入必須由使用者自行提供並確認授權。執行時間取決於硬體與資料；
HCC1395 內部具名 fixture 只保存結果 receipt，不提供一般 runtime promise。

```bash
# 將 placeholder 換成本機 profile；絕對路徑只維護於 untracked machine profile／registry
SITE_PROFILE="<SITE_PROFILE>"
"$REPO_ROOT/scripts/site/doctor" --profile "$SITE_PROFILE" --mode real-preflight
eval "$("$REPO_ROOT/scripts/site/site_profile.py" shell \
  --profile "$SITE_PROFILE" --sample HCC1395)"
ISM_DEMO_DIR="${TMPDIR:-/tmp}/ism_demo"
mkdir -p "$ISM_DEMO_DIR"

# shell contract 由 profile 提供 TUMOR_BAM、REFERENCE、SOMATIC_VCF
test -r "$TUMOR_BAM"
test -r "$REFERENCE" && test -r "${REFERENCE}.fai"
test -r "$SOMATIC_VCF"
samtools quickcheck "$TUMOR_BAM"

# 取輸入 VCF 的第一筆 biallelic SNV 作小型 smoke test
bcftools view -h "$SOMATIC_VCF" > "$ISM_DEMO_DIR/one_snv.vcf"
bcftools view -m2 -M2 -v snps -H "$SOMATIC_VCF" | head -n 1 >> "$ISM_DEMO_DIR/one_snv.vcf"
test "$(bcftools view -H "$ISM_DEMO_DIR/one_snv.vcf" | wc -l)" -eq 1

# 跑（只給三個必填參數）
"$BUILD_ROOT/bin/inter_sub_mod" \
  --tumor-bam  "$TUMOR_BAM" \
  --reference  "$REFERENCE" \
  --vcf        "$ISM_DEMO_DIR/one_snv.vcf" \
  --output-dir "$ISM_DEMO_DIR/out_min"
```

**內部 HCC1395 收據的實跑輸出**（不是任意輸入的 runtime／數值預期）：

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

![howto-six-steps](https://raw.githubusercontent.com/liaoyoyo/InterSubMod/ddd8909a838318d8a77969313e9561c8ff9d01c2/docs/images/howto-six-steps.png)

> **圖 1 · 六個步驟與各自的驗收點** —— 六個步驟依序為：編譯、跑測試、裝 Python 依賴、tiny synthetic DEMO、檢視輸出、驗證歷史數字；測試數須由該 commit 的 CTest/run_tests 動態輸出，驗收為 0 failure。Tiny fixture 是 public fresh-clone smoke，不是 real-data science。

### 每一步的驗收條件（不過就停下來，別往下走）

| 步驟 | 驗收條件 |
|---|---|
| 1 編譯 | immutable commit clean checkout；repo 外 build 完成；receipt 動態保存 executable inventory |
| 2 測試 | **動態取得 test/suite count；0 failure、退出碼 0** |
| 3 依賴 | `import pysam` 不報錯 |
| 4 試跑 | `TINY_E2E_RESULT` 含 `"all_pass": true`，且 receipt 標 `scope=DEMO` |
| 5 輸出 | region 目錄下有 methylation 與 distance |
| 6 驗證 | historical 35,332-site 指標加總打勾；不代表其他 claim family 已驗證 |

🔴 任何一步的退出碼不是 0，或數字對不上 —— 先解決再繼續，否則後面的結果都不可信。

---

## 08 · 常見狀況排除

| 症狀 | 原因與解法 |
|---|---|
| 程式直接停，說找不到參考基因組 | 參考 FASTA **缺 `.fai` 索引是硬性錯誤**。跑 `samtools faidx` 建索引。 |
| Python 腳本 import 就崩，錯誤指向 dataclass | 先確認 **Python ≥ 3.10**；低版本不是受支援環境。 |
| CI 判定 `--help` 失敗 | `--help` 的**退出碼是 1 不是 0**（與參數錯誤共用同一路徑）。改寫檢查方式。 |
| 某支腳本跑了但什麼都沒發生 | 可能是**函式庫不是入口**（如 `ism_heatmap_std.py`）。看它有沒有 argparse。 |
| 讀 CSV 後整列變成一欄 | `linkage_matrix.csv` **實際是 tab 分隔**。加 `sep='\t'`。 |
| 跨 run 比較時欄位對不上 | `significance_summary.csv` 欄數隨版本變動。Frozen release baseline `ddd8909a` 的 source header 為 199 欄，含 `VerificationSchemaVersion=2` 與 `RegionStratificationSchemaVersion=1` 兩個 component-level 欄，但仍無 single whole-file layout version；一律**用欄名並檢查 schema 欄**。歷史 `73afaeac-dirty` audit 不是 release source。 |
| 工作站 HTML 開很久／瀏覽器很卡 | 單檔可達 188 MB。示範時先開最小的那個（約 14.6 MB）。 |
| LongLineage 說 KernelBlocked | 🔴 **這是預期行為不是故障**。正式入口刻意 fail-closed，見 [LongLineage 分冊](https://github.com/liaoyoyo/InterSubMod/wiki/LongLineage-Engine)。 |

---

## 本頁的驗證方式

- **測試套件**：執行 `"$BUILD_ROOT/bin/run_tests"` 與 `ctest --test-dir "$BUILD_ROOT" --output-on-failure`；test/suite count 由本次輸出與 release receipt 動態產生，不在文件手抄固定數字。驗收為退出碼 0、0 failure。
- **Public synthetic smoke**：由 fresh clone 執行 tracked fixture，保存 exit code、schema 與 checksum receipt；scope 固定為 DEMO。
- **內部單點 receipt**：具名 HCC1395 inputs 的歷史輸出與退出碼 0；不作 runtime 或 public reproduction claim。
- **環境條件**：由 site doctor 在每次驗收動態記錄，不把舊機器 CPU、容量或 Python 狀態抄成讀者預期。
- **執行檔 inventory**：數量／名稱依 immutable commit 決定，由 clean-build receipt 動態列出；不把歷史「6 個」當成 current contract。

## ⚠️ 誠實標註

- 2026-08-05 的全新目錄 build 是**具名歷史 receipt**；handoff acceptance 仍須對 release manifest 指定 commit 另做 repo 外 clean build，舊 build 不可代替。
- **全基因組全量跑法未實跑** —— 屬長時間計算，且需先確認磁碟餘量。

---

**來源**：`InterSubMod/docs/explain/16_how-to-run.standalone.html` · 分支 `research/subclonal-reconstruction-202606` · 建立於 2026-08-06

[← 上一頁：分析與呈現層](https://github.com/liaoyoyo/InterSubMod/wiki/Analysis-and-Presentation) · [回系統全景](https://github.com/liaoyoyo/InterSubMod/wiki/System-Overview) · [回 Home](https://github.com/liaoyoyo/InterSubMod/wiki)
