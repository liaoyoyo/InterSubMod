# InterSubMod — 建置與自備資料執行指南

先讀 [2026-08-13 完整研究資料與軟體交接 snapshot](docs/handoff/20260813_完整研究資料與軟體交接_01/00_INDEX.md)：它是新人／AI、資料位置、軟體 I/O 與驗證 gate 的第一入口；科學 authority 日期仍為 2026-08-01，目前不得稱 release-ready。

本文檔先提供可由 public clone 重現的建置與測試，再說明如何以**使用者自備、已建立索引且帶 MM/ML tag** 的資料執行。

> [!IMPORTANT]
> 公開 repo 提供 **tiny synthetic DEMO fixture**，可驗證 clone→build→run→schema 接線，
> 但它不是生物資料、benchmark或science validation。`scripts/run_vcf_all_snv.sh`
> 仍含 `/big7`／`/big8` hard-coded paths，是實驗室 orchestration，不是 portable public quick start。

---

## 1. 環境準備與編譯 (Setup & Build)

### 系統需求

- **已驗證 OS**：Ubuntu 24.04 hosted runner；其他 Linux 組合須自行留下 receipt。
- **編譯器**：支援 C++17（hosted receipt：GCC 13.3）。
- **必要依賴**：CMake 3.14+、pkg-config、HTSlib、Boost headers、OpenMP；執行
  tiny fixture 另需 samtools，Python 驗收需 Python 3.10 與 `requirements-ci.lock`。

### 編譯步驟

```bash
# 1. Clone，並切到 handoff manifest 登錄的 immutable commit
git clone https://github.com/liaoyoyo/InterSubMod.git
cd InterSubMod
HANDOFF_COMMIT="<IMMUTABLE_HANDOFF_COMMIT_SHA>"
git checkout --detach "$HANDOFF_COMMIT"
test "$(git rev-parse HEAD)" = "$HANDOFF_COMMIT"
test -z "$(git status --porcelain)"

# 2. 在 repo 外建立 clean Release build
REPO_ROOT="$(pwd -P)"
BUILD_ROOT="$(mktemp -d "${TMPDIR:-/tmp}/ism-build.XXXXXXXX")"
cmake -S "$REPO_ROOT" -B "$BUILD_ROOT" -DCMAKE_BUILD_TYPE=Release
cmake --build "$BUILD_ROOT" -j$(nproc)
test -z "$(git -C "$REPO_ROOT" status --porcelain)"
```

編譯完成後，主程式位於 `$BUILD_ROOT/bin/inter_sub_mod`；build output 不進 Git。

---

## 2. 先驗證公開 clone 的建置

```bash
"$BUILD_ROOT/bin/run_tests"
ctest --test-dir "$BUILD_ROOT" --output-on-failure
```

test/suite count 必須由當下 commit 的 GoogleTest/CTest 輸出與release receipt動態產生；
驗收為exit 0、0 failure，不把歷史數字抄成未來版本的固定預期。

---

## 3. Public tiny synthetic E2E（DEMO）

```bash
# 建立 hash-locked Python 3.10 acceptance 環境
PYTHON_ENV="$(mktemp -d "${TMPDIR:-/tmp}/ism-python.XXXXXXXX")/venv"
python3.10 -m venv "$PYTHON_ENV"
"$PYTHON_ENV/bin/python" -m pip install --require-hashes \
  --requirement "$REPO_ROOT/requirements-ci.lock"
export PATH="$PYTHON_ENV/bin:$PATH"

# 會建立新的 /tmp work directory，不覆寫既有輸出
"$REPO_ROOT/scripts/handoff/run_tiny_public_e2e.sh" --repo-root "$REPO_ROOT" --jobs 4
```

`requirements.txt` 可作一般研究繪圖依賴或重建 lock 的輸入，但不屬 handoff acceptance 環境。

驗收片段應含 `TINY_E2E_RESULT`、`"all_pass": true` 與
`"tree_semantics": "read_dendrogram_from_methylation_distance_not_cellular_lineage"`。fixture source 位於
`tests/fixtures/tiny_public/`；BAM與index由builder在暫存目錄建立。此結果只證明
軟體接線與schema，scope固定為 **DEMO**，不得寫入science validation ledger。

---

## 4. 以自備資料執行 C++ 核心

```bash
SITE_PROFILE="<SITE_PROFILE_JSON>"
"$REPO_ROOT/scripts/site/doctor" --profile "$SITE_PROFILE" --mode real-preflight
eval "$("$REPO_ROOT/scripts/site/site_profile.py" shell --profile "$SITE_PROFILE" --sample HCC1395)"

"$BUILD_ROOT/bin/inter_sub_mod" \
    --tumor-bam "$TUMOR_BAM" \
    --reference "$REFERENCE" \
    --vcf "$SOMATIC_VCF" \
    --output-dir "${TMPDIR:-/tmp}/ism-real-smoke" \
    --threads 16 \
    --window-size 1000 \
    --distance-metric NHD
```

必要前置：

- BAM 有 `.bai` 並含 `MM`／`ML` methylation tags。
- FASTA 有 `.fai`。
- VCF 的 reference build 與 BAM／FASTA 相同。
- 無參數時 effective distance default 是 `NHD`；為了 provenance，仍建議明寫。

Exact-PS topology funnel 來自**獨立的 research solver pipeline**，不是 core executable 的 product。
`inter_sub_mod` 只寫出 per-region methylation／statistics 資料，也不會自動產生 PNG；兩條 provenance chain 不得混用。

Frozen release baseline `ddd8909a` 的 `significance_summary.csv` source header 是 **199 欄**；歷史
`73afaeac-dirty` audit 的 C++／CMake inputs 與它 byte-equivalent，但不是 release source。兩個 component-level
schema 欄不等於 whole-file layout version。下游必須用欄名並記錄 producer commit，不得套用歷史欄位位置。

---

## 5. 內部 orchestration（非 portable）

`scripts/run_vcf_all_snv.sh` 可在實驗室環境串接 C++ 與 Python 圖表，但目前依賴未納入 Git 的內部資料與絕對路徑。只有在那些依賴存在時才可執行：

```bash
# INTERNAL ONLY — 不是 fresh public clone acceptance command
./scripts/run_vcf_all_snv.sh --mode all-with-w1000 -o output/my_custom_run
```

它不能因 tiny synthetic fixture 存在就升格；要成為 public real-data workflow，仍需移除
`/big7`／`/big8` hard-coded paths，改成顯式 BAM／FASTA／VCF 參數與site profile。

---

## 6. 輸出結果檢視 (Output)

直接執行 core binary 時，若未指定 `--output-dir`，輸出根目錄是 `output/`。根層放 run/schema 摘要；region 結果依 `<VCF_STEM>/<chr>/<chr_pos>/<chr_start_end>/` 分層：

```text
output/
├── run_params.json
├── run_summary.json
├── significance_summary.csv
├── significance_statistics.txt
├── region_stratification_status.tsv
└── <VCF_STEM>/
    └── <chr>/
        └── <chr_pos>/
            └── <chr_start_end>/
                └── ... per-region CSV/TSV/Newick outputs
```

`full_execution_analysis.log`、日期命名目錄與 `filtered_snv_tp` 是特定內部 wrapper／VCF stem 的產物，不是 direct-core generic default。

> [!NOTE]
> 大型輸出應寫到 site profile 宣告的 local data plane；不要在 public clone 內建立絕對 symlink。
> 輸出內容不納入 Git，只在 registry／receipt 登錄 logical ID、位置、scope、producer 與 hash。

---

> [!NOTE]
> 直接執行 core 會產出數據與結構檔（JSON/TSV/CSV/Newick 等），不會自動生成教學網站的 HTML 視覺化。Python 只能在 validated data 上做呈現，不應靜默重算或改寫 science。
