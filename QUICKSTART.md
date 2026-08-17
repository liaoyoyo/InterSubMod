# InterSubMod — 建置與自備資料執行指南

本文檔先提供可由 public clone 重現的建置與測試，再說明如何以**使用者自備、已建立索引且帶 MM/ML tag** 的資料執行。

> [!IMPORTANT]
> 公開 repo 目前**沒有**可直接分析的 BAM／FASTA／VCF fixture；因此不能宣稱 fresh clone 後 10 分鐘內得到第一個生物分析結果。`scripts/run_vcf_all_snv.sh` 亦含 `/big7`／`/big8` 內部路徑，是實驗室 orchestration，不是 portable public quick start。

---

## 1. 環境準備與編譯 (Setup & Build)

### 系統需求
*   **OS**: Linux (Ubuntu 18.04+ Recommended)
*   **Compiler**: GCC 7+ 或 Clang 5+ (支援 C++17)
*   **Dependencies**: HTSlib 1.10+, CMake 3.14+, Python 3 (Matplotlib, Seaborn, Pandas)

### 編譯步驟

```bash
# 1. Clone 專案
git clone https://github.com/liaoyoyo/InterSubMod.git
cd InterSubMod

# 2. 建立並編譯 (使用最大並行數加速)
mkdir -p build && cd build
cmake ..
make -j$(nproc)
```

編譯完成後，主程式位於 `build/bin/inter_sub_mod`。

---

## 2. 先驗證公開 clone 的建置

```bash
# 在 repo 根目錄執行
./build/bin/run_tests
ctest --test-dir build --output-on-failure
```

2026-08-12 的 version-scoped audit 在 tracked core `73afaeac` 得到 GoogleTest **270 tests / 39 suites** 與 CTest **270/270**。未來 commit 的實際輸出才是當下 authority。

---

## 3. 以自備資料執行 C++ 核心

```bash
./build/bin/inter_sub_mod \
    --tumor-bam /path/to/tumor.mm_ml.bam \
    --reference /path/to/reference.fa \
    --vcf /path/to/candidates.vcf \
    --output-dir /path/to/results \
    --threads 16 \
    --window-size 1000 \
    --distance-metric NHD
```

必要前置：

- BAM 有 `.bai` 並含 `MM`／`ML` methylation tags。
- FASTA 有 `.fai`。
- VCF 的 reference build 與 BAM／FASTA 相同。
- 無參數時 effective distance default 是 `NHD`；為了 provenance，仍建議明寫。

`inter_sub_mod` 產生 per-region methylation／statistics 資料；它**不會**產生 exact-PS topology funnel，也不會自動產生 PNG。Exact-PS funnel 來自獨立 research solver pipeline。

---

## 4. 內部 orchestration（非 portable）

`scripts/run_vcf_all_snv.sh` 可在實驗室環境串接 C++ 與 Python 圖表，但目前依賴未納入 Git 的內部資料與絕對路徑。只有在那些依賴存在時才可執行：

```bash
# INTERNAL ONLY — 不是 fresh public clone acceptance command
./scripts/run_vcf_all_snv.sh --mode all-with-w1000 -o output/my_custom_run
```

要把它升格為 public workflow，需先移除 `/big7`／`/big8` hard-coded paths，改成顯式 BAM／FASTA／VCF 參數，並發布有授權的 tiny fixture、checksum 與預期輸出 receipt。

---

## 5. 輸出結果檢視 (Output)

執行完成後，請前往輸出目錄（預設為 `output/YYYYMMDD_vcf_*`）。
輸出包含根層摘要檔與依 VCF/染色體分層的結果目錄，示意如下：

```text
output/
├── full_execution_analysis.log  # 完整執行日誌
├── significance_summary.csv      # 各區域顯著性摘要
├── significance_statistics.txt   # 全域統計
└── filtered_snv_tp/
    ├── chr1/
    ├── chr2/
    └── ...
```

> [!NOTE]
> `output/` 建議在 repo 內保留為入口目錄，實際儲存可放在其他硬碟並以軟連結對接；`output` 內容不納入 Git 版本控管。

---

> [!NOTE]
> 手動執行僅產出數據檔案 (TSV/CSV)，不會自動生成熱圖。如需圖表，請接續執行 `tools/` 下的 Python 腳本。
