# InterSubMod (Under development)

InterSubMod - performs read-level integration of methylation profiles with haplotypes, somatic alleles, and tumor/normal labels, enabling somatic variant validation, subclone resolution, and single-molecule epigenomic clustering.

## 環境需求 (Prerequisites)

* **C++ Compiler**: 支援 C++17 (GCC 7+, Clang 5+, MSVC 2019+)
* **CMake**: >= 3.14
* **Libraries**:
  * htslib (需預先安裝, e.g., `sudo apt install libhts-dev`)
  * OpenMP (通常隨編譯器附帶)
  * Eigen3 (e.g., `sudo apt install libeigen3-dev`)
  * zlib

## 建置與編譯 (Build Instructions)

### 1. Clone 專案

```bash
git clone https://github.com/liaoyoyo/InterSubMod.git
cd InterSubMod
```

### 2. 建立並編譯 (Debug 模式)

Debug 模式包含偵錯符號 (`-g`) 並關閉優化 (`-O0`)，適合開發與測試。

```bash
mkdir -p build
cd build
cmake -DCMAKE_BUILD_TYPE=Debug ..
make -j$(nproc)
```

若需 Release 模式 (高效能)，請將 `-DCMAKE_BUILD_TYPE=Debug` 改為 `-DCMAKE_BUILD_TYPE=Release`。

## 執行與測試 (Usage & Testing)

### 執行主程式

編譯完成後，執行檔位於 `build/bin/inter_sub_mod`。

```bash
./bin/inter_sub_mod --help
```

範例指令 (需自行準備測試資料):

```bash
./bin/inter_sub_mod \
    --tumor-bam data/tumor.bam \
    --reference data/ref.fa \
    --vcf data/somatic.vcf \
    --output-dir results
```

### 執行單元測試

我們使用 GoogleTest 進行自動化測試。

```bash
./bin/run_tests
```

預期輸出應類似：

```text
[==========] Running 5 tests from 2 test suites.
[  PASSED  ] 5 tests.
```
