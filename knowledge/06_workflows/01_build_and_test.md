---
id: ism-kb-06-workflows-build-and-test
name: "Build & Test Workflow"
description: "編譯 ISM C++ 專案的標準流程：cmake + make + ctest；三種測試（quick/data/full）的用途與時長。"
status: active
last_verified: 2026-04-22
content_nature: runtime-fact
doc_type: howto
verified_scope: "build commands against .claude/rules/cpp-build.md"
related_ids:
  - ism-kb-06-workflows-index
  - ism-kb-06-workflows-full-vcf-analysis
  - ism-kb-06-workflows-cpp-change-pdd
  - ism-kb-00-governance-hooks-and-automation
tags: [workflow, build, test, cmake, make, ctest]
canonical_paths: [06_workflows/01_build_and_test.md]
alias_paths: []
---

# Build & Test Workflow

- 一句結論：`cmake .. && make -j$(nproc)` 編譯；測試分 quick (<30s) / data (~1min) / full (~5min)
- 適用對象：修改 C++ 後、初次使用專案
- 可直接執行命令（驗證日期：2026-04-22）：
  ```bash
  cd /big7_disk/liaoyoyo2001/InterSubMod/build && make -j$(nproc)
  ```

---

## 編譯流程

### 標準 Release build
```bash
cd /big7_disk/liaoyoyo2001/InterSubMod
mkdir -p build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
make -j$(nproc)
```

**產物**：`build/bin/inter_sub_mod`

### Debug build
```bash
cd build
cmake .. -DCMAKE_BUILD_TYPE=Debug
make -j$(nproc)
```

### 重新編譯單一變更
```bash
cd build && make -j$(nproc)
# make 會 incremental build，只重編改動的檔案
```

---

## 依賴（已列於 CMakeLists.txt）

| 依賴 | 用途 | 典型取得方式 |
|------|------|-------------|
| HTSlib | BAM 檔案處理 | `apt install libhts-dev` 或 conda |
| OpenMP | 平行運算 | 編譯器自帶（gcc 8+） |
| Eigen3 | 線性代數 | 頭檔-only，repo 內已包含 |
| GoogleTest | 單元測試 | FetchContent（CMake 自動下載） |
| jemalloc 5.3.0 | 記憶體分配 | 預編譯進 lib/ |
| CLI11 | CLI 解析 | 頭檔-only，repo 內已包含 |

**標準**：C++17

---

## 測試指令（按時長分三級）

### 1. `/test-quick` — chr19 驗證（<30 秒）
```bash
./scripts/run_vcf_all_snv.sh --mode chr19-verification
```
**目的**：確認主流程未破壞
**檢查**：chr19 上 10-20 個 SNV 的 ISM 分析成功率 ≥95%
**時機**：每次 C++ 修改後

### 2. `/test-data` — 輕量數據測試（~1 分鐘）
**目的**：端到端小規模驗證

### 3. `/test-full` — 全流程測試（~5 分鐘）
```bash
./scripts/run_batch_vcf_analysis.sh
```
**目的**：全量 HCC1395 TP+FP batch，F1 benchmark
**檢查**：F1 與上次 canonical 差異 < 0.01

### 4. `ctest` — GoogleTest 單元測試
```bash
cd build && ctest --output-on-failure
```

**主要 test**：
| Test | 涵蓋模組 |
|------|---------|
| `test_distance_matrix` | 6 種距離度量 |
| `test_hierarchical_clustering` | HierarchicalClustering |
| `test_read_parser` | ReadParser (MM/ML 解析) |
| `test_significance_stats` | GlobalTest / PERMANOVA |
| `test_label_test` | LabelTest (HP/Allele) |
| `test_per_cpg_asm` | Per-CpG ASM |
| `test_structure_bootstrap` | Bootstrap 穩定性 |

**單獨執行某 test**：
```bash
./build/tests/test_distance_matrix
```

---

## Post-Change Checklist（修改 C++ 後必做）

```bash
# 1. Compile
cd build && make -j$(nproc)

# 2. Quick verification (<30 sec)
cd .. && ./scripts/run_vcf_all_snv.sh --mode chr19-verification

# 3. 單元測試
cd build && ctest --output-on-failure

# 4. Full batch (若改動可能影響 F1)
cd .. && ./scripts/run_batch_vcf_analysis.sh

# 5. 比對 baseline F1 差異 < 0.01
```

**Hook 強制**：PreToolUse hook 會在 `git commit` 前阻擋未編譯的 C++ 變更（`exit 2`）

---

## 代碼規範

- **C++17**
- Google style，120 char line width，4 space indent（`.clang-format`）
- 格式化：`clang-format -i <file>`
- **Code comments 用英文**

---

## 新增模組流程

1. Header → `include/core/`
2. Implementation → `src/core/`
3. 更新 `CMakeLists.txt`（若需要）
4. 單元測試 → `tests/`

---

## Git Workflow（簡版）

專案使用 **feature branch 流**：

| 分支 | 用途 |
|------|------|
| `main` | 穩定發布；由 `develop` 合併 |
| `develop` | 開發主線；所有 feature branch 在此合併 |
| `refactor/phase1-safety` 等 | 長期 refactor 分支 |
| `feature/*` | 單一功能臨時分支 |

**常用 slash commands**：
| Command | 動作 |
|---------|------|
| `/git-start` | 從 develop 開新 feature branch |
| `/git-commit` | 含 Co-Author 簽名的 commit |
| `/git-finish` | 合併回 develop + 清理 |
| `/build` | cmake + make |
| `/test-quick` / `/test-data` / `/test-full` | 見上方測試指令章節 |

**commit message 慣例**（與 PDD 6 steps 對齊）：
- `chore: snapshot baseline ...`（Step 1）
- `test: ... unit test coverage`（Step 3）
- `feat: / fix: / refactor:`（Step 4）
- `docs: ... validation result`（Step 5）
- `chore: evidence_ledger record ...`（Step 6）

詳見 [07_cpp_change_pdd.md](07_cpp_change_pdd.md#commit-類型完整定義)

---

## 常見錯誤

| 錯誤 | 解法 |
|------|------|
| `htslib not found` | `apt install libhts-dev` 或設 `HTSLIB_DIR` |
| `openmp: command not found` | 確認 `gcc --version ≥ 8`；加 `-DCMAKE_CXX_COMPILER=g++-9` |
| `undefined reference to jemalloc` | 確認 `lib/jemalloc.a` 存在 |
| `Pre-commit hook blocked: build out of date` | 跑 `cd build && make -j$(nproc)` |

---

## 相關

- CLI 參數：[../04_parameters/01_cli_arguments.md](../04_parameters/01_cli_arguments.md)
- 專案規範：[../../.claude/rules/cpp-build.md](../../.claude/rules/cpp-build.md)
- Full VCF 分析：[02_full_vcf_analysis.md](02_full_vcf_analysis.md)
