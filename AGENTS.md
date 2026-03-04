# Repository Guidelines

## 語言
- **語言要求**：所有回覆、思考過程及任務清單必須使用**繁體中文**(zh-TW)
- **固定指令**： 'Implementation Plan, Task List and Thought in **Traditional Chinese**

## Project Structure & Module Organization
- `src/` holds the C++ core (split into `core/`, `io/`, `utils/`); headers live in `include/`.
- `tests/` contains GoogleTest unit tests; `src/test/` contains phase-specific test drivers.
- `tools/` houses Python analysis/plotting utilities; `scripts/` contains shell workflows.
- `data/` stores example inputs; `output/` is generated analysis output; `docs/` and `images/` support documentation.
- `build/` is the out-of-tree build output created by CMake.

## Build, Test, and Development Commands
- Build: `mkdir -p build && cd build && cmake .. && make -j$(nproc)`; binary at `build/bin/inter_sub_mod`.
- Run core manually: `./build/bin/inter_sub_mod --tumor-bam data/tumor.bam --reference data/ref.fa --vcf data/somatic.vcf --output-dir results`.
- Full pipeline script: `./scripts/run_vcf_all_snv.sh --mode all-with-w1000 --plot-type distance` (see `--help` for options).
- Output checks: `./scripts/verify_output.sh` validates expected files and matrix dimensions.
- Python deps for plotting: `pip install -r requirements.txt`.
- Optional container: `docker build -f Dockerfile.dev -t intersubmod:dev .` and `docker run -it --rm -v $(pwd):/workspace intersubmod:dev`.

## Coding Style & Naming Conventions
- C++17 code with `.hpp` headers and `.cpp` sources; namespace is `InterSubMod`.
- Formatting follows `.clang-format` (Google base, 4-space indent, 120 column limit); run `clang-format` on touched C++ files.
- Naming patterns: `CamelCase` classes (e.g., `BamReader`), `snake_case` methods and files.

## Testing Guidelines
- Unit tests live in `tests/test_*.cpp`; run with `ctest --test-dir build` or `./build/bin/run_tests`.
- Phase tests compile to `build/bin/test_phase*` from `src/test/`; `scripts/run_random_snv_test.sh` provides a quick smoke test.
- No explicit coverage target is enforced; add GTest coverage for new core logic when feasible.

## Commit & Pull Request Guidelines
- Commit messages in recent history use short imperative summaries like `Add ...` or `Refactor: ...`; keep to one line and add a prefix (`Fix:`, `Docs:`) when helpful.
- PRs should include a concise summary, commands run, and sample outputs/logs or plots when analysis or visualization changes.

## Data, Outputs, and Configuration Tips
- Many scripts default to absolute `/big8_disk/...` paths; prefer overriding via flags like `--vcf`, `--out`, `--threads`, and `--plot-type`.
- Keep generated artifacts in `output/` and avoid committing large datasets unless explicitly requested.

## 繼續研究前的必讀清單（每次對話開始時強制執行）

**每次開始研究/分析任務前，必須依序閱讀以下文件，不得省略：**

1. **`docs/CURRENT_FOCUS.md`** — 當前進行中的事項、阻塞點與風險
2. **`docs/experiments/INDEX.md`** — 過去所有研究方向的成功/失敗結論與建議後續
3. **`docs/README.md`** — 如需了解文件導航與查閱路徑

**目的**：
- 避免重複已失敗的方向
- 對齊當前最優先目標
- 了解哪些結論已驗證、哪些尚未解決

**觸發條件**：開始任何研究分析、實驗設計、程式改進、或延續前次工作時，此步驟為必要前置。

---

## AI Agent 預設操作政策（2026-03-01）
- `check_ai_agent_readiness.sh` 採「異常觸發」：僅在環境重建、路徑變更、腳本異常、或結果不一致時執行，不要求每次任務都先跑。
- `output/` 保持 repo 內入口；實體輸出放在 repo 外硬碟（以軟連結對接）屬預設建議策略。
- Agent 不可直接刪除檔案（包含 `rm`, `find -delete`, 覆寫式清空）。
- 若需移除內容，先搬移到 Archive 暫存區：`/big8_disk/liaoyoyo2001/InterSubMod_runs/Archive_pending_delete/`，並回報清單，待使用者手動最終刪除。
- 除非使用者明確要求，否則不做任何實際清除動作；若清理行為必須存在，需寫在可審核的執行腳本中。


## 實驗室知識庫 (Knowledge Base)

**路徑**：`/big8_disk/liaoyoyo2001/knowledge/`

當對話涉及以下主題時，**必須**先查閱知識庫對應文件確認細節，再進行回答或操作：

| 主題           | 查閱路徑                    | 觸發關鍵字                                                              |
| -------------- | --------------------------- | ----------------------------------------------------------------------- |
| 資料總覽與路徑 | `01_data_overview/`         | 資料位置、目錄結構、儲存空間                                            |
| 癌症樣本資訊   | `02_samples/`               | HCC1395, COLO829, H1437, H2009, HG002, purity, subsample                |
| 檔案格式規格   | `03_file_formats/`          | VCF, BAM, MM/ML, FILTER, phased VCF, modcall, HP tag                    |
| 資料庫與參考集 | `04_databases/`             | PON, gnomAD, dbSNP, CoLoRSdb, SEQC2, truth set, reference genome        |
| 工具使用與參數 | `05_tools/`                 | LongPhase, ClairS, ClairS-TO, DeepSomatic, InterSubMod                  |
| 分析流程       | `06_workflows/`             | somatic calling, phasing, haplotagging, methylation analysis, benchmark |
| 腳本操作說明   | `07_scripts/`               | auto_run.sh, benchmark script, 自動化腳本                               |
| 論文與參考資料 | `08_references/` + `paper/` | paper, 論文, server paths                                               |

### 查閱深度指引

| 情境                   | 查閱深度 | 動作                              |
| ---------------------- | -------- | --------------------------------- |
| 快速確認（路徑、名稱） | 淺層     | 讀 `README.md` 速查表             |
| 格式或參數細節         | 中層     | 讀對應子目錄的特定文件            |
| 完整流程或工具操作     | 深層     | 讀 workflow + tool 文件，交叉驗證 |
| 工具原始碼邏輯         | 最深層   | 讀 `codebase/` 目錄下的原始碼     |

### 查閱原則

- **不要憑記憶回答可以查證的事實**：檔案路徑、工具參數、VCF 欄位定義等務必查閱確認
- **引用來源**：回答時標註「根據 Knowledge/03_file_formats/vcf_clairs_to.md」
- **發現過時資訊時主動提醒使用者**
