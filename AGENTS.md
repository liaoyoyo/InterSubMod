# Repository Guidelines

## 語言
- **語言要求**：所有回覆、思考過程及任務清單必須使用**繁體中文**(zh-TW)
- **固定指令**： 'Implementation Plan, Task List and Thought in **Traditional Chinese**

## 與 big7 workspace root 規範的銜接
- 本 repo 之外，工作區另有一份 root 級規範：`/big7_disk/liaoyoyo2001/AGENTS.md`
- 分工原則：
  - `InterSubMod/AGENTS.md`：管理 `InterSubMod` repo 內的程式碼、研究文件、腳本、知識庫查閱與開發習慣
  - `/big7_disk/liaoyoyo2001/AGENTS.md`：管理整個 big7 工作區的目錄分工、輸出落點、runbook/meeting/canonical/synthesis 的角色與新建檔案規範
- 若任務只修改 `InterSubMod/` 內的程式碼、文件或腳本，先遵守本檔，再視需要參考 root 規範。
- 若任務涉及以下任一目錄，必須同時閱讀並遵守 root 規範：
  - `/big7_disk/liaoyoyo2001/big7_disk_output/`
  - `/big7_disk/liaoyoyo2001/InterSubMod_big7_runbook/`
  - `/big7_disk/liaoyoyo2001/Meeting/`
- 判斷原則：
  - 先用 root 規範判斷「檔案應該放在哪一類目錄」
  - 再用本檔判斷「若檔案放在 `InterSubMod/` 內，應遵守哪些 repo 內規則」
- 若兩份規範看似衝突，優先採用「路徑更具體者」：
  - `InterSubMod/` 內部細節以本檔為準
  - 跨目錄放置與 big7 輸出分流以 root 檔為準

## Project Structure & Module Organization
- `src/` holds the C++ core (split into `core/`, `io/`, `utils/`); headers live in `include/`.
- `tests/` contains GoogleTest unit tests; `src/test/` contains phase-specific test drivers.
- `tools/` houses Python analysis/plotting utilities; `scripts/` contains shell workflows.
- `data/` stores example inputs; `output/` is the repo-local symlink entry to `/big7_disk/liaoyoyo2001/big7_disk_output/`; `docs/` and `images/` support documentation.
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
- `output/` is the repo entry point to the current big7 output root: `/big7_disk/liaoyoyo2001/big7_disk_output/`.
- New formal sample/mode/run outputs belong under `output/canonical/{sample}/{canonical_mode}/{run_id}/`.
- New pure / tumor-only / exploratory rounds belong under `output/synthesis/research_rounds/`; cross-sample diagnostics and observation workspaces belong under `output/synthesis/observation_workspaces/`.
- `output/big8_output_archive/` and `output/bip8_output_archive/` are historical archive roots, not the default destination for new outputs.
- Some older scripts or docs still mention absolute `/big8_disk/...` output paths; treat those as legacy or archive-era references unless the task explicitly targets historical runs.
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
- `output/` 保持 repo 內入口；目前實體輸出根目錄固定為 `/big7_disk/liaoyoyo2001/big7_disk_output/`。
- 正式 sample/mode/run bundle 放在 `output/canonical/`；研究 round、觀察工作區與彙整產物放在 `output/synthesis/`。
- 除非任務明確指定歷史 archive，否則不要再把新的正式輸出寫到舊的 `/big8_disk/.../InterSubMod_runs/output` 類路徑。
- Agent 不可直接刪除檔案（包含 `rm`, `find -delete`, 覆寫式清空）。
- 若需移除內容，先搬移到與目前 `output/` 同卷、由 root 規範或使用者指定的 Archive 暫存區；若暫存區尚未建立，先回報並確認，不可直接刪除。
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
