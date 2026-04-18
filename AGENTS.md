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
- `data/` stores example inputs; `output/` is the repo-local symlink entry to `/big7_disk/liaoyoyo2001/big7_disk_output/`; `docs/` supports documentation; `research/` holds research workspaces with figures/data/scripts.
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

## 繼續研究前的必讀清單

見 `.claude/CLAUDE.md`「繼續研究前的必讀清單」區段（含速查表）。

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

**路徑**：`/big8_disk/liaoyoyo2001/knowledge/`（`Knowledge` 為 symlink，兩者等價）

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

---

## 研究輸出組織規範（2026-04-05）

### 目錄結構

| 目錄 | 用途 | 索引檔 |
|------|------|--------|
| `research/{study_name}/` | 單一研究主題的完整工作區 | `README.md` |
| `research/{study_name}/figures/` | 該研究的所有圖表（PNG/SVG） | — |
| `research/{study_name}/data/` | 中間數據（TSV/CSV） | — |
| `research/{study_name}/scripts/` | 分析腳本 | — |
| `research/{study_name}/reports/` | 研究報告 | — |
| `docs/reports/{topic}/` | 多檔案正式說明文件 | `00_INDEX.md` |
| `docs/reports/{topic}/figures/` | 說明文件圖表 | — |
| `docs/experiments/` | 實驗紀錄索引 | `INDEX.md` |

### 圖片存放規則

1. **統一命名**：純圖片子目錄一律叫 `figures/`（不用 `images/`、`plots/`）
   - **例外**：當子目錄包含混合類型資源（圖片 + JSON/TSV/PDF 等非圖片檔案）時，可使用 `assets/` 命名
2. **相對路徑**：.md 引用圖片必須用相對路徑 `figures/xxx.png`，禁止絕對路徑
3. **最大深度**：圖片相對路徑最多 2 層（`../figures/` 可，`../../../figures/` 禁止）
   - **例外**：引用 `output/` 或 `research/` 下的實驗圖片時，因目錄結構深度差異，允許超過 2 層（需確保目標檔案存在）
4. **命名格式**：`{NN}_{英文描述}.png`（如 `01_stability_overview.png`）

### 檔案命名格式

- Markdown：`{YYYYMMDD}_{中文說明目標}_{流水號}.md`
  - **例外**：`architecture/`（半永久結構文件）、`refactor_baseline/`（歷史基線快照）、`research_landscape/`（序號索引系列 `NN_名稱.md`）及特殊用途檔案（`INDEX_DETAIL_ARCHIVE.md`、`*_ARCHIVED.md`）不受此限
- 圖片：`{NN}_{英文描述}.png`
- 數據：`{YYYYMMDD}_{描述}.tsv`

### 多步驟研究專案目錄

多 Step 研究在 `plans/` 和 `architecture/` 下建專案子資料夾：
```
docs/plans/YYYY/MM/{YYYYMMDD}_{專案主題}/
  ├── 00_總覽與執行順序.md    # 索引
  ├── Step1_xxx.md            # Step 計劃
  └── Step2_xxx.md
docs/architecture/{YYYYMMDD}_{專案主題}/
  ├── 架構文件.md
  └── 資料追蹤表.md
```

### Git 追蹤規則

| 類型 | 追蹤 | 說明 |
|------|------|------|
| `docs/reports/*/figures/*.png` | 追蹤 | 正式說明文件的關鍵圖 |
| `docs/reports/*/*.md` | 追蹤 | 正式說明文件 |
| `research/*/figures/*.png` | gitignore | 研究中間圖表 |
| `research/*/data/` | gitignore | 中間數據 |
| `research/*/scripts/*.py` | 追蹤 | 可重現腳本 |
| `*.pdf`（根目錄） | gitignore | 簡報 PDF |

### 資訊分層與封存規則（2026-04-05）

| 層級 | 條件 | 位置 | 說明 |
|------|------|------|------|
| Active | 當月 + 進行中 | 原始目錄 | 完整內容 |
| Recent | 1-3 月內已完成 | 原始目錄 | 保留但可精簡 |
| Archive | >3 月 或 已被取代 | `docs/archive/YYYY/MM/` | 只搬移不刪除 |
| Deep Archive | 歷史快照/重複 | `docs/archive/deep/` | immutable |

**封存原則**：
1. **不刪除任何檔案**，一律搬移到 `docs/archive/YYYY/MM/`
2. 封存時建立 `SUMMARY.md` 提取重點結論
3. 原位置留 redirect notice（`ARCHIVED.md`）
4. 更新所有引用該檔案的索引連結
5. 大型檔案（>500 行）封存前提取精簡版保留在活躍目錄

### 元數據要求

每個 .md 檔案開頭必須有 HTML 註解元數據：

```markdown
<!--
建立時間: YYYY-MM-DD HH:MM
目標: [本檔案的目標或用途]
處理範圍: [涵蓋的工作範圍]
關聯檔案:
  - [相關檔案路徑 1]
  - [相關檔案路徑 2]
-->
```

---

## 主要查詢路徑與重點資訊（2026-04-05）

### 四層導航架構

| 層級 | 檔案 | 回答的問題 | 何時查閱 |
|------|------|-----------|---------|
| L1 入口 | `docs/README.md` | 文件在哪裡？怎麼導航？ | 首次接觸專案 |
| L2 焦點 | `docs/CURRENT_FOCUS.md` | 現在在做什麼？什麼阻塞？ | 每次對話開始 |
| L3 歷史 | `docs/experiments/INDEX.md` | 過去試過什麼？成功/失敗？ | 計劃新實驗前 |
| L4 深度 | `docs/reports/research_landscape/00_INDEX.md` | 完整研究推論鏈、證據與穩定性 | 需要完整理解 |

### 重點資訊速查表

| 我想知道... | 去哪裡找 |
|-------------|---------|
| TO FP 為什麼過濾不掉 | `docs/reports/research_landscape/01_TO_FP問題全貌.md` |
| Self-phasing 是什麼、影響多大 | `docs/reports/research_landscape/02_Self_Phasing根因.md` |
| ISM 哪些特徵可信、哪些不可信 | `docs/reports/research_landscape/03_ISM分析價值界定.md` |
| 哪些結論需要修正後重測 | `docs/reports/research_landscape/04_暫停判定與重評估.md` |
| 8 條證據鏈的完整推論 | `docs/reports/research_landscape/05_證據鏈總覽.md` |
| 14 個結論各自的穩定度評分 | `docs/reports/research_landscape/06_結論穩定性審查.md` |
| 當前研究策略與優先級 | `docs/CURRENT_FOCUS.md` |
| 過去實驗的成功/失敗記錄 | `docs/experiments/INDEX.md` |
| 樣本資訊、VCF 格式、工具參數 | Knowledge Base（見上方知識庫區段） |
| LOH 深度調查數據 | `research/loh_investigation/` |
| FP 來源追蹤分析 | `research/fp_provenance/` |

### Agent 查詢義務

- **開始研究前**：必讀 L2（CURRENT_FOCUS）+ L3（INDEX）
- **計劃新實驗前**：必讀 L3 確認未重複失敗方向
- **涉及 HP tag / LOH**：必讀 L4 的 02（Self-Phasing）和 04（暫停判定）
- **涉及 TO FP 過濾**：必讀 L4 的 01（FP 全貌）和 03（ISM 價值界定）
- **樣本/格式/工具問題**：必查 Knowledge Base
