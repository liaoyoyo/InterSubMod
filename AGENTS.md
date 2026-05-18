# AGENTS.md — InterSubMod Repository Guidelines

> **此檔職責**：跨 agent 共用 governance — 任何 agent 框架（OpenHands / Claude Code / Codex / Aider / Cline）都會讀。
> **Claude Code 特定行為（確認矩陣 / Skills / Hooks）→ `.claude/CLAUDE.md`**。

---

## §1 語言與輸出規範

- **回應語言**: 繁體中文（zh-TW）所有回覆、思考、任務清單
- **程式碼註解**: 英文
- **.md 路徑前綴**: `InterSubMod/...`（hook 強制檢核）；更高層級用全域絕對路徑
- **圖片嵌入**:
  - Markdown: `![標題](相對路徑)` 必須真實顯示，禁止只列路徑
  - HTML: 嵌入時保留**原始比例**，用相對路徑
- **設計準則**: 12 條 canonical（Tufte / CRAP / Assertion-Evidence / WCAG）→ `/doc-standards` skill

---

## §2 與 big7 workspace root 銜接

- Root 規範：`/big7_disk/liaoyoyo2001/AGENTS.md`
- **本檔權威範圍**：InterSubMod repo 內程式、研究文件、腳本、KB 查閱、開發習慣
- **Root 規範權威範圍**：跨目錄分工、輸出落點

**任務涉及以下任一目錄 → 必須同時閱讀並遵守 root 規範**：
- `/big7_disk/liaoyoyo2001/big7_disk_output/`
- `/big7_disk/liaoyoyo2001/InterSubMod_big7_runbook/`
- `/big7_disk/liaoyoyo2001/Meeting/`

衝突原則：路徑更具體者優先（InterSubMod 內以本檔 / 跨目錄以 root）。

**Agent 上下文 5 入口分工**（跨 agent 都需要知道）：
| 入口 | 權威範圍 |
|------|--------|
| `InterSubMod/AGENTS.md`（本檔）| 跨 agent governance |
| `InterSubMod/.claude/CLAUDE.md` | Claude Code 特定（其他 agent 不讀）|
| `InterSubMod/docs/references/manual/20260424_AI啟動壓縮上下文與研究索引_01.md` | 研究壓縮上下文、重要數據 |
| `InterSubMod/docs/CURRENT_FOCUS.md` | live 主軸、阻塞 |
| `InterSubMod/research/autoresearch/research_direction.md` | AutoResearch 候選 queue（**僅候選**，不作觸發）|

---

## §3 五大研究目標與當前主軸

**G1-G5**（每小任務必須回答「服務哪個 Gx」）：
- **G1**: longphase 家族 + ISM 整合突破（願景主軸）
- **G2**: longphase-s + longphase-to 協同 > 單獨任一
- **G3**: ISM read-level epigenetic 給領域突破
- **G4**: 多樣本一致性 + reproducibility
- **G5**: 業界級貢獻（可被外部驗證）

**當前主要 phase**: **Phase 2 Normal Methylation Reference (方向 A+D)** — 服務 G2/G3

**歷史觀察（非定論，仍應深入研究）**：
- Phase 1A paired-pure F1 delta = +0.0112（marginal）— **仍需更完整研究驗證**
- TO 模式甲基化整合無顯著 F1 增益（早期觀察）— **不應作為放棄方向理由**；可能機制：訊號被遮蔽 / 樣本特性 / metric 口徑不對
- **任何「F1 / TO mode / variant filter」方向仍是合理研究範圍** — 不可僅憑舊觀察直接放棄；應重新設計實驗驗證

**成功標準**: 每小任務完整詳細觀察 → 整合驗證 → 業界級可驗證突破

---

## §4 啟動研究任務 5 問

接到研究任務必先回答（防踩已知坑）：
1. **Thread D**（read-level epigenetic）相關？
2. **Thread B 撤回範圍**內？（已 NO-GO 方向）
3. 資料 **KDE-corrected** 否？
4. 需要 **VCF caller AF**（非 merged AF）否？
5. 觸及 **長計算 / C++ / 檔案搬移 / NO-GO gate** 否？

---

## §5 高影響場景識別

下列場景**典型落在高影響區**，需高信心理由才能自主推進，否則預設暫停：

| 場景 | 典型影響 | 危險 | 正確做法 |
|------|------|------|------|
| 研究重點排序 | 高（>1h 投入）| 浪費數天 | 列候選方向 + 收益/風險 → 等用戶決定 |
| 假說選擇 | 高（影響結論）| 隱含假設導致偏誤 | 寫出前提 + 可能 confound |
| 統計方法選擇 | 中-高 | 不同方法相反結論 | 說明為何選 + 替代方案 |
| 「改進」/「優化」模糊指令 | 高 | 方向無數 | 要求用戶定義成功標準 |
| 多檔案 / 多步驟重構 | 中-高 | 影響範圍不明 | 列受影響檔案 + 預期改動 |

---

## §6 假設陳述與 Step → Verify

**實作前必做**:
- 列關鍵假設；不確定標 `⚠ 待確認`
- 多種解讀時列所有合理選項 + 傾向理由（**不可默默選擇沉默執行**）
- 有更簡單替代方案要主動提

**Step → Verify 格式**（多步驟任務必用）：
```
1. [步驟] → 驗證: [具體可觀察]
2. [步驟] → 驗證: [具體可觀察]
```

**範例**（任務：加入 Normal BAM read 過濾邏輯）：
```
1. 讀取現有 ReadParser     → 驗證: 能指出 normal read 進入的函數位置
2. 新增 normal BAM 參數     → 驗證: make -j$(nproc) 無 warning
3. 實作過濾條件             → 驗證: 單元測試覆蓋 normal/tumor/mixed 三情境
4. 全量測試                  → 驗證: F1 差異 < 0.01
```

- **強驗證（要求）**: 數值範圍 / 檔案輸出含路徑 / 命令退出碼 / 預期輸出片段
- **弱驗證（禁止）**: 「看起來正確」「double-check」「讓它能動」「確認結果合理」

Opus 4.7 模型特性備註 → `.claude/CLAUDE.md` §2。

---

## §7 變動頻率規則

| 層級 | 上限 | 維護機制 |
|------|------|----------|
| Foundation（§3 目標）| 季度 | 手動 |
| Governance（§1-§8）| 月度 | 手動 |
| Working State（CURRENT_FOCUS）| **週級** ⚠ | SessionStart hook + 手動 |
| Memory Concluded 區 | grep-only | `/memory-consolidation` skill |

`InterSubMod/docs/CURRENT_FOCUS.md` > 7 天未更新 → hook 主動提醒。

> **業界對照**: `docs/CURRENT_FOCUS.md` 概念等同 [Cline Memory Bank `activeContext.md`](https://docs.cline.bot/features/memory-bank)（live working state，最高變動層）。命名保留 `CURRENT_FOCUS.md` 以維持既有 20+ 引用，不重新命名。
>
> **KB 鏡像**: `InterSubMod/knowledge/10_research_status/01_current_focus_snapshot.md` 為 2 週週期 KB 快照（每 14 天 freshness hook 警告需 verify）。

---

## §8 KB 義務（P-14 強制查詢順序）

**觸發**: 討論外部工具（longphase / clairs-to / claude-code）F1/AUC、論文 claim、業界口徑（bcftools isec / hap.py） → **不可從本專案 report 推論外部行為**

**強制 3 步**:
1. 查 KB：`mcp__knowledge__knowledge_search` 或 `Read /big8_disk/Knowledge/05_tools/<tool>.md`
2. 引用 paper §N 或 KB 段落（明示出處）
3. 才下結論，對照本專案 report 口徑

**反例 ❌**: 「我記得 longphase 不改 FILTER」（憑記憶）
**正例 ✅**: 「先 Read `Knowledge/05_tools/longphase-to.md` 確認論文 §4.3 F1 口徑為 V_H/V_L post-filter」

### KB 主題速查表

| 主題 | 路徑 | 觸發關鍵字 |
|------|------|-----------|
| 資料總覽 | `/big8_disk/Knowledge/01_data_overview/` | 資料位置、目錄結構 |
| 癌症樣本 | `02_samples/` | HCC1395, COLO829, H1437, H2009, HG002 |
| 檔案格式 | `03_file_formats/` | VCF, BAM, MM/ML, phased VCF, HP tag |
| 資料庫 | `04_databases/` | PON, gnomAD, dbSNP, SEQC2 |
| 工具 | `05_tools/` | LongPhase, ClairS, ClairS-TO, DeepSomatic, InterSubMod |
| 流程 | `06_workflows/` | somatic calling, phasing, methylation |
| 腳本 | `07_scripts/` | auto_run.sh, benchmark |
| 論文 | `08_references/` + `paper/` | paper, 論文 |

---

## §9 Project Structure

- `src/` C++ core（`core/` / `io/` / `utils/`）+ `include/` headers
- `tests/` GoogleTest unit tests；`src/test/` phase-specific drivers
- `tools/` Python 分析/繪圖；`scripts/` shell workflows
- `data/` 範例輸入；`output/` → symlink `/big7_disk/liaoyoyo2001/big7_disk_output/`
- `docs/` 文件；`research/` 研究工作區（含 `figures/` `data/` `scripts/` `reports/`）
- `build/` CMake out-of-tree 輸出

---

## §10 Build / Test / Dev 命令

```bash
# Build
mkdir -p build && cd build && cmake .. && make -j$(nproc)
# Binary: build/bin/inter_sub_mod

# Run core
./build/bin/inter_sub_mod --tumor-bam data/tumor.bam --reference data/ref.fa \
    --vcf data/somatic.vcf --output-dir results

# Full pipeline
./scripts/run_vcf_all_snv.sh --mode all-with-w1000 --plot-type distance

# Test
ctest --test-dir build  # 或 ./build/bin/run_tests

# Python plotting deps
pip install -r requirements.txt
```

Coding: C++17 / `.hpp` headers / namespace `InterSubMod` / `.clang-format`（Google 4-space 120-col）/ CamelCase class + snake_case method/file

---

## §11 執行檔案 IO 顯示規則

執行檔案、程式或命令時**必須**顯示：
1. **輸入路徑**（含完整路徑前綴 `InterSubMod/...` 或全域路徑）
2. **執行命令**（完整命令含參數）
3. **輸出路徑**（產出檔案位置）
4. **執行後**顯示**實際輸出片段**（供用戶確認結果與預期一致；如未跑過則先述「預期應產出」）

目的：用戶可逐步檢核，避免錯誤輸入 / 不知道輸出位置 / 結果未驗證。

---

## §12 Output 目錄結構

| 目錄 | 分類 | 說明 |
|------|------|------|
| `output/canonical/` | Canonical | 7 樣本 × 3 模式 ISM baseline（19 runs）|
| `output/synthesis/research_rounds/` | Research rounds | 新研究 round 落點 |
| `output/synthesis/observation_workspaces/` | Observation | 跨樣本診斷、觀察工作區 |
| `output/big8_output_archive/` / `bip8_output_archive/` | Archive | 歷史，**非**新輸出落點 |

**完整索引**: `InterSubMod/output/OUTPUT_INDEX.md`
**信任度規範**: `InterSubMod/docs/data_specs/20260414_output資料信任度與生命週期_01.md`
**Region 結構細節** → `.claude/rules/output-structure.md`（條件載入）

---

## §13 AI Agent 預設操作政策

- **不可直接刪除檔案**（含 `rm`, `find -delete`, 覆寫式清空）— 需先搬移到 Archive 區
- 不可直接寫到 `/big8_disk/.../InterSubMod_runs/output` 等舊路徑
- 環境異常 / 結果不一致 → 跑 `check_ai_agent_readiness.sh`（非每次必跑）
- 清理動作必須在可審核腳本中執行

---

## §14 任務切割與 Agent 啟動

**原則**：可切割且需要大量 context 處理的任務 → 啟動 agent 處理

**必要紀錄**（科學工程化）:
1. 子任務啟動前：列清楚輸入 + 預期輸出 + 驗證標準
2. 子任務完成後：**清楚回報主 agent**（不只是「完成了」，要列具體結果）
3. **紀錄成文件**（`InterSubMod/research/{topic}/` 或 `docs/experiments/`）供檢核驗證查詢
4. 文件必含：執行命令、輸入路徑、輸出路徑、驗證結果、邏輯鏈

---

## §15 回應分級機制

依任務內容選擇對應回應方式：
- **預設**: Markdown 對話
- **大量解釋 / 複雜架構**: 整理成 HTML 報告（`/html-report-build`）讓用戶離線細看
- **需要用戶理解學習**: 啟用一步步對話 skill（如 `/fast-learning-coach`）
- **需要用戶判斷與推進**: 預設使用 **AskUserQuestionTool**

### §15.1 30 秒法則 — Tier 3+ 任務首句強制 TL;DR（2026-05-18 P2 audit M5 落地）

**規則**: 中等以上複雜度的任務回應，**第一句必須是 TL;DR**（≤30 秒可讀完），讓用戶 5 秒內掌握重點。

**Tier 分級**:
| Tier | 場景 | TL;DR 強制 |
|------|------|-----------|
| **T1** 簡單 | 純命令 / 單檔查詢 / 簡短事實 | 不需 — 直接答 |
| **T2** 中等 | 1 個 file 修改 / 簡單分析 / 短 report | 建議但不強制 |
| **T3** 複雜 | 多 file change / 跨 skill 協作 / decision table | **強制首句 TL;DR** |
| **T4** 重大 | Hard Gate / NO-GO / 跨樣本驗證 / paper-scope | **強制 TL;DR + 影響/信心標註** |

**TL;DR 格式**:
```
[動詞] X — [關鍵結果 / 決策]（[影響: 低/中/高], [信心: 低/中/高]）
```

**範例**:
- ❌ Bad: 「我先做了 A 然後做了 B 接著做了 C 最後...」（埋頭分析無重點）
- ✅ Good: 「**P1 fix 完成 — 18 skills 修補，D2/D3 非工具類 ❌ 歸零（影響: 中, 信心: 高）。** 細節如下...」

**例外**: 純對話延續（「OK 繼續」「然後呢」）不需 TL;DR。

**業界對照**: Tufte data-ink + Garr Reynolds Zen「先講結論再展開」+ McKinsey Pyramid Principle SCQA。

---

## §16 文檔規範（精要 + 跳轉）

**識別性規則**（內聯）:
- 檔案命名：`{YYYYMMDD}_{中文說明}_{NN}.md`
- 圖片命名：`{NN}_{英文描述}.png`
- 每 .md 開頭必含 HTML 註解元數據（建立時間 / 目標 / 處理範圍 / 關聯檔案）
- **圖片必用 `![標題](相對路徑)`**，禁止只列路徑（除非明示列路徑）
- 圖片相對路徑最深 2 層（例外：`output/`、`research/` 因深度差異可超過）

**🔴 Hard Gate — 封存原則（不可逆操作保護）**:
- **不刪除任何檔案** — 一律搬移到 `InterSubMod/docs/archive/YYYY/MM/`
- 封存時建立 `SUMMARY.md` 提取重點結論
- 原位置留 redirect notice（`ARCHIVED.md`）
- 大型檔案（>500 行）封存前提取精簡版

**資訊分層**（4 層）:
- Active（當月+進行中，原目錄）/ Recent（1-3 月，原目錄精簡）/ Archive（>3 月，docs/archive/）/ Deep Archive（歷史快照，immutable）

**執行性細節**: 完整命名例外、多步驟專案目錄、Git tracking 表 → `/doc-standards` skill

---

## §17 主要查詢路徑（4 層導航）

| 層級 | 檔案 | 何時查 |
|------|------|------|
| L1 入口 | `InterSubMod/docs/README.md` | 首次接觸 |
| L2 焦點 | `InterSubMod/docs/CURRENT_FOCUS.md` | 每次對話開始 |
| L3 歷史 | `InterSubMod/docs/experiments/INDEX.md` | 計劃新實驗 |
| L4 深度 | `InterSubMod/docs/reports/research_landscape/00_INDEX.md` | 完整理解 |

**Agent 查詢義務**:
- 開始研究前必讀 L2 + L3
- HP tag / LOH 相關 → L4 02 (Self-Phasing) + 04 (暫停判定)
- TO FP 過濾相關 → L4 01 (FP 全貌) + 03 (ISM 價值)
- 樣本/格式/工具 → KB（§8）
