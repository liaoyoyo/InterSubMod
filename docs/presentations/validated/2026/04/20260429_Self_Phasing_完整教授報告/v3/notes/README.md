# v3/notes/ 目錄入口

> V3 是 v2 的 **12-slide 精簡演講版**（15-18 分鐘可講完）。本目錄只放 V3 自有的演講稿；所有深度資料、原始報告、業界對照、術語表、Q&A 預備都**繼續沿用 v2/notes/ 與 v2/source_materials/**，避免重複。

---

## V3 自有檔案

| 檔案 | 用途 | 適用時機 |
|------|------|---------|
| **`speaker_script_v3.md`** | 12-slide 演講稿（15-18 分鐘）；每 slide 含 開場/主敘事/過場/Q&A 預備；計時表與緊急加速指引 | 演講當天讀稿；演講前最後 30 分鐘 review |
| **`README.md`**（本檔） | v3/notes/ 入口；指向 v2 深度資料；V3 ↔ v2 對照速查 | 第一次接手 v3 時的導覽 |

---

## 深度資料（指向 v2/notes/）

V3 不重複任何 v2 已有的內容。下列 v2 文件**仍是 V3 的權威深度資料源**：

| v2 路徑 | 內容 | V3 何時需要 |
|---------|------|-----------|
| `v2/notes/speaker_script_v2.md` | 26-slide 完整講稿（30 分鐘）；嚴守 5 處口徑校準 | V3 跳過的細節需補完整脈絡時 |
| `v2/notes/key_metrics_table.md` | 所有重點數據 + 來源行號交叉驗證表 | 任何數字爭議；演講前 fact-check |
| `v2/notes/qa_11_questions.md` | Q1-Q13 教授預備問答（含 reviewer 補充 Q14） | Q&A 預備；現場深度提問即時查 |
| `v2/notes/glossary.md` | 30+ 術語表（工具、HP tag、PON、self-phasing、ISM 特徵、結構生物學、統計、內部命名、樣本、開放議題） | 教授不熟術語時備用 |
| `v2/notes/code_references.md` | 4 commit hash 與檔案行號對照（V2b / V3F / INDEL guard / V5） | 程式碼層細節提問 |
| `v2/notes/industry_references.md` | 業界對照（LongPhase 2022 / LongPhase-S 2025 / WhatsHap / DeepSomatic / MethPhaser） | S11 業界家族樹深度展開 |
| `v2/notes/00_background_context.md` | 6 部份背景說明：1 句總結 / 三系統切分 / HP tag 兩編碼 / PON 概念 / 三層 bug / V5 三層投票 / 4-commit / 量化驗證 / F1 因果鏈 / ISM 3-tier / LOH 兩層 / 跨樣本 / 五大目標 / R1-R9 caveat / F1-F8 行動 / Purity threshold | 全套背景脈絡，自包含閱讀 |

**v2 source_materials 三份原始報告**（深度數據最終源）：
- `v2/source_materials/01_V5_IGV_session_visual_audit.md` — IGV 視覺證據（SP1/2/3、V5max1/2/3）
- `v2/source_materials/02_Self_Phasing_Baseline_V5_Audit.md` — 一句結論審核版
- `v2/source_materials/03_longphase_TO_vs_V5_技術報告.md` — 13 段結構技術報告（V3 數據主要來源）

---

## 使用建議

| 場景 | 操作 |
|------|------|
| **首次接手 v3** | 先讀本 README → 再讀 `speaker_script_v3.md` → 必要時 fallback 到 `v2/notes/00_background_context.md` |
| **演講前 30 分鐘 review** | 只讀 `speaker_script_v3.md` 的「演講前準備」「計時表」「演講者準備檢核清單」 |
| **演講當天讀稿** | `speaker_script_v3.md` 主敘事段（每 slide）+ 過場句即可 |
| **Q&A 預備** | 先讀 `v2/notes/qa_11_questions.md` Q1-Q13；遇深度問題再 fallback 到 `v2/source_materials/` |
| **數據查證 / fact-check** | 所有數字以 `v2/notes/key_metrics_table.md` 為準（含來源行號） |
| **遇教授追問程式碼細節** | fallback 到 `v2/notes/code_references.md`（4 commit hash + 檔案行號） |
| **教授追問業界定位** | `v2/notes/industry_references.md` + `v2/notes/glossary.md` |

---

## V3 → v2 對照速查（12 slides ↔ 26 slides）

| V3 Slide | v2 Slides | 主訊息 | 處理 |
|---------|-----------|--------|------|
| **S1** 17.3:1 → 1:1 | v2 S1 | 衝擊 TL;DR + roadmap | 保留 |
| **S2** 工作流程一覽 | v2 S2+S3 | Pipeline 邊界 + InterSubMod 無 C++ 改動 | 合併 |
| **S3** HP tag 五值 + 三層證據 | v2 S3+S4 | HP tag 定義 + 全基因組 17.3:1 + 個別位點 | 合併 |
| **S4** Phasing OK，Tag 有問題 | v2 S5+S6+S7+S8 | Paired vs TO 對照 + LOH 兩層 + 拆鎖 | 合併 |
| **S5** 29/14/42 features | v2 S9+S10 | TO 根基 + ISM 3-tier 影響 | 合併 |
| **S6** ★ 根因樹 | v2 S11+S12+S13 | Purity 0.927 觸發 + 三 bug + PON 概念 | 保留（核心新發現） |
| **S7** ★ V5 三層投票 | v2 S14 | germline-first + Layer 1.5 + encode | 保留 |
| **S8** Sanity 15/15 + 5 證據鏈 | v2 S16 | 4 守恆律 + 多源驗證 | 保留 |
| **S9** 量化指標 | v2 S15+S17 | AMB / HP33 / concordance + 4-commit timeline | 合併 |
| **S10** 🎯 V5max1 climax | v2 S18+S19 | 39 reads 翻轉 + 3/3 SP-extreme + caveat | 合併 |
| **S11** 業界家族樹 + 五大目標 | v2 S20+S21+S22+S23 | LongPhase 世系 + 解鎖鏈 + 已有發現 | 合併 |
| **S12** Take-home + 3 P0 + Q&A | v2 S24+S25+S26 | 收束 + 行動承諾 | 合併 |

**移入 speaker note 不另開頁的 v2 內容**：
- v2 S13 程式碼 diff → 入 V3 S7 note + Q&A backup
- v2 S20 矛盾與盲點 → 入 V3 S11 note
- v2 S22 後續可動 → 入 V3 S11/S12 note

**結果**：26 → 12 slides（壓縮 54%），保留 4 個必講核心（S1 / S6 / S7 / S10），其餘合併或移入 speaker notes。

---

## 5 處口徑校準（V3 嚴守，與 v2 一致）

1. **V5 不修 self-phasing 本身**（V2b 已處理 phase scaffold）— S10 caveat 必說
2. 「**同實驗室相鄰工作**」非「業界共識」— S11 必說
3. **ISM 影響只列 count（29/14/42 共 85）**，來源報告誤寫 7%（實際 14/85=16.5%），slide 全部用 count — S5 必注意
4. `Util.h:20-25`（不是 21-25，HAPLOTYPE_UNDEFINED=-1 起於 line 20）— S6 Bug 2 必說
5. **+8.3 pp（clean PS blocks，全基因組）**，不是全基因組 raw，不是 cherry-picked — S1 + S9 必加 caveat

---

## V3 演講核心摘要（20 秒口袋版）

> 「Tumor-only 模式下 somatic read 在 LongPhase-TO 異常偏向 HP1，HCC1395 5kHz 全基因組 17.3:1。根因是 Purity 0.927 卡 0.95 閾值未觸發 Two-Pass，加上 getVote() priority bug 跟 enum vs int literal mismatch 兩層 tag 端 bug。V5 用三層投票（germline-first + somatic fallback confidence 0.6 + encode）修補，AMB% 17.5→8%，read-level concordance +8.3 pp（clean PS 全基因組）。InterSubMod 無 C++ 改動，修補在外部 fork longphase-to-mod。F1 不變是預期（V5 不改 caller），真實價值在解鎖 ISM 五大目標的可信度前提。」

---

## 文件版本

- **建立日期**：2026-04-30
- **對應 PPTX**：v3 12-slide 精簡版（產出於 `v3/output/`）
- **作者**：InterSubMod Research（演講稿由 AI agent 整合 v2 全套素材）
- **下次更新時機**：F1 V5 commit 完成後補 commit hash；F3 7 樣本擴展完成後補跨樣本 concordance
