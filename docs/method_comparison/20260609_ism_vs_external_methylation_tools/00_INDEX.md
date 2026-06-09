<!--
建立時間: 2026-06-10
狀態: in_progress (method-comparison study INDEX — Phase A 完成, Phase B 待執行)
報告類型: study_index
受眾: 任何讀者第一份 · PI · 外部協作者 · 未來的自己
framework: BLUF + 四層分層(INDEX 本身=L0+L1) + 導航表
data_sources: 01-06 各文件 + _assets/workflow_raw_result.json
provenance_note: 本 INDEX 為導航/摘要; 數字皆引自子文件(已各標來源)。
-->
<!-- provenance-verified: 摘要數字引自 01-06 子文件與 _assets; 子文件各自標一手來源。 -->

# 00 INDEX — ISM vs 外部甲基化方法比較研究

> **任何人先讀這份**。本研究回答：**網路上與 ISM 相關的甲基化軟體/論文有哪些、方法與我們的差異（含細節）、我們的方法可改進/學習什麼**。重點：ONT **modkit**、二代短讀癌症 DMR、甲基位點關聯方法、**Δβ 方法**、我們結果 × 外部論文數據交叉。

---

## L0 — 三句話結論

1. **方法定位**：外部甲基化方法分 **6 條軸**；ISM 主要站「**read-to-read 距離 + clustering（軸C）**」，而 modkit/DSS/methylKit 等主流站「**per-position 率差（軸A）**」——**兩者問的問題本質不同**（率差 vs 結構）。
2. **新穎性**：**無單一工具佔 ISM 完整組合**（read-read 距離矩陣 + **PERMANOVA 結構顯著性** + **normal-anchored somatic cis-test** + 5mC/5hmC 分軌 + LOH/CN 耦合），但**每個單獨元件都有更成熟對手**（cvlr 分群 / DSS 統計 / MHB within-read LD / qFDRP 同 kernel）。可防守新穎 = **這個組合 + supervised PERMANOVA + somatic cis 控制**。
3. **可改進**：最高 ROI = ① per-CpG Fisher → **beta-binomial + dispersion shrinkage**（DSS）；② 距離建在 **soft 機率**（cvlr/pycoMeth）；③ **停 5mC/5hmC max-collapse**。另發現**源碼引用錯誤**（DAMEfinder 被誤標「De Waele 2020」應為 Orjuela 2020）。

> ⚠ **誠實邊界**：「組合無人佔 ≠ 該組合更好/有用」—— 是否帶真訊號需 **Phase B 實機 benchmark** 驗證（06 已規劃，**尚未執行**）。所有 ASM 結果 = 單樣本 HCC1395 ⭐3 單 pipeline。

---

## L1 — 六文件導航（讀到夠就停）

| # | 文件 | 一句話 | 何時讀 |
|---|------|--------|--------|
| **01** | [`01_ism_method_spec_from_source.md`](01_ism_method_spec_from_source.md) | **我們的方法**：6 核心 + Δβ + cis-test + copy-partition，**源碼 file:line 確認** | 要懂 ISM 到底做什麼 |
| **02** | [`02_external_tools_survey.md`](02_external_tools_survey.md) | **外部地景**：82 工具依 6 軸分類，全 WebFetch 實證引用（modkit/DSS/MHB/cvlr…）| 要知道有哪些工具 |
| **03** | [`03_method_comparison_matrix.md`](03_method_comparison_matrix.md) | **對照矩陣**：ISM 能力 × 11 工具逐格 ✓/✗/部分 + point-by-point 細節差異 | 要答「細節哪裡不同」 |
| **04** | [`04_cross_analysis_data_vs_literature.md`](04_cross_analysis_data_vs_literature.md) | **結果×文獻**：我們 Δβ/ASM/結論 vs 外部論文數據（5 維 + alignment verdict）| 要答「我們數據 vs 外部觀察」 |
| **05** | [`05_improve_learn_recommendations.md`](05_improve_learn_recommendations.md) | **改進/學習**：10 條排序建議（FIX/LEARN/KEEP）+ 引用修正 | 要知道下一步改什麼 |
| **06** | [`06_benchmark_plan_phaseB.md`](06_benchmark_plan_phaseB.md) | **Phase B 方案**：實機下載運行規劃（modkit/cvlr/ASMS vs ISM）⚠未執行 | 要實機比較時 |
| HTML | [`report.standalone.html`](report.standalone.html) | 清楚整理的單頁互動版（同內容，PI 閱讀版）| 要好讀的視覺版 |
| 原始 | `_assets/workflow_raw_result.json` + `survey_digest.md` | 82 工具完整結構化資料 + 萃取 | 要查單一工具細節/引用 |

---

## L1 — 三個最像 ISM 的工具（投稿必正面對照）

| 工具 | 軸 | 與 ISM 的 delta | 狀態 |
|------|----|----------------|------|
| **cvlr** (Raineri 2023) | C | read 分群有，但**模型式 EM、無距離矩陣、無顯著性檢定、germline、固定 k** | ⚠ GitHub 404 |
| **ASMS** (Raineri 2024) | C | ONT no-phasing read-clustering ASM —— **單一最像**；ISM 多 PERMANOVA+normal-anchor | ⚠ UNVERIFIED 須親讀 |
| **qFDRP** (Scherer 2020) | C/E | **與 ISM NHD 同 normalized-Hamming kernel**，但塌成 scalar、無 label/PERMANOVA | ✓ |

---

## L1 — 五維交叉分析結論（我們結果 vs 外部論文）

| 維 | 結論 |
|----|------|
| Δβ 方法學 | 🟡 健全但缺 variance-stabilization，magnitude 需 M-value/beta-binomial 交叉 |
| ASM 量值/方向 | 🟡 與 Do2020「+5×/hypo」是**口徑差非矛盾**；Martin-Trujillo CN-confound 是最強對齊 |
| structure vs disorder | 🟢 **Landau 2014 明文排除 ASM 是 disorder 來源** = 最強背書 |
| locus-association 軸 | 🟢 三軸正交，ISM 組合無人佔（真空缺） |
| filter-NEG + phasing 白地 | 🟢 兩者皆被一手文獻支持（無 caller 用甲基當 filter；germline phaser 留 cancer/LOH future）|

---

## 狀態與下一步
- ✅ **Phase A 完成**（2026-06-10）：源碼確認 + 82 工具 survey + 對照矩陣 + 交叉分析 + 改進建議。產出 6 .md + HTML + 原始資料。
- ⏳ **Phase B 待用戶核准**：實機下載運行 benchmark（06 已規劃；涉及下載外部工具 + 長 compute，未自動啟動）。

## Provenance
- 外部工具事實：workflow **wf_9e64e51d-03d**（24 agents：10 survey + 8 deep-read + 5 cross + 1 synth-critic；全程 WebSearch→WebFetch 一手來源；~46 min，2.3M tokens）。原始：`_assets/workflow_raw_result.json`。
- ISM 方法：2026-06-09/10 直接讀 `src/core` + `include/core` + tsg `scripts/{03,34,37}` + `pipeline/lib/cis_asm_core.py`（file:line）。
- 內部結果數字：引自 `knowledge/11_external_literature/` validated 文件（各標來源）。
- 標 ✗UNVERIFIED 者為無法 fetch 確認；DAMEfinder/Sakamoto 引用修正見 05 #9/#10。
