<!--
建立時間: 2026-04-23
用途: PPT 資料夾自包含說明文件 — 任一 AI / 協作者從此入口即可理解全部背景
受眾: 未來需修改本 PPT 的 AI 或人類（新 session / 新成員）
-->

# 研究週報 PPT 資料夾 · 設計脈絡總入口（v2）

**主題**：LOH × AF × CN 下的 TP/FP 切分 — NG=2 LOH-constrained phasing 跨樣本確認 · CN KDE 校準 · PON 解 TO Self-Phasing
**涵蓋區間**：2026-04-16 ~ 2026-04-23（8 天，與 0421 週報重疊 6 天）
**受眾**：PI（熟悉 ONT、cancer genomics、somatic calling）
**總 slide 數**：26 頁（v2 採 7 段線性主軸，取代 v1 三幕偵探故事）

## v2 重要變更（2026-04-23 下午）

基於用戶 feedback「資料還是太雜亂 · 清楚的整理幾個報告主軸」，本 PPT 從「三幕偵探故事」改為**線性主題 7 段**：

| 段 | 涵蓋 Slides | 主題 |
|:-:|:---:|---|
| 段 0 | S1-S2 | 開場 + roadmap |
| 段 1 | S3-S6 | **背景觀察**（LOH.bed 佔比 / AF 分佈 / LOH×AF×CN 初步 / HCC1395 背景）|
| 段 2 | S7-S9 | **CN 用 KDE 重新定義**（舊 75× / KDE 驗收 / 各樣本新 CN）|
| 段 3 | S10-S13 | **TP rate 高的現象**（Phase 3 跨樣本 + Inner NG=2 + 生物學 + 統計）|
| 段 4 | S14-S16 | **LOH × AF × CN 切分框架**（S1-S7 scheme / 跨樣本熱圖 / S3+LOH_Noise）|
| 段 5 | S17-S18 | **新特徵有效性**（S4 secondary B4 / 下一步系統性計劃）|
| 段 6 | S19-S23 | **TO Self-Phasing 問題與解法**（HP tag 需要 / 根因 / 17.3:1 bias / PON 錨點）|
| 段 7 | S24-S25 | **未來工作**（P0/P1 下週 / 中長期目標）|
| 結語 | S26 | 小結 |

**同時整合 2026-04-23 下午新產 B1-B7 + Phase 3 Synthesis 結論**：
- Phase 3 Synthesis 跨 6 TO · S1 median 0.876 · S5 median 0.876 · 6/6 POSITIVE
- B1 Wilcoxon p=0.0156, W=21, median gap=0.365
- B2 HCC1954 outlier 源自 caller FP 背景 84%（非 phasing failure）
- **B3 Paired 對照 gap=0.00003, p=0.578 → TO-specific 完整驗證**
- B4 S4 內 AF+AlleleDelta LR AUC 0.717
- B5 COLO829 S1 fold 低為 small-n artifact
- B7 LOH_Noise 保留為獨立統計單元

**結論狀態更新**：Thread D 從「⭐5 TO-only（H-D4 懸掛）」升為「⭐5 TO-specific 確定性（H-D1/D2/D3/D4 全確認）」。

---

## 本資料夾內容

| 檔案 | 用途 | 優先閱讀順序 |
|------|------|:-----------:|
| `README.md` | **本文件** — 整個 PPT 的設計脈絡總入口 | **1** |
| `source_weekly_report.md` | **完整週報（1046 行、30 張圖）** — PPT 所有內容的權威來源 | **2** |
| `00_director_storyboard.md` | **分鏡稿** — 每頁 slide 的視覺/口頭/捨棄三層分配 | **3** |
| `build_pptx.py` | PPTX 生成腳本（python-pptx 1.0.2） | 修改時讀 |
| `20260423_NG2_LOH_constrained_phasing_pivot.pptx` | 輸出 PPTX（26 slides） | 最終產物 |
| `screenshots/slide-NN.png` | 每 slide 渲染 PNG（LibreOffice → PDF → pdftoppm）| 驗證對照 |

---

## 研究故事線（10 秒版本）

```
上週結尾 (0421)
    ↓  「LOH × AF × CN 切分 TP/FP 有雙極分佈」
本週進入
    ↓  想加入 HP / ISM 甲基 / ISM 聚類特徵進一步切分
Act I  · CN 方法學升級 (75× hardcoded → KDE per-sample 校準)
    ↓  S4 ambiguous bucket (75% TP + 76% FP) 卡關
Act II · --germline-hp-only flag 方案 (PON 當 germline 錨點)
    ↓  Phase 1 機制 ✓ 但 filter Gate FAIL；flag=on NG≥3=0 意外觀察
Act III · 用戶一句「NG=2 與甲基有關係嗎」觸發 C++ 原始碼回查
    ↓  HPFineNGroups 真機制 = 4-bucket phasing (非 methylation)
    ↓  obs18 跨 6 TO 樣本驗證 · Inner same-hap 93-99% · TP gap +0.37
TO-層論文主軸 pivot: methylation bimodality → LOH-constrained phasing
```

## 研究故事線（詳細版本 · 對應 slides）

### Act I · 定錨與方法學（Slides 1-9）

**目的**：建立「測量單位正確」+「filter 框架初步 POSITIVE」的信任

- **Slide 1-2**：封面 + Roadmap（定位本週在過去已關閉 × 下週優先之間）
- **Slide 3**：上週 0421 週報結論回顧（3 已解 + 2 blocker）
- **Slide 4**：本週提問 —— 四子問題 Q1-Q4（分別對應 Thread A/B/C/D）
- **Slide 5-7 · Thread A (CN KDE 校準)**：
  - 舊問題：`--expected-coverage 75.0` hardcoded 讓 CN tier 系統性偏移
  - 解方：KDE 雙 Pass 架構（Pass 1 並行 → Mid KDE peak detection → Pass 2 重算）
  - 驗收：HCC1395 pilot 53.0× bias −1.9% vs SEQC2 54×
- **Slide 8-9 · Thread B (LOH × AF × CN filter framework)**：
  - TP 熱圖雙極分佈（綠區 88-96% / 紅區 47-60%）
  - S1-S7 filter scheme：S3 Diploid Het TP 95.5%、S5 combo FP ↓99.37%、S4 ambiguous
  - Slide 9 伏筆：S4 含 75% TP + 76% FP → 需二級判別

### Act II · 卡關與方案（Slides 10-17）

**目的**：揭示「加特徵並未解決 S4」+「tag 方案機制對但 filter 無增益」的雙重卡關

- **Slide 10**：S4 bucket 5 維座標失效 → 引出 HP / ISM / 聚類特徵擴展
- **Slide 11**：HP 雜訊根因（0421 延續）—— LongPhase-TO 保留 somatic HP tags 污染 HP-derived 特徵
- **Slide 12 · Thread C 方案**：`--germline-hp-only` flag 設計 —— PON 當 germline 錨點，somatic demote
- **Slide 13**：Phase 1 機制驗證 ✓ —— Audit 獨立性（NHP_Somatic sum 兩 flag identical）
- **Slide 14**：Phase 1 AUC Gate FAIL —— 18 HP-related 特徵無一達 +0.02 Gate；AlleleDelta 對照 0 shift
- **Slide 15**：HPFineNGroups 分佈崩潰 —— flag=on 下 NG≥3=0（22,732 TP + 8,148 FP regions reassign）
- **Slide 16**：HP_Ratio shift 小 —— Plan R3 「HP_Ratio 歸零」條件未觸發 → bug 不在上游
- **Slide 17**：「tag 方案」定位 —— 按用戶決策 A·b，**技術可行但 filter 無增益**；Phase 2 判定前結論懸掛（避免「tag 可採用」肯定語）

### Act III · 機制揭露與 pivot（Slides 18-24）

**目的**：呈現本週偵探時刻與論文主軸 TO-層重定位

- **Slide 18**：偵探時刻 —— 用戶 04-22 晚間一句提問觸發 C++ 原始碼回查（時間軸 4h44m 轉折）
- **Slide 19**：HPFineNGroups 真機制（C++ 原始碼直接呈現）—— 4-bucket occupancy 純 phasing
- **Slide 20**：NG=2 四種組成（same_HP1/same_HP2/cross_het/cross_het_inv）
- **Slide 21 · obs18 跨樣本驗證**：Inner × NG=2 same-hap 93-99%（6/6 一致，median 97%）
- **Slide 22**：Inner vs Outer TP gap —— median +0.37（6/6 正向，range +0.05~+0.52）
- **Slide 23**：LOH-constrained phasing 機制圖（雙面板物理場景）+ TO-層論文主軸 pivot
- **Slide 24**：Thread C 為 Thread D 獨立 negative control —— somatic HP demote → NG≥3 消失 → NG 是 phasing 而非 methylation

### 結語（Slides 25-26）

- **Slide 25**：結論總表變動（6 項：新增 CL-LCP-001 ⭐5、CL-S3-001 ⭐4、CL-CN-KDE-001 ⭐4 / 降級 CL-016 ⭐4→⭐3 / 加註 CL-LAF-001 / 新增懸掛 CL-HP-ONLY-001 ⭐3）
- **Slide 26**：下週 P0/P1 行動（Paired obs18 對照 + Archive TO 6 樣本重跑 + HCC1954 outlier + Wilcoxon）

---

## 用戶決策約束（影響 PPT 敘事）

五個分歧點於 2026-04-23 plan 階段定案：

| 分歧點 | 決策 | PPT 反映 |
|------|------|---------|
| **A. 「tag 可採用」定位** | (b) 延後定論 | Slide 17 三張懸掛卡 + 結論 banner「技術可行但目前 filter 無增益」；**避免使用「tag 可採用」肯定語** |
| **B. 新 tag 重跑 ISM 觀察切分** | 發現即終點 | Slide 18-24 主軸建立在 flag=off 舊 tag obs17/obs18；**不等新 tag Phase 2 全樣本** |
| **C. 論文主軸 pivot 版本控制** | 以 TO 為本週重點（保留 paired 舊結論）| Slide 23 pivot banner 僅宣告 TO 層；**paired AF × NGroups POSITIVE 保留不撤回**，加註需獨立驗證 |
| **D. 週報涵蓋區間** | 2026-04-16 ~ 04-23（8 天）| 以 CN KDE 方法學為起點（Slide 5-7）；上週 Thread A 基礎設施三軌僅薄層引用（Slide 3） |
| **E. 執行模式** | 全自動（Phase 0-2/4-5 自動；3/6 Review）| 當前進入 Phase 3 Review + Phase 4 PPTX 生成 |

---

## 關鍵資料來源（所有數字 / 圖 / 引用的出處）

所有 slide 內的量化數字均可在 `source_weekly_report.md` 追溯至原實驗報告：

| 數字 / 結論 | 來源（相對於 repo root）|
|-----------|-------------------|
| HCC1395 KDE 53.0×、bias −1.9% | `docs/experiments/in_progress/2026/04/20260420_KDE_Fix_Acceptance_Validation_01.md` |
| CovM median 0.880 → 1.245 | 同上 |
| S1-S7 scheme 完整數字 | `docs/experiments/in_progress/2026/04/20260422_LOH_AF_KDE_TPFP_Discrimination_02.md` |
| 18 HP-related 特徵 AUC | `docs/experiments/in_progress/2026/04/20260421_ReadParser_GermlineHPOnly_HCC1395_01.md` |
| NG≥3 崩潰到 0 | 同上 |
| Audit 獨立性 NHP_Somatic sum | 同上 |
| HPFineN marker re-audit | `docs/experiments/in_progress/2026/04/20260421_HPFineNGroups_Marker_Reaudit_01.md` |
| obs18 Inner same-hap 93-99% | `research/tpfp_loh_af_kde_discrimination/09_TO_sample_af_lohside_ng.md` |
| TP gap +0.37 (6/6) | 同上 + `observations/step7_*.md` |
| C++ 原始碼引用 | `src/core/LabelTest.cpp:265-302`、`include/core/Stats.hpp:324` |
| Commits | `374fad4`、`12d9b3e`、`775036c`、`a61779c`、`4dc2d73`、`2e2df22`、`5abc659` |

---

## 圖片資產

所有 slide 使用的圖片均在 repo 內絕對路徑（`build_pptx.py` 使用 ROOT = `/big7_disk/liaoyoyo2001/InterSubMod`）：

| 類別 | 目錄 |
|------|------|
| Thread A · KDE | `docs/experiments/in_progress/2026/04/figures/20260420_kde_fix_acceptance/` |
| Thread B · Filter | `docs/experiments/in_progress/2026/04/figures/20260421_ng_kde_rescaled/` |
| Thread C · HP-only（本週新產）| `docs/experiments/in_progress/2026/04/figures/20260423_weekly_thread_c/` |
| Thread D · LCP | `research/tpfp_loh_af_kde_discrimination/figures/new/obs17/` 與 `.../obs18/` |

Thread C 4 張圖由 `scripts/analysis/20260423_weekly_thread_c_figures.py` 本週產出。

---

## PPT 設計常數（樣式一致性）

配色（ColorBrewer colorblind-safe）：
- 主色系：`dark=#1E2A44`（標題/文字）、`bg=#F7F3EC`（米色背景）、`accent=#A85540`（赤陶紅強調）
- 語意色：`positive=#009E73`（POSITIVE/PASS）、`negative=#D55E00`（NEGATIVE/FAIL）、`blue=#0072B2`（對照）
- 輔助色：`muted=#5E6572`、`en_color=#5E6572`、`line=#D6CCBF`、`card_bg=#FFFFFF`

字體：Arial（標題 Bold）· 中文主 + 英文副（EN 字號約為中文 60%）

版面：16:9（13.333 × 7.5 in）· 標題區 0-1.5in · 內容區 1.8-7.0in · footer 7.15-7.5in

每 slide 結構：
1. **頂部色條** + **Kicker**（章節標記）+ **中文標題** + **英文副標**
2. **主視覺區**（圖/表/機制圖）55-65%
3. **文字區**（≤4 bullet）≤15%
4. **Footer**：頁碼 + 週報名

---

## 未來 AI 修改 PPT 的指引

若未來 AI / 協作者需修改本 PPT：

1. **先讀 `README.md`（本文件）→ `00_director_storyboard.md` → `source_weekly_report.md`** 建立脈絡
2. **修改 `build_pptx.py`**：每 slide 一個 `slide_NN_XXX()` 函數；修改後執行 `python3 build_pptx.py` 重建
3. **渲染截圖驗證**：
   ```bash
   libreoffice --headless --convert-to pdf --outdir screenshots *.pptx
   pdftoppm -r 110 screenshots/*.pdf screenshots/slide -png
   ```
4. **避免違反用戶決策 A-E**（見上節「用戶決策約束」）：
   - 不用「tag 可採用」肯定語（決策 A）
   - Thread D 主軸不等新 tag 全樣本（決策 B）
   - paired 舊結論不撤回（決策 C）
5. **所有新增數字 / 結論 / 引用必須可在 `source_weekly_report.md` 找到對應**，避免 slide 與週報脫鉤
6. **修改視覺/口頭分配時，先更新 `00_director_storyboard.md` 對應 slide 章節，再修 `build_pptx.py`**

---

## 驗證清單（發布前檢核）

- [ ] 所有 26 slides 渲染 PNG 無缺圖 / 無文字溢出
- [ ] Slide 17 避用「tag 可採用」肯定語
- [ ] Slide 23 paired AF × NGroups 保留註明
- [ ] Slide 25 結論總表 tier 與週報 Layer 3 一致
- [ ] Slide 26 下週行動與週報 Layer 4 一致
- [ ] 所有引用 commit SHA 可在 `git log` 找到
- [ ] 所有圖片相對路徑可解析（build_pptx.py 使用絕對路徑，渲染無誤即通過）

---

**文件結束** · PPT 資料夾自包含性 ✓ · 可複製整個資料夾到新環境給新 AI 繼續設計修改
