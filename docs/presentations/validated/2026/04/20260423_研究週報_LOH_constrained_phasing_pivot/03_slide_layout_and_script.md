<!--
建立時間: 2026-04-23
用途: 佈局座標 + 雙語文案參考 · 供未來 AI 修改 PPT 時對照
備註: 實際座標與文字寫在 build_pptx.py 的各 slide_NN_* 函數，此文件為摘要對照
-->

# PPT Layout & Script Reference

## 設計常數（與 PPTX_PROTOCOL 對齊）

### 畫布
- Aspect ratio: **16:9** (13.333 × 7.5 inches)
- Background: `#F7F3EC` (米色)
- Top accent bar: 0.06" 高，`#A85540`（赤陶紅）

### 版面區（每 slide 通用）
| 區 | y 座標 | 高度 | 用途 |
|---|-----:|----:|------|
| Top accent bar | 0.00 | 0.06 | 章節色條 |
| Kicker | 0.18 | 0.30 | 章節標記（ACT I / II / III + Thread ID） |
| Title_zh | 0.50 | 0.70 | 中文主標題（26pt Bold） |
| Title_en | 1.15 | 0.35 | 英文副標（13pt 正常） |
| Content area | 1.80 | ~5.0 | 主視覺 + bullets + side card |
| Footer line | 7.15 | 0.02 | 分隔線（`#D6CCBF`） |
| Footer text | 7.18 | 0.25 | 頁名（左）+ 頁碼（右）9pt |

### 色彩（`THEME` dict）
| Key | Hex | 用途 |
|-----|-----|------|
| dark | `#1E2A44` | 標題/文字 |
| bg | `#F7F3EC` | 頁面背景 |
| bg_alt | `#EEE6DB` | 副區塊底色 |
| accent | `#A85540` | 強調色（赤陶紅） |
| positive | `#009E73` | POSITIVE / PASS |
| negative | `#D55E00` | NEGATIVE / FAIL |
| blue | `#0072B2` | 對照 / neutral |
| muted | `#5E6572` | 輔助文字 |
| en_color | `#5E6572` | 英文副標 |
| line | `#D6CCBF` | 分隔線 / 細框 |
| card_bg | `#FFFFFF` | 卡片底色 |

### 字體
- 中英文主字：`Arial` / `Arial Bold`（PPTX 內部；Wireframe 用 `DejaVuSans`）
- 程式碼：`Courier New`（PPTX）/ `DejaVuSansMono`（Wireframe）
- 字號範圍：標題 26-38pt / 主標 18-26pt / 內文 12-16pt / caption 9-11pt

### Wireframe 字型 fallback（見 `PPT_RENDERING_RULES.md`）
- Latin（code point < 0x2E80）→ `DejaVuSans.ttf`
- CJK（0x4E00-0x9FFF 等）→ `DroidSansFallbackFull.ttf`
- 程式碼 → `DejaVuSansMono.ttf`

## 雙語雙行規則

- **標題**：中文在上、英文在下；英文字號 60% + 縮排 0.25"
- **章節色條 / Kicker**：僅中文
- **小標題 / banner 文字**：可單語（視資訊量決定）
- **Speaker notes**：繁體中文為主，專有名詞保留英文

## 每 slide 關鍵文案（標題 zh + 標題 en）

| # | Title ZH | Title EN |
|---|---------|---------|
| 1 | NG=2 LOH-constrained phasing 發現與 TO-層論文主軸 pivot | Mechanism correction: from methylation bimodality to phasing signatures in tumor-only sequencing |
| 2 | ROADMAP · 研究脈絡 · 過去 → 本週 4 Thread → 下週優先 | Context · This-week deliverables · Next-week priorities |
| 3 | RECAP · 上週 0421 週報結論 · 3 已解 + 2 blocker | Settled / Outstanding from prior week |
| 4 | MOTIVATION · 從 LOH×AF×CN 切分出發 · 加入 HP/ISM 特徵能進一步區分？ | Extending the stratification beyond LOH × AF × CN |
| 5 | ACT I · Thread A · CN 方法學 · 舊做法共用 75× 讓 CN tier 系統偏移 | Old pipeline: hardcoded --expected-coverage 75.0 for all samples |
| 6 | ACT I · Thread A · 解方 KDE 雙 Pass · HCC1395 pilot 53.0× | Two-pass KDE calibration · bias −1.9% vs SEQC2 54× |
| 7 | ACT I · Thread A · KDE 修正量化驗收 | Distribution shift and category reclassification |
| 8 | ACT I · Thread B · TP 被切成雙極分佈 | TP rate heatmap reveals bipolar distribution |
| 9 | ACT I · Thread B · S3 Diploid Het TP 95.5% · S4 伏筆 | Biology-informed stratified filter schemes |
| 10 | ACT II · S4 卡關 → 擴展特徵 | Bridging to Thread C: HP / ISM / clustering features |
| 11 | ACT II · Thread C · TO self-phasing 污染 HP-derived 特徵 | Self-phasing contaminates HP bucket assignment |
| 12 | ACT II · Thread C · `--germline-hp-only` 方案 · PON 錨點 | Design: anchor HP assignment to germline phase-set blocks |
| 13 | ACT II · Thread C · Phase 1 機制驗證 Audit ✓ | Audit independence confirms mechanism correctness |
| 14 | ACT II · Thread C · Phase 1 AUC Gate FAIL | None of 18 HP-related features reach the +0.02 gate |
| 15 | ACT II · Thread C · HPFineNGroups 分佈崩潰 NG≥3=0 | NG≥3 categories → 0 in both TP and FP |
| 16 | ACT II · Thread C · HP_Ratio shift + R3 未觸發 | HP_Ratio shift minor; upstream hp_tag parsing is correct |
| 17 | ACT II · Thread C · Tag 方案定位 · 機制✓filter⏸ | Mechanism validated; filter hypothesis pending Phase 2 |
| 18 | ACT III · Thread D · 偵探時刻 · 用戶一句觸發 C++ 回查 | A single user question triggered the mechanism reinterpretation |
| 19 | ACT III · Thread D · HPFineNGroups 真機制 · 4-bucket count | C++ source code proves: pure phasing × variant-presence count |
| 20 | ACT III · Thread D · NG=2 的四種組成 | NG=2 composition determines biological meaning |
| 21 | ACT III · Thread D · obs18 跨 6 樣本 · Inner ≥93% same-hap | 6/6 TO samples: Inner NG=2 ≥93% same-hap |
| 22 | ACT III · Thread D · Inner vs Outer TP gap +0.37 · 6/6 正向 | TP rate gap: Inner same-hap vs Outer cross-het |
| 23 | ACT III · Thread D · LOH-constrained phasing 機制 + pivot | Two physical scenarios explain somatic vs germline leak |
| 24 | ACT III · Thread C 為 Thread D 獨立 negative control（H-D1） | flag=on NG≥3=0 as H-D1 negative control; H-D4 pending paired |
| 25 | 結論總表變動 · 3 新增 / 1 降級 / 1 加註 / 1 懸掛 | Conclusion tier table delta |
| 26 | 下週 P0/P1 · 4 項主要行動 · 2 週內定案 | Next-week P0/P1 priorities |

## 修改 PPTX 的流程（新 AI 入手時）

1. 先讀 `README.md` 與 `source_weekly_report.md` 建立脈絡
2. 找到要改的 slide 對應 `build_pptx.py` 的 `slide_NN_XXX()` 函數
3. 修改 shape 座標或文字
4. 執行 `python3 build_pptx.py` 重新 build
5. 執行 `python3 screenshot_all.py`，確認輸出 `[ISSUES] none detected`
6. 抽查 2-3 張關鍵 slide PNG
7. 更新本檔案（`03_slide_layout_and_script.md`）對應表格

## PPT 結構統計

- 總 slides: 26
- 帶圖片的 slides: 12（Thread A 3 + Thread B 2 + Thread C 4 + Thread D 3）
- 帶程式碼片段: 1（slide 19）
- 帶 banner：12（章節 + 結論 banner）
- 帶 speaker notes: 26（全部）
