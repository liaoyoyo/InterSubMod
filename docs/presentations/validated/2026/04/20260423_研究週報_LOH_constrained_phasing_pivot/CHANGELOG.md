# CHANGELOG — 20260423 PPT Pivot

## v3 · 2026-04-23 晚間（S5 + 7 頁深化修正 + 資料陷阱發現）

### 7 頁深化修正（依用戶要求）
- **S8** 換用新 KDE 4-step 視覺化圖（raw histogram → Gaussian smooth → 2nd-deriv peak → CN tier）
- **S10** S14 scheme 定義前移至此位置 + 每 scheme 一句中文意義解釋
- **S11** (原 S10) 改為 TO-only heatmap（6 TO 樣本 × S1-S7，移除 paired_full）
- **S12-S14** 順延並更新 footer 編號
- **S17** S4 Ambiguous 分析圖（B4 ROC + TP/FP in AF×AlleleDelta plane），確認數據 AUC 0.717
- **S18** 三類特徵類別加「預期數據意義」補充（HP / ISM 甲基 / ISM 聚類）
- **S20** Paired vs TO schematic（Chromosome + Haplotype + SNV 圖示對照）
- **S23** PON 錨點 before/after schematic（圖示案例 + 驗證策略）

### S5 AF 分佈與 HCC1395 資料陷阱發現與修正
**根因**：兩項 merged 合成檔資料陷阱被發現：
1. `merged_7samples_paired_full_plus_hcc1395_to.tsv.gz` 的 `AF` 欄位**非 caller_af**（5 樣本 p75<0.06，視覺全堆底）
2. HCC1395 `kde_status='phase1_new'` rows 的 LOH annotation 殘缺（Inner 5.7% vs 正常 58.8%）

**修正**：S5 腳本統一改從 stale master archive (`all_region_rows.tsv.gz`) 讀 6 樣本 TO，用真 `caller_af` 作 Y 軸、`NumReads/per-sample KDE bx` 作 X 軸。結果：
- HCC1395 vs HCC1395_DORADO 分佈高度吻合（同細胞株驗證）
- 6/6 樣本 Inner TP rate > baseline（Δ+0.6 至 +3.1pp，綠色標示）
- AF 分佈 0-1 合理分散（真 VAF）

**Panel 標題**：新增 `base=X.XX, Δ+Xpp` 讓 PI 看到相對 baseline 的 LOH 增益

**詳細 postmortem**（永久紀錄）：
- KB: `knowledge/05_data_formats/06_merged_dataset_pitfalls.md`
- KB: `knowledge/02_samples/01_hcc1395.md`（phase1_new 警告）
- KB: `knowledge/10_research_status/04_blockers_and_risks.md` B0

### 其他修正
- S19 `↓` 箭頭 text overflow 修正（box 高度 0.15 → 0.25）
- 全部 9 張使用的 figures 複製到 `figures/` 子資料夾
- 新增 `FIGURES_INDEX.md` Slide ↔ Figure 對照

---

## v2 · 2026-04-23 下午（線性主軸重構 + B1-B7 整合）

### 背景
用戶 feedback「資料還是太雜亂 · 圍繞在主軸敘述」，同時 2026-04-23 一天內新跑完 7 份 B1-B7 + Phase 3 Synthesis 驗證，解決了 v1 週報早上標示的多個「下週 P0」懸掛。

### 結構變更

| 維度 | v1 | v2 |
|------|----|----|
| 敘事結構 | 三幕偵探故事 (Act I/II/III) | 線性主題 7 段 (段 1-7) |
| 開場 | 4 問題 Q1-Q4 導入 | 背景觀察（LOH.bed 佔比 / AF 分佈） |
| CN KDE | 3 slides (Act I) | 3 slides (段 2 獨立) |
| 現象確認 | 分散（Act I S8 + Act III S21） | 集中段 3（S10-S13） |
| Self-phasing | Act II 中段（8 slides 含 Phase 1 filter） | 段 6 後段（5 slides 純機制說明） |
| NG=2 機制 | Act III 偵探高潮 | 併入段 3 現象確認 |
| pivot banner | Slide 23 有 | **移除** |

### 結論升級（B1-B7 + Phase 3 Synthesis）

| 原 v1 狀態 | v2 新狀態 | 依據 |
|-----------|-----------|------|
| Thread B HCC1395 pilot | **6/6 TO POSITIVE · S1/S5 median 0.876** | Phase 3 Synthesis |
| H-D4 「待下週 P0」懸掛 | **已驗證 · TO-specific 確認** | B3 Paired gap=0.00003, p=0.578 |
| Thread D ⭐5 TO-only | **⭐5 TO-specific 確定性** | B1 + B3 雙重證實 |
| HCC1954 outlier 「根因未知」 | **caller FP 背景 84% 主導** | B2 |
| S4 二級判別「下週 P2」 | **AF+AlleleDelta LR AUC 0.717 pilot ✓** | B4 |
| COLO829 S1 fold 低「待解」 | **small-n Wilson CI 寬 artifact 非 subclone 缺失** | B5 |
| LOH_Noise 分類「待決」 | **保留獨立單元 · HCC1395 TO TP 0.930** | B7 |

### PPT 新/改 slide 清單

完全新增：S2 Roadmap 線性 / S3 LOH.bed 佔比 / S4 AF 分佈 / S6 HCC1395 背景 / S9 各樣本新 CN / S10 Phase 3 heatmap / S13 B1/B2/B3 統計卡 / S17 B4 secondary / S18 下一步計劃 / S20 Paired vs TO 根因 / S21 HP tag 值表 / S22 三路徑影響 / S23 PON + 驗證 / S26 小結

完全改寫：S1 Cover 新標題 · S5 LOH×AF×CN 初步（改為主軸起點）· S8 KDE 驗收 · S11 obs18 · S12 機制圖 · S14 scheme 定義 · S15 Phase 3 雙 heatmap · S16 S3+LOH_Noise+COLO829 · S19 HP tag 需要 · S24 下週 P0/P1 · S25 其他目標

保留既有設計：helpers (add_rect / add_text / add_image fit-within / add_header / add_footer / Latin+CJK font fallback)

### 驗證結果

`python3 screenshot_all.py` → `[ISSUES] none detected · all text fits · all images loaded ✓`

修正：S19 `↓` 箭頭 box 高度從 0.15 增至 0.25 消除 text overflow。

---

## v1 · 2026-04-23 上午（初版）

三幕偵探故事結構 26 slides，Thread A/B/C/D 架構，pivot banner 於 Slide 23。詳見 git log。
