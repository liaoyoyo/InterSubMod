# CHECKLIST — 30 題決策速查

> **建立日期**: 2026-04-20
>
> **更新**: 2026-04-19 執行進度

> **用途**：一頁勾選全部 30 個決策點；填完後可直接轉為 TodoWrite 執行
>
> **格式**：每題列推薦選項；勾選完後回報「全選 A」或列出例外

---

## 執行進度（2026-04-21 更新）

### ✅ 已完成（免 Hard Gate）

| 項目 | 輸出 | 日期 |
|------|------|------|
| P1-A r=0.997 誤標修正 | 20260417_Zone_Aware:49,73 + 4 cards | 2026-04-20 |
| P2-A LOH 層次註解 | CURRENT_FOCUS + 02/04/05/07 補註 | 2026-04-19 |
| P3-A HPFineNGroups filter 更新 | CURRENT_FOCUS:101 新 filter | 2026-04-19 |
| P3-B 22 結論納入 06 穩定性審查 | 06 補「結論 17-22」章節 | 2026-04-19 |
| P3-C TODO 編號統一 | 06 章節頭加對照表 | 2026-04-19 |
| P3-D Archive 舊數據統一公告 | `docs/archive/OLD_DATA_CORRECTION_NOTICE.md` | 2026-04-19 |
| P3-E PRE-FIX 警告（umbrella） | 2× `output/*_output_archive/README.md` | 2026-04-19 |
| R-04 13 NO-GO 卡補註 | 9 cards 加 R-04 補註 | 2026-04-19 |
| R-05 Pooled OLS 全面盤點 | `cross_cutting/Pooled_OLS_Global_Inventory.md` | 2026-04-19 |
| P2-B AUC 0.58 門檻統一 | `docs/standards/auc_threshold_definition.md` | 2026-04-19 |
| **P0-A C++ 修正**（C1 logging + C2 audit column） | commits `ec0608b→5abc659` + config.sh typo fix `cdcd8a6` | **2026-04-19 ~ 2026-04-20** |
| **P0-C 方法論審查** | `docs/methodology/20260420_ReadParser_HP_PS_integration_01.md`（推薦 B 平行欄位） | **2026-04-20** |
| **P0-B C16 Within-Group Validation** | `scripts/analysis/methylation_cn_within_group_validation.py` + C16 card 補註；**REAL_SIGNAL 確認**（10/10 bin CI 排除 0, median \|r\|=0.178） | **2026-04-21** |
| **P0-B C15 Script Audit** | C15 card 補註：O15 scripts 已為 stratified (mode × loh_group × truth_label) 設計，無 Pooled OLS trap | **2026-04-21** |
| **P0-B C17 Script Audit** | C17 card 補註：step2 per-sample Mann-Whitney + step3 per-sample Spearman，無 LinearRegression/residualize 呼叫 | **2026-04-21** |

### 🟡 進行中

| 項目 | 狀態 | 預估 |
|------|------|------|
| **P0-A 7 樣本 rerun** | HCC1395 ✅ / HCC1395_DORADO 進行中（InterSubMod 階段）/ 餘 5 樣本 | 剩 ~1 天 |

### 🔴 待 Hard Gate 確認（需用戶明示）

| 項目 | 代價 | 阻塞 | 備註 |
|------|------|------|------|
| **P0-C /cpp-change 執行**（選 B 實作） | +0.5 天 | 29 HP-dependent 特徵 opt-in | 方法論已審查完成，等 P0-A rerun 結束避免 binary 版本混淆 |

### 🟡 P0 完成後可並行（3-7 天）

- P1-C: 1000× bootstrap CI + per-bin + 結合 P0-B
- P1-B: BH-FDR C16→C22→O11/O12/O13→C07
- P2-D: NG=3 NR=120 pilot
- P2-E: LongPhase-TO 源碼分析

### ⏳ Phase 2-dependent（延後）

- P2-C: ASM germline/somatic（綁 Phase 2 A+D）
- Q-05: Phase 2 A+D M1/M2/M3 成功標準定義
- R-06: F1-F2 核心結論文獻深度補齊（與論文同步）

---

## P0 Critical（3 題，必先決策）

### P0-A: CovM 75.0 hardcoded bug 修正策略
- [x] **A（推薦）** 立即 `/cpp-change` 修 KDE + 7 樣本全量重跑（3-5 天）
- [ ] B 僅修 code，延後重跑至 Phase 2 A+D 啟動時
- [ ] C 暫加警告標註，保留 bug
- [ ] Other: ___
- **詳細**：`decisions/01_P0_critical_decisions.md` § P0-A

### P0-B: Pooled OLS 重算範圍（C15/C16/C17）
- [ ] **A（推薦）** C17 優先（CovM bug 耦合）→ C16 次之 → C15 最後（3-5 天）
- [x] B 三者並行 within-group 驗證
- [ ] C 僅 C17
- [ ] Other: ___
- **詳細**：`decisions/01_P0_critical_decisions.md` § P0-B

### P0-C: Haplotag ReadParser 修正時機
- [x] **A（推薦）** CovM bug 修正後同一 rebuild 窗口修（+1 天）
- [ ] B 獨立 sprint
- [ ] C 延後至 Phase 2 A+D
- [ ] Other: ___
- **詳細**：`decisions/01_P0_critical_decisions.md` § P0-C

---

## P1 High（3 題）

### P1-A: r=0.997 誤標修正範圍
- [x] **A（推薦）** 修主源 + 4 卡引用 + grep 全 docs（30 分鐘）
- [ ] B 僅修主要文件
- [ ] C 打補丁加註保留原值
- [ ] Other: ___
- **詳細**：`decisions/02_P1_high_decisions.md` § P1-A

### P1-B: FDR 校正範圍
- [x] **A（推薦）** C16 優先 → C22 → O11/O12/O13 → C07（1 週）
- [ ] B 全 7 卡同時 BH-FDR
- [ ] C 僅 POSITIVE 結論
- [ ] Other: ___
- **詳細**：`decisions/02_P1_high_decisions.md` § P1-B

### P1-C: HPFineNGroups bootstrap CI
- [x] **A（推薦）** 1000× bootstrap + per-bin CI + 結合 P0-B（2-3 天）
- [ ] B 僅整體 CI
- [ ] C 省略
- [ ] Other: ___
- **詳細**：`decisions/02_P1_high_decisions.md` § P1-C

---

## P2 Precision（5 題）

### P2-A: LOH 層次混用精確化
- [x] **A（推薦）** 兩檔補註 + 全文 grep 62%（20 分鐘）
- [ ] B 新建術語表文件
- [ ] C 僅 C09 audit card 加註
- [ ] Other: ___
- **詳細**：`decisions/03_P2_precision_decisions.md` § P2-A

### P2-B: AUC 0.58 門檻定義統一
- [x] **A（推薦）** 建立 `docs/standards/auc_threshold_definition.md`（2 小時）
- [ ] B 只在新產出強制統一
- [ ] C 每結論保留原 cutoff 附註腳
- [ ] Other: ___
- **詳細**：`decisions/03_P2_precision_decisions.md` § P2-B

### P2-C: ASM germline vs somatic 區分
- [x] **A（推薦）** 綁 Phase 2 A+D 後區分；現階段註「pending」（2-3 週）
- [ ] B 用 imprinting region 先驗部分區分
- [ ] C 不做區分，論文改為「整體分佈」
- [ ] Other: ___
- **詳細**：`decisions/03_P2_precision_decisions.md` § P2-C

### P2-D: NG=3 非單調生物機制
- [x] **A（推薦）** 雙軌討論 + NR=120 pilot（1 週）
- [ ] B 僅加假說討論
- [ ] C 不解釋，僅 focus NG=4
- [ ] Other: ___
- **詳細**：`decisions/03_P2_precision_decisions.md` § P2-D

### P2-E: Self-Phasing 演算法層面根因
- [x] **A（推薦）** 閱讀 LongPhase-TO 源碼 + 補機制章節（3-5 天）
- [ ] B 僅引用文獻
- [ ] C 黑箱觀察
- [ ] Other: ___
- **詳細**：`decisions/03_P2_precision_decisions.md` § P2-E

---

## P3 Documentation Sync（5 題）

### P3-A: HPFineNGroups filter 版本更新
- [x] **A（推薦）** 修 `CURRENT_FOCUS.md:98-103` + grep 全 docs（20 分鐘）
- [ ] B 僅修 CURRENT_FOCUS
- [ ] C 不改，加 TODO
- [ ] Other: ___
- **詳細**：`decisions/04_P3_documentation_sync.md` § P3-A

### P3-B: 22 結論納入 06_結論穩定性審查
- [x] **A（推薦）** 新增「2026-04-18 補充結論 17-22」章節（1 小時）
- [ ] B 重寫 06_審查為統一主表
- [ ] C 不改，引用指向 audit card
- [ ] Other: ___
- **詳細**：`decisions/04_P3_documentation_sync.md` § P3-B

### P3-C: TODO 編號統一
- [x] **A（推薦）** 採 audit card 的 P0-A/P0-B/P0-C 為標準（30 分鐘）
- [ ] B 僅在 06_審查加註
- [ ] C 不統一
- [ ] Other: ___
- **詳細**：`decisions/04_P3_documentation_sync.md` § P3-C

### P3-D: Archive 舊數據 grep 修正
- [x] **A（推薦）** grep 全 docs + 逐處加「舊數據已矯正」註記（1 小時）
- [ ] B 僅修 CURRENT_FOCUS / INDEX
- [ ] C 不改，archive 視為歷史
- [ ] Other: ___
- **詳細**：`decisions/04_P3_documentation_sync.md` § P3-D

### P3-E: PRE-FIX 警告標示
- [x] **A（推薦）** 每 PRE-FIX 目錄新增 `_WARNING_PRE_FIX.md`（2 小時）
- [ ] B 建立 `bip8_output_archive/README.md`
- [ ] C 不標示
- [ ] Other: ___
- **詳細**：`decisions/04_P3_documentation_sync.md` § P3-E

---

## R Strategic（4 題）

### R-03: 期刊目標選擇
- [x] **A（推薦）** 暫不決定，Phase 2 完成後評估（3 個月後重議）
- [ ] B 鎖定 Genome Biology
- [ ] C 鎖定 Bioinformatics
- [ ] Other: ___
- **詳細**：`decisions/05_strategic_R_decisions.md` § R-03

### R-04: 13 NO-GO 結論是否回溯 within-group
- [x] **A（推薦）** 接受現有 NEGATIVE，每卡補註（1 小時）
- [ ] B 全 13 NO-GO within-group 驗證
- [ ] C 僅 C10/C14
- [ ] Other: ___
- **詳細**：`decisions/05_strategic_R_decisions.md` § R-04

### R-05: Pooled OLS 全面回溯範圍
- [x] **A（推薦）** grep 盤點所有使用點（1 天）
- [ ] B 盤點 + 全部 pooled-only 重跑
- [ ] C 僅 P0-B 處理
- [ ] Other: ___
- **詳細**：`decisions/05_strategic_R_decisions.md` § R-05

### R-06: 知識庫交叉驗證範圍
- [x] **A（推薦）** 漸進補齊 F1-F2 核心結論（1 週）
- [ ] B 全 22 系統化查證
- [ ] C 僅被 challenge 時補
- [ ] Other: ___
- **詳細**：`decisions/05_strategic_R_decisions.md` § R-06

---

## Q Open Questions（10 題）

### Q-01: E→A+D→B→C 優先序在 CovM bug 後是否成立？
- [x] **A（推薦）** 維持排序，P0-A 後設 rank correlation 決策點
- [ ] B 立即重排為 A+D 優先
- [ ] C 等 E 完全完成
- [ ] Other: ___

### Q-02: Zone 若失效是否重新定義？
- [x] **A（推薦）** P0-A 後依 cardinality 變化決定（≤30% 維持）
- [ ] B 預先定義新 zone
- [ ] C 不變 zone 定義
- [ ] Other: ___

### Q-03: ISM characterization 具體下游應用？
- [ ] **A（推薦）** 三方向並列（臨床 + 生物學 + 工具）
- [ ] B 僅臨床
- [ ] C 僅工具
- [x] Other: __ 生物學 + 工具 + 研究應用_

### Q-04: TO 甲基化增益為負是否獨立 negative result 章節？
- [x] **A（推薦）** 論文新增 NEGATIVE subsection
- [ ] B 保留散落於各卡
- [ ] C 放入 Discussion Limitations
- [ ] Other: ___

### Q-05: Phase 2 A+D 成功標準？
- [x] **A（推薦）** 3 層標準 M1/M2/M3
- [ ] B 不預定
- [ ] C 僅 M1
- [ ] Other: ___

### Q-06: 13 NO-GO 是否值得 within-group validation？
- [x] **A（推薦）** 僅 C08（F5 核心 NEGATIVE）
- [ ] B 全 13 NO-GO
- [ ] C 不補
- [ ] Other: ___

### Q-07: HCC1395 單樣本結論是否需其他 5 樣本驗證？
- [x] **A（推薦）** 核心 7/7 驗證，其他補 cross-sample context 敘述
- [ ] B 全單樣本結論在 COLO829 + H2009 驗證
- [ ] C 不做 generalizability 聲明
- [ ] Other: ___

### Q-08: ASM germline/somatic 區分是否依賴 Normal BAM？
- [x] **A（推薦）** Normal BAM 為主 + imprinting 補充
- [ ] B 僅 imprinting 先驗
- [ ] C 僅 Normal BAM
- [ ] Other: ___

### Q-09: NG=3 非單調：生物機制 vs artifact？
- [x] **A（推薦）** Pilot 後依數據決定（H1 branching / artifact）
- [ ] B 預先採 H1 敘事
- [ ] C 不解釋機制
- [ ] Other: ___

### Q-10: LOH Subclone 臨床 downstream cohort？
- [x] **A（推薦）** Cell line proof-of-concept + future work
- [ ] B 啟動 patient cohort pilot
- [ ] C 引用公開 WGBS 間接驗證
- [ ] Other: ___

**Q 題詳細**：`decisions/06_open_questions_Q.md`

---

## 快速填表（全選 A 範例）

若全部採納推薦：
```
P0: A, A, A
P1: A, A, A
P2: A, A, A, A, A
P3: A, A, A, A, A
R:  A, A, A, A
Q:  A, A, A, A, A, A, A, A, A, A
```

**回報格式**：
- 「全選 A」— 直接執行推薦方案
- 「除 X-Y 選 B 外其餘全選 A」— 列出例外
- 「X-Y 選 Other: ...」— 提出第四方案

---

## 決策執行優先序（選 A 後的時程）

### 立即可做（<1 天）
1. P1-A r=0.997 修正（30 分鐘）
2. P2-A LOH 層次註解（20 分鐘）
3. P3-A HPFineNGroups filter 更新（20 分鐘）
4. P3-C TODO 編號統一（30 分鐘）
5. P3-B 22 結論入 06_審查（1 小時）
6. R-04 NEGATIVE 卡補註（1 小時）
7. R-05 grep pooled OLS 盤點（1 天）

### 短期（1-7 天）
8. P0-A CovM fix + P0-C Haplotag（同 rebuild 窗口，3-5 天）
9. P0-B C17 within-group OLS（P0-A 後立即，2-3 天）
10. P1-C bootstrap CI（結合 P0-B，+1 天）
11. P2-E LongPhase-TO 源碼分析（3-5 天）
12. P3-D Archive 舊數據 grep（1 小時，可與其他並行）
13. P3-E PRE-FIX 警告（2 小時，可與其他並行）

### 中期（1-2 週）
14. P1-B FDR 校正（C16 優先，1 週）
15. P2-D NG=3 NR=120 pilot（1 週）
16. P2-B AUC 門檻文件（2 小時文件 + 逐卡修正）
17. Q-06 C08 within-group（2-3 天）

### 長期（綁 Phase 2）
18. P2-C ASM germline/somatic（綁 Phase 2 A+D）
19. Q-05 Phase 2 A+D 啟動（含 M1/M2/M3 標準）
20. R-06 知識庫深度補齊（與論文同步）

### 延後
- R-03 期刊目標（3 個月後）
- Q-01/Q-02 依 P0-A 結果觸發
- Q-03/Q-04 論文撰寫時決定
- Q-07/Q-10 論文 Discussion 決定

---

## 關聯文件

- 主決策文件：`EXECUTIVE_DECISION_BRIEF.md`
- 分類詳述：`decisions/00_research_path_tree.md` ~ `06_open_questions_Q.md`
- 跨題分析：`cross_cutting/*.md`
- Audit cards：`cards/C01-C22.md`

---

## 用戶回報欄

```
決策日期：YYYY-MM-DD
P0: [  ,  ,  ]
P1: [  ,  ,  ]
P2: [  ,  ,  ,  ,  ]
P3: [  ,  ,  ,  ,  ]
R:  [  ,  ,  ,  ]
Q:  [  ,  ,  ,  ,  ,  ,  ,  ,  ,  ]

備註：
```
