# 01 — P0 Critical Decisions（3 題）

> **建立日期**: 2026-04-19

> **優先等級**：P0（結論反轉風險 / 必先解決）
>
> **依賴**：本分類 3 題須在進 Phase 2 A+D 前全部決策並執行

---

## 問題 P0-A: CovM 75.0 hardcoded bug 修正策略

**問題描述**：
`expected_coverage=75.0` hardcoded default fallback 在 master dataset 全 7 樣本 × 2 mode（Paired/TO）共用。KDE-based per-sample baseline 存在但未啟用。導致 Coverage_Multiple 跨樣本比較時樣本間 coverage 差異被扁平化。R-01 已決定「立即修 + 重跑」，但**修正的具體策略與重跑範圍仍需決策**。

**影響範圍**：
- 受影響結論：**C17**（LOH Subclone step3 CN1）, **C18**（HCC1954 CovM=0.733）, **C19**（Z3 zone 邊界）, **C20**（CovM vs CN r=0.831 / z_extreme 0.15%）, **C22**（Zone Z1-Z5 定義）
- 受影響功能：F2 Subclone Marker, F3 Zone-Aware（CovM 為核心變量）
- 若不處理的風險：進入 Phase 2 A+D 時 baseline 汙染，新結論也會受影響；論文 reviewer 會質疑跨樣本 CovM 數據

**現況證據**：
- Audit card: `cards/C20_CovM_Non_Independent.md`, `cards/C22_Zone_Aware.md`
- 相關 cross_cutting: `cross_cutting/CovM_Bug_Impact.md`
- 原始報告：`docs/experiments/in_progress/2026/04/20260420_CovM_Baseline_Accuracy_Validation_01.md`
- 程式碼：`src/core/*.cpp` KDE baseline 路徑未啟用

**候選處理方向**：

| 選項 | 方法 | 代價 | 風險 |
|------|------|------|------|
| **A（推薦）** | 立即 `/cpp-change` 6 步流程修 KDE + 7 樣本全量重跑（Paired+TO+PON-only） | 3-5 天（含重跑）| 低 — 修正路徑已知（`--expected-coverage` CLI 已有） |
| B | 僅修 code、不重跑；延後重跑至 Phase 2 A+D 啟動時一併 | 1 天 | 中 — Phase 2 期間結論仍不可信；用戶已明確要求立即重跑 |
| C | 不修 code，暫加「bug affected」警告標註於所有受影響結論 | 2 小時 | 高 — 違背 R-01 決策；論文撰寫時仍需修 |

**驗證標準（無論選 A/B/C）**：
- **必達**：`--expected-coverage` 不傳時自動 KDE 估計啟用；7 樣本 expected_coverage 值**非全部 75.0**
- **必達**：修正前後 per-sample z-score rank 相關 >0.8（確保相對排序穩定）
- **觀察**：C20 CovM vs CN r 變動範圍（預期 0.5-0.9）
- **觀察**：C22 Z1-Z5 zone cardinality 變化（預期 ±30% 內）
- **驗收指令**：`./scripts/run_batch_vcf_analysis.sh` + `scripts/analysis/cross_sample_coverage_audit.py`

**依賴關係**：
- **前置依賴**：無（可立即啟動）
- **被阻塞項**：
  - P0-C Haplotag（建議同一 rebuild 窗口）
  - C17/C20/C22 重算
  - Phase 2 A+D 全量驗證
  - P1-B FDR 補算（C16/C22 部分需 post-fix 數據）

**用戶決策**：[ ] A（推薦）[ ] B [ ] C [ ] Other: ___

---

## 問題 P0-B: Pooled OLS 重算範圍（C15/C16/C17）

**問題描述**：
Pooled OLS residualization 在「confound 變量與分組標籤相關」時**無法移除分組信號**（見 `feedback_pooled_ols_residualization_trap.md`）。C15/C16/C17 三個結論在 residualize 步驟中使用 pooled OLS，可能保留 AF-bin 或 NR-bin 分組信號。**C17 r=+0.705 是最高風險位置**。

**影響範圍**：
- 受影響結論：**C15**（LOH methylation failure）, **C16**（HPFineNGroups marker）, **C17**（LOH Subclone AF×Methylation）
- 受影響功能：F2 Subclone Marker（C16+C17 為核心）
- 若不處理的風險：F2 POSITIVE 結論可能反轉為 CONDITIONAL 或 NEGATIVE，論文核心貢獻之一被動搖

**現況證據**：
- Audit card: `cards/C15_LOH_Methylation_Failure.md`, `cards/C16_HPFineNGroups.md`, `cards/C17_LOH_Subclone_AF.md`
- 相關 cross_cutting: `cross_cutting/Pooled_OLS_Audit.md`（含 SOP）
- Feedback：`feedback_pooled_ols_residualization_trap.md`

**候選處理方向**：

| 選項 | 方法 | 代價 | 風險 |
|------|------|------|------|
| **A（推薦）** | **C17 優先**（與 CovM bug 耦合，CovM fix 後一併做）→ **C16 次之**（結合 P1-C bootstrap）→ **C15 最後**（若 C17/C16 穩固則低優先） | 3-5 天 | 低 — 漸進式驗證，失敗也不傷害全局 |
| B | C15/C16/C17 三者並行 within-group 驗證 | 1 週+ | 中 — 平行投入分散，失敗難定位根因 |
| C | 僅 C17（C15 已 NEGATIVE 穩固，C16 filter-based 不受影響） | 2-3 天 | 中-高 — 省時但 C16 的 residualized AUC 聲明仍待驗證 |

**驗證標準（無論選 A/B/C）**：
- **必達**：C17 按 AF bins [0.05, 0.1, 0.2, 0.3, 0.4] within-group OLS，每組 r 值 + meta p value 寫回卡片
- **必達**：within-group r 方向一致（正相關）且 meta p<0.01
- **門檻**：若 within-group r ≥ 0.5（7/7 一致）→ 結論穩固；若 r<0.3 → 結論反轉
- **驗收指令**：`scripts/analysis/af_bin_within_group_ols.py`（需新建或擴展既有腳本）

**依賴關係**：
- **前置依賴**：P0-A CovM fix（C17 step3 CN1 依賴乾淨 CovM）
- **被阻塞項**：P1-B FDR（C16 90 組合 FDR 需 within-group AUC）, P1-C bootstrap CI

**用戶決策**：[ ] A（推薦）[ ] B [ ] C [ ] Other: ___

---

## 問題 P0-C: Haplotag ReadParser 修正時機

**問題描述**：
P0-2 haplotag ReadParser 修正（HP integer tag parsing）尚未完全落地。C03/C05/C08 目前處於「汙染條件下的觀測」狀態。雖然 HP tag 修正已完成初版（見 `project_hp_integer_tag_fix.md`），但完整重跑與驗證未補做。

**影響範圍**：
- 受影響結論：**C03**（TO AUC ceiling）, **C05**（O12 LOH Scenarios）, **C08**（Read-level FP）
- 受影響功能：無直接 F1-F5 影響，但 C08 是 F5 Variant Confidence 的負面證據
- 若不處理的風險：NEGATIVE 結論的 confound 歸因可能不精確；但結論方向（NEGATIVE）預期不變

**現況證據**：
- Audit card: `cards/C03_TO_AUC_Ceiling.md`, `cards/C05_O12_LOH_Scenarios.md`, `cards/C08_Read_Level_FP.md`
- Memory: `project_hp_integer_tag_fix.md`（已完成初版）
- Phase 1 Agent B 觀察：P0-2 haplotag ReadParser 修正未完成

**候選處理方向**：

| 選項 | 方法 | 代價 | 風險 |
|------|------|------|------|
| **A（推薦）** | CovM bug 修正後**同一 rebuild 窗口**修 Haplotag + 7 樣本重跑（P0-A、P0-C 合併一次執行） | +1 天（在 P0-A 重跑上加） | 低 — 共用 rebuild + 重跑，省時且乾淨 |
| B | 獨立 sprint，單獨修正 + 重跑 | 3-4 天（額外） | 中 — 重複 rebuild/重跑，浪費時間 |
| C | 延後至 Phase 2 A+D 啟動時一併 | 0 天（現在） | 中 — Phase 2 A+D 新結論仍受 Haplotag 汙染 |

**驗證標準（無論選 A/B/C）**：
- **必達**：HP integer tag 解析成功率 100%（log 檢查）
- **必達**：修正前後 AMB% 或 HP tag 分佈無異常變動
- **觀察**：C03/C05/C08 AUC 變動（預期 <0.02 差異，NEGATIVE 結論穩定）
- **驗收指令**：`./build/tests/test_read_parser` + VCF 分析後的 HP tag 分佈 TSV

**依賴關係**：
- **前置依賴**：無（程式碼修正部分）；建議與 P0-A 共用 rebuild 窗口
- **被阻塞項**：無（NEGATIVE 結論方向不依賴此修正）；Phase 2 A+D 若啟動會受此影響

**用戶決策**：[ ] A（推薦）[ ] B [ ] C [ ] Other: ___

---

## P0 區段總結

**推薦總動作**（全選 A）：
1. `/cpp-change` 觸發：同時修 **CovM KDE baseline** + **Haplotag ReadParser**
2. 一次 rebuild + 7 樣本全量重跑（Paired + TO + PON-only）
3. CovM 修正完成後，立即啟動 C17 within-group OLS 驗證
4. C17 穩固後推動 C16（結合 P1-C bootstrap）+ C15

**推薦理由**：
- R-01 已決定立即修，推薦選項全部對齊
- 共用 rebuild 窗口最省時
- 漸進式驗證（C17→C16→C15）失敗可及早停損
- 修正後 Phase 2 A+D 才有乾淨 baseline

**若全選 A 後的輸出**：
- `src/core/*.cpp` KDE baseline 啟用 patch
- `src/core/ReadParser.cpp` Haplotag fix
- `output/master_dataset_v6/`（7 樣本 × 3 mode 重跑結果）
- `cards/C17/C20/C22.md` 重算數字更新
- `00_INDEX.md` CovM bug 狀態 🔴 → 🟢

---

## 關聯文件

- EXECUTIVE_DECISION_BRIEF.md 第 4-5 節
- `cross_cutting/CovM_Bug_Impact.md`
- `cross_cutting/Pooled_OLS_Audit.md`
- CHECKLIST.md
