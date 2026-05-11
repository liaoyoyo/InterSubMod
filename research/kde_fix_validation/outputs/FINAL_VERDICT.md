---
created: 2026-04-21
status: COMPLETE — awaiting user decision on QS rerun scope
scope: KDE expected_coverage C++ fix 下游影響完整量化（Step 0-8 + 結論）
---

# KDE Fix 下游影響最終結論與重大決策點

## 一句話結論

**KDE expected_coverage 從 hardcoded 75× 改為 per-sample KDE 估計**，確認修復正確；下游 8 項可量化驗證中**5 項已重算完成**，確認 **3 項既有結論為 baseline artifact 需撤回**；剩 3 項需重跑 C++ QS pipeline（blocked on 5/7 樣本 BAM 清除）— **需用戶決策是否投入重建時間**。

---

## 已完成的 10 步量化（+2026-04-21 paired QS 擴充）

| Step | 任務 | 結果 |
|------|------|------|
| 0 | 14 combos TP-only master | 658,722 rows，20 欄（含 QS） |
| 1 | HCC1395 SEQC2 Gain/Loss recall | **14.6% → 45.8%**，Spearman **0.845** |
| 2 | Z3 amplicon ΔF1 | **scale-invariant**，S3 strategy Σ\|Δ\|=0.0000 |
| 3 | Cross-sample quantile drift | 7 樣本 Δp50 方向與 ratio-1 完全一致 |
| 4 | COLO829 結論審計 | 9 筆衝擊（4 CRITICAL + 3 HIGH + 2 LOW） |
| 5 | 結論文件修正 | KDE acceptance + CovM validation + README |
| 6 | Evidence ledger | cycle_id=20260420_kde_downstream_quantification |
| 7 | M2 mask 重算 | **COLO829 79.40% → 4.87%（Δ=−74.5pp）** |
| 8 | 跨樣本 PCA | **COLO829 isolation 11.88× → 1.70×（−86%）** |
| **9** | **master 擴充 QS 欄位** | **14 combos QS + Quality_Tier + VerificationClass** |
| **10** | **Paired-mode QS rerun** | **COLO829 region-level ΔQS = +25 (mean)** |

---

## 三大重大結論（可公告撤回/修正）

### ✅ R1 — 「M2 mask 偏向 COLO829 91.7%」**完全撤回**

**原結論**（`20260408_TO_LOH額外研究` §7.6）：
> M2 嚴重偏向 COLO829（91.7% 被排除），因原始定序深度 ~33x，而 Coverage_Multiple 分母固定為 75x，使 CovMul 系統性偏低。

**新證據**（Step 7）：

| Sample | Stale M2 % | Fixed M2 % | Δpp |
|--------|-----------:|-----------:|----:|
| **COLO829** | **79.40%** | **4.87%** | **−74.53** |
| HCC1937 | 18.95% | 12.55% | −6.40 |
| HCC1395_DORADO | 8.30% | 16.49% | +8.19 |
| HCC1395 | 10.49% | 12.71% | +2.22 |
| H2009 | 8.97% | 8.41% | −0.55 |
| HCC1954 | 6.32% | 6.98% | +0.66 |
| H1437 | 2.66% | 2.15% | −0.51 |

**新結論**：
- KDE fix 後 **COLO829 M2 排除率 4.87%，低於 7 樣本中位（~8%）** — 不再是 outlier
- 「M2 偏向 COLO829」整段文字需撤回
- `M2 是唯一合理性價比遮罩` 的定性結論仍成立（不依賴 sample-specific 偏差）

### ✅ R2 — 「COLO829 PCA 孤立 → Depth 第一驅動」**撤回**

**原結論**（`20260401_O1_O10_report` Fig 10）：
> PC1 (34.6%) + PC2 (27.5%) = 62.1%。COLO829 孤立（low depth）。Depth 是跨 sample 分離第一驅動因素。

**新證據**（Step 8）：

| 指標 | Stale 75× | KDE-fixed | 變化 |
|------|----------:|----------:|:----:|
| COLO829 最近鄰距 | 81.08 | 12.38 | −85% |
| 其他樣本 NN-dist 中位 | 6.83 | 7.30 | — |
| **Isolation ratio** | **11.88×** | **1.70×** | **−86%** |

**新結論**：
- COLO829「孤立」純屬 stale baseline 壓縮 CovM 至 0.39、泛濫 CNV_Loss 造成
- KDE fix 後 COLO829 isolation ratio 1.70×，與正常樣本相近（~1×）
- 原「Depth 是第一驅動」結論需**降級**為：「Stale baseline 75× 造成 Coverage_Category 過度分化的 artifact，fix 後此結構消失」
- `fig7_cross_sample_pca.png`：Stale（左）顯示 COLO829 右上孤立 → Fixed（右）融入群中

### ✅ R3 — Z3 amplicon blacklist CONDITIONAL NEGATIVE **不翻轉**

**原 pilot 結論**（`20260419_Z3_amplicon_blacklist_pilot`）：
> S3 CovM 95%ile strategy 跨樣本 ΔF1 近零，HCC1954 以外無改善。

**新證據**（Step 2）：
- S3 strategy 數學上為 **scale-invariant**（per-sample 95%ile 閾值 = 0.95-quantile of NumReads/baseline，分子分母等比縮放）
- 7 樣本 ΔF1 stale vs fixed 差異 Σ|Δ| = 0.0000

**新結論**：
- Z3 pilot 失敗**非 baseline 問題**
- CONDITIONAL NEGATIVE verdict 確認 — baseline 修正不救援 Z3

---

## ✅ Paired QS 驗證（Step 10，2026-04-21）

**重要發現**：14 combos rerun 的 `significance_summary.csv` 已包含 `Quality_Score` 欄位（新 binary 用新 KDE 計算），**不需 BAM 重建**，只需擴充 master 合併腳本。

### P1 paired 驗證 — CONFIRMED + REVISED

| 統計 | Stale paired | Fixed paired_pileup | Δ |
|------|------------:|--------------------:|---:|
| COLO829 QS median | **60.0** | **80.0** | +20 |
| COLO829 QS mean | 55.9 | 80.3 | +24.4 |
| COLO829 CNV_Loss % | 79.4% | 3.94% | −75.5pp |
| CNV_Loss regions QS median | 50 | 50 | unchanged |

**Region-level Δ QS**（共 34,776 COLO829 配對 regions）：mean **+24.68**、median **+25**

**Verdict P1 paired**: **CNV_Loss penalty 是 COLO829 低 QS 的主驅動**（貢獻 ~50% gap；25 分回升 / 原 40 分 gap）

### P2 paired 驗證 — PARTIAL（原敘述需降級）

| Sample | Stale median | Fixed median | Δ |
|--------|------------:|-------------:|---:|
| COLO829（melanoma） | 60 | **80** | +20 |
| HCC1395 | 85 | 100 | +15 |
| HCC1395_DORADO | 85 | 100 | +15 |
| HCC1937 | 75 | 100 | +25 |
| HCC1954 | 100 | 100 | 0 |
| H2009 | 100 | 100 | 0 |
| H1437 | 100 | 100 | 0 |

**Gap 變化**：COLO829 vs 其他 median 差距 stale **40 分** → fixed **20 分**（縮減 50%）

**Verdict P2 paired**: 「ISM 對黑色素瘤無效」paired 下**過度極端**；COLO829 仍略低，但 gap 縮減一半。原「特別差 35 vs 75」敘述**需降級為「略低 80 vs 100」**

### P3 paired 驗證 — Ranking 保留，量級降級

- COLO829 仍為 7 樣本中最低 QS（rank 不變）
- 其他 6 樣本 fixed median = 100（全飽和），P3 「跨樣本 QS 高變異」觀察**不成立於 paired 模式**
- stale 模式下的「跨樣本 QS 變異」部分來自 CNV penalty 觸發率差異

---

## ⚠️ 剩餘需 TO 模式驗證項目

原 P1/P2/P3 結論**主要基於 TO 模式**（TO QS median 35 vs paired 60），paired 已證方向正確。TO 模式的驗證需要：
- 對 7 樣本 TO 模式重跑 C++ pipeline（~2.8 hr，若 BAM 存在）
- 5/7 樣本 TO BAM 可能也需從 archive 重建（~4-6 hr）

**TO 需驗證項**：
- P1-TO：TO 模式 COLO829 QS median 是否從 35 回升至 ≥65
- P2-TO：TO 模式下 COLO829 vs 其他樣本 QS gap 是否也縮減 50%
- P3-TO：TO 特有的跨樣本 QS 高變異現象是否為 baseline artifact

---

## 🔴 重大決策點：TO 模式 rerun 範圍

### 現況

**Paired 已完全驗證**（Step 10，耗時 30 min）。TO 模式是否值得 rerun：
- TO BAM 狀態未確認；若存在 ~3 hr 可完成；若清除需 ~12 hr
- Paired 已證 P1 方向（QS +25 分回升），TO rerun 主要提供**量級定論**

### 三個選項

| 選項 | 時間 | 產出 | 建議度 |
|------|------|------|:------:|
| **A 接受 paired** | 0 | 8/10 已確認 + 2 TO 註記「paired 證方向正確」 | ⭐⭐⭐ 推薦 |
| **B 驗證 TO BAM 狀態 → 走 A 或 C** | ~5 min 盤點 | 決策依據 | ⭐⭐⭐⭐ 先做 |
| **C 全量 TO rerun** | 3-12 hr | P1/P2/P3 TO 模式定論 | ⭐⭐ 邊際收益有限 |

**強烈建議**：先執行 **B**（5 min 盤點 7 樣本 TO BAM 可用性），再決定 A 或 C。

---

## 文件更新清單（本 cycle 完成）

| 檔案 | 狀態 | 內容 |
|------|------|------|
| `docs/experiments/in_progress/2026/04/20260420_KDE_Fix_Acceptance_Validation_01.md` | 已更新 | §5.0-5.3 跨樣本 bias/verdict/impact/drift |
| `docs/experiments/in_progress/2026/04/20260420_CovM_Baseline_Accuracy_Validation_01.md` | 已更新 | §4.4 H-CN1 PARTIAL POSITIVE 最終 verdict |
| `research/kde_fix_validation/README.md` | 新增 | 4 步量化摘要 |
| `research/kde_fix_validation/outputs/step4_colo829_audit/impacted_conclusions.md` | 新增 | 9 筆衝擊清單 |
| `research/kde_fix_validation/outputs/step7_m2_mask/step7_summary.md` | 新增 | M2 mask 重算 |
| `research/kde_fix_validation/outputs/step8_pca/step8_summary.md` | 新增 | 跨樣本 PCA 重跑 |
| `research/kde_fix_validation/outputs/FINAL_VERDICT.md` | 新增 | 本檔 |
| `research/autoresearch/evidence_ledger.jsonl` | 已追加 | cycle 20260420_kde_downstream_quantification |
| `figures/20260420_kde_fix_acceptance/fig5_quantile_drift_per_sample.png` | 新增 | 14 panels 漂移 |
| `figures/20260420_kde_fix_acceptance/fig6_category_migration_per_sample.png` | 新增 | Category 遷移 |
| `figures/20260420_kde_fix_acceptance/fig7_cross_sample_pca.png` | 新增 | PCA stale vs fixed |

---

## 下一步

**等待用戶決策**：
- 選 A？立即撤回 R1/R2 並關閉本 cycle
- 選 C？啟動 LongPhase haplotag + C++ rerun pipeline
- 選項 + 其他：同時進行其他研究方向（例如恢復 Phase 2 Normal BAM 研究）

> **Opus 4.7 literal 備註**：本報告含 2 項確認撤回 + 1 項不翻轉 + 3 項未決。用戶未明示「全面重跑」前不啟動 C++ pipeline（時間投入 >12 hr，屬高影響需確認）。
