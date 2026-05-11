---
id: ism-kb-03-pipelines-f1-baseline-canonical
name: "F1 Baseline Canonical（Single Source of Truth）"
description: "ΔF1 locked 數字的唯一權威：Phase 1A paired-pure +0.0112、TO -0.0206、V5 F1=0.7154；含完整 provenance（運行條件、樣本、方法、驗證、限制、lock 日期）。其他 KB 文件應連結至此，不複製數字。"
status: active
last_verified: 2026-04-22
content_nature: frozen-decision
doc_type: reference
verified_scope: "F1 numbers against docs/reports/research_landscape/10_Research_Chain_Registry.md CL-002 + docs/CURRENT_FOCUS.md Phase 1A locked"
related_ids:
  - ism-kb-03-pipelines-index
  - ism-kb-03-pipelines-paired-full
  - ism-kb-03-pipelines-tumor-only
  - ism-kb-03-pipelines-pipeline-comparison
  - ism-kb-08-truth-and-benchmark-f1-calculation
  - ism-kb-10-research-status-current-focus-snapshot
  - ism-kb-09-conclusions-positive-findings
  - ism-kb-06-workflows-analysis-scripts-index
tags: [pipeline, f1, baseline, canonical, sot, provenance, locked]
canonical_paths: [03_pipelines/05_f1_baseline_canonical.md]
alias_paths: []
---

# F1 Baseline Canonical（Single Source of Truth）

- 一句結論：**所有 ΔF1 數字在 KB 內只在此頁定義**；其他位置必須連結至此表，不重複複製數字
- 適用對象：任何引用 F1/ΔF1 數字者、驗證這些數字來源者
- 可直接執行命令（驗證日期：2026-04-22）：
  ```bash
  # 回溯 CL-002 權威定義
  grep -A 15 "CL-002" /big7_disk/liaoyoyo2001/InterSubMod/docs/reports/research_landscape/10_Research_Chain_Registry.md
  ```

---

## 🎯 SoT 聲明

本文件為 **Single Source of Truth**（唯一權威）。規則：

1. ✅ **引用時**：寫「ΔF1=+0.0112（見 [03_pipelines/05_f1_baseline_canonical.md#cl-002](05_f1_baseline_canonical.md)）」
2. ❌ **禁止**：在其他文件複製 CI、method、樣本清單、lock 日期等 provenance 欄位
3. 🔄 **更新規則**：任何 ΔF1 數字變動**只改本頁**；引用頁不需動
4. ⚠️ **狀態翻轉**（locked → withdrawn 等）：同步 `docs/CURRENT_FOCUS.md` + 本頁 + CHANGELOG

---

## 📊 當前 Locked F1 數字（2026-04-22）

### CL-002：Paired-pure ΔF1 ⭐5
<a id="cl-002"></a>

| 欄位 | 值 |
|------|-----|
| **Claim** | ClairS Paired mode 下，ISM 作為 secondary filter 的 F1 delta 上限為 +0.0112 |
| **Sample** | HCC1395（ONT 5kHz simplex + 5mCG/5hmCG） |
| **Pipeline** | paired_full（[03_pipelines/01_paired_full.md](01_paired_full.md)） |
| **Mode variant** | paired-pure + methyl+context |
| **Baseline F1** | **0.9650**（ClairS 原始輸出，HC BED filtered） |
| **Post-ISM F1** | **0.9762** |
| **ΔF1** | **+0.0112** |
| **95% CI** | **[+0.0044, +0.0188]** |
| **Truth set** | SEQC2 high-confidence v1.2.1（39,447 variants）|
| **Benchmark 協議** | Internal `scripts/pipeline/steps/03_filter_analysis.py:229-234` + external som.py validation |
| **Method** | 60+ features × 748K regions, RF/XGBoost |
| **資料來源** | `docs/experiments/validated/2026/01/`（Phase 1A F1 優化）|
| **Data tier** | 🟢 A - Canonical（論文標準） |
| **Lock 日期** | 2026-04-17 |
| **Lock 理由** | Paired mode 已有 ClairS 內部 filter → 甲基化加持僅能篩出極少量殘留 FP → 上限小但統計顯著 |
| **Evidence rating** | ⭐5 |
| **Status** | `validated` / `locked` |
| **Last reviewed** | 2026-04-19（`10_Research_Chain_Registry.md` CL-002）|
| **Chain Registry ID** | CL-002 |

---

### CL-003 衍生：TO ΔF1 ⭐5 NEGATIVE

| 欄位 | 值 |
|------|-----|
| **Claim** | ClairS-TO + LongPhase-TO + ISM 的 ΔF1 為負（非 canonical filter） |
| **Sample** | 7 樣本平均（ONT 多平台混合）|
| **Pipeline** | tumor_only（[03_pipelines/03_tumor_only.md](03_tumor_only.md)）|
| **ΔF1** | **-0.0206**（負增益） |
| **95% CI** | 未 locked（跨樣本差異較大；H2009 ceiling effect 等 caveat） |
| **Benchmark 協議** | 同 CL-002；som.py 為主 |
| **Method** | 60+ 特徵 × 748K regions RF/XGBoost；HP-free dual path AUC=0.564；Fine-Pairwise 6 距離全無效 |
| **Lock 日期** | 2026-04-17（與 CL-002 一同凍結） |
| **Lock 理由** | 特徵空間耗盡（`project_beyond_auc_exhaustion_confirmed`）；TO 模式 HP tag self-phasing 導致甲基化訊號被偏差污染 |
| **Evidence rating** | ⭐5 |
| **Status** | `concluded` / NEGATIVE |
| **Chain Registry ID** | CL-003 |

---

### CL-005：V5 Somatic Fallback F1 ⭐4

| 欄位 | 值 |
|------|-----|
| **Claim** | V5 somatic fallback 機制將 ambiguous read 比例從 17.5% 降至 8.0% |
| **Sample** | HCC1395 |
| **Pipeline** | LongPhase-S V5（非 ISM 本身；上游 haplotag） |
| **F1（絕對值）** | **0.7154**（SEQC2 HCC1395） |
| **AMB%** | 17.5% → 8.0%（V3-Fixed → V5）|
| **Clean blocks** | 95% |
| **Method** | V3-fixed vs V5 per-sample F1 對照 |
| **Data tier** | 🟢 A（當前 canonical haplotag 版本） |
| **取代** | V3-Fixed 已被 V5 取代 |
| **Chain Registry ID** | CL-005 |
| **MEMORY** | `project_v5_somatic_fallback_verification` |

**注意**：此 F1=0.7154 是 **LongPhase-S 階段**的 F1（haplotag 正確率 × variant F1 混合），**非** ISM pipeline 整體 F1。勿與 CL-002 的 0.9762 混淆。

---

## 📐 不同 pipeline 變體的 F1 對照

（所有數字 anchored 到 CL-002/003/005，見上表）

| Pipeline | Baseline F1 | Post-ISM F1 | ΔF1 | CI | Status |
|----------|-------------|-------------|------|-----|--------|
| paired_full | 0.9650 | 0.9762 | **+0.0112** | [+0.0044, +0.0188] | 🟢 locked |
| paired_pileup | ~類似 paired_full | ~類似 paired_full | ~+0.01 範圍 | 無正式 CI | 🟡 未正式 lock |
| tumor_only (TO) | — | — | **-0.0206** | 無正式 CI | 🔴 NEGATIVE locked |

**為何 paired_pileup 無正式 CI**：Phase 1A F1 優化主表以 paired_full 為主；paired_pileup 作為模型變體對照，未做獨立 bootstrap CI。

---

## 🔬 運行條件完整規格（Reproducibility）

任何人想重現 +0.0112 的 Δ F1，必須使用以下組合：

### 輸入
- **Tumor BAM**：`/big8_disk/data/HCC1395/ONT_5khz_simplex_5mCG_5hmCG/HCC1395.bam`
- **Normal BAM**：`/big8_disk/data/HCC1395/ONT_5khz_simplex_5mCG_5hmCG/HCC1395BL.bam`
- **Reference**：`/big8_disk/ref/GRCh38_no_alt_analysis_set.fasta`
- **Tumor haplotagged BAM**：V5 Somatic Fallback 輸出
- **Truth VCF**：`/big8_disk/data/HCC1395/SEQC2/high-confidence_sSNV_in_HC_regions_v1.2.1.vcf.gz`
- **HC BED**：`/big8_disk/data/HCC1395/SEQC2/high-confidence_regions_v1.2.bed`

### ISM 執行參數
- Window size：`-w 5000`（Phase 1A canonical）
- Distance metric：`--distance-metric BERNOULLI`
- Min common coverage：`--min-common-coverage 3`
- Threads：`-j 120`（server-dependent）
- 其他：預設值（詳見 [04_parameters/01_cli_arguments.md](../04_parameters/01_cli_arguments.md)）

### 下游 ML Filter
- **Feature set**：60+ features（主要來自 significance_summary.csv 59 欄 + Python 衍生）
- **Model**：RF / XGBoost ensemble
- **Training set**：748K regions（7 樣本 × 2 modes 聚合）
- **Validation**：external，paired-pure + methyl+context split
- **Code**：`scripts/pipeline/steps/03_filter_analysis.py:229-234` 計算 F1

---

## ⚠️ 處理限制 / 已知 caveats

1. **僅 HCC1395 驗證**：CL-002 的 +0.0112 為 HCC1395 單樣本 external validation；跨樣本可能有變異
2. **HC BED 依賴**：F1 計算在 SEQC2 HC BED region 內；全基因體 F1 會較低
3. **ClairS 版本鎖定**：baseline 0.9650 對應特定 ClairS 版本；換 caller 版本需重算
4. **paired_pileup 僅近似**：無正式 CI，勿引用為主結論
5. **TO 為 NEGATIVE**：TO ΔF1=-0.0206 代表 ISM **不應**用於 TO mode 作 filter；僅 characterization 用途
6. **Phase 1A 凍結，不再追求**：ISM F1 提升已探索殆盡；重心轉向 Phase 2 characterization

---

## 🔍 驗證命令（如何重算）

```bash
# 1. 找 HCC1395 canonical paired_full run
CANON=$(ls -td /big7_disk/liaoyoyo2001/big7_disk_output/canonical/HCC1395/paired_full/*/ | head -1)
echo "Canonical: $CANON"

# 2. 取得 benchmark_comparison.tsv
cat ${CANON}/benchmark_comparison.tsv | head

# 3. 跨 sample 聚合
python3 scripts/analysis/compare_benchmark_f1.py \
  --run-dir /big7_disk/liaoyoyo2001/big7_disk_output/canonical/ \
  --output-tsv /tmp/all_f1.tsv

# 4. 檢查 Phase 1A external validation delta
grep -A 3 "paired-pure\|methyl+context" docs/experiments/validated/2026/01/*.md | head -20
```

**預期結果**：
- benchmark_comparison 顯示 `delta_f1 ≈ +0.0112`
- 跨 run 聚合後 HCC1395 paired_full 一致

---

## 🔗 Provenance 鏈（完整追溯）

```
[本頁 CL-002]
    │
    ├── docs/reports/research_landscape/10_Research_Chain_Registry.md#CL-002  (權威 chain)
    ├── docs/CURRENT_FOCUS.md line 57-59 (Phase 1A 已鎖定)
    ├── docs/reports/research_landscape/05_證據鏈總覽.md line 231 (external validation 表)
    ├── docs/reports/research_landscape/06_結論穩定性審查.md line 253 (S1 穩定性節點)
    ├── docs/experiments/validated/2026/01/ (Phase 1A 完整實驗紀錄)
    └── scripts/pipeline/steps/03_filter_analysis.py:229-234 (F1 公式)
```

---

## 📝 變動歷史

| 日期 | 事件 | 責任 |
|------|------|------|
| 2026-04-17 | Phase 1A locked，+0.0112 正式鎖定 | 研究 PI |
| 2026-04-19 | CL-002 加入 Research Chain Registry | 研究 PI |
| 2026-04-22 | KB SoT 頁面建立（本頁） | KB 維護者 |

**下次重檢**：Phase 2 A+D 完成後（預計 2026-06）— 若 Phase 2 有新 F1 加入需擴充本表

---

## 🔄 更新本頁的流程（重要）

若任何 ΔF1 數字有變動：

1. **先更新** `docs/CURRENT_FOCUS.md` 權威
2. **接著更新** `docs/reports/research_landscape/10_Research_Chain_Registry.md`
3. **最後更新本頁** + `CHANGELOG.md`
4. **不需** 改其他 21 處引用頁面（它們連結到本頁，自動取得最新）
5. 跑 3 驗證腳本 + freshness 更新
6. 若數字變動 > 0.001 → multi-agent 驗證（見 `00_governance/07_new_info_protocol.md`）

---

## 相關

- 完整 pipeline 規格：[01_paired_full.md](01_paired_full.md)
- F1 計算公式：[../08_truth_and_benchmark/02_f1_calculation.md](../08_truth_and_benchmark/02_f1_calculation.md)
- Truth set：[../08_truth_and_benchmark/01_truth_set_registry.md](../08_truth_and_benchmark/01_truth_set_registry.md)
- Current status：[../10_research_status/01_current_focus_snapshot.md](../10_research_status/01_current_focus_snapshot.md)
- New info protocol：[../00_governance/07_new_info_protocol.md](../00_governance/07_new_info_protocol.md)
- 權威 Chain Registry：[../../docs/reports/research_landscape/10_Research_Chain_Registry.md](../../docs/reports/research_landscape/10_Research_Chain_Registry.md)
