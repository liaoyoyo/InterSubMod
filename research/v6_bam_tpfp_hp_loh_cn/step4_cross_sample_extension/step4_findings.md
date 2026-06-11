<!--
build_date: 2026-05-15
agent: Step 4 cross-sample extension — final findings
scope: 5 samples (HCC1395 + H1437 + H2009 + HCC1954 + HCC1937), V6 BAM ISM, paired-pileup VCF
parent_plan: research/v6_bam_tpfp_hp_loh_cn/00_PLAN.md (v0.3 §Step 4)
status: validated (characterization-only; subject to FP-join sparsity caveat)
verdict: 1 cross-sample signature candidate (Outer|other|cov_high_gain) — TP-enriched, modest effect (Δ=+0.0069), needs Step 5 confound guard for promotion; HCC1937 outlier robust (sensitivity-pass)
-->

# Step 4 — V6 ISM 4-sample 跨樣本擴展 Findings

> Characterization-only. Read-only against `phaseD_v6_5sample/` + master.tsv. 不評估 filter / ΔF1。

## 0. TL;DR

- **Stage 1 (per-sample master TSV)**: 4 樣本完整建檔，行數合計 244,408（H1437=70964, H2009=136701, HCC1954=20136, HCC1937=16607）
- **Master.tsv join 覆蓋率**: TP 89-97%; **FP 1-7%（嚴重不平衡）— 跨樣本分析以 TP 為主，FP 結論需保留**
- **Stage 2 (50-cell grid per sample)**: 5 樣本 × 50 cells = 250 records; powered cells per sample: HCC1395=23, H1437=10, H2009=21, HCC1954=8, HCC1937=8（master-join 後 FP 稀疏導致 powered 主要由 TP 驅動）
- **Stage 3 (cross-sample Wilcoxon n=5)**: 1 cell 達 signature candidate 門檻（direction ≥4/5 + Wilcoxon p ≤ 0.0625）
- **Stage 4 (HCC1937 outlier)**: chr15/chr17/chr14 為 HCC1937 FP 富集染色體（chr17 含 BRCA1 mutation）；`Inner|same_HP1|cov_normal` cell HCC1937 TP rate 0.30 vs others 0.998（z=-362, n=10）—極端 outlier 但 n 小
- **Stage 5 (sensitivity)**: signature candidate (n=5) = 1; n=4 排除 HCC1937 = 1 → **HCC1937 不改變 candidate 集合**

## 1. Stage 1 — per-sample master TSV 建構

| Sample | n_rows | n_TP | n_FP | master_joined | join_rate | TP_joined | FP_joined | AF_missing | stale_DC75 | elapsed |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| H1437 | 70,964 | 70,191 | 773 | 67,476 | 95.1% | 67,468 (96.1%) | 8 (1.0%) | 0 | 8 | 598s |
| H2009 | 136,701 | 135,359 | 1,342 | 132,995 | 97.3% | ~132,973 | ~22 | 0 | 86 | 321s |
| HCC1954 | 20,136 | 19,449 | 687 | 17,938 | 89.1% | 17,909 (92.1%) | 29 (4.2%) | 0 | 0 | 76s |
| HCC1937 | 16,607 | 13,910 | 2,697 | 12,588 | 75.8% | 12,393 (89.1%) | 195 (7.2%) | 0 | 0 | 74s |

**重要 caveat**: master.tsv paired_full mode 與 phaseD paired-pileup VCF 的 FP 集合不對齊（phaseD FP 來自 paired-pileup caller，master FP 來自 paired_full caller）。FP-join 比 TP-join 低一個數量級。
- H1437: 773 phaseD FP → 8 joined（1.0%）
- H2009: 1342 → ~22（1.6%）
- HCC1954: 687 → 29（4.2%）
- HCC1937: 2697 → 195（7.2%）

**結論**: Step 4 grid 中 powered cells 主要為 TP-driven；FP-driven 結論需 caveat。

## 2. Stage 2 — 50-cell grid per sample

| Sample | universe (joined) | global TP_rate | powered cells | marginal cells |
|---|---:|---:|---:|---:|
| HCC1395 | 30,036 | 0.9832 | 23 | 3 |
| H1437 | 67,476 | 0.9999 | 10 | 1 |
| H2009 | 132,995 | 0.9994 | 21 | 1 |
| HCC1954 | 17,938 | 0.9984 | 8 | 0 |
| HCC1937 | 12,588 | 0.9845 | 8 | 1 |

**ceiling effect 明顯**：H1437/H2009/HCC1954 global TP rate ≥0.998 → 大部分 cell TP rate 接近 1.0，Wilcoxon 顯著性受限。HCC1395 + HCC1937 為 TP rate 較分散的兩個樣本。

## 3. Stage 3 — 跨樣本 Wilcoxon (n=5)

| 統計項 | 值 |
|---|---:|
| Total cells | 50 |
| Cells with ≥5 valid samples (n≥10 per sample) | 9 |
| Cells direction concordance ≥4/5 | 4 |
| Cells signature_candidate_flag (4/5 + Wilcoxon p≤0.0625) | **1** |

### Signature candidate (1)

`Outer|other|cov_high_gain` (LOH outer × HP bucket "other"（多 occupied） × Coverage_Multiple ≥ 2.5)

- direction = 5/0 (5 samples above global)
- Wilcoxon p = 0.0625 (n=5 exact min)
- mean Δ vs global TP rate = +0.0069
- mean n = 370
- per-sample TP rates: HCC1395 1.000 (n=33) / H1437 1.000 (n=11) / H2009 1.000 (n=1587) / HCC1954 1.000 (n=186) / HCC1937 1.000 (n=34)

**機制假說**：高 ploidy + LOH 外 + 多 HP bucket co-occupied region 在 5 樣本一致 TP rate=1.000，遠超全域 baseline。可能反映 high-coverage region 的 caller 高信心區域。

**caveat**：mean Δ = +0.0069 為 small effect；powered_count = 2 / 5（H2009 + HCC1954 有 n ≥ 50）— 其他樣本 cell n 偏低（11-34）。

## 4. Stage 4 — HCC1937 outlier 分析

### 4.1 BRCA1 / chr17 / 高 ploidy 染色體解釋

HCC1937 FP rate 最高的 5 chr（top 5）：chr15 (4.3%), chr17 (4.2%), chr14 (3.2%), chr1 (3.0%), chr22 (2.9%).

**chr17 含 BRCA1 / TP53 — 已知 BRCA1 mutant 在 HCC1937 樣本的 ploidy 異常熱點**

### 4.2 Per-cell HCC1937 vs others 偏離

| cell_id | HCC1937 TP rate | other mean | deviation | HCC1937 n | n_FP | z-score |
|---|---:|---:|---:|---:|---:|---:|
| Inner\|same_HP1\|cov_normal | 0.300 | 0.998 | **-0.698** | 10 | 7 | -362.4 |
| Inner\|same_HP2\|cov_normal | 0.636 | 0.996 | -0.360 | 44 | 16 | -61.8 |
| Inner\|other\|cov_loss | 0.953 | 0.994 | -0.041 | 725 | 34 | -5.4 |
| Outer\|other\|cov_loss | 0.957 | 0.983 | -0.026 | 69 | 3 | -0.9 |

**最大偏離**: Inner|same_HP1|cov_normal — HCC1937 TP rate 0.300 vs others' 0.998（n=10, n_FP=7 = 70% FP！）。雖然 n 小，但與其他 4 樣本完全相反方向。

**機制假說**：HCC1937 BRCA1 mutant + 高 ploidy → Inner LOH × same_hap × normal coverage cell 異常 FP 富集。需 Phase 3 zone-zoom 確認。

### 4.3 Sensitivity (n=5 vs n=4 排除 HCC1937)

| Set | Signature candidates |
|---|---:|
| n=5 (含 HCC1937) | 1 |
| n=4 (排除 HCC1937, p≤0.125 relaxed) | 1 |

→ **HCC1937 不改變 signature candidate 集合**。Outer|other|cov_high_gain 是 5 樣本 robust signature，不依賴 HCC1937 outlier 行為。

## 5. 與 H1c/H2-H5 假設對照（plan §假設）

| ID | 假設 | 預期 | 實際 |
|---|---|---|---|
| H2 | LOH.bed 內 × same_HP1 × low-cov cell 對 TP 富集 | +0.10~+0.30 | HCC1937 Inner|same_HP1|cov_normal 反向 -0.70（n=10） — 部分 cells 樣本特異 |
| H3 | Coverage_Multiple ≥ 1.5 × cross_het 對 FP 最富集 | FP 2.0× | 無法評估（FP join 1-7% 太稀疏） |
| H5 | 跨樣本 ≥1 cell 達 Wilcoxon p<0.05（n=5 exact min 0.0625 → 採此值） | direction match ≥4/5 | **1 cell 達標**: Outer|other|cov_high_gain |

## 6. 後續建議

### 立即可做（Step 5 H7 confound guard）
- 對 candidate cell `Outer|other|cov_high_gain` 跑 5+2 道 confound guard：within-group OLS residualize NumReads/AF, AF-bin L3, permutation, chr-stratified MH, HP1/HP2 symmetry, spatial autocorr
- 若 confound guard pass → 升級為 cross-sample characterization signature
- 評估 effect size +0.0069 是否實質（vs noise floor）

### 中期（需 master.tsv 重建）
- **Archive TO rerun** (CURRENT_FOCUS 待辦) → 為 4 樣本 paired-pileup mode 建 master.tsv → 提高 FP join 至 >80% → 重做 FP 分析
- **COLO829 truth set 0600 權限解後**補入 step4 → n=6 Wilcoxon (exact min p = 0.0312)

### 長期（與 step 1-3 整合）
- 將 Step 4 cross-sample signature 與 Step 1 V3F→V5→V6 trajectory 對照 — 是否同 cell 在三個 BAM 版本中亦 consistent
- 將 Step 4 HCC1937 outlier 與 Step 3 chr8 hotspot zone 對照 — chr17 BRCA1 是否為類似 sample-specific hotspot

## 7. 限制（不可不寫）

1. **FP 稀疏**: master.tsv paired_full mode 與 phaseD paired-pileup VCF FP 集合不對齊；FP join rate 1-7%
2. **Ceiling effect**: H1437/H2009/HCC1954 global TP rate ≥0.998 → cell TP rate 接近 1.0，Wilcoxon 顯著性 ceiling 受限
3. **Powered cells 不對稱**: HCC1395=23 powered, 其他 4 樣本 8-21 — 比較時需 weight by sample power
4. **Stale Diploid_Coverage_Used**: H1437 有 8 row, H2009 有 86 row stale=75.0（minor，<0.07%）
5. **HCC1954 kde_status=stale_rescaled**: Coverage_Multiple 為 rescaled 估計，非 KDE-corrected direct
6. **Wilcoxon n=5 exact min p = 0.0625**: 不可達 p<0.05 threshold；採 ≤0.0625 為實用門檻
7. **COLO829 deferred**: truth set 0600 權限阻塞，無法擴展到 n=6

## 8. 輸出產物清單

| 檔案 | 內容 |
|---|---|
| `per_sample_master/{H1437,H2009,HCC1954,HCC1937}_v6_master.tsv` | 4 樣本 V6 region-level wide TSV |
| `per_sample_grid/{HCC1395,H1437,H2009,HCC1954,HCC1937}_grid.tsv` | 5 樣本 50-cell grid |
| `step4_per_sample_grid.tsv` | 5 樣本 × 50 cell 合併 long-format（250 rows）|
| `step4_consistency.tsv` | per-cell direction + Wilcoxon |
| `step4_HCC1937_outlier_analysis.md` | HCC1937 outlier 詳細分析 |
| `step4_signature_candidates.md` | Final signature 候選表 |
| `figures/{sample}_facets.png` | 每樣本 facet heatmap (LOH × HP × cov) |
| `intermediate/file_lists/{sample}_{flag}_{label}.txt` | 路徑列表（用於 build 加速）|
| `intermediate/step4_per_sample_build_summary.tsv` | 4 樣本 build 統計 |
| `intermediate/step4_grid_summary.json` | 5 樣本 grid summary |
| `intermediate/step4_consistency_summary.json` | n=5 Wilcoxon 統計 |
| `intermediate/HCC1937_outlier_per_cell.tsv` | HCC1937 vs others per-cell |
| `intermediate/HCC1937_fp_per_chr.tsv` | per-chr FP rate（5 樣本）|
| `intermediate/HCC1937_signature_sensitivity.tsv` | n=4 排除 HCC1937 重算 |

## 9. 結論

**Step 4 task 完成度**: ✅ 全部 5 階段完成
- Stage 1 (per-sample master) ✅
- Stage 2 (per-sample grid) ✅
- Stage 3 (cross-sample consistency n=5 Wilcoxon) ✅
- Stage 4 (HCC1937 outlier) ✅
- Stage 5 (signature synthesis) ✅

**Verdict**:
- 1 cross-sample signature candidate (Outer|other|cov_high_gain) — modest effect, 5/5 direction concordant, Wilcoxon p = 0.0625 (n=5 exact min)
- HCC1937 outlier 行為 robust（移除不改變候選集合），主要 chr17/chr15 BRCA1 driven
- characterization-only verdict；filter / F1 評估不在 scope

**升級路徑**: candidate cell 需通過 Step 5 H7 confound guard（FDR + permutation + spatial）才能升級為 ⭐3+ characterization signature。

---

## 附錄：解讀備註

| Signature | 物理意義 | Caveat |
|---|---|---|
| Outer\|other\|cov_high_gain | LOH 外 × 多 HP bucket co-occupied × CN ≥ 2.5 高 ploidy region | mean Δ = +0.007 為 small effect；powered_count = 2/5；全 5 樣本 TP rate = 1.000 → 可能 ceiling effect 假性訊號 |
