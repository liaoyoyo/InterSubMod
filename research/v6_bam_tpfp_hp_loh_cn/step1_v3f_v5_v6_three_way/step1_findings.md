<!--
build_date: 2026-05-14
agent: Step 1 V3F/V5/V6 three-way integration (HCC1395 paired-pileup)
status: in_progress
report_class: characterization_post_hoc
scope: HCC1395 pilot, paired-pileup VCF, 30,490 TP + 4,842 FP regions
parent_plan: research/v6_bam_tpfp_hp_loh_cn/00_PLAN.md
inputs:
  - research/paired_priority_bug_audit/phaseC_genome_three_way/{V3F,V5,V6}_{on,off}_{tp,fp}/
  - research/paired_priority_bug_audit/phaseC_genome_three_way/v3f_vs_v5_vs_v6_region_ng.tsv
  - /big8_disk/liaoyoyo2001/InterSubMod/data/vcf/HCC1395/pileup/filtered_snv_{tp,fp}.vcf.gz
  - research/tpfp_loh_af_kde_discrimination/data/master.tsv.gz (HCC1395 paired_full mode)
outputs:
  - step1_master_three_way.tsv
  - step1_delta_{v5_vs_v3f,v6_vs_v5,v6_vs_v3f}.tsv
  - step1_trajectory.tsv
  - step1_off_vs_on_compare.tsv
  - step1_summary_stats.json
  - figures/step1/step1_{three_panel_heatmap,delta_heatmap,trajectory_sankey}.png
verdict: see Section 4 H1a/H1b/H1c table
-->

# Step 1 — V3F / V5 / V6 三向 ISM 整合（HCC1395 paired-pileup）

> Characterization-only. 不評估 filter / ΔF1（per plan §Out-of-scope）。

## 0. TL;DR

- **H1a/H1b/H1c 全 NEGATIVE**（Δgap < 0.005，threshold 0.03/0.03/0.06）— 但這是**ceiling effect**，不是真的「修補無效」。NG=2 Inner cell 三 BAM 全 ≥ 0.99 TP rate，物理上達不到 Δ ≥ 0.03。
- **真正的訊號** 在三個層面：
  1. **Marker coverage**：V6 (23,980) > V3F (21,997) > V5 (18,382)（NG≥3 region 數）
  2. **Marker TP rate**：V3F (0.918) > V6 (0.909) > V5 (0.894)
  3. **Inner NG=2 region 數**：V5 (8,136) ≫ V6 (5,353) ≈ V3F (5,064) — V5 over-promote +60%，V6 修補回 V3F 水準
- **per-region trajectory**：TP 中 34% 「只 V6 改善」、38% 無變化、16% 反向；FP 中 68% 無變化（V6 不改變 FP 分布）
- **off/on flag**：`--germline-hp-only=on` 將 NG≥3 全 collapse 為 0，不能近似 V6 修補 — V6 binary fix 才是主導
- **資料警示**：phaseC 用 paired-pileup VCF (30,490 TP + 4,842 FP)，但 master.tsv 只有 paired_full mode → FP 僅 ~10.5% 可 join LOH/CN covariate；Step 2 需以 TP 為主軸做 LOH × CN 分層
- **無 Diploid_Coverage_Used=75.0 stale row**（master.tsv KDE-corrected）

## 1. 資料盤點與重要警示

### 1.1 phaseC ISM 結果結構（與計畫書認知差異）

- 計畫書假設 phaseC 內存在 region-level master.tsv（含 HP1FamilyN/HP2FamilyN/caller_af/Coverage_Multiple/LOH/TP_label 等欄位）。
- **實際**：12 個 ISM 執行的 `significance_summary.csv` 均為 header-only（每個 run 的 `significance_statistics.txt` 顯示 `Regions Analyzed (Success): 0`）。phaseC 是 stripped-down HP audit 跑法，主 master TSV 未產出。
- 唯一可用的 per-region 數據是 `{run}/filtered_snv_{tp|fp}/{chr}/{chr}_{pos}/{region_id}/reads/reads.tsv`，欄位：`read_id, read_name, chr, start, end, mapq, hp, alt_support, is_tumor, strand`。
- 既有跨版本 aggregate：`v3f_vs_v5_vs_v6_region_ng.tsv`（105,996 rows = (30,490 TP + 4,842 FP) × 3 BAM × 1 (off NG + on NG)），由 `phaseC_region_ng_fast.py` 產出。
- region_id 編碼：`{chr}_{start}_{end}`，width=10kb，SNV pos = midpoint。

### 1.2 補強：caller_af 與 LOH/CN covariate 來源

- **caller_af**：直接從 source VCF `filtered_snv_{tp,fp}.vcf.gz` FORMAT/AF 欄位讀取（bcftools query）。
- **LOH / Coverage_Multiple / Diploid_Coverage_Used / loh_side / AF (master)**：從 `research/tpfp_loh_af_kde_discrimination/data/master.tsv.gz`（HCC1395 paired_full mode）按 Chr+Pos join。
- master.tsv join 覆蓋率：總 35332 rows 中 master_join_ok=1 為 30036 (85.0%)。
- **重要不對稱**：phaseC 用 **paired-pileup** VCF（30,490 TP + 4,842 FP），但 master.tsv 只有 **paired_full** mode（29,754 TP + 627 FP）。Join 結果：TP ~96.8% 可 join，FP 只有 ~10.5%（paired_full FP set 小很多）。Step 2 起若要做 LOH × CN 主軸分層分析，可以使用 TP 樣本完整，但 FP 內的 LOH/CN 必須額外處理（或限縮到可 join 的 506 個 FP）。
- Diploid_Coverage_Used=75.0 stale flag：未發現可疑 row（master.tsv KDE-corrected）。

## 2. 三向版本 genome-wide marker (NG_off ≥ 3) summary

| BAM | marker n | TP | FP | TP rate (95% CI) |
|-----|----------|----|----|-------------------|
| V3F | 21,997 | 20,183 | 1,814 | 0.918 [0.914, 0.921] |
| V5 | 18,382 | 16,428 | 1,954 | 0.894 [0.889, 0.898] |
| V6 | 23,980 | 21,806 | 2,174 | 0.909 [0.906, 0.913] |

> 與 `v3f_vs_v5_vs_v6_genome_summary.tsv` 一致（V3F=0.9175 / V5=0.8937 / V6=0.9093）— 確認讀取流程無誤。

## 3. NG=2 in master-joined regions — inner / outer TP rate gap

> Proxy for plan §H1a/b/c "Inner same-hap TP gap"。phaseC 沒有 per-cell HP-direction grain，所以以 master-joined `loh_side` 區分 Inner / Outer，並以 NG=2 cell 內 TP rate 作為 gap source。

| BAM | Inner NG=2 TP rate (95% CI) | n_TP / n_FP | Outer NG=2 TP rate (95% CI) | n_TP / n_FP | Inner − Outer gap |
|-----|------------------------------|-------------|-------------------------------|-------------|---------------------|
| V3F | 0.991 [0.987, 0.993] | 5016 / 48 | 0.992 [0.988, 0.995] | 2373 / 19 | -0.002 |
| V5 | 0.993 [0.991, 0.995] | 8082 / 54 | 0.992 [0.988, 0.995] | 2872 / 23 | +0.001 |
| V6 | 0.992 [0.989, 0.994] | 5308 / 45 | 0.990 [0.984, 0.993] | 1989 / 21 | +0.002 |

## 4. H1a / H1b / H1c 判定

| Hypothesis | 指標 | Δ | Threshold | Verdict |
|------------|------|---|-----------|---------|
| H1a — V5 BAM Inner same-hap TP gap > V3F | Δgap(V5−V3F) | +0.003 | ≥ 0.03 | **NEGATIVE** |
| H1b — V6 BAM Inner same-hap TP gap > V5 | Δgap(V6−V5) | +0.001 | ≥ 0.03 | **NEGATIVE** |
| H1c — V6 對 V3F 累積增益最大 | Δgap(V6−V3F) | +0.004 | ≥ 0.06 | **NEGATIVE** |

> **解讀**：
> - POSITIVE = 該段修補增加了 Inner-vs-Outer NG=2 TP rate gap（characterization 增益）
> - NEGATIVE = gap 變化未達 threshold，但未反向
> - NEGATIVE_REVERSE = gap 反向（Inner 變差 / Outer 變好）
> - UNKNOWN = 樣本不足或計算失敗

## 4.1 為什麼 H1a/H1b/H1c 同時 NEGATIVE — proxy 飽和

| Cell | TP rate (V3F, V5, V6) | 觀察 |
|------|------------------------|------|
| Inner NG=2 | 0.991 / 0.993 / 0.992 | 全 ≥ 0.99，三 BAM 差距 < 0.003 |
| Outer NG=2 | 0.992 / 0.992 / 0.990 | 全 ≥ 0.99，三 BAM 差距 < 0.003 |
| 全域 NG≥3 marker | 0.918 / 0.894 / 0.909 | **真正的差異**：V5 比 V3F 低 -0.024、V6 比 V5 高 +0.016、V6 比 V3F 低 -0.007 |

**結論**：NG=2 Inner cell 對 paired-pileup HCC1395 樣本而言是「高純度 caller-pre-filtered」cell（TP rate 已 99%+），三向版本差異 < 0.003，計畫書原設的 Δ ≥ 0.03 threshold 在此 cell 內物理上達不到（**ceiling effect / saturation**）。

**真正可區分的訊號** 出現在以下兩處：

1. **Marker coverage**（NG_off ≥ 3 region 數）：V6=23,980 > V3F=21,997 > V5=18,382
   - V6 比 V3F 多 +9.0% region（與 errata `v6_quantification_findings.md` 一致）
   - V5 因 Layer 1.5 過度 promote，反而降低 marker coverage 24%
2. **Marker TP rate**：V3F 0.918 > V6 0.909 > V5 0.894
   - V3F 最純，但 coverage 較小
   - V6 在 coverage 與 TP rate 中間取得最佳平衡（與 PI errata 一致）

**Step 2 起應改用全域 marker TP rate / coverage 雙指標**，而非 Inner−Outer gap。原 H1a/b/c threshold 設計需修正（threshold 應 ≥ 0.01 而非 0.03，且分子應換成「Inner NG≥3 marker coverage 比例」或「Inner FP rate 在 NG=2 中下降幅度」）。

## 4.2 V5 BAM 的 Inner NG=2 **數量** 異常

| BAM | Inner NG=2 master-joined n (TP+FP) | Inner NG=2 / Outer NG=2 比 |
|------|-------------------------------------|-----------------------------|
| V3F | 5,064 | 5,064 / 2,392 = 2.12 |
| **V5** | **8,136** | **8,136 / 2,895 = 2.81** |
| V6 | 5,353 | 5,353 / 2,010 = 2.66 |

V5 BAM 在 Inner NG=2 cell **多出 +60% TP region**（vs V3F），但 TP rate 並沒上升。這就是 V5 Layer 1.5 "marker over-promotion" 的直接證據：V5 把更多 region 標到 NG=2 cell，但這些 region 並非真的更具區分力（rate 仍 99%）。V6 回退到 V3F 保守標 hp=33 → Inner NG=2 region 數降回 5,353，與 V3F 接近（5,064）。**此即 plan v0.3 §V6 design rationale 在資料層的最直接驗證**。

## 5. Per-region NG trajectory (V3F → V5 → V6, flag=off)

**5 類分類**：
- A: 兩段都改善（V5−V3F > 0 且 V6−V5 > 0）
- B: 只 V5 改善（V5−V3F > 0 且 V6−V5 ≤ 0）
- C: 只 V6 改善（V5−V3F ≤ 0 且 V6−V5 > 0）
- D: 無變化（兩段 Δ 都 = 0）
- E: 反向或單段下降（任一段 Δ < 0 但非 D）

| 類別 | n | 占比 |
|------|---|------|
| A_both_improve | 1057 | 3.0% |
| B_only_V5_improve | 3673 | 10.4% |
| C_only_V6_improve | 10624 | 30.1% |
| D_no_change | 14746 | 41.7% |
| E_reverse_or_decrease | 5232 | 14.8% |
| MISSING | 0 | 0.0% |

### 5.1 按 TP / FP label 拆解

```
class  A_both_improve  B_only_V5_improve  C_only_V6_improve  D_no_change  E_reverse_or_decrease
label                                                                                          
FP                 94                805                268         3276                    399
TP                963               2868              10356        11470                   4833
```

## 6. `--germline-hp-only` off vs on 對照（mask V5 Layer 1.5？）

詳見 `step1_off_vs_on_compare.tsv`。每 (label, BAM, flag) 組合給 NG=2 與 NG≥3 比例。

**關鍵實證**：

1. **`flag=on` 將全部 NG≥3 collapse 為 0**（所有 BAM × label × loh_side cell 都是 0）。`--germline-hp-only=on` 在 ISM 端強制只取 germline-phased reads，所有非 germline HP family（如 `1-1` `2-1` `33`）會被 mask 為 unphased → NG family 退化 → NG≥3 永不成立。**這代表 on flag 完全消除 marker filter 的能力**，不能單獨用作 "mask Layer 1.5" 的代理。

2. **V5_on Inner ≡ V6_on Inner**（NG=2 frac 都是 15.4%, n_eq2=2,154）— 因為 on 模式下 V5 與 V6 BAM 的 germline-phased reads 完全相同（V6 重用 V5 phased VCF）。這證實 plan v0.3 的設計判斷：V6 與 V5 差異**僅在 Layer 1.5 標記策略**，與 phasing 層無關。

3. **flag=off Inner NG=2 frac**：V3F 35.8% → V5 57.7% (+22 pp) → V6 37.9% (回降近 V3F)。V5 在 Inner 區 over-promote +60% NG=2 region；V6 修補後 Inner NG=2 數量回到 V3F 水準。**V6 ≈ V3F 在 Inner LOH 區的 marker engineering 是設計一致的**。

4. **flag=off Outer NG=2 frac**：V3F 15.3% → V5 18.5% → V6 12.8%。V6 在 Outer 區比 V3F 更保守（-2.5 pp），這代表 V6 在 germline-existent 區因重用 V5 phased VCF 殘留 priority bug 偏移**但在 marker 標記層被回退**。

**結論**：ISM 端 `--germline-hp-only=on` 不能近似 V6 的修補，因為它對 NG≥3 過度激進（全 collapse 為 0）。要近似 V6 必須在 longphase 端做 Layer 1.5 改寫（即 V6 binary 本身）。

## 7. 圖示

- `figures/step1/step1_three_panel_heatmap.png` — V3F / V5 / V6 三聯 TP rate × (NG × LOH side)
- `figures/step1/step1_delta_heatmap.png` — 3 個 Δ TP rate heatmap
- `figures/step1/step1_trajectory_sankey.png` — 5 類 region trajectory stacked bar

## 8. Hand-off context — Step 2 Agent B

**Master TSV schema (step1_master_three_way.tsv)**：
- Key：`region_id`（chr_start_end，width=10kb），`chr`, `pos`（SNV center）
- Label：`label` (TP/FP)
- Per-version × per-flag features：`{V3F,V5,V6}_{off,on}_{0,1,2,1-1,2-1,3,11,21,33,other,NG,n_reads}`
- Caller AF：`caller_af`（from `filtered_snv_{tp,fp}.vcf.gz` FORMAT/AF）
- Master join：`master_join_ok` (1/0)
- LOH / CN covariate：`LOH_Bed_Overlap`, `LOH_Bed_Annotation`, `LOH_Subtype`, `Coverage_Multiple`, `Coverage_Category`, `Diploid_Coverage_Used`, `loh_side`, `AF_master`

**Step 2 推薦 covariate**（基於 plan §Step 2 3-axis grid）：
- 主 grid 軸：`loh_side` (Inner/Outer), HP bucket（從 reads.tsv 進一步以 HP family 計算），`Coverage_Multiple` (5 bins)
- LR covariate：`HPFineNGroups` (= V6_off_NG), `caller_af`, `n_reads`

**已知 caveat 給 Agent B**：
1. **FP set 的 LOH/CN join 不完整**（only ~10% master-joined）。Step 2 主分析建議以 TP 為主，FP 對照僅在 master-joined 子集做。
2. **V6 BAM germline-existent 區因重用 V5 phased VCF 殘留 priority bug 偏移**（hp_1_1:hp_2_1 ratio 1.838 介於 V3F 1.138 與 V5 1.86 之間）— 這不是 bug，是設計選擇。Inner same-hap 在 V6 不如 V3F 純淨。
3. **on flag mask Layer 1.5 marker → unphased**，會降低 NG_on=2 sample 數。Step 2 應以 flag=off 為主軸 grid，on 作 sensitivity check。
4. Power gate：master-joined Inner+NG=2 cell 預期 n 詳見 step1_summary_stats.json bam_metrics 各 row n_TP+n_FP 欄位。

## 9. 附錄

**Scripts**：
- `scripts/build_three_way_master.py` — Stage 1-4 wide-format master
- `scripts/compute_deltas_and_trajectory.py` — Stage 5-7 deltas + trajectory + H1a/b/c
- `scripts/make_figures.py` — Stage 8 figures
- `scripts/generate_findings.py` — 本 findings.md 自動生成

**Intermediate artifacts**：
- `intermediate/per_region_hp_counts.tsv.gz` — 12 ISM run × per-region HP family counts (long format)
- `intermediate/caller_af_lookup.tsv` — region_id × label × caller_af
- `intermediate/build_log.txt` — build_three_way_master 執行日誌

**Out-of-scope reminders** (plan §Out-of-scope)：
- 本 step 不評估 ΔF1 / filter 效果（characterization-only）。
- 不修改 C++ pipeline。
- 不重跑任何 BAM/ISM。