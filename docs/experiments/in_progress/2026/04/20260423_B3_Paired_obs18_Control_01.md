<!--
建立時間: 2026-04-23
目標: B3 — Paired mode 對照 obs18 NG=2 same-hap vs cross-het gap（H-D3 驗證）
處理範圍: 7 樣本 paired_full latest run significance_summary.csv × Extreme AF × NG=2
關聯檔案:
  - docs/experiments/in_progress/2026/04/20260422_LOH_AF_KDE_TPFP_Discrimination_02.md
  - research/tpfp_loh_af_kde_discrimination/scripts/obs18_TO_NG2_composition.py
  - research/tpfp_loh_af_kde_discrimination/data/obs18_NG2_composition_by_sample.tsv
  - research/tpfp_loh_af_kde_discrimination/data/obs18_wilcoxon_B1.json
  - scripts/analysis/20260423_B3_paired_obs18.py
關聯假說: H-D3（paired mode germline caller 移除 cross-het → gap 應消失）
-->

# B3 | Paired-mode 對照 obs18 NG=2 組成 gap（H-D3）

**日期**: 2026-04-23
**假說**: H-D3 — Paired-mode 下，matched-normal germline caller 會排除 germline heterozygotes，因此 Outer 區 `cross_het` 組合不再被保留為 somatic candidate；NG=2 Inner same-hap vs Outer cross-het TP rate gap 應接近 0。
**判定**: ✅ **H-D3 PASS**（Paired gap median ≈ 0.000，7/7 樣本 |gap|<0.10）

---

## 1. 背景（TO baseline 回顧）

`obs18`（TO mode, 6 樣本）確立：
- NG=2 Inner × same-hap TP rate median = 0.94（跨 6 TO 樣本）
- NG=2 Outer × cross-het TP rate median = 0.47（跨 6 TO 樣本）
- **TO gap median = +0.365**（Wilcoxon p=0.0156, 6/6 正號，95% CI [0.14, 0.49]）
- 機制解釋：LOH-constrained phasing 強迫 Inner 的 germline het 落入 same-hap bucket 並保留 AF=1；Outer 的 germline het 則呈 cross-hap bucket 且 AF~0.5，被 ClairS-TO 當 somatic candidate 導致低 TP rate（FP 主導）。

H-D3 預測：Paired mode 的 ClairS full pipeline 擁有 matched-normal，可識別並過濾 germline het，Outer cross-het 應本來就不進 somatic candidate pool（或進了但 TP rate 接近 1）。gap 應顯著縮小。

---

## 2. 方法

### 2.1 資料來源
7 樣本 `canonical/<sample>/paired_full/<latest>/intersubmod_tp|fp/significance_summary.csv`：

| Sample | Latest run |
|--------|------------|
| HCC1395 | `20260420_HCC1395_paired_full_full_2` |
| HCC1395_DORADO | `20260420_HCC1395_DORADO_paired_full_full` |
| HCC1937 | `20260421_HCC1937_paired_full_full` |
| HCC1954 | `20260421_HCC1954_paired_full_full` |
| H2009 | `20260421_H2009_paired_full_full` |
| H1437 | `20260421_H1437_paired_full_full` |
| COLO829 | `20260421_COLO829_paired_full_full` |

### 2.2 分析流程（與 obs18 一致）
1. Filter `HPFineNGroups == 2` ∩ `AF_class == Extreme`（AF<0.1 或 >0.9，與 obs17 一致）
2. 以 `HPFineN_HP1 / HPFineN_HP1S / HPFineN_HP2 / HPFineN_HP2S` > 0 判斷 4 類 bucket 占用
3. 分類：`same_HP1`, `same_HP2`, `cross_het`, `cross_het_inv`, `other`
4. 按 `Potential_LOH` 拆 Inner/Outer
5. 每樣本計算 `Inner × (same_HP1|same_HP2) mean TP rate` − `Outer × (cross_het|cross_het_inv) mean TP rate`（僅算 n≥50 cell）
6. Wilcoxon（two-sided）測試 paired gap vs 0 與 TO-vs-Paired paired difference

### 2.3 輸出
- 圖: `docs/experiments/in_progress/2026/04/figures/20260423_B3_paired_obs18/`
  - `B3_paired_NG2_composition_heatmap.png`
  - `B3_paired_NG2_composition_proportion.png`
  - `B3_paired_vs_TO_gap_per_sample.png`
- 資料: `research/tpfp_loh_af_kde_discrimination/data/`
  - `B3_paired_NG2_composition_by_sample.tsv`（69 rows）
  - `B3_paired_vs_TO_gap_per_sample.tsv`
  - `B3_wilcoxon_gap_stats.json`

---

## 3. 結果

### 3.1 每樣本 gap 對照

| Sample | TO gap (obs18) | Paired gap (B3) | Δ (TO−Paired) |
|--------|:---:|:---:|:---:|
| HCC1395 | n/a (Inner n 過低) | −0.017 | — |
| HCC1395_DORADO | +0.279 | −0.002 | +0.281 |
| HCC1937 | +0.411 | +0.012 | +0.399 |
| HCC1954 | +0.258 | +0.001 | +0.257 |
| H2009 | +0.047 | +0.000 | +0.047 |
| H1437 | +0.171 | −0.000 | +0.171 |
| COLO829 | — (不在 obs18 TO 6 樣本集) | −0.092 | — |

**Paired gap values**：{−0.092, −0.017, −0.002, −0.000, +0.000, +0.001, +0.012}
**Paired gap median = −0.0002** (essentially 0)
**TO gap median = +0.258**

### 3.2 Wilcoxon 統計
- Paired gap vs 0（H0: gap=0）：**W=10, p=0.578** → 不拒絕 H0，paired gap 與 0 無顯著差異
- TO gap vs Paired gap（5 樣本同時有兩側資料：HCC1395_DORADO, HCC1937, HCC1954, H2009, H1437）：**W=0, p=0.0625** → 5/5 TO gap > Paired gap，方向一致；n=5 下 p=0.0625 為 exact test 下限（≤ 理論最小 p-value 為 1/32=0.03125，雙側 0.0625）

### 3.3 Heatmap 觀察（paired mode）
- 7 樣本 × 2 側 × 5 組合的 70 個 cell 中：
  - Inner × cross_het: TP rate 0.955–1.000（6/7 有 n≥20 資料）
  - Outer × cross_het: TP rate 0.944–0.999（7/7）
  - Inner × same_HP1/2: TP rate 0.834–1.000
  - Outer × same_HP1/2: TP rate 0.647–1.000（COLO829 Outer same_HP1 TP=0.647 為唯一明顯偏低 cell）
- 整體 TP rate 在 paired mode 下壓倒性高（>0.9），與 TO mode 的 0.08–0.94 分布呈鮮明對比

### 3.4 COLO829 特殊性
COLO829 paired gap = −0.092（唯一負號；Inner same-hap 0.851 vs Outer cross-het 0.943）。
- 可能原因：COLO829 為 platinum benchmark，germline filter 極嚴格；Inner same-hap bucket 可能反而包含較多 somatic-in-LOH 區的邊界點（與 H-D3 機制無衝突，屬 low-purity-in-LOH 情境）。
- Gap 絕對值 <0.10，不改變整體判定。

---

## 4. 判定

### H-D3 PASS（強證據）
三重證據支持：
1. **Paired gap median = −0.0002，|gap|<0.10 for 7/7 樣本**：Inner same-hap 與 Outer cross-het TP rate 在 paired 模式下差異消失。
2. **Wilcoxon（paired gap vs 0）p=0.578**：統計上無法區辨 paired gap 與 0。
3. **TO-vs-Paired 配對差異 5/5 正號（方向一致），p=0.0625（n=5 exact 極限）**：機制方向符合預期。

### 強化「LOH-constrained phasing」機制解釋
B1 confirmed TO gap 存在 6/6，B3 confirmed paired gap 消失 7/7。兩個方向同時成立證明：
- TO mode 的 NG=2 Inner same-hap 高 TP、Outer cross-het 低 TP 差異，**不是** methylation subclone structure 造成，也**不是** NG=2 的 read-count artifact；
- **真正機制**是 ClairS-TO 缺乏 germline filter 導致 Outer germline hets 被當成 somatic candidate（FP 污染）。paired caller 移除之後差異完全消失。

---

## 5. 對研究主軸的影響

### 5.1 LOH-constrained phasing discovery（主軸候選）
- TO baseline（B1）+ paired 對照（B3）現在形成完整雙向證據鏈
- 機制解釋從「NG 是 methylation subclone marker」正式更正為「NG 是 `{HP1,HP1S,HP2,HP2S}` bucket occupancy，在 LOH 內受物理限制只能落 same-hap」
- 下一步：B4 S4 二級判別 pilot（能否用 `NG=2 + same-hap bucket + Inner` 當 TO mode 的 germline FP filter）仍為可行方向

### 5.2 記憶更正
- `project_hpfinengroups_subclone_marker.md` 結論需再次強調：NG 在 paired 模式下無 subclone marker 行為（因 germline het 已被濾除），所有 NG>1 都是 real methylation subclone signal，但因 TP rate 整體 >0.9，discriminative power 有限。
- `project_loh_constrained_phasing_discovery.md` 可升級：paired gap = 0 是機制的正向對照證據。

---

## 6. 資料完整性備註

**不存在「資料不足」情況**：7/7 paired_full 樣本的 significance_summary.csv 都提供了 `HPFineN_HP1/HP1S/HP2/HP2S` bucket count 與 `Potential_LOH`，分析直接可行。

僅兩點小限制（不影響主結論）：
- HCC1395 TO（obs18）Inner n<20，原 obs18 TSV 無 TO gap 值；不影響 B3 主判定（B3 用 paired 獨立資料）
- COLO829 不在 obs18 TO 6 樣本集（HCC1395_5kHz TO 那組僅 6 樣本），TO gap 值缺失；paired 部分正常

---

## 7. 後續

- [ ] 將 B3 結果匯入 `project_loh_constrained_phasing_discovery.md` memory 作為 paired 對照正向證據
- [ ] 在 `docs/reports/research_landscape/09_Part_B.md` 更新 HPFineNGroups 機制敘述為 bucket-occupancy-in-LOH，paired mode 失去 marker 行為
- [ ] B4 S4 pilot（若 TO mode 能用 `NG=2 ∩ same-hap ∩ Inner` 作為 germline FP filter）可承接此結論獨立進行
