<!--
建立時間: 2026-04-23
用途: PPT 內嵌圖片對照表（slide ↔ figure file ↔ 原始來源）
備註: figures/ 下均為複本；原始圖保留於 repo 各處以便重跑
-->

# Figures Index · v2（2026-04-23）

所有 PPT 使用的圖片均同時保存於 `figures/` 子資料夾（自包含）；下表也記錄原始路徑以便追溯重跑。

## 使用對照（Slide → Figure）

| Slide | 標題 | Figure (本地 `figures/`) | 原始路徑 |
|:-----:|------|--------------------------|----------|
| **S1** | Cover | `s1s7_heatmap_tp_rate.png` | `docs/experiments/in_progress/2026/04/figures/20260423_phase3_synthesis/s1s7_heatmap_tp_rate.png` |
| **S5** | LOH × AF × CN 初步觀察 ⭐ | `fig_s5_loh_inner_outer_af_cn_per_sample.png` | `docs/experiments/in_progress/2026/04/figures/20260423_s5_loh_af_cn_scatter/fig_s5_loh_inner_outer_af_cn_per_sample.png` |
| **S8** | KDE 方法 + 驗收 | `fig3_two_pass_architecture.png` | `docs/experiments/in_progress/2026/04/figures/20260420_kde_fix_acceptance/fig3_two_pass_architecture.png` |
| **S9** | 各樣本新 CN 狀況 (左) | `fig1_covm_density_per_sample.png` | `docs/experiments/in_progress/2026/04/figures/20260421_ng_kde_rescaled/fig1_covm_density_per_sample.png` |
| **S9** | 各樣本新 CN 狀況 (右) | `fig2_category_reclassification.png` | `docs/experiments/in_progress/2026/04/figures/20260420_kde_fix_acceptance/fig2_category_reclassification.png` |
| **S10** | Phase 3 跨樣本 TP rate heatmap | `s1s7_heatmap_tp_rate.png` | 同上 S1 |
| **S11** | Inner × NG=2 ≥93% same-hap | `obs18_NG2_composition_proportion.png` | `research/tpfp_loh_af_kde_discrimination/figures/new/obs18/obs18_NG2_composition_proportion.png` |
| **S15** | 跨樣本 scheme heatmap (左 TP rate) | `s1s7_heatmap_tp_rate.png` | 同上 |
| **S15** | 跨樣本 scheme heatmap (右 fold) | `s1s7_heatmap_fold.png` | `docs/experiments/in_progress/2026/04/figures/20260423_phase3_synthesis/s1s7_heatmap_fold.png` |
| **S16** | S3 最純 + LOH_Noise 保留 | `fig_v2_3_filter_scheme_bar.png` | `docs/experiments/in_progress/2026/04/figures/20260421_ng_kde_rescaled/fig_v2_3_filter_scheme_bar.png` |

其他 slides（S2, S3, S4, S6, S7, S12, S13, S14, S17-S26）均為純文字 / 圖表以 shape 手繪，不使用外部圖片。

## 圖片清單 (`figures/`, 9 張)

| 檔名 | 大小 | 類型 | 使用頁 |
|------|-----:|------|:------:|
| `s1s7_heatmap_tp_rate.png` | 178 KB | Phase 3 跨樣本 TP rate heatmap | S1, S10, S15 |
| `s1s7_heatmap_fold.png` | 167 KB | Phase 3 跨樣本 fold-improvement heatmap | S15 |
| `fig_s5_loh_inner_outer_af_cn_per_sample.png` | 644 KB | **S5 新圖**：6 樣本 × (Inner/Outer) × TP/FP 散點 + AF 三區 | S5 |
| `fig3_two_pass_architecture.png` | 176 KB | KDE 雙 Pass 架構 | S8 |
| `fig1_covm_density_per_sample.png` | 260 KB | Per-sample CovM 密度 | S9 |
| `fig2_category_reclassification.png` | 85 KB | Coverage_Category 重分類 | S9 |
| `fig_v2_1_to_tp_heatmap.png` | 140 KB | （舊 S5 用，已替換；保留備用） | — |
| `fig_v2_3_filter_scheme_bar.png` | 161 KB | S1-S7 scheme bar chart (HCC1395) | S16 |
| `obs18_NG2_composition_proportion.png` | 171 KB | NG=2 組成比例 6 樣本 Inner vs Outer | S11 |

## 重跑指令（如需更新）

### S5 新圖
```bash
python3 scripts/analysis/20260423_s5_loh_af_cn_scatter.py
cp docs/experiments/in_progress/2026/04/figures/20260423_s5_loh_af_cn_scatter/fig_s5_loh_inner_outer_af_cn_per_sample.png \
   docs/presentations/validated/2026/04/20260423_研究週報_LOH_constrained_phasing_pivot/figures/
```

### Phase 3 合成圖（s1s7_heatmap_*.png）
來源為 `scripts/analysis/20260423_phase3_synthesis.py`（B1-B7 合成）。

## 自包含規則

- **複製而非 symlink**：將來此 PPT 資料夾被複製到新環境時，`figures/` 內的圖片一併帶走。
- **原始保留**：上述原始路徑是 canonical 來源，有新結論需重跑時從原始路徑更新後再覆蓋 `figures/` 內副本。
- **`build_pptx.py` 路徑規則**：`FIG_*` 變數優先指向 `figures/`（本地），若本地無則 fallback 原始絕對路徑。
