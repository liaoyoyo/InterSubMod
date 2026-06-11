<!--
建立時間: 2026-05-30
類型: LOH-constrained phasing 主軸 Grade A 升級 — n=7 cross-sample Wilcoxon
任務類型: B Comprehensive validation (全 7 樣本)
樣本 scope: 7 sample-runs = 6 distinct cell lines (HCC1395/HCC1937/HCC1954/H2009/H1437/COLO829) + HCC1395_DORADO replicate
框架: Verdict-Pyramid
ledger: 20260530_loh_phasing_n7_cross_sample_wilcoxon
-->

# LOH-constrained phasing — n=7 跨樣本 Wilcoxon（Grade A 要求 #1）

## TL;DR（Verdict）

**NG=2 Inner same_HP1 vs Outer cross_het 的 TP-rate gap 在全 7 樣本一致正向（7/7），Wilcoxon W=28 p=0.0078（exact，= 1/2⁷ 理論最小）。** 加入 COLO829 後 cross-sample 軸達 Grade-A strength（從 n=6 W=21 p=0.0156 強化）。R-SELFREF 自參考軸仍部分（X3+V6C 機制已證，全 7-sample flag-on 負控待重跑）→ 整體 **Grade B+ 接近 A**。

## 1. 背景

LOH-constrained phasing（TO mode）：LOH 區物理上只保留單 haplotype，somatic SNV 必然形成 same-hap ref/alt 子族（NG=2 same-hap，TP-rich）；非 LOH 區 NG=2 多為 cross-hap germline-somatic 不可分（FP 源）。機制純 phasing（無需 methylation 實驗）。原 n=6 Grade B（W=21 p=0.0156）。

## 2. 方法

- 每樣本 gap = Inner same_HP1 TP rate − Outer cross_het TP rate（Extreme AF 子集，NG=2）
- 口徑同 obs18（`Potential_LOH` 欄定 Inner/Outer；4 HPFineN bucket 欄定 combo）
- Wilcoxon signed-rank（exact, one-sided greater）+ bootstrap median CI（seed 20260423）
- **COLO829 來源**：`big7_disk_output/synthesis/research_rounds/20260423_colo829_to_pilot/step05_intersubmod/intersubmod_{tp,fp}/significance_summary.csv`（2026-04-23 真實 KDE-corrected TO ISM，run.log "[TO pilot] Completed"）— 親自驗證存在（先前 memory `project_to_cross_sample_archive_data_exists` 與 1 個 scout 誤判 blocked，源於 2026-03-17 dry-run stub）
- **provenance fix**：原 B1 Wilcoxon 計算腳本從未 commit → 本 session 寫 `obs18b_wilcoxon_ng2_gap.py`（重現 B1 reproduces=True）+ `obs18c_add_colo829_n7.py`

## 3. 結果

| Sample | Inner same_HP1 TP | Outer cross_het TP | Gap |
|--------|:---:|:---:|:---:|
| HCC1937 | 0.76 | 0.24 | **+0.52** |
| HCC1395 | 0.96 | 0.50 | **+0.46** |
| HCC1395_DORADO | 0.94 | 0.55 | +0.39 |
| HCC1954 | 0.43 | 0.08 | +0.34 |
| H1437 | 0.92 | 0.69 | +0.23 |
| **COLO829** | **0.724** | **0.620** | **+0.104** |
| H2009 | 0.93 | 0.88 | +0.05（baseline 飽和）|

| 統計 | n=6 (Grade B) | **n=7 (本日)** |
|---|---|---|
| Wilcoxon W | 21 | **28** |
| p-value (exact, greater) | 0.015625 | **0.0078125** |
| 正向 | 6/6 | **7/7** |
| median gap | 0.365 | 0.345 |
| bootstrap median CI | [0.140, 0.491] | [0.104, 0.459] |

## 4. 結論判定：POSITIVE（cross-sample 軸 Grade-A strength）

- **Grade A 要求 #1（n=7 cross-sample Wilcoxon）✅ 達成** — 7/7 一致，p=0.0078。
- **Grade A 要求 #2（R-SELFREF 全 7-sample flag-on 負控）⚠ 部分** — Inner same_HP1 用 HPFineN_HP1S（HP1-1 count，C++ 標 circular dependency）；X3（HCC1395 flag-on/off 2026-04-24 demote HP1-1→germline 使 bucket 塌陷）+ V6C 4-sample chr19 已證機制；全 7-sample full-genome flag-on **需 C++ 重跑（~25-50hr，archive 無現成 flag-on 輸出）尚未做**。

## 5. Limitations

- COLO829 platform asymmetry：5kHz/DORADO ONT vs COLO829 ONT_PAO；truth set NYGC vs SEQC2（manifest H4 confound 已 flag）
- HCC1954 outlier：Inner same_HP1 TP 僅 0.43（baseline TP 全 TO 最低 + Potential_LOH polyploid 可靠性存疑）
- H2009 baseline TP 飽和（gap +0.05，phasing gap 空間被吃光）
- n=7 = 6 distinct cell lines + 1 platform replicate（HCC1395×2）

## 6. Next steps

1. **R-SELFREF 全 7-sample full-genome flag-on 負控**（~25-50hr C++ 重跑）→ 升 full Grade A 的唯一剩餘 gate
2. `/run-evaluator` tier 判定（cross-sample L2 已 A-strength；L3 mechanism + L4 orthogonal 對照）
3. publication：citation-verification + prior-art（TumorLens/Wakhan/SAVANA）+ framework coverage limitation

**Provenance**: ledger `20260530_loh_phasing_n7_cross_sample_wilcoxon`；data `obs18b_wilcoxon_ng2_gap_n7.json`；scripts `obs18b_wilcoxon_ng2_gap.py` + `obs18c_add_colo829_n7.py`；git HEAD 274f152。所有數字本 session 親自從 significance_summary.csv 重算。
