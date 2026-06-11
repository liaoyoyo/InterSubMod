---
title: cis-ASM 特徵矩陣試用 + 每特徵合適度記錄（BRCA2 + 9 對照）
date: 2026-06-02
sample: HCC1395 paired_full（單樣本）
status: in_progress / trial
scope: 10 位點（BRCA2 canonical + TP/het/FP/LOH 對照）× 38 特徵
artifacts: scripts/35_feature_matrix_trial.py + genome_survey_v2/feature_matrix_trial.{json,tsv}
figure: figures/feature_matrix_trial.png
commit: 859c55a
---

# cis-ASM 特徵矩陣試用 + 每特徵合適度記錄

> **目的**：先「試用」~30 個特徵對 BRCA2 + 對照做多角度交叉驗證，**記錄每個特徵是否合適**，再決定建完整 pipeline。

## TL;DR

**多角度交叉驗證法成立且實用**：BRCA2 在 **~10 個獨立特徵同時點亮**（多角度一致 = 高可信）；對照組要嘛無 cis、要嘛 untestable。但試用也**暴露 3 件事**：(1) 核心 cis 特徵只 **4/10 可測**（testability 根本限制）；(2) `dbeta_5mC` 口徑不一致須修；(3) `kstar/crosstag` small-n 不穩須加 min-n gate。

## BRCA2 多角度交叉驗證（~10 角度一致）

| 角度 | 特徵 | BRCA2 值 | 支持 cis? |
|------|------|---------|:---:|
| 訊號 | dbeta_HP / percpg_maxabs | −0.122 / 0.917（一個 CpG 幾乎全翻）| ✅ |
| 訊號位置 | maxabs_dist_to_var | −242（上游 242bp）| ✅ |
| 均勻度 | entropy HP1→HP1-1 | 0.211→0.158（somatic 更 homogeneous）| ✅ |
| cis vs drift | d_cis / d_drift / p_cis | −0.142 / −0.022 / 0.0 | ✅ |
| 內聚 | sil_HP11 / delta_cohesion | 0.313 / +0.194 | ✅ |
| 覆蓋/power | n_HP11 / n_paired_cpg / power_class | 38 / 193 / DETECTED | ✅ |
| 等位比例 | vaf_obs | 0.202 | ✅ |
| 調控脈絡 | dist_to_tss / gene / cpg_island | **235bp / ZAR1L / shelf** | ✅⭐ |
| 機械機制 | mechanical_cis | **NEUTRAL（非 CpG 破壞 → 調控性 cis）** | ⭐ |
| 拷貝 | total_cn | 5.0（gain）| 脈絡 |
| 非 error | mapq / nm / softclip | 60 / 0.025 / 0.026（乾淨）| ✅ |

→ **多角度一致 = 高可信**；任一角度矛盾會被 flag。對照如 het chr1:181698147 cohesion 高達 0.598 **但 d_cis≈−0.006（無 cis）** → 自動正確排除（cohesion≠cis 再驗證）。

## 試用帶出的新發現（trial 的價值）

1. **BRCA2 訊號在 ZAR1L TSS 235bp 內、CpG shelf** → 強調控脈絡支持 cis 合理性。
2. **mechanical_cis=NEUTRAL** → BRCA2 的 cis **不是「突變破壞 CpG」機械機制，而是調控性**（更有意思但更難證 causation）。
3. mechanical_cis 真為正交維度：對照中 het(DYTN)=**DESTROYS_CpG**、FP=**CREATES_CpG** 被自動標出。
4. somatic entropy < germline entropy（更 homogeneous）= 與 cohesion 獨立印證子克隆內一致。

## ★ 每特徵合適度記錄（核心交付）

| 判定 | 特徵 | 可測 | 記錄 |
|------|------|:---:|------|
| **✅ 合適·always-on** | vaf_obs, total_cn, **dist_to_tss**, **mechanical_cis**, cpg_island, gene, n_HP{1,11,2}, power_class | 10/10 | 全位點可算、有資訊；dist_to_tss + mechanical_cis 是試用新增的高價值正交角度 |
| **✅ 合適但 testability-gated** | dbeta_HP, d_cis, d_drift, d_somatic, p_cis, sil_HP11, delta_cohesion, entropy_HP11, percpg_maxabs, maxabs_dist | **4/10** | 核心 cis 量，但需 HP1-1 + 足量覆蓋；**~60% 位點無法測 HP-cis（VAF 低 / HP1-1 不足）= 根本 testability 限制**，必須先 gate |
| **⚠ 需修** | dbeta_5mC | 6/10 | BAM window 口徑（−0.027）≠ Level-1 paired any-mod（−0.122）；**須對同一 CpG set 做 paired 5mC 重算**才能比 |
| **⚠ small-n 不穩** | kstar, crosstag_chi2_p | 10/10, 7/10 | TP_hiARI n_HP11=2 卻 kstar=4（無意義）；BRCA2 crosstag p=0.25 不顯著卻有強 cohesion；**須加 min-n gate（如 n≥15 才算）** |
| **🔻 本 set 低資訊** | frac_untag, mapq, nm_rate, softclip, strand_balance | 10/10 | 全乾淨/全平衡 → 當 **error/品質 guard** 有用，但此 set 無區辨變異；保留作 flag 不作排序主軸 |

## 結論（合不合適）

1. **方法合適**：多角度特徵矩陣能交叉印證 cis 判斷、能自動排除偽陽（高 cohesion 但無 cis 的 het）、能帶出新角度（TSS 距離、機械機制）。
2. **最大限制 = testability**：HP-cis 核心特徵只 4/10 可測 → 找更多位點時，**多數位點只能用 always-on 特徵 + ALLELE-axis**，HP-cis 確認是少數高覆蓋位點的奢侈品。
3. **要修 2 個**：dbeta_5mC 口徑統一、kstar/crosstag 加 min-n gate。
4. **可丟/降級**：frac_untag（此 set 無資訊）；mapq/nm/softclip 降為 guard。

## 下一步建議

- 把這份「合適度記錄」併進 pipeline 計劃的特徵 schema：**always-on 特徵當主幹（全位點）+ HP-cis 特徵當 testable 子集的深度確認**。
- 修 dbeta_5mC（paired 同 CpG）+ kstar/crosstag min-n gate，重跑 trial 確認。
- 之後 P0 地基用這個**已試用、已記錄合適度**的 schema 建，避免建一堆無用欄位。

## 修正後確認（2026-06-02）— 修 2 個特徵 + 重跑

**FIX 1 — dbeta_5mC paired + window-restricted（逼出 cross-extractor 驗證）**
- 根因：BAM 抽取**未限窗**，ONT 長 read 的 `modified_bases` 回報整條 read（10-30kb）所有 mod 位置，把 27,456bp 範圍 CpG 全算進去稀釋焦點 ±1kb 訊號（−0.036 artifact；Level-1 焦點 197 CpG span 1,911bp）。修正 = CpG 限變異點 ±1kb。
- 結果：**BAM-windowed dbeta_any = −0.122 = Level-1 dbeta_HP（diff=0.000）** → 兩個獨立抽取器（MSA + pysam）完全收斂 → **canonical −0.122 經 cross-extractor 驗證 bulletproof**。
- 分解：BRCA2 **5mC=−0.121 / 5hmC=−0.001 → ASM definitively 純 5mC**。
- 合適度：dbeta_5mC ⚠需修 → ✅ **合適**（paired + windowed + 跨抽取器一致）。

**FIX 2 — kstar/crosstag gate on n_som + n_som 透明化**
- FP 位點 n_som=0 → kstar/crosstag 正確變 None ✓。
- 「TP_hiARI n_HP11=2 卻 kstar=4」澄清：**n_som=52**（somatic reads 在 HP2-1 軸，非 HP1）→ kstar=4 有真實 reads 支持，非雜訊。原「無意義」是把 n_HP11 誤當 somatic 總數。
- 合適度：kstar ⚠不穩 → ✅ **合適**（gated n_som≥10 + n_som 透明）；定位 secondary 結構特徵。

**意外收穫**：「試用 → 修正 → 重跑」迴圈不只修特徵，還 **cross-validate 了 canonical −0.122（2 獨立抽取器）+ 確認純 5mC** → 強化核心結果可信度。schema 修正完成，可進 P0。
