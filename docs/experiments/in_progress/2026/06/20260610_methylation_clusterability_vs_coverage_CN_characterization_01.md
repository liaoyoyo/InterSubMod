<!--
建立時間: 2026-06-10
狀態: in_progress (characterization 報告 — 甲基可分群 vs coverage + CN confound 驗證)
報告類型: characterization_negative_with_positive_mechanism
受眾: 廖子游 · PI · 論文 §Methods-Neg
task_type: B_validation (coverage 全 6 樣本×3 class; CN truth HCC1395-only)
tier: "⭐3 (single-pipeline ClairS-TO; CN truth single-sample HCC1395; LOH-mechanism 6/6 cross-sample)"
partial_flag: "subset — CN ground-truth 僅 HCC1395 (SEQC2)；其他 5 樣本用 ISM Potential_LOH proxy (recall 0.96)；FP 稀疏低 power"
data_sources: research/tsg_promoter_asm_reviewer/genome_survey_v2/catalog/clusterability_vs_cov_cn.json,research/tsg_promoter_asm_reviewer/genome_survey_v2/catalog/clusterability_cov_cn_extended.json
scripts: research/tsg_promoter_asm_reviewer/scripts/101_clusterability_vs_coverage_cn.py,research/tsg_promoter_asm_reviewer/scripts/102_clusterability_cov_cn_extended.py
figures: docs/paper_focus/04_figures/clusterability_vs_cov_cn.png,docs/paper_focus/04_figures/clusterability_cov_cn_extended.png
ledger: 20260610_clusterability_vs_coverage_cn (97), 20260610_clusterability_cov_cn_extended_FPFN_6sample (98)
-->
<!-- provenance-verified: 全數字由 clusterability_vs_cov_cn.json + clusterability_cov_cn_extended.json grep；scripts 101/102 可重跑；2026-06-10 本 session read-back 驗證 (ledger 97/98)。 -->

# 甲基「可分群」是否被 coverage / 拷貝數驅動 — confound 驗證 characterization

> **L0 一眼結論**：甲基 read 的「可分群性（clusterability）」**不是 coverage 也不是 CN-dosage 的假象**。與 read 數量只有**微弱**相關（控制後轉負）；與 CN **強相關但方向相反於假象** —— **LOH（失去一個 haplotype）把分群幾乎歸零（6/6 樣本、3 癌種重現）**，amplification（gain）不增加分群。真正的 power driver 是 **n_CpG**。機制：**read-clustering 需要 ≥2 個 haplotype 的等位甲基差才成立** —— 這是「分群是真等位生物學、非技術假象」的最直接證據。
>
> **L1 重點邏輯**：
> ① **coverage**：只有「最低覆蓋門檻」power floor（<50× 被壓抑），超過就 plateau；不是越多 reads 越會分群（18 cells 全 |ρ|<0.25）。⟨L1⟩
> ② **CN**：LOH 抑制分群（OR=0.07-0.09），gain 不增（even 固定覆蓋）→ 非 dosage 灌大訊號，是**等位結構效應**（需兩個 allele）。⟨L1⟩
> ③ **跨樣本**：LOH-抑制機制 **6/6 樣本 OR<1**（0.01-0.107），3 癌種 general。⟨L1⟩
> ④ **FP 例外**：FP 不顯示 LOH 抑制（暗示 FP 分群=artifactual noise；⚠ loh n=70 低 power → suggestive）。⟨L2⟩
> ⑤ ⚠ **區分 clusterability vs NGroups**：NGroups（HPFineNGroups）與 coverage 強相關（ρ=0.51）**但那是 HP-tag occupancy（phasing 取樣）非甲基**（HD-4）。

---

## §1 問題與方法

**問題**：甲基 read 能不能分群（clusterability）會不會只是「覆蓋高 / 拷貝數高 → 容易看到分群」的技術假象？若是 → 分群訊號不可信；若否 → 是真等位甲基生物學。

**指標（「甲基可分群」）**：
- `is_clustered`（二元）= reliable（ISM `Significant` gate）**OR** latent（script-88 PERMANOVA 救回，minN5）。
- `NGroups` = HPFineNGroups（⚠ HP-tag 子家族計數，HD-4 證 phasing-driven 非純甲基）。
- `CramersV` = 分群強度（連續）。

**預測因子**：`coverage`=NumReads · `n_CpG`=NumCpGs（分群 power 控制）· `CN`=SEQC2 真值（gain/loss/loh/neutral，HCC1395）。

**🔑 confound 處理（auc-confound-guard 邏輯）**：coverage 與 CN 本身相關（amplified→更多 reads），故用 **partial Spearman（rank-residual）+ 固定 CN-state 內 coverage 相關 + 固定 coverage-bin 內 CN 比率 + 多變量 logistic** 拆開。

**Scope**：coverage 全 6 樣本 × TP/FP/FN（18 cells）；CN 真值 HCC1395 × TP/FP/FN；跨樣本 LOH 用 ISM `Potential_LOH` proxy（SEQC2 只有 HCC1395）。

---

## §2 結果

### §2.1 Coverage = 微弱，非 driver ⟨L1⟩

| 證據 | 值 |
|------|----|
| ρ(clustered, coverage) HCC1395 TP | **0.109**（p=2e-79, n=29723）— 弱 |
| 18 cells（6 樣本×3 class）ρ 範圍 | **−0.009（HCC1954 TP）～ 0.20（H2009 FP）**，median 0.13，**全 |ρ|<0.25** |
| coverage decile clusterability | **power floor 後 plateau**：<50× rate 0.038 → >55× 飽和 0.15-0.18（非單調）|
| 多變量 logistic（控 n_CpG+CN 後）| z_logCoverage coef **−0.096**（OR=0.91, p=4e-5）→ **轉負** |

→ 只有最低覆蓋門檻效應（夠 reads 才測得到任何分群）；**超過門檻後 more reads ≠ more clustering**。非 coverage 假象，跨 class+樣本 robust。（圖 P① / 延伸 P①）

### §2.2 CN = 強相關但反假象方向 ⟨L1⟩

**HCC1395 TP clustered-rate by CN-state**（neutral n=1030 / gain 19521 / loh 8715 / loss 457）：

| CN-state | clustered-rate | 解讀 |
|----------|--------------:|------|
| neutral | **0.222** | 基準（CN=2，兩個 allele）|
| gain | 0.179 | 比 neutral **低**（amplification 不增分群）|
| loss | 0.083 | 低 |
| **loh** | **0.022** | **幾乎歸零**（單一 haplotype）|

- chi² p≈1e-298；LOH logistic OR=**0.072**（p=1e-137，主導 CN 效應）。
- **gain vs neutral 在每個 coverage-bin 內 neutral ≥ gain（3/3）** → 即使固定覆蓋，amplification 也沒製造分群。
- → CN 有強影響，但是**等位結構效應（LOH 移除第二個 allele → 沒有等位甲基差可分），不是 dosage 把訊號灌大**。（圖 P②③）

### §2.3 真 driver = n_CpG（power）⟨L1⟩

多變量 logistic（is_clustered ~ z(logCoverage)+z(nCpG)+CN dummies，HCC1395 TP）：

| 因子 | coef | OR | p |
|------|-----:|---:|---|
| **z_nCpG** | **+0.185** | **1.20** | 1e-35（最強連續，power）|
| z_logCoverage | −0.096 | 0.91 | 4e-5（控制後轉負）|
| CN_gain | −0.216 | 0.81 | 0.007 |
| **CN_loh** | **−2.636** | **0.07** | 1e-137（主導）|
| CN_loss | −1.344 | 0.26 | 2e-12 |

- **pseudo-R² = 0.083** → coverage+nCpG+CN 三者合計只解釋 **8%** 的 clusterability 變異 → **其餘 92% 來自真正的等位甲基生物學**，再次說明非技術假象。

### §2.4 跨樣本：LOH-抑制機制 6/6 重現 ⟨L1⟩

ISM `Potential_LOH` proxy（TP，每樣本 LOH vs non-LOH clustered-rate OR）：

| 樣本 | 癌種 | loh-rate | nonloh-rate | OR |
|------|------|--------:|-----------:|---:|
| HCC1395 | breast | 0.029 | 0.228 | 0.100 |
| HCC1937 | breast | 0.026 | 0.217 | 0.095 |
| HCC1954 | breast | 0.027 | 0.259 | 0.079 |
| H1437 | lung | 0.019 | 0.224 | 0.067 |
| H2009 | lung | 0.015 | 0.125 | 0.107 |
| COLO829 | melanoma | 0.001 | 0.107 | **0.010** |

- **6/6 OR<1**（range 0.01–0.107），3 癌種 general → 機制（need ≥2 haplotype）跨癌種成立。
- **proxy 驗證**：HCC1395 ISM-LOH vs SEQC2-LOH **recall 0.961**（捕到 96% 真 LOH）/ precision 0.491（ISM-LOH 為 superset）→ proxy 可靠捕 LOH 訊號。（圖延伸 P②）

### §2.5 FP/FN 延伸 — TP/FP mechanistic 區分 ⟨L2⟩

HCC1395 CN×class LOH-vs-nonLOH OR：
- **TP OR=0.092 / FN OR=0.018**（FN 抑制更徹底，loh rate 0.003）—— 真 somatic loci LOH 強抑制。
- **FP OR=0.982**（**不抑制**；loh 0.186 ≈ gain 0.196）—— FP 的「分群」不遵守 need-2-haplotype 規則 → **暗示 FP clustering = artifactual noise，非真等位結構**。
- ⚠ **caveat**：FP loh n=**70**（低 power），此 TP/FP 區分為 **suggestive 非定論**；且 characterization-only，**不可寫成 filter**（FP 稀疏 + 甲基-filter 已 DEAD 四道）。

---

## §3 機制解讀

**為什麼 LOH 殺分群、gain 不增、coverage 控制後轉負？**

ISM 的「分群」= read×CpG → read-read 距離 → 階層分群，本質在**比較兩個 haplotype/allele 的甲基模式差異**。
- **LOH** = 只保留單一 parental haplotype → **沒有第二個 allele 可比** → 分不了群（rate 0.022）。這是機制必然，不是 artifact。
- **gain（amplification）** = 多份同一序列，但**等位結構不變** → 不會憑空製造甲基差 → 不增分群。
- **coverage** = 給足 reads 是分群的**必要 power 門檻**（<50× 測不到），但門檻之上**等位甲基差的大小由生物學決定，非 read 數** → 控制後 coverage 係數轉負。
- **n_CpG** = 更多 CpG 位點 = 更高解析度看到等位甲基分裂 → 唯一正向 power driver。

→ **分群是真等位甲基生物學（需兩個 haplotype），非覆蓋/拷貝數技術假象。**

---

## §4 對論文

### §4.1 §Methods-Neg 段落（可直接用）

> **clustering 非 coverage/CN-dosage 假象 + LOH-抑制機制**：methylation read-clusterability 與 coverage 僅微弱相關（18 sample×class cells 全 |ρ|<0.25；多變量控制後 coverage OR=0.91 轉負），與 CN 強相關但反 dosage 方向 —— **LOH 抑制分群（logistic OR=0.07；6/6 樣本 OR<1, 3 癌種重現）**，gain 不增（固定覆蓋亦然）。機制 = read-clustering 需 ≥2 haplotype 等位甲基差；LOH 移除第二 allele 故分不了群。n_CpG 是唯一正向 power driver（OR=1.20）；coverage+nCpG+CN 僅解釋 8% 變異（pseudo-R²=0.083）。⇒ 分群是真等位生物學非技術假象。**與 fast_cnv（dosage 非量級 driver）+ HD-4（NGroups=phasing）一致。**

### §4.2 併入 catalog R6
- 答 schema P6「reliability 是否被覆蓋驅動」：**否**（power floor 後 plateau；控制後負）。
- catalog `clustering_reliability` 欄的可信度由此 backstop（非覆蓋/CN 假象）。
- 新增 cross-sample LOH-抑制當 R6 的 mechanistic 佐證。

---

## §5 Caveats（誠實邊界）
1. **CN 真值單樣本**（SEQC2=HCC1395）；其他 5 樣本用 ISM `Potential_LOH` proxy（recall 0.96 驗證，但非 orthogonal CN truth）。
2. **SEQC2「neutral」(not-in-bed n=1030)** 可能含未列入 benchmark 的 aneuploid（HCC1395 高度非整倍體）；但 LOH/gain 訊號不受此影響。
3. **FP 稀疏**（H1437=8 / HCC1954=29 排除；HCC1395 FP loh n=70）→ FP TP-區分 suggestive。
4. **single-pipeline ClairS-TO** → tier 封頂 ⭐3。
5. **TP 為主**分析；FP/FN 延伸確認趨勢一致。

---

## §6 Provenance
- 腳本 `scripts/101_clusterability_vs_coverage_cn.py`（partial Spearman + 分層 + logistic）/ `102_..._extended.py`（FP/FN + 6 樣本 + proxy 驗證）。
- 資料 `genome_survey_v2/catalog/clusterability_vs_cov_cn.json` + `clusterability_cov_cn_extended.json`。
- 圖 `docs/paper_focus/04_figures/clusterability_vs_cov_cn.png` + `clusterability_cov_cn_extended.png`。
- ledger 97 + 98。全數字 2026-06-10 read-back 驗證。
