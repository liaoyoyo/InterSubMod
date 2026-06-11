<!--
建立時間: 2026-05-29
作者: Claude (Opus 4.8) + InterSubMod Research
報告類型: 既有 artifact 驗證 (existing-artifact verification) + Phase 2 全基因組擴展
任務類型: B Comprehensive validation (subset → genome-wide)
樣本 scope: 單樣本 HCC1395 paired_full (longphase_s tagged BAM)；⚠ partial: 無跨樣本
驗證原則: feedback_existing_artifacts_must_verify — 每個數字 grep source + cross-check + 本 session 親自重算
框架: BLUF + Pyramid
ledger: 20260529_ism_zar1l_brca2_asm_verification (tier ⭐3)
-->

# ISM ZAR1L/BRCA2 ASM 驗證報告 — 位點確認 + read-level magnitude 調和 + 全基因組擴展

> ✅ **MAGNITUDE RESOLVED（2026-05-30，見 corrects ledger `20260530_ism_meth_hp_somatic_full_tp_characterization`）**：§3 的「reconciled −0.05~−0.07」**已被修正取代**。根因確認為 **MSA Level1 把 5mC+5hmC 雙修飾各寫一列、原 pivot 每列獨立計數把 beta 砍半**（非 CpG-set 差異）。改用 max-collapse any-modification 口徑全 39,447 TP 重跑後，**BRCA2 正式值 = HP-axis Δβ −0.122 / ALLELE-axis −0.099**（與 script03 5mC-only −0.122 完全吻合）。全基因組 A-characterization 正式結論見 `InterSubMod/research/tsg_promoter_asm_reviewer/output/phase2_dataverify_synthesis.md`。方向 hypo robust 不變；本報告 §3/§5 的舊 magnitude 數字僅留作方法學除錯軌跡。
>
> **partial-scope flag**: 全部基於單樣本 **HCC1395 paired_full**；跨樣本一致性 NOT verified。

---

## 0. TL;DR（BLUF）

針對使用者「確認 ISM 當時是否真的對此 TP 位點輸出圖+數據、統計現象是否真實、熱力圖是否可解釋、現流程能否重現、是否能用閥值定位甲基位點、哪些突變強關聯、哪些位點聚集」7 問，逐一驗證：

**核心裁決：BRCA2/ZAR1L 位點的 ASM 現象「真實、可重現、抗多重檢定」，但效應量小（Δβ≈−0.05~−0.07），且全基因組層級 ASM 不集中於 TSG promoter。BRCA2 是統計穩健的中等效應展示位點，非最強 ASM。**

| # | 問題 | 裁決 |
|---|------|------|
| Q1 | ISM 有用 TP 輸出此位點圖+數據？ | ✅ 是 |
| Q2 | 統計確實有此現象？ | ✅ 真實穩健（抗 Bonferroni）|
| Q3 | 各參數+熱力圖可解釋？ | ✅ 可解釋一致 |
| Q4 | 現流程能找到此位點+強關聯？ | ✅ 能（cross-check PASS，全基因組 top 0.14%）|
| Q5 | ISM 能用閥值法+定位甲基位點？ | ✅ 能（本就是三層閥值法）|
| Q6 | 哪些突變有強關聯？ | ⚠ 嚴格雙過濾僅 19 個；BRCA2 不在強效應名單 |
| Q7 | 哪些位點明顯聚集？ | ✅ 45 個 spatial cluster |

---

## 1. 位點確認（Q1）— TP variant + ISM 輸出

| 項目 | 值 | 來源 |
|------|-----|------|
| 位置 | chr13:32,315,128 G>A | `output/tp_in_tsg_promoter_nonLOH.tsv:2` |
| TVAF / NVAF | 0.189 / 0.002 → 真 somatic | 同上 |
| 基因歸屬 | **BRCA2/ZAR1L(ZAR2) bidirectional overlapping promoter**（+51bp BRCA2 TSS / −235bp ZAR1L TSS）| `output/literature_scout_ZAR1L_BRCA2.md`（L1 Misra 2010, PMID 20202217）|
| HP coverage | total=123, HP1=46, HP2=24, HP1-1=45, pass_strict=true | `output/step3_hp_coverage.json` |
| ISM 數據輸出 | `output/step4_ism_results.json`（per-CpG per-HP β 矩陣）| ✓ |
| ISM 圖片輸出 | IGV session `20260527_ASM_ZAR1L_ex.png` + `brca2_read_level_methylation.png` + `brca2_per_cpg_delta.png` | ✓ |

**結論**：ISM 確實對此 TP 位點輸出了完整數據與圖片。圖標「ZAR1L Gene」是正確標註（bidirectional promoter）。

---

## 2. 統計現象驗證（Q2）+ 熱力圖可解釋性（Q3）

ISM 方法（`scripts/03_step4_ism_methylation_diff.py`）：比較 **HP1（germline haplotype reads）vs HP1-1（帶 somatic allele 的 reconstructed reads）** 在同一 allele 上的 per-CpG 甲基化（Wilcoxon paired signed-rank + Cohen's d）。

- HP1 vs HP1-1：somatic 端 **低甲基化**，Wilcoxon **p=6.09e-11**，方向明確
- Negative control（5 random regions）Δβ 全 <0.06、p>0.1 → 對照無此訊號（`step4_negative_control.json`）
- 熱力圖可解釋：read-level matrix 肉眼可見 HP1-1 低甲基化；per-CpG Δβ 負值**聚集在變異上游 −500~0bp 的 bidirectional promoter 區**，與 CpG island（chr13:32315396-763, 距變異 268bp）位置吻合

---

## 3. ★ Caveat 2 — read-level magnitude 調和（核心方法學發現）

驗證初期發現兩條 pipeline 對同一位點 Δβ 不一致（per-locus script03 −0.122 vs genome survey script15 −0.063）。經 read-level 調和釘死根因：

### 3.1 閾值策略不是主因（red herring）

同一份 MSA Level1 數據，0.5-binary vs 0.8/0.2-3state 幾乎相同：

| 閾值 | germ β | som β | Δβ | p |
|------|--:|--:|--:|--:|
| 0.5-binary | 0.175 | 0.120 | −0.055 | 4.8e-13 |
| 0.8/0.2-3state | 0.170 | 0.116 | −0.054 | 8.3e-15 |

### 3.2 真因 = CpG 集差異（REF-anchored vs read-mod-call-anchored）

| 計算 | n | germ β | som β | Δβ |
|------|--:|--:|--:|--:|
| script03 @ 全 412 CpG | 197 | 0.215 | 0.093 | **−0.122** |
| script03 @ 共同 170 CpG | 119 | 0.164 | 0.078 | **−0.086** |
| MSA @ 全 524 CpG | 265 | 0.170 | 0.116 | **−0.054** |
| MSA @ 共同 170 CpG | 131 | 0.108 | 0.040 | **−0.068** |

- 兩 pipeline CpG 集只有 **170 共同**（script03 有 412、MSA 524）；script03 有 **242 個位置 MSA 不認**（非 REF CpG，read-mod-call 抓到的雜訊位置）
- script03 限縮到共同 CpG 後 Δβ 從 −0.122 **掉到 −0.086** → 那 242 個非標準 CpG **膨脹了 ~0.036**
- 殘差（共同 170 上 MSA −0.068 vs script03 −0.086，~0.018）來自 MSA 合併雙股 + REF 定義較嚴

### 3.3 裁決

> **方向（hypo）+ 顯著性（p≤1e-6）在所有切法下 100% robust；magnitude 方法依賴。可辯護單一數字 = MSA REF-anchored Δβ ≈ −0.05~−0.07。script03 的 −0.12 是高估（混入非標準 CpG）。**

---

## 4. 現流程可重現（Q4）+ 閥值定位甲基位點（Q5）

- **Q4 可重現**：5/28 的 4-agent MSA method audit（`output/msa_audit_synthesis.md`）結論 — C++ 工具 MethylHaploExtractor 認得 HP:Z:1-1 tag，BRCA2 單點重現成功，MSA 與 Python pivot 方向一致（cross-check PASS）。全基因組排名 **#25/18,149（top 0.14%）by p-value**。
- **Q5 閥值法定位**：ISM 本就是三層閥值法 — (a) per-CpG ML 閾值（0.8 meth / 0.2 unmeth）→ (b) per-CpG \|Δβ\| / max_abs_delta **定位驅動 ASM 的單個 CpG**（per_cpg_delta 圖中 −500~0bp 的負值 CpG）→ (c) locus-level \|Δβ\|≥0.1 & p<0.05 定義 sig ASM。✅ 能用閥值定位甲基位點。

---

## 5. Phase 2 全基因組擴展（Q6 / Q7）

genome survey **已完整跑完全 22 autosomes**（14,222 變異 / 18,149 variant-axis records；`genome_survey/genome_asm_summary.tsv`）。

### 5.1 多重檢定校正（必做）

| 篩選 | 通過數 |
|------|------:|
| 未校正 p<0.05 | 3,630（隨機期望 ~907 FP）|
| Bonferroni（p<2.75e-6）| **90** |
| BH-FDR 5% | 832 |
| **Bonferroni + \|Δβ\|≥0.1（強關聯雙過濾）** | **19**（8 hypo / 11 hyper）|

**BRCA2 抗 Bonferroni ✓**（p=1.84e-10）但 \|Δβ\|=0.063 < 0.1 → 不進 19 強名單。

### 5.2 Q6 — 強關聯突變（Bonferroni + \|Δβ\|≥0.1 的 top hypo）

| chr:pos | axis | nCpG | Δβ | p |
|---------|------|--:|--:|--:|
| chr1:200029019 | HP1_vs_HP1-1 | 50 | −0.350 | 2.8e-08 |
| chr20:50267392 | HP1_vs_HP1-1 | 74 | −0.238 | 3.1e-10 |
| chr2:94874876 | HP1_vs_HP1-1 | 86 | −0.214 | 2.9e-10 |
| chr15:67823497 | HP2_vs_HP2-1 | 261 | −0.155 | 3.3e-25 |
| chr22:34684238 | HP1_vs_HP1-1 | 154 | −0.131 | 1.2e-15 |

### 5.3 Q7 — 空間聚集（45 clusters, gap<10kb, ≥2 sig sites）

| cluster | sites | mean Δβ | best p |
|---------|--:|--:|--:|
| chr16:34637913-34640158 (2245bp) | 4 | −0.177 HYPO | 3.3e-04 |
| chr4:135127940-135135918 (7978bp) | 3 | −0.190 HYPO | 1.6e-03 |
| chr5:20193055-20193210 (155bp) | 3 | −0.165 HYPO | 8.4e-04 |

---

## 6. ⚠ Phase 2 三個結構性誠實發現

1. **BRCA2 統計穩健但效應量不在強者之列**：抗 Bonferroni ✓，但 \|Δβ\|=0.063 < 0.1；最強 hypo 在 chr1:200029019(−0.350) 等。BRCA2 = 中等效應展示位點，非最強 ASM。
2. **TSG promoter enrichment = 0.00×**：582 sig ASM 中 0 個落在 TSG promoter（主因：僅 3/14,222 變異落在狹窄 TSG promoter 窗 → underpowered，非「證明無關」）。全基因組 ASM 散佈，不集中 TSG promoter → BRCA2 是手挑展示例。
3. **方向以 hyper 為主**：強關聯 hyper(11) > hypo(8)；BRCA2 是 hypo（少數方向）。若 reviewer claim 為「TSG 低甲基化」，全基因組圖景 hyper 更常見。

---

## 7. Verdict + Tier

- **Verdict**: POSITIVE-but-modest — BRCA2/ZAR1L ASM 真實、可重現、抗 Bonferroni；但 effect size 小（−0.05~−0.07）、非 TSG-enriched、單樣本。
- **Tier**: ⭐3（單樣本 + reconciled magnitude marginal + cross-pipeline 方向一致 + 抗多重檢定；未達 ⭐4 因無跨樣本 + effect 小）。
- **4 caveat 必載**：(1) magnitude 方法依賴 −0.05~−0.12；(2) TSG enrichment=0；(3) hyper>hypo；(4) single-sample HCC1395。

---

## 8. Provenance（本 session 親自重算驗證）

| 數字 | source |
|------|--------|
| TP variant / HP coverage | `research/tsg_promoter_asm_reviewer/output/{tp_in_tsg_promoter_nonLOH.tsv, step3_hp_coverage.json}` |
| ISM stats (BRCA2 Δβ=−0.122) | `output/step4_ism_results.json` |
| MSA Level1 重算（Δβ=−0.054）| `msa_brca2/output/brca2_cloned/level1_raw_methylation_details.tsv.gz` |
| genome survey + 多重檢定 | `genome_survey/genome_asm_summary.tsv`（18,149 records）|
| 基因歸屬 | `output/literature_scout_ZAR1L_BRCA2.md`（L1 PMID 20202217）|
| method audit | `output/msa_audit_synthesis.md` |

**未驗證 / TBD**：跨樣本一致性（僅 HCC1395）；TSG enrichment underpowered（n=3）；strand-merge 殘差未逐 read 拆解。
