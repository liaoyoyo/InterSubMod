---
title: TP vs FP × 結構-標籤關聯狀況刻畫（HCC1395 paired_full，單樣本）
date: 2026-06-19
status: in_progress
task_type: B-validation (partial — 單樣本 subset)
tier: ⭐3 (single-sample single-pipeline 封頂)
partial_flag: true  # 僅 HCC1395；FN 維度未含（Phase 2）
data_sources: >
  docs/experiments/in_progress/2026/06/20260619_tp_fp_structure_association_HCC1395/_assets/crosstabs.json,
  docs/experiments/in_progress/2026/06/20260619_tp_fp_structure_association_HCC1395/_assets/confound_control.json,
  docs/experiments/in_progress/2026/06/20260619_tp_fp_structure_association_HCC1395/_assets/run_meta.json
scripts:
  - scripts/analysis/tp_fp_structure_label_association.py
  - scripts/analysis/tp_fp_structure_label_figures.py
inputs:
  - big7_disk_output/canonical/HCC1395/paired_full/20260420_HCC1395_paired_full_full_2/intersubmod_{tp,fp}/significance_summary.csv
related:
  - docs/methodology/20260617_structure_label_situation_inventory_01.md
  - docs/methodology/20260616_structure_label_verdict_false_negative_audit_01.md
  - memory project_ism_complete_tpfpfn_existence_cis
---

# TP vs FP × 結構-標籤關聯狀況刻畫（HCC1395）

> ⚠️ **partial**：僅 HCC1395 單樣本；**FN 維度未納入**（FN 無 region-level 結構資料，需 Phase 2 重跑 ISM）。
> 🔴 **定位**：本分析是 **描述性刻畫（characterization）**，**不是** 復活「結構/ASM 當 TP/FP filter」（該方向已 NEGATIVE/DEAD）。下述差異是**富集（enrichment）非可分離（separation）**，不可作為篩選器。

## BLUF（一句話）

TP 與 FP 在**甲基結構的「軸」上有清楚且覆蓋穩健的質性差異** — **TP 的結構偏 allele 軸（ALT/REF），FP 的結構偏 germline-HP 軸** — 但這是**富集而非排他**：沒有任何結構型態是真實變異獨有的，且兩類 90%+ 區域結構都很弱。〔L1 數字 + L3 詮釋〕

---

## 一、直接回答你的三個問題

**Q1「是否某種高關聯只在 TP（或 TP+FN）出現？」→ 否，無排他性。**〔L1〕
- 每一種結構型態在 TP 與 FP 都出現；最偏 TP 的 `allele_only` 軸在 FP 仍佔 **27.8%**。
- 僅有的「FP=0」格是極小的 CramérV 中段箱（`[.1,.2)`、`[.2,.3)`）與 cluster-PERMANOVA `nonsig`（TP 僅 2 個）。對這些格做 **rule-of-three**：FP n=627 下，「0 觀測」的真實 FP 率 95% 上界仍可達 **0.48%** → 是**抽樣不足非真不存在**，不構成「TP 獨有」。
- 結論：**不存在「只在真實變異出現」的結構簽名**（與 filter-DEAD 一致）。

**Q2「哪種狀況是 FP 常見的？」→ FP 的結構 = germline-haplotype（HP 軸）甲基差異，非 somatic-allele 驅動。**〔L1 富集 + L3 詮釋〕FP 相對富集於：
| FP 常見狀況 | TP% | FP% | 富集(FP/TP) | Fisher q(BH) | 覆蓋配對後存活？ |
|---|--:|--:|--:|--:|:--:|
| `hp_only`（只 HP 軸顯著）| 1.96 | **8.45** | **4.31×** | 2.9e-17 | ✅ OR 4.71 (p=3e-14) / CMH OR 4.10 |
| HP 軸 \|Δβ\|≥0.2 | 0.78 | **5.74** | **7.36×** | 2.4e-18 | （未單獨配對；CMH 見下文） |
| SampleASM \|Δβ\|≥0.2（tumor−normal）| 28.3 | **54.9** | 1.94× | 1.4e-41 | — |
| `neither`（無任何標籤軸結構）| 8.23 | **17.1** | 2.07× | 5.4e-12 | ✅ OR 2.17 |
| Noise verdict | 26.3 | **33.0** | 1.25× | 3.6e-4 | ✅ CMH OR 1.31 (但 Breslow-Day 異質) |
| within-HP 多群(proxy) | 60.1 | **76.7** | 1.28× | 1.6e-17 | — |

**Q3 反面 — TP 常見狀況 = allele 軸結構 + cluster×label 對齊 + Strong/Subclone verdict。**〔L1〕
| TP 偏多狀況 | TP% | FP% | OR(FP在該類) | 覆蓋配對後存活？ |
|---|--:|--:|--:|:--:|
| `allele_only`（只 allele 軸顯著）| **51.6** | 27.8 | 0.36 | ✅ OR 0.35 (p=1.6e-30) / CMH OR 0.44 |
| aligned（cluster×label 顯著對齊）| **8.34** | 2.39 | 0.27 | ✅ OR 0.34 / CMH OR 0.33 |
| Strong verdict | **6.52** | 2.07 | 0.30 | ✅ CMH OR 0.40 (p=7e-4) |
| allele 軸 \|Δβ\|≥0.2 | **7.31** | 3.19 | 0.42 | — |
| cluster-PERMANOVA 顯著 | **9.50** | 3.83 | 0.38 | — |

---

## 二、核心發現：差異在「結構的軸」，且**覆蓋穩健**

最強且最可解釋的維度是 **label-PERMANOVA 的 HP 軸 vs allele 軸 2×2**（omnibus χ²=260, p=4e-56）：

| label 軸 | TP%（全集）| FP%（全集）| TP%（配對）| FP%（配對）|
|---|--:|--:|--:|--:|
| `allele_only` | 51.6 | 27.8 | 52.5 | 27.8 |
| `both` | 38.2 | 46.7 | 36.9 | 46.7 |
| `hp_only` | 1.96 | 8.45 | 1.92 | 8.45 |
| `neither` | 8.23 | 17.1 | 8.65 | 17.1 |

→ 見 `fig1_label_axis_full_vs_matched.png`：**全集與覆蓋配對後幾乎重疊**，證明此差異**不是覆蓋假象**。

**為何覆蓋是頭號 confound、為何已被控制**〔L1〕：
- FP 的 NumCpGs 中位 **58** < TP **76**（NumReads 反而 FP 133 ≳ TP 124；LOH 平衡：TP 47.4% / FP 48.2%）→ FP 結構偵測力天生較低。
- 控制兩法並行且結論一致：(a) **配對子集** TP=3121/FP=627，NumReads/NumCpGs KS **p≈1**（完美平衡）；(b) **CMH** 跨 NumCpGs×NumReads 9 層校正。`hp_only` 的 FP 富集在 CMH 下 OR=**4.10**、Breslow-Day p=0.48（跨層同質）；`allele_only` 的 TP 富集 CMH OR=**0.44**。

**L3 生物學詮釋（推論，非 L1）**：真實 somatic SNV（TP）若有甲基結構，傾向**跟隨 allele 軸**（ALT/REF）— 即真實變異 cis-影響甲基或標記真實 (sub)clonal lineage；偽陽性（FP）位點非真 somatic，其 ALT「等位基因」是比對/定序假象，故無 allele-following 結構，但該位點常落在 **germline 等位特異甲基（germline ASM）** 區 → 表現為 HP 軸結構或 tumor-normal 差異（SampleASM）。**此詮釋待 cis-control 與跨樣本驗證**（單樣本不可形式排除 cis-ASM）。

---

## 三、與既有結論的關係

- **精煉 filter-DEAD（不推翻）**：既有結論「ASM/結構非 usable filter」成立 — 因為(1) 差異是富集非排他，(2) 兩類 90%+ 為 `no_cluster`/`Weak`/低 CramérV。本分析**新增**的是：TP/FP 的弱結構在**軸向**上系統不同（allele vs HP），且覆蓋穩健。此為刻畫，非篩選力。
- **新拆解**：`20260617` situation inventory 與 `20260616` FN audit 都 pool 了 TP+FP；本報告首次按 TP/FP 拆 verdict/軸/Δβ。
- **與 within-HP 多群的關係**：within-HP 多群 proxy 在 **FP 更常見**（76.7% vs 60.1%）→ 佐證「within-HP 碎裂不是乾淨 subclone 簽名」（記憶 `project_subcluster_cluster_count_determination`）。

---

## 四、限制（Limitations）

1. 🔴 **PERMDISP 退化**：本 run `ClusterDispersionWarn` 全 false、`ClusterDispersionP<0.05` 計數 = 0 → **無法分離 location vs dispersion**。所有「PERMANOVA 顯著」可能含 dispersion 成分。`hp_only` 的 FP 富集在覆蓋配對（平衡組大小）後仍存活，削弱（但未排除）純 dispersion 解釋。
2. 🔴 **VerificationClass 是 gate 截斷標籤**：90.5% region 從未跑 ClusterPERMANOVA（valid TP 9.5% / FP 3.8%）；verdict 軸與原始 PERMANOVA 軸答不同問題，故**分表**呈現、不合併。
3. **`HP_Residual_Delta` ≈ `HPMergedDelta`**：兩者分箱近乎一致 → 本 run 的「normal 扣除」殘差與原始 HP-Δβ 高度相關；其作為 somatic-controlled 軸的有效性**待驗**（不可僅憑此宣稱已扣 germline）。與記憶中「HP_Residual 顯著率 TP≈FP≈33%」的張力＝**顯著率 vs 量級**兩種視角差異，待對齊。
4. **單樣本 single-pipeline**：HCC1395 特有結構（chr8 hotspot、LOH 區）可能驅動，跨癌種未必複現 → tier 封頂 ⭐3。
5. **FP 小樣本**：n=627，Subclone(FP=2)/Strong(FP=13) 明確 underpowered；FP 率 <0.3% 的型態無偵測力。
6. **read-set 混合**：aligned/cluster 軸是 paired（含 normal reads，約半訊號來自 normal）；已將 somatic-controlled 軸（HP_Residual、SampleASM）分表，但見限制 3。

---

## 五、下一步

- **Phase 2（需用戶再確認）**：在 truth-missed 位點跑 ISM 產 `intersubmod_fn/` → 檢驗「TP=allele軸」簽名是否延伸到 FN（真實但被漏的變異）。**FN confound 必標**：FN 被漏常因低覆蓋 → 結構弱可能是品質非生物。屬 Hard-Gate 計算步驟。
- 跨樣本（5/6 有 normal）驗證 allele/HP 軸不對稱是否複現。
- cis-control：用 normal-anchored 檢定形式區分 cis-ASM vs germline-ASM（HCC1395 有 normal BAM）。

---

## Provenance

- 數字全部來自 `_assets/{crosstabs,confound_control,run_meta}.json`（腳本 `tp_fp_structure_label_association.py` 產出，本輪 Read 讀回驗證；§13 序列：跑→落檔→讀回→才撰寫，分析與本文 Write 分批）。
- 輸入：`20260420_HCC1395_paired_full_full_2`（注意：另有 `20260314_..._complete_matrix`；本分析用 20260420）。
- 規模核對：TP+FP 邊際和 = 29754 + 627 ✓；situation inventory 用 30490（差 109，疑不同 binary/date，本分析以實際檔案為準）。
- 圖：`fig1_label_axis_full_vs_matched.png`、`fig2_verdict_alignment.png`、`fig3_delta_magnitude.png`。
