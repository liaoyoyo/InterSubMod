<!--
建立時間: 2026-07-15T04:35:00+08:00
目標: 全量驗證 truth-labeled FP 單一 sSNV 的 focal-ALT reads 出現無監督甲基多群時，能否視為 subclone 或 linear evolution
處理範圍: 7 datasets / 6 biological samples；chr1-22 biallelic LongPhase-S recalibrated FILTER=PASS sSNV；FP 全量 7,745 sites + 1:1 matched TP controls
關聯檔案: InterSubMod/research/20260715_single_fp_alt_multicluster_subclone_validation/results/report_dataset_v1/report_dataset.json
狀態: validated · external Claude Code review complete · 2026-07-15
敘述框架: SCQA + hypothesis-falsification + claim-evidence tier
服務目標: G3 / G4 / G5
-->

# 單一 truth-FP sSNV 的 ALT reads 出現甲基多群，能否視為 subclone？

## TL;DR：不能直接判定為 subclone，也不能由單一 sSNV 推定 linear evolution

**直答結論：不可以把「單一 truth-labeled FP sSNV 的 focal-ALT reads 被 InterSubMod 無監督分成多群」直接解讀為高機率 subclone；更不能因為腫瘤可能呈 linear evolution，就從一個 sSNV 推定 linear ancestry。**

可被數據支持的較窄敘述是：

> 在最新 LongPhase-S recalibrated `FILTER=PASS` 的 `7,745` 個 truth-labeled FP sSNV 中，部分 focal-ALT readsets 呈現可重現的 read-level methylation heterogeneity；經 null、技術軸、REF 背景與 assignment 穩定性檢核後，剩 `10` sites / `9` unique ALT readsets / `6` layered regional components 可列為 **strict epigenetic follow-up candidates**，但沒有任何位點獲得 orthogonal cellular-subclone 或 focal linear-topology confirmation。

| 問題 | 結論 | 證據等級 | 影響 / 信心 |
|---|---|---|---|
| ALT reads 是否真的可分成甲基多群？ | 部分可；primary stable `583/4,967 = 11.74%` | L2：同 pipeline 直接量測 + null/seed sensitivity | 中 / 高 |
| 這些多群是否 FP-specific？ | 高門檻 endpoint 不具 specificity；FP `1.330%`、matched TP `1.265%`，`p=0.772` | L3：matched-control 反證 | 高 / 高 |
| 多群是否等於 subclone？ | 不支持；orthogonal confirmation `0` | L4：仍是假說 | 高 / 高 |
| 單一 sSNV 是否支持 linear evolution？ | 結構上不可識別；focal identifiable `0` | L5：缺第二 mutation state | 高 / 高 |
| 是否仍有後續價值？ | 有；strict `10/4,967 = 0.201%` 可優先做正交驗證 | L2：保守 sensitivity candidate | 中 / 中高 |

![證據鏈：由 7,745 個 FP 收斂到 10 個 strict follow-up candidates](figures/01_evidence_chain.png)

---

## SCQA：真正要驗證的不是「能否分群」，而是「分群是否代表同一細胞譜系」

**Situation**：ONT read 同時保留 focal allele、長距離 haplotype 與 MM/ML methylation；InterSubMod 可以在同一 focal ALT readset 內做無監督甲基分群。

**Complication**：read-level methylation heterogeneity 也可能由 coverage、missingness、read length/strand/start/end、HP family、mapping、normal/germline contamination、CN/LOH、區域性 epigenetic state 或 stochastic epimutation造成。Cluster 是統計物件，不自動等於細胞 clone。

**Question**：在最新 LongPhase-S PASS 的所有 truth FP 中，這種多群有多少？是否比 matched TP 多？是否經更嚴格 null、REF 背景、技術重複與 topology 後仍支持 subclone／linear evolution？

**Answer**：少數 readsets 有穩定甲基多群，但高門檻 endpoint 不具 FP specificity，strict candidates 僅 `10`；沒有第二 somatic haplotype tag、沒有 orthogonal clone marker、沒有 focal linear-topology identifiability。因此合理產物是 follow-up candidate list，不是 subclone callset 或 evolution tree。

---

## 輸入契約：主分析只使用同一次最新 LongPhase-S 的 recalibrated PASS sSNV

### 來源與 scope

| 項目 | 實際來源 / 數量 | 驗證 |
|---|---|---|
| 最新 canonical manifest | `/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260713_layered_reconstruction_v3_raw_all_lps_pass_v5/input_manifest.snapshot.json` | SHA-256 `16f2ef...c1a45` |
| LongPhase-S producer | `/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/research_rounds/20260711_longphase_s_raw_all_production_sidecars_v2/` | 同 run recalibrated VCF + exact HP/PS sidecar |
| 分析 scope | chr1-22、biallelic sSNV、LongPhase-S recalibrated `FILTER=PASS` | 7 datasets / 6 biological samples |
| truth split | 下游 post hoc label；不回饋 LongPhase-S 或 InterSubMod | FP `7,745`；TP `335,296` |
| FP overlap | latest ∩ legacy canonical `3,302` | latest-only / previously unmaterialized `4,443` |

### focal ALT 與 HP tag 是兩個不同欄位

- **focal ALT read**：只由每個 region 的 `reads/reads.tsv: alt_support=ALT` 定義；不能用 HP tag 代替。
- **LongPhase-S HP**：由最新 sidecar materialize 到 BAM，作 carrier / haplotype annotation。
- **methylation cluster label**：Python null-gated clustering 的統計群；報告統一改稱 `Methyl group A/B/C`。
- 即使原始 internal label 同樣長得像 `1-1`、`1-2`，它也不是 LongPhase-S HP `1-1`、`1-2`。

### 最新 tagged BAM 已輸出，但沒有覆寫原始或 canonical BAM

7 個 bounded tagged BAM 全部位於：

`/big7_disk/liaoyoyo2001/big7_disk_output/synthesis/observation_workspaces/20260715_single_fp_alt_multicluster_subclone_validation/latest_tagged_subset/`

實際 materialization receipt：

- unique alignments：`885,047`
- MM+ML 完整：`885,047/885,047`
- sidecar missing：`0`
- duplicate identities 依預註冊規則收斂：`86,291`
- output size：`14,507,454,186` bytes
- canonical/output 原路徑覆寫：`0`

HCC1937 的 `28,787` 個 duplicate identity 是「一份 MM/ML 完整、一份缺失」，已選完整 payload；HCC1395 僅 `1` 個 duplicate 有非分析 aux-tag 差異，沒有 sequence/quality/MAPQ/CIGAR 或 complete methylation payload 衝突。

---

## 完整流程：每一步都以可觀察輸出驗證

| Step | 方法 | 實際驗證 |
|---|---|---|
| 1. 鎖定 sSNV | 同一次 LongPhase-S recalibrated PASS，之後才 truth split | FP `7,745`、TP `335,296`；manifest hash 固定 |
| 2. materialize tags | raw methylation BAM + latest exact-coordinate HP/PS sidecar | 7 BAM pass；MM/ML `100%`；sidecar missing `0` |
| 3. 跑 InterSubMod | 每一 FP region 輸出 reads、methylation、distance | `7,745/7,745` per-site outputs complete；C++ exit `0` |
| 4. 選 focal ALT | `alt_support=ALT`；outlier peel 後至少 6 reads | evaluable `4,967/7,745 = 64.13%` |
| 5. primary multigroup | phylo-v4.1 column-permutation null95、10 seeds、modal K | stable `583`；assignment ARI gate `582/583` |
| 6. 排 measured axes | HP、strand、start/end/length、MAPQ、called CpG | residual `453`；高門檻 `103` |
| 7. matched TP | 同 dataset 以 logit AF / log DP / log GQ 最近鄰 1:1 配對 | `7,745` pairs；caliper sensitivity 另列 |
| 8. 背景與 strict null | tumor REF、row-circular null、empirical p、ARI | strict `10` sites / `9` readsets |
| 9. topology | 同一 latest layered v5 的 regional context | focal linear identifiable `0`；subclone confirmed `0` |
| 10. replication | HCC1395 與 Dorado exact-site 技術重複 | common `204`；both evaluable `128` |

`significance_summary` 少 `119` 個 site：其中 `46` 可評估、`6` stable、`5` residual。主分析直接讀完整 per-site reads/methylation/distance，不依賴 summary row，因此 `7,745` 分母未遺失；coverage omission 已另存 audit。

---

## 全量結果：11.74% 有 null-gated 多群，但只有 0.201% 通過 strict follow-up gate

### 兩種分母必須分開

| Endpoint | 數量 | evaluable 分母 | 全部 FP 分母 | 正確解讀 |
|---|---:|---:|---:|---|
| Evaluable ALT readset | 4,967 | 100% | 64.13% | 有足夠 ALT reads 與完整 distance |
| Forced silhouette multigroup | 891 | 17.94% | 11.50% | historical sensitivity；無 null gate |
| Stable null multigroup | 583 | 11.74%（95% Wilson 10.87-12.66%） | 7.53% | primary epigenetic heterogeneity endpoint |
| Residual multigroup | 453 | 9.12%（95% Wilson 8.35-9.95%） | 5.85% | 未由已量測技術軸解釋 |
| HP-tag-covered high threshold | 103 | 2.07% | 1.33% | 高門檻候選；不是 phase-confirmed clone |
| Strict follow-up | 10 | 0.201% | 0.129% | 保守 sensitivity 交集；不是 FDR / PPV |

forced silhouette 與 null-gated stable **不是巢狀流程**：兩者交集 `364`、forced-only `527`、stable-only `219`、Jaccard `0.328`。因此不可把 `891 -> 583` 描述成同一 funnel 的淘汰率。

### site 數不是獨立觀察數

- `583` stable sites 有 `14,619` 次 read assignment。
- 只有 `10,577` 個 unique dataset-read；`2,561` 條 read 被多個 sites 重用。
- 共享 read 連圖後形成 `380` connected components，其中 `105` 個含多個 sites，最大 component 有 `38` sites。
- 因此 site-level Wilson CI 只作描述；不能把每個 site 當完全獨立生物 replicate。

---

## 7 datasets 的比例差異大，而且偵測率受 ALT read depth 強烈影響

| Dataset | FP sites | Evaluable | Stable / evaluable | High threshold / evaluable | Strict | Matched stable RD | Matched high-threshold RD |
|---|---:|---:|---:|---:|---:|---:|---:|
| HCC1395 | 1,603 | 921 | 119 / 12.92% | 19 / 2.06% | 4 | +4.68 pp | +0.87 pp |
| HCC1395 Dorado | 904 | 540 | 92 / 17.04% | 28 / 5.19% | 1 | +1.55 pp | +0.44 pp |
| COLO829 | 3,224 | 2,466 | 230 / 9.33% | 32 / 1.30% | 4 | +0.47 pp | -0.22 pp |
| H1437 | 460 | 150 | 12 / 8.00% | 1 / 0.67% | 1 | -0.43 pp | 0.00 pp |
| H2009 | 735 | 328 | 60 / 18.29% | 7 / 2.13% | 0 | -3.13 pp | -1.09 pp |
| HCC1937 | 620 | 435 | 43 / 9.89% | 10 / 2.30% | 0 | +0.16 pp | +0.48 pp |
| HCC1954 | 199 | 127 | 27 / 21.26% | 6 / 4.72% | 0 | +0.50 pp | -0.50 pp |

HCC1395 與 HCC1395 Dorado 是同一 biological sample 的不同技術 dataset，不能算兩個 biological replication。COLO829 目前 truth split 沒有 confident-region BED，故其 `3,224` FP 的 label confidence 低於其他資料集，必須保留 caveat。

![各 dataset 的 FP-TP paired risk difference 方向不一致](figures/03_per_sample_risk_difference.png)

ALT read depth 對偵測率有明顯 power effect：stable 從 6-9 reads 的 `6.02%` 升至 >=40 reads 的 `19.23%`；高門檻 endpoint 從 `0%` 升至 `4.23%`。不同 dataset 的 coverage mixture 可直接改變表面 prevalence。

![ALT read depth 與 endpoint 發現率](figures/07_coverage_dependence.png)

---

## Matched TP 反證：高門檻多群沒有 truth-FP specificity

### 全部 7,745 matched pairs

| Endpoint | Truth FP | Matched truth TP | Paired RD | Exact McNemar p |
|---|---:|---:|---:|---:|
| Evaluable | 4,967 / 64.13% | 5,814 / 75.07% | -10.94 pp | `1.83e-131` |
| Stable | 583 / 7.53% | 502 / 6.48% | +1.05 pp | `0.00864` |
| Residual | 453 / 5.85% | 409 / 5.28% | +0.57 pp | `0.119` |
| High threshold | 103 / 1.330% | 98 / 1.265% | +0.065 pp | `0.772` |

![FP 與 matched TP 的 endpoint 比較](figures/02_matched_tp_specificity.png)

stable 的 pooled site-weighted 差異雖達名義顯著，但不能解讀為 FP-specific biological effect：

1. HCC1395 的 AF / DP / GQ matching imbalance 最大，且 stable RD `+4.68 pp`，主導 pooled 結果。
2. 6 biological samples 等權 stable RD 只剩 `+0.113 pp`；方向 4 正、2 負。
3. 移除 HCC1395 與 Dorado 後，stable RD 變成 `-0.153 pp`，`p=0.780`。
4. high-threshold 的 6-sample 等權 RD 為 `-0.111 pp`；方向 2 正、1 平、3 負。
5. 在雙方都 evaluable 的 `4,738` pairs 中，primary stable 為 FP `572/4,738 = 12.07%` vs TP `469/4,738 = 9.90%`，paired `p=0.000599`；此名義顯著結果必須完整揭露。
6. 對同一 both-evaluable estimand，6 biological samples 等權 RD 降為 `+0.419 pp`；排除 HCC1395 + Dorado 後為 `+0.120 pp`、`p=0.902`，顯示 site-weighted 顯著性不具跨 biological-sample 穩健性。
7. 在同一 `4,738` pairs 中，高門檻仍是 FP `2.132%` vs TP `2.026%`，`p=0.770`。

**反證結論**：無法把高門檻甲基多群當成 truth-FP 特有、或因此當成 FP 內隱藏 subclone 的高 PPV 指標。

---

## 更嚴格 null 與 allele background 將候選縮到 10 sites，但仍不是 subclone call

### 為何 column-permutation null 不足

column-permutation 會打破 read 內部 CpG correlation。為避免把自然的 row-level methylation block structure 當成群，本輪另跑：

- 每條 read 獨立 cyclic shift 的 row-circular null；保留該 read 的 values、missingness 與 covariance。
- `199` null replicates、empirical `p <= 0.05`、至少 `80%` valid null。
- 10 seeds 的 modal-K assignment；minimum pairwise ARI `>=0.8`。

實際存活：

- primary stable：`583`
- column empirical-p pass：`542`
- row-circular pass：`245`
- high threshold：`103`
- high threshold + row-circular + assignment：`39`
- high threshold + tumor-REF background：`30`
- 所有 strict gates 交集：`10`

![嚴格 null 與 allele background 的敏感度結果](figures/04_strict_sensitivity.png)

### tumor REF 背景不是空白

在 `583` 個 ALT-stable sites 中，tumor REF 有 `450` 個可評估，其中 `112/450 = 24.89%` 的 REF reads 本身也形成 stable multigroup。

在 `103` 個高門檻候選中：

- REF 可獨立評估 `76`，其中 `21/76 = 27.63%` REF stable。
- ALT+REF joint clustering 可評估 `78`，其中 `64` joint stable。
- `41` 個與 allele axis 對齊；`23` 個沒有 allele alignment。
- REF stable 或 joint-not-aligned 共削弱 `34/103` 個高門檻候選。

### normal sequence evidence 與 normal methylation background不同

- `43/103` 高門檻候選在 normal 有 ALT read 且 normal AF `>=1%`。
- `31/103` normal AF `>=5%`。
- `29/103` tumor caller AF `>=0.8`。
- 最新 InterSubMod batch **沒有跑 matched-normal methylation matrix**，所以 normal methylation background 仍是 `NOT_TESTED`。

### 10 個 strict epigenetic follow-up candidates

| Dataset | Site | ALT reads | Methyl group sizes | Tumor AF | Normal AF | LongPhase-S | Regional topology context |
|---|---|---:|---|---:|---:|---|---|
| COLO829 | chr4:92,183,865 C>T | 12 | 7 / 5 | 0.500 | 0 | PASS->PASS | linear regional |
| COLO829 | chr13:54,250,213 C>T | 15 | 9 / 5 | 0.484 | 0 | PASS->PASS | incomplete |
| COLO829 | chr20:52,761,564 C>T | 21 | 12 / 8 | 0.512 | 0 | PASS->PASS | linear regional |
| COLO829 | chr20:52,761,565 C>T | 21 | 12 / 8 | 0.512 | 0 | PASS->PASS | linear regional |
| H1437 | chr1:85,105,364 A>C | 35 | 21 / 14 | 0.468 | 0 | LowQual->PASS | branching / recurrence |
| HCC1395 | chr8:82,037,650 A>G | 14 | 6 / 5 | 0.149 | 0 | PASS->PASS | incomplete |
| HCC1395 | chr8:82,044,380 A>T | 12 | 7 / 5 | 0.180 | 0 | PASS->PASS | no regional tree |
| HCC1395 | chr8:82,160,837 T>C | 16 | 9 / 7 | 0.185 | 0 | PASS->PASS | no regional tree |
| HCC1395 | chr8:82,170,588 G>A | 17 | 12 / 5 | 0.212 | 0 | PASS->PASS | no regional tree |
| HCC1395 Dorado | chr5:750,311 A>G | 19 | 8 / 5 / 5 | 0.196 | 0 | LowQual->PASS | incomplete |

chr20 兩個相鄰 sites 使用相同 ALT readset，因此 10 sites 只有 9 unique ALT readsets。再依 latest layered reconstruction 的 `component_id` 聚合，只剩 **6 個 regional components**：HCC1395 四個 chr8 sites 同屬 `chr8:81975567-82178585`，COLO829 兩個 chr20 sites 同屬 `chr20:52680060-52761565`。這些位點的適當用途是 targeted follow-up，不是 10 個獨立事件，也不是直接輸入 clone tree。

---

## LongPhase-S tags 完整讀取，但沒有第二 somatic haplotype 證據

所有 `7,745` FP 的 focal ALT reads 共 `113,845` 條：

| Scope | ALT reads | HP 1-1 | HP 2-1 | HP3 ambiguous | Untagged | HP 1-2 / 2-2 |
|---|---:|---:|---:|---:|---:|---:|
| All FP | 113,845 | 53,552 | 53,444 | 6,560 | 289 | **0** |
| Stable | 15,412 | 6,897 | 7,857 | 620 | 38 | **0** |
| High threshold | 3,276 | 1,616 | 1,655 | 4 | 1 | **0** |
| Strict | 182 | 147 | 34 | 1 | 0 | **0** |

![LongPhase-S HP tags 在 focal ALT reads 的分布](figures/06_longphase_hp_tags.png)

`1-1/2-1` 都屬 LongPhase-S 第一 somatic haplotype family；`1-2/2-2` 才是明確第二 somatic haplotype。全體為 `0`，因此不能用 HP tag 說「ALT 內存在兩個 somatic haplotypes」。高門檻 `103/103` assignment ARI 穩定，只能說 methyl groups 穩定，不能把 methyl group 改名為 HP-defined clone。

---

## Topology 反證：相同 methylation endpoint 同時出現在 linear 與 branching regional context

### 全部高門檻 103 sites

| Regional context | 數量 |
|---|---:|
| No regional tree | 28 |
| Branching / recurrence | 25 |
| Candidate set incomplete | 22 |
| All stored regional trees linear | 22 |
| Trivial one-edge only | 6 |

- 高門檻候選相對其他 FP 並未富集 linear regional context：Fisher OR `1.217`、`p=0.441`。
- focal single-site linear topology identifiable：`0/7,745`。
- focal single-site subclone confirmed by layered tree：`0/7,745`。
- strict 10 sites 仍分成 linear `3`、branching `1`、incomplete `3`、no tree `3`。

![Methylation endpoint 橫跨 linear 與 branching regional context](figures/05_topology_context.png)

**結構性理由**：一個 focal sSNV 最多提供 `ROOT -> ALT` 的 trivial edge。要判斷 linear 或 branching，至少需要第二個可比較 mutation state，以及兩者在同一 reads／cells 的共現、互斥與 cellular prevalence。Same-pipeline regional tree 只能作 context，不是獨立正交確認。

![HCC1937：同一高門檻 endpoint 分別出現在 linear 與 branching context](figures/10_topology_contrast_heatmaps.png)

左例 HCC1937 chr9:23,323,391 是 regional linear context；右例 chr14:58,925,240 是 branching/recurrence context。兩者都可形成清楚兩群，但 topology 相反，直接否定「看到兩群就可假設 linear」的敘述。

---

## 技術重複：兩個 exact FP site 可重現，但 cluster 數與 HP orientation 不完全固定

HCC1395 與 HCC1395 Dorado 有 `204` 個共同 FP sites，`128` 個雙方可評估：

| Endpoint | HCC1395 positive | Dorado positive | Both | Either | Jaccard |
|---|---:|---:|---:|---:|---:|
| Stable | 7 | 12 | 2 | 17 | 0.118 |
| Residual | 6 | 9 | 2 | 13 | 0.154 |
| High threshold | 2 | 2 | 2 | 2 | 1.000 |

兩個高門檻 exact sites：

1. `chr1:175,089,892 C>T`：兩邊都是 2 methyl groups；HCC1395 focal ALT 全為 HP `1-1`，Dorado 全為 HP `2-1`。
2. `chr9:39,585,394 G>C`：HCC1395 是 2 groups、regional linear context；Dorado 是 3 groups、trivial one-edge context。

![HCC1395 與 Dorado 的實際 read-level methylation heatmaps](figures/09_technical_replication_heatmaps.png)

這支持「局部 epigenetic pattern 可跨 basecaller 技術資料重現」，但不支持「cluster label、cluster 數或 HP orientation 等同固定 clone identity」。而且兩個 dataset 仍是同一 biological sample。

---

## 個案反證：清楚的 heatmap 仍可能是 HP / 技術軸，strict candidate 也可能落在 branching context

![HP/技術混雜與 branching strict candidate 的實際 heatmaps](figures/11_confounder_heatmaps.png)

### HCC1937 chr2:143,925,487 G>A：HP-axis confound

- ALT reads `89`；群大小 `86/3`。
- ALT HP：untagged `19`、HP1-1 `18`、HP2-1 `39`、HP3 `13`。
- `hp_axis_confound=true`，所以 stable 但不是 residual/high-threshold。
- tumor AF `1.0`、normal AF `0.143`，也是不能稱為 tumor subclone 的強反證。

### H2009 chr1:244,024,747 T>C：technical-axis confound

- ALT reads `59`；4 groups `37/4/3/15`。
- `technical_axis_confound=true`，stable 但 residual=false。
- tumor AF `1.0`、normal AF `0.143`。

### H1437 chr1:85,105,364 A>C：strict candidate，但 regional context branching

- ALT reads `35`；兩群 `21/14`；tumor AF `0.468`、normal AF `0`。
- 通過 strict gates。
- regional context 是 branching/recurrence；它本身不能指定 linear topology。

**案例結論**：視覺上很乾淨的兩群或多群不等於生物 clone。必須把反例、背景與 topology 一起看。

---

## 先前 InterSubMod 數據與圖片可作歷史驗證，但舊比例不可直接沿用

latest 與 legacy canonical 共有 `3,302` FP sites：

- focal ALT readset：`3,302/3,302` 相同。
- `reads.tsv` 檔案 byte-level 完全相同：`3,016/3,302 = 91.34%`；其餘 `286` 個檔案含 HP/tag 欄位差異。亦即 ALT read **成員集合**相同，不代表整份 reads metadata 相同。
- methylation matrix：`3,302/3,302` 相同。
- Bernoulli distance matrix：`2,081/3,302 = 63.0%` 完全相同。
- Stable endpoint Jaccard：`0.548`。
- High-threshold endpoint Jaccard：`0.331`。

舊版 canonical stable `687` 中有 `374` sites 至少一個 passed node 的 `n_valid_null=0`；這是 fail-open 風險。新版採 fail-closed 與明確 missing-distance contract。

![Legacy 與 latest 的共同位點 sensitivity 比較](figures/08_legacy_compatibility.png)

因此：

- **可以參考**：先前 heatmap 是否呈現相同 read-level pattern、同一位點是否重現、相同 methylation matrix 的視覺個案。
- **不可合併**：legacy `687/3,188 evaluable = 21.55%` stable prevalence 與 current `583/4,967 = 11.74%`；legacy 分子分母來源為 `InterSubMod/research/20260715_single_fp_alt_multicluster_subclone_validation/results/focal_alt_multicluster/canonical_reference_v1/canonical_summary.json`，但兩版 distance/null semantics 不可交換。
- 歷史來源：`InterSubMod/output/canonical/` 與 `InterSubMod/research/20260715_single_fp_alt_multicluster_subclone_validation/results/focal_alt_multicluster/canonical_reference_v1/`。

---

## 為何文獻也不允許把單一 methylation multigroup 直接等同 subclone

1. [PhyloWGS 原始論文](https://link.springer.com/article/10.1186/s13059-015-0602-8)以大量 somatic mutations 的 VAF、copy-number 與 population frequency 重建 subclonal composition；subpopulation 通常由多個共同 mutations 定義，而非單一位點的表觀群。
2. [RETrace 原始論文](https://genome.cshlp.org/content/30/4/602)同時量測 lineage marker 與 DNA methylation，但 lineage 由 microsatellite mutations建立、methylation 用於 cell-state/type；兩者不是同一證據。
3. [Li et al., Nature Medicine 2016](https://pubmed.ncbi.nlm.nih.gov/27322744/)顯示 genetic 與 epigenetic heterogeneity 可有不同 kinetics，支持把 epigenetic multigroup 與 genetic clone 分開命名。
4. [Costabile et al., Biology Direct 2026](https://pubmed.ncbi.nlm.nih.gov/42231478/)觀察單一 clone expansion 後仍可重新形成多樣 epialleles；這提供一個直接替代機制：多 epialleles 不必來自多 genetic clones。
5. [LongPhase 官方 repository](https://github.com/twolinin/longphase)明確區分 somatic paired (`longphase-s`) 與 tumor-only (`longphase-to`) 工具；本報告的 tag semantics 另以本地 KB `05_tools/longphase-s.md` 核對。

因此，本研究最合理的新穎性不是「用一個 FP 發現 clone」，而是建立 **read-linked epigenetic heterogeneity candidate layer**，再與多 sSNV / CN / CCF / single-cell evidence 做整合。

---

## 哪些額外證據到位後，才可以把 candidate 升級為 subclone-supported

### 最低驗證階梯

| Gate | 必要證據 | 目前狀態 |
|---|---|---|
| A. Read-level pattern | stable null + assignment ARI + depth sensitivity | 已完成 |
| B. Allele specificity | tumor REF + matched-normal methylation background | tumor REF 已做；normal methylation 未做 |
| C. Genetic co-membership | 至少 2 個獨立 sSNV/CN markers 在同一 read group 共現 | 未做 |
| D. CCF consistency | purity/CN-adjusted cancer-cell fraction 與 group fraction相容 | 未做 |
| E. Technical replication | independent library/basecaller/platform 重現 group membership | HCC 技術重複只部分完成 |
| F. Cellular identity | single-cell DNA/methylation、colony、spatial 或 lineage barcode | 未做 |
| G. Topology | 多 mutation states 支持 ancestor/descendant 或 branching constraints | focal site `0` identifiable |

### 建議判定規則

- 通過 A+B：稱 **epigenetic heterogeneity candidate**。
- 再通過 C+D，且至少 E/F 之一：才稱 **subclone-supported candidate**。
- 再通過 G：才討論 **linear / branching topology**。

### 演化樹輸入應該使用什麼

樹的合理 genetic input 應繼續使用 **同一次最新 LongPhase-S 輸出的 recalibrated `FILTER=PASS` sSNV**，再整合 multi-locus read co-occurrence、CN/LOH、purity-adjusted CCF 與不確定性。單一 sSNV 的 methyl group 不能取代 sSNV callset，也不能直接當 tree node；strict candidates 可作 annotation 或安排 targeted linkage experiment。

---

## 多 agent 共識：四個獨立審查面向都拒絕 overclaim

完整紀錄：`InterSubMod/research/20260715_single_fp_alt_multicluster_subclone_validation/reviews/20260715_multiagent_consensus.md`

| Reviewer | 獨立結論 |
|---|---|
| Provenance | 最新 PASS、sidecar tags、MM/ML 與 no-overwrite 成立；HP tag 不等於 clone label |
| Statistics | cluster assignment 可穩定，但高門檻不具 FP specificity；strict 10 仍無 FDR/PPV |
| Topology | 一個 sSNV 無法識別 linear/branching；regional context 不是正交確認 |
| Legacy | 舊圖可驗證現象；舊 distance/null 口徑使舊率不可直接比較 |

四位 reviewer 無實質分歧：支持保留「epigenetic heterogeneity candidate」價值，但不支持「high-probability subclone」或「linear subclone」。外部 Claude Code 第二輪 read-only audit 的 verdict 為 **VALID WITH CORRECTIONS**；15 項獨立重算全部 PASS，未發現會翻轉 verdict 的數據錯誤。其三項透明度修正（both-evaluable stable、legacy 分母、readset vs reads.tsv）均已納入本版。完整紀錄：`InterSubMod/research/20260715_single_fp_alt_multicluster_subclone_validation/reviews/20260715_external_claude_code_review.md`。

---

## 限制與需更完整觀察之處

1. **Normal methylation background 未跑**：目前只有 normal AD/AF，不等於 matched-normal epigenetic control。
2. **CN / LOH / purity-adjusted CCF 未整合**：尤其 caller AF 接近 1 或 normal AF >0 的位點。
3. **沒有獨立 cellular identity**：orthogonal confirmation 結構上是 `0`，不是「已測且全陰性」。
4. **HCC 技術重複非獨立 biological sample**：不能當跨個體 replication。
5. **COLO829 無 confident-region BED**：其 truth-FP label confidence 必須降級。
6. **site dependency**：shared reads、相鄰 sSNVs、相同 readsets 使 site-level CI 偏樂觀。
7. **strict 不是 genome-wide FDR**：`10` 是 sensitivity intersection，不是宣稱 10 個真 subclones。
8. **cluster threshold 未完全校準**：min group size、modal fraction 與 null mode 都會改變 candidate 數。
9. **同 pipeline topology 非正交**：可作 context，不能獨立證明 methyl group 的細胞 lineage。

---

## 最終結論與下一個合理實驗

### 最終結論

**不支持原命題**：「單一 truth-FP sSNV 的 ALT reads 被分成多個 methylation groups」本身不是 high-probability subclone 證據，也不支持 linear evolution。

**支持縮限命題**：在所有最新 LongPhase-S PASS FP 中，約 `11.74%` evaluable sites 有 primary null-gated methylation multigroup；高門檻為 `2.07%`，但 matched TP 幾乎相同；經 strict sensitivity 後只剩 `0.201%` evaluable FP，適合作為正交驗證優先清單。

### 下一個最高資訊量實驗

以 10 strict sites 為 anchor，但不只分析單一位點：

1. 在每個 anchor 周圍擷取同一 ALT reads 跨越的其他 latest LongPhase-S PASS sSNVs、germline phase markers 與 CN/LOH segment。
2. 對 methyl group A/B 比較 multi-locus somatic co-occurrence；預註冊至少 2 個獨立 genetic markers。
3. 加入 matched-normal BAM 的同 region methylation matrix，使用相同 row-circular null。
4. 在 HCC1395 四個 chr8 strict sites 檢查是否形成區域共享 read component／相同 genetic state；避免把區域 epigenetic domain 當四個 clones。
5. 只有通過 genetic co-membership + CCF consistency 的 candidates 才送入 topology inference。

---

## 可重現命令、輸入與輸出

### 嚴格 null / REF survival

```bash
env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
    NUMEXPR_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 \
python3 research/20260715_single_fp_alt_multicluster_subclone_validation/scripts/audit_strict_null_assignment_survival.py \
  --site-results research/20260715_single_fp_alt_multicluster_subclone_validation/results/focal_alt_multicluster/latest_full_v3_frozen/latest_site_results.tsv \
  --ref-background research/20260715_single_fp_alt_multicluster_subclone_validation/results/focal_alt_multicluster/ref_background_v2/ref_background_site_results.tsv \
  --workers 42 --seeds 10 --rnull 199 --empirical-alpha 0.05 \
  --output-dir research/20260715_single_fp_alt_multicluster_subclone_validation/results/focal_alt_multicluster/strict_survival_v1
```

實際輸出：`processed=583/583 failures=0`；strict `10`；exit code `0`。

### Frozen topology annotation

```bash
python3 research/20260715_single_fp_alt_multicluster_subclone_validation/scripts/annotate_layered_topology_context.py \
  --site-results research/20260715_single_fp_alt_multicluster_subclone_validation/results/focal_alt_multicluster/latest_full_v3_frozen/latest_site_results.tsv \
  --output-tsv research/20260715_single_fp_alt_multicluster_subclone_validation/results/focal_alt_multicluster/latest_full_v3_frozen/latest_site_results_with_topology.tsv \
  --summary research/20260715_single_fp_alt_multicluster_subclone_validation/results/focal_alt_multicluster/latest_full_v3_frozen/latest_topology_context_summary.json
```

實際輸出：7 datasets / `7,745` sites；`pass=true`；exit code `0`。

### Matched-control robustness

```bash
python3 research/20260715_single_fp_alt_multicluster_subclone_validation/scripts/analyze_fp_tp_robustness.py \
  --comparison-summary research/20260715_single_fp_alt_multicluster_subclone_validation/results/focal_alt_multicluster/fp_matched_tp_comparison_v1/fp_matched_tp_comparison_summary.json \
  --output-dir research/20260715_single_fp_alt_multicluster_subclone_validation/results/focal_alt_multicluster/fp_matched_tp_robustness_v1
```

實際輸出：summary + leave-one-out TSV；`pass=true`；exit code `0`。

### Report SoT dataset 與 figures

```bash
python3 research/20260715_single_fp_alt_multicluster_subclone_validation/scripts/build_report_dataset.py \
  --output-dir research/20260715_single_fp_alt_multicluster_subclone_validation/results/report_dataset_v1

python3 research/20260715_single_fp_alt_multicluster_subclone_validation/scripts/generate_report_figures.py \
  --report-dataset research/20260715_single_fp_alt_multicluster_subclone_validation/results/report_dataset_v1/report_dataset.json \
  --assignments research/20260715_single_fp_alt_multicluster_subclone_validation/results/focal_alt_multicluster/latest_full_v3_frozen/latest_stable_multigroup_read_assignments.jsonl \
  --output-dir research/20260715_single_fp_alt_multicluster_subclone_validation/figures
```

實際輸出：report dataset `pass=true`、FP `7,745`、strict `10`；11 PNG manifest `pass=true`。

### Standalone HTML 與 browser QA

```bash
python3 research/20260715_single_fp_alt_multicluster_subclone_validation/scripts/qa_standalone_html.py \
  --html research/20260715_single_fp_alt_multicluster_subclone_validation/20260715_單一FP_sSNV_ALT_read甲基多群與subclone假說全量驗證_01.standalone.html \
  --output-dir research/20260715_single_fp_alt_multicluster_subclone_validation/qa

python3 -m unittest discover \
  -s research/20260715_single_fp_alt_multicluster_subclone_validation/tests \
  -p 'test_*.py' -v
```

實際輸出：desktop `1440x1000` 與 mobile `390x844` 均 `pass=true`；document width 等於 viewport、11/11 圖片載入、0 console/page errors、0 文字元件 overflow、目錄與 `<details>` 操作通過；14/14 unit tests PASS，exit code `0`。

---

## 關鍵資料位置與 provenance

| 類型 | 路徑 |
|---|---|
| Machine-readable report SoT | `InterSubMod/research/20260715_single_fp_alt_multicluster_subclone_validation/results/report_dataset_v1/report_dataset.json` |
| Per-sample metrics | `InterSubMod/research/20260715_single_fp_alt_multicluster_subclone_validation/results/report_dataset_v1/per_sample_metrics.tsv` |
| Strict candidate table | `InterSubMod/research/20260715_single_fp_alt_multicluster_subclone_validation/results/report_dataset_v1/strict_followup_candidates.tsv` |
| Actual case catalog | `InterSubMod/research/20260715_single_fp_alt_multicluster_subclone_validation/results/report_dataset_v1/case_catalog.tsv` |
| Latest frozen FP results | `InterSubMod/research/20260715_single_fp_alt_multicluster_subclone_validation/results/focal_alt_multicluster/latest_full_v3_frozen/` |
| Matched TP frozen results | `InterSubMod/research/20260715_single_fp_alt_multicluster_subclone_validation/results/focal_alt_multicluster/matched_tp_v2_frozen/` |
| Strict audit | `InterSubMod/research/20260715_single_fp_alt_multicluster_subclone_validation/results/focal_alt_multicluster/strict_survival_v1/` |
| Multi-agent consensus | `InterSubMod/research/20260715_single_fp_alt_multicluster_subclone_validation/reviews/20260715_multiagent_consensus.md` |
| External Claude Code audit | `InterSubMod/research/20260715_single_fp_alt_multicluster_subclone_validation/reviews/20260715_external_claude_code_review.md` |
| Figure manifest | `InterSubMod/research/20260715_single_fp_alt_multicluster_subclone_validation/figures/figure_manifest.json` |
| Standalone HTML | `InterSubMod/research/20260715_single_fp_alt_multicluster_subclone_validation/20260715_單一FP_sSNV_ALT_read甲基多群與subclone假說全量驗證_01.standalone.html` |
| Browser QA receipt | `InterSubMod/research/20260715_single_fp_alt_multicluster_subclone_validation/qa/browser_qa_summary.json` |
| Cross-session handoff | `InterSubMod/research/20260715_single_fp_alt_multicluster_subclone_validation/20260715_跨session研究交接與資料位置_01.md` |

Frozen hashes：

- `analyze_focal_alt_multicluster.py`: `c243ee9db03b4ccae4a97f0740b1994012c01d29be4701b1de01d987f5235125`
- `focal_alt_cluster_lib.py`: `da2aa67b4c8991fb7ffdd4ba77e14878c799b730d9ec83d1e909b8b2b248ee54`
- `latest_site_results_with_topology.tsv`: `479c5511745cb66337757401d96fe3a655233e0afeb7ee1f803f35a60d25c0f1`
- Git commit at report dataset build: `2bec873b40a668e6227cde8f68b1432796398e93`；working tree dirty，故以上 file hashes 是必要 provenance。

---

## 一句話帶走

**把這些結果稱為「少數、需正交驗證的 focal-ALT read-level epigenetic heterogeneity candidates」是可信的；把它們直接稱為 subclones 或 linear evolution，目前不可信。**
