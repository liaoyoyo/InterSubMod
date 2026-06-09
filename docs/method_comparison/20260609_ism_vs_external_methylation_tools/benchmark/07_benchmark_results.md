<!--
建立時間: 2026-06-10
狀態: in_progress (Phase B benchmark — 第一批 head-to-head 結果, 持續擴充)
報告類型: benchmark_results
受眾: PI · 方法驗證 · 論文 Methods
framework: 結果表 + 跨工具 concordance + 軸對齊 + caveat
data_sources:
  - runs/dmr_somatic.tsv, runs/dmr_family.tsv (modkit 0.6.3 實跑, 2026-06-10, 已 Read 讀回)
  - ISM 驗證值: knowledge/11_external_literature/07 (BRCA2 Δβ=-0.122, TBC1D16 d_within=+0.142)
provenance_note: modkit 數字皆自 runs/dmr_*.tsv 實跑讀回(非記憶); ISM 數字引 validated KB。本檔 Write 與 dmr 分析在不同 batch(§13.0)。單 subset(BRCA2+TBC1D16, 400 reads)第一批, 非全 panel。
-->
<!-- provenance-verified: modkit effect_size/score 引自 runs/dmr_somatic.tsv+dmr_family.tsv (modkit 0.6.3 實跑讀回); ISM Δβ 引 knowledge/11_external_literature/07 validated。產數字(dmr)與本 Write 不同 batch。 -->

# 07 — Phase B benchmark 結果（第一批 · modkit vs ISM）

> **狀態**：Phase B 進行中。本批 = **第一個 head-to-head 數據點**（modkit dmr 實跑 vs ISM 驗證值），位點 = BRCA2 + TBC1D16。DSS/ASMS 安裝中、cvlr repo 404、全 panel 待擴充。

---

## L0 — 一句話結果

**modkit（獨立工具，Bayesian pileup 率差）獨立重現了 ISM 的 BRCA2 somatic 軸低甲基化**：modkit `effect_size = −0.159`（5mC HP1 30.6% → HP1-1 14.6%）vs ISM 驗證 `Δβ = −0.122` —— **同方向、同量級** = 跨工具交叉驗證。**但兩者都只給「率差」，都未分離 copy**；只有 ISM 的 normal-anchored cis-test 進一步隔出邊際真 cis（d_within=−0.023）—— **這正是 modkit 沒有、ISM 獨有的能力**（對齊 03/05 結論）。

---

## L1 — 結果表（modkit 實跑 ↔ ISM 驗證值）

| 位點 | 軸 | modkit a%(HP1) | b%(HP1-1) | **modkit effect_size** | 5mC HP1→HP1-1 | modkit score(LR) | **ISM 驗證值** | concordance |
|------|----|---------------|-----------|----------------------|---------------|------------------|---------------|-------------|
| **BRCA2** | **somatic**(HP1 vs HP1-1) | 31.56% | 15.67% | **−0.159** | 30.56→14.64 | 145.6 | **Δβ=−0.122** | 🟢 同方向(hypo)+同量級 |
| **TBC1D16** | somatic | 77.81% | 86.57% | +0.0876 | 75.61→83.59 | 70.3 | d_within=+0.142 | 🟢 同方向(hyper) |
| BRCA2 | family(HP1 vs HP2) | 24.00% | 22.78% | −0.0122 | — | 1.69 | — | family 軸幾乎平 |
| TBC1D16 | family | 79.87% | 44.86% | −0.350 | — | 107.1 | — | ⚠HP2 cov 低(370)=cnLOH |

> `effect_size` = modkit 的 b_pct−a_pct（HP1-1 − HP1 的 C-modification 比例差）；`score` = Bayesian marginal-likelihood ratio（越高越 differential）；a/b% = 總 C-modification%（m+h），5mC 欄是分 code 的 m%。

### ISM 側 validated 分層值（同位點，取自 cis_scan_full.json + control_cohesion_cistest.json）

| 位點 | ISM dbeta_HP | dbeta_5mC | d_cis (p_cis) | d_drift | silhouette HP1-1 | ari_blind_allele |
|------|-------------|-----------|---------------|---------|-----------------|------------------|
| **BRCA2** | −0.122 | **−0.121** | −0.142 (p=0.0) | −0.021 | 0.313 | 0.79 |
| **TBC1D16** | +0.122 | **+0.109** | +0.231 (p=0.0) | +0.110 | — | — |

> 對齊 modkit 5mC：BRCA2 modkit **−0.159** ↔ ISM dbeta_5mC **−0.121**（同向同量級）；TBC1D16 modkit **+0.080** ↔ ISM **+0.109**。**率差層收斂**。
>
> 📐 **驗證新發現**：modkit dmr **原生輸出 Cohen's h + 95% CI**（dmr tsv 的 `cohen_h/low/high` 欄）：BRCA2 somatic **cohen_h=−0.38 [0.34–0.42]**、TBC1D16 **+0.23 [0.19–0.27]** —— 即改進 **#4「ISM 應加 Cohen's h(arcsine)」對手已原生提供**，可直接照搬。
>
> ✅ **2026-06-10 數據確認**：modkit effect_size 算術逐列核驗（=b_pct−a_pct）OK；ISM BRCA2 cis-test 基於 **n_shared_cpg=194**（穩健，非單點）；覆蓋充足，唯 TBC1D16 family HP2 b_total=370 偏低(cnLOH)。

---

## L2 — 三個解讀

### ① 跨工具交叉驗證 BRCA2 somatic 低甲基化 🟢
- modkit（pileup 計數率差 + Bayesian LR）與 ISM（per-CpG 配對 Wilcoxon Δβ）**獨立**得到 BRCA2 HP1-vs-HP1-1 **同方向、同量級**的低甲基化（−0.159 vs −0.122）。
- 兩者方法完全不同（modkit 聚合 per-position、ISM 保留 read×CpG 配對），卻收斂 → **這個 finding 不是單一工具 artifact**。
- 差異（−0.159 vs −0.122）來自：modkit = pooled count fraction（combine-strands、全 read 聚合、m+h）；ISM = per-CpG 配對均差（max-collapse any-mod，n=197 CpG）。量級一致。

### ② 「軸」決定性地影響結果 🔑
- BRCA2 **somatic 軸（HP1 vs HP1-1）= −0.159** vs **family 軸（HP1 vs HP2）= −0.012** —— 差 13×。
- 證實報告 03/04 的核心：**modkit 原生只能做 HP1 vs HP2（須先 pileup --phased）**；要得到 somatic 軸必須像我們一樣**重標 HP1-1**（本 benchmark 用 pysam 把 HP:Z `1-1`→HP:i `2` 才讓 modkit 跑出 somatic 軸）。ISM 原生就是 HP1 vs HP1-1。

### ③ modkit 缺 cis/copy 控制（ISM 獨有能力的實證）
- modkit effect_size −0.159 = **raw 率差，混入 copy**（無 normal-anchor、無 copy-partition）。
- ISM 在同位點：HP-axis Δβ=−0.122（~80% copy）→ normal-anchored cis-test 隔出 **d_within=−0.023（邊際真 cis）**。
- **modkit 無對應功能** → 若只看 modkit，會把 BRCA2 的 −0.159 當「強 ASM」；ISM 的 cis-control 才揭露大部分是 copy。**這正是 05 KEEP #2「normal-anchored cis-control 唯一」的實機證據。**

### ④ 分層分解 — modkit 只到第 1 層，ISM 再往下三層（BRCA2 實證）

| 層 | 問題 | modkit | ISM | BRCA2 值 |
|----|------|--------|-----|---------|
| **L1 率差** | HP1-1 vs HP1 甲基率差多少 | ✓ effect_size | ✓ dbeta | modkit −0.159 ≈ ISM −0.121（**收斂**）|
| **L2 normal-anchor cis vs drift** | 這率差是 cis-driven 還是 germline drift | ✗ | ✓ d_cis vs d_drift | d_cis −0.142 > d_drift −0.021 |
| **L3 copy-partition** | 扣掉 copy 後還剩多少真 cis | ✗ | ✓ d_within | **−0.023（~20% 真 cis，~80% copy）** |
| **結構** | reads 是否分成結構分離亞群 | ✗ | ✓ silhouette/ARI/PERMANOVA | HP1-1 silhouette 0.313；ari_blind 0.79 |

> 🔑 **一句話**：modkit 與 ISM 在 **L1 率差層收斂**（互相驗證 BRCA2 finding 為真），但 modkit **止步於此**；ISM 的 L2/L3/結構層揭露「−0.122 的率差其實 ~80% 是 copy、只 ~20% 真 cis」—— 這是「ISM 的組合是否帶來對手沒有的真資訊」的**正面實證**（synth gap 對這 2 個 concordant 位點的回答 = YES，ISM 多給了 cis/copy 分解；discordant case 仍待 panel 擴充）。

---

## L3 — 方法 + 可重現指令

**輸入**：HCC1395 longphase_s tagged BAM（HP:Z `1`/`1-1`/`2`）→ subset 到 loci（BRCA2 chr13:32,313,628-32,316,628；TBC1D16 chr17:79,989,620-79,992,620）= 400 reads。
**關鍵修正**（踩到的坑）：
1. modkit `--phased` 只認 `HP:i` 整數 → 用 `runs/retag_hp.py` 把 `HP:Z`→`HP:i`（somatic: 1→1,1-1→2；family: {1,1-1}→1,{2,2-1}→2）。
2. modkit 0.6.3 pileup 需 `--modified-bases C:m C:h`（不自動偵測）。
3. modkit dmr 指定 mod base 用 **`--base C`**（非 `--single-code`/`--assign-code`，後兩者報 "need to specify at least 1 modified base"）。
```bash
modkit pileup --phased --cpg --combine-strands --modified-bases C:m C:h -r $REF HCC1395_somatic_hp.bam pileup_somatic/
# bgzip+tabix hp1/hp2.bedmethyl
modkit dmr pair -a pileup_somatic/hp1.bedmethyl.gz -b pileup_somatic/hp2.bedmethyl.gz \
  --ref $REF -r loci.bed -o dmr_somatic.tsv --header --base C
```

---

## 持續驗證 — 待補（下一步）
- [ ] **discordant-case 測試（核心）**：在同 subset 重跑 **ISM PERMANOVA**，找「modkit 率差≈0 但 ISM PERMANOVA 顯著」的位點 → 回答 synth gap「ISM 結構訊號是真還噪」。
- [ ] **DSS**（安裝中）：用 beta-binomial+shrinkage 重算同位點 Δ，量化 ISM Fisher 是否在低覆蓋膨脹（05 #1）。
- [ ] **ASMS**（Raineri 2024，最像 ISM）：repo 待定位；no-phasing read-clustering 對照。
- [ ] **cvlr**：GitHub 404（作者 0 repo）→ 暫無法裝，紙上對照。
- [ ] 擴 panel（+imprinting ICR 正控 + control loci）+ runtime 比較。

## Provenance
- modkit 0.6.3（gh release binary `tools/dist_modkit_v0.6.3_26c3f9e/modkit`）；輸出 `runs/dmr_somatic.tsv` / `runs/dmr_family.tsv`（實跑讀回）。
- ISM 驗證值：`knowledge/11_external_literature/07`（BRCA2 Δβ=−0.122、d_within=−0.023；TBC1D16 d_within=+0.142）。
- ⚠ 單 subset 第一批；modkit=pileup 率差(無 copy/cis 控制)；ISM 值取自全基因組 validated 分析(非同 subset 重跑)→ 嚴格 apples-to-apples 需下一步同 subset 重跑 ISM。
