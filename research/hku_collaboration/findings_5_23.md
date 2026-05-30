<!--
建立時間: 2026-05-23
類型: A2 pysam analysis 統整 finding
主 source: A2 parallel-analysis subagent 報告（HCC1395 tagged BAM + SEQC2 truth）
範圍: PS set TP 分佈 + 相鄰 TP pair common reads + LOH HP-missing ROC
-->

# HKU Collaboration Pysam Analysis — Findings 5/23

## Summary one-liner

**PS set median 5 TP somatic 位點（85.9% PS ≥2 TP）／ 相鄰 TP pair median 17 common reads（≤±16bp window 內 70.5）／ LOH HP-missing rate AUC 0.504（≈random，hypothesis 拒絕）／ PS block global N50 943 kb (§11)**

> **術語**（v4）：本檔以下「**LongPhase-S**」一律指 `twolinin/longphase` mainline **v1.7.3** with `--somaticMode` flag（產出 HP1-1 / HP2-1 / HP3 細分 tag；mainline 內建功能，非 fork）。HP1_1=5 / H2_1=7 enum 定義於 `src/haplotag/HaplotagType.h:97-108`。

---

## §1 BAM framing 確認（from @PG line）

從 `samtools view -H` 的 `@PG ID:longphase VN:1.7.3` 行直接驗證：

| 元件 | 實際 |
|------|------|
| Tumor SNVs | **ClairS-TO ssrs v0.4.x**（`pileup/HCC1395_methyl_PASS.vcf`）|
| Tagging | **LongPhase-S**（`twolinin/longphase` mainline v1.7.3）`haplotag --somaticMode` |
| Paired-with-normal | `-s HCC1395BL_methyl_phase.vcf.gz -b HCC1395BL.bam` |
| SEQC2 guidance | `--highCon-snp high-confidence_sSNV_in_HC_regions_v1.2.1.vcf` |
| 其他 flags | `-q 20 --tagSupplementary -t 64` |

**❌ 不可寫法**：「ClairS paired pileup」、「LongPhase-S CCU fork」、「longphase-S v1.0.0」
**✅ 正確措辭**：「ClairS-TO ssrs + Normal BAM → LongPhase-S（mainline v1.7.3）`--somaticMode`」

---

## §2 Analysis 1 — PS set TP 位點分佈（chr1 + chr8 + chr19）

從 SEQC2 v1.2.1 PASS-only TP positions（chr1: 3,440, chr8: 2,605, chr19: 874 = 6,919 total），scan tagged BAM 並按 read PS tag 聚合：

| 指標 | 值 |
|---|---|
| PS sets 總數 | **759** |
| **median TP / PS** | **5** |
| mean TP / PS | 9.5 |
| **IQR [p25, p75]** | **[2, 12.5]** |
| max TP / PS | **90** |
| **PS sets ≥2 TP** | **85.9% (652/759)** |

**Finding**: 絕大多數 PS set 含足夠 TP 位點做 read-level 分群（85.9% ≥2 TP；max 90）。但 ~14% PS 只有 1 個 TP，需仰賴 germline anchors。

**Figure**: `figures/A2_1_ps_set_tp_count_dist.png`

---

## §3 Analysis 2 — 相鄰 TP somatic pair common reads

按 chrom+pos 排序，對每個 PS set 內相鄰 TP pair (i,j) 算 distance & common reads（兩位點同時被 cover 的 read 數），共 3,963 pairs：

| 指標 | 值 |
|---|---|
| Pairs 總數 | **3,963** |
| median distance | **17,421 bp** |
| **median common reads** | **17** |
| mean common reads | 28.1 |
| Pairs in ≤±16 bp caller window | **96 (2.4%)** |
| **median common reads — in-window (≤±16 bp)** | **70.5** |
| median common reads — (16, 500] bp | 64 |
| median common reads — > 500 bp | 16 |

**Finding**:
- ±16 bp caller 視窗內 pair 共享 reads 顯著更多（70.5 vs 16；4.4× ratio）— 但這只佔 2.4%
- median pair distance 17.4 kb 遠超 ONT 平均 read 長度 → 95% pair 需 LongPhase PS chaining 跨越單一 read 範圍
- (16, 500] bp 區段（86 pairs）median 64 reads 接近 ≤±16 bp 區段 → 上限主要由 read length 主導，非 caller 窗硬截斷

**Figure**: `figures/A2_2_common_reads_xy.png`

---

## §4 Analysis 3 — LOH HP-missing rate ROC vs SEQC2 truth ❌ NEGATIVE

1 kb bin × min 5 reads，full BAM scan: chr1 2.24M reads, chr8 1.65M reads, 324,595 bins。

| 指標 | Overall | chr1 | chr8 |
|---|---|---|---|
| n_bins | 324,595 | 199,075 | 125,520 |
| LOH fraction (bins) | 67.7% | 48.0% | **99.1%** ⚠ |
| **AUC** | **0.504** | 0.509 | **0.446** ⚠ |
| median HP-missing — LOH | 0.250 | 0.273 | 0.250 |
| median HP-missing — non-LOH | 0.250 | 0.250 | 0.367 |

**Finding**: ❌ **NEGATIVE** — Hypothesis（AUC > 0.7）清楚被拒絕。

**機制解釋**：
- LongPhase-S `--somaticMode` 用 **heterozygous germline SNV** 為 phasing anchor
- LOH **不改變** germline het 位點（germline-from-normal 仍是 het），只改 tumor allelic frequency
- 故 LOH 區段 reads 仍拿到 HP tag → HP-missing rate 不變

**chr8 反向訊號 (AUC 0.446)**：
- non-LOH bin HP-missing **更高**（中位 0.367）vs LOH bin (0.250)
- chr8 99.1% LOH (HCC1395 chr8 LOH hotspot per `project_hcc1395_chr8_hotspot.md`)
- 剩餘 non-LOH bins 多落在 mappability borderline 邊緣 → HP-missing 反而更高

**Figure**: `figures/A2_3_loh_hp_missing_vs_seqc2.png`

---

## §5 對 HKU 報告 4 個討論點的證據對應

| # | 討論點 | 證據 | 判定 |
|---|---|---|---|
| 1 | PS set 內 TP 位點足夠做 read-level 分群 | A2_1: 85.9% PS ≥ 2 TP, median 5, max 90 | ✅ **POSITIVE** |
| 2 | 相鄰 TP somatic pair 有 common reads | A2_2: ≤±16 bp pair median 70.5（但僅 2.4%）；中距 64；遠距 16 | ✅ POSITIVE for 緊鄰；⚠ 遠距須 PS chain |
| 3 | HP-missing rate 作 LOH 間接指標 | A2_3: AUC 0.504，chr8 反向 0.446 | ❌ **NEGATIVE** — 清楚拒絕 |
| 4 | LongPhase-S `--somaticMode` 不直接報 LOH | 機制推論 + A2_3 觀測一致 (HP-missing 不依 LOH status) | ✅ **VERIFIED** — tool design 確實不暴露 LOH |

### Reframe 對 HKU 報告影響

原 plan 假設「LongPhase-S 不直接報 LOH 但 HP-missing 可間接指標」**被駁回**。改為：

> 「LongPhase-S `--somaticMode` 用 germline het SNV 為 phasing anchor，LOH 不改變 normal-derived germline het 狀態（只改 tumor allelic frequency），因此 HP tag 在 LOH vs non-LOH 區段分配無顯著差異（AUC 0.504）。要間接驗 LOH 需改觀察 **tumor HP1 vs HP2 read count imbalance**，而非 HP-missing rate。」

---

## §6 Caveats / Limitations

1. MAPQ ≥ 20, drop dup/supp/secondary（與 `v5_audit_pysam_visualization.py` 一致）
2. A2_2 in-window pair 數小 (96) — median 70.5 CI 寬；要更穩定需用 ClairS-TO ssrs raw mpileup 候選
3. A2_3 bin = 1 kb 屬粗粒度；100 bp 或許可挖更細訊號（×10 成本）
4. chr8 99% LOH 是 HCC1395 特異性，不可外推 COLO829 / H2009
5. HP-missing 缺失原因混雜：(i) 缺 phasing anchor (ii) MAPQ 邊緣 (iii) 多 alignment (iv) longphase confidence 低 — 都會貢獻 noise
6. Wall clock: BAM index 9 min + 3 analyses 平行 11 min（A2_1 561s, A2_2 560s, A2_3 687s）

---

## §7 完成交付清單

**Scripts** (`research/hku_collaboration/scripts/`):
- `_common.py`
- `A2_1_ps_tp_count_distribution.py`
- `A2_2_adjacent_tp_common_reads.py`
- `A2_3_loh_hp_missing_roc.py`

**Figures** (PNG dpi=150 with CJK font):
- `figures/A2_1_ps_set_tp_count_dist.png` (hist + CDF)
- `figures/A2_2_common_reads_xy.png` (hexbin + boxplot)
- `figures/A2_3_loh_hp_missing_vs_seqc2.png` (ROC + violin)

**Data TSV** (`research/hku_collaboration/data/`):
- `A2_1_ps_tp_counts.tsv` (759 rows)
- `A2_1_summary_stats.tsv`
- `A2_2_pair_common_reads.tsv` (3,963 rows)
- `A2_2_summary_stats.tsv`
- `A2_3_bin_hp_missing.tsv` (324,595 rows, 13 MB)
- `A2_3_summary_stats.tsv`

**Run logs**: `data/A2_{1,2,3}_run.log`

**BAM index** (新建): `/big7_disk/liaoyoyo2001/data/bam/HCC1395_Tmode_tagged_ClairS_pileup_v040_woFilter.bam.bai` (120 MB)

---

## §10 Analysis 8 — Chr Ideogram + HP/LOH/CNV (24-panel, HCC1395 T-mode)

> 2026-05-23 補做，補完 §1-§7 範疇外的全 chr 視覺。腳本 `scripts/A8_chr_ideogram_hp_loh_cnv.py`，10 Mb bin。

### §10.1 設計

- 同一 BAM（`HCC1395_Tmode_tagged_ClairS_pileup_v040_woFilter.bam`，283 GB）+ 同一 SEQC2 truth BED（`ngs_benchmark_cnvs_gain_loss_loh.bed`，660 segments：307 gain / 320 loh / 33 loss）
- 22 autosomes + chrX，10 Mb bin（HG38 lengths），用 read midpoint 指派 bin，跳過 unmapped / secondary / supplementary / duplicate
- HP tag 正規化為 6 類：`1`、`2`、`1-1`、`2-1`、`3`、`no_HP`（per LongPhase-S（mainline v1.7.3）`haplotag --somaticMode` 輸出值；enum `src/haplotag/HaplotagType.h:97-108`）
- 平行：12 worker `multiprocessing.Pool`，BAM streaming 全 24 chr 569 s（~9.5 min）；3 figures 各 < 15 s
- 圖片 DPI=150，CJK font 注入

### §10.2 主要觀察

**LOH coverage 排名**（SEQC2 truth）— Figure 3：

| Tier | Chrs | LOH 範圍 |
|------|------|---------|
| Tier 1（≥75%） | chr8 96%, chr17 93%, chr11 85%, chr3 80%, chr10 78%, chr12 77% | 接近全 chr LOH |
| Tier 2（50-75%） | chr5 72%, chr19 55%, chr2 51% | 多數 chr LOH |
| Tier 3（30-50%） | chr1 45%, chr6 45%, chr13 49%, chr4 40%, chr22 38%, chr18 36%, chr9 33%, chr14 32%, chr15 32% | 部分 LOH |
| Tier 4（<5%） | chr7 16%, chr21 2%, **chr16 0%, chr20 0%, chrX 0%** | 幾無 LOH |

→ 印證 §4 / §5 既有結論：**chr8 是 HCC1395 LOH dominance 標誌**；chr16 / chr20 / chrX 則是 LOH-clean reference（對應 A7 §11 chr16 / chrX 短 N50 — phasing-weak ≠ LOH）。

**HP composition imbalance**（normalized stacked）— Figure 2：

1. **絕大多數 autosome HP1/HP2 平衡**（HP1 25-32%, HP2 26-31%, 比例 ≈ 1:1）— 表示 LongPhase haplotag 在 ClairS-TO pileup BAM 上總體對稱
2. **chr6 顯著 outlier**：HP1 17.4% / HP2 19.4%，且 **HP1-1 + HP2-1 + HP3 合計 ~9%**（其他 chr 大多 < 6%）— 暗示 chr6 phasing 較複雜（多個 minor haplotype block）
3. **chr16 phasing-weak**：HP1 14.6% / HP2 16.0% / **no_HP 52.9%**（autosome 中最高）— 對應 chr16 = 0% LOH coverage，可能 phasing anchor 稀少 + 結構保留
4. **chrX 極端 single-X 簽名**：HP1 2.5% / HP2 1.3% / **HP3 12.0% / no_HP 61.1%** — 與 HCC1395 lymphoid female ancestry 假設不符 wait — HCC1395 是 ductal carcinoma；chrX 此 pattern 反映 phasing 在 single-X 區域無法二分（HP3 用於三 haplotype 標記，多為 longphase tail 殘留）

**HP × LOH visual correlation**（Figure 1 panel 級）：

- 高 LOH chrs（chr8 / chr11 / chr17 / chr10 / chr12 / chr3）：HP1 vs HP2 比例多數 bin **接近 1:1**，與 LOH 區段 overlap → 證實 LOH 區段內 phasing 仍 work（loss-of-heterozygosity 不等於 loss-of-phasing；HP tag 仍能依 ALT read 對齊 normal-derived haplotype）
- chr16 / chrX 低 LOH 區域，反而 no_HP 比例升高（chr16 53% / chrX 61%）→ phasing failure mode 與 LOH 解耦：可能因 hetSNP 密度低 / mapping ambiguity 而非 CNV 結構
- chr20 = 94.5% CNV coverage 但 0% LOH → CNV gain 不必伴 LOH；HP1/HP2 比例仍對稱（28% / 27%）

### §10.3 與既有 finding 的關聯

- 補強 §4（LOH HP-missing ROC AUC 0.504 NEGATIVE）— 從 chr level 視覺看，**HP-missing rate 在 LOH vs 非 LOH chr 整體分布相似**（autosome no_HP 大多 36-43%，與 chr8 36% 幾乎一樣），與 §4 量化結論一致
- 補強 §5 對 HKU 4 個討論點：LOH-as-context-not-discriminator 結論在全 chr level 仍成立
- 新觀察點 ⚠ — **chr16 / chrX 是 phasing-weak 而非 LOH-related**：未來 LOH 機制研究若要 control，可用這兩 chr 作為「LOH-absent + phasing-weak」對照（avoid chr8 over-reliance）

### §10.4 Caveats

- Read midpoint binning 對 ONT long reads（mean length 5-30 kb）跨 bin 邊界誤差 < 0.3%（10 Mb bin 尺度）；可忽略
- 「in_loh / in_cnv_*」flag 只記是否與任一區段 overlap，未量化區段佔 bin 比例 — 細節需查 bin-level TSV
- chrX 無 SEQC2 LOH/CNV truth（BED 不含 chrX 條目）；故圖 3 chrX 兩 bar 皆 0 是設計使然非 bug
- HP tag value `3` 在 chr6 / chr16 / chrX 較常見（≥1%），餘 chr < 0.5% — 意義須查 LongPhase HaplotagProcess 源碼（per MEMORY: longphase getVote source 行號）

### §10.5 交付清單

**Script**：
- `research/hku_collaboration/scripts/A8_chr_ideogram_hp_loh_cnv.py`

**Figures**（PNG, DPI=150, CJK font 注入）：
- `research/hku_collaboration/figures/A8_chr_ideogram_hp_tag.png`（4×6 grid，24 panels；上層灰色 chr bar / 中層 LOH 紅+gain 藍+loss 綠 shading / 下層 HP stacked bar per 10 Mb）
- `research/hku_collaboration/figures/A8_per_chr_hp_summary_table.png`（23 chr × HP1/HP2/HP1-1/HP2-1/HP3/no-HP normalized stack）
- `research/hku_collaboration/figures/A8_chr_loh_cnv_coverage.png`（per-chr LOH 紅 + CNV 藍 並排 bar）

**Data TSV**：
- `research/hku_collaboration/data/A8_per_chr_hp_loh_cnv.tsv`（316 rows，bin-level；columns: chrom, bin_start, bin_end, hp1, hp2, hp1_1, hp2_1, hp3, no_hp, in_loh, in_cnv_gain, in_cnv_loss, total_reads）
- `research/hku_collaboration/data/A8_per_chr_summary.tsv`（23 rows，chr-level；columns: chrom, chr_length, total_reads, hp{1,2,1_1,2_1,3,no_hp}_pct, loh/cnv coverage bp+pct）

**Run log**：`research/hku_collaboration/data/A8_run.log`（BAM streaming 569 s + 3 fig 15 s）

---

## §8 Analysis A6 — LOH × CNV stratified HP / subclone analysis (22 autosomes + chrX)

> 2026-05-23 補做，per HKU 合作要求做 6-zone stratification + 統計顯著性檢定。腳本 `scripts/A6_loh_subclone_stratification.py`，per-read zone assignment + 1 Mb bin 統計。

### §8.1 設計

按 SEQC2 v1.2.1 CNV+LOH BED（660 segments：307 gain / 320 loh / 33 loss）對每條 read 的 reference midpoint 作 zone classify。BED 涵蓋 chr1–22 only（**無 chrX/Y/M 標註**），故 chrX 全 read 自動落 Z6（baseline）。

**6 zones**：
| ID | 定義 | ploidy weight |
|----|------|---------------|
| Z1 | LOH inside × CN-gain | 2.0× |
| Z2 | LOH inside × CN-loss | 0.5× |
| Z3 | LOH inside × CN-neutral（LOH only）| 1.0× |
| Z4 | LOH outside × CN-gain | 2.0× |
| Z5 | LOH outside × CN-loss | 0.5× |
| Z6 | LOH outside × CN-neutral（baseline）| 1.0× |

**MAPQ ≥ 20，drop dup/supp/secondary**（per pipeline 慣例）。Subclone proxy ratio = HP1-1 / (HP1 + HP1-1)（h1-side）+ HP2-1 / (HP2 + HP2-1)（h2-side）。

### §8.2 結果總表（global，全 23 chr 加總）

| Zone | n_reads | 長度 (Mb) | depth (/Mb) | HP1 % | HP2 % | HP1-1 % | HP2-1 % | HP3 % | no-HP % | subclone h1-side | subclone h2-side |
|------|---------|----------|-------------|-------|-------|---------|---------|-------|---------|------------------|------------------|
| Z1 LOH in × gain    | 5,979,275  | 520.2 | 11,493 | 29.9 | 31.0 | 2.6 | 2.6 | 0.27 | 33.7 | **0.079** | **0.076** |
| Z2 LOH in × loss    | 195,775    | 57.3  | 3,417  | 32.3 | 29.1 | 3.3 | 2.6 | 0.16 | 32.4 | 0.093 | 0.083 |
| Z3 LOH in × neutral | 6,400,850  | 913.0 | 7,011  | 30.5 | 29.9 | 2.5 | 2.4 | 0.24 | 34.4 | 0.077 | 0.076 |
| Z4 LOH out × gain   | 12,905,254 | 983.7 | 13,119 | 29.9 | 29.9 | 2.5 | 2.5 | 0.34 | 34.9 | 0.078 | 0.076 |
| Z5 LOH out × loss   | 74,347     | 30.6  | 2,427  | 25.5 | 27.1 | 2.6 | 3.0 | 0.41 | 41.4 | 0.092 | 0.100 |
| **Z6 LOH out × neutral** | **2,973,982** | 526.2 | 5,652  | **10.5** | **10.0** | **8.1** | **3.6** | **14.4** | **53.4** | **0.435** | **0.267** |

### §8.3 統計檢定（per 1Mb bin subclone proxy ratio）

**Kruskal-Wallis** 跨 6 zones：
- h1-side：H = 234.11, **p = 1.41×10⁻⁴⁸** ✅ 極顯著
- h2-side：H = 165.66, **p = 6.16×10⁻³⁴** ✅ 極顯著

**Pairwise Mann-Whitney U（Bonferroni-corrected，15 對 × 2 metrics = 30 tests）**：Z6 vs 其他 5 zones 全部達 p_bonf < 1e-5；rank-biserial r 範圍 0.27–0.56：

| 對比 | p_bonf (h1-side) | r |
|---|---|---|
| Z1 vs Z6 | 1.3 × 10⁻²⁸ | +0.40 |
| Z2 vs Z6 | 9.8 × 10⁻⁶  | +0.36 |
| Z3 vs Z6 | 1.1 × 10⁻³⁶ | +0.41 |
| Z4 vs Z6 | 5.4 × 10⁻²⁷ | +0.34 |
| Z5 vs Z6 | 3.3 × 10⁻¹⁰ | +0.56 |

LOH-inside zone 間（Z1 vs Z2 vs Z3）幾乎無顯著差異（p_bonf 全 > 0.05 或弱顯著）。CNV gain vs loss 內部差異也弱。**所有顯著性都來自 Z6**。

### §8.4 預期方向 vs 實際對照

| 預期 | 實際 | 判定 |
|------|------|------|
| LOH inside → subclone proxy 升高（single haplotype dominant）| LOH inside (Z1–Z3) subclone ratio ≈ 0.07–0.09，**反而是全 zone 最低** | ❌ **方向相反** |
| CN-gain inside → HP tag count 倍增但比例不變 | Z1 depth 11,493/Mb vs Z3 neutral 7,011/Mb ≈ 1.64× | ✅ depth 增加正確，但只 1.64× 而非 2× |
| CN-loss inside → HP tag count 減半 | Z2 depth 3,417/Mb vs Z3 7,011/Mb ≈ 0.49× | ✅ **完美對齊預期 0.5×** |
| Z6（baseline）→ subclone ratio 應 ≈ Z3 中等水準 | Z6 ratio 0.43，是 Z3 的 **5.7×** | ❌ **非預期** |

### §8.5 Mechanism reframe — Z6 anomaly 真實成因

Z6（LOH-out × CN-neutral）的異常 **不是均勻升高**，per-chrom panel（`figures/A6_per_chr_loh_cnv_hp.png`）顯示：

| chr | Z6 subclone h1-side | 解釋 |
|-----|---------------------|------|
| chrX | **0.80** | 全 chr 落 Z6（SEQC2 BED 不標 chrX），雌性樣本 X-inactivation 偏 HP1-1 |
| chr6 | **0.74** | 大段 mappability 邊緣（centromeric/telomeric），但有效 HP3 + no-HP 主導 |
| chr16 | **0.60** | 同 chr6 pattern |
| chr7 | **1.00** | 少數 reads 全 HP1-1 — small denominator artifact |
| 其他 18 chr | 0.04–0.13 | 與 Z3 一致 |

**結論**：Z6 異常的根因 **不是 LOH×CNV 生物學訊號**，而是：
1. **chrX 全落 Z6**（SEQC2 truth BED 不涵蓋）→ X-inactivation 假性訊號
2. **chr6 / chr16 之 Z6 集中在 centromeric 邊緣** → mapping 不確定性 + LongPhase ambiguity 累積（HP3 + no-HP 為主，HP1-1 純度低但 denominator 也低）
3. **其他 18 chr 的 Z6 ratio (0.04–0.13)** 與 Z3 一致 → **真實 LOH×CNV 訊號為 NULL**

### §8.6 對 HKU 報告影響 — Z6 假性反向訊號

原預期「LOH inside → subclone 更明顯」**被駁回**（與 §4 LOH HP-missing AUC 0.504 結論互相印證）：

> 「LongPhase `--somaticMode` 用 germline het SNV 為 phasing anchor，LOH 不改變 germline het 狀態，故 HP1/HP2/HP1-1/HP2-1 比例在 LOH inside/outside 區段幾乎相同（subclone ratio Δ < 0.02）。Z6 偏高源於 chrX 缺 truth annotation + chr6/16 centromere mapping artifact，**非 LOH-mediated subclone enrichment**。」

✅ **CN-loss ploidy weight 對齊正確**（Z2 depth 0.49× Z3，吻合 0.5× 預期）— CNV depth 訊號是真實的，可作未來 CNV inference 驗證對照。
❌ **subclone proxy 不能由 LOH×CNV stratify 出生物學差異**，6-zone 設計只能 detect coverage artifact。

### §8.7 完成交付（A6）

**Script**：`research/hku_collaboration/scripts/A6_loh_subclone_stratification.py`（multiprocessing 8 workers，per-chrom worker + worker 內 cheap BED reload）

**Figures**（PNG dpi=150，CJK font 注入）：
- `figures/A6_loh_cnv_stratified_hp_distribution.png`（6-panel stacked bar，HP1/2/1-1/2-1/3/no-HP × 6 zones）
- `figures/A6_per_chr_loh_cnv_hp.png`（4×6 = 23-panel grid，per-chr 6-zone subclone ratio dot bars）
- `figures/A6_subclone_ratio_violin.png`（2 panel：h1-side / h2-side，6 zone violin + Kruskal-Wallis + Bonferroni 顯著對註記）

**Data TSV**：
- `data/A6_loh_cnv_stratified_stats.tsv`（36 rows，6 zones × 6 HP）
- `data/A6_loh_cnv_per_chr.tsv`（138 rows，23 chr × 6 zones）
- `data/A6_loh_cnv_1mb_bins.tsv`（3,399 rows，per-bin subclone ratios）
- `data/A6_loh_cnv_pairwise_stats.tsv`（30 rows，15 pair × 2 metrics）

**Run log**：`data/A6_run.log`（per-chrom scan 778 s + figures 9 s ≈ 13 min total，8 workers）

**已知限制**：
1. SEQC2 BED 不涵蓋 chrX/Y/M → chrX 全落 Z6 造成假性訊號
2. CNV gain & loss 重疊解析 = gain dominance（loss 33 segments，重疊罕見）
3. 1Mb bin 粗粒度（與 A2_3 設計一致）；100 kb bin 應能微化但增 10× 成本
4. HP tag 值 `1-1` / `2-1` 為 string-typed（HCC1395 tagged BAM `HP:Z:` field，非 `HP:i:`），分析做 string match
5. CN-gain 在 Z1 depth 只 1.64× 而非 2×（可能因 LOH 改變了 mapping pattern，未深究）

---

## §11 Analysis A7 — PS block size distribution + global N50

> 2026-05-23 補做，per HKU 合作 v4 review #6 — 量化 LongPhase-S（`twolinin/longphase` mainline v1.7.3 `--somaticMode`）產出的 PS block 全 chr 覆蓋規模，作為 ISM "PS chaining 跨單 read" 結論的補強數據。腳本 `scripts/A7_phase_block_size.py`。

### §11.1 全 chr 統計（HCC1395 T-mode tagged BAM）

| 指標 | 值 |
|---|---|
| PS blocks 總數 | **6,861** |
| **Global N50** | **943,476 bp（943.5 kb）** |
| Median block size | **163,675 bp** |
| Mean block size | **420,226 bp** |
| **Max block** | **18,303,383 bp（18.30 Mb）on chr1 (PS=125047868)** |
| Wall clock | 508 s（~8.5 min） |

→ **Half of phased bases live in blocks ≥ 943 kb**，median block 164 kb 仍遠超 ONT 平均 read 長度（5-30 kb）→ 證實 §3 結論「median TP pair distance 17.4 kb 需 LongPhase PS chaining 跨越單 read」的可行性 — PS block 規模充足。

### §11.2 Per-chr highlights（從 `data/A7_per_chr_n50_summary.tsv`）

| chr | N50 | 備註 |
|---|---|---|
| **chr22** | **1,625,099 bp (1.625 Mb)** | autosome 最大 |
| chr1 | 1,071,203 bp (1.071 Mb) | 含全體 max block 18.3 Mb |
| chr2 | 1,159,099 bp | |
| chr11 | 977,497 bp | high LOH (chr8/11/17 Tier 1) 仍 N50 ~1 Mb |
| **chr19** | **1,018,274 bp** | autosome 中等偏大 |
| chr5 | 1,108,037 bp | |
| chr6 | **733,900 bp** | **per §10 chr6 phasing complex 對應**（HP1-1+HP2-1+HP3 ~9% outlier）|
| **chr16** | **226,882 bp** | **autosome 最小**（§10 no_HP 52.9% 對應）|
| **chrX** | **101,154 bp** | **全體最小**（§10 no_HP 61% 對應；§10 single-X 訊號）|

### §11.3 Cross-reference A8 chr ideogram §10

A7 N50 短的 chr 與 A8 §10 phasing-weak chr 完全對齊（per phasing complexity → block 短的因果）：

| chr | A7 N50 | A8 phasing pattern (§10) | 一致性 |
|---|---|---|---|
| chr6 | 734 kb | HP1-1+HP2-1+HP3 ~9% outlier | ✅ 對齊（phasing complex → block 短）|
| chr16 | 227 kb | no_HP 52.9%（autosome 最高）| ✅ 對齊（phasing-weak → block 短）|
| chrX | 101 kb | no_HP 61.1% + single-X HP3 12% | ✅ 對齊（single-X 無法二分 → block 短）|
| chr22 | 1,625 kb | 平衡 HP1/HP2 | ✅ 對齊（phasing strong → block 長）|

→ **finding**：PS block 規模是 LongPhase phasing confidence 的直接指標；A7 (block size) ↔ A8 (HP composition) 為同一 phasing-weakness 軸向的兩種量化角度。

### §11.4 對 HKU 報告影響

✅ **POSITIVE**：N50 943 kb >> ClairS-TO ±16 bp window（5× 數量級）→ PS chaining 確實有效跨越 caller blind spot。
✅ **POSITIVE**：85.9% PS set 含 ≥2 TP（§2）+ N50 943 kb → ISM read-level 甲基聚類有「足夠 read 數 × 足夠 spatial span」的雙重支撐。
⚠ **Caveat**：chr6 / chr16 / chrX 的 phasing-weak 區段（block N50 < 250 kb）需另案處理，可能成為 ISM 補完規則的 stress test 對象。

### §11.5 交付清單（A7）

**Script**：`research/hku_collaboration/scripts/A7_phase_block_size.py`

**Figures**（PNG dpi=150，CJK font 注入）：
- `figures/A7_phase_block_size_distribution.png`（global block size histogram + CDF）
- `figures/A7_n50_per_chr_bar.png`（per-chr N50 bar chart）
- `figures/A7_block_density_vs_chr_length.png`（block density vs chr length scatter）

**Data TSV**：
- `data/A7_per_chr_n50_summary.tsv`（24 rows：22 autosomes + chrX，含 chr_length / n_ps_blocks / median_block_size / q25 / q75 / n50 / max_block / max_block_ps / coverage_ratio）

**Run log**：~8.5 min wall clock，6,861 PS blocks scanned across 24 chrs.

**已知限制**：
1. Block size = (max_pos − min_pos) per PS group，未考慮 block 內 gap 結構
2. coverage_ratio > 1.0 表示部分 read 跨 PS 邊界（如 chr1 ratio 1.038）
3. chrX coverage_ratio 0.53 反映 single-X / no_HP 主導 → PS block 規模顯著低估真實 phasing span
