<!--
建立時間: 2026-05-24
類型: HKU handoff T10 — LongPhase-S Reviewer 3 Major #5 response (24-chr full extension)
主 source: pysam scan on HCC1395 Tmode tagged BAM + SEQC2 v1.2.1 truth (PASS)
範圍: chr1..chr22 + chrX 全 24 條染色體 HP3 read fraction + HP3 TP rate
延伸: T1 (chr1/chr8/chr19 pilot) → 全 24 chr
-->

# T10 — HP3 read fraction & TP rate per chromosome (HCC1395, 全 24 chr)

## Summary one-liner

**HP3 read 全基因組 28.53 M reads 中佔 503,095 (1.76 %)，落於 SEQC2 PASS truth 位點 71,952 reads (pooled TP rate 14.30 %)；然而扣除三個結構性 outlier (chr6 / chr16 / chrX) 後，20/23 常染色體呈現極為一致的 baseline：HP3 fraction 0.217–0.380 % (mean 0.283 %)、HP3 TP rate 78.6–95.4 % (mean 90.4 %)。outlier 三條 chr 揭露 LongPhase-S `--somaticMode` HP3 bucket 在 HLA region / repeat-rich segdup / 性染色體單套量區域會收斂為 phasing-failure fallback bucket（HP3 fraction ↑↑、TP rate ↓↓），與 reviewer 認為「HP3 = somatic-evidence-only」的單一解釋不符。**

> Reviewer 3 Major #5（LongPhase-S paper）詢問「HP3 reads 在總 reads 中佔比 + HP3 中 somatic TP 比例」。T1 已答覆 chr1/chr8/chr19 三條 chr。本檔 T10 將 scope 擴展至全 24 chr (chr1-22 autosomes + chrX；chrY 因 HCC1395 為女性樣本 + SEQC2 truth 無 chrY 條目而跳過)。

---

## §1 跑法

- **Tagged BAM**: `/big7_disk/liaoyoyo2001/data/bam/HCC1395_Tmode_tagged_ClairS_pileup_v040_woFilter.bam` (283 GB)
- **Truth VCF**: `data/answer/SEQC2/SEQC2_high-confidence_sSNV_in_HC_regions_v1.2.1.vcf.gz`（PASS-only；FILTER 包含 `PASS` 即收，bug fix vs naive `filt == ['PASS']`）
- **Read filter**: `MAPQ ≥ 20`、排除 secondary / supplementary / duplicate / unmapped（primary only）
- **HP3 判定**: `read.get_tag("HP")` 字串等於 `"3"`（LongPhase `--somaticMode` 寫 `HP:Z:3`）
- **HP3 TP 判定**: 該 HP3 read `get_reference_positions()` (0-based) +1 與該 chr 之 SEQC2 truth 位點 set 有交集 ≥1
- **Script**: `research/hku_collaboration/scripts/T10_hp3_tp_rate_24chr.py`（+ resume helper `T10_make_plot.py`）
- **Wall clock**: truth load 1.0 s + 全 24 chr scan ≈ 88 min（chr1 444s 最長，chr21 67s 最短）

**Sanity check（truth load）**: chr1 = 3,440 / chr8 = 2,605 / chr19 = 874 — 全數對齊 T1 baseline。chrX = 0 PASS（SEQC2 v1.2.1 truth set 不收 chrX 條目）。
**Sanity check（BAM count）**: chr1 / chr8 / chr19 三條 total_reads / hp3_reads / hp3_tp_reads 全與 T1 完全一致（`[T1 baseline OK]` × 3 在 stderr log 驗證）。

---

## §2 結果完整表

| chr | total_reads | hp3_reads | hp3_fraction | hp3_tp_reads | hp3_tp_rate | truth_PASS |
|---|---:|---:|---:|---:|---:|---:|
| chr1  | 2,236,562 | 6,419   | 0.287 %    | 5,658 | 88.14 % | 3,440 |
| chr2  | 2,241,428 | 7,482   | 0.334 %    | 6,845 | 91.49 % | 3,608 |
| chr3  | 2,030,234 | 4,642   | 0.229 %    | 4,153 | 89.47 % | 2,743 |
| chr4  | 1,693,079 | 4,845   | 0.286 %    | 4,446 | 91.76 % | 3,006 |
| chr5  | 1,728,890 | 4,901   | 0.283 %    | 4,539 | 92.61 % | 2,630 |
| **chr6**  | 1,649,102 | **155,184** | **9.410 %** ⚠ | 2,131 | **1.37 %** ⚠ | 1,599 |
| chr7  | 2,379,235 | 8,642   | 0.363 %    | 8,144 | 94.24 % | 3,057 |
| chr8  | 1,653,223 | 5,539   | 0.335 %    | 4,454 | 80.41 % | 2,605 |
| chr9  | 1,148,711 | 2,495   | 0.217 %    | 2,256 | 90.42 % | 1,656 |
| chr10 | 1,185,325 | 3,769   | 0.318 %    | 3,511 | 93.15 % | 1,712 |
| chr11 | 1,078,826 | 3,035   | 0.281 %    | 2,894 | 95.35 % | 1,594 |
| chr12 | 1,186,354 | 3,202   | 0.270 %    | 2,942 | 91.88 % | 1,747 |
| chr13 |   712,113 | 1,940   | 0.272 %    | 1,726 | 88.97 % | 1,144 |
| chr14 | 1,136,753 | 4,325   | 0.380 %    | 4,104 | 94.89 % | 1,574 |
| chr15 |   883,279 | 2,144   | 0.243 %    | 1,922 | 89.65 % | 1,250 |
| **chr16** | 1,091,872 | **93,638**  | **8.576 %** ⚠ | 3,306 | **3.53 %** ⚠ | 808   |
| chr17 |   772,342 | 2,209   | 0.286 %    | 1,870 | 84.65 % | 1,079 |
| chr18 |   777,662 | 2,048   | 0.263 %    | 1,902 | 92.87 % | 1,116 |
| chr19 |   535,588 | 1,264   | 0.236 %    | 1,184 | 93.67 % |   874 |
| chr20 |   765,439 | 2,485   | 0.325 %    | 2,325 | 93.56 % | 1,106 |
| chr21 |   282,365 |   616   | 0.218 %    |   484 | 78.57 % |   472 |
| chr22 |   522,931 | 1,259   | 0.241 %    | 1,156 | 91.82 % |   627 |
| **chrX**  |   838,170 | **181,012** | **21.596 %** ⚠ | **0** | **0.00 %** ⚠ |     0 |
| **TOTAL** | **28,529,483** | **503,095** | **1.763 %** (pooled) | **71,952** | **14.30 %** (pooled) | **39,447** |

完整精度 TSV：`findings_5_24/T10_HP3_TP_rate_24chr.tsv`
24-chr 雙軸條形圖：`figures/T10_HP3_TP_rate_24chr.png`（紅色背景 = HP3 fraction outlier >5%；藍色虛線框 = chr8 LOH 99% hotspot）

---

## §3 跨 chr 跨度與分布

### 3.1 Normal-chr baseline (20 條，排除 chr6 / chr16 / chrX)

| 統計 | HP3 fraction | HP3 TP rate |
|---|---:|---:|
| min  | 0.217 % (chr9)  | 78.57 % (chr21) |
| max  | 0.380 % (chr14) | 95.35 % (chr11) |
| mean | 0.283 %         | 90.38 %         |
| range / spread | 1.75×       | ±8.4 pp around mean |

→ **HP3 fraction 在 20 條正常 chr 跨度僅 1.75×**（vs. 三個 outlier 全部 > 22× normal mean），證實 `--somaticMode` HP3 分配在「典型雙套量、非 HLA 區」的 chr-invariance 假說在 T1 三條 pilot 之外仍站得住。
→ **HP3 TP rate 在 20 條正常 chr range 78.6–95.4 % (mean 90.4 %)**，跨度小且全部位於高 enrichment 區間（隨機背景 < 0.1 % per §3.3）。

### 3.2 Outlier 三條 chr（chr6 / chr16 / chrX）

| chr | HP3 frac | TP rate | 推測機制 |
|---|---:|---:|---|
| **chr6**  | **9.41 %** (33× normal) | **1.37 %** (66× ↓) | HLA / MHC region (chr6:28–34 Mb) 極高 polymorphism + segdup 密度 → LongPhase phasing 大量 failure → HP3 collapse 為 unphased fallback bucket；TP rate 崩潰因為這些 read 多數不覆蓋 somatic 位點 |
| **chr16** | **8.58 %** (30× normal) | **3.53 %** (26× ↓) | chr16 pericentromeric / acrocentric short-arm + 大量 segdup（HCC1395 reported chr16q amplicon），phasing 失效同 chr6；TP rate 略高於 chr6 但仍 1 個量級差距 |
| **chrX**  | **21.60 %** (76× normal) | **0.00 %** | HCC1395 為女性樣本（XX）但 chrX 易遭 LOH + skewed inactivation；SEQC2 v1.2.1 truth set **無 chrX PASS 條目** (truth_count=0) → TP rate 必為 0；HP3 fraction 大幅膨脹反映 LongPhase 在性染色體 phasing block 容易 fragment（即使是 XX）|

### 3.3 隨機背景對照

- 全 24 chr SEQC2 truth = 39,447 PASS / 約 3 × 10⁹ bp ≈ 1.3 × 10⁻⁵ /bp
- 一條 ~10 kb read 撞到任一 truth 位點之隨機機率 ≈ 0.013 % → normal-chr HP3 TP rate ~90 % 是隨機機率的 ~6,900× 富集
- chr6 / chr16 outlier 即使在最低 1.37 % 仍是隨機背景的 ~100×，**outlier 區的 HP3 並非「完全隨機 fallback」**，仍保留部分 somatic-evidence preference（只是被 phasing-failure inflation 稀釋）

---

## §4 對 R3 M5 reviewer 的回應（3 sentences）

> **In HCC1395 Tmode tagged BAM (LongPhase v1.7.3 `--somaticMode` over ClairS-TO ssrs v0.4.x), HP:Z:3 reads occupy 1.76 % of 28.5 M post-filter reads genome-wide, and 14.3 % of those HP3 reads cover at least one SEQC2 PASS somatic position (pooled HP3 TP rate). However, 20 of 23 chromosomes (excluding chr6, chr16, chrX) show a tight baseline of HP3 fraction 0.22–0.38 % (mean 0.28 %) and HP3 TP rate 78.6–95.4 % (mean 90.4 %), meaning HP3 strongly enriches for somatic-evidence reads (~6,900× over random background) outside three structural outliers. The outliers — chr6 (HP3 frac 9.41 %, TP 1.37 %), chr16 (8.58 %, 3.53 %), chrX (21.60 %, 0.00 % — SEQC2 v1.2.1 has no chrX PASS positions) — reveal that HP3 also functions as a fallback bucket when LongPhase phasing fails in HLA / MHC, segdup-rich, or single-copy sex chromosome regions, so HP3 should be described as somatic-evidence-enriched in well-phased regions rather than universally somatic-only.**

---

## §5 與 T1 (3-chr pilot) 對照

| 來源 | scope | HP3 reads | total reads | HP3 frac | HP3 TP rate |
|---|---|---:|---:|---:|---:|
| T1 (chr1+chr8+chr19) | 3 chr | 13,222   | 4.43 M | 0.299 % | 85.43 % (pooled) |
| **T10 (全 24 chr)**  | 24 chr | **503,095** | **28.53 M** | **1.763 %** | **14.30 %** (pooled) |
| T10 normal-only (20 chr) | 20 chr | 73,261 | 25.0 M | 0.293 % | ~90 % (mean) |

→ **全 24 chr pooled HP3 fraction 1.76 % 比 T1 3-chr 0.30 % 高 ~6×**，差異完全由 chr6 + chr16 + chrX 三個 outlier 拉高（contribute 429,834 HP3 reads = 全 HP3 的 85.4 %）。
→ **全 24 chr pooled TP rate 14.3 % 比 T1 85.4 % 低 ~6×** 同理由 outlier 主導（chrX 181 k HP3 全部 TP=0、chr6 + chr16 也都 < 4 %）。
→ **若以「normal-chr only」報告**：HP3 fraction 0.293 % / TP rate ≈ 90 %，與 T1 pooled 0.299 % / 85.4 % 完全一致。建議 reviewer response 同時呈現 raw pooled + outlier-excluded 兩個數字。

---

## §6 Per-chr 觀察重點

1. **chr8 LOH 99 % hotspot（藍色虛線框）但 HP3 fraction 仍正常 (0.34 %)**：與 chr6/chr16 outlier 區隔關鍵。LOH 高 ≠ phasing failure；LOH 區之 phasing block 仍有效，只是 het 變 hom。chr8 TP rate 80.4 % 是 normal range 下緣（vs chr11 95.4 %），可能反映 LOH 區 ClairS-TO 有較多 FP candidate（前述 T1 §4 已論述）。
2. **chr14 HP3 fraction 0.380 %（normal-chr max）**：略高於 mean，但 TP rate 94.9 % 也最高之一，無 outlier 含意。
3. **chr21 TP rate 78.6 %（normal-chr min）**：絕對數小（HP3 reads 僅 616），統計噪音較大；TP reads 484 仍合理。
4. **chr17 TP rate 84.7 %（中等偏低）**：HCC1395 chr17 TP53 區有 LOH，類似 chr8 邏輯但程度較輕。
5. **chr6 outlier 與 HLA position 的對應**：chr6:28–34 Mb 是 HLA class I/II 區，~6 Mb 佔 chr6 (170 Mb) 3.5 %；若 outlier 主導效應局限該 6 Mb，本地 HP3 fraction 估計可達 ~80 %（待 fine-grained per-window 驗證，T10 範圍外）。
6. **chr16 outlier 與 16p11/16q12 segdup 對應**：chr16 約 25 % 為 segdup（vs 全 genome 5 %），phasing failure 普遍。
7. **chrX 21.6 % outlier**：HCC1395 為 ATCC-derived 女性 lung-derived cell line。但 chrX 仍可能有 partial LOH 或 skewed XCI；外加 SEQC2 truth 完全無 chrX 條目，TP rate=0 為定義性結果，非品質問題。

---

## §7 Caveats

- 「TP read」定義 = read alignment 覆蓋任一 truth pos，**不要求該 read 實際支持 alt allele**；嚴格意義是「該 HP3 read 在 truth-rich locus 上 align」。allele-level concordance 需另跑 pysam pileup TSV（未在 T10 範圍）。
- chr6 / chr16 outlier 機制屬「推測（HLA / segdup）」；要嚴格驗證需 stratify by repeat-masker / segdup-track 重跑 HP3 fraction per 1 Mb window，T10 範圍外。
- chrX TP rate 0 是 truth set 限制（SEQC2 v1.2.1 chrX = 0 PASS），非 LongPhase 問題；reviewer 若追問 chrX 應澄清「truth coverage gap」而非「HP3 失效」。
- 全 24 chr scan 與 T1 結果完整對齊（chr1/chr8/chr19 三條 total/hp3/hp3_tp 數字逐筆一致）。
- BAM 為 `ClairS-TO ssrs v0.4.x + LongPhase v1.7.3 --somaticMode` 產出；換 caller 版本會變。
- MAPQ ≥ 20 / drop dup-sup-sec filter 與 T1 / deliverable A 一致。

---

## §8 交付清單

- **Script (main)**: `research/hku_collaboration/scripts/T10_hp3_tp_rate_24chr.py`
- **Script (plot helper)**: `research/hku_collaboration/scripts/T10_make_plot.py`
- **TSV**: `research/hku_collaboration/findings_5_24/T10_HP3_TP_rate_24chr.tsv`（24 rows + header；含 seqc2_truth_count 欄）
- **PNG**: `research/hku_collaboration/figures/T10_HP3_TP_rate_24chr.png`（雙軸條形圖，紅色 outlier highlight、藍色 chr8 LOH 框）
- **本檔**: `research/hku_collaboration/findings_5_24/T10_HP3_24chr_findings.md`

---

## §9 Provenance

- **BAM**: `/big7_disk/liaoyoyo2001/data/bam/HCC1395_Tmode_tagged_ClairS_pileup_v040_woFilter.bam` (283 GB, `ClairS-TO ssrs v0.4.x + LongPhase v1.7.3 --somaticMode`)
- **Truth**: `data/answer/SEQC2/SEQC2_high-confidence_sSNV_in_HC_regions_v1.2.1.vcf.gz` PASS-only（FILTER 含 PASS 即收，多 label 兼容）
- **Sample**: HCC1395（lung cell line，female，SEQC2 reference cancer sample）
- **Run date**: 2026-05-24
- **Wall clock**: 88 min total（chr1-8 第一次 run 41 min, killed; chr9-X resume run 44.7 min）
- **pysam**: 0.23.3 / matplotlib: Agg backend
- **Font**: Droid Sans Fallback (`/usr/share/fonts/truetype/droid/DroidSansFallbackFull.ttf`)
- **Baseline cross-ref**: T1 (`findings_5_24/T1_HP3_findings.md`) chr1/chr8/chr19 三條完全對齊；T10 是 T1 的 superset extension
