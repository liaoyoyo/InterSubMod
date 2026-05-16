# Thread D — TP-Enriched Phasing Signatures (LOH × cross_het) Tool Paper Draft

> **狀態**：initiated
> **建立日期**：2026-05-16
> **專案目錄**：`research/thread_d_paper/`
> **Parent plan**：`~/.claude/plans/tender-pondering-blossom.md`
> **主軸來源**：`InterSubMod/docs/reports/validated/2026/04/20260426_Thread_D_LOH_constrained_phasing_main_axis_01.md`

## 背景與動機

InterSubMod 14 個月研究經歷 4 次方向 pivot 後（filter → characterization → LOH-constrained → TP-enriched dual-pillar），於 2026-05-15 透過 V3F/V5/V6 三向 ISM 比較 + 4 樣本擴展 + paradigm reframe 確立 **「TP-enriched phasing signatures (LOH × cross_het) 雙重機制」** 為論文主軸。

**為什麼是 tool paper（非 specific discovery）**：
- 14 個月研究全部基於 cancer cell line，無 patient cohort（無法支撐 Genome Research 級 specific discovery 對 cross-population validation 的要求）
- InterSubMod 本身為首個整合 methylation × phasing × CN 的 read-level framework，contribution 在工具設計與機制可解釋性
- Tool paper 適合 Bioinformatics / NAR Genomics and Bioinformatics 級

**核心 contribution 一句話**：「InterSubMod is the first read-level framework integrating methylation × phasing × CN to characterize TP-enriched phasing signatures in tumor-only ONT long-read somatic calling.」

## 假說

### H1：Inner LOH NG=2 same-hap TP-enriched signature

**陳述**：LOH 區內 NG=2 同 haplotype bucket (`{HP1, HP1-1}` or `{HP2, HP2-1}`) 占用模式呈現顯著高 TP rate，相對 Outer cross_het bucket 形成 +0.302 median gap（n=6, Wilcoxon p=0.0156, X5 KDE-corrected）

**前提條件**：
- KDE-corrected master available (commit 374fad4)
- X3 `--germline-hp-only` flag-on negative control 確認 same-hap bucket 物理消失
- B3 paired-mode negative control 確認 gap 在有 matched normal 時消失

**已知 Confound**：
- HCC1954 caller FP background ceiling（Outer cross_het TP=0.084 vs 其他 6 樣本 0.55-0.88）
- HPFineN_HP1S self-reference（R-SELFREF）— bucket 定義本身包含 somatic attribution
- HCC1937 BRCA1 driver outlier chr15/17/14

**驗證標準**：
- **Positive**：n>=6 樣本 Wilcoxon p<0.05 + paired NC gap |median|<0.05 + flag-on bucket collapse
- **Negative**：n>=7 expansion (Archive TO rerun + COLO829) 下 <5/7 同方向

### H2：Outer cross_het TP-pure signature（2026-05-15 paradigm reframe）

**陳述**：原假設 "Outer cross_het = germline het FP" 不成立。實測 Z-OCH (Outer cross_het) FP rate 0.017 << global 0.137（Fisher p=3.8e-62 for TP enrichment），cross_het 模式實為 somatic-evidence marker。

**前提條件**：
- V3F/V5/V6 三向 ISM 比較 (HCC1395 全 chr) 已完成（2026-05-15）
- V5 over-promote 量化：Inner LOH NG=2 V5=8,136 vs V3F=5,064 (+60%)，V6=5,353 修補回 V3F 水準
- V6 binary 是 paradigm reframe 的乾淨基準

**已知 Confound**：
- V5 Layer 1.5 over-promote artifact（V5/V3F top cell ratio 5.95× 集中 cross_het_inv bucket）
- sample-specific CNV driver（HCC1937 BRCA1）
- COLO829 ONT R10 無 methylation → Outer cross_het 計算可能受影響

**驗證標準**：
- **Positive**：Cross-sample n>=4 同方向 + Fisher p<0.001 for TP enrichment + V6 重跑後仍保持
- **Negative**：V6 archive rerun 顯示 Z-OCH FP rate 回到 global level（~0.137）

### H3：Framework coverage 37% with interpretable mechanism

**陳述**：InterSubMod 三軸 framework（methylation × phasing × CN）解釋約 37% TP/FP 訊號，剩 63% 不被現有 framework cover（需 mappability/repeat/GC/SV 新軸）。此 limitation 為 paper Discussion 主動承認的 framework boundary，不影響核心 contribution。

**前提條件**：
- V6 production tag finalized（dependency on `selfphasing_v6_production`）
- Z-AUTO KDE 跨樣本 n>=4 擴展完成（plan Tier 2.1）

**已知 Confound**：
- Framework boundary 定義本身（axis selection bias）
- Sample-specific noise inflation

**驗證標準**：
- **Positive**：37% framework coverage 在 n>=4 樣本各自可重現
- **Negative**：37% 覆蓋率 HCC1395-only，其他樣本 framework coverage <20% 或 >60%（極不一致）

### H4：Strategy A — HCC1395 primary + 6 replication

**陳述**：論文敘事採 sample hierarchy Strategy A — HCC1395（SEQC2 truth set 詳細，2026-05-15 paradigm reframe primary discovery）作 §3 main result deep dive；其他 6 樣本作 §4 replication cohort 一致性驗證。HCC1954 / HCC1937 outlier 為 §4 subsection 解釋。

**前提條件**：
- HCC1395 全 chr deep panel 完成（2026-05-15 ✅）
- 6 cross-sample 經 V6 binary KDE-corrected Archive TO rerun（plan Tier 2 待執行）

**已知 Confound**：
- SEQC2 truth set vs 其他樣本 internal truth set 詳細度差異
- Platform diversity（5kHz vs DORADO vs PAO）可能模糊「樣本一致」概念

**驗證標準**：
- **Positive**：HCC1395 deep mechanism + 6 samples 對 H1/H2 direction match
- **Negative**：機制只在 HCC1395 成立，其他樣本 direction 散亂

## 方法

### 數據來源

| 數據集 | 路徑 | 描述 | 使用欄位 |
|--------|------|------|---------|
| HCC1395 V6 ISM | `output/canonical/HCC1395/V6/` | 全 chr V6 binary ISM | HPFineN_HP1/HP1S/HP2/HP2S, Inner_LOH, CN |
| V3F/V5/V6 master | `research/v6_bam_tpfp_hp_loh_cn/step1_master_three_way.tsv` | 三向 ISM 比較 (HCC1395) | 64 cols |
| 4 樣本 V6 ISM | `output/canonical/{H1437,H2009,HCC1954,HCC1937}/V6/` | 2026-05-15 擴展 | 同 HCC1395 |
| KDE-corrected master | (commit 374fad4) | 原 Thread D X5 evidence base | — |
| Knowledge sample profile | `/big8_disk/.../Knowledge/02_samples/cancer-samples.md` | 7 樣本 LOH/CN/purity | — |

### 分析步驟

```
Step 1: V6 binary archive TO rerun (剩餘 COLO829)
        → 驗證: ls output/canonical/COLO829/V6/ 存在; step1_master_three_way 加 COLO829 row
Step 2: HCC1395 primary discovery 章節骨架 (Strategy A §3)
        → 驗證: research/thread_d_paper/01_HCC1395_primary_chapter.md 含 mechanism figure list
Step 3: 6-sample replication cohort 章節骨架 (Strategy A §4)
        → 驗證: research/thread_d_paper/02_replication_cohort_chapter.md 含 HCC1954/HCC1937 subsections
Step 4: Z-AUTO KDE 跨樣本 (n>=4)，升 ⭐3 → ⭐4
        → 驗證: research/v6_bam_tpfp_hp_loh_cn/step5_zauto_kde_cross_sample.md 含 4 sample verdict
Step 5: Paper full outline (abstract + 6 章 + 6 主圖)
        → 驗證: research/thread_d_paper/00_outline.md 含 abstract 250 words
```

### 統計方法

- **Wilcoxon signed-rank exact test** for cross-sample direction consistency (H1, H2)
- **Fisher exact test** for TP enrichment per zone (H2 Z-OCH)
- **Bootstrap 95% CI** (n_iter=10000) for median gap (B1 already done)
- **Multiple testing correction**: BH-FDR for per-sample claims; main result single-test no correction
- **Effect size**: Cohen's d / r-rank-biserial for paired/non-paired comparisons

## 可行性評估

| 因素 | 評估 |
|------|------|
| 數據可用性 | △ 5/7 V6 已跑 + COLO829 待補；KDE-corrected master 完整 |
| 計算資源 | ✓ V6 archive rerun ~10 hr parallel；後續分析 <2 hr |
| 與已有結論衝突 | ✓ 不衝突 — Thread B 撤回 + Thread D 升級 + paradigm reframe 已 align |

## 已知風險

1. **R-SELFREF (HPFineN_HP1S self-reference)**：bucket 定義包含 somatic attribution，TP rate 可能循環論證；緩解：Phase 2B master + flag=on 重驗（已 listed P0-pending）
2. **HCC1954 caller ceiling**：Outer cross_het TP 0.084 為 caller FP 背景而非 mechanism artifact；緩解：B2 root cause 已建檔 + standalone case panel
3. **HCC1937 BRCA1 driver outlier**：chr15/17/14 sample-specific；緩解：§4 subsection 解釋 driver biology
4. **63% framework gap**：reviewer 可能要求補 GC/mappability/repeat 新軸；緩解：Tier 4.2 reactive pilot 預備
5. **Cancer-only limitation**：reviewer round 1 必問 non-cancer use case；緩解：Tier 4.4 HG002 pilot reactive；接受 limitation in Discussion

## Paper Structure（preliminary）

| § | 章節 | 主軸 | Figure |
|---|------|------|--------|
| §1 | Introduction | InterSubMod framework + tumor-only ONT challenge | F1: framework overview |
| §2 | Mechanism | 4-bucket HPFine + LOH × cross_het 雙重 signature | F2: 4-bucket physical model |
| §3 | HCC1395 primary discovery | SEQC2 truth + V3F/V5/V6 三向 + paradigm reframe | F3: HCC1395 TP/FP heatmap; F4: V5 over-promote |
| §4 | 6-sample replication cohort | Strategy A + HCC1954/HCC1937 outlier | F5: cross-sample direction |
| §5 | Reproducible release | Docker + benchmark suite + tutorial | F6: pipeline DAG |
| §6 | Discussion | 63% gap + cancer-only + Phase 2A future | — |

## 相關檔案

- `InterSubMod/docs/reports/validated/2026/04/20260426_Thread_D_LOH_constrained_phasing_main_axis_01.md` — 主軸報告（已 2026-05-15 加 banner + §2.5 paradigm reframe + 改名）
- `InterSubMod/docs/experiments/in_progress/2026/05/20260515_V6_TPFP_HP_LOH_CN_Characterization_01.md` — paradigm reframe primary evidence
- `InterSubMod/docs/reports/validated/2026/04/20260426_HCC1954_Outlier_CallerCeiling_CasePanel_01.md` — HCC1954 caveat
- `InterSubMod/docs/reports/validated/2026/04/20260426_Thread_B_Retraction_Whitelist_Cross_Sample_01.md` — Thread B 撤回對照
- `~/.claude/plans/tender-pondering-blossom.md` — parent plan (Tier 1-4 timeline)

## Plan tender-pondering-blossom Tier 對照

- **Tier 1**: T1.1 主軸正名 ✅ (2026-05-16) | T1.2 V6 production tag (dependency on `selfphasing_v6_production`) | T1.3 本 scaffolding ✅
- **Tier 2**: T2.1 Z-AUTO KDE cross-sample → step5 doc | T2.2 HCC1395 chapter skeleton | T2.3 replication cohort chapter skeleton
- **Tier 3**: T3.1 paper full outline | T3.2 GitHub public | T3.3 Docker | T3.4 benchmark suite | T3.5 Discussion
- **Tier 4** (reactive): T4.1 Phase 2A normal cross-sample | T4.2 GC/mappability axis | T4.4 HG002 pilot
