---
title: ClairS Caller-Augmentation 研究方向整理與計劃
date: 2026-05-24
author: InterSubMod Research (Claude Code conversation consolidation)
status: concept
tier: ⭐3 design-stage (PROBE verdict pre-cycle)
verified_scope: ClairS main HEAD 88f887a + longphase-s main, source code 直查
parent_project: research/methyl_augmented_filter_phase2/
related:
  - InterSubMod/research/methyl_augmented_filter_phase2/cycle1/
  - InterSubMod/research/methyl_augmented_filter_phase2/cycle4/loso_validation/
  - InterSubMod/docs/concepts/2026/05/20260523_HKU_collab_主軸與細節_01.md
sources_committed: 88f887a (ClairS), longphase-s main HEAD (2026-05-24 查驗)
---

# ClairS Caller-Augmentation 研究方向整理與計劃

> **Framework**：Verdict-Pyramid（結論在前）+ 每節 SCQA
> **目的**：把 8 round 對話中對 ClairS / longphase-s pipeline 的 verification、概念釐清、研究方向討論、PROBE verdict 全部固化為 canonical 文件，後續所有相關研究 cycle 引用此文件而非重複 verify。

---

## §0 Executive Summary（Verdict 在前）

### Top-line verdict

| 項目 | 結論 |
|---|---|
| **ClairS 預設流程** ✅ verified | NVTB phasing source + parallel pileup/full heads + germline-only HP integration |
| **longphase-s 後處理 ClairS 結果** ✅ verified | Downstream filter，ClairS 不知道 longphase-s 存在；HP:z 字串編碼與 ClairS HP:i integer 不相容 |
| **「直接用 somatic-aware HP 餵 ClairS pileup/full」這個 idea** ❌ NO-GO | 三重阻礙：pretrained model schema mismatch / pileup 無 HP channel / 程式碼 substring check `if hap in '12'` 把 1-1/2-1/3 silently drop |
| **「Caller-augmentation 研究方向」整體** 🟡 **PROBE** | 方向合理（gap 真實 + ISM 有資訊優勢），但需限定 framing 為 **post-hoc re-ranker**，不重訓 caller |
| **Phase 2 Cycle 1 既有 baseline** ⭐3 strong | HCC1395 ΔF1=+0.02236（9.24× v1.0），已證明 post-filter 路線可行；本計劃 = Cycle 2 自然延伸 |

### Immediate next action

啟動 **Phase 2 Cycle 2** with pre-registration：
- **H_C2A**：加 longphase-s HP:z annotation feature → HCC1395 ΔF1 > +0.002 over Cycle 1 baseline，5/5 LOSO 4+ POSITIVE
- **H_C2B**：加 ISM methylation feature → ΔF1 > +0.001 over Cycle 1 + HP-only baseline，3/5 POSITIVE
- **Stop criteria**：任一 H NEGATIVE → drop feature group；全 NEGATIVE → revert Cycle 1

---

## §1 概念理解：ClairS Pipeline 完整時序（verified）

### Situation
團隊在 InterSubMod 研究中需要明確 ClairS 在哪個 step 用 longphase HP tag、哪個 step 不用，以及為什麼。

### Complication
ClairS 文件不完整：README 沒講 STEP 流程細節、NVTB 預設不在 user-facing docs。前期理解（如「pileup 也用 HP」「normal BAM 被 tag」）多處錯誤。

### Question
ClairS v0.4.4 / main HEAD `88f887a` 實際做什麼？

### Answer：完整資料流圖（source-verified）

```
[Inputs] tumor.bam (raw) + normal.bam (raw) + ref.fa
   │
   ├─[Pre-A] Clair3 normal  → normal_merge_output.vcf.gz
   ├─[Pre-B] Clair3 tumor   → tumor_merge_output.vcf.gz (供 STEP 4 用)
   │
   ├─[Pre-C] Select Het SNP  ★ NVTB 預設: 從 normal Clair3 抽 het SNP ★
   │         → clair3_output/vcf/{chr}.vcf  (het SNP 位點 + GT)
   │
   ├─[Pre-D] Phase the Tumor BAM
   │         longphase phase -s {het_snp.vcf} -b tumor.bam   ← 用 tumor BAM 證據
   │         → phased_output/tumor_phased_{chr}.vcf.gz  (含 PS tag)
   │
   ├─[Pre-E] Haplotag the Tumor BAM
   │         longphase haplotag -o phased_output/tumor_{chr}.bam \
   │                            -s tumor_phased.vcf.gz -b tumor.bam
   │         ★ 產出 tagged tumor BAM, 含 HP:i:1/2 (integer) ★
   │
   │  (Normal BAM 整個 pre-phase 階段都不動 — phase_normal=False 預設)
   │
   ├─[STEP 1] Extract Variant Candidates
   │          + raw tumor.bam + raw normal.bam
   │          → CANDIDATES_FILES (位點 list)
   │
   ├─[STEP 2] Pileup Head ─── 並行 ─────────────┐
   │  ├─ create_pair_tensor_pileup              │
   │  │   + raw tumor.bam                       │
   │  │   + raw normal.bam                      │
   │  │   + CANDIDATES_FILES                    │
   │  │   (channel_size = 34, 無 HP)            │
   │  └─ predict (pileup model)                 │
   │     → pileup.vcf                           │
   │                                            │
   ├─[STEP 3] Full-Alignment Head ─── 並行 ─┐ │
   │  ├─ create_pair_tensor                  │ │
   │  │   + tagged tumor.bam (HP:i 1/2)      │ │  ★ 唯一吃 HP 的 model ★
   │  │   + raw normal.bam                   │ │
   │  │   + CANDIDATES_FILES                 │ │
   │  │   (HAP_TYPE = {1:30, 0:60, 2:90})    │ │
   │  └─ predict (full model)                │ │
   │     → full_alignment.vcf                │ │
   │                                          │ │
   ├─[STEP 5] Merge & Sort ←─────────────────┴─┘
   │          sort_vcf(pileup.vcf + full_alignment.vcf)
   │          → output.vcf.gz  (ClairS 最終輸出)
   │
   └─[Downstream optional, ClairS 不知道]
        longphase-S somatic_haplotag
        + output.vcf.gz + tumor.bam + normal.bam
        → 重 tag tumor BAM: HP:z:1/2/1-1/2-1/3 (string)
        → (供 ISM / 下游分析；不再回饋 ClairS)
```

### 關鍵 source 行號（main HEAD `88f887a`）

| 行為 | 檔案 | 行號 |
|---|---|---|
| NVTB 預設邏輯 | `run_clairs` | 707-712 |
| `use_normal_bam_for_intermediate_phasing=False` (NVTB 推導) | `run_clairs` | 712 |
| Phase the Tumor BAM (`-b args.tumor_bam_fn`) | `run_clairs` | 1064-1078 |
| Haplotag the Tumor BAM (`-o phased_output/tumor_{chr}`) | `run_clairs` | 1103-1112 |
| `phase_normal=False` 預設 (line 765) | `run_clairs` | 763-765 |
| STEP 2 pileup 用 `args.tumor_bam_fn` (原始) | `run_clairs` | 1176 |
| STEP 3 full 用條件 `phased_output/tumor_*.bam if args.phase_tumor` | `run_clairs` | 1226-1232 |
| Pileup channel_size = 34 (無 HP) | `shared/param.py` | 79 |
| `phase_tumor` ONT 預設 True | `shared/param.py` | 64 |
| HAP_TYPE 編碼 {1:30, 0:60, 2:90} | `src/create_pair_tensor.py` | (top of file) |
| `if hap in '12'` substring check | `src/create_pair_tensor_pileup.py` | (in tensor build fn) |
| longphase-s HP:z 字串寫入 | `longphase-s/src/somatic_haplotag/SomaticHaplotagProcess.cpp` | 533 |

---

## §2 概念理解：NVTB 預設的設計哲學

### NVTB = Normal Variant + Tumor BAM

| 環節 | 來源 | 為何 |
|---|---|---|
| Het SNP **位置 + GT 標注** | normal Clair3 VCF | normal 是 diploid，het 偵測 F1 接近 1.0；tumor 有 LOH / heterogeneity 干擾 |
| Phase 階段 **read evidence** | **tumor BAM** | 用 tumor 自己 reads 推 hap 結構 → 後續 tag 內部自洽 |
| Haplotag target | **tumor BAM** | 後續 full-alignment 看的就是 tumor reads；normal 不需 phase |

### 為何 v0.1.7 後預設 NVTB（不是 TVTB）

`run_clairs:707` 註解：
> `# By default, we use HET SNP in normal VCF for phasing after v0.1.7`

設計理由（L3 reasoning）：
1. **避免 tumor germline call 被 somatic 干擾** — tumor 的 germline GT 可能因 LOH 錯標
2. **避免 systematic bias** — tumor het SNP 可能 enrich 在 LOH 邊界
3. **與 ISM self-phasing pitfall 一致** — phasing 訊號必須來自獨立 germline source（你 [[project_self_phasing_causal_chain_confirmed]] 的核心 lesson）

### 為何 Normal BAM 預設不 phase + tag

`run_clairs:763-765`：
```python
if args.phase_tumor is False:
    args.phase_normal = False
```
1. **訊號需求不對稱**：full-alignment 主要看 tumor reads 的 ALT 在 hap 上的分佈
2. **節省 30-50% wall time**
3. **避免「double-phase artifact」**：tumor 與 normal 各自 phase 後 PS block 不一致

---

## §3 概念理解：為何 Germline HP 有用 / Somatic-aware HP 沒用

### Mechanism

```
真 somatic SNV (clonal, single chromatid origin)
   → ALT reads 偏向某一 germline HP (e.g. 8 reads HP1, 0 reads HP2)
   → full-alignment model 看到 HP channel 集中 → 高信心 PASS

隨機 sequencing / mapping error
   → ALT reads 均勻分佈兩 HP (e.g. 4 reads HP1, 4 reads HP2)
   → HP channel 分散 → 低信心 reject
```

→ Germline HP 是**獨立先驗坐標系**，讓 model 檢驗 ALT clustering。

### 為何 somatic-aware HP (1-1/2-1/3) 餵 ClairS 沒用

| 性質 | germline HP | somatic-aware HP |
|---|---|---|
| 計算時機 | somatic call **之前**（純 germline phasing） | somatic call **之後**（依賴 caller VCF） |
| Circular 風險 | ❌ 無 | ✅ 高 |
| 狀態數 | 3 (1/2/0) | 5+ (1/2/0/1-1/2-1/3) |
| Read 的 intrinsic property | ✅ 是 | ❌ 否（依賴 caller 推論） |
| Pretrained model 看過 | ✅ | ❌ schema OOD |
| 新增訊息量 | 高（獨立 prior） | 低（衍生 posterior） |
| 編碼相容性 | HP:i:1/2（integer） | HP:z:1-1（string），ClairS `if hap in '12'` substring check 把它 silently drop |

→ 直接餵 somatic-aware HP 進 ClairS pretrained model = information loss（既存的 HP channel 變稀疏，因為 1-1/2-1 reads 從 "12" 集合裡消失）。

---

## §4 概念理解：Circular Dependency 的 nuance

### 修正過度泛化的論點

前期討論中曾過度 generalize 「circular dependency = bad」。這是錯的。

**真實判斷準則**：error correlation structure。

| Circular 安全 | Circular 危險 |
|---|---|
| 上游訊號**與下游 task 弱相關** | 上游 error **與下游 task 強相關** |
| 上游 error 是**隨機 noise**（會 averaging out） | 上游 error 有 **systematic bias**（會 amplify） |
| 有 **independent ground truth** 校準 | 完全依賴自己（self-fulfilling） |

### ML 領域 circular 是常規

- EM algorithm（統計學基石）
- Pseudo-labeling / Noisy Student（ImageNet 漲 2%）
- DeepVariant + WhatsHap iteration（germline F1 漲 0.5%）
- ClairS 自己（germline phase → phased BAM → somatic call 本質是 circular）

### 對 ClairS 上下文的 empirical 結論

[[project_self_phasing_causal_chain_confirmed]] 已實證：
- ClairS FP **不是純隨機**，在 cnLOH region 有 systematic bias
- 這些 bias 與 phasing pattern **高度相關**
- 所以「用 caller 結果 inform phasing」會 amplify → 62% LOH 消失

→ **在 ClairS + cnLOH-heavy tumor 樣本上，circular 確實 high-risk**，但這是 empirical observation，不是先驗 deduction。
→ **在低 cnLOH / 純 random error 樣本上，circular 可能 self-correct**。

---

## §5 InterSubMod 整合分析

### ISM 在 ClairS 生態系的位置

```
ClairS (upstream caller)
   ↓ output: somatic.vcf.gz
longphase-S (read-level annotation)
   ↓ output: tagged BAM (HP:z:1/2/1-1/2-1/3)
InterSubMod (read-level feature extraction)
   ↓ output: per-candidate features (methyl, HPFineNGroups, LOH zone, NG, CN)
ISM Re-ranker (post-filter, your Phase 2 Cycle 1 Global FP Filter)
   ↓ output: refined VCF (PASS/LowQual flag adjusted)
```

### ClairS gap → ISM opportunity 對照表

| ClairS 不知道 | InterSubMod 對應目標 | 當前 status |
|---|---|---|
| Somatic-somatic phasing | 目標 1: Sub-clonal modification phasing | [[project_loh_constrained_phasing_discovery]] PARTIAL POSITIVE |
| Sub-clonal grouping | 目標 3: Two-hit order / driver-passenger | [[project_p4_target3_depends_on_target1]] pending 目標 1 |
| Sub-clonal methylation + SNV co-segregation | 目標 2: Methylation-SNV phasing | [[project_loh_subclone_af_methylation_positive]] paired POSITIVE |
| Read-level epigenetic context | 目標 4-5: Pan-cancer epigenetic landscape | future |
| Read groupings 跨 candidates | (新探索) | candidate for Cycle 3 |

### Phase 2 Cycle 1 既有 baseline（不可丟）

[[project_phase2_cycle1_global_fp_filter]]：
- 10 features global LR (drop NumReads VIF=217, L2 C=1.0, τ=0.39)
- HCC1395 ΔF1=+0.02236 (9.24× v1.0)
- 3/4 H PASS
- Step 5c lost TP 81% rescued
- methylation 5th-rank
- ⭐3 strong

→ **本研究計劃 = Cycle 2，建立在 Cycle 1 baseline 上**，不是從零開始。

---

## §6 研究方向評估：5-Dim Credibility

### 用戶提的方向

> 把 longphase-S tag 的 somatic-aware HP BAM + 甲基 sub-clone 特徵 + read 分群狀況 餵新 model → 取代或補強 caller TP/FP 辨識

### Credibility 評分

| 維度 | 分數 | 理由 |
|---|---|---|
| 理論合理性 | ⭐⭐⭐⭐ | 資訊論支持：更多 informative features → ML 能 leverage |
| 新穎性 | ⭐⭐⭐⭐⭐ | 業界無人做 somatic-aware HP + sub-clonal methyl 整合 caller |
| 可實作性 | ⭐⭐ | Schema 擴展 + retrain + ground truth 都高 cost |
| 可重現性 | ⭐⭐ | 缺 read-level sub-clonal truth set |
| 失敗風險 | ⭐⭐⭐⭐⭐ HIGH | Circular dependency + label leakage + confound 累積 |

**綜合**：⭐⭐⭐ PROBE（值得試但需嚴格 pre-reg + ablation）

### 7 條質疑（每條都可能讓專案掛掉）

| # | Concern | Mitigation |
|---|---|---|
| G1 | Circular dependency — longphase-S HP:z 標注源自 caller VCF；FP 也被標 1-1/2-1 → 無鑑別力 | 跑 `--disableFilter` 包含 LowQual 位點才有 FP-only HP 訊號 |
| G2 | Label leakage — feature (HP:z) 與 label (truth VCF) 都間接源自 ClairS training set | LOSO + permutation test (random shuffle HP) |
| G3 | Methylation sub-clone feature 成熟度不足 — ISM 五目標多數未 fully validated | 先固化已 validate 的 ISM features，不要混 untested marker |
| G4 | Ground truth 取得困難 — 無 read-level sub-clonal truth | Simulation 或退到 binary TP/FP（SEQC2 有 label） |
| G5 | Schema 爆炸 + retrain cost ~5-10K GPU-hours | 不 retrain caller，改 post-filter（lightweight re-ranker ~1 GPU-hour） |
| G6 | Confound 累積（HP×caller_af、methyl×coverage、subclone×LOH） | 每 feature 過 `/auc-confound-guard`；one-at-a-time ablation |
| G7 | 與 ClairS 職責 overlap — Replace 不切實際 | 限 post-filter framing，不替代 caller |

### 為何 Post-filter framing 是唯一可行路線

| 風險 | Replace caller | Post-filter |
|---|---|---|
| G1 Circular | 高（caller 內部 amplify） | 低（只調整邊界 QUAL） |
| G2 Leakage | 高（training set overlap 大） | 中（hold-out 可控） |
| G3 Untested feature | 高（feature 進 caller 影響全部 candidate） | 低（feature 只影響 re-rank 邊界） |
| G4 Ground truth | 需 sub-clonal truth | 只需 binary TP/FP（有） |
| G5 Cost | ~5-10K GPU-hours | ~1 GPU-hour |
| G6 Confound | 難 ablation（端到端訓練） | 易 ablation（10 features 規模） |
| G7 Overlap | 直接衝突 | 明確分工 |

---

## §7 研究計劃：Phase 2 Cycle 2 設計

### Architecture（Post-hoc Re-ranker）

```
Input per-candidate:
  - ClairS output (QUAL, AF, INFO fields)
  - longphase-S annotation (HP:z statistics: count_11, count_21, count_3, etc.)
  - ISM features (methyl AUC, HPFineNGroups, LOH zone, CN, NG, ...)
  - (Cycle 3) sub-clone grouping label (待目標 1 完成)
       ↓
Lightweight model (Gradient Boosting / small MLP, ~10K params)
       ↓
Output:
  - Per-candidate re-ranking score
  - → Adjust ClairS QUAL or set new PASS/LowQual flag
  - → 寫出 InterSubMod-augmented VCF
```

### Pre-registration（必須在 cycle 啟動前定）

| ID | Hypothesis | Effect size 預估 | Stop criteria |
|---|---|---|---|
| H_C2A | 加 longphase-S HP:z annotation → HCC1395 ΔF1 > Cycle 1 baseline +0.002 | Cohen d ≥ 0.3 | 5/5 LOSO 4+ POSITIVE 才 PASS；否則 drop HP feature |
| H_C2B | 加 ISM methylation feature → ΔF1 > Cycle 1+HP-only baseline +0.001 | Cohen d ≥ 0.2 | 3/5 POSITIVE 才 PASS；否則 drop methyl |
| H_C2C | (待 Cycle 3) 加 sub-clone grouping → ΔF1 > +0.001 | Cohen d ≥ 0.2 | 依賴目標 1 完成才 ablation |

**Negative control（必跑）**：
- NC1: Permutation HP annotation → expect ΔF1 ≈ 0（若 > +0.005 表示 label leakage）
- NC2: Random methyl feature → expect ΔF1 ≈ 0

### 4-Week Pilot Timeline

| Week | Task | Output artifact |
|---|---|---|
| 1 | Reproduce Cycle 1 baseline + 萃取 longphase-S HP:z 統計（含 `--disableFilter`）+ ISM features | features.tsv per candidate |
| 2 | Ablation: Model A (baseline) / B (+HP) / C (+HP+methyl) / D (+HP+methyl+subclone if ready) | ablation_results.tsv (per cell line) |
| 3 | Negative controls (NC1 permutation HP, NC2 random methyl) | nc_results.tsv |
| 4 | Cross-sample LOSO 5-7 cell line + statistics + P5 EVALUATE report | cycle2_evaluation.json + final report |

### 預期結果區間

| 場景 | ΔF1 over Cycle 1 | Tier | 後續行動 |
|---|---|---|---|
| Optimistic | +0.005 ~ +0.010 incremental | ⭐3-4 | 寫 PI report，考慮 publication framing |
| Realistic | +0.001 ~ +0.003 marginal | ⭐2 | 報告 marginal positive，方向 NEGATIVE 但 method 有價值 |
| Pessimistic | ≈ 0 或 negative | ⭐1 NEGATIVE | 學到 feature 不獨立 → 知識仍有價值；revert Cycle 1 |

---

## §8 未來研究與驗證方向

### Tier 1：本計劃內（Phase 2 Cycle 2，2026 Q2）

1. **Cycle 2 ablation pilot**（4 週）— 上述計劃
2. **HP:z annotation 提取工具開發** — 從 longphase-s tagged BAM 萃取 per-candidate HP statistics 的 Python script（可考慮加入 InterSubMod tools/）
3. **`--disableFilter` longphase-s 重跑** — 確保 LowQual 位點也有 HP tag 供 ablation

### Tier 2：自然延伸（Phase 2 Cycle 3-4，2026 Q3）

4. **Sub-clone grouping integration**（待目標 1 完成）— [[project_loh_constrained_phasing_discovery]] 固化後納入 feature set
5. **Cross-pipeline generalization** — ClairS-TO / DeepSomatic 上同樣 framework 驗證（不只 ClairS）
6. **Liquid biopsy / low-purity samples** — Cohen d 可能在低 purity 樣本更大（caller 在低 purity 表現差，re-ranker 改善空間大）

### Tier 3：戰略級（2026 Q4 +）

7. **Joint caller + re-ranker training** — 若 post-filter 證明有效，考慮把 features 整合進新 caller training（需大量 GPU resource）
8. **Open-source release** — InterSubMod-augmented ClairS post-filter 作為 community tool（提升專案能見度）
9. **Paper 投稿** — 主軸：「Sub-clonal-aware post-filter for somatic variant calling in long-read tumor sequencing」

### Tier 4：開放問題（research frontier）

10. **Somatic-somatic phasing 標準化** — 業界目前無 caller 做 chromosome-scale somatic phasing；ISM 可定義 "somatic phasing F1" 評估標準
11. **Sub-clonal truth set 建立** — 與 simulation 工具（pbsim2 / badread）合作建 read-level sub-clonal benchmark
12. **Methylation as somatic marker** — 5mCG pattern 作為 somatic SNV 同 sub-clone 證據的可行性 (你目標 2)

### 各 Tier 對應 risk / reward

| Tier | Time horizon | Reward | Risk | 建議資源分配 |
|---|---|---|---|---|
| 1 | 2026 Q2 | ⭐3-4 verification | 中（4 G NEGATIVE 概率 30%） | 70% effort |
| 2 | 2026 Q3 | ⭐4 generalization | 中 | 20% effort |
| 3 | 2026 Q4+ | ⭐5 paper-scope | 高 | 10% effort (exploratory) |
| 4 | 2027+ | open-ended | 高 | spike 才探索 |

---

## §9 對話歷史證據鏈（Provenance）

本文件結論來自 2026-05-23 到 2026-05-24 的 8-round conversation（議題：ClairS / longphase-s 細節確認）：

| Round | Topic | Verdict |
|---|---|---|
| 1 | ClairS full vs pileup 是否吃 longphase HP | full 吃，pileup 不吃 |
| 2 | 餵 phased BAM 給 pileup 會不會更好 | NO — pretrained model schema mismatch |
| 3 | ClairS GitHub 最新版本 | v0.4.4 stable，main HEAD `88f887a` |
| 4 | longphase-s filter LowQual 位點 read tagging 行為 | filter 後不 tag 1-1/2-1（除非 `--disableFilter`） |
| 5 | 迭代式 R1→R2 ClairS+longphase-s flow 可行性 | NO — HP:z vs HP:i 編碼不相容 + retrain 成本 |
| 6 | ClairS 是否只把 longphase-s 當 downstream filter | YES — ClairS 不知 longphase-s 存在 |
| 7 | Germline HP 為何有用 / somatic HP 為何沒用 | Independent prior vs derived posterior; 編碼簡潔 vs 爆炸 |
| 8 | Circular 不合理嗎？ | 修正過度泛化；circular 是中性，重點在 error structure |
| 9 | tumor BAM 何時 tag、何時進 full | Pre-E 階段 tag（longphase haplotag）；STEP 3 進 full |
| 10 | ClairS 流程的精確 4-claim 對證 | NVTB phasing + parallel pileup/full + germline-only HP gap |
| 11 | v0.4.4 是最新 stable | 確認；main HEAD 多 6 commits 為 docs/bugfix |
| 12 | Caller-augmentation 研究方向評估 | PROBE; post-filter framing only; 7 條質疑 G1-G7 |

---

## §10 Sources（最新 main HEAD `88f887a`，2026-05-24 查驗）

### ClairS

- Repository: https://github.com/HKU-BAL/ClairS
- Main HEAD commit: `88f887a` (2026-04-16) "add Clair-skills"
- Latest stable: v0.4.4 (2025-11-18)
- Latest prerelease: v0.4.4-Zenodo (2026-03-12) — Zenodo DOI archive
- Key files:
  - `run_clairs` — 主 driver script
  - `shared/param.py` — pileup_channel_size / HAP_TYPE / phase_tumor 預設
  - `src/create_pair_tensor.py` — full-alignment tensor + HP channel
  - `src/create_pair_tensor_pileup.py` — pileup tensor (無 HP)
  - `clairs/model.py` — pileup BiGRU + full ResNet

### longphase-s

- Repository: https://github.com/CCU-Bioinformatics-Lab/longphase-s
- Local install: `/big8_disk/liaoyoyo2001/longphase-s/`
- Key files:
  - `src/somatic_haplotag/SomaticHaplotagProcess.cpp:533` — HP:z 字串寫入
  - `src/somatic_haplotag/SomaticVarCaller.cpp:1062-1230` — somaticFeatureFilter
  - `src/somatic_haplotag/SomaticVarCaller.cpp:2402-2412` — getSomaticFlag (filter → tag gate)

### Knowledge Base

- `/big8_disk/liaoyoyo2001/Knowledge/05_tools/longphase.md`
- `/big8_disk/liaoyoyo2001/Knowledge/05_tools/longphase-s.md`
- `/big8_disk/liaoyoyo2001/Knowledge/05_tools/longphase-to.md`
- `/big8_disk/liaoyoyo2001/Knowledge/05_tools/variant-callers.md` ⚠ drift — 仍寫 v0.4.0/v0.4.1, 需更新到 v0.4.4

### ISM Internal Projects

- `InterSubMod/research/methyl_augmented_filter_phase2/cycle1/` — Phase 2 Cycle 1 baseline
- `InterSubMod/research/methyl_augmented_filter_phase2/cycle4/loso_validation/` — LOSO methodology
- `InterSubMod/docs/concepts/2026/05/20260523_HKU_collab_主軸與細節_01.md` — HKU collab context

### Memory references

- [[project_self_phasing_causal_chain_confirmed]] — 62% LOH 消失 amplification 證據
- [[project_phase2_cycle1_global_fp_filter]] — Cycle 1 ΔF1=+0.02236 baseline
- [[project_methyl_filter_pilot_marginal_positive]] — methyl single-feature marginal
- [[project_loh_constrained_phasing_discovery]] — 目標 1 PARTIAL POSITIVE
- [[project_loh_subclone_af_methylation_positive]] — 目標 2 paired POSITIVE
- [[project_getvote_per_read_concept]] — per-read framework foundation

### Literature 對照

- Cooke DP et al. 2021 *Nat Biotechnol* — Octopus haplotype-aware Bayesian
- Poplin R et al. 2018 *Nat Biotechnol* — DeepVariant + WhatsHap iteration
- Lee D-H. 2013 *ICML Workshop* — Pseudo-Label
- Xie Q et al. 2020 *CVPR* — Noisy Student
- ClairS paper: Zheng Z et al. 2023 *Nat Comp Sci* (待查 KB drift fix)

---

## §11 文件 metadata

- **Created**: 2026-05-24
- **Status**: concept (pre-cycle, awaiting `/cycle-init` for Phase 2 Cycle 2)
- **Tier**: ⭐3 design-stage
- **Next review**: After Phase 2 Cycle 2 Week 1 baseline reproduction
- **Drift triggers**:
  - ClairS 新版發布（v0.4.5+ → 需重查 phasing 邏輯是否改）
  - longphase-s 新版發布
  - InterSubMod 目標 1 fully validated（→ 升級到 Cycle 3 計劃）
  - Phase 2 Cycle 2 NEGATIVE/POSITIVE 出爐（→ 對應更新 verdict）
- **Owner**: InterSubMod Research
