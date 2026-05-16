# Self-Phasing V6 Binary Production Tag & PI Errata Package

> **狀態**：executing
> **建立日期**：2026-05-16
> **專案目錄**：`research/selfphasing_v6_production/`
> **Parent plan**：`~/.claude/plans/tender-pondering-blossom.md`
> **主軸報告**：`InterSubMod/docs/reports/validated/2026/05/20260508_Self_Phasing_完整觀察整合報告_01.md`

## 背景與動機

Self-Phasing 整合主軸 5/8-5/10 完成 5-commit 鏈，揭露 V5 Layer 1.5 在 germline-absent 區域繼承 priority bug 偏移（4.19:1 偏 HP1，與 baseline 完全相同）。V6 binary patch（V5 phasing + V3F-style hp=33 保守處理）已實作並在 chr19+全基因組三向驗證，marker coverage 超越 V3F 與 V5。

2026-05-15 進一步透過 HCC1395 全 chr V3F/V5/V6 三向 ISM 比較揭露 V5 over-promote 量化證據：Inner LOH NG=2 region V5=8,136 (+60% over V3F 5,064)，V6 修補回 V3F 水準 (5,353)。V5/V3F top cell ratio 達 5.95× 集中在 cross_het_inv bucket — Layer 1.5 機制只在 somatic-fallback heterozygous reads 作用。

**Tier 1 必須前置目標**：V6 production tag finalize（COLO829 ISM 補完 + binary commit hash 寫 manifest + git tag `v6-prod-*`）。完成後 unlock thread_d_paper Archive TO rerun + PI errata package 一併打包。

## 假說

### H1：V6 修補 V5 Layer 1.5 over-promote

**陳述**：V6 binary patch（V5 phasing + V3F-style hp=33）在 germline-absent + LOH 區修補 V5 過度推斷，同時保留 V5 在 phasing-weak 區的優勢，整體 marker coverage 不降反升。

**前提條件**：
- V3F/V5/V6 三向 ISM 比較完成（HCC1395 全 chr 2026-05-15 ✅）
- V5/V3F ratio 集中在 cross_het bucket（驗證 over-promote 機制定位）

**已知 Confound**：
- Binary commit hash drift（V5 baseline 從 938f0df 退到更早版本可能不一致）
- Sample-specific over-promote 幅度（5 樣本中 magnitude 可能差異）

**驗證標準**：
- **Positive**：Inner LOH NG=2 V6 ≈ V3F ±10% AND V5/V3F >1.5× 集中 cross_het bucket
- **Negative**：V6 region count > V5 OR V5/V3F ratio 分散在所有 bucket

**Verdict (2026-05-15)**：**supported** — HCC1395 V3F=5,064 / V5=8,136 / V6=5,353; ratio 5.95× isolated to cross_het_inv

### H2：V6 marker coverage > V3F/V5 with 0 critical regression

**陳述**：V6 binary 在 7 樣本 ISM 跑完後，marker coverage 平均提升 +5% 以上於 V3F，且 caller F1 vs SEQC2 truth set 全部 7 樣本無 regression。

**前提條件**：
- 7 樣本 V6 ISM 完整跑完（5/7 done; COLO829 pending）
- Caller F1 vs SEQC2 truth set computed for V3F/V5/V6 (HCC1395 已確認 0.7166 三版相同)

**已知 Confound**：
- COLO829 ONT R10 無 methylation（sample-specific limitation；marker coverage 計算可能受影響）
- Caller F1 metric 口徑（per longphase-to paper §4.3 V_H/V_L post-filter vs caller-level F1）

**驗證標準**：
- **Positive**：V6 marker coverage +5%+ over V3F，n>=5/7 同方向；caller F1 無樣本退步
- **Negative**：≥1 樣本 marker coverage 退步 OR caller F1 退步 >0.005

**Verdict (2026-05-15)**：**partial** — 5/7 done; HCC1395 +9% confirmed; 4 樣本 extending

### H3：V6 production tag unlocks downstream + PI errata package

**陳述**：V6 production tag finalize 後同時完成 (a) Archive TO 7 樣本 rerun unlock（thread_d_paper Tier 2）+ (b) PI Report 4-29 errata 5 條打包 + V6 sign-off 一併 written email。

**前提條件**：
- H1 + H2 supported
- Binary commit hash finalize + git tag created
- Manifest.yaml updated with v6 commit
- 5 條 errata content drift 已收斂（含 5/9 paired audit 新加 E5）

**已知 Confound**：
- PI ack timing（lab meeting vs written email 時程不同）
- Errata content drift（V6 sign-off 過程可能發現新 errata）

**驗證標準**：
- **Positive**：Git tag v6-prod-* 存在 + 5 errata + V6 sign-off email delivered
- **Negative**：PI 在 sign-off 前要求 new audit（T4.3 workload 擴大）

## 方法

### 數據來源

| 數據集 | 路徑 | 描述 | 使用欄位 |
|--------|------|------|---------|
| V3F/V5/V6 master | `research/v6_bam_tpfp_hp_loh_cn/step1_master_three_way.tsv` | 三向 ISM (HCC1395) | Region_ID, HPFineN_*, binary_version, TP/FP |
| 4 樣本 V6 ISM | `output/canonical/{H1437,H2009,HCC1954,HCC1937}/V6/` | 2026-05-15 擴展 | 同 HCC1395 |
| COLO829 V6 (pending) | `output/canonical/COLO829/V6/` | Archive TO rerun 待補 | — |
| Self-phasing 整合報告 | `docs/reports/validated/2026/05/20260508_*.md` | V6 patch motivation source | — |
| PI errata companion | `docs/reports/validated/2026/05/20260509_PI_Report_4_29_Errata_01/` | 5 條 errata candidate | — |

### 分析步驟

```
Step 1: COLO829 V6 ISM 補完 (Archive TO rerun + KDE-corrected)
        → 驗證: ls output/canonical/COLO829/V6/ 存在 + step1_master 加 COLO829 row
Step 2: 7-sample marker coverage 比較 (V3F vs V5 vs V6)
        → 驗證: research/v6_bam_tpfp_hp_loh_cn/step6_marker_coverage_7sample.tsv 含 7 樣本 ΔCoverage
Step 3: 7-sample caller F1 比較 (V3F vs V5 vs V6 vs SEQC2)
        → 驗證: research/v6_bam_tpfp_hp_loh_cn/step7_caller_f1_7sample.tsv F1 三版 within 0.005
Step 4: Binary commit hash 寫 manifest.yaml + git tag v6-prod-{YYYYMMDD}
        → 驗證: git tag --list 顯示 v6-prod-*; manifest.yaml v6_binary_commit 填入
Step 5: PI errata 5 條 + V6 sign-off written email 打包
        → 驗證: errata companion 5 條 final review + email draft saved
```

### 統計方法

- **Per-sample comparison**: paired t-test on marker coverage delta (V6-V3F)
- **Cross-sample direction**: sign test (n>=5/7 same direction)
- **Caller F1 equivalence**: |F1_V6 - F1_V3F| < 0.005 across all samples
- **No significance correction**: descriptive comparison, not hypothesis test

## 可行性評估

| 因素 | 評估 |
|------|------|
| 數據可用性 | △ 5/7 V6 已跑（HCC1395 全 chr + 4 樣本）；COLO829 待補 |
| 計算資源 | ✓ COLO829 V6 ISM ~2 hr；後續分析 <1 hr |
| 與已有結論衝突 | ✓ 不衝突 — Self-phasing 整合報告已 align |

## 已知風險

1. **COLO829 ONT R10 limitation**：無 methylation 訊號可能讓 V6 marker coverage 計算偏低；緩解：V6 sample-level gating 標註「優先 ONT_5mCG 或 5mCG+5hmCG」
2. **V5 baseline binary 不一致**：HEAD 938f0df 與更早 commit V5 行為可能差異；緩解：manifest 寫明 V5 baseline = 938f0df
3. **PI ack 時程**：errata + V6 sign-off 需 PI 同意才能 close；緩解：先 written email（24hr 內 ack），不開 lab meeting
4. **Errata content drift**：5/9 paired audit 已加 E5 候選；緩解：errata companion 已 patch（2026-05-09），週 3 finalize 前再 review

## Tier 1.2 V6 Production Tag Workflow

```
[Day 1-2] COLO829 V6 ISM 補完
   ↓
[Day 3] 7-sample marker coverage + caller F1 比較 (step6/7 docs)
   ↓
[Day 4] Binary commit hash 寫 manifest.yaml
   ↓
[Day 4] git tag v6-prod-{YYYYMMDD} (Hard Gate — 需 user 明確 ack)
   ↓
[Day 5] PI errata 5 條 + V6 sign-off written email draft
   ↓
[Day 5] User review email → send (Hard Gate — 需 user 確認)
   ↓
[Tier 1.2 完成] thread_d_paper Tier 2 unlock
```

## 相關檔案

- `InterSubMod/docs/reports/validated/2026/05/20260508_Self_Phasing_完整觀察整合報告_01.md` — 主軸 commit chain 報告
- `InterSubMod/docs/reports/validated/2026/05/20260505_self_phasing_V5_data_provenance_audit_01.md` — V5 audit
- `InterSubMod/docs/reports/validated/2026/05/20260509_PI_Report_4_29_Errata_01/` — errata companion
- `InterSubMod/research/v5_provenance_followup/` — T1/T1.2/T1.2-F1 audit material
- `InterSubMod/research/paired_priority_bug_audit/` — 5/9 paired Step A+C+D
- `InterSubMod/research/v6_bam_tpfp_hp_loh_cn/` — 2026-05-15 multi-agent fan-out
- `~/.claude/plans/tender-pondering-blossom.md` — parent plan (Tier 1.2 + T4.3)

## Plan tender-pondering-blossom Tier 對照

- **Tier 1.2**: V6 production tag finalize (本專案核心)
- **Tier 4.3**: PI errata sign-off (本專案 Step 5)
- **下游 dependency**: thread_d_paper Tier 2 Archive TO rerun（等本專案 git tag 完成）
