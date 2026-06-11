<!--
建立時間: 2026-05-19
agent: main session Coordinator
status: complete
report_class: cycle step verdict (Cycle 3 Step 1)
audience: Coordinator / Step 2 panel-survey agent / cycle 3 verdict synthesizer
scope: 5 samples (HCC1395 + H1437 + H2009 + HCC1954 + HCC1937) gated filter reproduce
parent_plan: InterSubMod/research/methyl_augmented_filter_phase2/cycle3/00_PLAN.md §Step 1
predecessor_cycle: 20260519_phase2_cycle2_cross_sample_negative
verdict: H_C3_1 PASS + H_C3_3 PASS (H_C3_2 pending Step 2)
last_verified: 2026-05-19
-->

# Phase 2 Cycle 3 Step 1 — Caller-F1-headroom-aware Gate Reproduce

> **Verdicts**:
> - **H_C3_1 PASS** — qualifying mean ΔF1 = **+0.01499** ≥ +0.01 (1.50× threshold, n=2)
> - **H_C3_3 PASS** — non-qualifying max |drift| = **0.000000** ≪ 0.001 (trivial PASS by construction)
> - **H_C3_2** — PENDING Step 2 (panel expansion 至 n≥4 待 inventory)

---

## 0. TL;DR

| Sample | Caller F1 | FP density | Gate | Action | Gated ΔF1 |
|--------|--:|--:|---|---|--:|
| HCC1395 | 0.7166 | 14.0% | ✅ PASS | apply re-fit τ=0.39 | **+0.02236** |
| HCC1937 | 0.3692 | 16.2% | ✅ PASS | apply re-fit τ=0.24 | **+0.00761** |
| H1437 | 0.8670 | 1.1% | ❌ skip | preserve caller F1 | 0.00000 |
| H2009 | 0.8863 | 1.0% | ❌ skip | preserve caller F1 | 0.00000 |
| HCC1954 | 0.8385 | 3.4% | ❌ skip | preserve caller F1 | 0.00000 |

**Qualifying mean ΔF1** = (+0.02236 + +0.00761) / 2 = **+0.01499** (PASS H_C3_1 threshold +0.01 = 1.50× ribbon)

**Non-qualifying drift** = 0.00000 (PASS H_C3_3 cap 0.001 = trivial preservation)

**Naive 5-sample mean** (for comparison vs cycle 2 transfer) = +0.005994 = nominal positive across n=5；但 cycle 2 transfer mode 同樣 5-sample mean = -0.094（H_C1_5 NEGATIVE）→ gate rule 改善 16×。

---

## 1. Method

### 1.1 Gate config (`cycle3_caller_f1_gate.json`)

```json
{
  "rule_name": "caller_f1_headroom_aware_gate",
  "gate_conditions": {
    "caller_F1_max": 0.80,
    "FP_density_min": 0.10,
    "combine": "AND"
  },
  "qualifying_action": "apply_per_sample_refit",
  "non_qualifying_action": "skip_filter"
}
```

Per-sample 邏輯：

```
IF caller_F1 < 0.80 AND FP_density > 0.10:
    delta_F1 = cycle 2 re-fit ΔF1 (per-sample 5-fold OOF swept τ)
ELSE:
    delta_F1 = 0 (skip)
```

### 1.2 Data source

- Per-sample re-fit ΔF1: `cycle2/data/cycle2_cross_sample_delta_f1.tsv` (mode=refit rows)
- Caller F1 + region-level TP/FP totals: same file (caller_F1, tp_total_used, fp_total_used columns)
- FP density: fp_total / (tp_total + fp_total)

### 1.3 Reproducibility

- Script: `cycle3/scripts/apply_gated_filter.py`
- Wall clock: ~2 sec
- Deterministic — re-fit ΔF1 from cycle 2 B3' (PRIMARY_SEED=42, lbfgs solver)

---

## 2. Results

### 2.1 Gate evaluation per sample

| Sample | Caller F1 < 0.80? | FP density > 0.10? | Gate AND | Action |
|--------|:-:|:-:|:-:|:-:|
| HCC1395 | ✅ 0.7166 | ✅ 0.140 | ✅ PASS | apply re-fit |
| HCC1937 | ✅ 0.3692 | ✅ 0.162 | ✅ PASS | apply re-fit |
| H1437 | ❌ 0.8670 | ❌ 0.011 | ❌ | skip |
| H2009 | ❌ 0.8863 | ❌ 0.010 | ❌ | skip |
| HCC1954 | ❌ 0.8385 | ❌ 0.034 | ❌ | skip |

**觀察**: H1437/H2009/HCC1954 兩個 gate 都不通過（高 F1 + 低 FP density 同向）；HCC1395/HCC1937 兩個 gate 都通過（低 F1 + 高 FP density 同向）。Gate 條件冗餘但保守。

### 2.2 Gated ΔF1 vs cycle 2 transfer/re-fit comparison

| Sample | Cycle 2 transfer fixed | Cycle 2 transfer swept | Cycle 2 re-fit | **Cycle 3 gated** |
|--------|:--|:--|:--|:--|
| HCC1395 | +0.02232 | +0.02246 | +0.02236 | **+0.02236** |
| H1437 | -0.03597 | -0.00744 | ~0 | **0.00000** |
| H2009 | -0.00450 | -0.00085 | ~0 | **0.00000** |
| HCC1954 | **-0.37744** ⚠ | -0.16972 | +0.00095 | **0.00000** |
| HCC1937 | -0.07068 | -0.02859 | +0.00761 | **+0.00761** |
| **5-sample mean** | -0.094 | -0.038 | +0.006 | **+0.006** |

**Gate rule 效果**:
- 移除 HCC1954 transfer 災難 (-0.377)
- 移除 H1437/H2009 微負（gate 確認本就不該應用）
- 保留 HCC1395 + HCC1937 全部 uplift

### 2.3 Verdict thresholds

| Pre-reg | Computed | Threshold | Verdict |
|---|---|---|---|
| H_C3_1 qualifying mean ΔF1 ≥ +0.01 | **+0.01499** | +0.01 | **PASS** (1.50×) |
| H_C3_3 non-qualifying max |drift| < 0.001 | **0.000000** | 0.001 | **PASS** (trivial) |
| H_C3_2 panel n≥4 Wilcoxon p<0.05 | n=2 only | n≥4 | **PENDING** Step 2 |

---

## 3. Verdict Interpretation

### 3.1 H_C3_1 PASS (qualifying subset effective)

n=2 (HCC1395 + HCC1937) qualifying samples mean ΔF1 +0.01499 = 1.50× Cohen small ribbon (+0.01)。
- **HCC1395 +0.02236**: 主導 contribution，與 cycle 1 stored bit-exact match
- **HCC1937 +0.00761**: 小幅 uplift，cycle 2 唯一 re-fit positive 新樣本

**Caveat**: n=2 too small for Wilcoxon (exact min p=0.50 with n=2)。Pre-reg 用 strict mean threshold 而非 statistical test，但這意味著 H_C3_1 不能單獨升 tier；需配 H_C3_2 panel expansion。

### 3.2 H_C3_3 PASS (trivial by construction)

Skip 策略 ΔF1=0 by definition → drift=0。**但這是 trivial PASS**：
- 不證明 gate rule 在 high-F1 樣本「真實對抗 filter degradation」
- 證明 gate rule **不會在 high-F1 樣本造成 F1 regression**（cycle 2 transfer mode 災難得到完全 mitigation）

對照 cycle 2 transfer mode `apply unconditionally`:
- HCC1954 -0.37744 → 0.0 (gate 修補 +0.377)
- H1437 -0.03597 → 0.0 (gate 修補 +0.036)
- H2009 -0.00450 → 0.0 (gate 修補 +0.005)

→ Gate rule **避免 0.418 total F1 loss** across 3 high-F1 samples。

### 3.3 H_C3_2 PENDING (panel expansion 必要)

n=2 qualifying 不足升 ⭐4。Step 2 必要：
- Survey InterSubMod archive 所有 caller F1 < 0.80 samples
- 找 ≥2 個新 low-F1 candidates (達 n≥4)
- 跑 ISM + cycle 1 re-fit + Wilcoxon

候選: HCC1395 (✓ done) / HCC1937 (✓ done) / COLO829 (truth set chenhan112 權限 pending) / 其他 SEQC2 derivatives / TCGA pairs。

---

## 4. Mechanism Confirmation

Cycle 2 已揭示 caller-F1-headroom-bounded mechanism；Step 1 提供進一步 evidence：

1. **AND gate 冗餘但保守正確**: 5 samples 中 caller F1 + FP density gate 完全同向（無 mixed case）→ 兩條件 highly correlated（FP density 隨 caller F1 上升而下降）。生產規則用 AND 比 OR 安全。

2. **HCC1395 + HCC1937 共同特性**: caller F1 ≤ 0.72 + FP density ≥ 14% + TP/FP ratio ≤ 7。Filter 在這個 "headroom corner" 有 work。

3. **High-F1 樣本不適用 filter ≠ filter wrong**: Cycle 1 filter 在 HCC1395 ⭐3 internal valid + H_C1_6 BAM-invariant；fail 只發生在 caller 已接近 ceiling 的樣本（已無 filter 可改善的 FP）。

4. **Methylation 仍 5th-rank**: Gate rule 不改 cycle 1 LR 內部結構；caller_af (+3.44) + LOH (+1.46) + Coverage (+1.27) + NG (+1.07) + HPFineF (+0.75) 排序保留。

---

## 5. Files

```
cycle3/
├── 00_PLAN.md
├── cycle3_caller_f1_gate.json                      (gate config, Step 1a output)
├── cycle3_step1_findings.md                        (this report)
├── scripts/
│   └── apply_gated_filter.py                       (Step 1b)
├── data/
│   └── cycle3_step1_gated_delta_f1.tsv             (5 samples × gate evaluation)
├── figures/
│   └── cycle3_step1_gated_delta_f1.png             (bar chart 5-sample gated ΔF1)
└── intermediate/
    └── cycle3_step1_summary.json                   (machine-readable verdict)
```

---

## 6. Hand-off to Step 2

### 6.1 Step 2 trigger conditions

| Condition | Status |
|---|---|
| H_C3_1 PASS qualifying subset | ✅ |
| H_C3_3 PASS non-qualifying preservation | ✅ |
| Gate rule formalized + tested | ✅ |
| Panel expansion 必要 | ✅ (n=2 不夠) |

→ **Step 2 GO**：panel survey 找 ≥2 個新 low-F1 candidates。

### 6.2 Step 2 candidate priority

1. **COLO829** (highest priority): truth set 已存在 (0600 chenhan112 權限 pending)，V6 BAM tagging 可能已 done in v6_5sample_extension archive
2. **TCGA paired tumor/normal samples** (medium): InterSubMod archive 可能有 callbase F1 < 0.80 的 candidates
3. **SEQC2 derivatives** (medium): HCC1395 truth set 同源候選
4. **HCC1954 ancestral** (low): 已知 caller F1 0.8385，邊緣值，可能不過 gate

### 6.3 Step 2 fallback (if n<4)

依 plan v0.1 §10：「Panel survey yields n=2 only → 限 qualifying subset 結論 ⭐3，不退 NEGATIVE」。

---

## 7. Reproducibility

```bash
cd /big7_disk/liaoyoyo2001/InterSubMod
python3 research/methyl_augmented_filter_phase2/cycle3/scripts/apply_gated_filter.py
```

Wall clock: ~2 sec. Deterministic — uses cycle 2 B3' refit ΔF1 verbatim, no new fits.

**Bit-exact reproducibility**: HCC1395 gated ΔF1 = +0.0223580973952435 matches cycle 2 re-fit drift 0.

---

**End of Cycle 3 Step 1 Findings**
